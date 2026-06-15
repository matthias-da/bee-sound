#!/usr/bin/env Rscript
# 03_mspb_modeling.R — core modelling on MSPB (the unblocked cohort)
#  C1: does predictive skill collapse from random k-fold -> leave-hive-out -> leave-apiary-out?
#  C2: does audio predict the bottom (mic) box better / more transferably than total FoB?
#  C3: robust hierarchical fit (rlmer vs lmer) + split-conformal coverage under apiary transfer.
#  Leave-apiary-out (Côté <-> Dubuc) is the cross-apiary transfer (C4) in miniature.

suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(ranger) })
have_rob <- requireNamespace("robustlmm", quietly = TRUE)
set.seed(1)
out_dir <- here::here("analysis", "output")
H <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))

bands <- paste0("f_", 1:16)
feat  <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation",
           bands, "t_in_mean", "rh_in_mean")
core  <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation",
           "t_in_mean", "rh_in_mean")                 # interpretable subset for the LMM
audio_only <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation", bands)

d <- H$mspb_tab %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes)) %>%
  dplyr::mutate(hive = factor(hive), yard = factor(yard)) %>%
  dplyr::select(hive, yard, eval, n_boxes, fob_total, brood1, all_of(feat)) %>%
  tidyr::drop_na(all_of(c(feat, "brood1")))

cat(sprintf("modelling rows: %d | hives: %d | apiaries: %s\n",
            nrow(d), nlevels(d$hive), paste(levels(d$yard), collapse = ", ")))
cat("fob_total by apiary:\n"); print(round(tapply(d$fob_total, d$yard, mean), 1))

metrics <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat); y <- y[ok]; yhat <- yhat[ok]
  c(RMSE = sqrt(mean((y - yhat)^2)), MAE = mean(abs(y - yhat)),
    R2 = 1 - sum((y - yhat)^2) / sum((y - mean(y))^2), corr = cor(y, yhat))
}
mkfolds <- function(d, scheme) switch(scheme,
  random = sample(rep_len(1:5, nrow(d))),
  hive   = as.integer(d$hive),
  apiary = as.integer(d$yard))

cv_predict <- function(d, scheme, engine, target = "fob_total", feats = feat) {
  fold <- mkfolds(d, scheme); pr <- rep(NA_real_, nrow(d))
  for (f in sort(unique(fold))) {
    tr <- d[fold != f, ]; te <- d[fold == f, ]
    if (nrow(te) == 0 || nrow(tr) < 10) next
    if (engine == "rf") {
      m <- ranger(reformulate(feats, target), data = tr, num.trees = 500, seed = 1)
      pr[fold == f] <- predict(m, te)$predictions
    } else if (engine == "lmm") {
      fit <- tryCatch(lmer(reformulate(c(core, "(1|hive)"), target), data = tr,
                           REML = FALSE, control = lmerControl(calc.derivs = FALSE)),
                      error = function(e) NULL)
      if (!is.null(fit))
        pr[fold == f] <- tryCatch(as.numeric(predict(fit, te, allow.new.levels = TRUE)),
                                  error = function(e) NA_real_)
    }
  }
  pr
}
row_metrics <- function(label_tbl, y, yhat) dplyr::bind_cols(label_tbl, as_tibble_row(metrics(y, yhat)))

## ---- Experiment 1: generalization gap (C1) + apiary transfer (C4) ----
schemes <- c("random", "hive", "apiary")
res <- map_dfr(schemes, function(s) dplyr::bind_rows(
  row_metrics(tibble(scheme = s, model = "RandomForest (audio+T/H)"),
              d$fob_total, cv_predict(d, s, "rf")),
  row_metrics(tibble(scheme = s, model = "Hier LMM (core + 1|hive)"),
              d$fob_total, cv_predict(d, s, "lmm"))))
cat("\n===== Experiment 1: fob_total prediction by CV regime =====\n")
print(as.data.frame(res %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

## ---- Experiment 2: bottom-box mediator (C2) ----
e2 <- map_dfr(c("hive", "apiary"), function(s) dplyr::bind_rows(
  row_metrics(tibble(scheme = s, target = "total FoB"),
              d$fob_total, cv_predict(d, s, "rf", "fob_total", audio_only)),
  row_metrics(tibble(scheme = s, target = "bottom box (mic)"),
              d$brood1,    cv_predict(d, s, "rf", "brood1",    audio_only))))
cat("\n===== Experiment 2: audio-only prediction — total vs bottom (mic) box =====\n")
print(as.data.frame(e2 %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

## ---- Experiment 3: robust fit + conformal under apiary transfer (C3) ----
cat("\n===== Experiment 3: robust hierarchical fit + conformal under transfer =====\n")
f_core <- reformulate(c(core, "(1|hive)"), "fob_total")
m_lmer <- lmer(f_core, data = d, REML = FALSE, control = lmerControl(calc.derivs = FALSE))
vc <- as.data.frame(VarCorr(m_lmer)); icc <- vc$vcov[vc$grp == "hive"] / sum(vc$vcov)
cat(sprintf("ICC(hive) for fob_total (LMM): %.3f\n", icc))
if (have_rob) {
  m_rob <- robustlmm::rlmer(f_core, data = d)
  cat("fixed effects — classical (lmer) vs robust (rlmer):\n")
  print(round(cbind(lmer = fixef(m_lmer), rlmer = fixef(m_rob)), 3))
} else cat("(robustlmm unavailable — skipping rlmer)\n")

# 90% split-conformal, calibrated on hive-disjoint subset of the training apiary, tested on the other
conf <- map_dfr(levels(d$yard), function(test_ap) {
  tr <- d[d$yard != test_ap, ]; te <- d[d$yard == test_ap, ]
  hv <- sample(unique(tr$hive)); cal_h <- hv[seq_len(max(2, floor(length(hv) * 0.3)))]
  fit <- ranger(reformulate(feat, "fob_total"), data = tr[!tr$hive %in% cal_h, ], num.trees = 500, seed = 1)
  cal <- tr[tr$hive %in% cal_h, ]
  q90 <- quantile(abs(cal$fob_total - predict(fit, cal)$predictions), 0.90, names = FALSE)
  pr  <- predict(fit, te)$predictions
  tibble(train_apiary = setdiff(levels(d$yard), test_ap), test_apiary = test_ap,
         n_test = nrow(te), nominal = 0.90,
         emp_coverage = round(mean(abs(te$fob_total - pr) <= q90), 3),
         mean_width = round(2 * q90, 2),
         RMSE_transfer = round(metrics(te$fob_total, pr)["RMSE"], 2))
})
cat("90% split-conformal (grouped calibration) under leave-apiary-out:\n")
print(as.data.frame(conf), row.names = FALSE)

## ---- save ----
theme_set(theme_minimal(base_size = 11))
p <- res %>%
  ggplot(aes(factor(scheme, levels = schemes), R2, fill = model)) +
  geom_col(position = "dodge") + geom_hline(yintercept = 0, linewidth = .3) +
  labs(title = "MSPB: predictive R² collapses under honest cross-validation",
       subtitle = "random k-fold (leaky) → leave-hive-out → leave-apiary-out (Côté↔Dubuc transfer)",
       x = "CV regime", y = expression(R^2), fill = NULL)
ggsave(file.path(out_dir, "fig_mspb_cv_gap.png"), p, width = 8, height = 5, dpi = 130)
write_csv(res, file.path(out_dir, "mspb_cv_results.csv"))
cat("\nWrote fig_mspb_cv_gap.png + mspb_cv_results.csv\n")

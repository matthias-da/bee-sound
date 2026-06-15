#!/usr/bin/env Rscript
# 04_mspb_refinements.R
#  C2 (structure-aware): does n_boxes (management-known, no inspection) carry the FoB signal that
#     audio cannot, and does it absorb the audio effect? -> tests the bottom-box-mediator mechanism
#     without the saturation confound of the naive test.
#  Conformal (firmed up): in-distribution vs apiary-transfer coverage, both directions, many seeds.

suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(ranger) })
out_dir <- here::here("analysis", "output")
H <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))

bands <- paste0("f_", 1:16)
audio_only <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation", bands)
feat_all   <- c(audio_only, "t_in_mean", "rh_in_mean")

d <- H$mspb_tab %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes)) %>%
  dplyr::mutate(hive = factor(hive), yard = factor(yard)) %>%
  dplyr::select(hive, yard, n_boxes, fob_total, brood1, all_of(feat_all)) %>%
  tidyr::drop_na(all_of(c(feat_all, "brood1")))

metrics <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat); y <- y[ok]; yhat <- yhat[ok]
  c(RMSE = sqrt(mean((y - yhat)^2)),
    R2 = 1 - sum((y - yhat)^2) / sum((y - mean(y))^2), corr = cor(y, yhat))
}
mkfolds <- function(d, scheme) switch(scheme,
  hive = as.integer(d$hive), apiary = as.integer(d$yard))
cv_rf <- function(d, scheme, feats, target = "fob_total") {
  fold <- mkfolds(d, scheme); pr <- rep(NA_real_, nrow(d))
  for (f in sort(unique(fold))) {
    tr <- d[fold != f, ]; te <- d[fold == f, ]
    if (nrow(te) == 0 || nrow(tr) < 10) next
    m <- ranger(reformulate(feats, target), data = tr, num.trees = 500, seed = 1)
    pr[fold == f] <- predict(m, te)$predictions
  }
  pr
}
row_m <- function(lbl, y, yhat) dplyr::bind_cols(lbl, as_tibble_row(metrics(y, yhat)))

## ===== Refinement 1: structure-aware C2 =====
set.seed(1)
fsets <- list("n_boxes only" = "n_boxes",
              "audio only"    = audio_only,
              "audio + n_boxes" = c(audio_only, "n_boxes"))
c2 <- map_dfr(c("hive", "apiary"), function(s)
  imap_dfr(fsets, ~ row_m(tibble(scheme = s, predictors = .y),
                          d$fob_total, cv_rf(d, s, .x))))
cat("===== C2: predicting total FoB (RF) — what carries the signal? =====\n")
print(as.data.frame(c2 %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

cat("\n--- mediation: does the audio effect shrink when structure (n_boxes) is added? ---\n")
m1 <- lmer(fob_total ~ scale(hive_power) + (1 | hive), d, REML = FALSE,
           control = lmerControl(calc.derivs = FALSE))
m2 <- lmer(fob_total ~ scale(hive_power) + n_boxes + (1 | hive), d, REML = FALSE,
           control = lmerControl(calc.derivs = FALSE))
cat(sprintf("hive_power coef (per SD):  audio-only model = %.2f  |  + n_boxes = %.2f  (%.0f%% shrink)\n",
            fixef(m1)["scale(hive_power)"], fixef(m2)["scale(hive_power)"],
            100 * (1 - fixef(m2)["scale(hive_power)"] / fixef(m1)["scale(hive_power)"])))

cat("\n--- within-colony-size association of audio (hive_power) with total FoB ---\n")
for (nb in 2:3) {
  sub <- d[d$n_boxes == nb, ]
  cat(sprintf("  n_boxes=%d (n=%d): cor(hive_power, FoB) = %.2f\n",
              nb, nrow(sub), cor(sub$hive_power, sub$fob_total)))
}

## ===== Refinement 2: conformal under transfer, firmed up =====
rf_fit <- function(tr) ranger(reformulate(feat_all, "fob_total"), data = tr, num.trees = 300, seed = 1)
q90    <- function(fit, cal) quantile(abs(cal$fob_total - predict(fit, cal)$predictions), .90, names = FALSE)
cw     <- function(fit, q, te) { pr <- predict(fit, te)$predictions
                                 c(cov = mean(abs(te$fob_total - pr) <= q), width = 2 * q) }
aps <- levels(d$yard)
one_rep <- function(seed) {
  set.seed(seed); out <- list()
  for (ap in aps) {                                  # in-distribution: hive-grouped 60/20/20 within apiary
    da <- d[d$yard == ap, ]; hv <- sample(unique(da$hive)); n <- length(hv)
    tr_h <- hv[seq_len(floor(.6 * n))]
    cal_h <- hv[(floor(.6 * n) + 1):floor(.8 * n)]
    te_h <- hv[(floor(.8 * n) + 1):n]
    if (length(te_h) < 1 || length(cal_h) < 1) next
    fit <- rf_fit(da[da$hive %in% tr_h, ]); r <- cw(fit, q90(fit, da[da$hive %in% cal_h, ]), da[da$hive %in% te_h, ])
    out[[length(out) + 1]] <- tibble(setting = paste0("in-dist: ", ap), cov = r["cov"], width = r["width"])
  }
  for (ap in aps) {                                  # transfer: train apiary -> test the other
    tr <- d[d$yard == ap, ]; te <- d[d$yard != ap, ]
    hv <- sample(unique(tr$hive)); n <- length(hv)
    tr_h <- hv[seq_len(floor(.7 * n))]; cal_h <- hv[(floor(.7 * n) + 1):n]
    fit <- rf_fit(tr[tr$hive %in% tr_h, ]); r <- cw(fit, q90(fit, tr[tr$hive %in% cal_h, ]), te)
    out[[length(out) + 1]] <- tibble(setting = paste0("transfer: ", ap, "->", setdiff(aps, ap)),
                                     cov = r["cov"], width = r["width"])
  }
  dplyr::bind_rows(out)
}
B <- 80
conf <- map_dfr(seq_len(B), one_rep)
conf_summary <- conf %>% dplyr::group_by(setting) %>%
  dplyr::summarise(reps = dplyr::n(), mean_cov = mean(cov), sd_cov = sd(cov), mean_width = mean(width), .groups = "drop")
cat("\n===== Conformal (90% nominal): in-distribution vs apiary transfer, over", B, "seeds =====\n")
print(as.data.frame(conf_summary %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

## ===== save =====
theme_set(theme_minimal(base_size = 11))
p <- conf %>%
  dplyr::mutate(kind = dplyr::if_else(grepl("^in-dist", setting), "in-distribution", "apiary transfer")) %>%
  ggplot(aes(setting, cov, fill = kind)) +
  geom_hline(yintercept = .90, linetype = 2) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.55) +
  geom_jitter(aes(colour = kind), width = 0.18, height = 0, size = 0.7, alpha = 0.5,
              show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white",
               colour = "black") +
  coord_flip() +
  scale_fill_manual(values = c("in-distribution" = "#0170B0", "apiary transfer" = "#D75E01")) +
  scale_colour_manual(values = c("in-distribution" = "#0170B0", "apiary transfer" = "#D75E01")) +
  labs(x = NULL, y = "empirical coverage", fill = NULL,
       caption = "points = 80 per-split values (jittered); diamond = mean; dashed = nominal 0.90")
ggsave(file.path(out_dir, "fig_mspb_conformal_transfer.png"), p, width = 8.5, height = 4.5, dpi = 130)
write_csv(c2, file.path(out_dir, "mspb_c2_structure.csv"))
write_csv(conf_summary, file.path(out_dir, "mspb_conformal_summary.csv"))
cat("\nWrote fig_mspb_conformal_transfer.png + mspb_c2_structure.csv + mspb_conformal_summary.csv\n")

#!/usr/bin/env Rscript
# 11_review_analyses.R
# Answers several of Jürgen's review points with concrete analyses; writes outputs + a log to
# analysis/output/ (supplementary material candidates). Input: analysis/output/mspb_harmonized.rds.
suppressPackageStartupMessages({
  library(tidyverse); library(lme4); library(performance); library(see)
  library(leaps); library(GGally); library(flextable); library(broom.mixed)
})
out_dir <- here::here("analysis", "output"); set.seed(1)
ctrl  <- lmerControl(calc.derivs = FALSE)
bands <- paste0("f_", 1:16)
core4 <- c("hive_power","audio_density","audio_density_ratio","density_variation")
audio20 <- c(core4, bands); th <- c("t_in_mean","rh_in_mean")

d <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))$mspb_tab %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes)) %>%
  dplyr::mutate(hive = factor(hive), yard = factor(yard),
                upper = fob_total - brood1) %>%
  tidyr::drop_na(all_of(c(audio20, th, "brood1")))

log <- file.path(out_dir, "review_analyses.txt")
con <- file(log, open = "wt"); sink(con); sink(con, type = "message")
cat("n =", nrow(d), "evaluations,", nlevels(droplevels(d$hive)), "hives\n")

## ---- (R citation) ----
cat("\n===== R citation() =====\n"); print(citation(), style = "text")
cat("\nR version:", R.version.string, "\n")

## ---- (Table 2 95% CIs, standardised) ----
cat("\n===== Table 2 standardised coefficients with 95% CIs (Wald) =====\n")
ds <- d %>% dplyr::mutate(dplyr::across(all_of(c(core4, th)), ~ as.numeric(scale(.))),
                          y = as.numeric(scale(fob_total)))
ms <- lmer(reformulate(c(core4, th, "(1|hive)"), "y"), ds, REML = TRUE, control = ctrl)
ci <- confint(ms, method = "Wald", parm = "beta_")
tab2 <- data.frame(term = rownames(summary(ms)$coef),
                   est = round(fixef(ms), 3),
                   lo  = round(ci[names(fixef(ms)), 1], 3),
                   hi  = round(ci[names(fixef(ms)), 2], 3),
                   t   = round(summary(ms)$coef[, "t value"], 2))
print(tab2, row.names = FALSE)
write_csv(tab2, file.path(out_dir, "table2_with_CIs.csv"))

## ---- (Multicollinearity: does it harm the TEMPERATURE estimate? McElreath leg-length Q) ----
cat("\n===== Collinearity & the temperature estimate =====\n")
cat("VIFs (raw-scale Table-2 model):\n")
m_raw <- lmer(reformulate(c(core4, th, "(1|hive)"), "fob_total"), d, REML = FALSE, control = ctrl)
print(performance::check_collinearity(m_raw))
cat("\nIs TEMPERATURE's estimate destabilised by adding the 4 collinear audio features?\n")
ds2 <- d %>% dplyr::mutate(dplyr::across(all_of(c(core4,th)), ~as.numeric(scale(.))))
mt_s  <- lmer(fob_total ~ t_in_mean + rh_in_mean + (1|hive), ds2, REML=TRUE, control=ctrl)
mta_s <- lmer(reformulate(c(core4,th,"(1|hive)"),"fob_total"), ds2, REML=TRUE, control=ctrl)
cmp <- data.frame(
  model = c("temp+humidity only", "+ 4 audio features"),
  temp_coef = c(fixef(mt_s)["t_in_mean"], fixef(mta_s)["t_in_mean"]),
  temp_SE   = c(sqrt(diag(vcov(mt_s)))["t_in_mean"], sqrt(diag(vcov(mta_s)))["t_in_mean"]))
print(cmp, row.names = FALSE)
cat("-> if temp_coef and temp_SE barely move, temperature (VIF~1.3) is NOT harmed by the audio collinearity.\n")

## ---- (Best-subset regression on the 20 audio features) ----
cat("\n===== Best-subset regression (leaps::regsubsets) on the 20 audio features =====\n")
rs <- regsubsets(reformulate(audio20, "fob_total"), data = d, nvmax = 20, really.big = TRUE)
sm <- summary(rs)
best <- data.frame(size = seq_along(sm$adjr2),
                   adjR2 = round(sm$adjr2, 3), BIC = round(sm$bic, 1), Cp = round(sm$cp, 1))
print(best, row.names = FALSE)
cat(sprintf("\nBest adjR2 over ALL audio subsets: %.3f (size %d) -- in-sample, no random effect.\n",
            max(sm$adjr2), which.max(sm$adjr2)))
# benchmark: box count alone (in-sample, plain lm for comparability)
cat(sprintf("Box count alone (lm, in-sample) R2 = %.3f ; adjR2 = %.3f\n",
            summary(lm(fob_total ~ n_boxes, d))$r.squared,
            summary(lm(fob_total ~ n_boxes, d))$adj.r.squared))
cat("-> even the best audio subset is compared against a single box-count predictor; out-of-sample\n")
cat("   (Table 3, leave-apiary-out) is the decisive comparison and already shows audio 0.17 vs box 0.71.\n")

## ---- (Localisation regression check_model + causal note) ----
cat("\n===== Localisation regression: hive_power ~ bottom + upper + (1|hive) =====\n")
# Pre-scale the variables into data columns rather than using scale() INSIDE the formula:
# with scale() in the formula, performance::check_predictions() returns the raw (unscaled)
# response while simulate() returns standardised values, so the posterior-predictive-check
# plots observed and simulated on different scales and looks broken. Pre-scaling fixes it; the
# fitted coefficients are identical either way.
dm <- d %>% dplyr::filter(n_boxes >= 2) %>%
  dplyr::mutate(z_hp = as.numeric(scale(hive_power)),
                z_bottom = as.numeric(scale(brood1)),
                z_upper  = as.numeric(scale(upper)))
ml <- lmer(z_hp ~ z_bottom + z_upper + (1 | hive), dm, REML = FALSE, control = ctrl)
print(round(summary(ml)$coefficients, 3))
cat("check_normality:"); print(performance::check_normality(ml))
cat("check_heteroscedasticity:"); print(performance::check_heteroscedasticity(ml))
sink(type = "message"); sink(); close(con)
cat("Wrote review_analyses.txt\n")

## ---- figures: ggpairs of the 4 audio features; check_model of localisation ----
gp <- GGally::ggpairs(d %>% dplyr::select(all_of(core4)),
                      title = "Pairwise relations of the four interpretable audio features (MSPB)")
ggsave(file.path(out_dir, "fig_supp_ggpairs_audio.png"), gp, width = 8, height = 7, dpi = 130)
cat("Wrote fig_supp_ggpairs_audio.png\n")

grDevices::png(file.path(out_dir, "fig_supp_checkmodel_localisation.png"), width=1500, height=1150, res=120)
print(plot(performance::check_model(ml, verbose = FALSE)))
grDevices::dev.off()
cat("Wrote fig_supp_checkmodel_localisation.png\n")

## ---- data-overview flextable (head + tail) for the supplement ----
show_cols <- c("hive","yard","eval","n_boxes","fob_total","brood1","upper",
               "hive_power","audio_density_ratio","t_in_mean","rh_in_mean")
ht <- d %>% dplyr::select(any_of(show_cols)) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 2)))
preview <- dplyr::bind_rows(utils::head(ht, 4),
                            ht[0,] %>% dplyr::add_row(),
                            utils::tail(ht, 4)) %>%
  dplyr::mutate(dplyr::across(everything(), ~ifelse(is.na(.x), "...", as.character(.x))))
ft <- flextable(preview) %>% theme_booktabs() %>% fontsize(size = 8, part = "all") %>%
  autofit() %>% set_caption("MSPB analysis table — first and last four rows (… = elision).")
flextable::save_as_html(ft, path = file.path(out_dir, "supp_data_overview.html"))
tryCatch(flextable::save_as_image(ft, path = file.path(out_dir, "supp_data_overview.png")),
         error = function(e) cat("(png export needs webshot2; html written)\n"))
cat("Wrote supp_data_overview.html\n")
cat("\nDone (11_review_analyses).\n")

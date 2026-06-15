#!/usr/bin/env Rscript
# 09_lmm_interpret.R
# Fit and INTERPRET the hierarchical linear mixed model used in the paper (the interpretable
# counterpart to the random forest), on the same MSPB modelling set as Table 1 (script 03).
# Reports: raw-scale fixed effects (FoB units) with SE/t; standardised effects (per SD, comparable);
# variance components (sigma_hive, sigma_resid) and conditional ICC; Nakagawa marginal/conditional R^2.
# Input: analysis/output/mspb_harmonized.rds
suppressPackageStartupMessages({ library(tidyverse); library(lme4) })
out_dir <- here::here("analysis", "output")
set.seed(1)
ctrl <- lmerControl(calc.derivs = FALSE)
bands <- paste0("f_", 1:16)
feat  <- c("hive_power","audio_density","audio_density_ratio","density_variation", bands, "t_in_mean","rh_in_mean")
core  <- c("hive_power","audio_density","audio_density_ratio","density_variation","t_in_mean","rh_in_mean")

d <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))$mspb_tab %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes)) %>%
  dplyr::mutate(hive = factor(hive)) %>%
  dplyr::select(hive, n_boxes, fob_total, brood1, all_of(feat)) %>%
  tidyr::drop_na(all_of(c(feat, "brood1")))
cat(sprintf("modelling set: n=%d evaluations, %d hives\n", nrow(d), nlevels(d$hive)))

## ---- raw-scale model (FoB units) ----
m <- lmer(reformulate(c(core, "(1|hive)"), "fob_total"), d, REML = TRUE, control = ctrl)
cat("\n--- RAW-scale fixed effects (response in frames of bees) ---\n")
print(round(summary(m)$coefficients, 3))
vc <- as.data.frame(VarCorr(m))
v_hive <- vc$vcov[vc$grp == "hive"]; v_res <- vc$vcov[vc$grp == "Residual"]
cat(sprintf("\nsigma_hive = %.2f frames | sigma_resid = %.2f frames | conditional ICC = %.3f\n",
            sqrt(v_hive), sqrt(v_res), v_hive / (v_hive + v_res)))
varF <- var(as.vector(model.matrix(m) %*% fixef(m)))
cat(sprintf("Nakagawa R2: marginal (fixed only) = %.3f | conditional (fixed+hive) = %.3f\n",
            varF / (varF + v_hive + v_res), (varF + v_hive) / (varF + v_hive + v_res)))

## ---- standardised model (predictors AND response z-scored: effects per SD, comparable) ----
ds <- d %>% dplyr::mutate(dplyr::across(all_of(core), ~ as.numeric(scale(.))), y = as.numeric(scale(fob_total)))
ms <- lmer(reformulate(c(core, "(1|hive)"), "y"), ds, REML = TRUE, control = ctrl)
cat("\n--- STANDARDISED fixed effects (per 1 SD of predictor; response in SD units) ---\n")
print(round(summary(ms)$coefficients, 3))
cat("\nDone.\n")

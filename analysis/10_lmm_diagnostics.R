#!/usr/bin/env Rscript
# 10_lmm_diagnostics.R
# Model-assumption checks for the hierarchical LMMs used in the paper (Table 2 model from script 09
# and the minimal hive-power mediator model). Produces performance::check_model() diagnostic panels
# plus numeric checks (normality, heteroscedasticity, collinearity/VIF, outliers, singularity) and
# saves everything to analysis/output/.
# Input: analysis/output/mspb_harmonized.rds (from 02).
suppressPackageStartupMessages({
  library(tidyverse); library(lme4); library(performance); library(see)
})
out_dir <- here::here("analysis", "output")
set.seed(1)
ctrl  <- lmerControl(calc.derivs = FALSE)
bands <- paste0("f_", 1:16)
feat  <- c("hive_power","audio_density","audio_density_ratio","density_variation", bands,
           "t_in_mean","rh_in_mean")
core  <- c("hive_power","audio_density","audio_density_ratio","density_variation",
           "t_in_mean","rh_in_mean")

d <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))$mspb_tab %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes)) %>%
  dplyr::mutate(hive = factor(hive)) %>%
  dplyr::select(hive, n_boxes, fob_total, brood1, all_of(feat)) %>%
  tidyr::drop_na(all_of(c(feat, "brood1")))
cat(sprintf("modelling set: n=%d evaluations, %d hives\n", nrow(d), nlevels(d$hive)))

## helper: save a performance::check_model() panel robustly to PNG
save_check <- function(model, file, title) {
  grDevices::png(file.path(out_dir, file), width = 1500, height = 1150, res = 120)
  on.exit(grDevices::dev.off())
  print(plot(performance::check_model(model, verbose = FALSE)) +
          patchwork::plot_annotation(title = title))
  invisible(NULL)
}

## numeric assumption checks dumped to a text log
sink_log <- file.path(out_dir, "lmm_diagnostics.txt")
con <- file(sink_log, open = "wt"); sink(con); sink(con, type = "message")
report_checks <- function(model, label) {
  cat("\n==================== ", label, " ====================\n", sep = "")
  cat("\n-- fixed effects --\n");        print(round(summary(model)$coefficients, 3))
  cat("\n-- R2 (Nakagawa) --\n");         print(performance::r2(model))
  cat("\n-- ICC --\n");                   print(performance::icc(model))
  cat("\n-- collinearity (VIF) --\n");    print(performance::check_collinearity(model))
  cat("\n-- normality of residuals --\n"); print(performance::check_normality(model))
  cat("\n-- normality of random effects --\n")
  print(tryCatch(performance::check_normality(model, effects = "random"),
                 error = function(e) conditionMessage(e)))
  cat("\n-- heteroscedasticity --\n");    print(performance::check_heteroscedasticity(model))
  cat("\n-- singularity --\n");           cat(performance::check_singularity(model), "\n")
  cat("\n-- influential observations (outliers) --\n")
  print(summary(performance::check_outliers(model)))
}

## ---- Model A: the Table-2 LMM (core acoustic + T/H, raw FoB scale) ----
mA <- lmer(reformulate(c(core, "(1|hive)"), "fob_total"), d, REML = TRUE, control = ctrl)
report_checks(mA, "Model A: fob_total ~ core acoustic + T/H + (1|hive)  [Table 2]")

## ---- Model B: minimal hive-power mediator model (Section 'structure') ----
mB <- lmer(fob_total ~ hive_power + (1 | hive), d, REML = TRUE, control = ctrl)
report_checks(mB, "Model B: fob_total ~ hive_power + (1|hive)  [minimal]")

sink(type = "message"); sink(); close(con)
cat("Wrote lmm_diagnostics.txt\n")

## ---- diagnostic panels ----
save_check(mA, "fig_lmm_checkmodel_tableA.png", "LMM diagnostics — Table 2 model (core acoustic + T/H + 1|hive)")
cat("Wrote fig_lmm_checkmodel_tableA.png\n")
save_check(mB, "fig_lmm_checkmodel_minimal.png", "LMM diagnostics — minimal hive-power model (+1|hive)")
cat("Wrote fig_lmm_checkmodel_minimal.png\n")

cat("\nDone. Diagnostics in analysis/output/: lmm_diagnostics.txt + 2 check_model panels.\n")

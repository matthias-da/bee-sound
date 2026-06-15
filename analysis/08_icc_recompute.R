#!/usr/bin/env Rscript
# 08_icc_recompute.R
# Like-for-like ICC of total frames of bees (FoB), replacing the apples-to-oranges
# "0.83 (UrBAN) vs 0.21 (MSPB)" headline. We report the UNCONDITIONAL ICC (intercept-only
# random-intercept model) of fob_total on the FULL inspection/evaluation set of each dataset,
# with a parametric-bootstrap 95% CI and the sample sizes. We also show (i) why the old UrBAN
# 0.83 was inflated (it used the n=17 late-season audio-matched subset) and (ii) the MSPB
# CONDITIONAL ICC (residual between-hive share after the acoustic + T/H fixed effects) = old 0.21.
# Inputs: urban/UrBAN/data/annotations/inspections_2021.csv ; analysis/output/mspb_harmonized.rds
#         (and analysis/output/urban_matched_table.rds from 07, for documentation).
suppressPackageStartupMessages({ library(tidyverse); library(lme4) })
set.seed(1)
out_dir <- here::here("analysis", "output")
ctrl <- lmerControl(calc.derivs = FALSE)

icc_of <- function(m) { v <- as.data.frame(VarCorr(m)); v$vcov[v$grp != "Residual"][1] / sum(v$vcov) }
report <- function(label, form, data, nsim = 1000) {
  m <- lmer(form, data = data, REML = TRUE, control = ctrl)
  icc <- icc_of(m)
  bb <- tryCatch(suppressWarnings(bootMer(m, icc_of, nsim = nsim, seed = 1, type = "parametric")),
                 error = function(e) NULL)
  ci <- if (!is.null(bb)) quantile(bb$t, c(.025, .975), na.rm = TRUE) else c(NA, NA)
  cat(sprintf("%-38s ICC=%.3f  95%% CI [%.2f, %.2f]  (n=%d obs, %d hives)\n",
              label, icc, ci[1], ci[2], nobs(m), ngrps(m)[1]))
  invisible(icc)
}

cat("===== UrBAN (full 2021 inspections) =====\n")
insp <- read_csv(here::here("urban", "UrBAN", "data", "annotations", "inspections_2021.csv"),
                 show_col_types = FALSE) %>%
  dplyr::transmute(hive = factor(`Tag number`), n_boxes = as.integer(`Colony Size`),
            fob1 = as.numeric(`Fob 1st`), fob2 = as.numeric(`Fob 2nd`), fob3 = as.numeric(`Fob 3rd`)) %>%
  dplyr::mutate(fob_total = rowSums(cbind(fob1, fob2, fob3), na.rm = TRUE),
         occ = fob_total / (10 * n_boxes)) %>%
  dplyr::filter(fob_total > 0, !is.na(n_boxes))
report("UrBAN full: FoB (unconditional)",      fob_total ~ 1 + (1|hive), insp)
report("UrBAN full: occupancy (unconditional)", occ ~ 1 + (1|hive), insp)
mt <- file.path(out_dir, "urban_matched_table.rds")
if (file.exists(mt))
  report("UrBAN audio-matched subset: FoB",     fob_total ~ 1 + (1|hive),
         readRDS(mt) %>% dplyr::mutate(hive = factor(hive)))

cat("\n===== MSPB (full evaluations) =====\n")
mspb <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))$mspb_tab %>%
  dplyr::filter(!is.na(fob_total), fob_total > 0) %>%
  dplyr::mutate(hive = factor(hive), occ = fob_total / (10 * n_boxes))
report("MSPB full: FoB (unconditional)",        fob_total ~ 1 + (1|hive), mspb)
report("MSPB full: occupancy (unconditional)",  occ ~ 1 + (1|hive), mspb)
core <- c("hive_power","audio_density","audio_density_ratio","density_variation","t_in_mean","rh_in_mean")
mc <- mspb %>% dplyr::filter(dplyr::if_all(all_of(core), is.finite))
report("MSPB conditional (resid. after audio+T/H)", reformulate(c(core,"(1|hive)"), "fob_total"), mc, nsim = 500)
cat("\nDone.\n")

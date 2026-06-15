#!/usr/bin/env Rscript
# 05_mspb_upperbox.R — localisation test for the bottom-box-vantage mechanism (completes lead C2).
# Does the in-hive microphone (bottom brood box) sense bees in its box but not those above it?

suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(ranger) })
out_dir <- here::here("analysis", "output")
set.seed(1)
H <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))
bands <- paste0("f_", 1:16)
audio_only <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation", bands)

d <- H$mspb_tab %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes), !is.na(brood1)) %>%
  dplyr::mutate(hive = factor(hive), yard = factor(yard),
         upper = fob_total - brood1) %>%        # frames above the microphone box
  tidyr::drop_na(all_of(audio_only)) %>%
  dplyr::filter(n_boxes >= 2)                          # an upper box must exist

cat(sprintf("rows (multi-box): %d | hives: %d\n", nrow(d), nlevels(droplevels(d$hive))))
cat(sprintf("target variances:  bottom(brood1)=%.2f  upper=%.2f  total=%.2f\n",
            var(d$brood1), var(d$upper), var(d$fob_total)))
cat(sprintf("mean frames:       bottom=%.1f  upper=%.1f  total=%.1f\n",
            mean(d$brood1), mean(d$upper), mean(d$fob_total)))

## 1. within-hive (group-mean-centered) correlation of audio energy with each compartment
wc <- d %>% dplyr::group_by(hive) %>%
  dplyr::mutate(dplyr::across(c(hive_power, brood1, upper, fob_total), ~ . - mean(.))) %>% dplyr::ungroup()
cat("\n[1] within-hive correlation of hive_power with:\n")
cat(sprintf("    bottom (brood1): %.3f\n    upper:           %.3f\n    total:           %.3f\n",
            cor(wc$hive_power, wc$brood1), cor(wc$hive_power, wc$upper), cor(wc$hive_power, wc$fob_total)))

## 2. mixed model: which compartment drives the measured audio energy?
m <- lmer(scale(hive_power) ~ scale(brood1) + scale(upper) + (1 | hive), d,
          REML = FALSE, control = lmerControl(calc.derivs = FALSE))
cat("\n[2] hive_power ~ bottom + upper (+ 1|hive), standardized fixed effects:\n")
print(round(summary(m)$coefficients, 3))

## 3. audio-only leave-apiary-out prediction of each compartment (corr is scale-free; R2 for context)
cv_metrics <- function(target) {
  fold <- as.integer(d$yard); pr <- rep(NA_real_, nrow(d))
  for (f in sort(unique(fold))) {
    tr <- d[fold != f, ]; te <- d[fold == f, ]
    mm <- ranger(reformulate(audio_only, target), data = tr, num.trees = 500, seed = 1)
    pr[fold == f] <- predict(mm, te)$predictions
  }
  y <- d[[target]]; ok <- is.finite(pr)
  c(corr = cor(y[ok], pr[ok]), R2 = 1 - sum((y[ok] - pr[ok])^2) / sum((y[ok] - mean(y[ok]))^2))
}
cat("\n[3] audio-only leave-apiary-out prediction:\n")
res <- map_dfr(c("brood1", "upper", "fob_total"), function(t) {
  r <- cv_metrics(t); tibble(target = t, corr = r["corr"], R2 = r["R2"])
})
print(as.data.frame(res %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
write_csv(res, file.path(out_dir, "mspb_upperbox.csv"))
cat("\nWrote mspb_upperbox.csv\n")

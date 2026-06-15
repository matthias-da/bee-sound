#!/usr/bin/env Rscript
# 07_urban_modeling.R
# UrBAN modelling arm (the high-ICC regime) + cross-dataset transfer with MSPB.
#  - reproduce the within-UrBAN random vs leave-hive-out collapse on UrBAN's own hand-crafted features
#  - structure-confound (box count vs audio) and ICC on UrBAN
#  - cross-dataset transfer: train on one dataset, test on the other, on the shared 20-feature space
#    (features z-scored WITHIN each dataset to remove sensor-calibration offset/scale)
# Inputs: analysis/output/urban_audio_features_2021.rds (from 06); inspections_2021.csv;
#         analysis/output/mspb_harmonized.rds (from 02).

suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(ranger); library(lubridate) })
set.seed(1)
out_dir <- here::here("analysis", "output")
feat20 <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation", paste0("f_", 1:16))

## ---- UrBAN labels (re-derived from inspections_2021.csv; self-contained) ----
insp <- read_csv(here::here("urban", "UrBAN", "data", "annotations", "inspections_2021.csv"),
                 show_col_types = FALSE) %>%
  dplyr::transmute(hive = as.character(`Tag number`), date = ymd(Date), n_boxes = as.integer(`Colony Size`),
            fob1 = as.numeric(`Fob 1st`), fob2 = as.numeric(`Fob 2nd`), fob3 = as.numeric(`Fob 3rd`)) %>%
  dplyr::mutate(fob_total = rowSums(cbind(fob1, fob2, fob3), na.rm = TRUE),
         share_bottom = fob1 / fob_total) %>%
  dplyr::filter(fob_total > 0, !is.na(n_boxes))

## ---- attach audio features: mean over recordings within +/-3 days of each inspection ----
feat <- readRDS(file.path(out_dir, "urban_audio_features_2021.rds")) %>%
  dplyr::mutate(hive = as.character(hive)) %>% dplyr::rename(rdate = date)
agg <- insp %>% dplyr::mutate(.iid = dplyr::row_number()) %>% dplyr::select(.iid, hive, idate = date) %>%
  dplyr::inner_join(feat %>% dplyr::select(hive, rdate, all_of(feat20)), by = "hive",
             relationship = "many-to-many") %>%
  dplyr::filter(abs(as.integer(rdate - idate)) <= 3) %>%
  dplyr::group_by(.iid) %>%
  dplyr::summarise(dplyr::across(all_of(feat20), ~ mean(.x, na.rm = TRUE)), n_rec = dplyr::n(), .groups = "drop")
urb <- insp %>% dplyr::mutate(.iid = dplyr::row_number()) %>% dplyr::left_join(agg, by = ".iid") %>%
  dplyr::filter(!is.na(hive_power)) %>% dplyr::mutate(hive = factor(hive))
cat(sprintf("UrBAN matched inspections with audio: %d / %d | hives: %d | median recs/cell: %d\n",
            nrow(urb), nrow(insp), nlevels(urb$hive), median(urb$n_rec)))

metrics <- function(y, yhat) { ok <- is.finite(y) & is.finite(yhat); y <- y[ok]; yhat <- yhat[ok]
  c(RMSE = sqrt(mean((y - yhat)^2)), R2 = 1 - sum((y - yhat)^2) / sum((y - mean(y))^2), corr = cor(y, yhat)) }
rowm <- function(lbl, y, yhat) dplyr::bind_cols(lbl, as_tibble_row(metrics(y, yhat)))

## ---- ICC on the matched UrBAN set ----
icc <- tryCatch({ m <- lmer(fob_total ~ 1 + (1 | hive), urb, REML = FALSE,
                            control = lmerControl(calc.derivs = FALSE))
  v <- as.data.frame(VarCorr(m)); v$vcov[v$grp == "hive"] / sum(v$vcov) }, error = function(e) NA)
cat(sprintf("ICC(hive) for fob_total on matched UrBAN set: %.3f\n", icc))

## ---- within-UrBAN CV gap (random vs leave-hive-out) ----
cv_rf <- function(d, scheme, feats, target = "fob_total") {
  fold <- if (scheme == "random") sample(rep_len(1:5, nrow(d))) else as.integer(d$hive)
  pr <- rep(NA_real_, nrow(d))
  for (f in sort(unique(fold))) { tr <- d[fold != f, ]; te <- d[fold == f, ]
    if (nrow(te) == 0 || nrow(tr) < 8) next
    m <- ranger(reformulate(feats, target), data = tr, num.trees = 500, seed = 1)
    pr[fold == f] <- predict(m, te)$predictions }
  pr
}
cat("\n===== UrBAN: fob_total prediction (random vs leave-hive-out) =====\n")
ucv <- dplyr::bind_rows(
  rowm(tibble(scheme = "random",          predictors = "audio (20)"), urb$fob_total, cv_rf(urb, "random", feat20)),
  rowm(tibble(scheme = "leave-hive-out",  predictors = "audio (20)"), urb$fob_total, cv_rf(urb, "hive",   feat20)),
  rowm(tibble(scheme = "random",          predictors = "box count"),  urb$fob_total, cv_rf(urb, "random", "n_boxes")),
  rowm(tibble(scheme = "leave-hive-out",  predictors = "box count"),  urb$fob_total, cv_rf(urb, "hive",   "n_boxes")))
print(as.data.frame(ucv %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

## ---- structure: coefficient shrink + within-stratum corr (UrBAN) ----
cat("\n----- UrBAN structure -----\n")
m1 <- lmer(fob_total ~ scale(hive_power) + (1 | hive), urb, REML = FALSE, control = lmerControl(calc.derivs = FALSE))
m2 <- lmer(fob_total ~ scale(hive_power) + n_boxes + (1 | hive), urb, REML = FALSE, control = lmerControl(calc.derivs = FALSE))
cat(sprintf("hive_power coef (per SD): audio-only %.2f  |  + n_boxes %.2f\n",
            fixef(m1)["scale(hive_power)"], fixef(m2)["scale(hive_power)"]))
for (nb in sort(unique(urb$n_boxes))) { s <- urb[urb$n_boxes == nb, ]
  if (nrow(s) > 3) cat(sprintf("  n_boxes=%d (n=%d): cor(hive_power, FoB)=%.2f\n", nb, nrow(s), cor(s$hive_power, s$fob_total))) }

## ---- cross-dataset transfer (shared features z-scored within each dataset) ----
mspb <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))$mspb_tab %>%
  dplyr::filter(!is.na(fob_total), fob_total > 0, !is.na(hive_power)) %>%
  dplyr::select(fob_total, n_boxes, all_of(feat20)) %>% tidyr::drop_na(all_of(feat20))
zscore <- function(df, cols) { df %>% dplyr::mutate(dplyr::across(all_of(cols), ~ as.numeric(scale(.x)))) }
mz <- zscore(mspb, feat20); uz <- zscore(urb %>% dplyr::select(fob_total, n_boxes, all_of(feat20)), feat20)

transfer <- function(train, test, feats) {
  m <- ranger(reformulate(feats, "fob_total"), data = train, num.trees = 500, seed = 1)
  metrics(test$fob_total, predict(m, test)$predictions)
}
cat("\n===== cross-dataset transfer (RF; features z-scored within dataset) =====\n")
tr <- dplyr::bind_rows(
  rowm(tibble(train = "MSPB", test = "UrBAN", predictors = "audio (20)"),  uz$fob_total, { m <- ranger(reformulate(feat20, "fob_total"), mz, num.trees = 500, seed = 1); predict(m, uz)$predictions }),
  rowm(tibble(train = "UrBAN", test = "MSPB", predictors = "audio (20)"),  mz$fob_total, { m <- ranger(reformulate(feat20, "fob_total"), uz, num.trees = 500, seed = 1); predict(m, mz)$predictions }),
  rowm(tibble(train = "MSPB", test = "UrBAN", predictors = "box count"),   urb$fob_total, { m <- ranger(fob_total ~ n_boxes, mspb, num.trees = 500, seed = 1); predict(m, urb)$predictions }),
  rowm(tibble(train = "UrBAN", test = "MSPB", predictors = "box count"),   mspb$fob_total,{ m <- ranger(fob_total ~ n_boxes, urb,  num.trees = 500, seed = 1); predict(m, mspb)$predictions }))
print(as.data.frame(tr %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

## ---- save ----
write_csv(ucv, file.path(out_dir, "urban_cv_results.csv"))
write_csv(tr,  file.path(out_dir, "crossdataset_transfer.csv"))
saveRDS(urb, file.path(out_dir, "urban_matched_table.rds"))
cat(sprintf("\nICC(UrBAN)=%.2f vs ICC(MSPB)=0.21; wrote urban_cv_results.csv + crossdataset_transfer.csv\n", icc))

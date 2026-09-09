#!/usr/bin/env Rscript
# 12_revision_analyses.R — analyses added for revision 1 of COMPAG-D-26-05686
#
# Each block answers a specific reviewer request:
#  (A) predictor-set ablation — acoustic / environment / structure and their combinations,
#      with R2 and RMSE and hive-cluster bootstrap CIs                    [R2-3, R1-3.4, R1-5.2]
#  (B) incremental value of acoustics over box count + environment        [R2-3]
#  (C) lagged (previous-visit) box count sensitivity                      [R2-2]
#  (D) temporal / doubly-blocked CV + fold sample allocation              [R1-2.1, R1-3.2]
#  (E) domain-shift diagnostics: calibration offsets vs. genuine failure  [R1-3.1]
#  (F) split-conformal for the structural model                          [R1-3.3]
#  (G) saturation of the frames-of-bees scale                            [R1-5.1]
#
# Input : analysis/output/mspb_harmonized.rds (from 02), urban_matched_table.rds (from 07)
# Output: analysis/output/revision_*.csv + revision_analyses.txt (full log)
# Run from the repo root: Rscript analysis/12_revision_analyses.R

suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(ranger) })
set.seed(1)
out_dir <- here::here("analysis", "output")

bands   <- paste0("f_", 1:16)
core4   <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation")
audio20 <- c(core4, bands)
th      <- c("t_in_mean", "rh_in_mean")
feat_all <- c(audio20, th)

H   <- readRDS(file.path(out_dir, "mspb_harmonized.rds"))
raw <- H$mspb_tab

d <- raw %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes)) %>%
  dplyr::mutate(hive = factor(hive), yard = factor(yard),
                round = as.integer(sub("Evaluation ", "", eval)),
                upper = fob_total - brood1) %>%
  tidyr::drop_na(all_of(c(feat_all, "brood1")))

log_file <- file.path(out_dir, "revision_analyses.txt")
con <- file(log_file, open = "wt"); sink(con, split = TRUE)

cat("=========================================================================\n")
cat(" 12_revision_analyses.R —", R.version.string, "\n")
cat("=========================================================================\n")
cat(sprintf("MSPB modelling set: n = %d evaluations, %d colonies, apiaries: %s\n",
            nrow(d), nlevels(droplevels(d$hive)), paste(levels(d$yard), collapse = ", ")))

metrics <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat); y <- y[ok]; yhat <- yhat[ok]
  c(RMSE = sqrt(mean((y - yhat)^2)),
    MAE  = mean(abs(y - yhat)),
    R2   = 1 - sum((y - yhat)^2) / sum((y - mean(y))^2),
    corr = suppressWarnings(cor(y, yhat)))
}

## fold constructors -------------------------------------------------------
mkfolds <- function(df, scheme) switch(scheme,
  random = sample(rep_len(1:5, nrow(df))),
  hive   = as.integer(droplevels(df$hive)),
  apiary = as.integer(droplevels(df$yard)),
  round  = df$round)

# out-of-fold predictions for one predictor set under one scheme
cv_rf <- function(df, scheme, feats, target = "fob_total", ntree = 500) {
  fold <- mkfolds(df, scheme); pr <- rep(NA_real_, nrow(df))
  for (f in sort(unique(fold))) {
    tr <- df[fold != f, , drop = FALSE]; te <- df[fold == f, , drop = FALSE]
    if (nrow(te) == 0 || nrow(tr) < 10) next
    m <- ranger(reformulate(feats, target), data = tr, num.trees = ntree, seed = 1)
    pr[fold == f] <- predict(m, te)$predictions
  }
  pr
}

# strict double blocking: test cells share neither colony group nor evaluation round with training
cv_rf_blocked <- function(df, feats, target = "fob_total", n_hive_groups = 5, ntree = 500) {
  hv  <- levels(droplevels(df$hive))
  grp <- setNames(rep_len(seq_len(n_hive_groups), length(hv)), sample(hv))
  hg  <- grp[as.character(df$hive)]
  pr  <- rep(NA_real_, nrow(df))
  for (g in seq_len(n_hive_groups)) for (r in sort(unique(df$round))) {
    te_i <- which(hg == g & df$round == r)
    tr_i <- which(hg != g & df$round != r)
    if (!length(te_i) || length(tr_i) < 10) next
    m <- ranger(reformulate(feats, target), data = df[tr_i, , drop = FALSE],
                num.trees = ntree, seed = 1)
    pr[te_i] <- predict(m, df[te_i, , drop = FALSE])$predictions
  }
  pr
}

## =========================================================================
## (A) predictor-set ablation, with hive-cluster bootstrap CIs      [R2-3, R1-3.4]
## =========================================================================
cat("\n\n=== (A) Predictor-set ablation: what actually carries the signal? ===\n")
cat("Predictor sets are reported separately so that 'audio model' is never ambiguous.\n")

sets <- list(
  "acoustic only (20)"                  = audio20,
  "temperature + humidity only (2)"     = th,
  "acoustic + temp/humidity (22)"       = c(audio20, th),
  "box count only (1)"                  = "n_boxes",
  "box count + temp/humidity (3)"       = c("n_boxes", th),
  "box count + acoustic (21)"           = c("n_boxes", audio20),
  "box count + acoustic + temp/hum (23)" = c("n_boxes", audio20, th))

schemes_ab <- c("random", "hive", "apiary")
oof <- list()
for (s in schemes_ab) for (nm in names(sets)) {
  set.seed(1)
  oof[[paste(s, nm, sep = "|")]] <- cv_rf(d, s, sets[[nm]])
}

# one joint hive-cluster bootstrap so that differences between sets are paired
boot_all <- function(df, oof, keys, B = 2000, seed = 7) {
  set.seed(seed)
  hv  <- levels(droplevels(df$hive))
  idx <- split(seq_len(nrow(df)), droplevels(df$hive))
  reps <- vector("list", B)
  for (b in seq_len(B)) {
    hb   <- sample(hv, length(hv), replace = TRUE)
    rows <- unlist(idx[hb], use.names = FALSE)
    y <- df$fob_total[rows]
    reps[[b]] <- vapply(keys, function(k) {
      yh <- oof[[k]][rows]; ok <- is.finite(yh); yk <- y[ok]; yhk <- yh[ok]
      c(RMSE = sqrt(mean((yk - yhk)^2)),
        R2   = 1 - sum((yk - yhk)^2) / sum((yk - mean(yk))^2))
    }, numeric(2))
  }
  reps
}

keys_ab <- names(oof)
reps_ab <- boot_all(d, oof, keys_ab)

ci_of <- function(reps, key, stat) {
  v <- vapply(reps, function(m) m[stat, key], numeric(1))
  quantile(v, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
}

abl <- purrr::map_dfr(schemes_ab, function(s) purrr::map_dfr(names(sets), function(nm) {
  k <- paste(s, nm, sep = "|"); mm <- metrics(d$fob_total, oof[[k]])
  ci_r <- ci_of(reps_ab, k, "RMSE"); ci_q <- ci_of(reps_ab, k, "R2")
  tibble(scheme = s, predictors = nm,
         R2 = mm["R2"], R2_lo = ci_q[1], R2_hi = ci_q[2],
         RMSE = mm["RMSE"], RMSE_lo = ci_r[1], RMSE_hi = ci_r[2],
         MAE = mm["MAE"], corr = mm["corr"])
}))
print(as.data.frame(abl %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(abl, file.path(out_dir, "revision_ablation.csv"))

## =========================================================================
## (B) incremental contribution of acoustics                             [R2-3]
## =========================================================================
cat("\n\n=== (B) Incremental value of acoustics ON TOP OF box count (+ environment) ===\n")
cat("Paired hive-cluster bootstrap of the RMSE/R2 difference (same resampled colonies).\n")

contrasts <- tribble(
  ~label,                                            ~base,                            ~full,
  "add acoustics to box count",                      "box count only (1)",             "box count + acoustic (21)",
  "add acoustics to box count + temp/hum",           "box count + temp/humidity (3)",  "box count + acoustic + temp/hum (23)",
  "add temp/hum to box count",                       "box count only (1)",             "box count + temp/humidity (3)",
  "add box count to acoustic + temp/hum",            "acoustic + temp/humidity (22)",  "box count + acoustic + temp/hum (23)",
  "add acoustics to temp/hum only",                  "temperature + humidity only (2)", "acoustic + temp/humidity (22)")

inc <- purrr::map_dfr(c("hive", "apiary"), function(s) purrr::pmap_dfr(contrasts, function(label, base, full) {
  kb <- paste(s, base, sep = "|"); kf <- paste(s, full, sep = "|")
  dR <- vapply(reps_ab, function(m) m["RMSE", kf] - m["RMSE", kb], numeric(1))
  dQ <- vapply(reps_ab, function(m) m["R2",   kf] - m["R2",   kb], numeric(1))
  tibble(scheme = s, contrast = label,
         dRMSE = metrics(d$fob_total, oof[[kf]])["RMSE"] - metrics(d$fob_total, oof[[kb]])["RMSE"],
         dRMSE_lo = quantile(dR, .025, names = FALSE), dRMSE_hi = quantile(dR, .975, names = FALSE),
         dR2 = metrics(d$fob_total, oof[[kf]])["R2"] - metrics(d$fob_total, oof[[kb]])["R2"],
         dR2_lo = quantile(dQ, .025, names = FALSE), dR2_hi = quantile(dQ, .975, names = FALSE))
}))
print(as.data.frame(inc %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(inc, file.path(out_dir, "revision_incremental.csv"))
cat("\n(dRMSE > 0 and dR2 < 0 mean the added block HURTS out-of-sample accuracy.)\n")

## =========================================================================
## (C) lagged (previous-visit) box count                                  [R2-2]
## =========================================================================
cat("\n\n=== (C) Sensitivity: box count from the PREVIOUS visit (management-endogeneity check) ===\n")
cat("Reviewer 2: the same-visit configuration may encode the inspector's contemporaneous judgement.\n")
cat("Lags are built from the FULL inspection record (before the sensor filter) so no lag is lost.\n")

lagtab <- raw %>%
  dplyr::mutate(hive = as.character(hive),
                round = as.integer(sub("Evaluation ", "", eval))) %>%
  dplyr::arrange(hive, round) %>%
  dplyr::group_by(hive) %>%
  dplyr::mutate(n_boxes_prev = dplyr::lag(n_boxes),
                round_prev   = dplyr::lag(round),
                days_since_prev = as.numeric(date - dplyr::lag(date))) %>%
  dplyr::ungroup() %>%
  dplyr::select(hive, round, n_boxes_prev, round_prev, days_since_prev)

dl <- d %>%
  dplyr::mutate(hive = as.character(hive)) %>%
  dplyr::left_join(lagtab, by = c("hive", "round")) %>%
  dplyr::filter(!is.na(n_boxes_prev), round_prev == round - 1) %>%
  dplyr::mutate(hive = factor(hive), yard = factor(yard))

cat(sprintf("\nlagged subset: n = %d evaluations, %d colonies (rounds 2-6); median gap %.0f days\n",
            nrow(dl), nlevels(droplevels(dl$hive)), median(dl$days_since_prev, na.rm = TRUE)))
cat(sprintf("agreement between current and previous box count: %.1f%% identical; cor = %.2f\n",
            100 * mean(dl$n_boxes == dl$n_boxes_prev), cor(dl$n_boxes, dl$n_boxes_prev)))
cat("box count changes between consecutive visits:\n")
print(table(`change in boxes` = dl$n_boxes - dl$n_boxes_prev))

sets_lag <- list("box count (current visit)"                  = "n_boxes",
                 "box count (previous visit)"                 = "n_boxes_prev",
                 "box count (previous visit) + temp/humidity" = c("n_boxes_prev", th),
                 "box count (previous visit) + acoustic"      = c("n_boxes_prev", audio20),
                 "box count (previous visit) + acoustic + temp/hum" = c("n_boxes_prev", audio20, th),
                 "acoustic only (20)"                         = audio20,
                 "acoustic + temp/humidity (22)"              = c(audio20, th))
# Out-of-fold predictions are kept so one colony resample serves every predictor set, making the
# contrasts paired — the same construction as block (A), applied to the lagged subset.
oof_lag <- list()
for (s in c("hive", "apiary")) for (nm in names(sets_lag)) {
  set.seed(1)
  oof_lag[[paste(s, nm, sep = "|")]] <- cv_rf(dl, s, sets_lag[[nm]])
}
reps_lag <- boot_all(dl, oof_lag, names(oof_lag))

lagres <- purrr::map_dfr(c("hive", "apiary"), function(s) purrr::map_dfr(names(sets_lag), function(nm) {
  k <- paste(s, nm, sep = "|"); mm <- metrics(dl$fob_total, oof_lag[[k]])
  ci_r <- ci_of(reps_lag, k, "RMSE"); ci_q <- ci_of(reps_lag, k, "R2")
  tibble(scheme = s, predictors = nm, n = nrow(dl),
         R2 = mm["R2"], R2_lo = ci_q[1], R2_hi = ci_q[2],
         RMSE = mm["RMSE"], RMSE_lo = ci_r[1], RMSE_hi = ci_r[2])
}))
print(as.data.frame(lagres %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(lagres, file.path(out_dir, "revision_lagged_boxcount.csv"))

cat("\n-- paired contrast: previous-visit box count minus acoustic + temp/humidity --\n")
cat("(the text leans on 'matches or beats'; this attaches an interval to it)\n")
lag_pairs <- tribble(
  ~label,                                              ~base,                                        ~full,
  "previous-visit box count vs acoustic + temp/hum",   "acoustic + temp/humidity (22)",              "box count (previous visit)",
  "+ acoustic, to previous-visit box count",           "box count (previous visit)",                 "box count (previous visit) + acoustic",
  "+ acoustic, to previous-visit box count + temp/hum","box count (previous visit) + temp/humidity", "box count (previous visit) + acoustic + temp/hum")
lag_contrast <- purrr::map_dfr(c("hive", "apiary"), function(s) purrr::pmap_dfr(lag_pairs,
  function(label, base, full) {
    kb <- paste(s, base, sep = "|"); kf <- paste(s, full, sep = "|")
    dQ <- vapply(reps_lag, function(m) m["R2", kf] - m["R2", kb], numeric(1))
    tibble(scheme = s, contrast = label,
           dR2 = metrics(dl$fob_total, oof_lag[[kf]])["R2"] - metrics(dl$fob_total, oof_lag[[kb]])["R2"],
           dR2_lo = quantile(dQ, .025, names = FALSE), dR2_hi = quantile(dQ, .975, names = FALSE))
  }))
print(as.data.frame(lag_contrast %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(lag_contrast, file.path(out_dir, "revision_lagged_contrast.csv"))

cat("\n-- hive-power coefficient on the lagged subset (frames of bees per SD of hive power) --\n")
cat("No percentage shrink is quoted here: on this subset (rounds 2-6, early season dropped) the\n")
cat("marginal hive-power effect is already near zero, so a ratio would not be meaningful.\n")
hp_row <- function(fit, lbl) {
  s <- summary(fit)$coefficients["scale(hive_power)", ]
  tibble(model = lbl, est = s[["Estimate"]], se = s[["Std. Error"]], t = s[["t value"]])
}
att <- dplyr::bind_rows(
  hp_row(lmer(fob_total ~ scale(hive_power) + (1 | hive), dl, REML = FALSE,
              control = lmerControl(calc.derivs = FALSE)), "hive power alone"),
  hp_row(lmer(fob_total ~ scale(hive_power) + n_boxes_prev + (1 | hive), dl, REML = FALSE,
              control = lmerControl(calc.derivs = FALSE)), "+ previous-visit box count"),
  hp_row(lmer(fob_total ~ scale(hive_power) + n_boxes + (1 | hive), dl, REML = FALSE,
              control = lmerControl(calc.derivs = FALSE)), "+ current-visit box count"))
print(as.data.frame(att %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(att, file.path(out_dir, "revision_lagged_attenuation.csv"))

## =========================================================================
## (D) temporal / doubly-blocked CV + fold allocation             [R1-2.1, R1-3.2]
## =========================================================================
cat("\n\n=== (D1) Fold sample allocation for every cross-validation scheme ===\n")
alloc <- purrr::map_dfr(c("random", "hive", "apiary", "round"), function(s) {
  set.seed(1); fold <- mkfolds(d, s)
  sz <- as.integer(table(fold))
  hv <- vapply(split(as.character(d$hive), fold), function(x) length(unique(x)), integer(1))
  tibble(scheme = s, n_folds = length(sz),
         test_n_min = min(sz), test_n_median = median(sz), test_n_max = max(sz),
         train_n_median = nrow(d) - median(sz),
         test_colonies_median = median(hv),
         train_colonies_median = nlevels(droplevels(d$hive)) - median(hv))
})
print(as.data.frame(alloc), row.names = FALSE)
readr::write_csv(alloc, file.path(out_dir, "revision_fold_allocation.csv"))

cat("\n=== (D2) Temporal blocking: does residual temporal correlation inflate our scores? ===\n")
cat("leave-one-evaluation-round-out holds out a whole seasonal round;\n")
cat("hive+round blocked holds out cells sharing neither colony nor round with training.\n")
tsets <- list("acoustic + temp/humidity (22)" = c(audio20, th), "box count only (1)" = "n_boxes")
tcv <- dplyr::bind_rows(
  purrr::imap_dfr(tsets, function(f, nm) {
    set.seed(1); pr <- cv_rf(d, "round", f)
    dplyr::bind_cols(tibble(scheme = "leave-round-out", predictors = nm),
                     as_tibble_row(metrics(d$fob_total, pr))) }),
  purrr::imap_dfr(tsets, function(f, nm) {
    set.seed(1); pr <- cv_rf_blocked(d, f)
    dplyr::bind_cols(tibble(scheme = "hive+round blocked", predictors = nm),
                     as_tibble_row(metrics(d$fob_total, pr))) }))
print(as.data.frame(tcv %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(tcv, file.path(out_dir, "revision_temporal_cv.csv"))

cat("\n-- within-colony temporal autocorrelation of the seasonal evaluations --\n")
safe_lag_cor <- function(x) {
  y <- dplyr::lag(x); ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}
ac <- d %>% dplyr::arrange(hive, round) %>% dplyr::group_by(hive) %>%
  dplyr::summarise(r_fob = safe_lag_cor(fob_total),
                   r_hp  = safe_lag_cor(hive_power), .groups = "drop")
cat(sprintf("median lag-1 within-colony correlation: FoB %.2f | hive power %.2f (n = %d colonies)\n",
            median(ac$r_fob, na.rm = TRUE), median(ac$r_hp, na.rm = TRUE), sum(!is.na(ac$r_fob))))

## =========================================================================
## (E) domain-shift diagnostics: calibration vs. genuine failure          [R1-3.1]
## =========================================================================
cat("\n\n=== (E) Is failed transfer sensor-calibration drift, or a genuine limit? ===\n")

urb <- readRDS(file.path(out_dir, "urban_matched_table.rds"))
mspb_x <- raw %>%
  dplyr::filter(!is.na(fob_total), fob_total > 0, !is.na(hive_power)) %>%
  dplyr::select(fob_total, n_boxes, all_of(audio20)) %>% tidyr::drop_na(all_of(audio20))
urb_x <- urb %>% dplyr::select(fob_total, n_boxes, all_of(audio20))
z <- function(df, cols) dplyr::mutate(df, dplyr::across(all_of(cols), ~ as.numeric(scale(.x))))

cat("\n-- (E1) cross-dataset transfer WITH vs WITHOUT within-dataset z-scoring --\n")
cat("If failure were only an additive/multiplicative calibration offset, z-scoring would repair it.\n")
tr_res <- purrr::map_dfr(c(TRUE, FALSE), function(zs) {
  M <- if (zs) z(mspb_x, audio20) else mspb_x
  U <- if (zs) z(urb_x,  audio20) else urb_x
  set.seed(1)
  m1 <- ranger(reformulate(audio20, "fob_total"), M, num.trees = 500, seed = 1)
  m2 <- ranger(reformulate(audio20, "fob_total"), U, num.trees = 500, seed = 1)
  dplyr::bind_rows(
    dplyr::bind_cols(tibble(zscored = zs, train = "MSPB", test = "UrBAN", predictors = "acoustic (20)"),
                     as_tibble_row(metrics(U$fob_total, predict(m1, U)$predictions))),
    dplyr::bind_cols(tibble(zscored = zs, train = "UrBAN", test = "MSPB", predictors = "acoustic (20)"),
                     as_tibble_row(metrics(M$fob_total, predict(m2, M)$predictions))))
})
print(as.data.frame(tr_res %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(tr_res, file.path(out_dir, "revision_zscore_transfer.csv"))

cat("\n-- (E2) how large is the feature shift between sites, and between datasets? --\n")
cat("Standardised mean differences of the 20 shared acoustic features.\n")
smd <- function(a, b) {
  vapply(audio20, function(f) {
    x <- a[[f]]; y <- b[[f]]
    (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) /
      sqrt((var(x, na.rm = TRUE) + var(y, na.rm = TRUE)) / 2)
  }, numeric(1))
}
smd_site <- smd(dplyr::filter(d, yard == levels(d$yard)[1]), dplyr::filter(d, yard == levels(d$yard)[2]))
smd_data <- smd(mspb_x, urb_x)
shift <- tibble(comparison = c(paste(levels(d$yard), collapse = " vs "), "MSPB vs UrBAN"),
                hardware   = c("identical (same study, same season)", "different (different study, year, city)"),
                median_abs_SMD = c(median(abs(smd_site)), median(abs(smd_data))),
                max_abs_SMD    = c(max(abs(smd_site)),    max(abs(smd_data))),
                n_features_abs_SMD_gt_0.5 = c(sum(abs(smd_site) > 0.5), sum(abs(smd_data) > 0.5)))
print(as.data.frame(shift %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 2)))),
      row.names = FALSE)
readr::write_csv(shift, file.path(out_dir, "revision_domain_shift.csv"))
cat("\nKey point: the leave-apiary-out failure happens WITHIN MSPB, where the two apiaries share\n")
cat("the same sensor hardware, firmware, feature pipeline and season, so hardware/calibration\n")
cat("differences cannot be the explanation there.\n")

cat("\n-- (E3) cross-dataset transfer: acoustics vs. the box count, side by side --\n")
env_tr <- purrr::map_dfr(c("acoustic (20)", "box count"), function(p) {
  f <- if (p == "box count") "n_boxes" else audio20
  M <- if (p == "box count") mspb_x else z(mspb_x, audio20)
  U <- if (p == "box count") urb_x  else z(urb_x,  audio20)
  set.seed(1)
  m <- ranger(reformulate(f, "fob_total"), M, num.trees = 500, seed = 1)
  dplyr::bind_cols(tibble(train = "MSPB", test = "UrBAN", predictors = p),
                   as_tibble_row(metrics(U$fob_total, predict(m, U)$predictions)))
})
print(as.data.frame(env_tr %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)

## =========================================================================
## (F) split-conformal for the structural model                          [R1-3.3]
## =========================================================================
cat("\n\n=== (F) Split-conformal intervals: acoustic model vs. structural model ===\n")
cat("Nonconformity score: absolute residual |y - yhat|; nominal level 90% (alpha = 0.10);\n")
cat("calibration set always colony-disjoint from training; 80 random calibration splits.\n")

conf_sets <- list("acoustic + temp/humidity (22)" = c(audio20, th), "box count only (1)" = "n_boxes")
aps <- levels(d$yard)

# Finite-sample-valid split-conformal quantile: the ceil((n+1)(1-alpha))-th order statistic of the
# calibration scores, NOT the type-7 empirical quantile. Only this gives the "at least 1-alpha"
# marginal guarantee; with the small calibration sets used here the difference is material
# (type-7 targets ~0.87 expected coverage at n~28, the valid rule ~0.93).
conf_q <- function(scores, alpha = 0.10) {
  s <- sort(scores[is.finite(scores)]); n <- length(s)
  k <- ceiling((n + 1) * (1 - alpha))
  if (k > n) Inf else s[k]
}
one_rep <- function(seed, feats) {
  set.seed(seed); out <- list()
  for (ap in aps) {                       # in-distribution: colony-grouped 60/20/20 within apiary
    da <- d[d$yard == ap, ]; hv <- sample(unique(as.character(da$hive))); n <- length(hv)
    tr_h  <- hv[seq_len(floor(.6 * n))]
    cal_h <- hv[(floor(.6 * n) + 1):floor(.8 * n)]
    te_h  <- hv[(floor(.8 * n) + 1):n]
    if (!length(te_h) || !length(cal_h)) next
    fit <- ranger(reformulate(feats, "fob_total"), da[as.character(da$hive) %in% tr_h, ],
                  num.trees = 500, seed = 1)
    cal <- da[as.character(da$hive) %in% cal_h, ]; te <- da[as.character(da$hive) %in% te_h, ]
    q <- conf_q(abs(cal$fob_total - predict(fit, cal)$predictions))
    pr <- predict(fit, te)$predictions
    out[[length(out) + 1]] <- tibble(setting = paste0("in-dist: ", ap),
                                     cov = mean(abs(te$fob_total - pr) <= q), width = 2 * q)
  }
  for (ap in aps) {                       # transfer: calibrate on one apiary, test on the other
    tr <- d[d$yard == ap, ]; te <- d[d$yard != ap, ]
    hv <- sample(unique(as.character(tr$hive))); n <- length(hv)
    tr_h <- hv[seq_len(floor(.7 * n))]; cal_h <- hv[(floor(.7 * n) + 1):n]
    fit <- ranger(reformulate(feats, "fob_total"), tr[as.character(tr$hive) %in% tr_h, ],
                  num.trees = 500, seed = 1)
    cal <- tr[as.character(tr$hive) %in% cal_h, ]
    q <- conf_q(abs(cal$fob_total - predict(fit, cal)$predictions))
    pr <- predict(fit, te)$predictions
    out[[length(out) + 1]] <- tibble(setting = paste0("transfer: ", ap, "->", setdiff(aps, ap)),
                                     cov = mean(abs(te$fob_total - pr) <= q), width = 2 * q)
  }
  dplyr::bind_rows(out)
}
conf <- purrr::imap_dfr(conf_sets, function(f, nm)
  purrr::map_dfr(seq_len(80), ~ one_rep(.x, f)) %>% dplyr::mutate(model = nm))
conf_sum <- conf %>% dplyr::group_by(model, setting) %>%
  dplyr::summarise(reps = dplyr::n(), mean_cov = mean(cov), sd_cov = sd(cov),
                   mean_width = mean(width), .groups = "drop")
print(as.data.frame(conf_sum %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(conf_sum, file.path(out_dir, "revision_conformal_bymodel.csv"))

## =========================================================================
## (G) saturation of the frames-of-bees scale                            [R1-5.1]
## =========================================================================
cat("\n\n=== (G) Does the frames-of-bees scale saturate? ===\n")
cat("Occupancy = FoB / (10 x boxes): the fraction of installed frame capacity covered by bees.\n")
sat <- d %>% dplyr::group_by(n_boxes) %>%
  dplyr::summarise(n = dplyr::n(), mean_fob = mean(fob_total), mean_occ = mean(occupancy),
                   pct_occ_ge_0.9 = 100 * mean(occupancy >= 0.9),
                   pct_occ_ge_1.0 = 100 * mean(occupancy >= 1.0), .groups = "drop")
print(as.data.frame(sat %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 2)))),
      row.names = FALSE)
cat(sprintf("\nMSPB overall: %.0f%% of colony-evaluations at >= 90%% of nominal capacity; %.0f%% sit exactly\n",
            100 * mean(d$occupancy >= 0.9), 100 * mean(abs(d$occupancy - 1) < 1e-9)))
cat(sprintf("at the ten-frames-per-box ceiling (median occupancy %.2f); only %.1f%% exceed it.\n",
            median(d$occupancy), 100 * mean(d$occupancy > 1)))
cat(sprintf("bottom brood box full (>= 9.5 of 10 frames) in %.0f%% of multi-box evaluations.\n",
            100 * mean(d$brood1[d$n_boxes >= 2] >= 9.5)))

u21 <- readRDS(file.path(out_dir, "urban_harmonized.rds"))$analysis21 %>%
  dplyr::filter(!is.na(fob_total), fob_total > 0, !is.na(occupancy))
cat(sprintf("UrBAN (n = %d): median occupancy %.2f; %.0f%% at >= 90%% of nominal capacity; %.0f%% exactly at the ceiling.\n",
            nrow(u21), median(u21$occupancy), 100 * mean(u21$occupancy >= 0.9),
            100 * mean(abs(u21$occupancy - 1) < 1e-9)))
cat("Interpretation: within a fixed hive configuration the scale is close to its ceiling, so further\n")
cat("population growth is expressed mainly by adding a box rather than by a higher frame count.\n")
readr::write_csv(sat, file.path(out_dir, "revision_saturation.csv"))

## ---- figure: the ablation, for the manuscript ----
theme_set(theme_minimal(base_size = 11))
lvl <- rev(names(sets))
pf <- abl %>%
  dplyr::filter(scheme %in% c("hive", "apiary")) %>%
  dplyr::mutate(predictors = factor(predictors, levels = lvl),
                scheme = factor(scheme, levels = c("hive", "apiary"),
                                labels = c("leave-one-colony-out", "leave-one-apiary-out")),
                block = dplyr::case_when(grepl("^box count only", predictors) ~ "structure",
                                         grepl("^box count", predictors) ~ "structure + sensors",
                                         grepl("^acoustic only", predictors) ~ "acoustic",
                                         grepl("^temperature", predictors) ~ "environment",
                                         TRUE ~ "acoustic + environment"))
# The row labels already name every predictor set, so the colour legend is redundant; dropping it
# frees the vertical space needed to keep the axis type legible at the printed \linewidth (~6.7 in).
p <- ggplot(pf, aes(R2, predictors, colour = block)) +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey60") +
  geom_errorbar(aes(xmin = R2_lo, xmax = R2_hi), orientation = "y", width = .22, linewidth = .5) +
  geom_point(size = 2.1) +
  facet_wrap(~ scheme) +
  scale_colour_manual(values = c("structure" = "#0170B0", "structure + sensors" = "#6BAED6",
                                 "acoustic" = "#D75E01", "environment" = "#009E73",
                                 "acoustic + environment" = "#E69F00"), guide = "none") +
  labs(x = expression(paste("out-of-sample ", R^2, " (95% colony-cluster bootstrap CI)")), y = NULL) +
  theme_minimal(base_size = 10)
ggsave(file.path(out_dir, "fig_mspb_ablation.png"), p, width = 6.7, height = 2.9, dpi = 400)
cat("Wrote fig_mspb_ablation.png\n")

## ---- figure: conformal coverage, BOTH models (replaces the single-model version from 04) ----
pc <- conf %>%
  dplyr::mutate(kind = dplyr::if_else(grepl("^in-dist", setting), "in-distribution", "apiary transfer"),
                setting = sub("->", "→", setting, fixed = TRUE),
                model = factor(model, levels = names(conf_sets),
                               labels = c("acoustic + temperature/humidity", "box count only")))
pcf <- ggplot(pc, aes(setting, cov, fill = kind)) +
  geom_hline(yintercept = .90, linetype = 2, linewidth = .4) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.55, linewidth = 0.3) +
  geom_jitter(aes(colour = kind), width = 0.18, height = 0, size = 0.5, alpha = 0.45,
              show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 1.8, fill = "white", colour = "black") +
  coord_flip(ylim = c(0.35, 1.02)) +
  facet_wrap(~ model) +
  scale_fill_manual(values = c("in-distribution" = "#0170B0", "apiary transfer" = "#D75E01")) +
  scale_colour_manual(values = c("in-distribution" = "#0170B0", "apiary transfer" = "#D75E01")) +
  labs(x = NULL, y = "empirical coverage of nominal 90% intervals", fill = NULL) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig_mspb_conformal_transfer.png"), pcf,
       width = 6.4, height = 2.9, dpi = 400)
cat("Wrote fig_mspb_conformal_transfer.png (both models)\n")

## =========================================================================
## (H) robustness of the incremental-acoustics result to the mtry rule
## =========================================================================
cat("\n\n=== (H) Is the 'acoustics do not help' result an artefact of ranger's default mtry? ===\n")
cat("ranger's default mtry = floor(sqrt(p)) CHANGES with the predictor count, so the box count is a\n")
cat("split candidate in every split when alone but only ~1 in 5 once 20 acoustic features are present.\n")
cat("We therefore repeat the ablation under three mtry rules.\n")

cv_rf_mtry <- function(df, scheme, feats, rule, ntree = 500) {
  p <- length(feats)
  mt <- switch(rule, sqrt = max(1, floor(sqrt(p))), third = max(1, ceiling(p / 3)), all = p)
  fold <- mkfolds(df, scheme); pr <- rep(NA_real_, nrow(df))
  for (f in sort(unique(fold))) {
    tr <- df[fold != f, , drop = FALSE]; te <- df[fold == f, , drop = FALSE]
    if (nrow(te) == 0 || nrow(tr) < 10) next
    m <- ranger(reformulate(feats, "fob_total"), data = tr, num.trees = ntree,
                mtry = mt, seed = 1)
    pr[fold == f] <- predict(m, te)$predictions
  }
  pr
}

mtry_sets <- c("acoustic only (20)", "temperature + humidity only (2)", "acoustic + temp/humidity (22)",
               "box count only (1)", "box count + temp/humidity (3)", "box count + acoustic (21)",
               "box count + acoustic + temp/hum (23)")
oof_mtry <- list()
for (rule in c("sqrt", "third", "all")) for (s in c("hive", "apiary")) for (nm in mtry_sets) {
  set.seed(1)
  oof_mtry[[paste(rule, s, nm, sep = "|")]] <- cv_rf_mtry(d, s, sets[[nm]], rule)
}
reps_mtry <- boot_all(d, oof_mtry, names(oof_mtry))
mtry_res <- purrr::map_dfr(c("sqrt", "third", "all"), function(rule)
  purrr::map_dfr(c("hive", "apiary"), function(s)
    purrr::map_dfr(mtry_sets, function(nm) {
      k <- paste(rule, s, nm, sep = "|")
      dplyr::bind_cols(tibble(mtry_rule = rule, scheme = s, predictors = nm),
                       as_tibble_row(metrics(d$fob_total, oof_mtry[[k]]))) })))
print(as.data.frame(mtry_res %>% dplyr::select(mtry_rule, scheme, predictors, R2, RMSE) %>%
                      dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
readr::write_csv(mtry_res, file.path(out_dir, "revision_mtry_sensitivity.csv"))

get_r2 <- function(rule, s, nm) mtry_res$R2[mtry_res$mtry_rule == rule & mtry_res$scheme == s &
                                              mtry_res$predictors == nm]
cat("\n-- incremental effect of the acoustic block under each mtry rule (delta R2) --\n")
# paired intervals per rule, so the +0.024 under bagging can be judged against zero
mtry_ci <- function(rule, s, base, full) {
  kb <- paste(rule, s, base, sep = "|"); kf <- paste(rule, s, full, sep = "|")
  v <- vapply(reps_mtry, function(m) m["R2", kf] - m["R2", kb], numeric(1))
  c(quantile(v, .025, names = FALSE), quantile(v, .975, names = FALSE))
}
inc_mtry <- purrr::map_dfr(c("sqrt", "third", "all"), function(rule)
  purrr::map_dfr(c("hive", "apiary"), function(s) {
    ci_ab  <- mtry_ci(rule, s, "box count only (1)", "box count + acoustic (21)")
    ci_abe <- mtry_ci(rule, s, "box count + temp/humidity (3)", "box count + acoustic + temp/hum (23)")
    ci_b   <- mtry_ci(rule, s, "acoustic + temp/humidity (22)", "box count + acoustic + temp/hum (23)")
    tibble(mtry_rule = rule, scheme = s,
      d_add_A_to_B  = get_r2(rule, s, "box count + acoustic (21)") - get_r2(rule, s, "box count only (1)"),
      A_to_B_lo = ci_ab[1], A_to_B_hi = ci_ab[2],
      d_add_A_to_BE = get_r2(rule, s, "box count + acoustic + temp/hum (23)") -
                      get_r2(rule, s, "box count + temp/humidity (3)"),
      A_to_BE_lo = ci_abe[1], A_to_BE_hi = ci_abe[2],
      d_add_B_to_AE = get_r2(rule, s, "box count + acoustic + temp/hum (23)") -
                      get_r2(rule, s, "acoustic + temp/humidity (22)"),
      B_to_AE_lo = ci_b[1], B_to_AE_hi = ci_b[2],
      R2_box_only = get_r2(rule, s, "box count only (1)"),
      R2_best_acoustic = max(get_r2(rule, s, "acoustic only (20)"),
                             get_r2(rule, s, "acoustic + temp/humidity (22)")))}))
print(as.data.frame(inc_mtry %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))),
      row.names = FALSE)
readr::write_csv(inc_mtry, file.path(out_dir, "revision_mtry_increments.csv"))
cat(sprintf("\nRANGE of the acoustic increment on top of box count + environment: %.3f to %.3f\n",
            min(inc_mtry$d_add_A_to_BE), max(inc_mtry$d_add_A_to_BE)))
cat(sprintf("RANGE of box-count-only R2: %.3f to %.3f | best acoustic set: %.3f to %.3f\n",
            min(inc_mtry$R2_box_only), max(inc_mtry$R2_box_only),
            min(inc_mtry$R2_best_acoustic), max(inc_mtry$R2_best_acoustic)))
cat("-> the box-count advantage is invariant to the mtry rule; the SIGN of the acoustic increment is not.\n")

## =========================================================================
## (I) supporting diagnostics requested during internal review
## =========================================================================
cat("\n\n=== (I) Supporting diagnostics ===\n")

cat("\n-- (I1) scale of the cross-dataset transfer TEST set (UrBAN matched cells) --\n")
cat(sprintf("UrBAN matched test set: n = %d, colonies = %d, FoB range %.0f-%.0f, sd %.2f; box counts: %s\n",
            nrow(urb), length(unique(as.character(urb$hive))), min(urb$fob_total), max(urb$fob_total),
            sd(urb$fob_total), paste(names(table(urb$n_boxes)), table(urb$n_boxes),
                                     sep = "x", collapse = ", ")))
set.seed(11)
m_bx <- ranger(fob_total ~ n_boxes, mspb_x, num.trees = 500, seed = 1)
pr_bx <- predict(m_bx, urb)$predictions
hv_u <- unique(as.character(urb$hive)); idx_u <- split(seq_len(nrow(urb)), as.character(urb$hive))
bs <- replicate(2000, {
  hb <- sample(hv_u, length(hv_u), replace = TRUE)
  r <- unlist(idx_u[hb], use.names = FALSE); y <- urb$fob_total[r]
  1 - sum((y - pr_bx[r])^2) / sum((y - mean(y))^2)
})
cat(sprintf("box-count transfer R2 = %.3f, colony-cluster bootstrap 95%% CI [%.2f, %.2f]\n",
            metrics(urb$fob_total, pr_bx)["R2"],
            quantile(bs, .025, names = FALSE), quantile(bs, .975, names = FALSE)))

cat("\n-- (I2) per-round bias under leave-one-round-out --\n")
set.seed(1); pr_ae <- cv_rf(d, "round", c(audio20, th)); set.seed(1); pr_bx2 <- cv_rf(d, "round", "n_boxes")
rb <- d %>% dplyr::mutate(e_ae = pr_ae - fob_total, e_bx = pr_bx2 - fob_total) %>%
  dplyr::group_by(round) %>%
  dplyr::summarise(n = dplyr::n(), mean_fob = mean(fob_total), mean_hp = mean(hive_power),
                   bias_acoustic_env = mean(e_ae), rmse_acoustic_env = sqrt(mean(e_ae^2)),
                   bias_box = mean(e_bx), .groups = "drop")
print(as.data.frame(rb %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 2)))), row.names = FALSE)
int <- d$round %in% 2:5
cat(sprintf("interior rounds 2-5 only: acoustic+env R2 = %.3f | box count R2 = %.3f\n",
            metrics(d$fob_total[int], pr_ae[int])["R2"], metrics(d$fob_total[int], pr_bx2[int])["R2"]))

cat("\n-- (I3) within-stratum acoustic-population correlation, ALL box-count strata --\n")
ws <- d %>% dplyr::group_by(n_boxes) %>%
  dplyr::summarise(n = dplyr::n(), cor_hivepower_fob = cor(hive_power, fob_total), .groups = "drop")
print(as.data.frame(ws %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

cat("\n-- (I4) why the marginal hive-power effect vanishes on the rounds-2-6 subset --\n")
cat(sprintf("marginal cor(hive power, FoB): all rounds %.3f | round 1 only %.3f | rounds 2-6 %.3f\n",
            cor(d$hive_power, d$fob_total),
            cor(d$hive_power[d$round == 1], d$fob_total[d$round == 1]),
            cor(d$hive_power[d$round != 1], d$fob_total[d$round != 1])))
cat("round means (FoB | hive power):\n")
print(as.data.frame(d %>% dplyr::group_by(round) %>%
                      dplyr::summarise(mean_fob = round(mean(fob_total), 1),
                                       mean_hive_power = round(mean(hive_power), 2), .groups = "drop")),
      row.names = FALSE)

## =========================================================================
## (J) does "acoustics add nothing" depend on the acoustic block, the learner, or the window?
## =========================================================================
cat("\n\n=== (J) Robustness of the acoustic-increment result beyond the mtry rule ===\n")

cat("\n-- (J1) the 4-feature interpretable acoustic block instead of all 20 --\n")
cat("The 16 band powers are collinear with hive power, so 20 features may simply dilute the splits.\n")
sets_j1 <- list("acoustic4 only"        = core4,
                "acoustic4 + temp/hum"  = c(core4, th),
                "box count only"        = "n_boxes",
                "box count + temp/hum"  = c("n_boxes", th),
                "box count + acoustic4 + temp/hum" = c("n_boxes", core4, th))
j1 <- purrr::map_dfr(c("hive", "apiary"), function(s) purrr::imap_dfr(sets_j1, function(f, nm) {
  set.seed(1); pr <- cv_rf(d, s, f)
  dplyr::bind_cols(tibble(scheme = s, predictors = nm), as_tibble_row(metrics(d$fob_total, pr)))
}))
print(as.data.frame(j1 %>% dplyr::select(scheme, predictors, R2, RMSE) %>%
                      dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
g1 <- function(s, nm) j1$R2[j1$scheme == s & j1$predictors == nm]
for (s in c("hive", "apiary"))
  cat(sprintf("  %s: dR2 of adding acoustic4 to box count + temp/hum = %+.3f\n",
              s, g1(s, "box count + acoustic4 + temp/hum") - g1(s, "box count + temp/hum")))
readr::write_csv(j1, file.path(out_dir, "revision_acoustic4.csv"))

cat("\n-- (J2) a non-tree learner: the same contrasts under a linear mixed model --\n")
cat("Answers 'is it the information or the learner?' using the interpretable 4-feature block.\n")
cv_lmm <- function(df, scheme, feats) {
  fold <- mkfolds(df, scheme); pr <- rep(NA_real_, nrow(df))
  for (f in sort(unique(fold))) {
    tr <- df[fold != f, , drop = FALSE]; te <- df[fold == f, , drop = FALSE]
    if (nrow(te) == 0 || nrow(tr) < 10) next
    fit <- tryCatch(lmer(reformulate(c(feats, "(1|hive)"), "fob_total"), tr, REML = FALSE,
                         control = lmerControl(calc.derivs = FALSE)), error = function(e) NULL)
    if (!is.null(fit))
      pr[fold == f] <- tryCatch(as.numeric(predict(fit, te, allow.new.levels = TRUE)),
                                error = function(e) NA_real_)
  }
  pr
}
j2 <- purrr::map_dfr(c("hive", "apiary"), function(s) purrr::imap_dfr(sets_j1, function(f, nm) {
  pr <- cv_lmm(d, s, f)
  dplyr::bind_cols(tibble(scheme = s, predictors = nm), as_tibble_row(metrics(d$fob_total, pr)))
}))
print(as.data.frame(j2 %>% dplyr::select(scheme, predictors, R2, RMSE) %>%
                      dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
g2 <- function(s, nm) j2$R2[j2$scheme == s & j2$predictors == nm]
for (s in c("hive", "apiary"))
  cat(sprintf("  %s (LMM): dR2 of adding acoustic4 to box count + temp/hum = %+.3f | box count to acoustic4+temp/hum = %+.3f\n",
              s, g2(s, "box count + acoustic4 + temp/hum") - g2(s, "box count + temp/hum"),
              g2(s, "box count + acoustic4 + temp/hum") - g2(s, "acoustic4 + temp/hum")))
readr::write_csv(j2, file.path(out_dir, "revision_lmm_ablation.csv"))

cat("\n-- (J3) sensitivity to the sensor aggregation window (3, 7, 14 days before the inspection) --\n")
cat("The box count is read at the visit; the acoustic features are a window average, so the window is\n")
cat("a design choice that could disadvantage the acoustics. We vary it.\n")
sens <- H$sensor_daily %>% dplyr::mutate(hive = as.character(hive))
pop_dates <- raw %>% dplyr::transmute(hive = as.character(hive), eval, date, fob_total, n_boxes,
                                      brood1) %>%
  dplyr::filter(!is.na(fob_total), fob_total > 0, !is.na(n_boxes))
build_window <- function(k) {
  pop_dates %>%
    dplyr::inner_join(sens, by = "hive", relationship = "many-to-many",
                      suffix = c("", ".s")) %>%
    dplyr::filter(date.s <= date, date.s > date - k) %>%
    dplyr::group_by(hive, eval) %>%
    dplyr::summarise(dplyr::across(all_of(feat_all), ~ mean(.x, na.rm = TRUE)),
                     fob_total = dplyr::first(fob_total), n_boxes = dplyr::first(n_boxes),
                     .groups = "drop") %>%
    dplyr::left_join(raw %>% dplyr::transmute(hive = as.character(hive), eval, yard), by = c("hive","eval")) %>%
    dplyr::mutate(hive = factor(hive), yard = factor(yard)) %>%
    tidyr::drop_na(all_of(feat_all))
}
j3 <- purrr::map_dfr(c(3, 7, 14), function(k) {
  dk <- build_window(k)
  purrr::map_dfr(c("hive", "apiary"), function(s) purrr::imap_dfr(
    list("acoustic only (20)" = audio20, "acoustic + temp/humidity (22)" = c(audio20, th),
         "box count only (1)" = "n_boxes"),
    function(f, nm) {
      set.seed(1); pr <- cv_rf(dk, s, f)
      dplyr::bind_cols(tibble(window_days = k, n = nrow(dk), scheme = s, predictors = nm),
                       as_tibble_row(metrics(dk$fob_total, pr)))
    }))
})
print(as.data.frame(j3 %>% dplyr::select(window_days, n, scheme, predictors, R2, RMSE) %>%
                      dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
readr::write_csv(j3, file.path(out_dir, "revision_window_sensitivity.csv"))

cat("\n-- (J4) distribution of the target, for interpreting RMSE --\n")
cat(sprintf("MSPB modelling set: FoB range %.0f-%.0f, mean %.1f, sd %.2f, IQR %.0f-%.0f\n",
            min(d$fob_total), max(d$fob_total), mean(d$fob_total), sd(d$fob_total),
            quantile(d$fob_total, .25), quantile(d$fob_total, .75)))

## =========================================================================
## (K) bottom-vs-upper contrast: difference of the standardised coefficients
##     in scale(hive_power) ~ scale(brood1) + scale(upper) + (1|hive), with a
##     colony-cluster bootstrap CI on the DIFFERENCE (multi-box subset, as 05)
## =========================================================================
cat("\n\n=== (K) bottom-vs-upper coefficient difference, cluster-bootstrap CI ===\n")
dk <- H$mspb_tab %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes), !is.na(brood1)) %>%
  dplyr::mutate(hive = factor(hive), upper = fob_total - brood1) %>%
  dplyr::filter(n_boxes >= 2)
fit_bu <- function(dat) {
  m <- suppressWarnings(suppressMessages(
    lme4::lmer(scale(hive_power) ~ scale(brood1) + scale(upper) + (1 | hive), dat,
               REML = FALSE, control = lme4::lmerControl(calc.derivs = FALSE))))
  b <- lme4::fixef(m); unname(b["scale(brood1)"] - b["scale(upper)"])
}
est_bu <- fit_bu(dk)
set.seed(1)
hv <- levels(droplevels(dk$hive))
boot_bu <- replicate(2000, {
  smp <- sample(hv, length(hv), replace = TRUE)
  db <- dplyr::bind_rows(lapply(seq_along(smp), function(i)
    dk[dk$hive == smp[i], , drop = FALSE] %>% dplyr::mutate(hive = paste0(smp[i], "_", i))))
  tryCatch(fit_bu(db), error = function(e) NA_real_)
})
boot_bu <- boot_bu[is.finite(boot_bu)]
ci_bu <- quantile(boot_bu, c(.025, .975), names = FALSE)
cat(sprintf("beta_bottom - beta_upper = %.3f | 95%% colony-cluster bootstrap CI [%.3f, %.3f] (%d/2000 fits)\n",
            est_bu, ci_bu[1], ci_bu[2], length(boot_bu)))
readr::write_csv(tibble(estimate = est_bu, lo = ci_bu[1], hi = ci_bu[2], n_boot = length(boot_bu)),
                 file.path(out_dir, "revision_bottom_upper_contrast.csv"))

cat("\n\n=== done (12_revision_analyses) ===\n")
sink()
close(con)
cat("Wrote", log_file, "and revision_*.csv\n")

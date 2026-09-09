#!/usr/bin/env Rscript
# 13_hierarchical_conformal.R — analysis added for revision 2 of COMPAG-D-26-05686
#
# Reviewer 1 (R2, point 1): the split-conformal order-statistic quantile is finite-sample valid only
# for exchangeable scores. MSPB calibration residuals are clustered within colonies (up to six
# assessments per colony), so assessment-level residuals are not exchangeable with the score of a
# new observation from a new colony. This script repeats the conformal analysis of script 12 (F)
# with a colony-level (hierarchical) quantile that carries the finite-sample guarantee under
# *hierarchical* exchangeability: hierarchical conformal prediction (HCP) of Lee, Barber & Willett
# (2026, ACM J. Data Sci., doi:10.1145/3786352), eq. (6):
#
#   T = Q_{1-alpha}( sum_k sum_i  1/((K1+1) N_k) * delta_{s_ki}  +  1/(K1+1) * delta_{+inf} )
#
# where K1 is the number of calibration colonies and N_k the number of assessments of colony k.
# Each colony contributes total weight 1/(K1+1) regardless of how many assessments it has.
# With alpha = 0.10 the threshold is finite only if K1/(K1+1) >= 0.90, i.e. K1 >= 9 calibration
# colonies. The splits of script 12 (20 % of ~27 colonies in-distribution, 30 % under transfer)
# give K1 = 5 and 8, so the colony-level quantile is +inf there. We therefore (i) document that
# fact and (ii) re-run both the assessment-level rule and the colony-level rule on the SAME
# enlarged splits (in-distribution 40/40/20 colonies, transfer 60/40), so the two rules are
# compared on identical training and calibration data.
#
# Input : analysis/output/mspb_harmonized.rds (from 02)
# Output: analysis/output/hierarchical_conformal_*.csv + hierarchical_conformal.txt (full log)
# Run from the repo root, after 02: Rscript analysis/13_hierarchical_conformal.R
#
# What to compare with the paper (revision 2):
#   Design 1 ("manuscript splits")  -> the cov_assess / width_assess columns reproduce Table 9
#                                      (in-distribution 0.90/0.91 sensor, 0.91/0.92 box count;
#                                      transfer 0.70/0.99 sensor, 0.91/0.89 box count)
#   Design 2 ("enlarged calibration") -> Table 10, both rules on the same splits
# A copy of the log this script produced for the paper is kept in analysis/reference_output/.

suppressPackageStartupMessages({ library(tidyverse); library(ranger) })
set.seed(1)
in_dir  <- here::here("analysis", "output")
out_dir <- here::here("analysis", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bands   <- paste0("f_", 1:16)
core4   <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation")
audio20 <- c(core4, bands)
th      <- c("t_in_mean", "rh_in_mean")
feat_all <- c(audio20, th)

H   <- readRDS(file.path(in_dir, "mspb_harmonized.rds"))
raw <- H$mspb_tab
d <- raw %>%
  dplyr::filter(has_sensor, !is.na(fob_total), fob_total > 0, !is.na(n_boxes)) %>%
  dplyr::mutate(hive = factor(hive), yard = factor(yard),
                round = as.integer(sub("Evaluation ", "", eval)),
                upper = fob_total - brood1) %>%
  tidyr::drop_na(all_of(c(feat_all, "brood1")))

log_file <- file.path(out_dir, "hierarchical_conformal.txt")
con <- file(log_file, open = "wt"); sink(con, split = TRUE)
cat("=========================================================================\n")
cat(" 13_hierarchical_conformal.R —", R.version.string, "\n")
cat("=========================================================================\n")
cat(sprintf("MSPB modelling set: n = %d assessments, %d colonies, apiaries: %s\n",
            nrow(d), nlevels(droplevels(d$hive)), paste(levels(d$yard), collapse = ", ")))
cat(sprintf("Assessments per colony: median %d, range %d-%d\n",
            as.integer(median(table(d$hive))), min(table(d$hive)), max(table(d$hive))))

conf_sets <- list("sensor (22)" = c(audio20, th), "box count (1)" = "n_boxes")
aps   <- levels(d$yard)
alpha <- 0.10
B     <- 80

## ---- two quantile rules -------------------------------------------------------------------
# (a) assessment-level split conformal: the ceil((n+1)(1-alpha))-th smallest of n scores.
q_assess <- function(scores, alpha = 0.10) {
  s <- sort(scores[is.finite(scores)]); n <- length(s)
  k <- ceiling((n + 1) * (1 - alpha))
  if (k > n) Inf else s[k]
}
# (b) colony-level HCP (Lee, Barber & Willett 2026, eq. 6): weighted quantile with weight
#     1/((K1+1) N_k) per score of colony k and weight 1/(K1+1) on +inf.
q_hcp <- function(scores, group, alpha = 0.10) {
  ok <- is.finite(scores); scores <- scores[ok]; group <- as.character(group)[ok]
  K1 <- length(unique(group))
  Nk <- table(group)
  w  <- 1 / ((K1 + 1) * as.numeric(Nk[group]))
  o  <- order(scores); s <- scores[o]; w <- w[o]
  cw <- cumsum(w)                           # mass on finite scores; the remaining 1/(K1+1) sits at +inf
  j  <- which(cw >= 1 - alpha - 1e-12)[1]
  if (is.na(j)) Inf else s[j]
}
# sanity: HCP reduces to split conformal when every group has one observation
stopifnot(isTRUE(all.equal(q_hcp(1:20, 1:20), q_assess(1:20))))

## ---- one replicate: given colony splits, fit once, apply both rules ----------------------
one_rep <- function(seed, feats, frac_cal_id = 0.20, frac_tr_id = 0.60, frac_cal_tr = 0.30) {
  set.seed(seed); out <- list()
  eval_both <- function(fit, cal, te, setting) {
    sc <- abs(cal$fob_total - predict(fit, cal)$predictions)
    qa <- q_assess(sc, alpha); qh <- q_hcp(sc, cal$hive, alpha)
    pr <- predict(fit, te)$predictions; err <- abs(te$fob_total - pr)
    tibble(setting = setting,
           n_cal = nrow(cal), K1 = length(unique(as.character(cal$hive))), n_train = NA_integer_,
           cov_assess = mean(err <= qa), width_assess = 2 * qa,
           cov_hcp = mean(err <= qh), width_hcp = 2 * qh)
  }
  for (ap in aps) {                       # in-distribution: colony-grouped split within apiary
    da <- d[d$yard == ap, ]; hv <- sample(unique(as.character(da$hive))); n <- length(hv)
    n_tr  <- floor(frac_tr_id * n); n_cal <- floor((frac_tr_id + frac_cal_id) * n) - n_tr
    tr_h  <- hv[seq_len(n_tr)]
    cal_h <- hv[n_tr + seq_len(n_cal)]
    te_h  <- hv[(n_tr + n_cal + 1):n]
    if (!length(te_h) || !length(cal_h)) next
    trd <- da[as.character(da$hive) %in% tr_h, ]
    fit <- ranger(reformulate(feats, "fob_total"), trd, num.trees = 500, seed = 1)
    cal <- da[as.character(da$hive) %in% cal_h, ]; te <- da[as.character(da$hive) %in% te_h, ]
    r <- eval_both(fit, cal, te, paste0("in-dist: ", ap)); r$n_train <- nrow(trd)
    out[[length(out) + 1]] <- r
  }
  for (ap in aps) {                       # transfer: train + calibrate on one apiary, test on the other
    tr <- d[d$yard == ap, ]; te <- d[d$yard != ap, ]
    hv <- sample(unique(as.character(tr$hive))); n <- length(hv)
    n_tr <- floor((1 - frac_cal_tr) * n)
    tr_h <- hv[seq_len(n_tr)]; cal_h <- hv[(n_tr + 1):n]
    trd <- tr[as.character(tr$hive) %in% tr_h, ]
    fit <- ranger(reformulate(feats, "fob_total"), trd, num.trees = 500, seed = 1)
    cal <- tr[as.character(tr$hive) %in% cal_h, ]
    r <- eval_both(fit, cal, te, paste0("transfer: ", ap, "->", setdiff(aps, ap))); r$n_train <- nrow(trd)
    out[[length(out) + 1]] <- r
  }
  dplyr::bind_rows(out)
}

run_design <- function(label, frac_cal_id, frac_tr_id, frac_cal_tr) {
  res <- purrr::imap_dfr(conf_sets, function(f, nm)
    purrr::map_dfr(seq_len(B), function(s)
      one_rep(s, f, frac_cal_id = frac_cal_id, frac_tr_id = frac_tr_id, frac_cal_tr = frac_cal_tr)) %>%
      dplyr::mutate(model = nm))
  res$design <- label
  stopifnot(!anyNA(res$cov_assess), !anyNA(res$cov_hcp))
  res
}

## ---- design 1: the splits of the manuscript (in-dist 60/20/20, transfer 70/30) -----------
cat("\n\n=== Design 1: manuscript splits (in-distribution 60/20/20 colonies; transfer 70/30) ===\n")
cat("Colony-level (HCP) quantile at alpha = 0.10 needs K1 >= 9 calibration colonies to be finite.\n")
res1 <- run_design("manuscript splits", frac_cal_id = 0.20, frac_tr_id = 0.60, frac_cal_tr = 0.30)
sum1 <- res1 %>% dplyr::group_by(design, model, setting) %>%
  dplyr::summarise(reps = dplyr::n(), n_train = mean(n_train), n_cal = mean(n_cal), K1 = mean(K1),
                   sd_assess = sd(cov_assess), cov_assess = mean(cov_assess), width_assess = mean(width_assess),
                   sd_hcp = sd(cov_hcp), cov_hcp = mean(cov_hcp),
                   width_hcp = mean(width_hcp), share_hcp_finite = mean(is.finite(width_hcp)),
                   .groups = "drop")
print(as.data.frame(sum1 %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

## ---- design 2: enlarged calibration sets (in-dist 40/40/20, transfer 60/40) --------------
cat("\n\n=== Design 2: enlarged calibration (in-distribution 40/40/20 colonies; transfer 60/40) ===\n")
cat("Both rules applied to the SAME fitted model and the SAME calibration colonies.\n")
res2 <- run_design("enlarged calibration", frac_cal_id = 0.40, frac_tr_id = 0.40, frac_cal_tr = 0.40)
sum2 <- res2 %>% dplyr::group_by(design, model, setting) %>%
  dplyr::summarise(reps = dplyr::n(), n_train = mean(n_train), n_cal = mean(n_cal), K1 = mean(K1),
                   sd_assess = sd(cov_assess), cov_assess = mean(cov_assess), width_assess = mean(width_assess),
                   sd_hcp = sd(cov_hcp), cov_hcp = mean(cov_hcp),
                   width_hcp = mean(width_hcp), share_hcp_finite = mean(is.finite(width_hcp)),
                   .groups = "drop")
print(as.data.frame(sum2 %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

## ---- paired difference between the two rules, per split ----------------------------------
cat("\n\n=== Design 2: paired difference (HCP - assessment-level) in coverage and width, per split ===\n")
diff2 <- res2 %>% dplyr::mutate(d_cov = cov_hcp - cov_assess, d_width = width_hcp - width_assess) %>%
  dplyr::group_by(model, setting) %>%
  dplyr::summarise(d_cov_mean = mean(d_cov), d_cov_q05 = quantile(d_cov, .05), d_cov_q95 = quantile(d_cov, .95),
                   d_width_mean = mean(d_width), d_width_q05 = quantile(d_width, .05), d_width_q95 = quantile(d_width, .95),
                   .groups = "drop")
print(as.data.frame(diff2 %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))), row.names = FALSE)

readr::write_csv(dplyr::bind_rows(sum1, sum2), file.path(out_dir, "hierarchical_conformal_summary.csv"))
readr::write_csv(dplyr::bind_rows(res1, res2), file.path(out_dir, "hierarchical_conformal_splits.csv"))
readr::write_csv(diff2, file.path(out_dir, "hierarchical_conformal_paired.csv"))
cat("\nWrote hierarchical_conformal_summary.csv, hierarchical_conformal_splits.csv, hierarchical_conformal_paired.csv\n")
sink(); close(con)

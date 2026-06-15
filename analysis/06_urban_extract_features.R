#!/usr/bin/env Rscript
# 06_urban_extract_features.R
# Extract UrBAN hand-crafted acoustic features from the 2021 raw audio, reproducing the upstream
# `nectar_hand_crafted` definition (UrBAN feature_extraction.py) in R so the features match MSPB's
# stored 20-feature set. No Python/librosa dependency.
#
# Per recording: load WAV (16 kHz) -> resample to 15625 Hz (MSPB grid) -> three 1-min windows (start at
# minutes 1, 6, 11), each 29 non-overlapping 512-sample frames -> Hann window -> rFFT power (29 x 257).
# From bins 4:18 (hive band, R 5:18) and 4:257 (full, R 5:257) compute hive_power, audio_density,
# audio_density_ratio, density_variation and 16 band coefficients f_1..f_16 (bins 4:20, R 5:20);
# average over the 3 windows. Filenames: DD-MM-YYYY_HHhMM_HIVE-<id>.(wav|WAV).

suppressPackageStartupMessages({ library(tuneR); library(signal); library(tidyverse); library(lubridate) })
set.seed(1)

audio_dir <- here::here("urban", "UrBAN", "data", "audio", "beehives_2021")
insp_csv  <- here::here("urban", "UrBAN", "data", "annotations", "inspections_2021.csv")
out_dir   <- here::here("analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SR  <- 15625L; NF <- 512L; NFR <- 29L
hann <- 0.5 - 0.5 * cos(2 * pi * (0:(NF - 1)) / NF)            # periodic Hann (librosa fftbins=TRUE)
win_starts <- c(1, 6, 11) * 60L * SR
need <- max(win_starts) + NFR * NF
b_hive <- 5:18; b_full <- 5:257; b_coef <- 5:20
gcd <- function(a, b) { while (b) { t <- b; b <- a %% b; a <- t }; a }

pow_spec <- function(sig) {
  m <- matrix(sig, nrow = NF, ncol = NFR) * hann
  X <- mvfft(m)[1:(NF / 2 + 1), , drop = FALSE]
  t(Mod(X)^2)
}
nectar_features <- function(path) {
  w <- tryCatch(readWave(path), error = function(e) NULL); if (is.null(w)) return(NULL)
  x <- if (w@stereo) (w@left + w@right) / 2 else w@left
  x <- x / 2^(w@bit - 1)
  if (w@samp.rate != SR) { d <- gcd(SR, w@samp.rate); x <- signal::resample(x, SR %/% d, w@samp.rate %/% d) }
  if (length(x) < need) return(NULL)
  hp <- ad <- adr <- dv <- numeric(length(win_starts)); fb <- matrix(NA_real_, length(win_starts), 16)
  for (i in seq_along(win_starts)) {
    seg <- x[(win_starts[i] + 1):(win_starts[i] + NFR * NF)]
    P <- pow_spec(seg)
    band_sum <- rowSums(P[, b_hive, drop = FALSE]); full_sum <- rowSums(P[, b_full, drop = FALSE])
    hp[i] <- 10 * log10(mean(band_sum)); ad[i] <- 10 * log10(sum(band_sum) / 420)
    adr[i] <- sum(band_sum / full_sum) / 30; dv[i] <- 10 * log10(max(full_sum) / min(full_sum))
    fb[i, ] <- 10 * log10(colSums(P[, b_coef, drop = FALSE]) / 30)
  }
  c(hive_power = mean(hp), audio_density = mean(ad), audio_density_ratio = mean(adr),
    density_variation = mean(dv), setNames(colMeans(fb), paste0("f_", 1:16)))
}

## ---- file list, pre-filtered to within +/-3 days of an inspection ----
files <- list.files(audio_dir, pattern = "\\.(wav|WAV)$", full.names = TRUE)
if (length(files) == 0L)
  stop(sprintf(paste0("No .wav files found in %s\n",
                      "  The UrBAN raw audio (~8.3 GB, FRDR/Globus) is not included in the repo.\n",
                      "  Download it into that directory before running 06/07; the other scripts do not need it."),
               audio_dir), call. = FALSE)
insp_dates <- read_csv(insp_csv, show_col_types = FALSE) %>% dplyr::pull(Date) %>% ymd() %>% unique()
m <- str_match(basename(files), "^(\\d{2})-(\\d{2})-(\\d{4})_(\\d{2})h(\\d{2})_HIVE-(\\d+)")
fdate <- make_date(as.integer(m[, 4]), as.integer(m[, 3]), as.integer(m[, 2]))
keep  <- !is.na(fdate) & vapply(fdate, function(d) any(abs(as.integer(d - insp_dates)) <= 3), logical(1))
files <- files[keep]; mk <- m[keep, , drop = FALSE]; fdate <- fdate[keep]
cat(sprintf("recordings total in subset within +/-3d of inspections: %d of %d\n",
            length(files), length(keep)))

## ---- extract ----
rows <- vector("list", length(files)); n_ok <- 0L
for (k in seq_along(files)) {
  fe <- nectar_features(files[k]); if (is.null(fe) || any(!is.finite(fe))) next
  rows[[k]] <- dplyr::bind_cols(
    tibble(hive = mk[k, 7],
           datetime = make_datetime(as.integer(mk[k, 4]), as.integer(mk[k, 3]), as.integer(mk[k, 2]),
                                    as.integer(mk[k, 5]), as.integer(mk[k, 6])),
           date = fdate[k]),
    as_tibble_row(fe))
  n_ok <- n_ok + 1L
  if (k %% 50 == 0) cat(sprintf("  %d/%d (%d ok)\n", k, length(files), n_ok))
}
feat <- dplyr::bind_rows(rows)
if (nrow(feat) == 0L)
  stop("Found audio files but extracted 0 usable recordings (no finite features). ",
       "Check the audio files and the date/inspection matching window.", call. = FALSE)
cat(sprintf("\nextracted %d recordings | hives: %s | %s to %s\n", nrow(feat),
            paste(sort(unique(feat$hive)), collapse = ","), format(min(feat$date)), format(max(feat$date))))
cat("recordings per hive:\n"); print(table(feat$hive))
saveRDS(feat, file.path(out_dir, "urban_audio_features_2021.rds"))
write_csv(feat, file.path(out_dir, "urban_audio_features_2021.csv"))
cat("Wrote urban_audio_features_2021.{rds,csv}\n")

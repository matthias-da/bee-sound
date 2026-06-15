#!/usr/bin/env Rscript
# 02_build_mspb_table.R
# Harmonize MSPB (Zenodo 11398835) into the same hive x time structure as UrBAN:
#  - parse per-evaluation population (per-box frames -> total FoB, occupancy, composition)
#  - join the SHARED hand-crafted audio features + temp/humidity from the D1 sensor file
# This establishes the cross-dataset (UrBAN <-> MSPB) shared feature space for the C4 transfer arm.
#
# Hive-id linkage (verified): evaluation "Hive ID" == Nectar colony number ("02056");
# sensor tag_number == 200000 + Nectar number  =>  hive = sprintf("%05d", tag_number - 200000).

suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(lubridate); library(data.table)
})

mspb    <- here::here("mspb")
out_dir <- here::here("analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
hr <- function(t) cat("\n========== ", t, " ==========\n", sep = "")

# readxl returns the Dates column as a serial number OR as POSIXct depending on the sheet;
# handle both.
to_date <- function(x) {
  if (inherits(x, c("POSIXct", "POSIXt", "Date"))) as.Date(x)
  else as.Date(suppressWarnings(as.numeric(x)), origin = "1899-12-30")
}

box_names <- c("brood1", "brood2", "super1", "super2", "super3", "super4")  # bottom -> top

## ---- 1. Population evaluations 1-6 (per-box frames) ----
read_eval <- function(sheet) {
  d <- suppressMessages(read_excel(file.path(mspb, "D1_ant.xlsx"),
                                   sheet = sheet, skip = 2, col_names = FALSE))
  d <- d[, 1:min(10, ncol(d))]                       # 4 id cols + up to 6 box cols
  nb <- ncol(d) - 4L                                 # box columns actually present in this sheet
  names(d)[1:4] <- c("date_serial", "yard", "hive", "n_boxes")
  if (nb >= 1) names(d)[5:ncol(d)] <- box_names[seq_len(nb)]
  for (b in box_names) if (!b %in% names(d)) d[[b]] <- NA_real_   # pad missing box columns
  d %>%
    dplyr::filter(!is.na(hive)) %>%
    dplyr::mutate(
      eval    = sheet,
      date    = to_date(date_serial),
      hive    = sprintf("%05d", suppressWarnings(as.integer(hive))),   # normalise to Nectar id
      n_boxes = suppressWarnings(as.integer(n_boxes)),
      dplyr::across(all_of(box_names), ~ suppressWarnings(as.numeric(.)))
    ) %>%
    dplyr::select(eval, date, yard, hive, n_boxes, all_of(box_names))
}
pop <- map_dfr(paste("Evaluation", 1:6), read_eval) %>%
  dplyr::mutate(
    fob_total    = rowSums(dplyr::across(all_of(box_names)), na.rm = TRUE),
    n_boxes_obs  = rowSums(!is.na(dplyr::across(all_of(box_names)))),
    capacity     = 10L * n_boxes,                        # 10-frame boxes (assumption, as UrBAN)
    occupancy    = fob_total / capacity,
    share_bottom = brood1 / fob_total
  )

## ---- 2. D1 sensor -> daily per-hive aggregates (shared features + T/H) ----
feat_named <- c("hive_power", "audio_density", "audio_density_ratio", "density_variation")
sd1 <- fread(file.path(mspb, "D1_sensor_data.csv"),
             showProgress = FALSE)
hz_cols <- grep("^hz_", names(sd1), value = TRUE)
hz_cols <- hz_cols[order(as.numeric(sub("hz_", "", hz_cols)))]   # ascending frequency = f_1..f_16
sd1[, hive := sprintf("%05d", as.integer(tag_number) - 200000L)]
sd1[, date := as.Date(substr(published_at, 1, 10))]
sensor_daily <- sd1[, c(lapply(.SD, mean, na.rm = TRUE), .(n_readings = .N)),
                    by = .(hive, date),
                    .SDcols = c("temperature", "humidity", feat_named, hz_cols)]
setnames(sensor_daily, hz_cols, paste0("f_", seq_along(hz_cols)))   # harmonise band names to UrBAN
sensor_daily <- as_tibble(sensor_daily) %>%
  dplyr::rename(t_in_mean = temperature, rh_in_mean = humidity)

## ---- 3. Join sensor features over a 7-day window up to each evaluation ----
window_days <- 7
feat_all <- c("t_in_mean", "rh_in_mean", feat_named, paste0("f_", seq_along(hz_cols)))
sd2 <- sensor_daily %>% dplyr::rename(sdate = date)
sens_agg <- pop %>%
  dplyr::mutate(.eid = dplyr::row_number()) %>%
  dplyr::select(.eid, hive, edate = date) %>%
  dplyr::inner_join(sd2, by = "hive", relationship = "many-to-many") %>%
  dplyr::filter(sdate >= edate - (window_days - 1), sdate <= edate) %>%
  dplyr::group_by(.eid) %>%
  dplyr::summarise(dplyr::across(all_of(feat_all), ~ mean(.x, na.rm = TRUE)),
            n_sensor_days = dplyr::n(), .groups = "drop")
mspb_tab <- pop %>%
  dplyr::mutate(.eid = dplyr::row_number()) %>%
  dplyr::left_join(sens_agg, by = ".eid") %>%
  dplyr::mutate(has_sensor = !is.na(hive_power))

## ---- 4. EDA (mirror the UrBAN script) ----
hr("Coverage")
cat("hives (evaluated):", dplyr::n_distinct(pop$hive),
    "| apiaries:", paste(sort(unique(pop$yard)), collapse = ", "),
    "| evaluations:", dplyr::n_distinct(pop$eval), "\n")
cat("population rows:", nrow(pop),
    "| date range:", format(min(pop$date)), "to", format(max(pop$date)), "\n")
cat("sensor daily rows:", nrow(sensor_daily),
    "| sensor hives:", dplyr::n_distinct(sensor_daily$hive), "\n")

hr("Outcome: total FoB & occupancy")
cat("fob_total summary:\n"); print(summary(pop$fob_total))
cat("occupancy summary:\n"); print(summary(pop$occupancy))
cat("n_boxes (stated) distribution:\n"); print(table(pop$n_boxes))
cat("rows where stated n_boxes != observed non-NA boxes:",
    sum(pop$n_boxes != pop$n_boxes_obs, na.rm = TRUE), "\n")

hr("Bottom-box share by colony size (C2 mediator — compare to UrBAN)")
print(pop %>% dplyr::filter(!is.na(share_bottom)) %>%
        dplyr::group_by(n_boxes) %>%
        dplyr::summarise(n = dplyr::n(),
                  mean_share_bottom = mean(share_bottom, na.rm = TRUE),
                  mean_fob = mean(fob_total), .groups = "drop"))

hr("Between-hive variance teaser (ICC of occupancy) — compare to UrBAN 0.51")
icc_txt <- tryCatch({
  m  <- lme4::lmer(occupancy ~ 1 + (1 | hive), data = pop)
  vc <- as.data.frame(lme4::VarCorr(m))
  icc <- vc$vcov[vc$grp == "hive"] / sum(vc$vcov)
  sprintf("ICC(hive) for occupancy = %.3f", icc)
}, error = function(e) paste("ICC failed:", conditionMessage(e)))
cat(icc_txt, "\n")

hr("Shared feature space + join coverage")
shared <- c(feat_named, paste0("f_", seq_along(hz_cols)))
cat("shared hand-crafted features (n=", length(shared), "): ",
    paste(shared, collapse = ", "), "\n", sep = "")
cat("population rows with same-day sensor features:",
    sum(mspb_tab$has_sensor), "/", nrow(mspb_tab), "\n")

## ---- 5. Save + figures ----
write_csv(mspb_tab, file.path(out_dir, "mspb_analysis_table.csv"))
saveRDS(list(pop = pop, sensor_daily = sensor_daily, mspb_tab = mspb_tab),
        file.path(out_dir, "mspb_harmonized.rds"))

theme_set(theme_minimal(base_size = 11))
share_df <- pop %>% dplyr::filter(!is.na(share_bottom))
share_n  <- share_df %>% dplyr::count(n_boxes)
p1 <- ggplot(share_df, aes(factor(n_boxes), share_bottom)) +
  geom_boxplot(fill = "#0170B0", alpha = 0.30, colour = "grey30", outlier.size = 0.7) +
  geom_text(data = share_n, aes(factor(n_boxes), 1.06, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3, colour = "grey30") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
  labs(x = "number of boxes", y = "share in bottom brood box")
ggsave(file.path(out_dir, "fig_mspb_share_bottom.png"), p1, width = 8, height = 5, dpi = 130)

p2 <- ggplot(pop, aes(date, fob_total, group = hive, colour = yard)) +
  geom_line(alpha = .4) +
  labs(title = "MSPB: total frames of bees over the summer, by apiary",
       x = NULL, y = "frames of bees (total)", colour = "apiary")
ggsave(file.path(out_dir, "fig_mspb_fob.png"), p2, width = 9, height = 5, dpi = 130)

hr("DONE")
cat("Wrote:", file.path(out_dir, "mspb_analysis_table.csv"), "\n")
cat("Wrote:", file.path(out_dir, "mspb_harmonized.rds"), "\n")
cat("Figures: fig_mspb_share_bottom.png, fig_mspb_fob.png\n")

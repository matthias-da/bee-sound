#!/usr/bin/env Rscript
# 01_build_urban_table.R
# Step 0 of the bee colony-strength paper: harmonize the UrBAN inspection + in-hive
# sensor + external weather CSVs into a tidy hive x time analysis table, derive the
# outcome candidates (total FoB, occupancy, per-box composition), and run first EDA.
# Uses ONLY the CSVs already in the repo. Raw audio (FRDR) and MSPB (Zenodo) are not
# needed for this step and are not touched here.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

data_dir <- here::here("urban", "UrBAN", "data")
out_dir  <- here::here("analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hr <- function(t) cat("\n========== ", t, " ==========\n", sep = "")

## ---- 1. Inspections 2021 (wide: per-box FoB, colony size, queen) ----
insp21 <- read_csv(file.path(data_dir, "annotations/inspections_2021.csv"),
                   show_col_types = FALSE) %>%
  dplyr::transmute(
    hive    = factor(`Tag number`),
    date    = ymd(Date),
    n_boxes = as.integer(`Colony Size`),
    fob1    = as.numeric(`Fob 1st`),
    fob2    = as.numeric(`Fob 2nd`),
    fob3    = as.numeric(`Fob 3rd`),
    brood   = as.numeric(FoBrood),
    honey   = as.numeric(`Frames of Honey`),
    queen   = dplyr::na_if(toupper(str_trim(`Queen status`)), "")   # fixes "Qr" -> "QR", "" -> NA
  ) %>%
  dplyr::mutate(
    fob_total    = rowSums(cbind(fob1, fob2, fob3), na.rm = TRUE),
    capacity     = 10L * n_boxes,                 # 10-frame Langstroth boxes
    occupancy    = fob_total / capacity,          # in [0,1]; the scale-invariant target
    share_bottom = fob1 / fob_total,              # C2 mediator: mic sits in the bottom brood box
    share_mid    = fob2 / fob_total,
    share_top    = fob3 / fob_total,
    year         = 2021L
  )

## ---- 2. Inspections 2022 (long/event: total FoB only, sparse) ----
insp22 <- read_csv(file.path(data_dir, "annotations/inspections_2022.csv"),
                   show_col_types = FALSE) %>%
  dplyr::filter(str_to_lower(str_trim(Category)) == "frames of bees") %>%
  dplyr::transmute(
    hive      = factor(`Tag number`),
    datetime  = ymd_hms(Date, tz = "UTC"),
    date      = as_date(datetime),
    fob_total = suppressWarnings(as.numeric(`Action detail`)),
    queen     = str_to_lower(str_trim(`Queen status`)),
    is_alive  = as.integer(`Is alive`),
    year      = 2022L
  ) %>%
  dplyr::filter(!is.na(fob_total))

## ---- 3. In-hive sensor 2021 -> daily per-hive aggregates ----
sensor_daily <- read_csv(file.path(data_dir, "temperature_humidity/sensor_2021.csv"),
                         show_col_types = FALSE) %>%
  dplyr::transmute(hive = factor(`Tag number`),
            date = as_date(ymd_hms(Date, tz = "UTC")),
            temperature, humidity) %>%
  dplyr::group_by(hive, date) %>%
  dplyr::summarise(t_in_mean = mean(temperature, na.rm = TRUE),
            t_in_min  = min(temperature,  na.rm = TRUE),
            t_in_max  = max(temperature,  na.rm = TRUE),
            t_in_sd   = sd(temperature,   na.rm = TRUE),
            rh_in_mean = mean(humidity,   na.rm = TRUE),
            n_readings = dplyr::n(), .groups = "drop")

## ---- 4. External weather (global, LST-hourly) -> daily aggregates ----
weather_raw <- read_csv(file.path(data_dir, "weather_info/weather_2021_2022.csv"),
                        show_col_types = FALSE)
names(weather_raw) <- c("dt", "temp_ext", "rh_ext", "precip")
weather_daily <- weather_raw %>%
  dplyr::mutate(date = as_date(ymd_hms(dt))) %>%      # LST clock time -> local date
  dplyr::group_by(date) %>%
  dplyr::summarise(t_ext_mean = mean(temp_ext, na.rm = TRUE),
            t_ext_max  = max(temp_ext,  na.rm = TRUE),
            rh_ext_mean = mean(rh_ext,  na.rm = TRUE),
            precip_sum = sum(precip,    na.rm = TRUE), .groups = "drop")

## ---- 5. Build the 2021 labeled analysis table (same-day join) ----
analysis21 <- insp21 %>%
  dplyr::left_join(sensor_daily, by = c("hive", "date")) %>%
  dplyr::left_join(weather_daily, by = "date") %>%
  dplyr::mutate(has_sensor = !is.na(t_in_mean))

## ---- 6. EDA / sanity ----
hr("Hives: inspections vs in-hive sensor coverage")
hives_insp   <- levels(droplevels(insp21$hive))
hives_sensor <- sort(unique(as.character(sensor_daily$hive)))
cat("inspected hives (n=", length(hives_insp), "): ", paste(hives_insp, collapse = ", "), "\n", sep = "")
cat("sensor hives   (n=", length(hives_sensor), "): ", paste(hives_sensor, collapse = ", "), "\n", sep = "")
cat("inspected but NO in-hive sensor: ",
    paste(setdiff(hives_insp, hives_sensor), collapse = ", "), "\n", sep = "")

hr("Inspection counts & date ranges")
cat("2021 inspections:", nrow(insp21), "rows;",
    format(min(insp21$date)), "to", format(max(insp21$date)), "\n")
cat("2022 FoB inspections:", nrow(insp22), "rows;",
    format(min(insp22$date)), "to", format(max(insp22$date)), "\n")
cat("inspections per hive (2021):\n"); print(table(droplevels(insp21$hive)))
cat("distinct inspection dates (2022):", length(unique(insp22$date)), "->",
    paste(format(sort(unique(insp22$date))), collapse = ", "), "\n")

hr("Outcome: total FoB and occupancy")
cat("2021 fob_total summary:\n"); print(summary(insp21$fob_total))
cat("2022 fob_total summary:\n"); print(summary(insp22$fob_total))
cat("2021 occupancy summary (fob_total / (10*n_boxes)):\n"); print(summary(insp21$occupancy))
cat("n_boxes (colony size) distribution 2021:\n"); print(table(insp21$n_boxes))

hr("Structural vs count zeros (2021 per-box)")
# A box that exists (index <= n_boxes) but holds 0 frames is a COUNT zero;
# a non-existent box (index > n_boxes) is a STRUCTURAL zero (encoded NA).
count_zero_existing <- with(insp21,
  sum(fob1 == 0, na.rm = TRUE) +
  sum(fob2 == 0 & n_boxes >= 2, na.rm = TRUE) +
  sum(fob3 == 0 & n_boxes >= 3, na.rm = TRUE))
struct_zero <- with(insp21, sum(n_boxes < 2) + sum(n_boxes < 3))
cat("count zeros among existing boxes:", count_zero_existing, "\n")
cat("structural zeros (non-existent boxes):", struct_zero, "\n")

hr("Bottom-box share (C2 mediator candidate)")
cat("share_bottom overall summary:\n"); print(summary(insp21$share_bottom))
cat("mean share_bottom by colony size:\n")
print(insp21 %>% dplyr::group_by(n_boxes) %>%
        dplyr::summarise(n = dplyr::n(), mean_share_bottom = mean(share_bottom, na.rm = TRUE),
                  mean_fob = mean(fob_total), .groups = "drop"))

hr("Between-hive variance teaser (ICC of occupancy)")
icc_txt <- tryCatch({
  m <- lme4::lmer(occupancy ~ 1 + (1 | hive), data = insp21)
  vc <- as.data.frame(lme4::VarCorr(m))
  icc <- vc$vcov[vc$grp == "hive"] / sum(vc$vcov)
  sprintf("ICC(hive) for occupancy = %.3f  (share of variance that is between-hive)", icc)
}, error = function(e) paste("ICC model failed:", conditionMessage(e)))
cat(icc_txt, "\n")

hr("Join coverage (2021 analysis table)")
cat("rows:", nrow(analysis21),
    "| with in-hive sensor (same day):", sum(analysis21$has_sensor),
    "| with weather:", sum(!is.na(analysis21$t_ext_mean)), "\n")

## ---- 7. Save derived data + figures ----
write_csv(analysis21, file.path(out_dir, "urban_2021_analysis_table.csv"))
saveRDS(list(insp21 = insp21, insp22 = insp22, sensor_daily = sensor_daily,
             weather_daily = weather_daily, analysis21 = analysis21),
        file.path(out_dir, "urban_harmonized.rds"))

theme_set(theme_minimal(base_size = 11))

p1 <- ggplot(insp21, aes(date, fob_total, colour = hive, group = hive)) +
  geom_line(alpha = .6) + geom_point(size = 1.3) +
  labs(title = "UrBAN 2021: total frames of bees over the season",
       x = NULL, y = "frames of bees (total)", colour = "hive")
ggsave(file.path(out_dir, "fig_fob_2021.png"), p1, width = 9, height = 5, dpi = 130)

p2 <- insp21 %>% dplyr::filter(!is.na(share_bottom)) %>%
  ggplot(aes(date, share_bottom, colour = factor(n_boxes), group = hive)) +
  geom_line(alpha = .5) + geom_point(size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Bottom-brood-box share over the season (the mic's-eye view)",
       subtitle = "Fraction of total FoB in box 1, where the microphone sits (C2 mediator)",
       x = NULL, y = "share in bottom box", colour = "n boxes")
ggsave(file.path(out_dir, "fig_share_bottom_2021.png"), p2, width = 9, height = 5, dpi = 130)

p3 <- ggplot(insp21, aes(occupancy)) +
  geom_histogram(bins = 20, fill = "steelblue", colour = "white") +
  labs(title = "Occupancy = total FoB / (10 x n_boxes)", x = "occupancy", y = "count")
ggsave(file.path(out_dir, "fig_occupancy_2021.png"), p3, width = 7, height = 4, dpi = 130)

hr("DONE")
cat("Wrote:", file.path(out_dir, "urban_2021_analysis_table.csv"), "\n")
cat("Wrote:", file.path(out_dir, "urban_harmonized.rds"), "\n")
cat("Figures: fig_fob_2021.png, fig_share_bottom_2021.png, fig_occupancy_2021.png\n")

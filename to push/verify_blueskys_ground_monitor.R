library(tidyverse)
library(readxl)
library(lubridate)
library(terra)
library(sf)
library(exactextractr)
library(ggplot2)
library(ggtext)

raw <- read_excel("2021-2025_more_locations.xlsx", col_names = FALSE)
station_names <- as.character(raw[2, 3:ncol(raw)])

data_rows <- raw[5:nrow(raw), ]
colnames(data_rows) <- c("date_raw", "time_raw", station_names)

data_rows <- data_rows %>%
  filter(!is.na(date_raw), !date_raw %in% c("Minimum", "Maximum", "Median",
                                            "Mean", "Num", "Data[%]", "STD")) %>%
  mutate(date = as.Date(date_raw, "%m/%d/%Y")) %>%
  select(-date_raw, -time_raw)

monitors_daily <- data_rows %>%
  pivot_longer(-date, names_to = "station", values_to = "pm25_raw") %>%
  mutate(obs_pm25 = suppressWarnings(as.numeric(pm25_raw))) %>%
  select(-pm25_raw) %>%
  filter(!is.na(obs_pm25)) %>%
  mutate(
    month = month(date),
    year  = year(date),
    doy   = yday(date)
  )

coords_raw <- read_excel("air quality monitoring locations.xlsx", col_names = TRUE)
station_coords <- coords_raw %>%
  mutate(station = gsub("\u00a0", " ", station))

stations_sf <- station_coords %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Background removal
background_monthly <- monitors_daily %>%
  filter(month %in% c(11, 12, 1, 2, 3, 4)) %>%
  group_by(station, month) %>%
  summarise(bg_monthly_median = median(obs_pm25, na.rm = TRUE),
            n_obs = n(), .groups = "drop")

# Per-station baseline
station_bg_baseline <- monitors_daily %>%
  filter(month %in% c(11, 12, 1, 2, 3, 4)) %>%
  group_by(station) %>%
  summarise(bg_baseline = median(obs_pm25, na.rm = TRUE), .groups = "drop")

monitors_daily <- monitors_daily %>%
  left_join(background_monthly, by = c("station", "month")) %>%
  left_join(station_bg_baseline, by = "station") %>%
  mutate(
    bg_final            = dplyr::coalesce(bg_monthly_median, bg_baseline),
    obs_pm25_bg_removed = pmax(obs_pm25 - bg_final, 0)
  )

# check on a heavy smoke day
monitors_daily %>%
  filter(date == "2021-08-15") %>%
  select(station, obs_pm25, bg_final, obs_pm25_bg_removed) %>%
  arrange(desc(obs_pm25))

bluesky_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/University /Year 5/COGS 402/Data/All_years_blueskys_tif"
bluesky_files <- list.files(bluesky_dir, pattern = "\\.tif$", full.names = TRUE)

FILE_START_UTC <- 9  

get_noon_layer <- function(file_date, target_bc_hour = 12) {
  month_val   <- month(file_date)
  utc_offset  <- if (month_val %in% 3:11) -7L else -8L
  target_utc  <- (target_bc_hour - utc_offset) %% 24  
  layer_index <- target_utc - FILE_START_UTC + 1
  layer_index
}

extract_bluesky_noon <- function(filepath) {
  date_str  <- str_extract(basename(filepath), "\\d{8}")
  file_date <- as.Date(date_str, "%Y%m%d")
  
  if (is.na(file_date)) {
    warning("Could not parse date from: ", basename(filepath))
    return(NULL)
  }
  
  r           <- rast(filepath)
  layer_index <- get_noon_layer(file_date)
  
  if (layer_index < 1 || layer_index > nlyr(r)) {
    warning(basename(filepath), ": computed layer ", layer_index,
            " out of range (nlyr=", nlyr(r), "). Falling back to nearest layer.")
    layer_index <- pmin(pmax(layer_index, 1), nlyr(r))
  }
  
  r_noon <- r[[layer_index]]
  
  stations_reproj <- st_transform(stations_sf, crs(r_noon))
  vals <- terra::extract(r_noon, vect(stations_reproj), method = "bilinear")
  
  station_coords %>%
    mutate(
      station            = gsub("\u00a0", " ", station),
      date               = file_date,
      bluesky_pm25       = vals[, 2],
      bluesky_layer_used = layer_index
    ) %>%
    select(station, date, bluesky_pm25, bluesky_layer_used)
}

bluesky_extracted <- map_dfr(bluesky_files, extract_bluesky_noon)

# Quick check
bluesky_extracted %>%
  group_by(month = month(date)) %>%
  summarise(n = n(), mean_pm25 = round(mean(bluesky_pm25, na.rm = TRUE), 2)) %>%
  arrange(month)

paired <- monitors_daily %>%
  inner_join(bluesky_extracted, by = c("station", "date")) %>%
  filter(!is.na(obs_pm25), !is.na(bluesky_pm25)) %>%
  mutate(
    season = case_when(
      month %in% 5:10 ~ "Fire Season (May-Oct)",
      TRUE            ~ "Non-Fire Season (Nov-Apr)"
    )
  )

das_sf <- st_read("bc_das_population_for_raster.gpkg")

station_buffers <- st_buffer(
  st_transform(stations_sf, st_crs(das_sf)),
  dist = 10000  # 10 km radius
)

station_pop <- st_join(station_buffers, das_sf %>% select(pop_total_2021)) %>%
  st_drop_geometry() %>%
  group_by(station) %>%
  summarise(pop_total = sum(pop_total_2021, na.rm = TRUE), .groups = "drop")

station_weights <- station_pop %>%
  mutate(pop_weight = pop_total / sum(pop_total, na.rm = TRUE))

paired <- paired %>%
  left_join(station_weights, by = "station")

paired %>%
  select(station, pop_weight) %>%
  distinct() %>%
  arrange(desc(pop_weight)) %>%
  print(n = 10)

compute_stats <- function(obs, mod) {
  complete <- !is.na(obs) & !is.na(mod)
  obs <- obs[complete]
  mod <- mod[complete]
  
  n <- sum(complete)
  
  if (n < 3 || sd(obs) == 0 || sd(mod) == 0) {
    return(tibble(n = n, r = NA_real_, r2 = NA_real_,
                  MB = NA_real_, RMSE = NA_real_, NMB_pct = NA_real_))
  }
  
  mb   <- mean(mod - obs)
  rmse <- sqrt(mean((mod - obs)^2))
  nmb  <- mb / mean(obs) * 100
  r    <- cor(obs, mod)
  r2   <- r^2
  
  tibble(n = n, r = round(r, 3), r2 = round(r2, 3),
         MB = round(mb, 2), RMSE = round(rmse, 2), NMB_pct = round(nmb, 1))
}

compute_weighted_stats <- function(df, obs_col = "obs_pm25_bg_removed") {
  obs <- df[[obs_col]]
  mod <- df$bluesky_pm25
  w   <- df$pop_weight
  
  complete <- !is.na(obs) & !is.na(mod) & !is.na(w)
  obs <- obs[complete]; mod <- mod[complete]; w <- w[complete]
  w   <- w / sum(w)
  
  w_mean_obs <- weighted.mean(obs, w)
  w_mean_mod <- weighted.mean(mod, w)
  
  cov_wt    <- sum(w * (obs - w_mean_obs) * (mod - w_mean_mod))
  sd_obs_wt <- sqrt(sum(w * (obs - w_mean_obs)^2))
  sd_mod_wt <- sqrt(sum(w * (mod - w_mean_mod)^2))
  
  r_w    <- cov_wt / (sd_obs_wt * sd_mod_wt)
  r2_w   <- r_w^2
  mb_w   <- weighted.mean(mod - obs, w)
  rmse_w <- sqrt(weighted.mean((mod - obs)^2, w))
  nmb_w  <- mb_w / w_mean_obs * 100
  
  tibble(n = sum(complete),
         r_wtd = round(r_w, 3), r2_wtd = round(r2_w, 3),
         MB_wtd = round(mb_w, 2), RMSE_wtd = round(rmse_w, 2),
         NMB_pct_wtd = round(nmb_w, 1))
}

smoke_stations <- c("Castlegar Zinio Park", "Grand Forks City Hall",
                    "Kamloops Federal Building", "Kelowna KLO Road",
                    "Vernon Science Centre", "Cranbrook Muriel Baxter",
                    "Williams Lake Columneetza School", "Quesnel Johnston Avenue",
                    "Prince George Plaza 400", "Golden Helipad")

problem_stations <- c("Valemount", "Fort St John Key Learning Centre",
                      "Willow Creek Mine")

# --- Main validation summary ---
validation_summary <- bind_rows(
  paired %>%
    summarise(compute_stats(obs_pm25_bg_removed, bluesky_pm25)) %>%
    mutate(subset = "All stations, all seasons"),
  
  paired %>%
    filter(season == "Fire Season (May-Oct)") %>%
    summarise(compute_stats(obs_pm25_bg_removed, bluesky_pm25)) %>%
    mutate(subset = "All stations, fire season"),
  
  paired %>%
    filter(season == "Fire Season (May-Oct)", station %in% smoke_stations) %>%
    summarise(compute_stats(obs_pm25_bg_removed, bluesky_pm25)) %>%
    mutate(subset = "Interior stations, fire season"),
  
  paired %>%
    filter(season == "Fire Season (May-Oct)", !station %in% problem_stations) %>%
    summarise(compute_stats(obs_pm25_bg_removed, bluesky_pm25)) %>%
    mutate(subset = "Excl. industrial stations, fire season"),
  
  paired %>%
    filter(season == "Fire Season (May-Oct)") %>%
    summarise(compute_stats(obs_pm25, bluesky_pm25)) %>%
    mutate(subset = "All stations, fire season [RAW — no bg removal]")
  
) %>%
  select(subset, n, r, r2, MB, RMSE, NMB_pct)

print(validation_summary)

# --- Population-weighted summary (fire season) ---
paired %>%
  filter(season == "Fire Season (May-Oct)") %>%
  compute_weighted_stats() %>%
  print()

# --- By year (fire season) ---
paired %>%
  filter(season == "Fire Season (May-Oct)") %>%
  group_by(year) %>%
  reframe(compute_stats(obs_pm25_bg_removed, bluesky_pm25)) %>%
  print()

# --- By station ---
paired %>%
  group_by(station) %>%
  reframe(compute_stats(obs_pm25_bg_removed, bluesky_pm25)) %>%
  print()

# --- By season ---
paired %>%
  group_by(season) %>%
  reframe(compute_stats(obs_pm25_bg_removed, bluesky_pm25)) %>%
  print()

# --- Monthly R² breakdown (fire season) ---
monthly_r2 <- paired %>%
  filter(season == "Fire Season (May-Oct)") %>%
  mutate(month_label = month(date, label = TRUE, abbr = FALSE)) %>%
  group_by(month_label) %>%
  summarise(
    n              = n(),
    r2             = cor(obs_pm25_bg_removed, bluesky_pm25, use = "complete.obs")^2,
    r2_raw         = cor(obs_pm25, bluesky_pm25, use = "complete.obs")^2,
    mean_obs_bgrem = round(mean(obs_pm25_bg_removed, na.rm = TRUE), 2),
    mean_obs_raw   = round(mean(obs_pm25, na.rm = TRUE), 2),
    mean_bluesky   = round(mean(bluesky_pm25, na.rm = TRUE), 2),
    .groups        = "drop"
  ) %>%
  mutate(month_label = factor(month_label,
                              levels = c("May", "June", "July", "August",
                                         "September", "October")))

print(monthly_r2)

# --- Peak smoke days (top 10% of background-removed observed PM2.5) ---
thresh <- quantile(paired$obs_pm25_bg_removed, 0.90, na.rm = TRUE)
message("\n=== Peak smoke days (top 10% obs, threshold = ", round(thresh, 1), " µg/m³) ===")
paired %>%
  filter(obs_pm25_bg_removed >= thresh) %>%
  summarise(
    r2_peak   = round(cor(obs_pm25_bg_removed, bluesky_pm25, use = "complete.obs")^2, 3),
    n         = n(),
    threshold = round(thresh, 1)
  ) %>%
  print()

paired %>%
  group_by(season) %>%
  filter(station %in% smoke_stations) %>%
  summarise(n_both_nonNA = sum(!is.na(obs_pm25_bg_removed) & !is.na(bluesky_pm25)))

paired %>%
  filter(season == "Fire Season (May-Oct)") %>%
  summarise(
    mean_bgrem = mean(obs_pm25_bg_removed),
    sd_bgrem   = sd(obs_pm25_bg_removed),
    n_nonzero  = sum(obs_pm25_bg_removed > 0)
  )

# timeseries figure
ts_data <- paired %>%
  group_by(date) %>%
  summarise(
    obs   = mean(obs_pm25_bg_removed, na.rm = TRUE),
    model = mean(bluesky_pm25,        na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(obs, model),
               names_to  = "source",
               values_to = "pm25") %>%
  mutate(
    source = factor(source,
                    levels = c("obs", "model"),
                    labels = c("Observed (monitor, bg-removed)", "BlueSky modelled"))
  )

fire_bands <- tibble(
  xmin = as.Date(paste0(2021:2026, "-05-01")),
  xmax = as.Date(paste0(2021:2026, "-10-31"))
)

ggplot(ts_data, aes(x = date, y = pm25, colour = source)) +
  
  geom_rect(data = fire_bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "#F5A623", alpha = 0.08) +
  
  geom_line(linewidth = 0.5, alpha = 0.9) +
  
  facet_wrap(~ source, ncol = 1, scales = "free_y",
             strip.position = "left") +
  
  scale_colour_manual(
    values = c("Observed (monitor, bg-removed)" = "#2171B5",
               "BlueSky modelled"               = "#238B45"),
    guide = "none"
  ) +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               expand      = expansion(add = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  
  labs(
    title    = paste0("BlueSky Canada smoke forecast system vs. observed PM~2.5~ across British Columbia's Ministry of Environment<br>Air Quality Monitoring Stations (2021\u20132026)"),
    subtitle = "Station-mean daily values; shaded bands = May\u2013Oct fire season",
    x        = NULL,
    y        = "PM~2.5~ (µg/m³)"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    strip.placement      = "outside",
    strip.text           = element_markdown(size = 10, colour = "grey30"),
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    panel.spacing        = unit(1, "lines"),
    plot.title           = element_markdown(size = 12, face = "plain"),
    plot.subtitle        = element_markdown(size = 10, colour = "grey50"), 
    axis.title.y         = element_markdown(size = 9, colour = "grey30"),
    axis.title.x         = element_markdown(size = 9, colour = "grey30")
  )

ggsave("bluesky_vs_obs_timeseries_faceted.png", width = 10, height = 5, dpi = 300)


bind_rows(
  validation_summary,
  monthly_r2 %>%
    select(subset = month_label, n, r2) %>%
    mutate(subset = paste("Month:", subset))
) %>%
  write_csv("validation_stats_v3.csv")

# Counterfactual determination 
region_lookup <- tibble::tribble(
  ~station,                                        ~air_zone,
  
  # --- Northeast ---
  "Fort St John Key Learning Centre",              "Northeast",
  "Fort St John 85th Avenue_60",                   "Northeast",
  "Fort St John North Camp C_60",                  "Northeast",
  "Fort St John Old Fort_60",                        "Northeast",
  "Farmington Community Hall",                     "Northeast",
  "Peace Valley Attachie Flat Upper Terrace_60",   "Northeast",
  "Willow Creek Mine",                             "Northeast",
  
  # --- Central Interior ---
  "Prince George Plaza 400",                       "Central Interior",
  "Quesnel Johnston Avenue",                       "Central Interior",
  "Quesnel Kinchant St MAML",                      "Central Interior",
  "Williams Lake Columneetza School",              "Central Interior",
  "Valemount",                                     "Central Interior",
  "Burns Lake Fire Centre",                        "Central Interior",
  "Houston Firehall",                              "Central Interior",
  "Smithers Muheim Memorial",                      "Central Interior",
  
  # --- Coastal ---
  "Prince Rupert Fairview",                        "Coastal",
  "Terrace Skeena Middle School",                  "Coastal",
  "Kitimat Riverlodge",                            "Coastal",
  "Kitimat Whitesail",                             "Coastal",
  
  # --- Southern Interior ---
  "Castlegar Zinio Park",                          "Southern Interior",
  "Cranbrook Muriel Baxter",                       "Southern Interior",
  "Golden Helipad",                                "Southern Interior",
  "Grand Forks City Hall",                         "Southern Interior",
  "Kamloops Federal Building",                     "Southern Interior",
  "Kelowna KLO Road",                              "Southern Interior",
  "Vernon Science Centre",                         "Southern Interior",
  "Penticton Debeck Road",                         "Southern Interior",
  "Penticton Industrial Place",                    "Southern Interior",
  "Elkford Rocky Mountain Elementary School",      "Southern Interior",
  "Sparwood Centennial Square",                    "Southern Interior",
  "Skookumchuck Farstad Way",                      "Southern Interior",
  
  # --- Georgia Strait ---
  "Victoria Topaz",                                "Georgia Strait",
  "Colwood City Hall",                             "Georgia Strait",
  "Nanaimo Labieux Road",                          "Georgia Strait",
  "Courtenay Elementary School",                   "Georgia Strait",
  "Elk Falls Dogwood",                             "Georgia Strait",
  "Powell River James Thomson School",             "Georgia Strait",
  "Langdale Elementary",                           "Georgia Strait",
  "Squamish Elementary",                           "Georgia Strait",
  "Whistler Meadow Park",                          "Georgia Strait",
  "Crofton Substation",                            "Georgia Strait",
  "Duncan College Street",                         "Georgia Strait",
  "Duncan Deykin Avenue",                          "Georgia Strait",
  
  # --- Lower Fraser Valley ---
  "Abbotsford A Columbia Street",                  "Lower Fraser Valley",
  "Abbotsford Central",                            "Lower Fraser Valley",
  "Agassiz Municipal Hall",                        "Lower Fraser Valley",
  "Chilliwack Airport",                            "Lower Fraser Valley",
  "Mission School Works Yard",                     "Lower Fraser Valley",
  "Pitt Meadows Meadowlands School",               "Lower Fraser Valley",
  "Port Moody Rocky Point Park",                   "Lower Fraser Valley",
  "Burnaby Kensington Park",                       "Lower Fraser Valley",
  "Burnaby South",                                 "Lower Fraser Valley",
  "New Westminster Sapperton Park",                "Lower Fraser Valley",
  "North Delta",                                   "Lower Fraser Valley",
  "North Vancouver Mahon Park",                    "Lower Fraser Valley",
  "North Vancouver Second Narrows",                "Lower Fraser Valley",
  "Richmond South",                                "Lower Fraser Valley",
  "Surrey East",                                   "Lower Fraser Valley",
  "Tsawwassen",                                    "Lower Fraser Valley",
  "Vancouver Clark Drive",                         "Lower Fraser Valley",
  "Vancouver International Airport #2",            "Lower Fraser Valley",
  "Langley Central",                               "Lower Fraser Valley"
  
) %>%
  mutate(station = gsub("\u00a0", " ", station))

non_fire_months <- c(11, 12, 1, 2, 3, 4)

# check assignment works!
station_coords %>%
  mutate(station = gsub("\u00a0", " ", station)) %>%
  left_join(region_lookup, by = "station") %>%
  select(station, air_zone) %>%
  arrange(air_zone, station) %>%
  print(n = Inf)

non_fire_regional <- monitors_daily %>%
  filter(month %in% non_fire_months) %>%
  left_join(region_lookup, by = "station") %>%
  mutate(
    fire_year = case_when(
      month %in% c(11, 12) ~ year + 1L,
      TRUE                 ~ year
    )
  ) %>%
  filter(fire_year >= 2021)

# Air-zone-specific CF
cf_annual_regional <- non_fire_regional %>%
  group_by(air_zone, fire_year) %>%
  summarise(
    cf_mean    = round(mean(obs_pm25,   na.rm = TRUE), 3),
    cf_median  = round(median(obs_pm25, na.rm = TRUE), 3),
    cf_p25     = round(quantile(obs_pm25, 0.25, na.rm = TRUE), 3),
    cf_p75     = round(quantile(obs_pm25, 0.75, na.rm = TRUE), 3),
    n_obs      = n(),
    n_stations = n_distinct(station),
    .groups    = "drop"
  )

print(cf_annual_regional, n = Inf)
write_csv(cf_annual_regional, "cf_annual_regional.csv")

# Wide format for easy reading
cf_annual_regional %>%
  select(air_zone, fire_year, cf_median) %>%
  pivot_wider(names_from = fire_year, values_from = cf_median) %>%
  arrange(air_zone) %>%
  print()

# cap threshold analysis
fire_obs_all <- monitors_daily %>%
  filter(month %in% 5:10, !is.na(obs_pm25))

# Province-wide percentile profile
extreme_summary <- fire_obs_all %>%
  summarise(
    overall_max = round(max(obs_pm25, na.rm = TRUE), 1),
    p999        = round(quantile(obs_pm25, 0.999, na.rm = TRUE), 1),
    p99         = round(quantile(obs_pm25, 0.99,  na.rm = TRUE), 1),
    p975        = round(quantile(obs_pm25, 0.975, na.rm = TRUE), 1),
    p95         = round(quantile(obs_pm25, 0.95,  na.rm = TRUE), 1),
    n_obs_total = n()
  )

print(extreme_summary)

candidate_caps <- c(150, 200, 250, 300, 500)

obs_cap_exceedance <- map_dfr(candidate_caps, function(cap) {
  fire_obs_all %>%
    summarise(
      cap             = cap,
      n_days_exceed   = sum(obs_pm25 > cap, na.rm = TRUE),
      pct_days_exceed = round(sum(obs_pm25 > cap, na.rm = TRUE) / n() * 100, 3),
      n_stations_ever = n_distinct(station[obs_pm25 > cap]),
      max_obs         = round(max(obs_pm25, na.rm = TRUE), 1)
    )
})

print(obs_cap_exceedance)
write_csv(obs_cap_exceedance, "obs_cap_exceedance_monitor.csv")

fire_obs_all %>%
  left_join(region_lookup, by = "station") %>%
  arrange(desc(obs_pm25)) %>%
  select(date, station, air_zone, obs_pm25) %>%
  head(20) %>%
  print()


# tail histogram
p_cap_tail <- ggplot(filter(fire_obs_all, obs_pm25 > 100), aes(x = obs_pm25)) +
  geom_histogram(binwidth = 10, fill = "#D94F3D", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = c(150, 200, 250, 300),
             linetype = "dashed", colour = "black", linewidth = 0.6) +
  annotate("text", x = c(150, 200, 250, 300),
           y = Inf, label = c("150", "200", "250", "300"),
           vjust = 1.5, hjust = -0.15, size = 3.5, colour = "black") +
  labs(
    title    = "Tail of fire-season PM~2.5~ distribution (>100 \u00b5g/m\u00b3)",
    subtitle = "Used to assess appropriate cap threshold",
    x        = "Observed PM~2.5~ (\u00b5g/m\u00b3)",
    y        = "Count of station-days"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_markdown(size = 12),
        plot.subtitle = element_markdown(size = 8))

print(p_cap_tail)
ggsave("obs_pm25_tail_distribution.png", p_cap_tail,
       width = 9, height = 5, dpi = 300, bg = "white")


# Province-wide CF comparision
cf_province <- monitors_daily %>%
  filter(month %in% non_fire_months) %>%
  mutate(fire_year = case_when(
    month %in% c(11, 12) ~ year + 1L,
    TRUE                 ~ year
  )) %>%
  filter(fire_year >= 2021) %>%
  group_by(fire_year) %>%
  summarise(cf_median = round(median(obs_pm25, na.rm = TRUE), 3), .groups = "drop")

p_cf_region <- ggplot(cf_annual_regional,
                      aes(x = fire_year, y = cf_median,
                          colour = air_zone, group = air_zone)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_line(data = cf_province,
            aes(x = fire_year, y = cf_median),
            inherit.aes = FALSE,
            linetype = "dashed", colour = "black", linewidth = 0.8) +
  scale_colour_brewer(palette = "Set1", name = "Air Zone") +
  scale_x_continuous(breaks = min(cf_annual_regional$fire_year):
                       max(cf_annual_regional$fire_year)) +
  labs(
    title    = "Regional vs. province-wide background PM~2.5~ counterfactual (non-fire-season median)",
    subtitle = "Dashed line = province-wide background; coloured lines = regional background",
    x        = "Fire year",
    y        = "Background PM~2.5~ (\u00b5g/m\u00b3)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_markdown(size = 12),
        plot.subtitle = element_markdown(size = 8))

print(p_cf_region)
ggsave("regional_cf_comparison.png", p_cf_region,
       width = 9, height = 5, dpi = 300, bg = "white")

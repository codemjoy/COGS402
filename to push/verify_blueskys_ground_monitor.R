library(tidyverse)
library(readxl)
library(lubridate)
library(terra)
library(sf)

raw <- read_excel("2021-2026 24 hr mean air quality archieve.xlsx", col_names = FALSE)

station_names <- as.character(raw[2, 3:ncol(raw)])

data_rows <- raw[5:nrow(raw), ]
colnames(data_rows) <- c("date_raw", "time_raw", station_names)

# Strip summary footer rows 
data_rows <- data_rows %>%
  filter(!is.na(date_raw), !date_raw %in% c("Minimum", "Maximum", "Median",
                                            "Mean", "Num", "Data[%]", "STD"))

data_rows <- data_rows %>%
  mutate(date = as.Date(date_raw, "%m/%d/%Y")) %>%
  select(-date_raw, -time_raw)

# Pivot to long format
monitors_daily <- data_rows %>%
  pivot_longer(-date, names_to = "station", values_to = "pm25_raw") %>%
  mutate(obs_pm25 = suppressWarnings(as.numeric(pm25_raw))) %>%
  select(-pm25_raw) %>%
  filter(!is.na(obs_pm25))

# STEP 2: Load station coordinates
coords_raw <- read_excel("air quality monitoring locations.xlsx", col_names = TRUE)

# The first column is the row label ("Lat" / "Long"), remaining cols are stations
station_coords <- coords_raw %>%
  rename(coord_type = 1) %>%
  pivot_longer(-coord_type, names_to = "station", values_to = "value") %>%
  pivot_wider(names_from = coord_type, values_from = value) %>%
  rename(lat = Lat, lon = Long)

# Convert to spatial object (WGS84)
stations_sf <- station_coords %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# STEP 3: Extract BlueSky values at monitor locations
bluesky_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/University /Year 5/COGS 402/Data/All_years_blueskys_tif"
bluesky_files <- list.files(bluesky_dir, pattern = "\\.tif$", full.names = TRUE)

extract_bluesky_day <- function(filepath) {
  
  r <- rast(filepath)
  
  # Average all hourly disepersion layer to single mean
  r_daily <- mean(r, na.rm = TRUE) * 1e9
  
  # Check CRS and reproject stations if needed
  stations_reproj <- st_transform(stations_sf, crs(r_daily))
  
  # Extract PM2.5 value at each station location
  vals <- terra::extract(r_daily, vect(stations_reproj), method = "bilinear")

  # Parse date from filename 
  date_str <- str_extract(basename(filepath), "\\d{8}")
  parsed_date <- as.Date(date_str, "%Y%m%d")
  
  if (is.na(parsed_date)) {
    warning("Could not parse date from: ", basename(filepath))
    return(NULL)
  }
  
  station_coords %>%
    mutate(
      date         = parsed_date,
      bluesky_pm25 = vals[[2]] * 1e9   
    ) %>%
    select(station, date, bluesky_pm25)
}

bluesky_extracted <- map_dfr(bluesky_files, extract_bluesky_day)

# STEP 4: Build paired dataset
paired <- monitors_daily %>%
  inner_join(bluesky_extracted, by = c("station", "date")) %>%
  filter(!is.na(obs_pm25), !is.na(bluesky_pm25)) %>%
  mutate(
    month  = month(date),
    season = case_when(
      month %in% 5:10 ~ "Fire Season (May-Oct)",
      TRUE            ~ "Non-Fire Season (Nov-Apr)"
    ),
    year = year(date)
  )

# verify
paired %>%
  filter(date == "2021-08-02") %>%
  select(station, obs_pm25, bluesky_pm25) %>%
  arrange(desc(obs_pm25))

# STEP 5: Compute validation statistics
compute_stats <- function(df) {
  obs <- df$obs_pm25
  mod <- df$bluesky_pm25
  n   <- sum(!is.na(obs) & !is.na(mod))
  
  mb   <- mean(mod - obs, na.rm = TRUE)                        # Mean Bias
  rmse <- sqrt(mean((mod - obs)^2, na.rm = TRUE))              # RMSE
  nmb  <- mb / mean(obs, na.rm = TRUE) * 100                   # Normalized Mean Bias (%)
  r    <- cor(obs, mod, use = "complete.obs")                  # Pearson correlation
  r2   <- r^2                                                   # R-squared
  
  tibble(n = n, r = round(r, 3), r2 = round(r2, 3),
         MB = round(mb, 2), RMSE = round(rmse, 2), NMB_pct = round(nmb, 1))
}

# Overall performance
stats_overall <- paired %>%
  summarise(compute_stats(cur_data())) %>%
  mutate(group = "Overall")

# By season
stats_season <- paired %>%
  group_by(season) %>%
  group_modify(~ compute_stats(.x)) %>%
  rename(group = season) %>%
  ungroup()

# Fire Season specifically
cor(paired$obs_pm25, paired$bluesky_pm25, use = "complete.obs")^2
paired %>%
  filter(season == "Fire Season (May-Oct)") %>%
  summarise(r2 = cor(obs_pm25, bluesky_pm25, use = "complete.obs")^2)

# By station
stats_station <- paired %>%
  group_by(station) %>%
  group_modify(~ compute_stats(.x)) %>%
  rename(group = station) %>%
  ungroup()

# For heavy fire stations
smoke_stations <- c("Castlegar Zinio Park", "Grand Forks City Hall", 
                    "Kamloops Federal Building", "Kelowna KLO Road",
                    "Vernon Science Centre", "Cranbrook Muriel Baxter",
                    "Williams Lake Columneetza School", "Quesnel Johnston Avenue",
                    "Prince George Plaza 400", "Golden Helipad")

paired %>%
  filter(season == "Fire Season (May-Oct)",
         station %in% smoke_stations) %>%
  summarise(r2 = cor(obs_pm25, bluesky_pm25, use = "complete.obs")^2)

# For heavy smoke days
paired %>%
  filter(season == "Fire Season (May-Oct)",
         obs_pm25 > 15) %>%
  summarise(r2 = cor(obs_pm25, bluesky_pm25, use = "complete.obs")^2)

# By year
stats_year <- paired %>%
  group_by(year) %>%
  group_modify(~ compute_stats(.x)) %>%
  mutate(group = as.character(year)) %>%
  select(-year) %>%
  ungroup()

paired %>%
  filter(season == "Fire Season (May-Oct)") %>%
  group_by(year) %>%
  summarise(r2 = cor(obs_pm25, bluesky_pm25, use = "complete.obs")^2,
            n  = n(),
            mean_obs = mean(obs_pm25))

problem_stations <- c("Valemount", "Fort St John Key Learning Centre", 
                      "Willow Creek Mine")

# summary table
validation_summary <- bind_rows(
  paired %>%
    summarise(compute_stats(cur_data())) %>%
    mutate(subset = "All stations, all seasons"),
  
  paired %>%
    filter(season == "Fire Season (May-Oct)") %>%
    summarise(compute_stats(cur_data())) %>%
    mutate(subset = "All stations, fire season"),
  
  paired %>%
    filter(season == "Fire Season (May-Oct)",
           station %in% smoke_stations) %>%
    summarise(compute_stats(cur_data())) %>%
    mutate(subset = "Interior stations, fire season"),
  
  paired %>%
    filter(season == "Fire Season (May-Oct)",
           !station %in% problem_stations) %>%
    summarise(compute_stats(cur_data())) %>%
    mutate(subset = "Excl. industrial stations, fire season")
) %>%
  select(subset, n, r, r2, MB, RMSE, NMB_pct)

print(validation_summary)

# Save stats to CSV
bind_rows(stats_overall, stats_season, stats_station, stats_year) %>%
  write_csv("validation_stats.csv")

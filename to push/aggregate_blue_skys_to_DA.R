library(sf)
library(dplyr)
library(terra)
library(exactextractr)
library(lubridate)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyr)

year        <- 2025
bluesky_dir <- "/Users/minalander/Library/Mobile Documents/com~apple~CloudDocs/Documents/University /Year 5/COGS 402/Data/Blue skys tif /2025"

# BlueSky Canada BSC00CA12-01: 51-band hourly output, initialized 00:00 UTC
# UTC-7 = PDT (summer); band 1 = 00:00 UTC = 17:00 PDT prior day
# Band 8  = 07:00 UTC = 00:00 PDT (start of filename calendar day)
# Band 31 = 30:00 UTC = 23:00 PDT (end of filename calendar day)
# Bands 32-51 = 00:00-19:00 PDT following day (forecast overhang, excluded)
#
# PRIMARY:     bands 8-31  = full 24-hr PDT mean (aligns with Chen et al. CRF)
# SENSITIVITY: bands 20-26 = 12:00-18:00 PDT afternoon mean
#              (avoids nocturnal BlueSky accumulation artifact;
#               produces higher estimates — reported as sensitivity)

BAND_PRIMARY     <- 8:31  
BAND_AFTERNOON   <- 20:26  

census_das             <- st_read("bc_das_population_for_raster.gpkg")
census_das_transformed <- st_transform(census_das, 4326)

bluesky_files <- list.files(bluesky_dir, pattern = "\\.tif$", full.names = TRUE)
cat("Found", length(bluesky_files), "files for year", year, "\n")

extract_date <- function(filename) {
  date_str <- str_extract(basename(filename), "(?<=_)\\d{8}(?=_)")
  ymd(date_str)
}

file_df <- data.frame(filepath = bluesky_files) %>%
  mutate(date = extract_date(filepath)) %>%
  filter(!is.na(date)) %>%
  arrange(date)

# extraction loop
daily_results <- list()

for(i in seq_along(file_df$filepath)) {
  
  current_date <- file_df$date[i]
  current_file <- file_df$filepath[i]
  
  r_full <- rast(current_file)
  
  
  if(nlyr(r_full) != 51) {
    warning(sprintf(
      "Unexpected band count %d in %s — band mapping may be wrong, check before proceeding",
      nlyr(r_full), basename(current_file)
    ))
  }
  
  r_24hr       <- mean(r_full[[BAND_PRIMARY]],   na.rm = TRUE)
  r_afternoon  <- mean(r_full[[BAND_AFTERNOON]], na.rm = TRUE)
  
  # --- Extract to DAs ---
  da_24hr      <- exact_extract(r_24hr,      census_das_transformed,
                                fun = 'mean', progress = FALSE)
  da_afternoon <- exact_extract(r_afternoon, census_das_transformed,
                                fun = 'mean', progress = FALSE)
  
  daily_results[[i]] <- data.frame(
    DGUID              = census_das_transformed$DGUID,
    date               = current_date,
    pm25_mean          = da_24hr,        
    pm25_mean_afternoon = da_afternoon,  
    n_bands_primary    = length(BAND_PRIMARY),
    n_bands_afternoon  = length(BAND_AFTERNOON)
  )
  
  if(i %% 50 == 0) {
    cat("Processed", i, "of", nrow(file_df), "— Date:", format(current_date), "\n")
  }
}

da_daily_exposure <- bind_rows(daily_results) %>%
  mutate(
    above_15   = pm25_mean > 15,
    above_37.5 = pm25_mean > 37.5,
    above_50   = pm25_mean > 50
  )


# SAVE OUTPUTS

saveRDS(da_daily_exposure, paste0("da_daily_pm25_exposure_", year, ".rds"))
write.csv(da_daily_exposure, paste0("da_daily_pm25_exposure_", year, ".csv"),
          row.names = FALSE)

# ============================================================
# ANNUAL SUMMARY GPKG
# Required downstream for population data and spatial geometry
# ============================================================

da_annual_summary <- da_daily_exposure %>%
  group_by(DGUID) %>%
  summarise(
    annual_mean_pm25       = mean(pm25_mean,             na.rm = TRUE),
    annual_mean_pm25_afternoon = mean(pm25_mean_afternoon, na.rm = TRUE),
    days_above_15          = sum(above_15,               na.rm = TRUE),
    days_above_37.5        = sum(above_37.5,             na.rm = TRUE),
    days_above_50          = sum(above_50,               na.rm = TRUE),
    n_days                 = n(),
    .groups = "drop"
  ) %>%
  mutate(year = year)

# Join population and geometry from census file
da_annual_spatial <- census_das %>%
  left_join(da_annual_summary, by = "DGUID")

st_write(da_annual_spatial,
         paste0("da_annual_pm25_summary_", year, ".gpkg"),
         delete_layer = TRUE)
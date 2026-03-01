library(tidyverse)
library(sf)

DATA_DIR       <- "PM25data"         
FILE_PATTERN   <- "da_annual_pm25_summary_{year}.gpkg"
YEARS          <- 2021:2025
DGUID_COL      <- "DGUID" 
PM25_COL       <- "annual_mean_pm25"
PM25_FILE    <- "pm25_da.csv" 
OUTPUT_PM25    <- "pm25_da.csv"  
OUTPUT_CENSUS <- "census_da.csv"
BCSTATS_FILE <- "bcstats_proj.csv"

POP_COLS <- c(
  "pop_total_45_49", "pop_total_50_54",   
  "pop_total_55_59", "pop_total_60_64",   
  "pop_total_65_69", "pop_total_70_74",   
  "pop_total_75_79", "pop_total_80_84",   
  "pop_total_65plus"                      
)

all_years_pm25 <- map_dfr(YEARS, function(yr) {
  fpath <- file.path(DATA_DIR, gsub("\\{year\\}", yr, FILE_PATTERN))
  
  if (!file.exists(fpath)) {
    warning("File not found, skipping: ", fpath)
    return(NULL)
  }
  
  layer <- st_read(fpath, quiet = TRUE) %>%
    st_drop_geometry() %>%
    select(all_of(c(DGUID_COL, PM25_COL))) %>%
    rename(DGUID = all_of(DGUID_COL),
           pm25  = all_of(PM25_COL)) %>%
    mutate(year = yr)
})
# =============================================================================
# COMPUTE PER-DA MEAN PM2.5 ACROSS 2021–2025
year_counts <- all_years_pm25 %>%
  group_by(DGUID) %>%
  summarise(n_years = n_distinct(year), .groups = "drop")

incomplete <- filter(year_counts, n_years < length(YEARS))
if (nrow(incomplete) > 0) {
  cat("\nNote:", nrow(incomplete), "DAs have fewer than", length(YEARS),
      "years of data — included using available years.\n")
}

pm25_mean <- all_years_pm25 %>%
  group_by(DGUID) %>%
  summarise(
    mean_pm25    = mean(pm25, na.rm = TRUE),
    n_years_used = sum(!is.na(pm25)),
    min_annual   = min(pm25, na.rm = TRUE),
    max_annual   = max(pm25, na.rm = TRUE),
    sd_annual    = sd(pm25, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(dguid = DGUID)

cat("\n=== PM2.5 Exposure Summary (2021–2025 mean per DA) ===\n")
cat("DAs:               ", nrow(pm25_mean), "\n")
cat("Mean PM2.5:        ", round(mean(pm25_mean$mean_pm25), 2), "µg/m³\n")
cat("Median PM2.5:      ", round(median(pm25_mean$mean_pm25), 2), "µg/m³\n")
cat("Max PM2.5:         ", round(max(pm25_mean$mean_pm25), 2), "µg/m³\n")
cat("DAs above 8 µg/m³: ", sum(pm25_mean$mean_pm25 > 8),
    "(", round(100 * mean(pm25_mean$mean_pm25 > 8), 1), "%)\n")
cat("DAs above 15 µg/m³:", sum(pm25_mean$mean_pm25 > 15),
    "(", round(100 * mean(pm25_mean$mean_pm25 > 15), 1), "%)\n")

# Year-by-year summary (handy for your methods section)
cat("\nYear-by-year BC-wide mean PM2.5:\n")
all_years_pm25 %>%
  group_by(year) %>%
  summarise(bc_mean = round(mean(pm25, na.rm = TRUE), 2), .groups = "drop") %>%
  print()

# =============================================================================
# EXTRACT POPULATION COHORT DATA FROM 2021 FILE

census_file_2021 <- file.path(DATA_DIR, gsub("\\{year\\}", 2021, FILE_PATTERN))

census_da <- st_read(census_file_2021, quiet = TRUE) %>%
  st_drop_geometry() %>%
  select(all_of(c(DGUID_COL, POP_COLS))) %>%
  rename(dguid = all_of(DGUID_COL)) %>%
  mutate(
    # Pre-aggregate into the four cohorts the projection script expects
    pop_45_54 = pop_total_45_49 + pop_total_50_54,
    pop_55_64 = pop_total_55_59 + pop_total_60_64,
    pop_65_74 = pop_total_65_69 + pop_total_70_74,
    pop_75_84 = pop_total_75_79 + pop_total_80_84
  ) %>%
  select(dguid, pop_45_54, pop_55_64, pop_65_74, pop_75_84, pop_total_65plus)

cat("Population cohort data extracted for", nrow(census_da), "DAs\n")
cat("\nCohort sizes (BC total):\n")
census_da %>%
  summarise(across(starts_with("pop_"), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "cohort", values_to = "total_population") %>%
  print()

write_csv(pm25_mean, OUTPUT_PM25)
write_csv(census_da, OUTPUT_CENSUS)


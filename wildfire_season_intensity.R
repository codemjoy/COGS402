library(dplyr)
library(ggplot2)

# Read all years
all_years <- bind_rows(
  read.csv("bluesky_comprehensive_2021.csv"),  
  read.csv("bluesky_comprehensive_2022.csv"),  
  read.csv("bluesky_comprehensive_2023.csv"),  
  read.csv("bluesky_comprehensive_2024.csv"), 
  read.csv("bluesky_comprehensive_2025.csv"),  
  read.csv("bluesky_comprehensive_2026.csv")  
)

all_years$date <- as.Date(all_years$date)
all_years <- all_years %>%
  filter(!is.na(year), success == TRUE)


fire_season <- all_years %>%
  mutate(month = month(date), 
         month_name = month(date, label = TRUE)) %>%
  filter(month %in% c(5, 6, 7, 8, 9))


average_fire_season <- fire_season %>%
  group_by(year, month, month_name) %>%
  summarise(
    mean_daily_peak_24hr = mean(peak_max_24hr, na.rm = TRUE),
    max_daily_peak_24hr = max(peak_max_24hr, na.rm = TRUE),
    p95_peak_24hr = quantile(peak_max_24hr, 0.95, na.rm = TRUE),
    mean_daily_mean_24hr = mean(mean_mean_24hr, na.rm = TRUE),
    
    mean_area_over_15_km2 = mean(area_over_15_km2, na.rm = TRUE),
    max_area_over_15_km2 = max(area_over_15_km2, na.rm = TRUE),
    mean_area_over_37.5_km2 = mean(area_over_37.5_km2, na.rm = TRUE),
    max_area_over_37.5_km2 = max(area_over_37.5_km2, na.rm = TRUE),
    mean_area_over_50_km2 = mean(area_over_50_km2, na.rm = TRUE),
    max_area_over_50_km2 = max(area_over_50_km2, na.rm = TRUE),
    
    mean_pct_over_15 = mean(pct_bc_over_15, na.rm = TRUE),
    max_pct_over_15 = max(pct_bc_over_15, na.rm = TRUE),
    mean_pct_over_37.5 = mean(pct_bc_over_37.5, na.rm = TRUE),
    max_pct_over_37.5 = max(pct_bc_over_37.5, na.rm = TRUE),
    mean_pct_over_50 = mean(pct_bc_over_50, na.rm = TRUE),
    max_pct_over_50 = max(pct_bc_over_50, na.rm = TRUE),
    
    mean_consecutive_15 = mean(consecutive_hrs_over_15, na.rm = TRUE),
    max_consecutive_15 = max(consecutive_hrs_over_15),
    mean_consecutive_37.5 = mean(consecutive_hrs_over_37.5, na.rm = TRUE),
    max_consecutive_37.5 = max(consecutive_hrs_over_37.5),
    mean_consecutive_50 = mean(consecutive_hrs_over_50, na.rm = TRUE),
    max_consecutive_50 = max(consecutive_hrs_over_50),
    
    total_days = n(),
    .groups = "drop"
  ) %>%
  arrange(year)

print(average_fire_season)
wildfire_season_intensity <- paste0("wildfire_season_intensity.csv")
write.csv(average_fire_season, wildfire_season_intensity, row.names = FALSE)


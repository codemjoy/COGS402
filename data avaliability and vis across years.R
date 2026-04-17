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

# Resolution distribution

resolution_check <- all_years %>%
  group_by(year, resolution_km, data_source) %>%
  summarise(
    n_days = n(),
    .groups = "drop"
  )

print(resolution_check)

# 1. Check data completeness
data_coverage <- all_years %>%
  group_by(year) %>%
  summarise(
    total_files = n(),
    date_range = paste(min(date, na.rm = TRUE), "to", max(date, na.rm = TRUE)),
    unique_dates = n_distinct(date, na.rm = TRUE),
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE),
    .groups = "drop"
  )

print(data_coverage)

# Check monthly distribution
monthly_coverage <- all_years %>%
  mutate(
    month = month(date, label = TRUE),
    year_month = format(date, "%Y-%m")
  ) %>%
  group_by(year, month) %>%
  summarise(
    n_days = n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = month, values_from = n_days, values_fill = 0)

print(monthly_coverage)

# Visualize coverage
coverage_plot <- all_years %>%
  mutate(
    month = month(date, label = TRUE), 
    year_fac = as.factor(year)
  ) %>%
  group_by(year_fac, month) %>%
  summarize(n_days = n(), .groups = "drop")

p_coverage <- ggplot(coverage_plot, aes(x = month, y = year_fac, fill = n_days)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = n_days), color = "white", fontface = "bold") +
  scale_fill_gradient(low = "lightblue", high = "brown4", name = "Days") +
  labs(title = "BlueSky Data Coverage by Month and Year",
       x = "Month", y = "Year") +
  theme_minimal() +
  theme(panel.grid = element_blank())

print(p_coverage)
ggsave("data_coverage_heatmap.png", p_coverage, width = 12, height = 6)


# 2. Fire season timing 

summer_check <-all_years %>%
  mutate(month = month(date)) %>%
  filter(month %in% c(7, 8)) %>%
  group_by(year, month) %>%
  summarise(
    n_days = n(),
    mean_peak_24hr = mean(peak_max_24hr, na.rm = TRUE),
    max_peak_24hr = max(peak_max_24hr, na.rm = TRUE),
    mean_peak_51hr = max(peak_max_51hr, na.rm = TRUE),
    max_peak_51hr = max(peak_max_51hr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year, month)

print(summer_check)

summer_data <- all_years %>%
  mutate(month = month(date)) %>%
  filter(month %in% c(6, 7, 8, 9))

p_summer_fire_season <- ggplot(summer_data, aes(x = date, y = peak_max_24hr)) +
  geom_line(alpha = 0.3, color = "gray50") +
  geom_point(aes(color = peak_max_24hr > 37.5), alpha = 0.6, size = 2) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "orange",
             alpha = 0.7, linewidth = 0.8) +
  geom_hline(yintercept = 37.5, linetype = "dashed", color = "red",
             alpha = 0.7, linewidth = 0.8) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "darkred",
             alpha = 0.7, linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                     labels = c("Below 37.5 µg/m³", "Above 37.5 µg/m³"), 
                     name = NULL) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000),
    labels = c("10", "50", "100", "500", "1000", "5000", "10000",
               "50,000", "100,000")
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") + 
  facet_wrap(~year, scales = "free_x", ncol = 1) +
  labs(title = "Daily Peak PM2.5 (24-hour) Across BC Summer (June- September, 2021-2025)",
       subtitle = "Reference lines: Orange = WHO AQG (15 µg/m³) | Red = WHO IT-3 (37.5 µg/m³) | Dark Red = WHO IT-2 (50 µg/m³)",
       x = "Date", 
       y = "Peak PM2.5 (µg/m³, log scale)") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"),
        plot.subtitle = element_text(size = 20, color = "gray30", hjust = 0.5),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        strip.text = element_text(size = 15, face = "bold")
        )

print(p_summer_fire_season)
ggsave("summer_fire_season_pm25.png", p_summer_fire_season, width = 20, height = 18)


# 3. Spatial Coverage vs intensity

spatial_comparison <- all_years %>%
  group_by(year) %>%
  summarise(
    mean_daily_peak_24hr = mean(peak_max_24hr, na.rm = TRUE),
    max_daily_peak_24hr = max(peak_max_24hr, na.rm = TRUE),
    p95_peak_24hr = quantile(peak_max_24hr, 0.95, na.rm = TRUE),
    mean_daily_mean_24hr = mean(mean_mean_24hr, na.rm = TRUE),
    
    mean_area_over_15_km2 = mean(area_over_15_km2, na.rm = TRUE),
    mean_area_over_15_km2 = mean(area_over_15_km2, na.rm = TRUE),
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

print(spatial_comparison)

p_intensity_extent <- ggplot(spatial_comparison, aes(x = mean_pct_over_15, y = mean_daily_peak_24hr)) +
  geom_point(aes(size = max_daily_peak_24hr, color = as.factor(year)), alpha = 0.7) +
  geom_text(aes(label = year), vjust = -1.5, size = 4, fontface = "bold") +
  scale_size_continuous(name = "Max Peak\nPM2.5", range = c(5, 15)) +
  scale_color_viridis_d(name = "Year") +
  labs(title = "Spatial Extent vs Intensity of PM2.5 Pollution",
       subtitle = "Size of points = maximum 24-hour peak observed",
       x = "Mean % of BC Exceeding 15 µg/m³",
       y = "Mean Daily Peak PM2.5 (µg/m³, 24-hr)") +
  theme_minimal()

print(p_intensity_extent)
ggsave("intensity_vs_extent_15.png", p_intensity_extent, width = 10, height = 7)

p_intensity_extent <- ggplot(spatial_comparison, aes(x = mean_pct_over_37.5, y = mean_daily_peak_24hr)) +
  geom_point(aes(size = max_daily_peak_24hr, color = as.factor(year)), alpha = 0.7) +
  geom_text(aes(label = year), vjust = -1.5, size = 4, fontface = "bold") +
  scale_size_continuous(name = "Max Peak\nPM2.5", range = c(5, 15)) +
  scale_color_viridis_d(name = "Year") +
  labs(title = "Spatial Extent vs Intensity of PM2.5 Pollution",
       subtitle = "Size of points = maximum 24-hour peak observed",
       x = "Mean % of BC Exceeding 37.5 µg/m³",
       y = "Mean Daily Peak PM2.5 (µg/m³, 24-hr)") +
  theme_minimal()

print(p_intensity_extent)
ggsave("intensity_vs_extent_37.5.png", p_intensity_extent, width = 10, height = 7)

p_intensity_extent <- ggplot(spatial_comparison, aes(x = mean_pct_over_50, y = mean_daily_peak_24hr)) +
  geom_point(aes(size = max_daily_peak_24hr, color = as.factor(year)), alpha = 0.7) +
  geom_text(aes(label = year), vjust = -1.5, size = 4, fontface = "bold") +
  scale_size_continuous(name = "Max Peak\nPM2.5", range = c(5, 15)) +
  scale_color_viridis_d(name = "Year") +
  labs(title = "Spatial Extent vs Intensity of PM2.5 Pollution",
       subtitle = "Size of points = maximum 24-hour peak observed",
       x = "Mean % of BC Exceeding 50 µg/m³",
       y = "Mean Daily Peak PM2.5 (µg/m³, 24-hr)") +
  theme_minimal()

print(p_intensity_extent)
ggsave("intensity_vs_extent_50.png", p_intensity_extent, width = 10, height = 7)

# 4. Check Peak Fire Days
worst_days_all_15 <- all_years %>%
  arrange(desc(peak_max_24hr)) %>%
  select(date, year, peak_max_24hr, mean_mean_24hr, consecutive_hrs_over_15,
         area_over_15_km2, pct_bc_over_15) %>%
  head(20)

print(worst_days_all_15)

worst_days_file_15 <- paste0("worst_days_all_15", year, ".csv")
write.csv(worst_days_all_15, worst_days_file_15, row.names = FALSE)


worst_by_year_15 <- all_years %>%
  group_by(year) %>%
  slice_max(peak_max_24hr, n = 5) %>%
  select(year, date, peak_max_24hr, mean_mean_24hr, consecutive_hrs_over_15,
         pct_bc_over_15, area_over_15_km2) %>%
  arrange(year, desc(peak_max_24hr))

print(worst_by_year_15)
worst_days_by_year_file_15 <- paste0("worst_days_by_year_15", year, ".csv")
write.csv(worst_by_year_15, worst_days_by_year_file_15, row.names = FALSE)

extreme_events_15 <- all_years %>%
  filter(peak_max_24hr > 15) %>%
  group_by(year) %>%
  summarise(
    n_extreme_days = n(),
    mean_peak_on_extreme_days = mean(peak_max_24hr, na.rm = TRUE),
    max_peak = max(peak_max_24hr, na.rm = TRUE),
    mean_duration = mean(consecutive_hrs_over_15, na.rm = TRUE),
    mean_area_affected = mean(area_over_15_km2, na.rm = TRUE),
    mean_pct_bc_affected = mean(pct_bc_over_15, na.rm = TRUE),
    .groups = "drop"
  )

print(extreme_events_15)
extreme_events_file_15 <- paste0("extreme_events_file_15", year, ".csv")
write.csv(extreme_events_15, extreme_events_file_15, row.names = FALSE)

worst_days_all_37.5 <- all_years %>%
  arrange(desc(peak_max_24hr)) %>%
  select(date, year, peak_max_24hr, mean_mean_24hr, consecutive_hrs_over_37.5,
         area_over_37.5_km2, pct_bc_over_37.5) %>%
  head(20)

print(worst_days_all_37.5)
worst_days_file_37.5 <- paste0("worst_days_all_37.5", year, ".csv")
write.csv(worst_days_all_37.5, worst_days_file_37.5, row.names = FALSE)

worst_by_year_37.5 <- all_years %>%
  group_by(year) %>%
  slice_max(peak_max_24hr, n = 5) %>%
  select(year, date, peak_max_24hr, mean_mean_24hr, consecutive_hrs_over_37.5,
         pct_bc_over_37.5, area_over_37.5_km2) %>%
  arrange(year, desc(peak_max_24hr))

print(worst_by_year_37.5)
worst_days_by_year_file_37.5 <- paste0("worst_days_by_year_37.5", year, ".csv")
write.csv(worst_by_year_37.5, worst_days_by_year_file_37.5, row.names = FALSE)

extreme_events_37.5 <- all_years %>%
  filter(peak_max_24hr > 37.5) %>%
  group_by(year) %>%
  summarise(
    n_extreme_days = n(),
    mean_peak_on_extreme_days = mean(peak_max_24hr, na.rm = TRUE),
    max_peak = max(peak_max_24hr, na.rm = TRUE),
    mean_duration = mean(consecutive_hrs_over_37.5, na.rm = TRUE),
    mean_area_affected = mean(area_over_37.5_km2, na.rm = TRUE),
    mean_pct_bc_affected = mean(pct_bc_over_37.5, na.rm = TRUE),
    .groups = "drop"
  )

print(extreme_events_37.5)
extreme_events_file_37.5 <- paste0("extreme_events_file_37.5", year, ".csv")
write.csv(extreme_events_37.5, extreme_events_file_37.5, row.names = FALSE)

worst_days_all_50 <- all_years %>%
  arrange(desc(peak_max_24hr)) %>%
  select(date, year, peak_max_24hr, mean_mean_24hr, consecutive_hrs_over_50,
         area_over_50_km2, pct_bc_over_50) %>%
  head(20)

print(worst_days_all_50)
worst_days_file_15 <- paste0("worst_days_all_15", year, ".csv")
write.csv(worst_days_all_15, worst_days_file_15, row.names = FALSE)

worst_by_year_50 <- all_years %>%
  group_by(year) %>%
  slice_max(peak_max_24hr, n = 5) %>%
  select(year, date, peak_max_24hr, mean_mean_24hr, consecutive_hrs_over_50,
         pct_bc_over_50, area_over_50_km2) %>%
  arrange(year, desc(peak_max_24hr))

print(worst_by_year_50)
worst_days_by_year_file_50 <- paste0("worst_days_by_year_50", year, ".csv")
write.csv(worst_by_year_50, worst_days_by_year_file_50, row.names = FALSE)


extreme_events_50 <- all_years %>%
  filter(peak_max_24hr > 50) %>%
  group_by(year) %>%
  summarise(
    n_extreme_days = n(),
    mean_peak_on_extreme_days = mean(peak_max_24hr, na.rm = TRUE),
    max_peak = max(peak_max_24hr, na.rm = TRUE),
    mean_duration = mean(consecutive_hrs_over_50, na.rm = TRUE),
    mean_area_affected = mean(area_over_50_km2, na.rm = TRUE),
    mean_pct_bc_affected = mean(pct_bc_over_50, na.rm = TRUE),
    .groups = "drop"
  )

print(extreme_events_50)

extreme_events_file_50 <- paste0("extreme_events_file_50", year, ".csv")
write.csv(extreme_events_50, extreme_events_file_50, row.names = FALSE)


# 5. TIME SERIES VISUALIZATION

# Full time series - Peak PM2.5 37.5 
p_timeseries_peak <- ggplot(all_years, aes(x = date, y = peak_max_24hr)) +
  geom_line(alpha = 0.3, color = "gray50") +
  geom_point(aes(color = peak_max_24hr > 37.5), alpha = 0.6, size = 1) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "orange", alpha = 0.7) +
  geom_hline(yintercept = 37.5, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "darkred", alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                     name = "Exceeds\n37.5 µg/m³") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 10, 50, 100, 500, 1000, 5000, 10000),
    labels = c("0", "1", "10", "50", "100", "500", "1,000", "5,000", "10,000")
                                     ) +
  facet_wrap(~year, scales = "free_x", ncol = 1) +
  labs(title = "Daily Peak PM2.5 (24-hour) Across BC (2021-2026)",
       subtitle = "Dashed lines: WHO AQG (15), IT-3 (37.5), IT-2 (50)",
       x = "Date", 
       y = "Peak PM2.5 (µg/m³, log scale)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"),
        plot.subtitle = element_text(size = 20, color = "gray30", hjust = 0.5),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        strip.text = element_text(size = 15, face = "bold"))

print(p_timeseries_peak)
ggsave("timeseries_peak_pm25_24hr.png", p_timeseries_peak, width = 20, height = 16)

# Full time series - Mean PM2.5 (24-hr equivalent)
p_timeseries_mean <- ggplot(all_years, aes(x = date, y = mean_mean_24hr)) +
  geom_line(alpha = 0.3, color = "gray50") +
  geom_point(aes(color = mean_mean_24hr > 15), alpha = 0.6, size = 1) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "orange", alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                     name = "Exceeds\nWHO AQG") +
  facet_wrap(~year, scales = "free_x", ncol = 1) +
  labs(title = "Daily Mean PM2.5 (24-hour average) Across BC (2021-2026)",
       subtitle = "Dashed line: WHO AQG (15 µg/m³)",
       x = "Date", 
       y = "Mean PM2.5 (µg/m³)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"),
        plot.subtitle = element_text(size = 20, color = "gray30", hjust = 0.5),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        strip.text = element_text(size = 15, face = "bold"))

print(p_timeseries_mean)
ggsave("timeseries_mean_pm25_24hr.png", p_timeseries_mean, width = 20, height = 16)

library(lubridate)
library(viridis)

# Duration time series



p_timeseries_duration <- ggplot(all_years, aes(x = date, y = consecutive_hrs_over_37.5)) +
  geom_line(alpha = 0.3, color = "gray50") +
  geom_point(aes(color = consecutive_hrs_over_37.5, size = consecutive_hrs_over_37.5),
             alpha = 0.7) +
  scale_color_gradient2(
    low = "chartreuse4",
    mid = "goldenrod1",
    high = "darkred",
  ) +
  scale_size_continuous(range = c(0.5, 4), name = "Duration \n(hours)", 
                        breaks = c(0, 12, 25, 40, 50) )+
  scale_x_date(date_labels = "%b", date_breaks = "2 months", limits = function(x) {
    c(floor_date(min(x), "year"), ceiling_date(max(x), "year") - days(1))
  }) +
  facet_wrap(~year, scales = "free", ncol = 1) +
  labs(title = "Duration of Elevated PM2.5 Events (>37.5 µg/m³)",
       subtitle = "Consecutive hours above WHO IT-3 threshold",
       x = "Date", 
       y = "Consecutive Hours > 37.5 µg/m³") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(1.5, "cm"),
        plot.subtitle = element_text(size = 20, color = "gray30", hjust = 0.5),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 18),
        strip.text = element_text(size = 15, face = "bold"),
        panel.spacing = unit(1, "lines")
        )

print(p_timeseries_duration)
ggsave("timeseries_duration.png", p_timeseries_duration, width = 20, height = 17)

# Spatial extent over time
p_timeseries_spatial <- ggplot(all_years, aes(x = date, y = pct_bc_over_37.5)) +
  geom_line(alpha = 0.3, color = "gray50") +
  geom_point(aes(color = pct_bc_over_37.5 > 10), alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                     name = "> 10% of BC\nAffected") +
  facet_wrap(~year, scales = "free_x", ncol = 2) +
  labs(title = "Spatial Extent of Unhealthy Air Quality",
       subtitle = "Percentage of BC exceeding 37.5 µg/m³ (WHO IT-3)",
       x = "Date", 
       y = "% of BC Exceeding 37.5 µg/m³") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_timeseries_spatial)
ggsave("timeseries_spatial_extent.png", p_timeseries_spatial, width = 14, height = 10)

p_timeseries_area <- ggplot(all_years, aes(x = date, y = area_over_37.5_km2)) +
  geom_line(alpha = 0.3, color = "gray50") +
  geom_point(aes(color = area_over_37.5_km2 > 50000), alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                     name = "> 50,000 km²\nAffected") +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~year, scales = "free_x", ncol = 2) +
  labs(title = "Area of BC with Unhealthy Air Quality",
       subtitle = "Area (km²) exceeding 37.5 µg/m³ - resolution-normalized",
       x = "Date", 
       y = "Area (km²) Exceeding 37.5 µg/m³") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_timeseries_area)
ggsave("timeseries_area_extent.png", p_timeseries_area, width = 14, height = 10)


# 6. FIRE SEASON COMPARISON (MAY-SEPTEMBER)

fire_season_data <- all_years %>%
  mutate(month = month(date)) %>%
  filter(month >= 5 & month <= 9)

# Monthly patterns during fire season
fire_season_monthly <- fire_season_data %>%
  mutate(month_name = month(date, label = TRUE)) %>%
  group_by(year, month_name) %>%
  summarise(
    n_days = n(),
    mean_peak_24hr = mean(peak_max_24hr, na.rm = TRUE),
    max_peak_24hr = max(peak_max_24hr, na.rm = TRUE),
    days_over_37.5 = sum(peak_max_24hr > 37.5, na.rm = TRUE),
    mean_consecutive = mean(consecutive_hrs_over_37.5, na.rm = TRUE),
    mean_area_affected = mean(area_over_37.5_km2, na.rm = TRUE),
    .groups = "drop"
  )

print("Fire Season Monthly Summary:")
print(fire_season_monthly)

# Visualize fire season progression
p_fire_season <- ggplot(fire_season_monthly, aes(x = month_name, y = mean_peak_24hr, 
                                                 group = year, color = as.factor(year))) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_viridis_d(name = "Year") +
  labs(title = "Fire Season Progression (May-September)",
       subtitle = "Mean daily peak PM2.5 (24-hour) by month",
       x = "Month",
       y = "Mean Peak PM2.5 (µg/m³)") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_fire_season)
ggsave("fire_season_progression.png", p_fire_season, width = 10, height = 6)

# Peak intensity by month
p_fire_season_max <- ggplot(fire_season_monthly, aes(x = month_name, y = max_peak_24hr, 
                                                     fill = as.factor(year))) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(name = "Year") +
  labs(title = "Peak PM2.5 Events by Month (Fire Season)",
       subtitle = "Maximum 24-hour peak observed each month",
       x = "Month",
       y = "Maximum Peak PM2.5 (µg/m³)") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_fire_season_max)
ggsave("fire_season_max_by_month.png", p_fire_season_max, width = 10, height = 6)

# Days exceeding thresholds by month
p_fire_season_days <- ggplot(fire_season_monthly, aes(x = month_name, y = days_over_37.5, 
                                                      fill = as.factor(year))) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(name = "Year") +
  labs(title = "Days Exceeding WHO IT-3 (37.5 µg/m³) During Fire Season",
       x = "Month",
       y = "Number of Days > 37.5 µg/m³") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_fire_season_days)
ggsave("fire_season_exceedance_days.png", p_fire_season_days, width = 10, height = 6)

# 7. COMPREHENSIVE ANNUAL COMPARISON

annual_summary <- all_years %>%
  group_by(year) %>%
  summarise(
    total_days = n(),
    resolution_mix = paste(unique(sort(resolution_km)), collapse = ", "),
    
    # Peak metrics
    highest_peak_24hr = max(peak_max_24hr, na.rm = TRUE),
    mean_daily_peak_24hr = mean(peak_max_24hr, na.rm = TRUE),
    median_daily_peak_24hr = median(peak_max_24hr, na.rm = TRUE),
    p95_peak_24hr = quantile(peak_max_24hr, 0.95, na.rm = TRUE),
    
    # Mean metrics (24-hr)
    highest_mean_24hr = max(mean_mean_24hr, na.rm = TRUE),
    mean_daily_mean_24hr = mean(mean_mean_24hr, na.rm = TRUE),
    
    # Duration
    max_consecutive_37.5 = max(consecutive_hrs_over_37.5, na.rm = TRUE),
    mean_consecutive_37.5 = mean(consecutive_hrs_over_37.5, na.rm = TRUE),
    
    # Frequency - Peak
    days_peak_over_15 = sum(peak_max_24hr > 15, na.rm = TRUE),
    days_peak_over_35 = sum(peak_max_24hr > 37.5, na.rm = TRUE),
    days_peak_over_50 = sum(peak_max_24hr > 50, na.rm = TRUE),
    
    # Frequency - Mean (24-hr)
    days_mean_over_15 = sum(mean_mean_24hr > 15, na.rm = TRUE),
    days_mean_over_37.5 = sum(mean_mean_24hr > 37.5, na.rm = TRUE),
    
    # Spatial extent
    mean_pct_over_15 = mean(pct_bc_over_15, na.rm = TRUE),
    mean_pct_over_37.5 = mean(pct_bc_over_37.5, na.rm = TRUE),
    max_pct_over_37.5 = max(pct_bc_over_37.5, na.rm = TRUE),
    
    mean_area_over_15_km2 = mean(area_over_15_km2, na.rm = TRUE),
    max_area_over_15_km2 = max(area_over_15_km2, na.rm = TRUE),
    mean_area_over_37.5_km2 = mean(area_over_37.5_km2, na.rm = TRUE),
    max_area_over_37.5_km2 = max(area_over_37.5_km2, na.rm = TRUE),
    
    # Prolonged exposure
    mean_area_24hrs_over_37.5_km2 = mean(area_24hrs_over_37.5_km2, na.rm = TRUE),
    days_with_prolonged_exposure = sum(area_24hrs_over_37.5_km2 > 10000, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(year)

print(annual_summary)
write.csv(annual_summary, "annual_comprehensive_summary.csv", row.names = FALSE)

# ============================================
# 8. FINAL COMPARISON VISUALIZATIONS
# ============================================

# Peak comparison
p_annual_peaks <- ggplot(annual_summary, aes(x = as.factor(year))) +
  geom_col(aes(y = mean_daily_peak_24hr), fill = "steelblue", alpha = 0.7) +
  geom_point(aes(y = highest_peak_24hr), color = "red", size = 5) +
  geom_text(aes(y = highest_peak_24hr, label = round(highest_peak_24hr, 0)), 
            vjust = -0.7, fontface = "bold") +
  labs(title = "Annual PM2.5 Peak Comparison (24-hour)",
       subtitle = "Blue bars = mean daily peak | Red points = highest peak observed",
       x = "Year",
       y = "PM2.5 (µg/m³)") +
  theme_minimal()

print(p_annual_peaks)
ggsave("annual_peaks_comparison_24hr.png", p_annual_peaks, width = 10, height = 6)

# Area comparison (resolution-normalized)
p_annual_area <- ggplot(annual_summary, aes(x = as.factor(year))) +
  geom_col(aes(y = mean_area_over_37.5_km2), fill = "coral", alpha = 0.7) +
  geom_point(aes(y = max_area_over_37.5_km2), color = "darkred", size = 5) +
  geom_text(aes(y = max_area_over_37.5_km2, 
                label = scales::comma(round(max_area_over_37.5_km2, 0))), 
            vjust = -0.7, fontface = "bold", size = 3) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Annual Area Exceeding WHO IT-3 (37.5 µg/m³)",
       subtitle = "Coral bars = mean daily area | Red points = maximum area (resolution-normalized)",
       x = "Year",
       y = "Area (km²)") +
  theme_minimal()

print(p_annual_area)
ggsave("annual_area_comparison.png", p_annual_area, width = 10, height = 6)

# Exceedance days comparison
exceedance_long <- annual_summary %>%
  select(year, days_peak_over_15, days_peak_over_35, days_peak_over_50) %>%
  tidyr::pivot_longer(cols = starts_with("days"), 
                      names_to = "threshold",
                      values_to = "count") %>%
  mutate(threshold = factor(threshold,
                            levels = c("days_peak_over_100", "days_peak_over_35", "days_peak_over_15"),
                            labels = c(">50 µg/m³", ">37.5 µg/m³", ">15 µg/m³")))

p_exceedances <- ggplot(exceedance_long, aes(x = as.factor(year), y = count, fill = threshold)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(">15 µg/m³" = "orange", 
                               ">35 µg/m³" = "red", 
                               ">50 µg/m³" = "darkred")) +
  labs(title = "Days Exceeding PM2.5 Thresholds",
       x = "Year",
       y = "Number of Days",
       fill = "Threshold") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_exceedances)
ggsave("annual_exceedances_comparison.png", p_exceedances, width = 10, height = 6)

print("\n========================================")
print("ANALYSIS COMPLETE!")
print("========================================")

# Convert date
all_years$date <- as.Date(all_years$date)

# Annual summary with correct column names
annual_summary <- all_years %>%
  filter(success == TRUE) %>%
  group_by(year) %>%
  summarise(
    # Peak metrics
    highest_peak = max(peak_max_24hr, na.rm = TRUE),
    mean_daily_peak = mean(peak_max_24hr, na.rm = TRUE),
    median_daily_peak = median(peak_max_24hr, na.rm = TRUE),
    
    # Mean metrics (24-hr average)
    highest_mean = max(mean_mean_24hr, na.rm = TRUE),
    mean_daily_mean = mean(mean_mean_24hr, na.rm = TRUE),
    
    # Duration metrics
    mean_consecutive_37.5 = mean(consecutive_hrs_over_37.5, na.rm = TRUE),
    max_consecutive_37.5 = max(consecutive_hrs_over_37.5, na.rm = TRUE),
    mean_consecutive_15 = mean(consecutive_hrs_over_15, na.rm = TRUE),
    max_consecutive_15 = max(consecutive_hrs_over_15, na.rm = TRUE),
    
    # Frequency metrics - PEAK
    days_peak_over_15 = sum(peak_max_24hr > 15, na.rm = TRUE),
    days_peak_over_37.5 = sum(peak_max_24hr > 37.5, na.rm = TRUE),
    days_peak_over_50 = sum(peak_max_24hr > 50, na.rm = TRUE),
    
    # Frequency metrics - MEAN (24-hr)
    days_mean_over_15 = sum(mean_mean_24hr > 15, na.rm = TRUE),
    days_mean_over_37.5 = sum(mean_mean_24hr > 37.5, na.rm = TRUE),
    
    # Spatial extent
    mean_pct_over_15 = mean(pct_bc_over_15, na.rm = TRUE),
    mean_pct_over_37.5 = mean(pct_bc_over_37.5, na.rm = TRUE),
    
    # Prolonged exposure
    days_with_24hrs_over_37.5 = sum(cells_24hrs_over_37.5 > 0, na.rm = TRUE),
    
    total_days = n()
  ) %>%
  arrange(year)

print("========== ANNUAL COMPARISON ==========")
print(annual_summary)

# Save
write.csv(annual_summary, "annual_pm25_comprehensive.csv", row.names = FALSE)

# Visualizations
p1 <- ggplot(annual_summary, aes(x = as.factor(year))) +
  geom_col(aes(y = mean_daily_peak), fill = "steelblue", alpha = 0.7) +
  geom_point(aes(y = highest_peak), color = "red", size = 4) +
  geom_text(aes(y = highest_peak, label = round(highest_peak, 0)), 
            vjust = -0.5, size = 3) +
  labs(title = "Annual PM2.5 Peaks in BC",
       subtitle = "Blue bars = mean daily peak | Red points = highest peak",
       x = "Year",
       y = "PM2.5 (µg/m³)") +
  theme_minimal()

print(p1)
ggsave("annual_peaks.png", p1, width = 10, height = 6)



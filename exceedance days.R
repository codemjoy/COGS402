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

exceedance_days_annually <- all_years %>%
  group_by (year) %>%
  summarise (
    total_days = n(),
    
    # using 24-hour peak
    days_peak_over_15_AQG = sum(peak_max_24hr > 15, na.rm = TRUE),
    days_peak_over_37.5_IT3 = sum(peak_max_24hr > 37.5, na.rm = TRUE),
    days_peak_over_50_IT2 = sum(peak_max_24hr > 50, na.rm = TRUE),
    
    # using 24-hour mean 
    days_mean_over_15_AQG = sum(mean_mean_24hr > 15, na.rm = TRUE),
    days_mean_over_37.5_IT3 = sum(mean_mean_24hr > 37.5, na.rm = TRUE),
    days_mean_over_50_IT2 = sum(mean_mean_24hr > 50, na.rm = TRUE),
    
    # percentages
    pct_days_peak_over_15 = (sum(peak_max_24hr > 15, na.rm = TRUE) / n ()) * 100, 
    pct_days_peak_over_37.5 = (sum(peak_max_24hr > 37.5, na.rm = TRUE) / n()) * 100, 
    pct_days_peak_over_50 = (sum(peak_max_24hr > 50, na.rm = TRUE) / n()) * 100, 
    
    .groups = "drop"
  ) %>%
  arrange(year)

print(exceedance_days_annually)
write.csv(exceedance_days_annually, "exceedance_days_annually.csv", row.names = FALSE)


exceedance_days_monthly <- all_years %>%
  mutate(month = month(date),
         month_name = month(date, label = TRUE)) %>%
  group_by (year, month, month_name) %>%
  summarise (
    total_days = n(),
    
    # using 24-hour peak
    days_peak_over_15_AQG = sum(peak_max_24hr > 15, na.rm = TRUE),
    days_peak_over_37.5_IT3 = sum(peak_max_24hr > 37.5, na.rm = TRUE),
    days_peak_over_50_IT2 = sum(peak_max_24hr > 50, na.rm = TRUE),
    
    # using 24-hour mean 
    days_mean_over_15_AQG = sum(mean_mean_24hr > 15, na.rm = TRUE),
    days_mean_over_37.5_IT3 = sum(mean_mean_24hr > 37.5, na.rm = TRUE),
    days_mean_over_50_IT2 = sum(mean_mean_24hr > 50, na.rm = TRUE),
    
    # percentage of month
    pct_month_peak_over_15 = (sum(peak_max_24hr > 15, na.rm = TRUE) / n ()) * 100, 
    pct_month_peak_over_37.5 = (sum(peak_max_24hr > 37.5, na.rm = TRUE) / n()) * 100, 
    pct_month_peak_over_50 = (sum(peak_max_24hr > 50, na.rm = TRUE) / n()) * 100, 
    
    .groups = "drop"
  ) %>%
  arrange(year)

print(exceedance_days_monthly)
write.csv(exceedance_days_monthly, "exceedance_days_monthly.csv", row.names = FALSE)


# For the fire season only: 

exceedance_days_fire_season <- all_years %>%
  mutate(month = month(date)) %>%
  filter(month %in% c(5, 6, 7, 8, 9)) %>%
  group_by (year) %>%
  summarise (
    fire_season_days = n(),
    
    # using 24-hour peak
    days_peak_over_15_AQG = sum(peak_max_24hr > 15, na.rm = TRUE),
    days_peak_over_37.5_IT3 = sum(peak_max_24hr > 37.5, na.rm = TRUE),
    days_peak_over_50_IT2 = sum(peak_max_24hr > 50, na.rm = TRUE),
    
    # using 24-hour mean 
    days_mean_over_15_AQG = sum(mean_mean_24hr > 15, na.rm = TRUE),
    days_mean_over_37.5_IT3 = sum(mean_mean_24hr > 37.5, na.rm = TRUE),
    days_mean_over_50_IT2 = sum(mean_mean_24hr > 50, na.rm = TRUE),
    
    # percentage of month
    pct_fire_season_over_15 = (sum(peak_max_24hr > 15, na.rm = TRUE) / n ()) * 100, 
    pct_fire_season_over_37.5 = (sum(peak_max_24hr > 37.5, na.rm = TRUE) / n()) * 100, 
    pct_fire_season_over_50 = (sum(peak_max_24hr > 50, na.rm = TRUE) / n()) * 100, 
    
    .groups = "drop"
  ) %>%
  arrange(year)

print(exceedance_days_fire_season)
write.csv(exceedance_days_fire_season, "exceedance_days_fire_season.csv", row.names = FALSE)

# Comparison of three measures

comparision_summary <- exceedance_days_annually %>%
  select(year, total_days, days_peak_over_15_AQG, days_peak_over_37.5_IT3,
         days_peak_over_50_IT2) %>%
  pivot_longer(
    cols = starts_with("days"),
    names_to = "threshold",
    values_to = "count"
  ) %>%
  mutate(
    threshold = case_when(
      threshold == "days_peak_over_15_AQG" ~ "WHO AQG (15 µg/m³)",
      threshold == "days_peak_over_37.5_IT3" ~ "WHO IT-3 (37.5 µg/m³)",
      threshold == "days_peak_over_50_IT2" ~ "WHO IT-2 (50 µg/m³)",
      TRUE ~ threshold
    ),
    threshold = factor(threshold, 
                       levels = c("WHO AQG (15 µg/m³)",
                                  "WHO IT-3 (37.5 µg/m³)",
                                  "WHO IT-2 (50 µg/m³)"))
  )

print (comparision_summary)
write.csv(comparision_summary, "comparision_summary.csv", row.names = FALSE)


p_annual <- ggplot(comparision_summary, aes(x = as.factor(year), y = count, 
                                            fill = threshold)) +
  geom_col(position = "dodge") + 
  geom_text(aes(label = count), position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_manual(
    values = c("WHO AQG (15 µg/m³)" = "darkorange2", 
               "WHO IT-3 (37.5 µg/m³)" = "firebrick2",
               "WHO IT-2 (50 µg/m³)" = "darkred"),
    name = "Threshold"
  ) +
  labs(
    title = "Days Exceeding WHO Air Quality Guidelines",
    subtitle = "Based on 24-hour peak PM2.5",
    x = "Year",
    y = "Number of Days"
  ) + 
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "bottom",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )


print(p_annual)
ggsave("exceedance_days_annual_plot.png", p_annual, width = 12, height = 7, dpi = 300)

# plot for only fire season

fire_season_long <- exceedance_days_fire_season %>%
  select(year, days_peak_over_15_AQG, days_peak_over_37.5_IT3, days_peak_over_50_IT2) %>%
  pivot_longer(
    cols = starts_with("days_"),
    names_to = "threshold",
    values_to = "count"
  ) %>%
  mutate(
    threshold = case_when(
      threshold == "days_peak_over_15_AQG" ~ "WHO AQG (15 µg/m³)",
      threshold == "days_peak_over_37.5_IT3" ~ "WHO IT-3 (37.5 µg/m³)",
      threshold == "days_peak_over_50_IT2" ~ "WHO IT-2 (50 µg/m³)",
      TRUE ~ threshold
    ), 
    threshold = factor(threshold,
                       levels = c("WHO AQG (15 µg/m³)",
                                  "WHO IT-3 (37.5 µg/m³)",
                                  "WHO IT-2 (50 µg/m³)"))
  )


p_fire_season <- ggplot(fire_season_long, aes(x = as.factor(year), y = count, fill = threshold)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(
    values = c("WHO AQG (15 µg/m³)" = "chocolate1", 
               "WHO IT-3 (37.5 µg/m³)" = "firebrick2",
               "WHO IT-2 (50 µg/m³)" = "darkred"),
    name = "Threshold"
  ) +
  labs(
    title = "Fire Season Days Exceeding WHO Guidelines (May-September)",
    subtitle = "Based on 24-hour peak PM2.5",
    x = "Year",
    y = "Number of Days"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "bottom",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
  )

print(p_fire_season)
ggsave("exceedance_days_fire_season_plot.png", p_fire_season, width = 12, height = 7, dpi = 300)


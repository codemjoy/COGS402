library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggh4x)  
library(ggtext)
library(sf)

years <- 2021:2025

da_pop <- st_read("bc_das_population_for_raster.gpkg") %>%
  st_drop_geometry() %>%
  select(DGUID, population = pop_total_2021) %>%
  distinct()

cf_annual <- read_csv("annual_cf_stats.csv") %>%
  select(fire_year, cf_applied = cf_median)

ts_data <- map_dfr(years, function(yr) {
  readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds")) %>%
    mutate(pm25_mean = pmin(pm25_mean, 150),
           month_num = month(date),
           fire_year = yr
           ) %>%  
    left_join(cf_annual, by = "fire_year") %>%
    mutate(
      pm25_cf_removed = pmax(pm25_mean - cf_applied, 0)
    ) %>%
    left_join(da_pop, by = "DGUID") %>%
    group_by(date) %>%
    summarise(
      bc_daily_mean    = weighted.mean(pm25_mean,        w = population, na.rm = TRUE),
      bc_daily_mean_cf = weighted.mean(pm25_cf_removed,  w = population, na.rm = TRUE),
      cf_used          = first(cf_applied),   # <-- was cf_median, now cf_applied
      .groups = "drop"
    ) %>%
    mutate(
      year      = yr,
      plot_date = as.Date(paste0("2000-", format(date, "%m-%d")))
    ) %>%
    filter(month(date) >= 5 & month(date) <= 10)
}) %>%
  mutate(
    above_threshold = bc_daily_mean >= 15,
    year_label = as.character(year)
  ) 

ts_data %>%
  mutate(month_num = month(date)) %>%
  group_by(year, month_num) %>%
  summarise(
    cf_used       = first(cf_used),
    mean_raw      = round(mean(bc_daily_mean,    na.rm = TRUE), 2),
    mean_cf_remov = round(mean(bc_daily_mean_cf, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print(n = 40)

p <- ggplot(ts_data, aes(x = plot_date, y = bc_daily_mean)) +
  
  geom_line(color = "grey75", linewidth = 0.35) +
  
  geom_point(aes(color = above_threshold), size = 1.6, alpha = 0.85) +
  
  scale_color_manual(
    values = c("FALSE" = "#6baed6", "TRUE" = "#CC0000"),
    labels = c("FALSE" = "Below 15 µg/m³", "TRUE" = "Above 15 µg/m³"),
    name   = NULL
  ) +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::comma
  ) +
  
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  
  facet_wrap2(~ year_label, ncol = 2, strip.position = "top", axes = "all", remove_labels = "none") +  # x-axis on all panels
  
  labs(
    title    = "Population-Weighted Mean PM<sub>2.5</sub> During BC Fire Season (May - October)",
    subtitle = paste0("Daily Exposure Averaged Across BC Residents by Dissemination Area Population"),
    x        = "Month",
    y = "Population-Weighted Mean PM<sub>2.5</sub> (\u00b5g/m\u00b3)",
    caption  = "BlueSky Canada smoke forecast system (Noon PST observed snapshot). Statistics Canada 2021 Census."
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_markdown(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle      = element_markdown(color = "grey40", size = 16, hjust = 0.5,
                                      margin = margin(b = 8)),
    plot.caption       = element_text(color = "grey50", size = 15, hjust = 0.5,
                                      lineheight = 1.4, margin = margin(t = 10)),
    strip.text         = element_text(face = "bold", size = 15),
    legend.position    = "bottom",
    legend.text        = element_text(size = 15),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text          = element_text(color = "grey30", size = 15),
    axis.title.x = ggtext::element_markdown(color = "grey30", size = 15),
    axis.title.y = ggtext::element_markdown(color = "grey30", size = 15),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

ts_data %>% arrange(desc(bc_daily_mean)) %>% head(10)

print(p)
ggsave(
  "timeseries_fire_season_2021_2025.png",
  plot   = p,
  width  = 11,
  height = 9,
  dpi    = 300,
  bg     = "white"
)
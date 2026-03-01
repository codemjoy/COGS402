library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggh4x)  

years <- 2021:2025

ts_data <- map_dfr(years, function(yr) {
  readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds")) %>%
    group_by(date) %>%
    summarise(
      bc_daily_max  = max(pm25_mean,  na.rm = TRUE),
      bc_daily_mean = mean(pm25_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      year      = yr,
      plot_date = as.Date(paste0("2000-", format(date, "%m-%d")))
    ) %>%
    filter(month(date) >= 5 & month(date) <= 10)   # fire season only
}) %>%
  mutate(
    above_threshold = bc_daily_max >= 15,
    year_label = as.character(year)
  ) %>%
  mutate(bc_daily_max = pmax(bc_daily_max, 0.1))

p <- ggplot(ts_data, aes(x = plot_date, y = bc_daily_max)) +
  
  geom_line(color = "grey75", linewidth = 0.35) +
  
  geom_point(aes(color = above_threshold), size = 1.6, alpha = 0.85) +
  
  scale_color_manual(
    values = c("FALSE" = "#6baed6", "TRUE" = "#CC0000"),
    labels = c("FALSE" = "Below 15 µg/m³", "TRUE" = "Above 15 µg/m³"),
    name   = NULL
  ) +
  
  scale_y_log10(
    breaks = c(1, 5, 10, 50, 100, 500, 1000),
    labels = scales::comma
  ) +
  
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  
  facet_wrap2(~ year_label, ncol = 2, strip.position = "top", axes = "all", remove_labels = "none") +  # x-axis on all panels
  
  labs(
    title    = "Daily Peak PM\u2082.\u2085 Across BC During Fire Season, 2021\u20132025",
    subtitle = "Maximum observed PM\u2082.\u2085 across all BC dissemination areas on each day \u00b7 BlueSky band 1\n (9am UTC observed snapshot) Fire season: May\u2013October",
    x        = "Month",
    y        = "Peak PM\u2082.\u2085 (\u00b5g/m\u00b3, log scale)",
    caption  = "BlueSky fire weather model. Statistics Canada 2021 Census dissemination areas."
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle      = element_text(color = "grey40", size = 14, hjust = 0.5,
                                      margin = margin(b = 8)),
    plot.caption       = element_text(color = "grey50", size = 13, hjust = 0.5,
                                      lineheight = 1.4, margin = margin(t = 10)),
    strip.text         = element_text(face = "bold", size = 15),
    legend.position    = "bottom",
    legend.text        = element_text(size = 15),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text          = element_text(color = "grey30", size = 8),
    axis.title         = element_text(color = "grey30", size = 9),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

print(p)
ggsave(
  "timeseries_fire_season_2021_2025.png",
  plot   = p,
  width  = 11,
  height = 9,
  dpi    = 300,
  bg     = "white"
)
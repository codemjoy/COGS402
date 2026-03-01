library(sf)
library(ggplot2)
library(dplyr)
library(rmapshaper)
library(scales)
library(lubridate)

if (!dir.exists("figures")) dir.create("figures")

# ---- Set year here ----
year <- 2023

daily_file <- paste0("da_daily_pm25_exposure_", year, ".rds")

da_daily <- readRDS(daily_file) %>%
  filter(month(date) %in% 5:9)

da_fire_season <- da_daily %>%
  group_by(DGUID) %>%
  summarise(
    fire_season_mean_pm25 = mean(pm25_mean, na.rm = TRUE),
    fire_season_max_pm25  = max(pm25_max, na.rm = TRUE),
    days_above_15         = sum(pm25_mean > 15, na.rm = TRUE),
    days_above_37.5       = sum(pm25_mean > 37.5, na.rm = TRUE),
    days_above_50         = sum(pm25_mean > 50, na.rm = TRUE),
    n_days                = n(),
    .groups = "drop"
  )

gpkg_file <- paste0("da_annual_pm25_summary_", year, ".gpkg")

da_spatial <- st_read(gpkg_file, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(3005) %>%   
  ms_simplify(keep = 0.05, keep_shapes = TRUE) %>%
  select(DGUID) %>%
  left_join(da_fire_season, by = "DGUID")

rm(da_spatial)

names(da_spatial)

# ---- City reference points ----
cities <- data.frame(
  name = c("Vancouver", "Victoria", "Kelowna", "Prince George", "Fort St. John"),
  lon  = c(-123.12, -123.37, -119.50, -122.75, -120.85),
  lat  = c(  49.28,   48.43,   49.88,   53.92,   56.25)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3005)  

# ---- Summary stats for annotation ----
stats <- da_spatial %>%
  sf::st_drop_geometry() %>%
  summarise(
    mean_val   = mean(fire_season_mean_pm25, na.rm = TRUE),
    median_val = median(fire_season_mean_pm25, na.rm = TRUE),
    p95        = quantile(fire_season_mean_pm25, 0.95, na.rm = TRUE),
    n_da       = n()
  )

# ---- Build map ----
p <- ggplot(da_spatial) +
  geom_sf(aes(fill = fire_season_mean_pm25), color = NA, linewidth = 0) +
  geom_sf(data = cities, shape = 21, fill = "grey40", color = "grey40",
          size = 2, stroke = 0.5) +
  geom_sf_text(data = cities, aes(label = name),
               size = 4, color = "grey40", fontface = "bold",
               nudge_y = 30000,   # nudge in metres (BC Albers)
               check_overlap = TRUE) +
  scale_fill_viridis_c(
    option   = "magma",
    direction = -1, 
    limits = c(0, 75),
    oob = scales::squish,
    name     = "PM2.5\n(µg/m³)",
    na.value = "grey50",
    labels   = label_number(accuracy = 0.1)
  ) +
  labs(
    title    = paste0("Wildfire PM2.5 Exposure by Dissemination Area — BC ", year),
    subtitle = paste0("Mean PM2.5 across the 2023 Wildfire Season (May-September)| BlueSky Canada\n",
    "fire emissions model"),
    caption  = "Data: BlueSky Canada | StatCan DA boundaries | Projection: BC Albers (EPSG:3005)"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "#ffffffC9", color = NA),
    panel.background = element_rect(fill = "#ffffffC9", color = NA),
    plot.title       = element_text(color = "grey40", size = 20, face = "bold",
                                    margin = margin(b = 4, t = 8, l = 8)),
    plot.subtitle    = element_text(color = "grey70", size = 18,
                                    margin = margin(b = 6, l = 8)),
    plot.caption     = element_text(color = "grey50", size = 12,
                                    margin = margin(t = 6, b = 6, r = 8), hjust = 1),
    legend.position   = "right",
    legend.title      = element_text(color = "grey40", size = 9),
    legend.text       = element_text(color = "grey80", size = 8),
    legend.key.height = unit(3, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    plot.margin       = margin(10, 10, 10, 10)
  )

print(p)


# ---- Save ----
out_file <- paste0("figures/choropleth_pm25_DA_", year, ".png")
ggsave(out_file, plot = p, width = 12, height = 9, dpi = 300, bg = "#ffffffC9")
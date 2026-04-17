library(sf)
library(ggplot2)
library(dplyr)
library(rmapshaper)
library(scales)
library(lubridate)
library(patchwork)
library(ggtext)
library(tidyr)

YEARS <- 2021:2025

cities <- data.frame(
  name = c("Vancouver", "Victoria", "Kelowna", "Prince George", "Fort St. John"),
  lon  = c(-123.12, -123.37, -119.50, -122.75, -120.85),
  lat  = c(  49.28,   48.43,   49.88,   53.92,   56.25)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3005)

da_pop <- st_read("bc_das_population_for_raster.gpkg", quiet = TRUE) %>%
  st_drop_geometry() %>%
  select(DGUID, pop_total = pop_total_2021)

da_list <- lapply(YEARS, function(yr) {
  
  da_daily <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds")) %>%
    filter(month(date) %in% 5:10) %>%
    mutate(pm25_mean = pmin(pm25_mean, 250))   
  
  da_fire_season <- da_daily %>%
    group_by(DGUID) %>%
    summarise(
      fire_season_mean_pm25 = mean(pm25_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(da_pop, by = "DGUID") %>%
    mutate(
      pop_total = replace_na(pop_total, 0),
      pw_pm25 = fire_season_mean_pm25 * pop_total
    )
  
  st_read(paste0("da_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE) %>%
    st_make_valid() %>%
    st_transform(3005) %>%
    ms_simplify(keep = 0.05, keep_shapes = TRUE) %>%
    select(DGUID) %>%
    left_join(da_fire_season, by = "DGUID") %>%
    mutate(year = yr)
})

all_vals <- sapply(da_list, function(d) {
  d %>% st_drop_geometry() %>% pull(fire_season_mean_pm25)
}) %>% unlist()

scale_max <- quantile(all_vals, 0.99, na.rm = TRUE)  

make_panel <- function(da_sf, yr, show_legend = FALSE) {
  ggplot(da_sf) +
    geom_sf(aes(fill = fire_season_mean_pm25), color = NA, linewidth = 0) +
    geom_sf(data = cities, shape = 21, fill = "grey40", color = "grey40",
            size = 1.5, stroke = 0.5) +
    geom_sf_text(data = cities, aes(label = name),
                 size = 3.5, color = "grey30", fontface = "bold",
                 nudge_y = 30000, check_overlap = TRUE) +
    scale_fill_viridis_c(
      option    = "magma",
      direction = -1,
      limits    = c(0, scale_max),
      oob       = scales::squish,
      name      = "PM2.5\n(µg/m³)",
      na.value  = "grey50",
      labels    = label_number(accuracy = 0.1)
    ) +
    labs(title = as.character(yr)) +
    theme_void(base_size = 10) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(face = "bold", size = 14, hjust = 0.5,
                                      margin = margin(b = 4, t = 6)),
      legend.position  = if (show_legend) "right" else "none",
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.key.height = unit(3, "cm"),
      legend.key.width  = unit(0.5, "cm"),
      plot.margin      = margin(6, 6, 6, 6)
    )
}

#Build panels
panels <- lapply(seq_along(YEARS), function(i) {
  make_panel(da_list[[i]], YEARS[i], show_legend = (i == 2))  # legend on 2022
})

top_row    <- panels[[1]] + panels[[2]]
bottom_row <- panels[[3]] + panels[[4]] + panels[[5]]

p_combined <- top_row / bottom_row +
  plot_annotation(
    title    = "Wildfire PM~2.5~ Exposure by Dissemination Area — BC, 2021–2025",
    subtitle = "Mean fire-season (May–Sep) PM~2.5~ per DA; BlueSky values capped at 250 µg/m³",
    caption  = "Data: BlueSky Canada smoke forecast system | StatCan DA boundaries | Projection: BC Albers (EPSG:3005)",
    theme = theme(
      plot.title    = element_markdown(face = "bold", size = 20, margin = margin(b = 4, t = 8)),
      plot.subtitle = element_markdown(colour = "grey50", size = 16, margin = margin(b = 8)),
      plot.caption  = element_text(colour = "grey50", size = 15, hjust = 1)
    )
  )
print(p_combined)
ggsave("fig_pm25_choropleth_combined.png", p_combined,
       width = 18, height = 14, dpi = 300, bg = "white")


make_dot_panel <- function(da_sf, yr, one_dot_per = 200,
                           show_legend = FALSE, scale_max = 75) {
  
  da_valid <- da_sf %>%
    filter(
      !is.na(pop_total), pop_total > 0,
      !is.na(fire_season_mean_pm25)
    ) %>%
    mutate(n_dots = pmax(1L, as.integer(round(pop_total / one_dot_per))))
  
  n_vec  <- da_valid$n_dots
  pm_vec <- da_valid$fire_season_mean_pm25
  
  set.seed(42)
  # st_sample with a size vector returns one sfc_POINT per DA row
  # We must pass by= to ensure one geometry list per DA
  sampled <- st_sample(da_valid, size = n_vec, type = "random", by_polygon = FALSE)
  
  # sampled is a flat sfc — its length should equal sum(n_vec)
  # If it doesn't, something went wrong upstream
  message("  Expected pts: ", sum(n_vec), "  Got: ", length(sampled))
  
  if (length(sampled) != sum(n_vec)) {
    # Fall back: sample per-row to guarantee alignment
    pts_list <- lapply(seq_len(nrow(da_valid)), function(j) {
      st_sample(da_valid[j, ], size = n_vec[j], type = "random")
    })
    sampled <- do.call(c, pts_list)
    # actual_n from row-wise sampling
    actual_n <- vapply(pts_list, length, integer(1L))
  } else {
    actual_n <- n_vec   
  }
  
  # Defensive check
  stopifnot(
    all(is.finite(actual_n)),
    all(actual_n >= 0L),
    length(actual_n) == length(pm_vec)
  )
  
  dots <- st_sf(
    fire_season_mean_pm25 = rep(pm_vec, actual_n),
    geometry              = sampled,
    crs                   = 3005
  )
  
  ggplot() +
    geom_sf(data = da_sf,  fill = "grey92", color = NA) +
    geom_sf(data = dots,
            aes(color = fire_season_mean_pm25),
            size = 0.15, alpha = 0.6, shape = 16) +
    geom_sf(data = cities,
            shape = 21, fill = "white", color = "grey30",
            size = 1.5, stroke = 0.5) +
    geom_sf_text(data = cities, aes(label = name),
                 size = 9.0, color = "grey20", fontface = "bold",
                 nudge_y = 30000, check_overlap = TRUE) +
    scale_color_viridis_c(
      option    = "magma",
      direction = -1,
      limits    = c(0, scale_max),
      oob       = scales::squish,
      name      = "PM\u2082.\u2085\n(\u03bcg/m\u00b3)",
      labels    = label_number(accuracy = 0.1)
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      plot.title        = element_markdown(face = "bold", size = 20, hjust = 0.5,
                                       margin = margin(b = 4, t = 6)),
      plot.caption      = element_markdown(colour = "grey50", size = 13, hjust = 0.5,
                                       margin = margin(t = 4)),
      legend.position   = if (show_legend) "right" else "none",
      legend.title      = element_text(size = 12),
      legend.text       = element_text(size = 12),
      legend.key.height = unit(3, "cm"),
      legend.key.width  = unit(0.5, "cm"),
      plot.margin       = margin(6, 6, 6, 6)
    )
}

scale_max <- quantile(all_vals, 0.99, na.rm = TRUE)
message("Colour scale max (99th pctile): ", round(scale_max, 1), " µg/m³")

# ---- Build panels ----
dot_panels <- lapply(seq_along(YEARS), function(i) {
  message("  Panel ", i, ": ", YEARS[i])
  make_dot_panel(
    da_sf       = da_list[[i]],
    yr          = YEARS[i],
    one_dot_per = 200,
    show_legend = (i == 2),
    scale_max   = scale_max
  )
})

# ---- Assemble with patchwork ----
dot_top    <- dot_panels[[1]] + dot_panels[[2]]
dot_bottom <- dot_panels[[3]] + dot_panels[[4]] + dot_panels[[5]]

p_dot_combined <- dot_top / dot_bottom +
  plot_annotation(
    title    = "Population Exposure to Wildfire PM~2.5~ — BC, 2021\u20132025",
    subtitle = paste0("Each dot = 200 people; colour = DA mean fire-season (May\u2013Oct) ",
                      "PM~2.5~ capped at 250 \u00b5g/m\u00b3"),
    caption  = paste0("Data: BlueSky Canada smoke forecast system | ",
                      "StatCan DA boundaries | Projection: BC Albers (EPSG:3005)"),
    theme = theme(
      plot.title    = element_markdown(face = "bold", size = 21,
                                   margin = margin(b = 4, t = 8)),
      plot.subtitle = element_markdown(colour = "grey50", size = 16,
                                   margin = margin(b = 8)),
      plot.caption  = element_text(colour = "grey50", size = 14, hjust = 1)
    )
  )


print(p_dot_combined)
ggsave("fig_pm25_dot_density_combined.png", p_dot_combined,
       width = 18, height = 14, dpi = 300, bg = "white")

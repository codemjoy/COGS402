library(sf)
library(tidyverse)
library(ggtext)
library(scales)
library(patchwork)
library(bcmaps)

bc_boundary <- bc_bound() %>% st_transform(3005)

fire_centres <- st_read("DRPMFFRCNT_polygon.shp", quiet = TRUE) %>%
  st_transform(3005)

fire_centre_region_map <- tribble(
  ~FIRE_CENTRE_NAME,           ~region,
  "Cariboo Fire Centre",       "Interior",
  "Kamloops Fire Centre",      "Interior",
  "Southeast Fire Centre",     "Interior",
  "Coastal Fire Centre",       "Coast",
  "Prince George Fire Centre", "North",
  "Northwest Fire Centre",     "North"
)

nbac_files <- list(
  "2021" = "NBAC_2021_20250506.shp",
  "2022" = "nbac_2022_20250506.shp",
  "2023" = "nbac_2023_20250506.shp",
  "2024" = "nbac_2024_20250506.shp"
)

nbac_all <- nbac_files %>%
  imap(~ {
    st_read(.x, quiet = TRUE) %>%
      st_transform(3005) %>%
      st_filter(bc_boundary, .predicate = st_intersects) %>%
      mutate(year = as.integer(.y))
  }) %>%
  bind_rows()

burned_province <- nbac_all %>%
  group_by(year) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  mutate(total_ha = as.numeric(st_area(geometry)) / 10000) %>%
  st_drop_geometry() %>%
  select(year, total_ha)

message("Province-wide burned area (NBAC):")
print(burned_province)

nbac_centroids <- nbac_all %>%
  st_centroid() %>%
  st_join(
    fire_centres %>% select(FIRE_CENTRE_NAME = MFFRCNTRNM),
    join = st_within
  )

# Assign unmatched polygons to nearest fire centre
unmatched_idx <- which(is.na(nbac_centroids$FIRE_CENTRE_NAME))
if (length(unmatched_idx) > 0) {
  nearest_idx <- st_nearest_feature(nbac_centroids[unmatched_idx, ], fire_centres)
  nbac_centroids$FIRE_CENTRE_NAME[unmatched_idx] <-
    fire_centres$MFFRCNTRNM[nearest_idx]
}

message("NBAC polygons matched to a fire centre: ",
        sum(!is.na(nbac_centroids$FIRE_CENTRE_NAME)),
        " / ", nrow(nbac_centroids))

fc_assignment <- nbac_centroids %>%
  st_drop_geometry() %>%
  select(year, FIRE_CENTRE_NAME) %>%
  bind_cols(st_drop_geometry(nbac_all) %>% select(-year))

burned_area <- nbac_all %>%
  mutate(FIRE_CENTRE_NAME = fc_assignment$FIRE_CENTRE_NAME) %>%
  left_join(fire_centre_region_map, by = "FIRE_CENTRE_NAME") %>%
  group_by(year, FIRE_CENTRE_NAME, region) %>%
  summarise(
    geometry   = st_union(geometry),
    n_polygons = n(),
    .groups    = "drop"
  ) %>%
  st_as_sf() %>%
  mutate(total_ha = as.numeric(st_area(geometry)) / 10000) %>%
  st_drop_geometry() %>%
  select(year, FIRE_CENTRE_NAME, region, total_ha, n_polygons)

print(burned_area %>% arrange(year, FIRE_CENTRE_NAME), n = Inf)
write_csv(burned_area %>% arrange(year, FIRE_CENTRE_NAME),
          "annual_burned_area_by_fire_centre.csv")

burned_area_regional <- burned_area %>%
  group_by(year, region) %>%
  summarise(
    total_ha   = sum(total_ha),
    n_polygons = sum(n_polygons),
    .groups    = "drop"
  )

message("\nBurned area by region and year:")
print(burned_area_regional %>% arrange(region, year), n = Inf)
write_csv(burned_area_regional, "burned_area_by_region.csv")

# ── IER dementia cases ────────────────────────────────────────────────────────
mc_annual_summary <- readRDS("mc_annual_summary.rds")

ier_cases <- mc_annual_summary %>%
  filter(Scenario == "IER GBD 2019, CRA") %>%
  transmute(
    year         = as.integer(Year),
    cases_median = Median,
    ci_lo        = CI_lo,
    ci_hi        = CI_hi
  ) %>%
  filter(year %in% burned_province$year)  # align to NBAC years only

# ── Figure 1: Burned area by fire centre ─────────────────────────────────────
p_fc <- ggplot(burned_area,
               aes(x = year, y = total_ha / 1000, fill = FIRE_CENTRE_NAME)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Spectral", name = "Fire Centre") +
  scale_x_continuous(breaks = 2021:2024) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "BC Wildfire Burned Area by Fire Centre, 2021\u20132024",
    subtitle = "Canadian National Burned Area Composite (NBAC), clipped to BC",
    x        = "Fire season year",
    y        = "Burned area (thousand ha)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )

print(p_fc)
ggsave("burned_area_by_fire_centre.png", p_fc,
       width = 9, height = 5, dpi = 300, bg = "white")

# ── Figure 2: Burned area by region ──────────────────────────────────────────
p_region <- ggplot(burned_area_regional,
                   aes(x = year, y = total_ha / 1000, fill = region)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Set1", name = "Region") +
  scale_x_continuous(breaks = 2021:2024) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "BC Wildfire Burned Area by Analysis Region, 2021\u20132024",
    subtitle = "Canadian National Burned Area Composite (NBAC), clipped to BC",
    x        = "Fire season year",
    y        = "Burned area (thousand ha)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_region)
ggsave("burned_area_by_region.png", p_region,
       width = 9, height = 5, dpi = 300, bg = "white")

# ── Figure 3: Two-panel comparison ───────────────────────────────────────────
p_burned <- ggplot(burned_province, aes(x = year, y = total_ha)) +
  geom_col(fill = "#998731", alpha = 0.6, width = 0.6) +
  scale_x_continuous(breaks = 2021:2024) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Provincial burned area, 2021\u20132024",
    x     = "Fire season year",
    y     = "Burned area (ha)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_cases <- ggplot(ier_cases, aes(x = year)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.15, colour = "#998731", linewidth = 0.7) +
  geom_point(aes(y = cases_median), colour = "#998731", size = 3) +
  geom_line(aes(y = cases_median), colour = "#998731", linewidth = 0.8) +
  scale_x_continuous(breaks = 2021:2024) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "IER-attributable dementia cases, 2021\u20132024",
    x     = "Fire season year",
    y     = "Attributable cases (95% CI)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_compare <- p_burned / p_cases +
  plot_annotation(
    title    = "Wildfire burned area and PM~2.5~-attributable Alzheimer\u2019s burden, BC 2021\u20132024",
    subtitle = "Top: NBAC provincial burned area \u00b7 Bottom: IER attributable cases (point = median, bars = 95% CI)",
    theme    = theme(plot.title = element_markdown(size = 13, face = "bold"))
  )

print(p_compare)
ggsave("burned_area_vs_cases_comparison.png", p_compare,
       width = 7, height = 7, dpi = 300, bg = "white")

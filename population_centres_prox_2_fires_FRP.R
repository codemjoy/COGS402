library(sf)
library(tidyverse)
library(ggtext)
library(scales)
library(patchwork)
library(bcmaps)

YEARS         <- 2021:2024   
FIRE_MONTHS   <- 5:10
CRS_BC        <- 3005        

bc_boundary <- bc_bound() %>% st_transform(CRS_BC)

cf_regional <- read_csv("cf_annual_regional.csv") %>%
  filter(fire_year %in% YEARS) %>%
  mutate(fire_year = as.integer(fire_year))

# Fire centre boundaries
fire_centres <- st_read("DRPMFFRCNT_polygon.shp", quiet = TRUE) %>%
  st_transform(CRS_BC)

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
      st_transform(CRS_BC) %>%
      st_filter(bc_boundary, .predicate = st_intersects) %>%
      mutate(year = as.integer(.y))
  }) %>%
  bind_rows()

# VIIRS fire hotspots 
viirs_raw <- st_read("fire_archive_J1V-C2_735471.shp", quiet = TRUE) %>%
  st_transform(CRS_BC) %>%
  filter(
    month(ACQ_DATE) %in% FIRE_MONTHS,
    year(ACQ_DATE)  %in% YEARS,
    CONFIDENCE %in% c("h", "n"),   # high and nominal confidence only
    TYPE == 0,                      # vegetation fires only 
    FRP >= 10                       # try to remove industrial activity?
  ) %>%
  st_filter(bc_boundary, .predicate = st_intersects) %>%
  mutate(year = year(ACQ_DATE))

message("VIIRS hotspots in BC fire season: ", nrow(viirs_raw))

# Major BC cities
cities <- tribble(
  ~city_label,      ~pop_total, ~longitude,  ~latitude,
  "Vancouver",        662248,   -123.1207,   49.2827,
  "Surrey",           568322,   -122.8490,   49.1913,
  "Burnaby",          249125,   -122.9774,   49.2488,
  "Richmond",         209937,   -123.1376,   49.1666,
  "Kelowna",          144576,   -119.4960,   49.8880,
  "Abbotsford",       153524,   -122.3017,   49.0504,
  "Coquitlam",        148625,   -122.7932,   49.2838,
  "Kamloops",         97902,    -120.3273,   50.6745,
  "Nanaimo",          99863,    -123.9401,   49.1659,
  "Prince George",    74003,    -122.7497,   53.9171,
  "Chilliwack",       93203,    -121.9543,   49.1579,
  "Vernon",           65282,    -119.2720,   50.2674,
  "Victoria",         92141,    -123.3656,   48.4284,
  "Penticton",        48583,    -119.5885,   49.4991,
  "Langley",          132603,   -122.6604,   49.1044,
  "Terrace",          16283,    -128.6032,   54.5164,
  "Prince Rupert",    12220,    -130.3208,   54.3150,
  "Williams Lake",    12859,    -122.1417,   52.1418,
  "Fort St. John",    23817,    -120.8487,   56.2527,
  "Dawson Creek",     13474,    -120.2327,   55.7596,
  "Quesnel",          11916,    -122.4924,   52.9784,
  "Smithers",         5401,     -127.1669,   54.7806,
  "Merritt",          7139,     -120.7867,   50.1133,
  "Cranbrook",        26083,    -115.7694,   49.5097,
  "Castlegar",        8540,     -117.6561,   49.3253,
  "Campbell River",   37571,    -125.2733,   50.0163,
  "Courtenay",        28420,    -124.9897,   49.6878
) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(CRS_BC)

message("Cities identified: ", nrow(cities))
print(cities %>% st_drop_geometry() %>%
        select(city_label, pop_total,) %>%
        arrange(desc(pop_total)))

# days with active fires proximate to cities
hotspot_daily <- map_dfr(YEARS, function(yr) {
  viirs_yr <- viirs_raw %>% filter(year == yr)
  if (nrow(viirs_yr) == 0) return(NULL)
  message("  Year: ", yr)
  
  map_dfr(seq_len(nrow(cities)), function(i) {
    city  <- cities[i, ]
    dists <- st_distance(viirs_yr, city) %>% as.numeric()
    
    # All unique dates with a hotspot within each band
    dates_10  <- sort(unique(viirs_yr$ACQ_DATE[dists <= 10000]))
    dates_50  <- sort(unique(viirs_yr$ACQ_DATE[dists <= 50000]))
    dates_100 <- sort(unique(viirs_yr$ACQ_DATE[dists <= 100000]))
    
    # Helper: longest consecutive run of days
    longest_streak <- function(dates) {
      if (length(dates) == 0) return(0L)
      gaps <- c(1, diff(as.integer(dates)))
      runs <- rle(gaps == 1)
      # streak length = run length of TRUE gaps + 1 for first day
      if (!any(runs$values)) return(1L)
      max(runs$lengths[runs$values]) + 1L
    }
    
    tibble(
      year                = yr,
      city_label          = cities$city_label[i],
      pop_total           = cities$pop_total[i],
      # Total days
      total_days_10km     = length(dates_10),
      total_days_50km     = length(dates_50),
      total_days_100km    = length(dates_100),
      # Longest consecutive streak
      max_streak_10km     = longest_streak(dates_10),
      max_streak_50km     = longest_streak(dates_50),
      max_streak_100km    = longest_streak(dates_100),
      # First and last date of exposure within 100 km
      first_date_100km    = if (length(dates_100) > 0) min(dates_100) else NA_Date_,
      last_date_100km     = if (length(dates_100) > 0) max(dates_100) else NA_Date_,
      # Peak FRP within 100 km
      max_frp_100km       = if (any(dists <= 100000))
        max(viirs_yr$FRP[dists <= 100000], na.rm = TRUE) else NA_real_
    )
  })
})

write_csv(hotspot_daily, "city_hotspot_proximity_temporal.csv")

hotspot_summary <- hotspot_daily %>%
  mutate(
    exposure_window_100km = case_when(
      is.na(first_date_100km) ~ "None",
      TRUE ~ paste0(format(first_date_100km, "%b %d"),
                    " to ",
                    format(last_date_100km, "%b %d"))
    )
  ) %>%
  select(year, city_label, pop_total,
         total_days_10km, max_streak_10km,
         total_days_50km, max_streak_50km,
         total_days_100km, max_streak_100km,
         exposure_window_100km, max_frp_100km) %>%
  arrange(year, desc(total_days_100km))

print(hotspot_summary, n = Inf)
write.csv(hotspot_summary, "city_fire_exposure_summary.csv")

# age comp of populations near fires
da_sf_bc <- st_read("bc_das.shp", quiet = TRUE) %>%
  st_transform(CRS_BC) %>%
  filter(PRUID == "59")

da_census <- st_read("bc_das_population_for_raster.gpkg", quiet = TRUE)

da_census_df <- if (inherits(da_census, "sf")) {
  as.data.frame(da_census)
} else {
  da_census
}

da_joined <- da_sf_bc %>%
  left_join(
    da_census_df %>%
      select(DGUID, pop_total_2021, pop_total_65_74,
             pop_total_75_84, pop_total_85plus),
    by = "DGUID"
  )

age_near_fire <- map_dfr(YEARS, function(yr) {
  nbac_yr <- nbac_all %>% filter(year == yr)
  
  da_cents <- da_joined %>%
    st_centroid()
  
  # Calculate distance from each DA centroid to nearest individual fire polygon
  dist_matrix <- st_distance(da_cents, nbac_yr)  
  min_dist    <- apply(dist_matrix, 1, min)       
  
  da_cents <- da_cents %>%
    mutate(dist_to_nearest_fire = as.numeric(min_dist))
  
  map_dfr(c(10, 50, 100), function(km) {
    da_cents %>%
      st_drop_geometry() %>%
      filter(dist_to_nearest_fire <= km * 1000) %>%
      summarise(
        year       = yr,
        dist_band  = paste0("Within ", km, " km"),
        n_das      = n(),
        pop_total  = sum(pop_total_2021,  na.rm = TRUE),
        pop_65_74  = sum(pop_total_65_74, na.rm = TRUE),
        pop_75_84  = sum(pop_total_75_84, na.rm = TRUE),
        pop_85plus = sum(pop_total_85plus, na.rm = TRUE),
        pop_65plus = pop_65_74 + pop_75_84 + pop_85plus,
        pct_65plus = round(pop_65plus / pop_total * 100, 1),
        pct_85plus = round(pop_85plus / pop_total * 100, 1)
      )
  })
})

message("\nAge composition near fires:")
age_near_fire |> as.data.frame() |> View()
write_csv(age_near_fire, "age_composition_near_fires.csv")


# Fig1: Provincial NBAC perimeters + VIIRS FRP
viirs_plot <- viirs_raw %>%
  mutate(
    frp_class = cut(FRP,
                    breaks = c(0, 50, 150, 300, 500, Inf),
                    labels = c("<50", "50\u2013150", "150\u2013300",
                               "300\u2013500", ">500"),
                    right  = FALSE)
  )

frp_colours <- c(
  "<50"          = "#FFEDA0",
  "50\u2013150"  = "#FEB24C",
  "150\u2013300" = "#FC4E2A",
  "300\u2013500" = "#BD0026",
  ">500"         = "#67000D"
)

cities_plot <- cities %>% 
  filter(city_label %in% c(
    "Vancouver", "Kamloops", "Kelowna", "Penticton", "Prince George", "Fort St. John", "Prince Rupert", 
    "Victoria", "Cranbrook", "Fort Nelson", "100 Mile House"
  ))

nbac_all <- st_intersection(nbac_all, bc_boundary)

p_maps <- ggplot() +
  geom_sf(data = bc_boundary,
          fill = "grey95", colour = "grey70", linewidth = 0.3) +
  geom_sf(data = viirs_plot,
          aes(colour = frp_class), size = 0.25, alpha = 0.6) +
  geom_sf(data = cities_plot,
          shape = 21, fill = "white", colour = "black",
          size = 1.6, stroke = 0.5) +
  geom_sf_text(
    data = cities_plot,
    aes(label = city_label),
    size = 2.0,
    colour = "black",
    nudge_y = 15000
  ) +
  facet_wrap(~ year, ncol = 2) +
  scale_colour_manual(
    values = frp_colours,
    name   = "FRP (MW)",
    guide  = guide_legend(override.aes = list(size = 2.5, alpha = 1))
  ) +
  coord_sf(crs = CRS_BC) +
  labs(
    title    = "BC Wildfire Activity and Fire Radiative Power, 2021\u20132024",
    subtitle = paste0("Points = VIIRS active fire hotspots coloured by FRP (MW)\n",
                      "Fire season May\u2013Oct, high/nominal confidence, vegetation fires only"),
    caption  = "Sources: Natural Resources Canada NBAC; NASA FIRMS VIIRS NOAA-20 Collection 2"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = "grey40", size = 8.5),
    plot.caption     = element_text(colour = "grey50", size = 7.5),
    strip.text       = element_text(face = "bold", size = 11),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.margin      = margin(10, 10, 10, 10)
  )


print(p_maps)
ggsave("fig_fire_maps_frp.png", p_maps,
       width = 12, height = 10, dpi = 300, bg = "white")

# Fig 2: Hotspot proximity heatmap 
hotspot_days <- hotspot_daily %>%
  filter(city_label %in% cities_plot$city_label) %>%
  select(year, city_label, pop_total,
         days_hotspot_10km  = total_days_10km,
         days_hotspot_50km  = total_days_50km,
         days_hotspot_100km = total_days_100km)

hotspot_heat <- hotspot_days %>%
  mutate(city_label = fct_reorder(city_label, pop_total, .desc = TRUE)) %>%
  pivot_longer(cols = starts_with("days_hotspot"),
               names_to = "band", values_to = "days") %>%
  mutate(band = recode(band,
                       "days_hotspot_10km"  = "Within 10 km",
                       "days_hotspot_50km"  = "Within 50 km",
                       "days_hotspot_100km" = "Within 100 km"),
         band = factor(band, levels = c("Within 10 km",
                                        "Within 50 km",
                                        "Within 100 km")))

p_heat <- ggplot(hotspot_heat,
                 aes(x = factor(year), y = city_label, fill = days)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(days == 0, "", days)),
            size = 2.6, colour = "grey15") +
  facet_wrap(~ band, ncol = 3) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    name      = "Days with\nactive fire",
    limits    = c(0, NA)
  ) +
  labs(
    title    = "Days with Active Fire Within Distance Bands of BC Cities, 2021\u20132024",
    subtitle = paste0("Count of fire-season days (May\u2013Oct) with \u22651 VIIRS hotspot ",
                      "within each distance band of city centre"),
    x        = "Year",
    y        = NULL,
    caption  = "Source: NASA FIRMS VIIRS NOAA-20 Collection 2; Statistics Canada 2021 Census"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(colour = "grey40", size = 9),
    plot.caption    = element_text(colour = "grey50", size = 8),
    axis.text.y     = element_text(size = 8.5),
    strip.text      = element_text(face = "bold"),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

print(p_heat)
ggsave("fig_hotspot_proximity_heatmap.png", p_heat,
       width = 14, height = 7, dpi = 300, bg = "white")

# Fig 3: Age comp near fires
age_long <- age_near_fire %>%
  pivot_longer(
    cols      = c(pop_65_74, pop_75_84, pop_85plus),
    names_to  = "age_group",
    values_to = "population"
  ) %>%
  mutate(
    age_group = recode(age_group,
                       "pop_65_74"  = "65\u201374",
                       "pop_75_84"  = "75\u201384",
                       "pop_85plus" = "85+"),
    age_group = factor(age_group,
                       levels = c("85+", "75\u201384", "65\u201374")),
    dist_band = factor(dist_band,
                       levels = c("Within 10 km",
                                  "Within 50 km",
                                  "Within 100 km"))
  )

p_age <- ggplot(age_long,
                aes(x = factor(year), y = population / 1000,
                    fill = age_group)) +
  geom_col(position = "stack") +
  facet_wrap(~ dist_band, ncol = 3) +
  scale_fill_manual(
    values = c("65\u201374" = "#9ECAE1",
               "75\u201384" = "#3182BD",
               "85+"        = "#08306B"),
    name = "Age group"
  ) +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Older Adult Population Near BC Wildfires by Age Group, 2021\u20132024",
    subtitle = "Population aged 65+ residing within each distance band of NBAC fire perimeters",
    x        = "Year",
    y        = "Population (thousands)",
    caption  = "Source: Statistics Canada 2021 Census; Natural Resources Canada NBAC"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(colour = "grey40", size = 9),
    plot.caption       = element_text(colour = "grey50", size = 8),
    strip.text         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )
print(p_age)
ggsave("fig_age_composition_near_fires.png", p_age,
       width = 12, height = 5, dpi = 300, bg = "white")

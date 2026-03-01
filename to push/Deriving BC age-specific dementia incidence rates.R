library(sf)
library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)

# PARAMETERS 
beta_central <- log(1.04) / 4.8
beta_low     <- log(1.03) / 4.8
beta_high    <- log(1.05) / 4.8
counterfactual <- 15

primary_rates <- list(
  rate_65_74  =  6.10 / 1000,
  rate_75_84  = 34.50 / 1000,
  rate_85plus = 44.68 / 1000
)

years <- 2021:2026
# ============================================================
# CRF FUNCTION 
calculate_crf <- function(data, beta, counterfactual, rates) {
  data %>%
    st_drop_geometry() %>%
    mutate(
      delta_pm25      = pmax(fireseason_mean_pm25 - counterfactual, 0),
      RR              = exp(beta * delta_pm25),
      PAF             = (RR - 1) / RR,
      expected_65_74  = pop_total_65_74  * rates$rate_65_74,
      expected_75_84  = pop_total_75_84  * rates$rate_75_84,
      expected_85plus = pop_total_85plus * rates$rate_85plus,
      expected_total  = expected_65_74 + expected_75_84 + expected_85plus,
      attrib_65_74    = expected_65_74  * PAF,
      attrib_75_84    = expected_75_84  * PAF,
      attrib_85plus   = expected_85plus * PAF,
      attrib_total    = expected_total  * PAF
    )
}
# ============================================================
# MULTI-YEAR LOOP
all_years_results <- map_dfr(years, function(yr) {
  
  cat("Processing", yr, "...\n")
  
  # --- Load daily file and calculate fire-season exposure ---
  daily <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds"))
  
  fireseason <- daily %>%
    filter(
      date >= as.Date(paste0(yr, "-05-01")),
      date <= as.Date(paste0(yr, "-09-30"))
    ) %>%
    group_by(DGUID) %>%
    summarise(
      fireseason_mean_pm25       = mean(pm25_mean,    na.rm = TRUE),
      fireseason_days_above_15   = sum(above_15,      na.rm = TRUE),
      fireseason_days_above_37.5 = sum(above_37.5,    na.rm = TRUE),
      fireseason_days_above_50   = sum(above_50,      na.rm = TRUE),
      n_fireseason_days          = n(),
      .groups = "drop"
    )
  
  # --- Load annual gpkg (carries population data) ---
  annual <- st_read(
    paste0("da_annual_pm25_summary_", yr, ".gpkg"),
    quiet = TRUE
  )
  
  # --- Join ---
  da_joined <- annual %>%
    left_join(fireseason, by = "DGUID")
  
  # --- Run CRF ---
  calculate_crf(da_joined, beta_central, counterfactual, primary_rates) %>%
    mutate(year = yr)
})

cat("Done. Total rows:", nrow(all_years_results), "\n")

# ============================================================
# EXPOSURE SUMMARY STATISTICS BY YEAR
exposure_summary <- map_dfr(years, function(yr) {
  
  # Get the fire-season exposure for this year
  daily <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds"))
  
  fireseason <- daily %>%
    filter(
      date >= as.Date(paste0(yr, "-05-01")),
      date <= as.Date(paste0(yr, "-09-30"))
    ) %>%
    group_by(DGUID) %>%
    summarise(
      fireseason_mean_pm25 = mean(pm25_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  annual <- st_read(
    paste0("DA_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE
  ) %>%
    st_drop_geometry() %>%
    select(DGUID, DAUID, pop_total_65plus)
  
  joined <- annual %>%
    left_join(fireseason, by = "DGUID")
  
  # Calculate metrics
  total_das        <- nrow(joined)
  exposed_das      <- sum(joined$fireseason_mean_pm25 > 15, na.rm = TRUE)
  pct_das_exposed  <- (exposed_das / total_das) * 100
  
  total_pop_65plus   <- sum(joined$pop_total_65plus, na.rm = TRUE)
  exposed_pop_65plus <- sum(
    joined$pop_total_65plus[joined$fireseason_mean_pm25 > 15], na.rm = TRUE
  )
  pct_pop_exposed  <- (exposed_pop_65plus / total_pop_65plus) * 100
  
  tibble(
    year               = yr,
    total_das          = total_das,
    exposed_das        = exposed_das,
    pct_das_exposed    = round(pct_das_exposed, 1),
    total_pop_65plus   = total_pop_65plus,
    exposed_pop_65plus = exposed_pop_65plus,
    pct_pop_exposed    = round(pct_pop_exposed, 1)
  )
})

cat("\n=== EXPOSURE SUMMARY BY YEAR ===\n")
print(exposure_summary)

# Summarise across all years (mean and range)
cat("\n=== EXPOSURE SUMMARY ACROSS 2021-2026 ===\n")
cat("Mean % DAs exceeding WHO guideline:   ",
    round(mean(exposure_summary$pct_das_exposed), 1), "%",
    "(range:", min(exposure_summary$pct_das_exposed), "–",
    max(exposure_summary$pct_das_exposed), "%)\n")

cat("Mean exposed pop 65+ per year:        ",
    round(mean(exposure_summary$exposed_pop_65plus), 0), "people",
    "(range:", min(exposure_summary$exposed_pop_65plus), "–",
    max(exposure_summary$exposed_pop_65plus), ")\n")

cat("Mean % of BC pop 65+ exposed per year:",
    round(mean(exposure_summary$pct_pop_exposed), 1), "%",
    "(range:", min(exposure_summary$pct_pop_exposed), "–",
    max(exposure_summary$pct_pop_exposed), "%)\n")

# Active fire years only
cat("\n=== EXPOSURE SUMMARY — ACTIVE FIRE YEARS (2021-2023) ONLY ===\n")
exposure_summary %>%
  summarise(
    mean_pct_das_exposed    = round(mean(pct_das_exposed), 1),
    range_pct_das           = paste0(min(pct_das_exposed), "–", max(pct_das_exposed), "%"),
    mean_exposed_pop        = round(mean(exposed_pop_65plus), 0),
    range_exposed_pop       = paste0(min(exposed_pop_65plus), "–", max(exposed_pop_65plus)),
    mean_pct_pop_exposed    = round(mean(pct_pop_exposed), 1),
    range_pct_pop           = paste0(min(pct_pop_exposed), "–", max(pct_pop_exposed), "%")
  ) %>%
  print()

# ============================================================
# YEAR-BY-YEAR SUMMARY TABLE
summary_by_year <- all_years_results %>%
  group_by(year) %>%
  summarise(
    pop_65plus         = sum(pop_total_65plus,  na.rm = TRUE),
    expected_cases     = sum(expected_total,    na.rm = TRUE),
    attributable_cases = sum(attrib_total,      na.rm = TRUE),
    pct_attributable   = (attributable_cases / expected_cases) * 100,
    mean_pm25_wtd      = weighted.mean(
      fireseason_mean_pm25, pop_total_65plus, na.rm = TRUE),
    n_das_exposed      = sum(delta_pm25 > 0, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_by_year)

# Cumulative total
cat("\n=== CUMULATIVE 2021-2026 ===\n")
cat("Total attributable cases (central):",
    round(sum(summary_by_year$attributable_cases), 1), "\n")
# ============================================================
# REGIONAL BREAKDOWN BY YEAR
regional_by_year <- all_years_results %>%
  filter(delta_pm25 > 0) %>%
  mutate(
    cd_code = substr(DAUID, 1, 4),
    region = case_when(
      cd_code == "5901" ~ "Capital / Greater Victoria",
      cd_code == "5903" ~ "Cowichan Valley",
      cd_code == "5905" ~ "Nanaimo",
      cd_code == "5907" ~ "Alberni-Clayoquot / Central Island",
      cd_code == "5909" ~ "Comox Valley",
      cd_code == "5911" ~ "Strathcona",
      cd_code == "5915" ~ "Metro Vancouver",
      cd_code == "5919" ~ "Squamish-Lillooet",
      cd_code == "5921" ~ "Fraser Valley",
      cd_code == "5933" ~ "Thompson-Nicola / Kamloops",
      cd_code == "5935" ~ "Central Okanagan / Kootenay",
      cd_code == "5937" ~ "Okanagan-Similkameen",
      cd_code == "5939" ~ "Peace River / Northeast BC",
      cd_code == "5941" ~ "Cariboo",
      cd_code == "5943" ~ "Skeena",
      cd_code == "5945" ~ "Kitimat-Stikine",
      cd_code == "5947" ~ "Bulkley-Nechako",
      cd_code == "5949" ~ "Fraser-Fort George / Prince George",
      cd_code == "5953" ~ "Central Okanagan",
      cd_code == "5955" ~ "North Okanagan",
      cd_code == "5957" ~ "Columbia-Shuswap",
      cd_code == "5959" ~ "East Kootenay",
      TRUE ~ paste0("CD ", cd_code)
    )
  ) %>%
  group_by(year, region) %>%
  summarise(
    n_das              = n(),
    pop_65plus         = sum(pop_total_65plus,   na.rm = TRUE),
    mean_pm25          = weighted.mean(
      fireseason_mean_pm25, pop_total_65plus, na.rm = TRUE),
    attributable_cases = sum(attrib_total,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year, desc(attributable_cases))

# Top 5 regions per year
top5_by_year <- regional_by_year %>%
  group_by(year) %>%
  slice_max(attributable_cases, n = 5) %>%
  ungroup()

print(top5_by_year, n = 30)
# ============================================================
# CUMULATIVE REGIONAL BURDEN (all years combined)
regional_cumulative <- regional_by_year %>%
  group_by(region) %>%
  summarise(
    total_attributable = sum(attributable_cases, na.rm = TRUE),
    mean_annual        = mean(attributable_cases, na.rm = TRUE),
    n_years_exposed    = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(total_attributable))

print(regional_cumulative, n = 20)
# ============================================================
# SENSITIVITY: run low and high beta across all years
sensitivity_low <- map_dfr(years, function(yr) {
  daily <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds"))
  fireseason <- daily %>%
    filter(date >= as.Date(paste0(yr, "-05-01")),
           date <= as.Date(paste0(yr, "-09-30"))) %>%
    group_by(DGUID) %>%
    summarise(fireseason_mean_pm25 = mean(pm25_mean, na.rm = TRUE),
              .groups = "drop")
  annual <- st_read(paste0("DA_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE)
  da_joined <- annual %>% left_join(fireseason, by = "DGUID")
  calculate_crf(da_joined, beta_low, counterfactual, primary_rates) %>%
    mutate(year = yr)
})

sensitivity_high <- map_dfr(years, function(yr) {
  daily <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds"))
  fireseason <- daily %>%
    filter(date >= as.Date(paste0(yr, "-05-01")),
           date <= as.Date(paste0(yr, "-09-30"))) %>%
    group_by(DGUID) %>%
    summarise(fireseason_mean_pm25 = mean(pm25_mean, na.rm = TRUE),
              .groups = "drop")
  annual <- st_read(paste0("DA_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE)
  da_joined <- annual %>% left_join(fireseason, by = "DGUID")
  calculate_crf(da_joined, beta_high, counterfactual, primary_rates) %>%
    mutate(year = yr)
})

# Uncertainty summary by year
uncertainty_by_year <- bind_rows(
  all_years_results %>% mutate(scenario = "central"),
  sensitivity_low   %>% mutate(scenario = "low"),
  sensitivity_high  %>% mutate(scenario = "high")
) %>%
  group_by(year, scenario) %>%
  summarise(
    attributable_cases = sum(attrib_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = scenario, values_from = attributable_cases) %>%
  mutate(across(c(low, central, high), ~round(., 1)))

print(uncertainty_by_year)

cat("\nCumulative low: ",  round(sum(uncertainty_by_year$low),     1), "\n")
cat("Cumulative central:", round(sum(uncertainty_by_year$central), 1), "\n")
cat("Cumulative high:",    round(sum(uncertainty_by_year$high),    1), "\n")
# ============================================================
# SENSITIVITY ANALYSIS: BC Ministry of Health LHA-specific rates

lha_boundaries <- st_read("BCHA_LOCAL_HEALTH_AREA_SP.geojson", quiet = TRUE) %>%
  st_transform(3005) %>%
  mutate(lha_code = as.character(LOCAL_HLTH_AREA_CODE))

moh_rates_combined <- read_csv("BC Community Health Data Download- age standardized dementia.csv") %>%
  rename(
    lha_code      = `Jurisdiction Code`,
    fiscal_year   = `Year(s)`,
    sex           = Sex,
    rate_per_1000 = `Indicator Value`
  ) %>%
  filter(sex == "Total", !is.na(lha_code), !is.na(rate_per_1000)) %>%
  filter(fiscal_year %in% c(
    "FY 2016.2017", "FY 2017.2018", "FY 2018.2019",
    "FY 2019.2020", "FY 2020.2021"
  )) %>%
  group_by(lha_code) %>%
  summarise(rate_per_1000 = mean(rate_per_1000, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    lha_code        = as.character(lha_code),
    rate_per_person = rate_per_1000 / 1000
  )

# Spatial join — DA centroids to LHA
census_das_albers <- census_das %>% st_transform(3005)

da_centroids <- census_das_albers %>%
  st_centroid() %>%
  select(DGUID, DAUID)

da_lha_lookup <- st_join(
  da_centroids,
  lha_boundaries %>% select(lha_code),
  join = st_within
) %>%
  st_drop_geometry()

# Diagnostic
cat("Total DAs:", nrow(da_lha_lookup), "\n")
cat("DAs matched to an LHA:", sum(!is.na(da_lha_lookup$lha_code)), "\n")
cat("DAs unmatched:", sum(is.na(da_lha_lookup$lha_code)), "\n")

# Get the unmatched DAs
unmatched_das <- da_lha_lookup %>% filter(is.na(lha_code))

# Find nearest LHA for each unmatched DA centroid
unmatched_centroids <- da_centroids %>%
  filter(DGUID %in% unmatched_das$DGUID)

nearest_lha <- st_join(
  unmatched_centroids,
  lha_boundaries %>% select(lha_code),
  join = st_nearest_feature
) %>%
  st_drop_geometry()

# Patch into lookup
da_lha_lookup <- da_lha_lookup %>%
  rows_update(nearest_lha, by = "DGUID")

# Verify
cat("DAs unmatched after fix:", sum(is.na(da_lha_lookup$lha_code)), "\n")

# ------------------------------------------------------------
# Run MoH sensitivity across all years

# Age band columns to sum for pop_total_40plus
age_bands_40plus <- c(
  "pop_total_40_44", "pop_total_45_49", "pop_total_50_54",
  "pop_total_55_59", "pop_total_60_64", "pop_total_65_69",
  "pop_total_70_74", "pop_total_75_79", "pop_total_80_84",
  "pop_total_85_89", "pop_total_90_94", "pop_total_95_99",
  "pop_total_100plus"
)

sensitivity_moh <- map_dfr(years, function(yr) {
  
  daily <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds"))
  fireseason <- daily %>%
    filter(
      date >= as.Date(paste0(yr, "-05-01")),
      date <= as.Date(paste0(yr, "-09-30"))
    ) %>%
    group_by(DGUID) %>%
    summarise(
      fireseason_mean_pm25 = mean(pm25_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  annual <- st_read(
    paste0("da_annual_pm25_summary_", yr, ".gpkg"),
    quiet = TRUE
  ) %>%
    st_drop_geometry()
  
  da_moh <- annual %>%
    # Sum age bands to get pop_total_40plus
    mutate(
      pop_total_40plus = rowSums(
        across(all_of(age_bands_40plus)),
        na.rm = TRUE
      )
    ) %>%
    left_join(fireseason,      by = "DGUID") %>%
    left_join(da_lha_lookup,   by = c("DGUID", "DAUID")) %>%
    left_join(moh_rates_combined, by = "lha_code") %>%
    mutate(
      delta_pm25     = pmax(fireseason_mean_pm25 - counterfactual, 0),
      RR             = exp(beta_central * delta_pm25),
      PAF            = (RR - 1) / RR,
      expected_total = pop_total_40plus * rate_per_person,
      attrib_total   = expected_total * PAF,
      year           = yr
    )
  
  da_moh
})

# ------------------------------------------------------------
# Summarise MoH sensitivity by year

sensitivity_moh_summary <- sensitivity_moh %>%
  group_by(year) %>%
  summarise(
    attributable_cases = sum(attrib_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(scenario = "MoH LHA (40+)")

print(sensitivity_moh_summary)

# ------------------------------------------------------------
# Update uncertainty table to include MoH scenario

uncertainty_by_year_extended <- bind_rows(
  all_years_results %>%
    mutate(scenario = "central") %>%
    group_by(year, scenario) %>%
    summarise(attributable_cases = sum(attrib_total, na.rm = TRUE), .groups = "drop"),
  sensitivity_low %>%
    mutate(scenario = "low") %>%
    group_by(year, scenario) %>%
    summarise(attributable_cases = sum(attrib_total, na.rm = TRUE), .groups = "drop"),
  sensitivity_high %>%
    mutate(scenario = "high") %>%
    group_by(year, scenario) %>%
    summarise(attributable_cases = sum(attrib_total, na.rm = TRUE), .groups = "drop"),
  sensitivity_moh_summary %>%
    rename(scenario_label = scenario) %>%
    mutate(scenario = "moh") %>%
    select(year, scenario, attributable_cases)
) %>%
  group_by(year, scenario) %>%
  summarise(attributable_cases = sum(attributable_cases, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = scenario, values_from = attributable_cases) %>%
  mutate(across(c(low, central, high, moh), ~round(., 1)))

cat("\n=== EXTENDED UNCERTAINTY TABLE (WITH MOH SENSITIVITY) ===\n")
print(uncertainty_by_year_extended)

cat("\nCumulative low:     ", round(sum(uncertainty_by_year_extended$low),     1), "\n")
cat("Cumulative central: ", round(sum(uncertainty_by_year_extended$central), 1), "\n")
cat("Cumulative high:    ", round(sum(uncertainty_by_year_extended$high),    1), "\n")
cat("Cumulative MoH:     ", round(sum(uncertainty_by_year_extended$moh),     1), "\n")

# ------------------------------------------------------------
#  Plot 1
plot_data_extended <- uncertainty_by_year_extended %>%
  mutate(
    is_small            = central < 50,
    uncertainty_label_y = ifelse(is_small, high + 25, high + 8)
  )

p1_extended <- ggplot(plot_data_extended, aes(x = factor(year))) +
  
  # Uncertainty range (low-high beta)
  geom_rect(
    aes(
      xmin = as.numeric(factor(year)) - 0.35,
      xmax = as.numeric(factor(year)) + 0.35,
      ymin = low,
      ymax = high
    ),
    fill = "#ca706a", alpha = 0.25
  ) +
  
  # Central estimate bar
  geom_col(
    aes(y = central),
    fill = "#ca706a", alpha = 0.85, width = 0.5
  ) +
  
  # Central estimate label
  geom_text(
    aes(y = central, label = round(central, 0)),
    vjust = -0.3, size = 3.5, fontface = "bold", color = "#1a1f2e"
  ) +
  
  # Uncertainty range label
  geom_text(
    aes(y = uncertainty_label_y,
        label = paste0("(", round(low, 0), "–", round(high, 0), ")")),
    vjust = 0, size = 2.8, color = "#1a1f2e"
  ) +
  
  # Dotted connector for small bars
  geom_segment(
    data = plot_data_extended %>% filter(is_small),
    aes(
      x = as.numeric(factor(year)), xend = as.numeric(factor(year)),
      y = high + 3, yend = high + 20
    ),
    color = "grey60", linewidth = 0.3, linetype = "dotted"
  ) +
  
  # MoH sensitivity — diamond points connected by line
  geom_line(
    aes(y = moh, group = 1),
    color = "#2A9D8F", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_point(
    aes(y = moh, shape = "BC MoH LHA rates (40+)"),
    color = "#2A9D8F", size = 3, fill = "#2A9D8F"
  ) +
  geom_text(
    aes(y = moh, label = round(moh, 0)),
    vjust = -0.7, size = 2.8, color = "#2A9D8F"
  ) +
  
  scale_shape_manual(
    values = c("BC MoH LHA rates (40+)" = 18),
    name   = "Sensitivity analysis"
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25)),
    labels = scales::comma
  ) +
  
  labs(
    title    = "Estimated Alzheimer's Cases Attributable to Wildfire PM\u2082.\u2085",
    subtitle = "British Columbia, 2021\u20132026 Fire Seasons (May\u2013September) \u00b7 Population Aged 65+",
    x        = "Year",
    y        = "Attributable Alzheimer's cases",
    caption  = paste0(
      "Bars: central estimate (Fasoro et al. 2025 / CCDSS rates, 65+ age bands). ",
      "Shaded range: \u03b2 sensitivity bounds (Chen et al. 2017 HR 1.03\u20131.05 per 4.8 \u00b5g/m\u00b3 IQR).\n",
      "Teal diamonds: BC MoH LHA-specific sensitivity (single 40+ rate; population base differs from primary analysis).\n",
      "Counterfactual: 15 \u00b5g/m\u00b3 (WHO AQG). Differences between primary and MoH estimates reflect both rate source and age-group definition."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 20, hjust = 0),
    plot.subtitle      = element_text(color = "grey40", size = 15, hjust = 0),
    plot.caption       = element_text(color = "grey50", size = 9, hjust = 0,
                                      lineheight = 1.4),
    legend.position    = c(0.85, 0.85),
    legend.background  = element_rect(fill = "white", color = "grey90"),
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "grey30"),
    axis.title         = element_text(color = "grey30", size = 10),
    plot.background    = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

print(p1_extended)

ggsave(
  "attributable_cases_by_year_with_moh_sensitivity.png",
  plot   = p1_extended,
  width  = 11,
  height = 7,
  dpi    = 300,
  bg     = "#ffffffC9"
)
# ============================================================
# PLOT 2: Choropleth — cumulative attributable cases per DA
top5_cumulative <- regional_by_year %>%
  group_by(region) %>%
  summarise(
    total_attributable = sum(attributable_cases, na.rm = TRUE),
    mean_pm25          = weighted.mean(mean_pm25, attributable_cases, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_attributable)) %>%
  slice_head(n = 5) %>%
  mutate(
    region_short = case_when(
      region == "Alberni-Clayoquot / Central Island" ~ "Alberni-Clayoquot /\nCentral Island",
      region == "Thompson-Nicola / Kamloops"         ~ "Thompson-Nicola /\nKamloops",
      region == "Peace River / Northeast BC"         ~ "Peace River /\nNortheast BC",
      region == "Central Okanagan / Kootenay"        ~ "Central Okanagan /\nKootenay",
      region == "Capital / Greater Victoria"         ~ "Capital /\nGreater Victoria",
      TRUE ~ region
    ),
    # Scale bar width relative to max
    bar_pct = total_attributable / max(total_attributable)
  )

# Total cumulative burden
total_burden <- sum(summary_by_year$attributable_cases)

# Top 5 cumulative regions
top5_share <- regional_cumulative %>%
  slice_head(n = 5) %>%
  summarise(
    top5_cases = sum(total_attributable),
    pct_of_total = top5_cases / total_burden * 100
  )

cat("Top 5 regions cumulative cases:", round(top5_share$top5_cases, 1), "\n")
cat("% of total burden in top 5:    ", round(top5_share$pct_of_total, 1), "%\n")

print(regional_cumulative %>% 
        slice_head(n = 5) %>% 
        select(region, total_attributable))

# ============================================================
# MAIN MAP
map_plot <- ggplot() +
  
  geom_sf(data = bc_outline, fill = "grey92", color = "grey70", linewidth = 0.3) +
  geom_sf(data = da_map %>% filter(cumulative_attrib == 0 | is.na(cumulative_attrib)),
          fill = "grey88", color = NA) +
  geom_sf(data = da_map %>% filter(cumulative_attrib > 0),
          aes(fill = attrib_per_1000_capped), color = NA) +
  
  scale_fill_viridis_c(
    option    = "magma",
    direction = -1,
    name      = "Attributable cases\nper 1,000 pop. 65+\n(2021–2026 cumulative)",
    na.value  = "grey88",
    labels    = scales::number_format(accuracy = 0.1),
    guide     = guide_colorbar(
      title.position = "top",
      barwidth       = unit(0.8, "cm"),
      barheight      = unit(5, "cm"),
      ticks          = TRUE
    )
  ) +
  
  labs(
    title    = paste0("Wildfire Smoke Exposure Concentrates Alzheimer's Risk in BC's Interior\n",
    "and South Coast"),
    subtitle = "British Columbia dissemination areas · 2021–2026 fire seasons · Population aged 65+"
  ) +
  
  theme_void(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 20, hjust = 0,
                                      margin = margin(b = 4)),
    plot.subtitle      = element_text(color = "grey40", size = 14, hjust = 0,
                                      margin = margin(b = 8)),
    legend.position    = c(0.03, 0.72),
    legend.justification = c(0, 1),
    legend.background      = element_rect(fill = "white", color = NA, linewidth = 0),
    legend.title           = element_text(size = 9, lineheight = 1.4),
    legend.text            = element_text(size = 8.5),
    legend.margin          = margin(4, 6, 4, 6),
    plot.background        = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin            = margin(12, 4, 4, 4)   # reduced bottom margin
  ) +
  
  coord_sf(crs = 3005)

# ============================================================
# RIGHT PANEL: Top 5 regions bar chart
panel_plot <- ggplot(top5_cumulative,
                     aes(y = reorder(region_short, total_attributable))) +
  
  # Background track
  geom_col(
    aes(x = max_cases * 1.15),
    fill = "grey92", width = 0.6
  ) +
  
  # Attributable cases bar
  geom_col(
    aes(x = total_attributable),
    fill = "#ca706a", alpha = 0.85, width = 0.6
  ) +
  
  # Cases label inside/end of bar
  geom_text(
    aes(
      x     = total_attributable,
      label = round(total_attributable, 0)
    ),
    hjust = -0.2, size = 3.2, fontface = "bold", color = "#1a1f2e"
  ) +
  
  # Mean PM2.5 label on right
  geom_text(
    aes(
      x     = max_cases * 1.35,
      label = paste0(round(mean_pm25, 1), " µg/m³")
    ),
    hjust = 0.5, size = 2.9, color = "grey40"
  ) +
  
  # Column headers
  annotate("text", x = max_cases * 0.5,  y = 5.7,
           label = "Attributable cases", size = 2.8,
           color = "grey30", fontface = "bold", hjust = 0.5) +
  annotate("text", x = max_cases * 1.35, y = 5.7,
           label = "Mean PM₂.₅", size = 2.8,
           color = "grey30", fontface = "bold", hjust = 0.5) +
  
  scale_x_continuous(
    limits = c(0, max_cases * 1.55),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_discrete(expand = expansion(add = c(0.5, 0.8))) +
  
  labs(
    title    = "Top 5 Regions",
    subtitle = "Cumulative 2021–2026"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle    = element_text(color = "grey40", size = 8.5, hjust = 0,
                                    margin = margin(b = 8)),
    plot.caption     = element_text(color = "grey50", size = 7, hjust = 0,
                                    lineheight = 1.4, margin = margin(t = 10)),
    axis.text.y      = element_text(size = 8, color = "grey30", lineheight = 1.3),
    axis.text.x      = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    plot.background  = element_rect(fill = "#ffffffC9", color = NA),
    panel.background = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin      = margin(12, 12, 4, 4)
  )

# ============================================================
# COMBINE 
combined <- map_plot + panel_plot +
  plot_layout(widths = c(2.2, 1)) +
  plot_annotation(
    caption  = paste0(
      "Attributable cases per 1,000 pop. 65+, summed across 2021–2026 fire seasons.Grey areas: fire-season mean PM₂.₅ below 15 µg/m³ WHO counterfactual in all years.\n",
      "Chen et al. (2017) CRF. BlueSky fire weather model. Statistics Canada 2021 Census."
    ),
    theme = theme(
      plot.caption    = element_text(size = 7.5, color = "grey50", hjust = 0,
                                     lineheight = 1.5, margin = margin(t = 6)),
      plot.background = element_rect(fill = "#ffffffC9", color = NA),
      panel.background = element_rect(fill = "#ffffffC9", color = NA),
      plot.margin     = margin(0, 0, 8, 0)
    )
  )

ggsave(
  "bc_choropleth_with_panel.png",
  plot   = combined,
  width  = 14,
  height = 9.5,
  dpi    = 300,
  bg     = "#ffffffC9"
)

# ============================================================
# PLOT 3: Age group breakdown of attributable cases
age_breakdown <- all_years_results %>%
  summarise(
    `65–74` = sum(attrib_65_74,  na.rm = TRUE),
    `75–84` = sum(attrib_75_84,  na.rm = TRUE),
    `85+`   = sum(attrib_85plus, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "age_group", values_to = "cases") %>%
  mutate(
    age_group = factor(age_group, levels = c("65–74", "75–84", "85+")),
    pct       = round((cases / sum(cases)) * 100, 1),
    label     = paste0(round(cases, 0), "\n(", pct, "%)")
  )

p3 <- ggplot(age_breakdown, aes(x = age_group, y = cases, fill = age_group)) +
  
  geom_col(width = 0.55, alpha = 0.88) +
  
  geom_text(
    aes(label = label),
    vjust = -0.4,
    size  = 3.4,
    color = "grey20",
    lineheight = 1.3
  ) +
  
  scale_fill_manual(
    values = c(
      "65–74" = "#f2b49b",
      "75–84" = "#ca706a",
      "85+"   = "#7a2c2c"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2)),
    labels = scales::comma
  ) +
  
  labs(
    title    = "Attributable Cases by Age Group",
    subtitle = "2021–2026 · cumulative",
    x        = "Age group",
    y        = "Attributable dementia cases",
    caption  = "Cases summed across 2021–2026 fire seasons (May–September).\nChen et al. (2017) CRF. Fasoro et al. (2025) / CCDSS baseline rates."
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle      = element_text(color = "grey40", size = 9, hjust = 0,
                                      margin = margin(b = 10)),
    plot.caption       = element_text(color = "grey50", size = 7.5, hjust = 0,
                                      lineheight = 1.4, margin = margin(t = 10)),
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "grey30", size = 10),
    axis.title         = element_text(color = "grey30", size = 9),
    plot.background    = element_rect(fill = "#ffffffC9", color = NA),
    panel.background   = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

ggsave(
  "age_group_breakdown.png",
  plot = p3, width = 6, height = 5, dpi = 300,
  bg = "#ffffffC9"
)

# ============================================================
# SEX BREAKDOWN ANALYSIS
moh_rates_sex <- read_csv(
  "BC Community Health Data Download- age standardized dementia.csv",
  show_col_types = FALSE
) %>%
  rename(
    lha_code      = `Jurisdiction Code`,
    fiscal_year   = `Year(s)`,
    sex           = Sex,
    rate_per_1000 = `Indicator Value`
  ) %>%
  # Keep only Female and Male rows (drop Total and blanks)
  filter(sex %in% c("Female", "Male")) %>%
  filter(!is.na(lha_code), !is.na(rate_per_1000)) %>%
  filter(fiscal_year %in% c(
    "FY 2016.2017", "FY 2017.2018", "FY 2018.2019",
    "FY 2019.2020", "FY 2020.2021"
  )) %>%
  group_by(lha_code, sex) %>%
  summarise(
    rate_per_1000 = mean(rate_per_1000, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lha_code        = as.character(lha_code),
    rate_per_person = rate_per_1000 / 1000
  ) %>%
  # Pivot wide so each LHA has one row with male and female rates
  pivot_wider(
    names_from  = sex,
    values_from = c(rate_per_1000, rate_per_person),
    names_glue  = "{tolower(sex)}_{.value}"
  )

# Age band columns for male/female 40+ populations

age_bands_male_40plus <- c(
  "male40.to.44.years", "male45.to.49.years", "male50.to.54.years",
  "male55.to.59.years", "male60.to.64.years", "male65.to.69.years",
  "male70.to.74.years", "male75.to.79.years", "male80.to.84.years",
  "male85.to.89.years", "male90.to.94.years", "male95.to.99.years",
  "male100.years.and.over"
)

age_bands_female_40plus <- c(
  "female40.to.44.years", "female45.to.49.years", "female50.to.54.years",
  "female55.to.59.years", "female60.to.64.years", "female65.to.69.years",
  "female70.to.74.years", "female75.to.79.years", "female80.to.84.years",
  "female85.to.89.years", "female90.to.94.years", "female95.to.99.years",
  "female100.years.and.over"
)

# ------------------------------------------------------------
# Run sex-stratified CRF across all years
sex_results <- map_dfr(years, function(yr) {
  
  cat("Sex breakdown processing:", yr, "...\n")
  
  daily <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds"))
  
  fireseason <- daily %>%
    filter(
      date >= as.Date(paste0(yr, "-05-01")),
      date <= as.Date(paste0(yr, "-09-30"))
    ) %>%
    group_by(DGUID) %>%
    summarise(
      fireseason_mean_pm25 = mean(pm25_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  annual <- st_read(
    paste0("da_annual_pm25_summary_", yr, ".gpkg"),
    quiet = TRUE
  ) %>%
    st_drop_geometry()
  
  annual %>%
    # Sum age bands separately for each sex
    mutate(
      pop_male_40plus   = rowSums(across(all_of(age_bands_male_40plus)),   na.rm = TRUE),
      pop_female_40plus = rowSums(across(all_of(age_bands_female_40plus)), na.rm = TRUE)
    ) %>%
    left_join(fireseason,       by = "DGUID") %>%
    left_join(da_lha_lookup,    by = c("DGUID", "DAUID")) %>%
    left_join(moh_rates_sex,    by = "lha_code") %>%
    mutate(
      delta_pm25      = pmax(fireseason_mean_pm25 - counterfactual, 0),
      RR              = exp(beta_central * delta_pm25),
      PAF             = (RR - 1) / RR,
      
      # Expected cases by sex
      expected_male   = pop_male_40plus   * male_rate_per_person,
      expected_female = pop_female_40plus * female_rate_per_person,
      
      # Attributable cases by sex
      attrib_male     = expected_male   * PAF,
      attrib_female   = expected_female * PAF,
      attrib_total    = attrib_male + attrib_female,
      
      # Female share of burden in this DA
      female_fraction = ifelse(attrib_total > 0, attrib_female / attrib_total, NA),
      
      year = yr
    )
})

# ------------------------------------------------------------
# BC-wide summary by sex and year

bc_sex_by_year <- sex_results %>%
  group_by(year) %>%
  summarise(
    attrib_male         = sum(attrib_male,   na.rm = TRUE),
    attrib_female       = sum(attrib_female, na.rm = TRUE),
    attrib_total        = sum(attrib_total,  na.rm = TRUE),
    pct_female          = attrib_female / attrib_total * 100,
    pct_male            = attrib_male   / attrib_total * 100,
    .groups = "drop"
  )

cat("\n=== BC-WIDE ATTRIBUTABLE CASES BY SEX AND YEAR ===\n")
print(bc_sex_by_year %>%
        mutate(across(c(attrib_male, attrib_female, attrib_total), ~round(., 1)),
               across(c(pct_male, pct_female), ~round(., 1))))

# Cumulative across all years
cumulative_sex <- sex_results %>%
  summarise(
    attrib_male   = sum(attrib_male,   na.rm = TRUE),
    attrib_female = sum(attrib_female, na.rm = TRUE),
    attrib_total  = sum(attrib_total,  na.rm = TRUE)
  ) %>%
  mutate(
    pct_male   = attrib_male   / attrib_total * 100,
    pct_female = attrib_female / attrib_total * 100
  )

cat("\n=== CUMULATIVE SEX BREAKDOWN 2021-2026 ===\n")
cat("Female attributable cases:", round(cumulative_sex$attrib_female, 1),
    paste0("(", round(cumulative_sex$pct_female, 1), "%)"), "\n")
cat("Male attributable cases:  ", round(cumulative_sex$attrib_male,   1),
    paste0("(", round(cumulative_sex$pct_male,   1), "%)"), "\n")
cat("Total:                    ", round(cumulative_sex$attrib_total,  1), "\n")

# ------------------------------------------------------------
# DA-level sex variance 
da_sex_summary <- sex_results %>%
  group_by(DGUID, DAUID) %>%
  summarise(
    attrib_male     = sum(attrib_male,   na.rm = TRUE),
    attrib_female   = sum(attrib_female, na.rm = TRUE),
    attrib_total    = sum(attrib_total,  na.rm = TRUE),
    female_fraction = attrib_female / attrib_total,
    sex_gap         = abs(attrib_female - attrib_male),
    dominant_sex    = if_else(attrib_female >= attrib_male, "Female", "Male"),
    .groups = "drop"
  ) %>%
  filter(attrib_total > 0)

cat("\n=== DA-LEVEL SEX DISTRIBUTION ===\n")
cat("Median female share of burden:", round(median(da_sex_summary$female_fraction, na.rm = TRUE) * 100, 1), "%\n")
cat("Mean female share of burden:  ", round(mean(da_sex_summary$female_fraction,   na.rm = TRUE) * 100, 1), "%\n")
cat("DAs where female burden > 60%:", sum(da_sex_summary$female_fraction > 0.6, na.rm = TRUE), "\n")
cat("DAs where male burden > 60%:  ", sum(da_sex_summary$female_fraction < 0.4, na.rm = TRUE), "\n")

# Regional sex breakdown — which regions have highest female vs male burden
regional_sex <- sex_results %>%
  filter(delta_pm25 > 0) %>%
  mutate(
    cd_code = substr(DAUID, 1, 4),
    region = case_when(
      cd_code == "5901" ~ "Capital / Greater Victoria",
      cd_code == "5903" ~ "Cowichan Valley",
      cd_code == "5905" ~ "Nanaimo",
      cd_code == "5907" ~ "Alberni-Clayoquot / Central Island",
      cd_code == "5909" ~ "Comox Valley",
      cd_code == "5911" ~ "Strathcona",
      cd_code == "5915" ~ "Metro Vancouver",
      cd_code == "5919" ~ "Squamish-Lillooet",
      cd_code == "5921" ~ "Fraser Valley",
      cd_code == "5933" ~ "Thompson-Nicola / Kamloops",
      cd_code == "5935" ~ "Central Okanagan / Kootenay",
      cd_code == "5937" ~ "Okanagan-Similkameen",
      cd_code == "5939" ~ "Peace River / Northeast BC",
      cd_code == "5941" ~ "Cariboo",
      cd_code == "5943" ~ "Skeena",
      cd_code == "5945" ~ "Kitimat-Stikine",
      cd_code == "5947" ~ "Bulkley-Nechako",
      cd_code == "5949" ~ "Fraser-Fort George / Prince George",
      cd_code == "5953" ~ "Central Okanagan",
      cd_code == "5955" ~ "North Okanagan",
      cd_code == "5957" ~ "Columbia-Shuswap",
      cd_code == "5959" ~ "East Kootenay",
      TRUE ~ paste0("CD ", cd_code)
    )
  ) %>%
  group_by(region) %>%
  summarise(
    attrib_male     = sum(attrib_male,   na.rm = TRUE),
    attrib_female   = sum(attrib_female, na.rm = TRUE),
    attrib_total    = sum(attrib_total,  na.rm = TRUE),
    pct_female      = attrib_female / attrib_total * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(attrib_total))

cat("\n=== REGIONAL SEX BREAKDOWN (CUMULATIVE 2021-2026) ===\n")
print(regional_sex %>%
        mutate(across(c(attrib_male, attrib_female, attrib_total), ~round(., 1)),
               pct_female = round(pct_female, 1)),
      n = 20)

# ------------------------------------------------------------
# Plot — BC-wide attributable cases by sex per year

plot_sex <- bc_sex_by_year %>%
  select(year, attrib_male, attrib_female) %>%
  pivot_longer(
    cols      = c(attrib_male, attrib_female),
    names_to  = "sex",
    values_to = "cases"
  ) %>%
  mutate(
    sex = if_else(sex == "attrib_male", "Male", "Female"),
    sex = factor(sex, levels = c("Female", "Male"))
  )

p_sex <- ggplot(plot_sex, aes(x = factor(year), y = cases, fill = sex)) +
  
  geom_col(position = "dodge", width = 0.6, alpha = 0.88) +
  
  geom_text(
    aes(label = round(cases, 0)),
    position = position_dodge(width = 0.6),
    vjust    = -0.4,
    size     = 3,
    color    = "grey20"
  ) +
  
  scale_fill_manual(
    values = c("Female" = "#ca706a", "Male" = "#6baed6"),
    name   = NULL
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2)),
    labels = scales::comma
  ) +
  
  labs(
    title    = "Wildfire PM\u2082.\u2085-Attributable Alzheimer's Cases by Sex",
    subtitle = "British Columbia, 2021\u20132026 Fire Seasons \u00b7 Population Aged 40+",
    x        = "Year",
    y        = "Attributable cases",
    caption  = paste0(
      "BC Ministry of Health LHA-specific age-standardised dementia incidence rates (40+) by sex,\n",
      "averaged across FY 2016/17\u20132020/21. Chen et al. (2017) CRF applied equally to both sexes\n",
      "(HR = 1.04 per 4.8 \u00b5g/m\u00b3 IQR). Sex categories reflect Statistics Canada census definitions.\n",
      "Differences reflect population size and baseline incidence, not differential PM\u2082.\u2085 susceptibility."
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle      = element_text(color = "grey40", size = 9, hjust = 0,
                                      margin = margin(b = 8)),
    plot.caption       = element_text(color = "grey50", size = 7.5, hjust = 0,
                                      lineheight = 1.4, margin = margin(t = 10)),
    legend.position    = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "grey30"),
    axis.title         = element_text(color = "grey30", size = 9),
    plot.background    = element_rect(fill = "#ffffffC9", color = NA),
    panel.background   = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

print(p_sex)

ggsave(
  "attributable_cases_by_sex.png",
  plot   = p_sex,
  width  = 9,
  height = 6,
  dpi    = 300,
  bg     = "#ffffffC9"
)

# Save DA-level sex summary
write_csv(da_sex_summary, "da_sex_burden_summary_2021_2026.csv")

# ============================================================
# PLOT 4: Concentration-response curve

# Generate a range of PM2.5 values to plot the curve across
pm25_range <- seq(0, 100, by = 0.5)

crf_curve <- tibble(pm25 = pm25_range) %>%
  mutate(
    delta_pm25   = pmax(pm25 - counterfactual, 0),
    RR_central   = exp(beta_central * delta_pm25),
    RR_low       = exp(beta_low     * delta_pm25),
    RR_high      = exp(beta_high    * delta_pm25),
    PAF_central  = (RR_central - 1) / RR_central,
    PAF_low      = (RR_low     - 1) / RR_low,
    PAF_high     = (RR_high    - 1) / RR_high
  )

# Pull actual DA-level fire-season exposures from 2021 and 2023
# for the rug/distribution overlay
exposure_dist <- all_years_results %>%
  filter(fireseason_mean_pm25 > 0) %>%
  select(year, fireseason_mean_pm25, pop_total_65plus) %>%
  mutate(year = factor(year))

p4 <- ggplot() +
  
  # Uncertainty ribbon
  geom_ribbon(
    data = crf_curve,
    aes(x = pm25, ymin = PAF_low, ymax = PAF_high, fill = "95% uncertainty range"),
    alpha = 0.2
  ) +
  
  # Central CRF line
  geom_line(
    data = crf_curve,
    aes(x = pm25, y = PAF_central, linetype = "Central estimate"),
    color = "#ca706a", linewidth = 1.2
  ) +
  
  scale_fill_manual(
    values = c("95% uncertainty range" = "#ca706a"),
    name = NULL,
    guide = guide_legend(order = 1)
  ) +
  
  scale_linetype_manual(
    values = c("Central estimate" = "solid"),
    name = NULL,
    guide = guide_legend(order = 1)
  ) +
  
  # WHO counterfactual vertical line
  geom_vline(
    xintercept = 15,
    linetype = "dashed", color = "grey40", linewidth = 0.7
  ) +
  
  # Counterfactual label
  annotate(
    "text", x = 15.8, y = 0.72,
    label = "WHO AQG\n15 µg/m³",
    hjust = 0, size = 2.8, color = "grey40", lineheight = 1.3
  ) +
  
  # Rug layers — color mapped inside aes() so it feeds scale_color_manual
  geom_rug(
    data = exposure_dist %>% filter(year == "2021"),
    aes(x = fireseason_mean_pm25, color = "2021 DA exposures"),
    alpha = 0.15, sides = "b", length = unit(0.04, "npc")
  ) +
  geom_rug(
    data = exposure_dist %>% filter(year == "2022"),
    aes(x = fireseason_mean_pm25, color = "2022 DA exposures"),
    alpha = 0.15, sides = "b", length = unit(0.04, "npc")
  ) +
  geom_rug(
    data = exposure_dist %>% filter(year == "2023"),
    aes(x = fireseason_mean_pm25, color = "2023 DA exposures"),
    alpha = 0.15, sides = "b", length = unit(0.04, "npc")
  ) +
  geom_rug(
    data = exposure_dist %>% filter(year == "2024"),
    aes(x = fireseason_mean_pm25, color = "2024 DA exposures"),
    alpha = 0.15, sides = "b", length = unit(0.04, "npc")
  ) +
  geom_rug(
    data = exposure_dist %>% filter(year == "2025"),
    aes(x = fireseason_mean_pm25, color = "2025 DA exposures"),
    alpha = 0.15, sides = "b", length = unit(0.04, "npc")
  ) +
  geom_rug(
    data = exposure_dist %>% filter(year == "2026"),
    aes(x = fireseason_mean_pm25, color = "2026 DA exposures"),
    alpha = 0.15, sides = "b", length = unit(0.04, "npc")
  ) +

  
  scale_color_manual(
    name = "DA exposure distribution",
    values = c(
      "2021 DA exposures" = "#7a2c2c",
      "2022 DA exposures" = "#ca706a",
      "2023 DA exposures" = "#e8a598",
      "2024 DA exposures" = "rosybrown2",
      "2025 DA exposures" = "rosybrown1",
      "2026 DA exposures" = "snow3"
    ),
    guide = guide_legend(
      order = 2,
      override.aes = list(alpha = 1)
    )
  ) +
  
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 15, 25, 50, 75, 100),
    labels = c("0", "15", "25", "50", "75", "100")
  ) +
  
  scale_y_continuous(
    limits = c(0, 0.85),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  coord_cartesian(xlim = c(0, 100), clip = "off") +
  
  labs(
    title    = "Concentration–Response Function",
    subtitle = "Population attributable fraction (PAF) by fire-season PM₂.₅ · Chen et al. (2017)",
    x        = "Fire-season mean PM₂.₅ (µg/m³)",
    y        = "Population attributable fraction",
    caption  = paste0(
      "Shaded band: 95% uncertainty range across β sensitivity bounds (HR 1.03–1.05 per 4.8 µg/m³ IQR).\n",
      "Rug marks show distribution of DA-level fire-season mean PM₂.₅ for 2021–2025. (x-axis truncated at 100 µg/m³; 94 DAs omitted)\n",
      "PAF = 0 below counterfactual of 15 µg/m³ (WHO 24-hr AQG)."
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "right",
    legend.box         = "vertical",
    legend.title       = element_text(size = 8, color = "grey30"),
    legend.text        = element_text(size = 7.5, color = "grey30"),
    legend.key.height  = unit(0.4, "cm"),
    legend.key.width  = unit(1.2, "cm"),   
    legend.spacing.y   = unit(0.2, "cm"),
    legend.background  = element_rect(fill = "#ffffffC9", color = NA),
    plot.title         = element_text(face = "bold", size = 20, hjust = 0),
    plot.subtitle      = element_text(color = "grey40", size = 9, hjust = 0,
                                      margin = margin(b = 10)),
    plot.caption       = element_text(color = "grey50", size = 7.5, hjust = 0,
                                      lineheight = 1.4, margin = margin(t = 10)),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "grey30", size = 9),
    axis.title         = element_text(color = "grey30", size = 9),
    plot.background    = element_rect(fill = "#ffffffC9", color = NA),
    panel.background   = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

ggsave(
  "concentration_response_curve.png",
  plot = p4, width = 7, height = 5, dpi = 300,
  bg = "#ffffffC9"
)
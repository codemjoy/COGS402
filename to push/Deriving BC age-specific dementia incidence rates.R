library(sf)
library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)
library(ggtext)
library(bcmaps)
library(purrr)
library(dplyr)
library(scales)

# PARAMETERS 
# primary crf: Chen et al 2017
beta_central <- log(1.04) / 4.8
beta_low     <- log(1.03) / 4.8
beta_high    <- log(1.05) / 4.8

# exploratory upper bound - Elser et al. 2025
beta_elser         <- log(1.12)
beta_elser_low     <- log(0.98)   
beta_elser_high    <- log(1.28) 

# High-end cohort estimate Peters et al. 2019, Carey et al. 2018
beta_highcohort <- log(1.09) / 10

# counterfactuals 
cf_annual_emp <- read_csv("annual_cf_stats.csv") %>%
  select(fire_year, cf_median, cf_p25, cf_p75) %>%
  filter(fire_year %in% years) %>%
  mutate(cf_se = (cf_p75 - cf_p25) / 2)

cf_central_emp <- setNames(cf_annual_emp$cf_median, cf_annual_emp$fire_year)

#sensitivity comparison 
CF_WHO  <- 15
CF_5    <- 5
CF_NONE <- 0

PM25_DAILY_CAP <- 150

# Incidence rates (Fasoro et al. 2025 / CCDSS)
primary_rates <- list(
  rate_65_74  =  6.10 / 1000,
  rate_75_84  = 34.50 / 1000,
  rate_85plus = 44.68 / 1000
)

years <- 2021:2025

get_fireseason <- function(yr, cap = PM25_DAILY_CAP,
                           pm25_col = "pm25_mean") {
  daily <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds"))
  
  if (!pm25_col %in% names(daily)) {
    stop(paste0(
      "Column '", pm25_col, "' not found in RDS for ", yr,
      ". Re-run rasterization script to generate afternoon bands."
    ))
  }
  
  daily %>%
    filter(
      date >= as.Date(paste0(yr, "-05-01")),
      date <= as.Date(paste0(yr, "-09-30"))
    ) %>%
    mutate(
      pm25_val = pmin(.data[[pm25_col]], cap)
    ) %>%
    group_by(DGUID) %>%
    summarise(
      fireseason_mean_pm25       = mean(pm25_val,   na.rm = TRUE),
      fireseason_days_above_15   = sum(above_15,    na.rm = TRUE),
      fireseason_days_above_37.5 = sum(above_37.5,  na.rm = TRUE),
      fireseason_days_above_50   = sum(above_50,    na.rm = TRUE),
      n_fireseason_days          = n(),
      n_days_capped              = sum(.data[[pm25_col]] > cap, na.rm = TRUE),
      .groups = "drop"
    )
}

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

# determine how many days affected by cap 
cap_audit <- map_dfr(years, function(yr) {
  readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds")) %>%
    filter(
      date >= as.Date(paste0(yr, "-05-01")),
      date <= as.Date(paste0(yr, "-09-30"))
    ) %>%
    summarise(
      total_da_days = n(),
      days_capped   = sum(pm25_mean > PM25_DAILY_CAP, na.rm = TRUE),
      das_affected  = n_distinct(DGUID[pm25_mean > PM25_DAILY_CAP]),
      max_raw_value = round(max(pm25_mean, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    mutate(
      year       = yr,
      pct_capped = round(days_capped / total_da_days * 100, 3)
    )
}) %>%
  select(year, total_da_days, days_capped, pct_capped, das_affected, max_raw_value)

print(cap_audit)


# Chen central beta, 24-hr mean
all_years_results <- map_dfr(years, function(yr) {
  cf_emp <- cf_central_emp[[as.character(yr)]]
  
  fireseason <- get_fireseason(yr, cap = PM25_DAILY_CAP, pm25_col = "pm25_mean")
  annual     <- st_read(paste0("da_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE)
  
  annual %>%
    left_join(fireseason, by = "DGUID") %>%
    calculate_crf(beta_central, cf_emp, primary_rates) %>%
    mutate(year = yr, cf_used = cf_emp)
})

# multiple sensitivity scenarios loop
scenarios <- list(
  # Primary CRF, CF=15 (beta bounds)
  chen_central_cf15 = list(beta = beta_central, cf = CF_WHO,  col = "pm25_mean",
                           label = "Chen central (CF=15 WHO AQG)"),
  chen_low_cf15     = list(beta = beta_low,     cf = CF_WHO,  col = "pm25_mean",
                           label = "Chen low β (CF=15)"),
  chen_high_cf15    = list(beta = beta_high,    cf = CF_WHO,  col = "pm25_mean",
                           label = "Chen high β (CF=15)"),
  
  # Dual counterfactual — Chen central
  chen_central_cf5  = list(beta = beta_central, cf = CF_5,    col = "pm25_mean",
                           label = "Chen central (CF=5)"),
  chen_low_cf5      = list(beta = beta_low,     cf = CF_5,    col = "pm25_mean",
                           label = "Chen low β (CF=5)"),
  chen_high_cf5     = list(beta = beta_high,    cf = CF_5,    col = "pm25_mean",
                           label = "Chen high β (CF=5)"),
  
  # No-threshold (CF=0): "total burden above background"
  chen_nothreshold  = list(beta = beta_central, cf = CF_NONE, col = "pm25_mean",
                           label = "No-threshold linear (CF=0)"),
  
  # Elser corrected — EXPLORATORY, non-significant
  elser_cf15        = list(beta = beta_elser,   cf = CF_WHO,  col = "pm25_mean",
                           label = "Elser corrected (CF=15)*"),
  elser_cf5         = list(beta = beta_elser,   cf = CF_5,    col = "pm25_mean",
                           label = "Elser corrected (CF=5)*"),
  # High cohort HR
  highcohort_cf15   = list(beta = beta_highcohort, cf = CF_WHO, col = "pm25_mean",
                           label = "High cohort HR ~1.09/10µg/m³ (CF=15)"),
  
  # Afternoon exposure sensitivity (requires re-run of rasterization)
  afternoon_cf15    = list(beta = beta_central, cf = CF_WHO,  col = "pm25_mean_afternoon",
                           label = "Afternoon exposure 12-18h PDT (CF=15)"),
  afternoon_cf5     = list(beta = beta_central, cf = CF_5,    col = "pm25_mean_afternoon",
                           label = "Afternoon exposure 12-18h PDT (CF=5)")
)

all_scenarios_results <- map_dfr(years, function(yr) {
  cf_emp <- cf_central_emp[[as.character(yr)]]
  
  annual <- st_read(paste0("da_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE)
  
  # Cache both exposure columns for this year to avoid re-reading
  fs_24hr      <- get_fireseason(yr, PM25_DAILY_CAP, "pm25_mean")
  fs_afternoon <- tryCatch(
    get_fireseason(yr, PM25_DAILY_CAP, "pm25_mean_afternoon"),
    error = function(e) {
      message("  NOTE: pm25_mean_afternoon not found for ", yr,
              " — afternoon scenarios skipped. Re-run rasterization to enable.")
      NULL
    }
  )
  
  fixed_scenarios <- list(
    chen_central_cf_emp = list(beta = beta_central, cf = cf_emp,  
                               label = paste0("Chen central (CF=", round(cf_emp,1), " empirical)")),
    chen_low_cf_emp     = list(beta = beta_low,     cf = cf_emp,  
                               label = "Chen low β (empirical CF)"),
    chen_high_cf_emp    = list(beta = beta_high,    cf = cf_emp,  
                               label = "Chen high β (empirical CF)"),
    chen_central_cf15   = list(beta = beta_central, cf = CF_WHO,  
                               label = "Chen central (CF=15 WHO AQG)"),
    chen_low_cf15       = list(beta = beta_low,     cf = CF_WHO,  
                               label = "Chen low β (CF=15)"),
    chen_high_cf15      = list(beta = beta_high,    cf = CF_WHO,  
                               label = "Chen high β (CF=15)"),
    chen_central_cf5    = list(beta = beta_central, cf = CF_5,    
                               label = "Chen central (CF=5)"),
    chen_nothreshold    = list(beta = beta_central, cf = CF_NONE, 
                               label = "No-threshold linear (CF=0)"),
    elser_cf_emp        = list(beta = beta_elser,   cf = cf_emp,  
                               label = "Elser corrected (empirical CF)*"),
    highcohort_cf_emp   = list(beta = beta_highcohort, cf = cf_emp, 
                               label = "High cohort HR ~1.09/10µg/m³ (empirical CF)")
  )
  
  map_dfr(names(scenarios), function(sc_name) {
    sc <- scenarios[[sc_name]]
    
    
    fs <- if (sc$col == "pm25_mean_afternoon") {
      if (is.null(fs_afternoon)) return(NULL)   
      fs_afternoon
    } else {
      fs_24hr
    }
    
    annual %>%
      left_join(fs, by = "DGUID") %>%
      calculate_crf(sc$beta, sc$cf, primary_rates) %>%
      mutate(year = yr, scenario = sc_name, scenario_label = sc$label)
  })
})

# Generate empirical CF scenarios from existing all_scenarios_results
emp_cf_scenarios <- map_dfr(years, function(yr) {
  cf_emp <- cf_central_emp[[as.character(yr)]]
  
  annual <- st_read(paste0("da_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE)
  fs     <- get_fireseason(yr, PM25_DAILY_CAP, "pm25_mean")
  
  bind_rows(
    annual %>% left_join(fs, by = "DGUID") %>%
      calculate_crf(beta_central, cf_emp, primary_rates) %>%
      mutate(year = yr, scenario = "chen_central_cf_emp",
             scenario_label = paste0("Chen central (CF=", round(cf_emp,1), " empirical)")),
    annual %>% left_join(fs, by = "DGUID") %>%
      calculate_crf(beta_low, cf_emp, primary_rates) %>%
      mutate(year = yr, scenario = "chen_low_cf_emp",
             scenario_label = "Chen low β (empirical CF)"),
    annual %>% left_join(fs, by = "DGUID") %>%
      calculate_crf(beta_high, cf_emp, primary_rates) %>%
      mutate(year = yr, scenario = "chen_high_cf_emp",
             scenario_label = "Chen high β (empirical CF)")
  )
})

# Summarise and append to existing sensitivity_summary
emp_cf_summary <- emp_cf_scenarios %>%
  group_by(year, scenario, scenario_label) %>%
  summarise(
    attributable_cases = round(sum(attrib_total, na.rm = TRUE), 1),
    .groups = "drop"
  )

sensitivity_summary <- bind_rows(sensitivity_summary, emp_cf_summary)

# Verify
sensitivity_summary %>% distinct(scenario) %>% print(n = 30)

# Summary
# --- Year-by-year primary results ---
summary_by_year <- all_years_results %>%
  group_by(year) %>%
  summarise(
    pop_65plus         = sum(pop_total_65plus,  na.rm = TRUE),
    expected_cases     = sum(expected_total,    na.rm = TRUE),
    attributable_cases = sum(attrib_total,      na.rm = TRUE),
    pct_attributable   = attributable_cases / expected_cases * 100,
    mean_pm25_wtd      = weighted.mean(fireseason_mean_pm25,
                                       pop_total_65plus, na.rm = TRUE),
    n_das_exposed      = sum(delta_pm25 > 0, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_by_year %>% mutate(across(where(is.numeric), ~round(., 1))))

# --- All scenarios summary ---
sensitivity_summary <- all_scenarios_results %>%
  group_by(year, scenario, scenario_label) %>%
  summarise(
    attributable_cases = round(sum(attrib_total, na.rm = TRUE), 1),
    .groups = "drop"
  )

# --- Dual counterfactual comparison table ---
cf_comparison <- sensitivity_summary %>%
  filter(scenario %in% c("chen_central_cf15", "chen_central_cf5")) %>%
  pivot_wider(
    id_cols     = year,
    names_from  = scenario,
    values_from = attributable_cases
  ) %>%
  rename(cases_cf15 = chen_central_cf15, cases_cf5 = chen_central_cf5) %>%
  mutate(
    additional_cases = round(cases_cf5 - cases_cf15, 1),
    pct_increase     = round((cases_cf5 - cases_cf15) / cases_cf15 * 100, 1)
  )

print(cf_comparison)

# Add absolute context to the comparison
cf_comparison %>%
  mutate(
    pct_increase_label = case_when(
      pct_increase > 500 ~ paste0(round(pct_increase, 0), "% — low base effect"),
      TRUE               ~ paste0(round(pct_increase, 1), "%")
    ),
    interpretation = case_when(
      cases_cf15 < 5  ~ "CF=15 near zero; pct increase uninformative",
      pct_increase > 200 ~ "Large relative increase; moderate absolute increase",
      TRUE            ~ "Meaningful increase in both absolute and relative terms"
    )
  ) %>%
  print()

# --- Afternoon vs 24-hr comparison ---
afternoon_comparison <- sensitivity_summary %>%
  filter(scenario %in% c("chen_central_cf15", "afternoon_cf15")) %>%
  pivot_wider(
    id_cols     = year,
    names_from  = scenario,
    values_from = attributable_cases
  ) %>%
  rename(
    cases_24hr      = chen_central_cf15,
    cases_afternoon = afternoon_cf15
  ) %>%
  mutate(
    pct_difference = round((cases_afternoon - cases_24hr) / cases_24hr * 100, 1)
  )

print(afternoon_comparison)

# --- Full scenario cumulative table ---
sensitivity_summary %>%
  group_by(scenario, scenario_label) %>%
  summarise(cumulative = round(sum(attributable_cases), 1), .groups = "drop") %>%
  arrange(cumulative) %>%
  mutate(note = if_else(str_detect(scenario, "elser"),
                        "* Non-significant after correction (OR 0.98-1.28)",
                        "")) %>%
  print(n = 20)

# --- Wide uncertainty table (primary + beta bounds + counterfactuals) ---
uncertainty_wide <- sensitivity_summary %>%
  filter(scenario %in% c("chen_low_cf_emp", "chen_central_cf_emp",
                         "chen_high_cf_emp", "chen_central_cf15",
                         "chen_central_cf5")) %>%
  pivot_wider(
    id_cols     = year,
    names_from  = scenario,
    values_from = attributable_cases
  ) %>%
  rename(
    low_emp     = chen_low_cf_emp,
    central_emp = chen_central_cf_emp,
    high_emp    = chen_high_cf_emp,
    central_cf15 = chen_central_cf15,
    central_cf5  = chen_central_cf5
  )

print(uncertainty_wide)

# exposure summary by year
exposure_summary <- map_dfr(years, function(yr) {
  fs <- get_fireseason(yr, PM25_DAILY_CAP, "pm25_mean")
  
  annual <- st_read(paste0("da_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE) %>%
    st_drop_geometry() %>%
    select(DGUID, DAUID, pop_total_65plus)
  
  joined <- annual %>% left_join(fs, by = "DGUID")
  
  total_das          <- nrow(joined)
  exposed_das        <- sum(joined$fireseason_mean_pm25 > 15, na.rm = TRUE)
  total_pop_65plus   <- sum(joined$pop_total_65plus, na.rm = TRUE)
  exposed_pop_65plus <- sum(
    joined$pop_total_65plus[joined$fireseason_mean_pm25 > 15], na.rm = TRUE
  )
  
  tibble(
    year               = yr,
    total_das          = total_das,
    exposed_das        = exposed_das,
    pct_das_exposed    = round(exposed_das / total_das * 100, 1),
    total_pop_65plus   = total_pop_65plus,
    exposed_pop_65plus = exposed_pop_65plus,
    pct_pop_exposed    = round(exposed_pop_65plus / total_pop_65plus * 100, 1)
  )
})

print(exposure_summary)

# REGIONAL BREAKDOWN 
region_lookup <- function(cd_code) {
  case_when(
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
} 

regional_by_year <- all_years_results %>%
  filter(delta_pm25 > 0) %>%
  mutate(
    cd_code = substr(DAUID, 1, 4),
    region  = region_lookup(cd_code)
  ) %>%
  group_by(year, region) %>%
  summarise(
    n_das              = n(),
    pop_65plus         = sum(pop_total_65plus, na.rm = TRUE),
    mean_pm25          = weighted.mean(fireseason_mean_pm25,
                                       pop_total_65plus, na.rm = TRUE),
    attributable_cases = sum(attrib_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year, desc(attributable_cases))

regional_cumulative <- regional_by_year %>%
  group_by(region) %>%
  summarise(
    total_attributable = sum(attributable_cases, na.rm = TRUE),
    mean_annual        = mean(attributable_cases, na.rm = TRUE),
    n_years_exposed    = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(total_attributable))

regional_by_year %>%
  group_by(year) %>%
  slice_max(attributable_cases, n = 5) %>%
  print(n = 30)

print(regional_cumulative, n = 20)


# ------------------------------------------------------------
#  Plot 1: Attribuable cases by year 
plot_data_main <- uncertainty_wide %>%
  mutate(
    is_small            = central_emp < 50,
    ci_label_y = ifelse(is_small, high_emp + 25, high_emp + 8)
  )

# Pull afternoon estimates for overlay
afternoon_line <- sensitivity_summary %>%
  filter(scenario == "afternoon_cf15") %>%
  select(year, afternoon = attributable_cases)

plot_data_main <- plot_data_main %>%
  left_join(afternoon_line, by = "year")

# primary atttribuable cases 
p1 <- ggplot(plot_data_main, aes(x = factor(year))) +
  
  # Beta CI shading
  geom_rect(
    aes(
      xmin = as.numeric(factor(year)) - 0.35,
      xmax = as.numeric(factor(year)) + 0.35,
      ymin = low_cf15,
      ymax = high_cf15
    ),
    fill = "#ca706a", alpha = 0.25
  ) +
  
  # Primary bar
  geom_col(aes(y = central_cf15),
           fill = "#ca706a", alpha = 0.85, width = 0.5) +
  
  # Afternoon overlay only — stays close to primary scale
  geom_point(aes(y = afternoon,
                 shape = "Afternoon exposure (12\u201318h PDT)"),
             color = "#457B9D", size = 3.5, na.rm = TRUE) +
  geom_line(aes(y = afternoon, group = 1),
            color = "#457B9D", linewidth = 0.7,
            linetype = "dashed", na.rm = TRUE) +
  
  # Value labels
  geom_text(aes(y = central_cf15, label = round(central_cf15, 0)),
            vjust = -0.3, size = 3.3, fontface = "bold",
            color = "#1a1f2e") +
  geom_text(aes(y = ci_label_y,
                label = paste0("(", round(low_cf15, 0),
                               "\u2013", round(high_cf15, 0), ")")),
            vjust = 0, size = 2.6, color = "#1a1f2e") +
  
  geom_segment(
    data = plot_data_main %>% filter(is_small),
    aes(
      x    = as.numeric(factor(year)),
      xend = as.numeric(factor(year)),
      y    = high_cf15 + 3,
      yend = high_cf15 + 20
    ),
    color = "grey60", linewidth = 0.3, linetype = "dotted"
  ) +
  
  scale_shape_manual(
    name   = "Sensitivity",
    values = c("Afternoon exposure (12\u201318h PDT)" = 18)
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.28)),
    labels = scales::comma
  ) +
  
  labs(
    title    = "Estimated Alzheimer's Cases Attributable to Wildfire PM~2.5~",
    subtitle = paste0(
      "British Columbia, 2021\u20132025 Fire Seasons (May\u2013September)",
      " \u00b7 Population Aged 65+ \u00b7 Empirical CF (BC MoE non-fire-season median)"
    ),
    x       = "Year",
    y       = "Attributable Alzheimer's cases",
    caption = paste0(
      "Bars: Chen et al. (2017) central estimate, empirical counterfactual",
      " ((BC MoE monitoring network non-fire-season median, range 4.8–5.3 µg/m³),.\n",
      "Shaded range: β sensitivity (HR 1.03–1.05 per 4.8 µg/m³ IQR). ",
      "Daily cap: 150 µg/m³. CCDSS age-specific rates (65+)"
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_markdown(face = "bold", size = 18, hjust = 0),
    plot.subtitle      = element_text(color = "grey40", size = 12, hjust = 0),
    plot.caption       = element_text(color = "grey50", size = 8.5, hjust = 0,
                                      lineheight = 1.4),
    legend.position    = c(0.80, 0.85),
    legend.background  = element_rect(fill = "white", color = "grey90"),
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "grey30"),
    axis.title         = element_text(color = "grey30", size = 10),
    plot.background    = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

print(p1)

ggsave("attributable_cases_primary.png",
       plot = p1, width = 11, height = 7,
       dpi = 300, bg = "#ffffffC9")

# plot 1.5 CF=5 comparison
cf5_plot_data <- cf_comparison %>%
  select(year, `WHO AQG (CF=15)` = cases_cf15,
         `Near-background (CF=5)` = cases_cf5) %>%
  pivot_longer(-year, names_to = "counterfactual", values_to = "cases") %>%
  mutate(
    counterfactual = factor(counterfactual,
                            levels = c("WHO AQG (CF=15)",
                                       "Near-background (CF=5)"))
  )

p1b <- ggplot(cf5_plot_data,
              aes(x = factor(year), y = cases, fill = counterfactual)) +
  
  geom_col(position = "dodge", width = 0.6, alpha = 0.88) +
  
  geom_text(
    aes(label = if_else(cases < 1, "<1", as.character(round(cases, 0)))),
    position = position_dodge(width = 0.6),
    vjust = -0.4, size = 3, color = "grey20"
  ) +
  
  scale_fill_manual(
    name   = "Counterfactual",
    values = c(
      "WHO AQG (CF=15)"        = "#ca706a",
      "Near-background (CF=5)" = "#7a2c2c"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.2)),
    labels = scales::comma
  ) +
  
  labs(
    title    = "Counterfactual Comparison: WHO AQG vs Near-Background",
    subtitle = paste0(
      "British Columbia, 2021\u20132025 \u00b7 Chen et al. (2017) CRF \u00b7",
      " Population Aged 65+"
    ),
    x       = "Year",
    y       = "Attributable Alzheimer's cases",
    caption = paste0(
      "CF=15: cases attributable to smoke above WHO 24-hr AQG (15 \u00b5g/m\u00b3). ",
      "CF=5: cases attributable to smoke above near-background (5 \u00b5g/m\u00b3).\n",
      "Difference represents burden in populations exposed below the WHO guideline. ",
      "Daily cap: 150 \u00b5g/m\u00b3. CCDSS age-specific rates (65+)."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle      = element_text(color = "grey40", size = 11, hjust = 0),
    plot.caption       = element_text(color = "grey50", size = 8.5, hjust = 0,
                                      lineheight = 1.4),
    legend.position    = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "grey30"),
    axis.title         = element_text(color = "grey30", size = 10),
    plot.background    = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )
print(p1b)

ggsave("attributable_cases_counterfactual_comparison.png",
       plot = p1b, width = 10, height = 6,
       dpi = 300, bg = "#ffffffC9")

# PLOT 2: CRF comparison - all scenarios 
cf_emp_min <- min(cf_annual_emp$cf_median)  # 4.8
cf_emp_max <- max(cf_annual_emp$cf_median)  # 5.3
cf_emp_mid <- mean(cf_annual_emp$cf_median) # ~5.1

pm25_range <- seq(0, 100, by = 0.5)

crf_curves_multi <- bind_rows(
  tibble(pm25 = pm25_range, beta = beta_central, cf = cf_emp_mid,
         crf = "Chen et al. (2017), empirical CF [primary]"),
  tibble(pm25 = pm25_range, beta = beta_low,     cf = cf_emp_mid,
         crf = "Chen low \u03b2, empirical CF"),
  tibble(pm25 = pm25_range, beta = beta_high,    cf = cf_emp_mid,
         crf = "Chen high \u03b2, empirical CF"),
  tibble(pm25 = pm25_range, beta = beta_central, cf = CF_WHO,
         crf = "Chen et al. (2017), CF=15 (sensitivity)"),
  tibble(pm25 = pm25_range, beta = beta_highcohort, cf = cf_emp_mid,
         crf = "High cohort HR ~1.09/10\u00b5g/m\u00b3, empirical CF"),
  tibble(pm25 = pm25_range, beta = beta_central, cf = CF_NONE,
         crf = "No-threshold (CF=0)")
) %>%
  mutate(
    delta_pm25 = pmax(pm25 - cf, 0),
    PAF        = (exp(beta * delta_pm25) - 1) / exp(beta * delta_pm25)
  )

p_crf <- ggplot(crf_curves_multi,
                aes(x = pm25, y = PAF, color = crf, linetype = crf)) +
  
  annotate("rect",
           xmin = cf_emp_min, xmax = cf_emp_max,
           ymin = 0, ymax = 0.85,
           fill = "#2d6a4f", alpha = 0.08) +
  geom_line(linewidth = 1.0) +
  geom_vline(xintercept = cf_emp_mid,
             linetype = "dashed", color = "#2d6a4f", linewidth = 0.7) +
  geom_vline(xintercept = CF_WHO,
             linetype = "dotted", color = "grey50", linewidth = 0.6) +
  
  annotate("text", x = cf_emp_mid + 0.4, y = 0.78,
           label = paste0("Empirical CF\n(", round(cf_emp_min,1),
                          "\u2013", round(cf_emp_max,1), " \u00b5g/m\u00b3)"),
           size = 2.8, hjust = 0, color = "#2d6a4f", lineheight = 1.2) +
  annotate("text", x = CF_WHO + 0.4, y = 0.68,
           label = "CF=15\n(sensitivity)",
           size = 2.8, hjust = 0, color = "grey50", lineheight = 1.2) +
  
  scale_color_manual(name = NULL, values = c(
    "Chen et al. (2017), empirical CF [primary]"     = "#ca706a",
    "Chen low \u03b2, empirical CF"                  = "#e8a89e",
    "Chen high \u03b2, empirical CF"                 = "#9e3a2e",
    "Chen et al. (2017), CF=15 (sensitivity)"        = "#ca706a",
    "High cohort HR ~1.09/10\u00b5g/m\u00b3, empirical CF" = "#457B9D",
    "No-threshold (CF=0)"                            = "#E9C46A"
  )) +
  scale_linetype_manual(name = NULL, values = c(
    "Chen et al. (2017), empirical CF [primary]"     = "solid",
    "Chen low \u03b2, empirical CF"                  = "dotted",
    "Chen high \u03b2, empirical CF"                 = "dotted",
    "Chen et al. (2017), CF=15 (sensitivity)"        = "longdash",
    "High cohort HR ~1.09/10\u00b5g/m\u00b3, empirical CF" = "longdash",
    "No-threshold (CF=0)"                            = "dotdash"
  )) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = c(0, cf_emp_mid, 15, 25, 50, 75, 100),
                     labels = c("0", paste0(round(cf_emp_mid,1), "\nCF"), 
                                "15", "25", "50", "75", "100")) +
  scale_y_continuous(limits = c(0, 0.85),
                     labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  
  labs(
    title    = "Concentration\u2013Response Function",
    subtitle = "Population attributable fraction (PAF) by fire-season PM~2.5~",
    x        = "Fire-season mean PM~2.5~ (\u00b5g/m\u00b3)",
    y        = "Population attributable fraction",
    caption  = paste0(
      "Shaded band: 95% uncertainty range across \u03b2 sensitivity bounds ",
      "(HR 1.03\u20131.05 per 4.8 \u00b5g/m\u00b3 IQR).\n",
      "Green band: empirical CF range across years (", round(cf_emp_min,1),
      "\u2013", round(cf_emp_max,1), " \u00b5g/m\u00b3, BC MoE non-fire-season median). ",
      "CF=15 shown as sensitivity only.\n",
      "High cohort HR: upper range of non-wildfire ambient PM\u2082.\u2085\u2013dementia literature. ",
      "No-threshold: Chen \u03b2 from 0 \u00b5g/m\u00b3."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "right",
    legend.text      = element_markdown(size = 9),
    plot.title       = element_markdown(face = "bold", size = 16, hjust = 0),
    plot.subtitle    = element_markdown(color = "grey40", size = 11, hjust = 0,
                                    margin = margin(b = 10)),
    axis.title.x     = element_markdown(color = "grey30", size = 11),
    axis.title.y     = element_markdown(color = "grey30", size = 11),
    plot.caption     = element_markdown(color = "grey50", size = 8.5, hjust = 0,
                                    lineheight = 1.4, margin = margin(t = 10)),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin      = margin(16, 16, 16, 16)
  )

print(p_crf)
ggsave("crf_comparison_scenarios.png",
       plot = p_crf, width = 9, height = 6, dpi = 300, bg = "#ffffffC9")
    
# Plot 3: Dual counterfactual 

plot_cf_data <- sensitivity_summary %>%
  filter(scenario %in% c(
    "chen_central_cf15", "chen_low_cf15", "chen_high_cf15",
    "chen_central_cf5",  "chen_low_cf5",  "chen_high_cf5"
  )) %>%
  mutate(
    cf    = if_else(str_detect(scenario, "cf15"),
                    "WHO 24-hr AQG counterfactual (15 \u00b5g/m\u00b3)",
                    "WHO annual mean AQG (5 \u00b5g/m\u00b3)"),
    cf    = factor(cf, levels = c("WHO 24-hr AQG counterfactual (15 \u00b5g/m\u00b3)",
                                  "WHO annual mean AQG (5 \u00b5g/m\u00b3)")),
    bound = case_when(
      scenario %in% c("chen_low_cf15",     "chen_low_cf5")     ~ "low",
      scenario %in% c("chen_high_cf15",    "chen_high_cf5")    ~ "high",
      scenario %in% c("chen_central_cf15", "chen_central_cf5") ~ "central",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(bound)) %>%
  select(year, cf, bound, attributable_cases) %>%
  group_by(year, cf, bound) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_wider(names_from = bound, values_from = attributable_cases) %>%
  mutate(
    central = as.numeric(central),
    low     = as.numeric(low),
    high    = as.numeric(high)
  )

# Verify
cat("Rows in plot_cf_data:", nrow(plot_cf_data), "\n")  # should be 10 (5 years × 2 CFs)
cat("Any NAs?\n")
plot_cf_data %>%
  summarise(across(c(central, low, high), ~sum(is.na(.)))) %>%
  print()
glimpse(plot_cf_data)

p_cf <- ggplot(plot_cf_data, aes(x = factor(year))) +
  
  geom_rect(
    aes(xmin = as.numeric(factor(year)) - 0.35,
        xmax = as.numeric(factor(year)) + 0.35,
        ymin = low, ymax = high),
    fill = "#ca706a", alpha = 0.25
  ) +
  geom_col(aes(y = central),
           fill = "#ca706a", alpha = 0.85, width = 0.5) +
  geom_text(aes(y = central, label = round(central, 0)),
            vjust = -0.3, size = 3.2, fontface = "bold", color = "#1a1f2e") +
  geom_text(aes(y = high * 1.08,
                label = paste0("(", round(low, 0), "–", round(high, 0), ")")),
            vjust = 0, size = 2.6, color = "#1a1f2e") +
  
  facet_wrap(~ cf, ncol = 1, scales = "free_y") +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.3)),
                     labels = comma) +
  
  labs(
    title    = "Attributable Alzheimer's Cases: Two Counterfactual Scenarios",
    subtitle = "British Columbia, 2021\u20132025 \u00b7 Chen et al. (2017) CRF \u00b7 Population 65+",
    x        = "Year",
    y        = "Attributable cases",
    caption  = paste0(
      "Top: WHO 24-hr AQG (15 \u00b5g/m\u00b3) — burden above health guideline. ",
      "Bottom: WHO annual mean AQG (5 \u00b5g/m\u00b3) — representing burden near-background.\n",
      "Shaded band: \u03b2 sensitivity (HR 1.03\u20131.05). Daily cap 150 \u00b5g/m\u00b3. CCDSS rates."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle      = element_text(color = "grey40", size = 11, hjust = 0),
    plot.caption       = element_text(color = "grey50", size = 8.5, hjust = 0,
                                      lineheight = 1.4),
    strip.text         = element_text(face = "bold", size = 10, hjust = 0),
    strip.background   = element_rect(fill = "grey95", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.background    = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

print(p_cf)
ggsave("attributable_cases_dual_counterfactual.png",
       plot = p_cf, width = 10, height = 9, dpi = 300, bg = "#ffffffC9")


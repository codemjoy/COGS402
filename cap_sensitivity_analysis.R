library(tidyverse)
library(sf)
library(purrr)

years             <- 2021:2025
PM25_DAILY_CAP    <- 150          
CF_PRIMARY            <- 15           
cap_thresholds    <- c(100, 150, 200, 250, 500)   
fire_season_start <- "-05-01"
fire_season_end   <- "-10-31"

cra_cases_chen <- function(pm25_vec, pop_65_74, pop_75_84, pop_85plus,
                           beta, inc_65_74, inc_75_84, inc_85plus, cf) {
  rr_obs <- exp(beta * pm25_vec)
  rr_cf  <- exp(beta * cf)
  paf    <- pmax((rr_obs - rr_cf) / rr_obs, 0)
  sum(paf * (pop_65_74 * inc_65_74 +
               pop_75_84 * inc_75_84 +
               pop_85plus * inc_85plus),
      na.rm = TRUE)
}

build_da_cap <- function(yr, cap) {
  
  # Daily exposure
  fireseason <- readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds")) %>%
    filter(
      date >= as.Date(paste0(yr, fire_season_start)),
      date <= as.Date(paste0(yr, fire_season_end))
    ) %>%
    mutate(pm25_capped = pmin(pm25_mean, cap)) %>%
    group_by(DGUID) %>%
    summarise(pm25_mean = mean(pm25_capped, na.rm = TRUE), .groups = "drop")
  
  annual <- st_read(paste0("da_annual_pm25_summary_", yr, ".gpkg"), quiet = TRUE) %>%
    st_drop_geometry() %>%
    select(DGUID, pop_total_65_74, pop_total_75_84, pop_total_85plus, pop_total_65plus)
  
  annual %>%
    left_join(fireseason, by = "DGUID") %>%
    rename(
      pop_65_74  = pop_total_65_74,
      pop_75_84  = pop_total_75_84,
      pop_85plus = pop_total_85plus,
      pop_65plus = pop_total_65plus
    )
}

get_fireseason_raw <- function(yr) {
  readRDS(paste0("da_daily_pm25_exposure_", yr, ".rds")) %>%
    filter(
      date >= as.Date(paste0(yr, fire_season_start)),
      date <= as.Date(paste0(yr, fire_season_end))
    )
}

# attributable cases by year/cap threshold
cap_results <- map_dfr(cap_thresholds, function(cap) {
  map_dfr(years, function(yr) {
    d <- build_da_cap(yr, cap)
    tibble(
      year         = yr,
      cap          = paste0(cap, " µg/m³"),
      attrib_cases = round(
        cra_cases_chen(
          d$pm25_mean, d$pop_65_74, d$pop_75_84, d$pop_85plus,
          beta       = CHEN_BETA_CENTRAL,
          inc_65_74  = INCIDENCE_65_74,
          inc_75_84  = INCIDENCE_75_84,
          inc_85plus = INCIDENCE_85PLUS,
          cf         = CF_PRIMARY
        ), 1)
    )
  })
})

cap_wide <- cap_results %>%
  mutate(cap = factor(cap, levels = paste0(cap_thresholds, " µg/m³"))) %>%
  pivot_wider(names_from = cap, values_from = attrib_cases) %>%
  arrange(year)

cap_wide_with_total <- cap_wide %>%
  mutate(year = as.character(year)) %>%
  bind_rows(
    cap_wide %>%
      summarise(across(where(is.numeric), ~round(sum(.), 1))) %>%
      mutate(year = "Cumulative")
  )

print(cap_wide_with_total)

# mean fire season pm2.5 by year and cap threshold
pm25_sensitivity <- map_dfr(cap_thresholds, function(cap) {
  map_dfr(years, function(yr) {
    d <- build_da_cap(yr, cap) %>%
      filter(!is.na(pm25_mean), pop_65plus > 0)
    
    wtd_mean   <- sum(d$pm25_mean * d$pop_65plus) / sum(d$pop_65plus)
    unwtd_mean <- mean(d$pm25_mean, na.rm = TRUE)
    
    tibble(
      year            = yr,
      cap             = paste0(cap, " µg/m³"),
      mean_pm25_wtd   = round(wtd_mean,   3),
      mean_pm25_unwtd = round(unwtd_mean, 3)
    )
  })
})

cap_levels <- paste0(cap_thresholds, " µg/m³")

pm25_sens_wide <- pm25_sensitivity %>%
  mutate(cap = factor(cap, levels = cap_levels)) %>%
  select(year, cap, mean_pm25_wtd) %>%
  pivot_wider(names_from = cap, values_from = mean_pm25_wtd) %>%
  arrange(year)

print(pm25_sens_wide)

pm25_sens_unwtd <- pm25_sensitivity %>%
  mutate(cap = factor(cap, levels = cap_levels)) %>%
  select(year, cap, mean_pm25_unwtd) %>%
  pivot_wider(names_from = cap, values_from = mean_pm25_unwtd) %>%
  arrange(year)

print(pm25_sens_unwtd)

# frequency of daily DA exceedance
exceedance_thresholds <- cap_thresholds 

exceedance_table <- map_dfr(years, function(yr) {
  
  raw <- get_fireseason_raw(yr)
  
  total_da_days <- nrow(raw)
  total_das     <- n_distinct(raw$DGUID)
  
  thresh_stats <- map_dfr(exceedance_thresholds, function(thresh) {
    exceed <- raw %>% filter(pm25_mean > thresh)
    tibble(
      threshold_ug_m3  = thresh,
      da_days_exceeded = nrow(exceed),
      das_affected     = n_distinct(exceed$DGUID),
      pct_da_days      = round(nrow(exceed) / total_da_days * 100, 3),
      pct_das          = round(n_distinct(exceed$DGUID) / total_das * 100, 2),
      max_raw_pm25     = round(max(raw$pm25_mean, na.rm = TRUE), 1)
    )
  })
  
  thresh_stats %>% mutate(year = yr, total_da_days = total_da_days, total_das = total_das)
})

exceedance_wide <- exceedance_table %>%
  select(year, threshold_ug_m3, da_days_exceeded, pct_da_days) %>%
  pivot_wider(
    names_from  = threshold_ug_m3,
    values_from = c(da_days_exceeded, pct_da_days),
    names_glue  = "{.value}_{threshold_ug_m3}"
  ) %>%
  arrange(year)

cat("DA-days exceeding each threshold (count | % of fire-season DA-days):\n")
print(exceedance_wide)

exceedance_table %>%
  select(year, threshold_ug_m3, da_days_exceeded, pct_da_days,
         das_affected, pct_das, max_raw_pm25, total_da_days) %>%
  arrange(year, threshold_ug_m3) %>%
  print(n = Inf)


#DA-days affected by 150 cap
cap_audit <- map_dfr(years, function(yr) {
  raw <- get_fireseason_raw(yr)
  
  raw %>%
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

write_csv(cap_wide_with_total,  "cap_sensitivity_attributable_cases.csv")
write_csv(pm25_sens_wide,       "cap_sensitivity_mean_pm25_weighted.csv")
write_csv(pm25_sens_unwtd,      "cap_sensitivity_mean_pm25_unweighted.csv")
write_csv(exceedance_table %>%
            select(year, threshold_ug_m3, da_days_exceeded, pct_da_days,
                   das_affected, pct_das, max_raw_pm25, total_da_days) %>%
            arrange(year, threshold_ug_m3),
          "cap_sensitivity_exceedance_supplemental.csv")
write_csv(cap_audit,            "cap_sensitivity_audit.csv")


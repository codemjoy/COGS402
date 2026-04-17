library(tidyverse)
library(sf)
library(terra)
library(ggtext)

set.seed(42)

YEARS          <- 2021:2025
FIRE_MONTHS    <- 5:10
ARTIFACT_CAP   <- 250
N_SIM          <- 10000

# GEMM PARAMETERS 
GEMM_THETA_CENTRAL <- 0.12
GEMM_THETA_SE      <- 0.03
GEMM_MU            <- -1.3      
GEMM_TAU           <- 0.1       
GEMM_R             <- 143.3     
GEMM_XCF_LO        <- 2.7       
GEMM_XCF_HI        <- 7.6       


gemm_rr <- function(z, theta, mu, tau, r, xcf) {
  x     <- pmax(z - xcf, 0)
  f_x   <- log(x + 1)
  omega <- 1 / (1 + exp(-(x - mu) / (tau * r)))
  exp(theta * f_x * omega)
}

LAG_WINDOWS <- c(10, 15, 20)

cf_regional <- read_csv("cf_annual_regional.csv") %>%
  filter(fire_year %in% YEARS) %>%
  mutate(
    fire_year = as.integer(fire_year),
    cf_se     = (cf_p75 - cf_p25) / 2
  )

da_census <- st_read("bc_das_population_for_raster.gpkg")

da_census <- da_census %>%
  mutate(
    pop_45_64 = female45.to.49.years + female50.to.54.years +
      female55.to.59.years + female60.to.64.years +
      male45.to.49.years   + male50.to.54.years   +
      male55.to.59.years   + male60.to.64.years,
    
    pop_65_74 = female65.to.69.years + female70.to.74.years +
      male65.to.69.years   + male70.to.74.years,
    
    pop_75_84 = female75.to.79.years + female80.to.84.years +
      male75.to.79.years   + male80.to.84.years,
    
    pop_85plus = female85.to.89.years + female90.to.94.years +
      female95.to.99.years + female100.years.and.over +
      male85.to.89.years   + male90.to.94.years   +
      male95.to.99.years   + male100.years.and.over,
    
    pop_65plus    = pop_65_74 + pop_75_84 + pop_85plus,
    cohort_midage = 55
  )

# fire exposure
exposure_list <- lapply(YEARS, function(yr) {
  
  f <- paste0("da_daily_pm25_exposure_", yr, ".rds")
  if (!file.exists(f)) {
    warning("File not found: ", f, " — skipping year ", yr)
    return(NULL)
  }
  
  daily <- readRDS(f)
  
  daily %>%
    filter(month(date) %in% FIRE_MONTHS) %>%
    mutate(pm25_capped = pmin(pm25_mean, ARTIFACT_CAP)) %>%
    group_by(DGUID) %>%
    summarise(
      pm25_fireseason = mean(pm25_capped, na.rm = TRUE),
      n_days          = n(),
      .groups         = "drop"
    ) %>%
    mutate(year = yr)
})

exposure_all <- bind_rows(exposure_list)

exposure_with_pop <- exposure_all %>%
  left_join(
    st_drop_geometry(da_census) %>%
      select(DGUID, pop_45_64, pop_65_74, pop_75_84, pop_85plus, pop_65plus),
    by = "DGUID"
  ) %>%
  filter(!is.na(pop_45_64), pop_45_64 > 0) %>%
  mutate(
    cd_code  = substr(DGUID, 10, 13),
    air_zone = case_when(
      cd_code %in% c("5955", "5957", "5959", "5939",
                     "5943", "5945", "5947", "5951") ~ "Northeast",
      cd_code %in% c("5941", "5949", "5953")         ~ "Central Interior",
      cd_code %in% c("5901", "5903", "5905", "5907",
                     "5909", "5911")                 ~ "Coastal",
      cd_code %in% c("5933", "5935", "5937")         ~ "Southern Interior",
      cd_code %in% c("5913", "5915", "5917", "5919",
                     "5921", "5923")                 ~ "Georgia Strait",
      cd_code %in% c("5925", "5927", "5929", "5931") ~ "Lower Fraser Valley",
      TRUE ~ "Southern Interior"
    )
  )

# life table & survival 
life_raw <- read_csv(
  "Stat_Can_Life_Expentancy.csv",
  skip      = 8,
  col_names = c("age_label", "e_2018_2020", "e_2019_2021",
                "e_2020_2022", "e_2021_2023", "e_2022_2024"),
  col_types = cols(.default = "c"),
  n_max     = 111
)

life_table <- life_raw %>%
  filter(!is.na(age_label), age_label != "") %>%
  mutate(
    age_x = as.integer(str_extract(age_label, "^\\d+")),
    ex    = as.numeric(e_2021_2023)
  ) %>%
  filter(!is.na(age_x), !is.na(ex)) %>%
  select(age_x, ex)

compute_survival_from_ex <- function(from_age, lag, lt = life_table) {
  ex_from <- lt %>% filter(age_x == from_age)       %>% pull(ex)
  ex_to   <- lt %>% filter(age_x == from_age + lag) %>% pull(ex)
  if (length(ex_from) == 0 | length(ex_to) == 0) {
    stop("Age not found in life table: ", from_age, " or ", from_age + lag)
  }
  (ex_to + lag) / (ex_from + lag)
}

surv_lag10 <- compute_survival_from_ex(from_age = 55, lag = 10)
surv_lag15 <- compute_survival_from_ex(from_age = 55, lag = 15)
surv_lag20      <- compute_survival_from_ex(from_age = 55, lag = 20)

life_table_prelim <- life_raw %>%
  filter(!is.na(age_label), age_label != "") %>%
  mutate(
    age_x = as.integer(str_extract(age_label, "^\\d+")),
    ex    = as.numeric(e_2022_2024)
  ) %>%
  filter(!is.na(age_x), !is.na(ex)) %>%
  select(age_x, ex)

surv_lag10_sens <- compute_survival_from_ex(55, 10, lt = life_table_prelim)
surv_lag15_sens <- compute_survival_from_ex(55, 15, lt = life_table_prelim)
surv_lag20_sens <- compute_survival_from_ex(55, 20, lt = life_table_prelim)

message(sprintf("Survival lag-10 (55\u219265): %.4f  [sensitivity 2022-2024: %.4f]",
                surv_lag10, surv_lag10_sens))
message(sprintf("Survival lag-15 (55\u219270): %.4f  [sensitivity 2022-2024: %.4f]",
                surv_lag15, surv_lag15_sens))
message(sprintf("Survival lag-20 (55→75): %.4f  [sensitivity 2022-2024: %.4f]",
                surv_lag20, surv_lag20_sens))


# INCIDENCE RATES
INCIDENCE_65_74  <- 2.80 / 1000
INCIDENCE_75_84  <- 34.5 / 1000
INCIDENCE_85PLUS <- 44.6 / 1000

INCIDENCE_CV_YOUNG <- (617  - 604)  / (2 * 1.96 * 610)
INCIDENCE_CV_OLD   <- (3697 - 3640) / (2 * 1.96 * 3669)

INCIDENCE_65_74_SE  <- INCIDENCE_65_74  * INCIDENCE_CV_YOUNG
INCIDENCE_75_84_SE  <- INCIDENCE_75_84  * INCIDENCE_CV_OLD
INCIDENCE_85PLUS_SE <- INCIDENCE_85PLUS * INCIDENCE_CV_OLD

summarise_mc <- function(mc_results, cf_label) {
  
  per_year <- mc_results %>%
    group_by(year) %>%
    summarise(
      lag10_median = median(attr_lag10_total),
      lag10_lo95   = quantile(attr_lag10_total, 0.025),
      lag10_hi95   = quantile(attr_lag10_total, 0.975),
      lag15_median = median(attr_lag15_total),
      lag15_lo95   = quantile(attr_lag15_total, 0.025),
      lag15_hi95   = quantile(attr_lag15_total, 0.975),
      lag20_median = median(attr_lag20_total),           
      lag20_lo95   = quantile(attr_lag20_total, 0.025),  
      lag20_hi95   = quantile(attr_lag20_total, 0.975),  
      .groups      = "drop"
    ) %>%
    mutate(counterfactual = cf_label)
  
  cumulative <- mc_results %>%
    group_by(sim_id) %>%
    summarise(
      cumul_lag10 = sum(attr_lag10_total),
      cumul_lag15 = sum(attr_lag15_total),
      cumul_lag20 = sum(attr_lag20_total),              
      .groups     = "drop"
    ) %>%
    summarise(
      lag10_cumul_median = median(cumul_lag10),
      lag10_cumul_lo95   = quantile(cumul_lag10, 0.025),
      lag10_cumul_hi95   = quantile(cumul_lag10, 0.975),
      lag15_cumul_median = median(cumul_lag15),
      lag15_cumul_lo95   = quantile(cumul_lag15, 0.025),
      lag15_cumul_hi95   = quantile(cumul_lag15, 0.975),
      lag20_cumul_median = median(cumul_lag20),          
      lag20_cumul_lo95   = quantile(cumul_lag20, 0.025), 
      lag20_cumul_hi95   = quantile(cumul_lag20, 0.975) 
    ) %>%
    mutate(counterfactual = cf_label)
  
  list(per_year = per_year, cumulative = cumulative)
}

#monte carlo
run_lag_simulation_gemm <- function(cf_mode = c("empirical", "fixed"),
                                    fixed_cf = 5) {
  cf_mode <- match.arg(cf_mode)
  
  results <- lapply(seq_len(N_SIM), function(i) {
    

    theta_i <- max(rnorm(1, GEMM_THETA_CENTRAL, GEMM_THETA_SE), 0.01)
    xcf_i   <- runif(1, GEMM_XCF_LO, GEMM_XCF_HI)
    
    # Sample incidence uncertainty
    inc_65_74_i  <- rlnorm(1, log(INCIDENCE_65_74),
                           INCIDENCE_65_74_SE / INCIDENCE_65_74)
    inc_75_84_i  <- rlnorm(1, log(INCIDENCE_75_84),
                           INCIDENCE_75_84_SE / INCIDENCE_75_84)
    inc_85plus_i <- rlnorm(1, log(INCIDENCE_85PLUS),
                           INCIDENCE_85PLUS_SE / INCIDENCE_85PLUS)
    
    yearly <- lapply(YEARS, function(yr) {
      
      d <- exposure_with_pop %>% filter(year == yr)
      
      if (cf_mode == "empirical") {
        cf_draws <- cf_regional %>%
          filter(fire_year == yr) %>%
          mutate(cf_i = pmax(rnorm(n(), mean = cf_median, sd = cf_se), 0.5)) %>%
          select(air_zone, cf_i)
        
        d <- d %>% left_join(cf_draws, by = "air_zone")
      } else {
        d <- d %>% mutate(cf_i = fixed_cf)
      }
      
      d %>%
        mutate(
          rr_obs = gemm_rr(pm25_fireseason, theta_i,
                           GEMM_MU, GEMM_TAU, GEMM_R, xcf_i),
          rr_cf  = gemm_rr(cf_i, theta_i,
                           GEMM_MU, GEMM_TAU, GEMM_R, xcf_i),
          PAF    = pmax((rr_obs - rr_cf) / rr_obs, 0),
          
          # Lag-10: 45–64 cohort projected to 65–74 incidence band (~2031–2035)
          cohort_lag10   = pop_45_64 * surv_lag10,
          expected_lag10 = cohort_lag10 * inc_65_74_i,
          attr_lag10     = expected_lag10 * PAF,
          
          # Lag-15: 45–64 cohort projected to 75–84 incidence band (~2036–2040)
          cohort_lag15   = pop_45_64 * surv_lag15,
          expected_lag15 = cohort_lag15 * inc_75_84_i,
          attr_lag15     = expected_lag15 * PAF,
          
          # Lag-20: 45–64 cohort projected to 85+ incidence band (~2041–2045)
          cohort_lag20   = pop_45_64 * surv_lag20,
          expected_lag20 = cohort_lag20 * inc_85plus_i,
          attr_lag20     = expected_lag20 * PAF
        ) %>%
        summarise(
          year             = yr,
          attr_lag10_total = sum(attr_lag10, na.rm = TRUE),
          attr_lag15_total = sum(attr_lag15, na.rm = TRUE),
          attr_lag20_total = sum(attr_lag20, na.rm = TRUE)
        )
    })
    
    bind_rows(yearly)
  })
  
  bind_rows(results, .id = "sim_id")
}


mc_emp <- run_lag_simulation_gemm(cf_mode = "empirical")
mc_cf5 <- run_lag_simulation_gemm(cf_mode = "fixed", fixed_cf = 5)

summary_emp <- summarise_mc(mc_emp,  cf_label = "Empirical CF")
summary_cf5 <- summarise_mc(mc_cf5,  cf_label = "CF=5 \u03bcg/m\u00b3")

per_year_combined   <- bind_rows(summary_emp$per_year,  summary_cf5$per_year)
cumulative_combined <- bind_rows(summary_emp$cumulative, summary_cf5$cumulative)

print(cumulative_combined %>% print(width = Inf))

# DA-LEVEL BURDEN

da_boundaries <- st_read("bc_das.shp")

map_data <- exposure_with_pop %>%
  left_join(
    cf_regional %>% select(air_zone, fire_year, cf_median) %>%
      rename(year = fire_year),
    by = c("air_zone", "year")
  ) %>%
  mutate(
    rr_obs = gemm_rr(pm25_fireseason, GEMM_THETA_CENTRAL,
                     GEMM_MU, GEMM_TAU, GEMM_R, GEMM_XCF_LO),
    rr_cf  = gemm_rr(cf_median, GEMM_THETA_CENTRAL,
                     GEMM_MU, GEMM_TAU, GEMM_R, GEMM_XCF_LO),
    PAF    = pmax((rr_obs - rr_cf) / rr_obs, 0),
    attr_lag10_da = pop_45_64 * surv_lag10 * INCIDENCE_65_74  * PAF,
    attr_lag15_da = pop_45_64 * surv_lag15 * INCIDENCE_75_84  * PAF,
    attr_lag20_da = pop_45_64 * surv_lag20 * INCIDENCE_85PLUS * PAF
  ) %>%
  group_by(DGUID) %>%
  summarise(
    cumul_attr_lag10 = sum(attr_lag10_da, na.rm = TRUE),
    cumul_attr_lag15 = sum(attr_lag15_da, na.rm = TRUE),
    cumul_attr_lag20 = sum(attr_lag20_da, na.rm = TRUE),
    mean_pm25        = mean(pm25_fireseason, na.rm = TRUE),
    .groups          = "drop"
  )

map_sf <- da_boundaries %>%
  select(DGUID, geometry) %>%
  left_join(map_data, by = "DGUID")


# PLOT 1: cases by year and lag window

plot_data <- per_year_combined %>%
  pivot_longer(
    cols      = c(lag10_median, lag15_median, lag20_median),
    names_to  = "lag",
    values_to = "projected_cases"
  ) %>%
  mutate(
    lo95 = case_when(
      lag == "lag10_median" ~ lag10_lo95,
      lag == "lag15_median" ~ lag15_lo95,
      lag == "lag20_median" ~ lag20_lo95
    ),
    hi95 = case_when(
      lag == "lag10_median" ~ lag10_hi95,
      lag == "lag15_median" ~ lag15_hi95,
      lag == "lag20_median" ~ lag20_hi95
    ),
    lag_label = case_when(
      lag == "lag10_median" ~ "Lag-10 (projected ~2031\u20132035)",
      lag == "lag15_median" ~ "Lag-15 (projected ~2036\u20132040)",
      lag == "lag20_median" ~ "Lag-20 (projected ~2041\u20132045)"
    )
  ) %>%
  select(-lag10_lo95, -lag10_hi95, -lag15_lo95, 
         -lag15_hi95, -lag20_lo95, -lag20_hi95)

p_lag <- ggplot(plot_data,
                aes(x      = year,
                    y      = projected_cases,
                    colour = lag_label,
                    fill   = lag_label)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ counterfactual, ncol = 2) +
  scale_colour_manual(
    values = c("Lag-10 (projected ~2031\u20132035)" = "#E07B39",
               "Lag-15 (projected ~2036\u20132040)" = "#3B6E8F",
               "Lag-20 (projected ~2041\u20132045)" = "#6A3D9A")
  ) +
  scale_fill_manual(
    values = c("Lag-10 (projected ~2031\u20132035)" = "#E07B39",
               "Lag-15 (projected ~2036\u20132040)" = "#3B6E8F",
               "Lag-20 (projected ~2041\u20132045)" = "#6A3D9A")
  ) +
  labs(
    title    = "Projected Future Dementia Burden from Wildfire PM~2.5~ Exposure",
    subtitle = "GEMM nonlinear CRA (Ru et al. 2021; Burnett et al. 2018) \u00b7 Cases attributable to 45\u201364 cohort exposure, projected to peak incidence window",
    x        = "Fire Season Year (Exposure)",
    y        = "Projected Attributable Cases",
    colour   = "Lag Window",
    fill     = "Lag Window",
    caption  = paste0(
      "GEMM nonlinear CRA (Ru et al. 2021; Burnett et al. 2018), AAP+SHS model ",
      "(\u03b8 = 0.12, SE = 0.03; xcf ~ Uniform(2.7\u20137.6 \u03bcg/m\u00b3)).\n",
      "Empirical CF derived from non-fire-season (Nov\u2013Apr) BC MoE station medians per year.\n",
      "Survival correction from BC life tables (Statistics Canada 13-10-0114-01).\n",
      "Shaded bands: 95% uncertainty interval from Monte Carlo (n=", N_SIM, " simulations)." )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    plot.caption     = element_text(size = 8, colour = "grey40"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_markdown(size = 12, face = "plain")
  )

print(p_lag)
ggsave("lag_gemm_projection_plot.png", p_lag, width = 10, height = 5.5,
       dpi = 300, bg = "white")


mc_cumul_summary <- readRDS("mc_cumul_summary.rds")

cumul_concurrent <- mc_cumul_summary %>%
  filter(Scenario == "GEMM CRA") %>%
  transmute(
    frame        = Scenario,
    cases_median = Median,
    lo95         = CI_lo,
    hi95         = CI_hi
  )

cumul_lag <- cumulative_combined %>%
  filter(counterfactual == "Empirical CF") %>%
  pivot_longer(
    cols      = c(lag10_cumul_median, lag15_cumul_median,
                  lag20_cumul_median,
                  lag10_cumul_lo95,   lag15_cumul_lo95,
                  lag20_cumul_lo95,
                  lag10_cumul_hi95,   lag15_cumul_hi95,
                  lag20_cumul_hi95
                  ),
    names_to  = c("lag", ".value"),
    names_pattern = "(lag\\d+)_cumul_(.*)"
  ) %>%
  transmute(
    frame = case_when(
      lag == "lag10" ~ "Lag-10 projection (45\u201364 \u2192 65\u201374)",
      lag == "lag15" ~ "Lag-15 projection (45\u201364 \u2192 75\u201384)",
      lag == "lag20" ~ "Lag-20 projection (45\u201364 \u2192 85+)"
    ),
    cases_median = median,
    lo95         = lo95,
    hi95         = hi95
  )

cumul_comparison <- bind_rows(cumul_concurrent, cumul_lag) %>%
  mutate(ci = sprintf("%.0f (95%% CI: %.0f\u2013%.0f)",
                      cases_median, lo95, hi95))

print(cumul_comparison %>% select(frame, ci))

# Forest plot Lag analysis
forest_annual <- per_year_combined %>%
  filter(counterfactual == "Empirical CF") %>%
  pivot_longer(
    cols      = c(lag10_median, lag15_median, lag20_median,
                  lag10_lo95,   lag15_lo95,   lag20_lo95,
                  lag10_hi95,   lag15_hi95,   lag20_hi95),
    names_to  = c("lag", ".value"),
    names_pattern = "(lag\\d+)_(.*)"
  ) %>%
  rename(Median = median, CI_lo = lo95, CI_hi = hi95) %>%
  mutate(
    Year = factor(as.character(year),
                  levels = rev(as.character(YEARS))),
    lag_label = case_when(
      lag == "lag10" ~ "Lag-10 (~2031\u20132035)",
      lag == "lag15" ~ "Lag-15 (~2036\u20132040)",
      lag == "lag20" ~ "Lag-20 (~2041\u20132045)"
    ),
    label_txt = sprintf("%.0f (%.0f\u2013%.0f)", Median, CI_lo, CI_hi)
  )

# --- Cumulative data, empirical CF only ---
forest_cumul <- cumulative_combined %>%
  filter(counterfactual == "Empirical CF") %>%
  pivot_longer(
    cols      = c(lag10_cumul_median, lag15_cumul_median, lag20_cumul_median,
                  lag10_cumul_lo95,   lag15_cumul_lo95,   lag20_cumul_lo95,
                  lag10_cumul_hi95,   lag15_cumul_hi95,   lag20_cumul_hi95),
    names_to  = c("lag", ".value"),
    names_pattern = "(lag\\d+)_cumul_(.*)"
  ) %>%
  rename(Median = median, CI_lo = lo95, CI_hi = hi95) %>%
  mutate(
    lag_label = case_when(
      lag == "lag10" ~ "Lag-10 (~2031\u20132035)",
      lag == "lag15" ~ "Lag-15 (~2036\u20132040)",
      lag == "lag20" ~ "Lag-20 (~2041\u20132045)"
    ),
    label_txt = sprintf("%.0f (%.0f\u2013%.0f)", Median, CI_lo, CI_hi)
  ) %>%
  mutate(
    lag_label = factor(lag_label, levels = c(
      "Lag-10 (~2031\u20132035)",
      "Lag-15 (~2036\u20132040)",
      "Lag-20 (~2041\u20132045)"
    ))
  )

# --- Shared colour + shape scales ---
lag_colours <- c(
  "Lag-10 (~2031\u20132035)" = "#E07B39",
  "Lag-15 (~2036\u20132040)" = "#3B6E8F",
  "Lag-20 (~2041\u20132045)" = "#6A3D9A"
)
lag_shapes <- c(
  "Lag-10 (~2031\u20132035)" = 16,
  "Lag-15 (~2036\u20132040)" = 17,
  "Lag-20 (~2041\u20132045)" = 15
)

# --- Left panel: per-year ---
p_left <- ggplot(
  forest_annual,
  aes(x      = Median,
      y      = Year,
      colour = lag_label,
      shape  = lag_label)
) +
  geom_errorbar(
    aes(xmin = CI_lo, xmax = CI_hi),
    width       = 0.3,
    linewidth   = 0.75,
    orientation = "y",
    position    = position_dodge(width = 0.6)
  ) +
  geom_point(
    size     = 2.8,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    aes(x = CI_hi, label = label_txt),
    hjust    = -0.12,
    size     = 2.5,
    colour   = "grey35",
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  scale_colour_manual(values = lag_colours, name = "Lag window") +
  scale_shape_manual(values  = lag_shapes,  name = "Lag window") +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.02, 0.40))
  ) +
  labs(
    title    = "Annual burden by lag window",
    x        = "Projected attributable cases",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8.5),
    plot.title         = element_text(size = 10, face = "bold", colour = "grey30")
  )

# --- Right panel: cumulative ---
p_right <- ggplot(
  forest_cumul,
  aes(x      = Median,
      y      = lag_label,
      colour = lag_label,
      shape  = lag_label)
) +

  geom_errorbar(
    aes(xmin = CI_lo, xmax = CI_hi),
    width       = 0.25,
    linewidth   = 0.9,
    orientation = "y"
  ) +
  geom_point(size = 3.5) +
  geom_text(
    aes(x = CI_hi, label = label_txt),
    hjust  = -0.12,
    size   = 2.6,
    colour = "grey35"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  scale_colour_manual(values = lag_colours, guide = "none") +
  scale_shape_manual(values  = lag_shapes,  guide = "none") +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.02, 0.42))
  ) +
  labs(
    title = "Cumulative 2021\u20132025",
    x     = "Projected attributable cases",
    y     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    plot.title         = element_text(size = 12, face = "bold", colour = "grey30")
  )

# --- Combine with patchwork ---
p_forest_lag <- (p_left | p_right) +
  plot_layout(widths = c(1.6, 1), guides = "collect") &
  theme(legend.position = "bottom")

p_forest_lag <- p_forest_lag +
  plot_annotation(
    title    = "Projected Dementia Burden from Wildfire PM~2.5~ — GEMM CRA, BC 2021\u20132025",
    subtitle = "Attributable cases (median, 95% Monte Carlo CI) \u00b7 45\u201364 cohort projected to peak incidence windows \u00b7 Empirical counterfactual",
    caption  = paste0(
      "GEMM nonlinear CRA (Ru et al. 2021; Burnett et al. 2018), AAP+SHS model (\u03b8 = 0.12, SE = 0.03; xcf ~ Uniform(2.7\u20137.6 \u03bcg/m\u00b3)).\n",
      "Empirical CF: non-fire-season (Nov\u2013Apr) BC MoE station medians per year. ",
      "Lag-10: 65\u201374 incidence (2.80/1,000); lag-15: 75\u201384 (34.5/1,000); lag-20: 85+ (44.6/1,000).\n",
      "Survival correction from BC life tables (Statistics Canada 13-10-0114-01). Monte Carlo n = 10,000."
    ),
    theme = theme(
      plot.title    = element_markdown(face = "bold", size = 20),
      plot.subtitle = element_text(colour = "grey40", size = 15),
      plot.caption  = element_text(colour = "grey45", size = 11,
                                   lineheight = 1.4)
    )
  )

print(p_forest_lag)
ggsave("lag_gemm_forest_plot.png", p_forest_lag,
       width = 13, height = 5.5, dpi = 300, bg = "white")

write_csv(per_year_combined,   "lag_gemm_summary_by_year.csv")
write_csv(cumulative_combined, "lag_gemm_summary_cumulative.csv")
saveRDS(map_sf,                "lag_gemm_map_data.rds")
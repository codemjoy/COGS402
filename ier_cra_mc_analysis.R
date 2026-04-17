library(tidyverse)
library(sf)
library(ggtext)
library(scales)
library(ggplot2)
library(patchwork)
library(cowplot)
library(viridis)

set.seed(42)

years <- 2021:2025

PM25_CAP <- 250

da_by_year <- map(
  set_names(years),
  \(yr) {
    st_read(sprintf("da_annual_pm25_summary_%d.gpkg", yr), quiet = TRUE) |>
      st_drop_geometry() |>
      filter(!is.na(annual_mean_pm25), annual_mean_pm25 <= PM25_CAP) |>
      transmute(
        DGUID,
        pm25_mean      = annual_mean_pm25,
        pop_total_2021,
        pop_65_74      = pop_total_65_74,
        pop_75_84      = pop_total_75_84,
        pop_85plus     = pop_total_85plus,
        cd_code  = substr(DGUID, 10, 13),
        air_zone = case_when(
          cd_code %in% c("5955", "5957", "5959", "5939") ~ "Northeast",
          cd_code %in% c("5941", "5943", "5945", "5947",
                         "5949", "5951", "5953")          ~ "Central Interior",
          cd_code %in% c("5901", "5903", "5905", "5907",
                         "5909", "5911")                  ~ "Coastal",
          cd_code %in% c("5933", "5935", "5937")          ~ "Southern Interior",
          cd_code %in% c("5913", "5915", "5917", "5919",
                         "5921", "5923")                  ~ "Georgia Strait",
          cd_code %in% c("5925", "5927", "5929", "5931")  ~ "Lower Fraser Valley",
          TRUE ~ "Southern Interior"
        )
      )
  }
)

# Chen paramters
CHEN_HR_CENTRAL   <- 1.04   
CHEN_HR_LO        <- 1.02    
CHEN_HR_HI        <- 1.06
CHEN_IQR          <- 4.8   
CHEN_BETA_CENTRAL <- log(CHEN_HR_CENTRAL) / CHEN_IQR
CHEN_BETA_SE      <- (log(CHEN_HR_HI) - log(CHEN_HR_LO)) / (2 * 1.96) / CHEN_IQR

# Incidence rates 
INCIDENCE_65_74  <- 2.80 / 1000    # 65–74
INCIDENCE_75_84  <- 34.5 / 1000    # 75–84
INCIDENCE_85PLUS <- 44.6 / 1000    # 85+

# Fasaro Uncertainty
INCIDENCE_CV_YOUNG <- (617 - 604) / (2 * 1.96 * 610)  
INCIDENCE_CV_OLD   <- (3697 - 3640) / (2 * 1.96 * 3669) 

INCIDENCE_65_74_SE  <- INCIDENCE_65_74  * INCIDENCE_CV_YOUNG
INCIDENCE_75_84_SE  <- INCIDENCE_75_84  * INCIDENCE_CV_YOUNG
INCIDENCE_85PLUS_SE <- INCIDENCE_85PLUS * INCIDENCE_CV_OLD

# Counterfactuals
cf_regional <- read_csv("cf_annual_regional.csv") %>%
  filter(fire_year %in% years) %>%
  mutate(
    fire_year = as.integer(fire_year),
    cf_se     = (cf_p75 - cf_p25) / 2
  )

# Ru et al. (2021) GEMM parameters 
GEMM_THETA_CENTRAL <- 0.12
GEMM_THETA_SE      <- 0.03
GEMM_MU            <- -1.3      
GEMM_TAU           <- 0.1      
GEMM_R             <- 143.3     
GEMM_XCF_LO        <- 2.7     
GEMM_XCF_HI        <- 7.6       

# GEMM 
gemm_rr <- function(z, theta, mu, tau, r, xcf) {
  x     <- pmax(z - xcf, 0)
  f_x   <- log(x + 1)
  omega <- 1 / (1 + exp(-(x - mu) / (tau * r)))
  exp(theta * f_x * omega)
}

# CRA total cases — GEMM model
cra_cases_gemm <- function(pm25_vec, pop_65_74, pop_75_84, pop_85plus,
                           theta, mu, tau, r, xcf,
                           inc_65_74, inc_75_84, inc_85plus, cf) {
  rr_obs <- gemm_rr(pm25_vec, theta, mu, tau, r, xcf)
  rr_cf  <- gemm_rr(cf,       theta, mu, tau, r, xcf)
  paf    <- pmax((rr_obs - rr_cf) / rr_obs, 0)
  sum(paf * (pop_65_74  * inc_65_74 +
               pop_75_84  * inc_75_84 +
               pop_85plus * inc_85plus))
}

# CRA total cases — Chen log-linear model
cra_cases_chen <- function(pm25_vec, pop_65_74, pop_75_84, pop_85plus,
                           beta, inc_65_74, inc_75_84, inc_85plus, cf) {
  rr_obs <- exp(beta * pm25_vec)
  rr_cf  <- exp(beta * cf)
  paf    <- pmax((rr_obs - rr_cf) / rr_obs, 0)
  sum(paf * (pop_65_74 * inc_65_74 +
               pop_75_84 * inc_75_84 +
               pop_85plus * inc_85plus))
}

point_estimates <- map_dfr(years, \(yr) {
  d <- da_by_year[[as.character(yr)]] %>%
    left_join(
      cf_regional %>% filter(fire_year == yr) %>% select(air_zone, cf_median),
      by = "air_zone"
    )
  
  tibble(
    year = yr,
    point_chen_cf_regional = cra_cases_chen(
      d$pm25_mean, d$pop_65_74, d$pop_75_84, d$pop_85plus,
      beta       = CHEN_BETA_CENTRAL,
      inc_65_74  = INCIDENCE_65_74,
      inc_75_84  = INCIDENCE_75_84,
      inc_85plus = INCIDENCE_85PLUS,
      cf         = d$cf_median
    ),
    point_gemm_regional = cra_cases_gemm(
      d$pm25_mean, d$pop_65_74, d$pop_75_84, d$pop_85plus,
      theta      = GEMM_THETA_CENTRAL,
      mu         = GEMM_MU, tau = GEMM_TAU, r = GEMM_R,
      xcf        = GEMM_XCF_LO,           # use lower bound for point estimate
      inc_65_74  = INCIDENCE_65_74,
      inc_75_84  = INCIDENCE_75_84,
      inc_85plus = INCIDENCE_85PLUS,
      cf         = d$cf_median
    )
  )
})

print(point_estimates)
cat(sprintf("Cumulative Chen: %.0f  |  Cumulative GEMM: %.0f\n",
            sum(point_estimates$point_chen_cf_regional),
            sum(point_estimates$point_gemm_regional)))


# DA-days between 146 and 250 µg/m³ ?
da_all_years <- map_dfr(years, \(yr) {
  da_by_year[[as.character(yr)]] %>%
    mutate(year = yr)
})

da_all_years %>%
  summarise(
    total_da_days  = n(),
    n_146_to_250   = sum(pm25_mean > 146 & pm25_mean <= 250, na.rm = TRUE),
    pct_146_to_250 = round(mean(pm25_mean > 146 & pm25_mean <= 250, na.rm = TRUE) * 100, 2),
    n_above_146    = sum(pm25_mean > 146, na.rm = TRUE),
    pct_above_146  = round(mean(pm25_mean > 146, na.rm = TRUE) * 100, 2)
  ) %>%
  print()

# Broken down by year
da_all_years %>%
  group_by(year) %>%
  summarise(
    total_das      = n(),
    n_146_to_250   = sum(pm25_mean > 146 & pm25_mean <= 250, na.rm = TRUE),
    pct_146_to_250 = round(mean(pm25_mean > 146 & pm25_mean <= 250, na.rm = TRUE) * 100, 2),
    n_above_146    = sum(pm25_mean > 146, na.rm = TRUE),
    pct_above_146  = round(mean(pm25_mean > 146, na.rm = TRUE) * 100, 2),
    .groups = "drop"
  ) %>%
  print()

# MONTE CARLO
N_SIM <- 10000

# Results matrices: rows = simulations, columns = years
mc_chen_mat <- matrix(NA_real_, nrow = N_SIM, ncol = length(years),
                      dimnames = list(NULL, as.character(years)))
mc_gemm_mat <- matrix(NA_real_, nrow = N_SIM, ncol = length(years),
                      dimnames = list(NULL, as.character(years)))

# Parameter store for variance decomposition
mc_params <- tibble(
  sim        = 1:N_SIM,
  beta       = NA_real_,
  inc_65_74  = NA_real_,
  inc_75_84  = NA_real_,
  inc_85plus = NA_real_,
  gemm_theta = NA_real_,
  gemm_xcf   = NA_real_
)

pb_interval <- N_SIM / 10

mc_gemm_mat <- matrix(NA_real_, nrow = N_SIM, ncol = length(years),
                      dimnames = list(NULL, as.character(years)))

for (i in seq_len(N_SIM)) {
  if (i %% pb_interval == 0)
    cat(sprintf("  sim %d / %d\n", i, N_SIM))
  
  beta_i       <- rnorm(1, CHEN_BETA_CENTRAL, CHEN_BETA_SE)
  inc_65_74_i  <- rlnorm(1, log(INCIDENCE_65_74),  INCIDENCE_65_74_SE  / INCIDENCE_65_74)
  inc_75_84_i  <- rlnorm(1, log(INCIDENCE_75_84),  INCIDENCE_75_84_SE  / INCIDENCE_75_84)
  inc_85plus_i <- rlnorm(1, log(INCIDENCE_85PLUS), INCIDENCE_85PLUS_SE / INCIDENCE_85PLUS)
  
  # GEMM: only theta and xcf are uncertain; mu and tau are fixed modal values
  theta_i <- rnorm(1, GEMM_THETA_CENTRAL, GEMM_THETA_SE)
  xcf_i   <- runif(1, GEMM_XCF_LO, GEMM_XCF_HI)
  
  mc_params$beta[i]       <- beta_i
  mc_params$inc_65_74[i]  <- inc_65_74_i
  mc_params$inc_75_84[i]  <- inc_75_84_i
  mc_params$inc_85plus[i] <- inc_85plus_i
  mc_params$gemm_theta[i] <- theta_i
  mc_params$gemm_xcf[i]   <- xcf_i
  
  for (yr in years) {
    yr_chr <- as.character(yr)
    d      <- da_by_year[[yr_chr]]
    
    cf_draws <- cf_regional %>%
      filter(fire_year == as.integer(yr_chr)) %>%
      mutate(cf_i = pmax(rnorm(n(), mean = cf_median, sd = cf_se), 0.5)) %>%
      select(air_zone, cf_i)
    
    d_with_cf <- d %>% left_join(cf_draws, by = "air_zone")
    
    mc_chen_mat[i, yr_chr] <- cra_cases_chen(
      d_with_cf$pm25_mean, d_with_cf$pop_65_74,
      d_with_cf$pop_75_84, d_with_cf$pop_85plus,
      beta       = beta_i,
      inc_65_74  = inc_65_74_i, inc_75_84 = inc_75_84_i, inc_85plus = inc_85plus_i,
      cf         = d_with_cf$cf_i
    )
    
    mc_gemm_mat[i, yr_chr] <- cra_cases_gemm(
      d_with_cf$pm25_mean, d_with_cf$pop_65_74,
      d_with_cf$pop_75_84, d_with_cf$pop_85plus,
      theta      = theta_i,
      mu         = GEMM_MU, tau = GEMM_TAU, r = GEMM_R,
      xcf        = xcf_i,
      inc_65_74  = inc_65_74_i, inc_75_84 = inc_75_84_i, inc_85plus = inc_85plus_i,
      cf         = d_with_cf$cf_i
    )
  }
}

colSums(is.na(mc_chen_mat))
colSums(is.na(mc_gemm_mat))

# Cumulative (2021–2025) totals per simulation
mc_chen_cumul <- rowSums(mc_chen_mat)
mc_gemm_cumul  <- rowSums(mc_gemm_mat)

# --- Summarise MC results ---

# Annual summaries
mc_annual_summary <- bind_rows(
  map_dfr(as.character(years), \(yr) {
    sims <- mc_chen_mat[, yr]
    tibble(
      Year     = yr,
      Scenario = "Log-linear CRA",
      Median   = median(sims),
      CI_lo    = unname(quantile(sims, 0.025)),
      CI_hi    = unname(quantile(sims, 0.975))
    )
  }),
  map_dfr(as.character(years), \(yr) {
    sims <- mc_gemm_mat[, yr]
    tibble(
      Year     = yr,
      Scenario = "GEMM CRA",
      Median   = median(sims),
      CI_lo    = unname(quantile(sims, 0.025)),
      CI_hi    = unname(quantile(sims, 0.975))
    )
  })
)

# Cumulative summaries
mc_cumul_summary <- tibble(
  Scenario = c("Log-linear CRA", "GEMM CRA"),
  Median   = c(median(mc_chen_cumul),             median(mc_gemm_cumul)),
  CI_lo    = c(unname(quantile(mc_chen_cumul, 0.025)), unname(quantile(mc_gemm_cumul, 0.025))),
  CI_hi    = c(unname(quantile(mc_chen_cumul, 0.975)), unname(quantile(mc_gemm_cumul, 0.975)))
)

cat(sprintf("  Log-linear — median: %.0f  95%% UI: %.0f–%.0f\n",
            mc_cumul_summary$Median[1], mc_cumul_summary$CI_lo[1], mc_cumul_summary$CI_hi[1]))
cat(sprintf("  GEMM  — median: %.0f  95%% UI: %.0f–%.0f\n",
            mc_cumul_summary$Median[2], mc_cumul_summary$CI_lo[2], mc_cumul_summary$CI_hi[2]))

# Annual per-year summaries )
saveRDS(mc_annual_summary, "mc_annual_summary.rds")

# Cumulative summaries
saveRDS(mc_cumul_summary,  "mc_cumul_summary.rds")

# Full simulation matrices 
saveRDS(mc_chen_mat,       "mc_chen_mat.rds")
saveRDS(mc_gemm_mat,        "mc_gemm_mat.rds")

# --- Log-linear vs Gemm annual comparison figure ---
p_model_comparison <- mc_annual_summary %>%
  mutate(Year = as.integer(Year)) %>%
  ggplot(aes(x = Year, y = Median, colour = Scenario, fill = Scenario)) +
  geom_ribbon(aes(ymin = CI_lo, ymax = CI_hi), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(
    "Log-linear CRA" = "#C05040",
    "GEMM CRA"       = "#2A6AB0"
  )) +
  scale_fill_manual(values = c(
    "Log-linear CRA" = "#E07B6A",
    "GEMM CRA"       = "#4A90D9"
  )) +
  scale_x_continuous(breaks = 2021:2025) +
  scale_y_continuous(labels = comma, limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Annual Attributable Dementia Cases: Log-Linear vs GEMM CRA",
    subtitle = "BC, 2021–2025 · Empirical CF · Median and 95% CI from Monte Carlo (n=10,000)",
    x        = "Year",
    y        = "Attributable cases",
    colour   = NULL, fill = NULL,
    caption  = paste0(
      "Log-linear CRA: HR 1.04 per 4.8 µg/m³ IQR (Chen et al., 2017). ",
      "Nonlinear CRA: GEMM (Ru et al., 2021; Burnett et al., 2018), AAP+SHS model ",
      "(θ = 0.12, SE = 0.03; xcf ~ Uniform(2.7–7.6 µg/m³)).\n",
      "Empirical counterfactual: BC MoE non-fire-season annual median. ",
      "Daily cap: 250 µg/m³. Monte Carlo n = 10,000."
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(colour = "grey40", size = 10),
    plot.caption     = element_text(colour = "grey50", size = 8.5,
                                    lineheight = 1.4, hjust = 0),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

print(p_model_comparison)
ggsave("fig_chen_vs_gemm_annual.png",
       p_model_comparison, width = 9, height = 5.5, dpi = 300)


# --- Colours ---
scenario_colours <- c(
  "Log-linear CRA" = "#C05040",
  "GEMM CRA"       = "#2A6AB0"
)

# --- Prepare data ---
annual_data <- mc_annual_summary %>%
  mutate(
    Year   = as.integer(Year),
    CI_lo  = unname(CI_lo),
    CI_hi  = unname(CI_hi),
    Scenario = factor(Scenario, levels = c("Log-linear CRA","GEMM CRA"))
  )

cumul_data <- mc_cumul_summary %>%
  mutate(
    CI_lo     = unname(CI_lo),
    CI_hi     = unname(CI_hi),
    plot_year = "2021–2025",
      Scenario  = factor(Scenario,levels = c("Log-linear CRA", "GEMM CRA"))
    )

cat("Annual rows:", nrow(annual_data), "\n")
cat("Cumul rows:",  nrow(cumul_data),  "\n")

# --- Annual panel ---
p_annual <- ggplot(
  annual_data,
  aes(x     = Median,
      y     = factor(Year, levels = rev(c(2021,2022,2023,2024,2025))),
      colour = Scenario)
) +
  geom_vline(xintercept = seq(0, 2000, by = 500),
             colour = "grey90", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = CI_lo, xmax = CI_hi),
    width       = 0.25,
    linewidth   = 0.8,
    orientation = "y",
    position    = position_dodge(width = 0.6)
  ) +
  geom_point(
    size     = 2.8,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = annual_data,
    aes(x      = CI_hi,
        y      = factor(Year, levels = rev(c(2021,2022,2023,2024,2025))),
        colour = Scenario,
        label  = paste0(round(Median), " (",
                        round(CI_lo),  "–",
                        round(CI_hi),  ")")),
    hjust    = -0.08,
    size     = 2.8,
    position = position_dodge(width = 0.6)
  ) +
  scale_colour_manual(values = scenario_colours, guide = "none") +
  scale_x_continuous(
    limits = c(0, 2000),
    breaks = seq(0, 2000, by = 500),
    expand = expansion(mult = c(0, 0.02)),
    labels = comma
  ) +
  labs(title = "**Annual**", x = "Attributable cases", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title         = element_markdown(size = 20, hjust = 0.5),
    panel.grid         = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
    axis.text.y        = element_text(size = 13, colour = "grey20"),
    axis.text.x        = element_text(size = 13, colour = "grey40"),
    plot.margin        = margin(8, 32, 8, 8)
  )

# --- Cumulative panel ---
p_cumul <- ggplot(
  cumul_data,
  aes(x      = Median,
      y      = plot_year,
      colour = Scenario)
) +
  geom_vline(xintercept = seq(0, 3500, by = 500),
             colour = "grey90", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = CI_lo, xmax = CI_hi),
    width       = 0.15,
    linewidth   = 0.8,
    orientation = "y",
    position    = position_dodge(width = 0.25)
  ) +
  geom_point(
    size     = 2.8,
    position = position_dodge(width = 0.25)
  ) +
  geom_text(
    aes(x     = CI_hi,
        label = paste0(round(Median), " (",
                       round(CI_lo),  "–",
                       round(CI_hi),  ")")),
    hjust    = -0.08,
    size     = 2.8,
    position = position_dodge(width = 0.25)
  ) +
  scale_colour_manual(values = scenario_colours, guide = "none") +
  scale_x_continuous(
    limits = c(0, 4500),
    breaks = seq(0, 3500, by = 500),
    expand = expansion(mult = c(0, 0.02)),
    labels = comma
  ) +
  labs(title = "**Cumulative 2021–2025**", x = "Attributable cases", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title         = element_markdown(size = 20, hjust = 0.5),
    panel.grid         = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
    axis.text.y        = element_text(size = 14, colour = "grey20"),
    axis.text.x        = element_text(size = 14, colour = "grey40"),
    plot.margin        = margin(8, 8, 8, 32)
  )

# --- Shared legend ---
legend_plot <- ggplot(
  annual_data,
  aes(x = Median, y = factor(Year), colour = Scenario)
) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = scenario_colours, name = NULL) +
  theme_void() +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 11)
  )

shared_legend <- cowplot::get_legend(legend_plot)

# --- Combine ---
p_combined <- wrap_plots(p_annual, p_cumul, ncol = 2, widths = c(1.6, 1)) +
  plot_annotation(
    title    = "Attributable Dementia Cases: Annual and Cumulative",
    subtitle = "Median (95% Monte Carlo CI) · Log-Linear (Chen et al., 2017) vs GEMM (Ru et al., 2021)",
    caption  = paste0(
      "Log-linear CRA: HR 1.04 per 4.8 µg/m³ IQR (Chen et al., 2017). ",
      "Nonlinear CRA: GEMM (Ru et al., 2021; Burnett et al., 2018), AAP+SHS model ",
      "(θ = 0.12, SE = 0.03;\n",
      "xcf ~ Uniform(2.7–7.6 µg/m³)). Empirical counterfactual: BC MoE stations non-fire-season annual median. ",
      "Monte Carlo n = 10,000. Daily cap: 250 µg/m³."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(colour = "grey40", size = 15),
      plot.caption  = element_text(colour = "grey50", size = 11,
                                   lineheight = 1.4, hjust = 0)
    )
  )

final_plot <- cowplot::plot_grid(
  shared_legend, p_combined,
  ncol        = 1,
  rel_heights = c(0.05, 1)
)


print(final_plot)
ggsave("fig_forest_annual_cumul.png",
       final_plot, width = 16, height = 6, dpi = 300, bg = "white")

# FIGURE: GEMM Burden Map (2021 vs 2023) + Top 5 Bar Charts
# --- Per-DA GEMM point estimates ---
da_burden_map <- map_dfr(c(2021, 2023), \(yr) {
  d <- da_by_year[[as.character(yr)]] %>%
    left_join(
      cf_regional %>% filter(fire_year == yr) %>% select(air_zone, cf_median),
      by = "air_zone"
    )
  
  rr_obs <- gemm_rr(d$pm25_mean, GEMM_THETA_CENTRAL, GEMM_MU, GEMM_TAU, GEMM_R, GEMM_XCF_LO)
  rr_cf  <- gemm_rr(d$cf_median, GEMM_THETA_CENTRAL, GEMM_MU, GEMM_TAU, GEMM_R, GEMM_XCF_LO)
  paf <- pmax((rr_obs - rr_cf) / rr_obs, 0)
  d |>
    mutate(
      year               = yr,
      attributable_cases = paf * (pop_65_74  * INCIDENCE_65_74 +
                                    pop_75_84  * INCIDENCE_75_84 +
                                    pop_85plus * INCIDENCE_85PLUS)
    ) |>
    select(DGUID, year, attributable_cases, pm25_mean)
})

# --- Spatial join ---
bc_geom <- st_read("da_annual_pm25_summary_2021.gpkg", quiet = TRUE) |>
  select(DGUID)

burden_sf <- bc_geom |>
  st_simplify(preserveTopology = TRUE, dTolerance = 100) |>
  left_join(da_burden_map, by = "DGUID") |>
  filter(!is.na(year)) |>
  mutate(
    year_label    = paste0(year, " Fire Season"),
    cases_display = sqrt(pmax(attributable_cases, 0))
  )

# Colour scale parameters (shared between map and bars)
real_breaks  <- c(0, 0.1, 0.5, 1, 2, 5)
sqrt_breaks  <- sqrt(real_breaks)
break_labels <- c("0", "0.1", "0.5", "1", "2", "5+")
bar_fill     <- "#C05040"   

# --- MAP ---
p_map <- ggplot(burden_sf) +
  geom_sf(aes(fill = cases_display), colour = NA, linewidth = 0) +
  facet_wrap(~ year_label, ncol = 2) +
  coord_sf(crs = st_crs(3005)) +
  scale_fill_viridis_c(
    option    = "magma",
    direction = -1,
    name      = "Attributable\ncases",
    na.value  = "grey85",
    begin     = 0.15,
    end       = 1,
    breaks    = sqrt_breaks,
    labels    = break_labels,
    limits    = c(0, sqrt(6)),
    oob       = scales::squish
  ) +
  labs(
    title   = "Wildfire PM~2.5~-Attributable Alzheimer's Burden by British Columbian Dissemination Areas",
    subtitle = "Nonlinear CRA: GEMM (Ru et al., 2021; Burnett et al., 2018), AAP+SHS model (θ = 0.12, SE = 0.03)"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    plot.title        = element_markdown(face = "bold", size = 12, hjust = 0),
    plot.subtitle     = element_text(colour = "grey40", size = 10, hjust = 0),
    strip.text        = element_text(face = "bold", size = 11,
                                     margin = margin(b = 6)),
    legend.position   = "bottom",
    legend.key.width  = unit(3, "cm"),
    legend.key.height = unit(0.5, "cm"),
    plot.margin       = margin(10, 10, 10, 10),
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 9)
  )

add_region <- function(df) {
  df |>
    mutate(
      cd_code = substr(DGUID, 10, 13),   
      region  = case_when(
        cd_code == "5901" ~ "Capital / Greater Victoria",
        cd_code == "5903" ~ "Cowichan Valley",
        cd_code == "5905" ~ "Nanaimo",
        cd_code == "5907" ~ "Alberni-Clayoquot / Central Island",
        cd_code == "5909" ~ "Comox Valley",
        cd_code == "5911" ~ "Strathcona",
        cd_code == "5913" ~ "Powell River",
        cd_code == "5915" ~ "Metro Vancouver",
        cd_code == "5917" ~ "Sunshine Coast",
        cd_code == "5919" ~ "Squamish-Lillooet",
        cd_code == "5921" ~ "Fraser Valley",
        cd_code == "5923" ~ "Chilliwack",
        cd_code == "5925" ~ "Hope",
        cd_code == "5927" ~ "Central Fraser Valley",
        cd_code == "5929" ~ "Langley / Surrey",
        cd_code == "5931" ~ "Greater Vancouver South",
        cd_code == "5933" ~ "Thompson-Nicola / Kamloops",
        cd_code == "5935" ~ "Central Okanagan / Kootenay",
        cd_code == "5937" ~ "Okanagan-Similkameen",
        cd_code == "5939" ~ "Peace River / Northeast BC",
        cd_code == "5941" ~ "Cariboo",
        cd_code == "5943" ~ "Skeena",
        cd_code == "5945" ~ "Kitimat-Stikine",
        cd_code == "5947" ~ "Bulkley-Nechako",
        cd_code == "5949" ~ "Fraser-Fort George / Prince George",
        cd_code == "5951" ~ "Stikine",
        cd_code == "5953" ~ "Central Okanagan",
        cd_code == "5955" ~ "North Okanagan",
        cd_code == "5957" ~ "Columbia-Shuswap",
        cd_code == "5959" ~ "East Kootenay",
        TRUE ~ paste0("CD ", cd_code)
      )
    )
}

# --- TOP 5 BAR CHART helper ---
make_top5_bar <- function(yr, data) {
  top5 <- data |>
    filter(year == yr) |>
    slice_max(attributable_cases, n = 5) |>
    add_region()
  
  top5 <- top5 |>
    arrange(attributable_cases) |>
    mutate(
      da_label = paste0(region, " (", substr(DGUID, 14, 17), ")"),
      da_label = factor(da_label, levels = da_label),
      label    = sprintf("%.2f", attributable_cases)
    )
  
  ggplot(top5, aes(x = attributable_cases, y = da_label)) +
    geom_col(fill = bar_fill, width = 0.6) +
    geom_text(aes(label = label),
              hjust = -0.15, size = 3.2, colour = "grey20") +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.25)),
      labels = comma
    ) +
    labs(
      title = paste0("**Top 5 DAs — ", yr, "**"),
      x     = "Attributable cases",
      y     = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title        = element_markdown(size = 11, hjust = 0.5),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.grid        = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
      axis.text.y       = element_text(size = 9, colour = "grey20"),
      axis.text.x       = element_text(size = 9, colour = "grey40"),
      plot.margin       = margin(4, 16, 4, 4)
    )
}

p_bar_2021 <- make_top5_bar(2021, da_burden_map)
p_bar_2023 <- make_top5_bar(2023, da_burden_map)

# --- COMBINE with patchwork ---
p_bars <- p_bar_2021 / p_bar_2023 +
  plot_layout(heights = c(1, 1))

p_final <- p_map + p_bars +
  plot_layout(widths = c(2.5, 1)) +
  plot_annotation(
    caption = paste0(
      "Nonlinear CRA: GEMM (Ru et al., 2021; Burnett et al., 2018), AAP+SHS model (θ = 0.12, SE = 0.03). ",
      "Empirical counterfactual: BC MoE non-fire-season annual median.\n",
      "Yellow areas: no exposure above counterfactual."
    ),
    theme = theme(
      plot.caption     = element_text(colour = "grey50", size = 8.5,
                                      lineheight = 1.4, hjust = 0),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  )

print(p_final)
ggsave("fig_gemm_burden_map_top5.png",
       p_final, width = 16, height = 8, dpi = 300, bg = "white")



lower_mainland_cds <- c("5915", "5917", "5919", "5921", "5923",
                        "5925", "5927", "5929", "5931")

lm_burden <- map_dfr(c(2021, 2023), \(yr) {
  d <- da_by_year[[as.character(yr)]] %>%
    left_join(
      cf_regional %>% filter(fire_year == yr) %>% select(air_zone, cf_median),
      by = "air_zone"
    ) %>%
    filter(cd_code %in% lower_mainland_cds)
  
  rr_obs <- gemm_rr(d$pm25_mean, GEMM_THETA_CENTRAL, GEMM_MU, GEMM_TAU, GEMM_R, mean(c(GEMM_XCF_LO, GEMM_XCF_HI)))
  rr_cf  <- gemm_rr(d$cf_median, GEMM_THETA_CENTRAL, GEMM_MU, GEMM_TAU, GEMM_R, mean(c(GEMM_XCF_LO, GEMM_XCF_HI)))
  paf    <- pmax((rr_obs - rr_cf) / rr_obs, 0)
  
  tibble(
    year               = yr,
    attributable_cases = round(sum(paf * (d$pop_65_74 * INCIDENCE_65_74 +
                                            d$pop_75_84 * INCIDENCE_75_84 +
                                            d$pop_85plus * INCIDENCE_85PLUS),
                                   na.rm = TRUE), 2),
    pop_65plus         = sum(d$pop_65_74 + d$pop_75_84 + d$pop_85plus, na.rm = TRUE),
    n_das_exposed      = sum(d$pm25_mean > d$cf_median, na.rm = TRUE),
    mean_pm25_wtd      = weighted.mean(d$pm25_mean,
                                       d$pop_65_74 + d$pop_75_84 + d$pop_85plus,
                                       na.rm = TRUE)
  )
})

print(lm_burden)
library(tidyverse)
library(sf)
library(ggtext)

# ---- 1. PARAMETERS ------------------------------------------

# Chen et al. (2017) HR and 95% CI on the log scale, per unit PM2.5
hr_central <- 1.04
hr_low     <- 1.03
hr_high    <- 1.05
iqr        <- 4.8          # IQR used in Chen et al.

# Convert to beta (ln(HR) per µg/m³)
beta_mean <- log(hr_central) / iqr

beta_se   <- (log(hr_high) - log(hr_low)) / (2 * 1.96 * iqr)

cat("Beta mean:", round(beta_mean, 6), "\n")
cat("Beta SE:  ", round(beta_se,   6), "\n")
cat("Implied 95% CI on beta: [",
    round(beta_mean - 1.96 * beta_se, 6), ",",
    round(beta_mean + 1.96 * beta_se, 6), "]\n\n")

# Monte Carlo settings
n_sims       <- 10000
set.seed(2402)  

counterfactual <- 15   # WHO AQG µg/m³

primary_rates <- list(
  rate_65_74  =  6.10 / 1000,
  rate_75_84  = 34.50 / 1000,
  rate_85plus = 44.68 / 1000
)

years <- 2021:2026

# ---- SAMPLE BETA DISTRIBUTION ----------------------------

beta_samples <- rnorm(n_sims, mean = beta_mean, sd = beta_se)

cat("Sampled beta distribution:\n")
cat("  Mean:", round(mean(beta_samples), 6), "\n")
cat("  2.5%:", round(quantile(beta_samples, 0.025), 6), "\n")
cat(" 97.5%:", round(quantile(beta_samples, 0.975), 6), "\n\n")

# ---- PRELOAD EXPOSURE DATA -------------------------------.
year_data <- map(years, function(yr) {
  
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
    st_drop_geometry() %>%
    select(DGUID, DAUID,
           pop_total_65_74, pop_total_75_84, pop_total_85plus)
  
  annual %>%
    left_join(fireseason, by = "DGUID") %>%
    mutate(
      delta_pm25      = pmax(fireseason_mean_pm25 - counterfactual, 0),
      expected_65_74  = pop_total_65_74  * primary_rates$rate_65_74,
      expected_75_84  = pop_total_75_84  * primary_rates$rate_75_84,
      expected_85plus = pop_total_85plus * primary_rates$rate_85plus,
      expected_total  = expected_65_74 + expected_75_84 + expected_85plus,
      year            = yr
    ) %>%
    filter(!is.na(delta_pm25))
})

names(year_data) <- years

# ---- MONTE CARLO LOOP ------------------------------------

# Pre-extract the vectors 
prepped <- map(year_data, function(df) {
  list(
    delta    = df$delta_pm25,
    expected = df$expected_total,
    yr       = df$year[1]
  )
})

# Matrix approach: rows = simulations, cols = years
mc_results <- matrix(NA_real_,
                     nrow = n_sims,
                     ncol = length(years),
                     dimnames = list(NULL, as.character(years)))

pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

for (i in seq_len(n_sims)) {
  b <- beta_samples[i]
  for (j in seq_along(years)) {
    d <- prepped[[j]]
    RR  <- exp(b * d$delta)
    PAF <- (RR - 1) / RR
    mc_results[i, j] <- sum(d$expected * PAF, na.rm = TRUE)
  }
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\nMonte Carlo complete.\n\n")

# ---- SUMMARISE RESULTS -----------------------------------

# Per-year summary
mc_year_summary <- map_dfr(seq_along(years), function(j) {
  sims <- mc_results[, j]
  tibble(
    year           = years[j],
    central        = median(sims),          
    mean_sim       = mean(sims),
    ci_lower       = quantile(sims, 0.025),
    ci_upper       = quantile(sims, 0.975),
    ci_width       = quantile(sims, 0.975) - quantile(sims, 0.025)
  )
}) %>%
  mutate(across(where(is.numeric), ~round(., 1)))

cat("=== PER-YEAR MONTE CARLO RESULTS ===\n")
print(mc_year_summary)

# Cumulative (sum across years for each simulation)
mc_cumulative <- rowSums(mc_results)

cat("\n=== CUMULATIVE 2021-2026 ===\n")
cat("Median:      ", round(median(mc_cumulative), 1), "\n")
cat("Mean:        ", round(mean(mc_cumulative),   1), "\n")
cat("95% CI:      [",
    round(quantile(mc_cumulative, 0.025), 1), ",",
    round(quantile(mc_cumulative, 0.975), 1), "]\n\n")

# ---- 6. COMPARE TO ORIGINAL THREE-POINT SENSITIVITY ---------
# Your existing low/central/high for reference
beta_central_orig <- log(1.04) / 4.8
beta_low_orig     <- log(1.03) / 4.8
beta_high_orig    <- log(1.05) / 4.8

three_point_check <- map_dfr(seq_along(years), function(j) {
  d <- prepped[[j]]
  tibble(
    year   = years[j],
    low    = sum(d$expected * ((exp(beta_low_orig  * d$delta) - 1) /
                                exp(beta_low_orig  * d$delta)), na.rm = TRUE),
    central = sum(d$expected * ((exp(beta_central_orig * d$delta) - 1) /
                                 exp(beta_central_orig * d$delta)), na.rm = TRUE),
    high   = sum(d$expected * ((exp(beta_high_orig * d$delta) - 1) /
                                exp(beta_high_orig * d$delta)), na.rm = TRUE)
  )
}) %>%
  mutate(across(where(is.numeric), ~round(., 1)))

cat("=== COMPARISON: MC 95% CI vs THREE-POINT SENSITIVITY ===\n")
comparison <- mc_year_summary %>%
  select(year, mc_lower = ci_lower, mc_median = central, mc_upper = ci_upper) %>%
  left_join(three_point_check %>% select(year, tp_low = low,
                                          tp_central = central,
                                          tp_high = high),
            by = "year")

print(comparison)

write_csv(mc_year_summary, "mc_uncertainty_by_year.csv")

mc_cumulative_df <- tibble(
  sim_id     = seq_len(n_sims),
  total_cases = mc_cumulative
)
write_csv(mc_cumulative_df, "mc_cumulative_distribution.csv")


# ---- 8. PLOT: Uncertainty distribution (cumulative) ---------
p_mc <- ggplot(mc_cumulative_df, aes(x = total_cases)) +
  
  geom_histogram(
    bins  = 80,
    fill  = "#ca706a",
    alpha = 0.75,
    color = NA
  ) +
  
  geom_vline(
    xintercept = median(mc_cumulative),
    color = "#7a2c2c", linewidth = 1, linetype = "solid"
  ) +
  geom_vline(
    xintercept = c(quantile(mc_cumulative, 0.025),
                   quantile(mc_cumulative, 0.975)),
    color = "#7a2c2c", linewidth = 0.8, linetype = "dashed"
  ) +
  
  annotate("text",
           x = median(mc_cumulative),
           y = Inf, vjust = 1.5, hjust = -0.1,
           label = paste0("Median: ", round(median(mc_cumulative), 0)),
           size = 3.5, color = "#7a2c2c", fontface = "bold") +
  
  annotate("text",
           x = quantile(mc_cumulative, 0.975),
           y = Inf, vjust = 3, hjust = -0.1,
           label = paste0("97.5%: ", round(quantile(mc_cumulative, 0.975), 0)),
           size = 3, color = "#7a2c2c") +
  
  annotate("text",
           x = quantile(mc_cumulative, 0.025),
           y = Inf, vjust = 3, hjust = 1.1,
           label = paste0("2.5%: ", round(quantile(mc_cumulative, 0.025), 0)),
           size = 3, color = "#7a2c2c") +
  
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  labs(
    title    = "Monte Carlo Distribution of Cumulative Attributable Cases",
    subtitle = paste0("2021–2026 · ", scales::comma(n_sims),
                      " simulations · β ~ N(mean, SE) from Chen et al. (2017) 95% CI"),
    x        = "Total attributable Alzheimer's cases (2021–2026)",
    y        = "Frequency",
    caption  = paste0(
      "β sampled from normal distribution with mean = log(1.04)/4.8 and\n",
      "SE = (log(1.05) − log(1.03)) / (2 × 1.96 × 4.8). Counterfactual: 15 µg/m³."
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle    = element_text(color = "grey40", size = 9, hjust = 0,
                                    margin = margin(b = 10)),
    plot.caption     = element_text(color = "grey50", size = 8, hjust = 0,
                                    lineheight = 1.4, margin = margin(t = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background  = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin      = margin(16, 16, 16, 16)
  )

ggsave(
  "mc_cumulative_distribution.png",
  plot = p_mc, width = 8, height = 5, dpi = 300,
  bg = "#ffffffC9"
)

# ---- UPDATED BAR PLOT WITH MC CIs -----------------------
plot_data_mc <- mc_year_summary %>%
  mutate(
    is_small            = central < 50,
    ci_label_y          = ifelse(is_small, ci_upper + 25, ci_upper + 8)
  )

p_mc_bars <- ggplot(plot_data_mc, aes(x = factor(year))) +
  # MC uncertainty range
  geom_rect(
    aes(
      xmin = as.numeric(factor(year)) - 0.35,
      xmax = as.numeric(factor(year)) + 0.35,
      ymin = ci_lower,
      ymax = ci_upper
    ),
    fill = "#ca706a", alpha = 0.25
  ) +
  
  # Central (median) bar
  geom_col(
    aes(y = central),
    fill = "#ca706a", alpha = 0.85, width = 0.5
  ) +
  
  # Central label
  geom_text(
    aes(y = central, label = round(central, 0)),
    vjust = -0.3, size = 4, fontface = "bold", color = "#1a1f2e"
  ) +
  
  # CI label
  geom_text(
    aes(y = ci_label_y,
        label = paste0("(", round(ci_lower, 0), "–", round(ci_upper, 0), ")")),
    vjust = 0, size = 4, color = "#1a1f2e"
  ) +
  
  # Dotted connector for small bars
  geom_segment(
    data = plot_data_mc %>% filter(is_small),
    aes(
      x = as.numeric(factor(year)), xend = as.numeric(factor(year)),
      y = ci_upper + 3, yend = ci_upper + 20
    ),
    color = "grey60", linewidth = 0.3, linetype = "dotted"
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25)),
    labels = scales::comma
  ) +
  
  labs(
    title    = "Estimated Alzheimer's Cases Attributable to Wildfire PM<sub>2.5</sub>",
    subtitle = "British Columbia, 2021\u20132026 Fire Seasons (May\u2013September) \u00b7 Population Aged 65+",
    x        = "Year",
    y        = "Attributable Alzheimer's cases",
    caption  = paste0(
      "Bars: median of Monte Carlo distribution (10,000 simulations).\n",
      "Shaded range: 95% CI from probabilistic sampling of \u03b2 ~ N(log(1.04)/4.8, SE),\n",
      "where SE derived from Chen et al. (2017) HR 95% CI (1.03\u20131.05 per 4.8 \u00b5g/m\u00b3 IQR).\n",
      "Counterfactual: 15 \u00b5g/m\u00b3 (WHO AQG). BlueSky Canada smoke forecast system. Statistics Canada 2021 Census."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_markdown(face = "bold", size = 20, hjust = 0),
    plot.subtitle      = element_markdown(color = "grey40", size = 16, hjust = 0),
    plot.caption       = element_text(color = "grey50", size = 13, hjust = 0,
                                      lineheight = 1.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "grey30"),
    axis.title         = element_text(color = "grey30", size = 13),
    plot.background    = element_rect(fill = "#ffffffC9", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

ggsave(
  "attributable_cases_mc_ci.png",
  plot   = p_mc_bars,
  width  = 11,
  height = 7,
  dpi    = 300,
  bg     = "#ffffffC9"
)


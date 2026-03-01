library(sf)
library(dplyr)
library(terra)
library(exactextractr)
library(lubridate)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyr)

# Set year parameter (remember to cd at same time!!)
year <- 2021
bluesky_dir <- "/Users/minalander/Library/Mobile Documents/com~apple~CloudDocs/Documents/University /Year 5/COGS 402/Data/Blue skys tif /2021"

census_das <- st_read("bc_das_population_for_raster.gpkg") 
census_das_transformed <- st_transform(census_das, 4326)
bluesky_files <- list.files(bluesky_dir,
                            pattern = "\\.tif$",
                            full.names = TRUE)

print(paste("Found", length(bluesky_files), "files for year", year))

extract_date <- function(filename) {
  date_str <- str_extract(basename(filename), "(?<=_)\\d{8}(?=_)")
  ymd(date_str)
}

file_df <- data.frame(
  filepath = bluesky_files
) %>%
  mutate(date = extract_date(filepath)) %>%
  filter(!is.na(date)) %>%
  arrange(date)

daily_results <- list()

for(i in seq_along(file_df$filepath)) {
  current_date <- file_df$date[i]
  current_file <- file_df$filepath[i]
  pm25_raster_full <- rast(current_file)
  pm25_raster <- pm25_raster_full[[1]]
  da_mean <- exact_extract(pm25_raster, census_das_transformed,fun = 'mean',
    progress = FALSE)
  daily_results[[i]] <- data.frame(DGUID = census_das_transformed$DGUID,
    date = current_date,pm25_mean = da_mean)
  if(i %% 50 == 0) {
    print(paste("Processed", i, "of", nrow(file_df),
                "- Date:", current_date))
  }
}

da_daily_exposure <- bind_rows(daily_results)

da_daily_exposure <- da_daily_exposure %>%
  mutate(
    above_15 = pm25_mean > 15,
    above_37.5 = pm25_mean > 37.5,
    above_50 = pm25_mean > 50
  )

output_file_rds <- paste0("da_daily_pm25_exposure_", year, ".rds")
output_file_csv <- paste0("da_daily_pm25_exposure_", year, ".csv")

saveRDS(da_daily_exposure, output_file_rds)
write.csv(da_daily_exposure, output_file_csv, row.names = FALSE)

print(paste("\n✓ Daily exposure data for", year, "saved"))

da_annual_summary <- da_daily_exposure %>%
  summarise(
    .by = DGUID,
    annual_mean_pm25 = mean(pm25_mean, na.rm = TRUE),
    days_above_15 = sum(above_15, na.rm = TRUE),
    days_above_37.5 = sum(above_37.5, na.rm = TRUE), 
    days_above_50 = sum(above_50, na.rm = TRUE),
    n_days = n()
  ) %>%
  mutate(year = year)

da_annual_spatial <- census_das %>%
  left_join(da_annual_summary, by = "DGUID")

output_gpkg <- paste0("da_annual_pm25_summary_", year, ".gpkg")
st_write(da_annual_spatial, output_gpkg, delete_layer = TRUE)

print(paste("✓ Annual summary for", year, "saved"))

# Age-stratified exposure profiles:
da_exposure_pop <- da_annual_summary %>%
  left_join(
    census_das %>% st_drop_geometry() %>%
      select(DGUID, starts_with("pop_total_"), starts_with("pop_male_"), 
             starts_with("pop_female_")),
    by = "DGUID"
  )

# calculate population-weighted exposures by age-group

age_stratified <- da_exposure_pop %>%
  mutate(
    exposure_0_4 = pop_total_0_4 * annual_mean_pm25,
    exposure_5_9 = pop_total_5_9 * annual_mean_pm25,
    exposure_10_14 = pop_total_10_14 * annual_mean_pm25,
    exposure_15_19 = pop_total_15_19 * annual_mean_pm25,
    exposure_20_24 = pop_total_20_24 * annual_mean_pm25,
    exposure_25_29 = pop_total_25_29 * annual_mean_pm25,
    exposure_30_34 = pop_total_30_34 * annual_mean_pm25,
    exposure_35_39 = pop_total_35_39 * annual_mean_pm25,
    exposure_40_44 = pop_total_40_44 * annual_mean_pm25,
    exposure_45_49 = pop_total_45_49 * annual_mean_pm25,
    exposure_50_54 = pop_total_50_54 * annual_mean_pm25,
    exposure_55_59 = pop_total_55_59 * annual_mean_pm25,
    exposure_60_64 = pop_total_60_64 * annual_mean_pm25,
    exposure_65_69 = pop_total_65_69 * annual_mean_pm25,
    exposure_70_74 = pop_total_70_74 * annual_mean_pm25,
    exposure_75_79 = pop_total_75_79 * annual_mean_pm25,
    exposure_80_84 = pop_total_80_84 * annual_mean_pm25,
    exposure_85_89 = pop_total_85_89 * annual_mean_pm25,
    exposure_90_94 = pop_total_90_94 * annual_mean_pm25,
    exposure_95_99 = pop_total_95_99 * annual_mean_pm25,
    exposure_100plus = pop_total_100plus * annual_mean_pm25,
    
    exposure_65plus = pop_total_65plus * annual_mean_pm25,
    exposure_75plus = pop_total_75plus * annual_mean_pm25
  )

write.csv(age_stratified, paste0("da_age_stratified_exposure_", year, ".csv"), row.names = FALSE)
saveRDS(age_stratified, paste0("da_age_stratified_exposure_", year, ".rds"))

# Calculate BC-wide exposure by age group
bc_age_summary <- age_stratified %>%
  summarise(
    # Total populations
    pop_0_4 = sum(pop_total_0_4, na.rm = TRUE),
    pop_5_9 = sum(pop_total_5_9, na.rm = TRUE),
    pop_10_14 = sum(pop_total_10_14, na.rm = TRUE),
    pop_15_19 = sum(pop_total_15_19, na.rm = TRUE),
    pop_20_24 = sum(pop_total_20_24, na.rm = TRUE),
    pop_25_29 = sum(pop_total_25_29, na.rm = TRUE),
    pop_30_34 = sum(pop_total_30_34, na.rm = TRUE),
    pop_35_39 = sum(pop_total_35_39, na.rm = TRUE),
    pop_40_44 = sum(pop_total_40_44, na.rm = TRUE),
    pop_45_49 = sum(pop_total_45_49, na.rm = TRUE),
    pop_50_54 = sum(pop_total_50_54, na.rm = TRUE),
    pop_55_59 = sum(pop_total_55_59, na.rm = TRUE),
    pop_60_64 = sum(pop_total_60_64, na.rm = TRUE),
    pop_65_69 = sum(pop_total_65_69, na.rm = TRUE),
    pop_70_74 = sum(pop_total_70_74, na.rm = TRUE),
    pop_75_79 = sum(pop_total_75_79, na.rm = TRUE),
    pop_80_84 = sum(pop_total_80_84, na.rm = TRUE),
    pop_85_89 = sum(pop_total_85_89, na.rm = TRUE),
    pop_90_94 = sum(pop_total_90_94, na.rm = TRUE),
    pop_95_99 = sum(pop_total_95_99, na.rm = TRUE),
    pop_100plus = sum(pop_total_100plus, na.rm = TRUE),
    
    # Total exposures
    total_exposure_0_4 = sum(exposure_0_4, na.rm = TRUE),
    total_exposure_5_9 = sum(exposure_5_9, na.rm = TRUE),
    total_exposure_10_14 = sum(exposure_10_14, na.rm = TRUE),
    total_exposure_15_19 = sum(exposure_15_19, na.rm = TRUE),
    total_exposure_20_24 = sum(exposure_20_24, na.rm = TRUE),
    total_exposure_25_29 = sum(exposure_25_29, na.rm = TRUE),
    total_exposure_30_34 = sum(exposure_30_34, na.rm = TRUE),
    total_exposure_35_39 = sum(exposure_35_39, na.rm = TRUE),
    total_exposure_40_44 = sum(exposure_40_44, na.rm = TRUE),
    total_exposure_45_49 = sum(exposure_45_49, na.rm = TRUE),
    total_exposure_50_54 = sum(exposure_50_54, na.rm = TRUE),
    total_exposure_55_59 = sum(exposure_55_59, na.rm = TRUE),
    total_exposure_60_64 = sum(exposure_60_64, na.rm = TRUE),
    total_exposure_65_69 = sum(exposure_65_69, na.rm = TRUE),
    total_exposure_70_74 = sum(exposure_70_74, na.rm = TRUE),
    total_exposure_75_79 = sum(exposure_75_79, na.rm = TRUE),
    total_exposure_80_84 = sum(exposure_80_84, na.rm = TRUE),
    total_exposure_85_89 = sum(exposure_85_89, na.rm = TRUE),
    total_exposure_90_94 = sum(exposure_90_94, na.rm = TRUE),
    total_exposure_95_99 = sum(exposure_95_99, na.rm = TRUE),
    total_exposure_100plus = sum(exposure_100plus, na.rm = TRUE)
  ) %>%
  mutate(
    # Calculate population-weighted mean exposure for each age group
    mean_pm25_0_4 = total_exposure_0_4 / pop_0_4,
    mean_pm25_5_9 = total_exposure_5_9 / pop_5_9,
    mean_pm25_10_14 = total_exposure_10_14 / pop_10_14,
    mean_pm25_15_19 = total_exposure_15_19 / pop_15_19,
    mean_pm25_20_24 = total_exposure_20_24 / pop_20_24,
    mean_pm25_25_29 = total_exposure_25_29 / pop_25_29,
    mean_pm25_30_34 = total_exposure_30_34 / pop_30_34,
    mean_pm25_35_39 = total_exposure_35_39 / pop_35_39,
    mean_pm25_40_44 = total_exposure_40_44 / pop_40_44,
    mean_pm25_45_49 = total_exposure_45_49 / pop_45_49,
    mean_pm25_50_54 = total_exposure_50_54 / pop_50_54,
    mean_pm25_55_59 = total_exposure_55_59 / pop_55_59,
    mean_pm25_60_64 = total_exposure_60_64 / pop_60_64,
    mean_pm25_65_69 = total_exposure_65_69 / pop_65_69,
    mean_pm25_70_74 = total_exposure_70_74 / pop_70_74,
    mean_pm25_75_79 = total_exposure_75_79 / pop_75_79,
    mean_pm25_80_84 = total_exposure_80_84 / pop_80_84,
    mean_pm25_85_89 = total_exposure_85_89 / pop_85_89,
    mean_pm25_90_94 = total_exposure_90_94 / pop_90_94,
    mean_pm25_95_99 = total_exposure_95_99 / pop_95_99,
    mean_pm25_100plus = total_exposure_100plus / pop_100plus,
    
    year = year
  )

# Reshape for plotting
bc_age_long <- bc_age_summary %>%
  select(year, starts_with("mean_pm25_")) %>%
  pivot_longer(
    cols = starts_with("mean_pm25_"),
    names_to = "age_group",
    values_to = "mean_pm25"
  ) %>%
  mutate(
    age_group = gsub("mean_pm25_", "", age_group),
    age_group = gsub("_", "-", age_group)
  )

# Create visualization
if(!dir.exists("figures")) {
  dir.create("figures")
}

p_age <- ggplot(bc_age_long, aes(x = age_group, y = mean_pm25)) +
  geom_col(fill = "#2A9D8F") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = paste("Population-Weighted Mean PM2.5 Exposure by Age Group -", year),
    subtitle = "British Columbia",
    x = "Age Group",
    y = "Mean PM2.5 (µg/m³)"
  )

print(p_age)

ggsave(paste0("figures/age_stratified_pm25_", year, ".png"), 
       plot = p_age, width = 10, height = 8, dpi = 300)

# Save summary
write.csv(bc_age_summary, paste0("bc_age_stratified_summary_", year, ".csv"), row.names = FALSE)

# Fix invalid geometries
da_plot_clean <- st_make_valid(da_plot)

# Verify
sum(!st_is_valid(da_plot_clean))

# Annual mean PM2.5
p1 <- ggplot(data = da_plot_clean) +   
  geom_sf(aes(fill = annual_mean_pm25), color = NA, linewidth = 0) +
  scale_fill_viridis_c(option = "inferno", name = "PM2.5 (µg/m³)", na.value = "grey50") + 
  coord_sf(xlim = c(-140, -114), ylim = c(48, 60)) + 
  theme_minimal() +
  labs(title = paste("Annual Mean PM2.5 by Dissemination Area -", year),
       subtitle = "British Columbia") +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

print(p1)

ggsave(paste0("map_annual_mean_pm25_.png", year, ".png"),
       width = 10, height = 8, dpi = 300)

# Monthly fire season mean PM2.5

# Calculate monthly means for May- September
months <- c("May", "June", "July", "August", "September")
month_nums <- 5:9

calculate_monthly_mean <- function(file_df, month_num) {
  monthly_files <- file_df %>%
    filter(month(date) == month_num) %>%
    pull(filepath)
  
  if (length(monthly_files) == 0) {
    warning(paste("No files found for", year, month_num))
    return(NULL)
  }
  
  print(paste("Processing", length(monthly_files), "files for", month.name[month_num]))
  
  daily_means <- list()
  reference_rast <- NULL
  
  for (i in seq_along(monthly_files)) {
    r_full <- rast(monthly_files[i])
    r_mean <- r_full[[1]]
    
    if (is.null(reference_rast)) {
      reference_rast <- r_mean
    } else if (!all(as.vector(ext(r_mean)) == as.vector(ext(reference_rast)))) {
      message(paste("Extent mismatch:", basename(monthly_files[i]), "- resampling"))
      r_mean <- resample(r_mean, reference_rast, method = "bilinear")
    }
    
    daily_means[[i]] <- r_mean
  }
  
  daily_stack <- rast(daily_means)
  return(mean(daily_stack, na.rm = TRUE))
}

monthly_rasters <- list()
for (i in seq_along(month_nums)) {
  monthly_rasters[[months[i]]] <- calculate_monthly_mean(file_df, month_nums[i])
}

da_vect <- vect(da_plot)

raster_stack <- rast(lapply(months, function(m) monthly_rasters[[m]]))
names(raster_stack) <- tolower(months)

extracted <- terra::extract(
  raster_stack,
  da_vect, 
  fun = mean, 
  na.rm = TRUE, 
  weights = TRUE
)
for (i in seq_along(months)) {
  col_name <- paste0(tolower(months[i]), "_pm25")
  da_plot[[col_name]] <- extracted [, i + 1]
}

month_cols <- paste0(tolower(months), "_pm25")
all_values <- unlist(da_plot[month_cols])
scale_limits <- c(0, max(all_values, na.rm = TRUE))

plot_list <- lapply(seq_along(months), function(i) {
  ggplot() + 
    geom_sf(data = da_plot, aes(fill = .data[[month_cols[i]]]),
            color = NA, linewidth = 0) +
    scale_fill_viridis_c(option = "inferno",
                         name = "PM2.5(µg/m³)",
                         na.value = "grey50",
                         limits = scale_limits) +
    theme_minimal() +
    labs(title = months[i]) +
    theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 20, hjust = 0.5),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
})

combined_plot <- wrap_plots(plot_list, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

final_plot <- combined_plot +
  plot_annotation(
    title = paste("Fire Season PM2.5 by Month -", year),
    subtitle = "British Columbia"
  )

print(final_plot)

ggsave(paste0("montly_pm25_maps.png", year, ".png"),
       width = 20, height = 5, dpi = 300)

#Map 3: Days above 37.5 threshold
p3 <- ggplot() +
  geom_sf(data = da_plot, aes(fill = days_above_37.5), color = NA, linewidth = 0) +
  scale_fill_viridis_c(option = "inferno", name = "Days", na.value = "grey50") +
  theme_minimal() +
  labs(title = paste("Days Above 37.5 µg/m³ by Dissemination Area -", year),
       subtitle = "British Columbia") +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank())

print(p3)

ggsave(paste0("map_days_above_37.5_", year, ".png"),
       width = 10, height = 8, dpi = 300)

# Map 4: Days above 15 threshold
p4 <- ggplot() +
  geom_sf(data = da_plot, aes(fill = days_above_15), color = NA, linewidth = 0) +
  scale_fill_viridis_c(option = "inferno", name = "Days", na.value = "grey50") +
  theme_minimal() +
  labs(title = paste("Days Above 15 µg/m³ by Dissemination Area -", year),
       subtitle = "British Columbia") +
  theme(legend.position = "bottom",
        axis.text = element_blank(),
        axis.ticks = element_blank())

print(p4)

ggsave(paste0("map_days_above_15_", year, ".png"),
       width = 10, height = 8, dpi = 300)

# Ranked DAs by Exposure 

da_summary <- da_annual_summary %>%
  st_drop_geometry()

top_20_mean <- da_annual_summary %>%
  arrange(desc(annual_mean_pm25)) %>%
  head(20) %>%
  mutate(rank = row_number())

ggplot(top_20_mean, aes(x = reorder(DGUID, annual_mean_pm25), y = annual_mean_pm25)) +
  geom_col(fill = "brown3") +
  coord_flip() +
  theme_minimal() + 
  labs(title = paste("Top 20 Dissemination Areas by Annual Mean PM2.5 -", year),
       x = "Dissemination Area",
       y = "Annual Mean PM2.5 (µg/m³)")

ggsave(paste0("top20_mean_pm25_", year, ".png"),
       width = 10, height = 8, dpi = 300)

# Rank by Peak exposure 
top_20_max <- da_summary %>%
  arrange(desc(annual_max_pm25)) %>%
  head(20) %>%
  mutate(rank = row_number())

# Create bar chart - Peak
ggplot(top_20_max, aes(x = reorder(DGUID, annual_max_pm25), y = annual_max_pm25)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  theme_minimal() +
  labs(title = paste("Top 20 Dissemination Areas by Peak PM2.5 -", year),
       x = "Dissemination Area",
       y = "Peak PM2.5 (µg/m³)")

ggsave(paste0("top20_peak_pm25_", year, ".png"),
       width = 10, height = 8, dpi = 300)

# Rank by days above threshold (37.5)
top_20_days <- da_summary %>%
  arrange(desc(days_above_37.5)) %>%
  head(20) %>%
  mutate(rank = row_number())

# Create bar chart - days above theshold
ggplot(top_20_days, aes(x = reorder(DGUID, days_above_37.5), y = days_above_37.5)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  theme_minimal() +
  labs(title = paste("Top 20 Dissemination Areas by Days Above 37.5 µg/m³ -", year),
       x = "Dissemination Area",
       y = "Number of Days")

ggsave(paste0("top20_days_above_37.5_", year, ".png"),
       width = 10, height = 8, dpi = 300)
  
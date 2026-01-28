library(terra)
library(sf)
library(dplyr)
library(lubridate)

extract_pm25_comprehensive <- function(file_path, bc_boundary, year) {  
  tryCatch({
    # Load raster (51 hourly layers)
    raster <- rast(file_path)
    bc_transformed <- st_transform(bc_boundary, crs = crs(raster))
    raster_bc_all <- mask(crop(raster, bc_transformed), bc_transformed)
    
    # PEAK EXPOSURE (max across all 51 hours) 
    raster_max_51hr <- max(raster_bc_all, na.rm = TRUE)
    vals_max_51hr <- values(raster_max_51hr, na.rm = FALSE)
    vals_max_51hr_clean <- vals_max_51hr[!is.na(vals_max_51hr)]  
    
    # MEAN ACROSS 51 HOURS 
    raster_mean_51hr <- mean(raster_bc_all, na.rm = TRUE)
    vals_mean_51hr <- values(raster_mean_51hr, na.rm = FALSE)
    vals_mean_51hr_clean <- vals_mean_51hr[!is.na(vals_mean_51hr)]  
    
    # 24-HOUR METRICS 
    raster_24hr <- raster_bc_all[[1:24]]
    raster_max_24hr <- max(raster_24hr, na.rm = TRUE)
    raster_mean_24hr <- mean(raster_24hr, na.rm = TRUE)
    
    vals_max_24hr <- values(raster_max_24hr, na.rm = FALSE)
    vals_max_24hr_clean <- vals_max_24hr[!is.na(vals_max_24hr)]
    
    vals_mean_24hr <- values(raster_mean_24hr, na.rm = FALSE)
    vals_mean_24hr_clean <- vals_mean_24hr[!is.na(vals_mean_24hr)]
    
    # DURATION: Count hours exceeded thresholds 
    hours_over_15 <- app(raster_bc_all, fun = function(x) sum(x > 15, na.rm = TRUE))
    hours_over_37.5 <- app(raster_bc_all, fun = function(x) sum(x > 37.5, na.rm = TRUE))
    hours_over_50 <- app(raster_bc_all, fun = function(x) sum(x > 50, na.rm = TRUE))
    
    # BC-WIDE HOURLY EVOLUTION
    hourly_maxes <- global(raster_bc_all, "max", na.rm = TRUE)[,1]
    hourly_means <- global(raster_bc_all, "mean", na.rm = TRUE)[,1]
    
    # When peak occurred
    peak_hour <- which.max(hourly_maxes)
    
    # CONSECUTIVE HOURS OVER THRESHOLDS 
    # WHO AQG (15 µg/m³)
    event_15 <- hourly_maxes > 15 
    consecutive_15 <- if(any(event_15)) {
      rle_result <- rle(event_15)
      max(rle_result$lengths[rle_result$values], na.rm = TRUE)
    } else {0}
    
    # WHO IT-3 (37.5 µg/m³)
    event_37.5 <- hourly_maxes > 37.5 
    consecutive_37.5 <- if(any(event_37.5)) {
      rle_result <- rle(event_37.5)
      max(rle_result$lengths[rle_result$values], na.rm = TRUE)
    } else {0}
    
    # WHO IT-2 (50 µg/m³)
    event_50 <- hourly_maxes > 50
    consecutive_50 <- if(any(event_50)) {
      rle_result <- rle(event_50)
      max(rle_result$lengths[rle_result$values], na.rm = TRUE)
    } else {0}
    
    # SPATIAL EXTENT - PEAK values (51-hour max) 
    non_na_cells <- sum(!is.na(vals_max_51hr))
    cells_over_0 <- sum(vals_max_51hr_clean > 0, na.rm = TRUE)
    cells_over_10 <- sum(vals_max_51hr_clean > 10, na.rm = TRUE)
    cells_over_15 <- sum(vals_max_51hr_clean > 15, na.rm = TRUE)
    cells_over_25 <- sum(vals_max_51hr_clean > 25, na.rm = TRUE)
    cells_over_37.5 <- sum(vals_max_51hr_clean > 37.5, na.rm = TRUE)
    cells_over_50 <- sum(vals_max_51hr_clean > 50, na.rm = TRUE)
    cells_over_75 <- sum(vals_max_51hr_clean > 75, na.rm = TRUE)
    cells_over_100 <- sum(vals_max_51hr_clean > 100, na.rm = TRUE)
    cells_over_250 <- sum(vals_max_51hr_clean > 250, na.rm = TRUE)
    
    # SPATIAL EXTENT - MEAN values
    # 51-hour mean
    cells_mean_51hr_over_15 <- sum(vals_mean_51hr_clean > 15, na.rm = TRUE)  # FIXED
    cells_mean_51hr_over_37.5 <- sum(vals_mean_51hr_clean > 37.5, na.rm = TRUE)  # FIXED
    
    # 24-hour mean
    cells_mean_24hr_over_15 <- sum(vals_mean_24hr_clean > 15, na.rm = TRUE)  # FIXED
    cells_mean_24hr_over_37.5 <- sum(vals_mean_24hr_clean > 37.5, na.rm = TRUE)  # FIXED
    
    # DURATION-BASED SPATIAL METRICS
    cells_6hrs_over_15 <- sum(values(hours_over_15) >= 6, na.rm = TRUE)
    cells_12hrs_over_15 <- sum(values(hours_over_15) >= 12, na.rm = TRUE)
    cells_24hrs_over_15 <- sum(values(hours_over_15) >= 24, na.rm = TRUE)
    
    cells_6hrs_over_37.5 <- sum(values(hours_over_37.5) >= 6, na.rm = TRUE)
    cells_12hrs_over_37.5 <- sum(values(hours_over_37.5) >= 12, na.rm = TRUE)
    cells_24hrs_over_37.5 <- sum(values(hours_over_37.5) >= 24, na.rm = TRUE)
    
    # MEAN DURATION across affected cells 
    mean_hours_over_15 <- mean(values(hours_over_15)[values(hours_over_15) > 0], na.rm = TRUE)
    mean_hours_over_37.5 <- mean(values(hours_over_37.5)[values(hours_over_37.5) > 0], na.rm = TRUE)
    
    # DETERMINE RESOLUTION AND CALCULATE CELL AREA
    date_str <- gsub(".*_(\\d{8})_.*", "\\1", basename(file_path))  # FIXED regex
    file_date <- as.Date(date_str, format = "%Y%m%d")
    file_month <- month(file_date)
    
    # Assign resolution
    resolution_km <- if(year %in% c(2021, 2022, 2023)) {
      12
    } else if(year == 2024 && file_month %in% c(5, 6)) {
      12  # May-June 2024 exception
    } else if(year %in% c(2024, 2025, 2026)) {
      4
    } else {
      NA_real_
    }
    
    cell_area_km2 <- resolution_km^2
    
    # Data source label
    data_source <- if(year %in% c(2021, 2022, 2023)) {
      "12km (standard)"
    } else if(year == 2024 && file_month %in% c(5, 6)) {
      "12km (May-June gap-fill)"
    } else if(year %in% c(2024, 2025, 2026)) {
      "4km (standard)"
    } else {
      NA_character_
    }
    
    # === CONVERT CELL COUNTS TO AREA (km²) ===
    total_area_km2 <- non_na_cells * cell_area_km2
    
    area_over_0_km2 <- cells_over_0 * cell_area_km2
    area_over_10_km2 <- cells_over_10 * cell_area_km2
    area_over_15_km2 <- cells_over_15 * cell_area_km2
    area_over_25_km2 <- cells_over_25 * cell_area_km2
    area_over_37.5_km2 <- cells_over_37.5 * cell_area_km2
    area_over_50_km2 <- cells_over_50 * cell_area_km2
    area_over_75_km2 <- cells_over_75 * cell_area_km2
    area_over_100_km2 <- cells_over_100 * cell_area_km2
    area_over_250_km2 <- cells_over_250 * cell_area_km2
    
    area_mean_51hr_over_15_km2 <- cells_mean_51hr_over_15 * cell_area_km2
    area_mean_51hr_over_37.5_km2 <- cells_mean_51hr_over_37.5 * cell_area_km2
    area_mean_24hr_over_15_km2 <- cells_mean_24hr_over_15 * cell_area_km2
    area_mean_24hr_over_37.5_km2 <- cells_mean_24hr_over_37.5 * cell_area_km2
    
    area_6hrs_over_15_km2 <- cells_6hrs_over_15 * cell_area_km2
    area_12hrs_over_15_km2 <- cells_12hrs_over_15 * cell_area_km2
    area_24hrs_over_15_km2 <- cells_24hrs_over_15 * cell_area_km2
    area_6hrs_over_37.5_km2 <- cells_6hrs_over_37.5 * cell_area_km2
    area_12hrs_over_37.5_km2 <- cells_12hrs_over_37.5 * cell_area_km2
    area_24hrs_over_37.5_km2 <- cells_24hrs_over_37.5 * cell_area_km2
    
    # === COMPILE RESULTS ===
    data.frame(
      filename = basename(file_path),
      n_hours = as.integer(nlyr(raster)),
      
      # Resolution info
      resolution_km = as.numeric(resolution_km),
      cell_area_km2 = as.numeric(cell_area_km2),
      data_source = as.character(data_source),
      
      # PEAK METRICS (51-hour)
      peak_max_51hr = as.numeric(max(vals_max_51hr_clean, na.rm = TRUE)),
      peak_mean_51hr = as.numeric(mean(vals_max_51hr_clean, na.rm = TRUE)),
      peak_median_51hr = as.numeric(median(vals_max_51hr_clean, na.rm = TRUE)),
      peak_p95_51hr = as.numeric(quantile(vals_max_51hr_clean, 0.95, na.rm = TRUE)),
      
      # MEAN METRICS (51-hour)
      mean_mean_51hr = as.numeric(mean(vals_mean_51hr_clean, na.rm = TRUE)),
      mean_median_51hr = as.numeric(median(vals_mean_51hr_clean, na.rm = TRUE)),
      mean_max_51hr = as.numeric(max(vals_mean_51hr_clean, na.rm = TRUE)),
      
      # PEAK METRICS (24-hour)
      peak_max_24hr = as.numeric(max(vals_max_24hr_clean, na.rm = TRUE)),
      peak_mean_24hr = as.numeric(mean(vals_max_24hr_clean, na.rm = TRUE)),
      peak_median_24hr = as.numeric(median(vals_max_24hr_clean, na.rm = TRUE)),
      
      # MEAN METRICS (24-hour)
      mean_mean_24hr = as.numeric(mean(vals_mean_24hr_clean, na.rm = TRUE)),
      mean_median_24hr = as.numeric(median(vals_mean_24hr_clean, na.rm = TRUE)),
      mean_max_24hr = as.numeric(max(vals_mean_24hr_clean, na.rm = TRUE)),
      
      # TEMPORAL METRICS
      hour_of_peak = as.integer(peak_hour),
      consecutive_hrs_over_15 = as.integer(consecutive_15), 
      consecutive_hrs_over_37.5 = as.integer(consecutive_37.5),
      consecutive_hrs_over_50 = as.integer(consecutive_50),
      bc_wide_peak = as.numeric(max(hourly_maxes, na.rm = TRUE)),
      
      # SPATIAL EXTENT - CELLS
      total_cells = as.integer(non_na_cells),
      cells_any_pollution = as.integer(cells_over_0),
      cells_over_10 = as.integer(cells_over_10),
      cells_over_15_WHO_AQG = as.integer(cells_over_15),
      cells_over_25_WHO_IT4 = as.integer(cells_over_25),
      cells_over_37.5_WHO_IT3 = as.integer(cells_over_37.5),
      cells_over_50_WHO_IT2 = as.integer(cells_over_50),
      cells_over_75_WHO_IT1 = as.integer(cells_over_75),
      cells_over_100 = as.integer(cells_over_100),  # ADDED
      cells_over_250 = as.integer(cells_over_250),
      
      # SPATIAL EXTENT - AREA
      total_area_km2 = as.numeric(total_area_km2),
      area_over_0_km2 = as.numeric(area_over_0_km2),
      area_over_10_km2 = as.numeric(area_over_10_km2),
      area_over_15_km2 = as.numeric(area_over_15_km2),
      area_over_25_km2 = as.numeric(area_over_25_km2),
      area_over_37.5_km2 = as.numeric(area_over_37.5_km2),
      area_over_50_km2 = as.numeric(area_over_50_km2),
      area_over_75_km2 = as.numeric(area_over_75_km2),
      area_over_100_km2 = as.numeric(area_over_100_km2),
      area_over_250_km2 = as.numeric(area_over_250_km2),
      
      # MEAN CELLS
      cells_mean_51hr_over_15 = as.integer(cells_mean_51hr_over_15),
      cells_mean_51hr_over_37.5 = as.integer(cells_mean_51hr_over_37.5),
      cells_mean_24hr_over_15 = as.integer(cells_mean_24hr_over_15),
      cells_mean_24hr_over_37.5 = as.integer(cells_mean_24hr_over_37.5),
      
      # MEAN AREA
      area_mean_51hr_over_15_km2 = as.numeric(area_mean_51hr_over_15_km2),
      area_mean_51hr_over_37.5_km2 = as.numeric(area_mean_51hr_over_37.5_km2),
      area_mean_24hr_over_15_km2 = as.numeric(area_mean_24hr_over_15_km2),
      area_mean_24hr_over_37.5_km2 = as.numeric(area_mean_24hr_over_37.5_km2),
      
      # DURATION CELLS
      cells_6hrs_over_15 = as.integer(cells_6hrs_over_15),
      cells_12hrs_over_15 = as.integer(cells_12hrs_over_15),
      cells_24hrs_over_15 = as.integer(cells_24hrs_over_15),
      cells_6hrs_over_37.5 = as.integer(cells_6hrs_over_37.5),
      cells_12hrs_over_37.5 = as.integer(cells_12hrs_over_37.5),
      cells_24hrs_over_37.5 = as.integer(cells_24hrs_over_37.5),
      
      # DURATION AREA
      area_6hrs_over_15_km2 = as.numeric(area_6hrs_over_15_km2),
      area_12hrs_over_15_km2 = as.numeric(area_12hrs_over_15_km2),
      area_24hrs_over_15_km2 = as.numeric(area_24hrs_over_15_km2),
      area_6hrs_over_37.5_km2 = as.numeric(area_6hrs_over_37.5_km2),
      area_12hrs_over_37.5_km2 = as.numeric(area_12hrs_over_37.5_km2),
      area_24hrs_over_37.5_km2 = as.numeric(area_24hrs_over_37.5_km2),
      
      # MEAN DURATION 
      mean_hrs_over_15 = as.numeric(mean_hours_over_15),
      mean_hrs_over_37.5 = as.numeric(mean_hours_over_37.5),
      
      # PERCENTAGES
      pct_bc_over_15 = as.numeric((cells_over_15 / non_na_cells) * 100),
      pct_bc_over_37.5 = as.numeric((cells_over_37.5 / non_na_cells) * 100),
      pct_bc_over_50 = as.numeric((cells_over_50 / non_na_cells) * 100),
      pct_bc_over_100 = as.numeric((cells_over_100 / non_na_cells) * 100),
      
      success = TRUE,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    # Error handler - must match success output exactly
    data.frame(
      filename = basename(file_path),
      n_hours = NA_integer_,
      resolution_km = NA_real_,
      cell_area_km2 = NA_real_,
      data_source = NA_character_,
      
      # Peak 51hr
      peak_max_51hr = NA_real_, peak_mean_51hr = NA_real_, 
      peak_median_51hr = NA_real_, peak_p95_51hr = NA_real_,
      
      # Mean 51hr
      mean_mean_51hr = NA_real_, mean_median_51hr = NA_real_, mean_max_51hr = NA_real_,
      
      # Peak 24hr
      peak_max_24hr = NA_real_, peak_mean_24hr = NA_real_, peak_median_24hr = NA_real_,
      
      # Mean 24hr
      mean_mean_24hr = NA_real_, mean_median_24hr = NA_real_, mean_max_24hr = NA_real_,
      
      # Temporal
      hour_of_peak = NA_integer_,
      consecutive_hrs_over_15 = NA_integer_, 
      consecutive_hrs_over_37.5 = NA_integer_, 
      consecutive_hrs_over_50 = NA_integer_,
      bc_wide_peak = NA_real_, 
      
      # Cells
      total_cells = NA_integer_, 
      cells_any_pollution = NA_integer_,
      cells_over_10 = NA_integer_, 
      cells_over_15_WHO_AQG = NA_integer_, 
      cells_over_25_WHO_IT4 = NA_integer_,
      cells_over_37.5_WHO_IT3 = NA_integer_, 
      cells_over_50_WHO_IT2 = NA_integer_, 
      cells_over_75_WHO_IT1 = NA_integer_, 
      cells_over_100 = NA_integer_, 
      cells_over_250 = NA_integer_,
      
      # Area
      total_area_km2 = NA_real_,
      area_over_0_km2 = NA_real_,
      area_over_10_km2 = NA_real_,
      area_over_15_km2 = NA_real_,
      area_over_25_km2 = NA_real_,
      area_over_37.5_km2 = NA_real_,
      area_over_50_km2 = NA_real_,
      area_over_75_km2 = NA_real_,
      area_over_100_km2 = NA_real_,
      area_over_250_km2 = NA_real_,
      
      # Mean cells
      cells_mean_51hr_over_15 = NA_integer_, 
      cells_mean_51hr_over_37.5 = NA_integer_,
      cells_mean_24hr_over_15 = NA_integer_,
      cells_mean_24hr_over_37.5 = NA_integer_,
      
      # Mean area
      area_mean_51hr_over_15_km2 = NA_real_,
      area_mean_51hr_over_37.5_km2 = NA_real_,
      area_mean_24hr_over_15_km2 = NA_real_,
      area_mean_24hr_over_37.5_km2 = NA_real_,
      
      # Duration cells
      cells_6hrs_over_15 = NA_integer_, 
      cells_12hrs_over_15 = NA_integer_, 
      cells_24hrs_over_15 = NA_integer_,
      cells_6hrs_over_37.5 = NA_integer_, 
      cells_12hrs_over_37.5 = NA_integer_, 
      cells_24hrs_over_37.5 = NA_integer_,
      
      # Duration area
      area_6hrs_over_15_km2 = NA_real_,
      area_12hrs_over_15_km2 = NA_real_,
      area_24hrs_over_15_km2 = NA_real_,
      area_6hrs_over_37.5_km2 = NA_real_,
      area_12hrs_over_37.5_km2 = NA_real_,
      area_24hrs_over_37.5_km2 = NA_real_,
      
      # Mean duration
      mean_hrs_over_15 = NA_real_,
      mean_hrs_over_37.5 = NA_real_,
      
      # Percentages
      pct_bc_over_15 = NA_real_, 
      pct_bc_over_37.5 = NA_real_, 
      pct_bc_over_50 = NA_real_,
      pct_bc_over_100 = NA_real_,
      
      success = FALSE, 
      error = as.character(e$message),
      stringsAsFactors = FALSE
    )
  })
}

# ============================================
# MAIN PROCESSING SCRIPT
# ============================================

# SET YOUR YEAR
year <- 2025

# Get all .tif files in current directory
files <- list.files(pattern = "\\.tif$", full.names = TRUE)

print(paste("Processing year:", year))
print(paste("Found", length(files), "files"))

# Initialize results list
results_list <- list()

# Process each file with progress updates
start_time <- Sys.time()

for (i in seq_along(files)) {
  if (i %% 20 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    print(paste("Processed", i, "of", length(files), "files |",
                "Elapsed:", round(elapsed, 1), "mins"))
  }
  
  results_list[[i]] <- extract_pm25_comprehensive(files[i], bc_boundary, year)  
}

print("Combining results...")
results_df <- bind_rows(results_list)

# Add year and date
results_df$year <- year
results_df$date <- as.Date(gsub(".*_(\\d{8})_.*", "\\1", results_df$filename), 
                           format = "%Y%m%d")

# Save results
output_file <- paste0("bluesky_comprehensive_", year, ".csv")
write.csv(results_df, output_file, row.names = FALSE)

total_time <- difftime(Sys.time(), start_time, units = "mins")
print(paste("\nResults saved to:", output_file))
print(paste("Total processing time:", round(total_time, 1), "minutes"))

# summary statistics

print("\n========================================")
print("COMPREHENSIVE SUMMARY")
print("========================================")
print(paste("Total files processed:", nrow(results_df)))
print(paste("Successful:", sum(results_df$success)))
print(paste("Failed:", sum(!results_df$success)))

if(sum(results_df$success) > 0) {
  successful <- results_df %>% filter(success == TRUE)
  
  print("\n--- RESOLUTION INFO ---")
  resolution_summary <- successful %>%
    group_by(resolution_km, data_source) %>%
    summarise(
      n_days = n(),
      date_range = paste(min(date), "to", max(date)),
      .groups = "drop"
    )
  print(resolution_summary)
  
  print("\n--- PEAK EXPOSURE (24-hr) ---")
  print(paste("Highest PM2.5 peak:", round(max(successful$peak_max_24hr, na.rm = TRUE), 2), "µg/m³"))
  print(paste("Mean of daily peaks:", round(mean(successful$peak_max_24hr, na.rm = TRUE), 2), "µg/m³"))
  print(paste("Median of daily peaks:", round(median(successful$peak_max_24hr, na.rm = TRUE), 2), "µg/m³"))
  
  print("\n--- MEAN EXPOSURE (24-hr) ---")
  print(paste("Highest daily mean:", round(max(successful$mean_mean_24hr, na.rm = TRUE), 2), "µg/m³"))
  print(paste("Average daily mean:", round(mean(successful$mean_mean_24hr, na.rm = TRUE), 2), "µg/m³"))
  
  print("\n--- DURATION METRICS ---")
  print(paste("Mean consecutive hours > 15:", 
              round(mean(successful$consecutive_hrs_over_15, na.rm = TRUE), 1), "hours"))
  print(paste("Max consecutive hours > 15:", 
              max(successful$consecutive_hrs_over_15, na.rm = TRUE), "hours"))
  
  print(paste("Mean consecutive hours > 37.5:", 
              round(mean(successful$consecutive_hrs_over_37.5, na.rm = TRUE), 1), "hours"))
  print(paste("Max consecutive hours > 37.5:", 
              max(successful$consecutive_hrs_over_37.5, na.rm = TRUE), "hours"))
  
  print("\n--- EXPOSURE DAYS ---")
  print(paste("Days with peak > 15 (WHO AQG):", sum(successful$peak_max_24hr > 15, na.rm = TRUE)))
  print(paste("Days with peak > 37.5 (WHO IT-3):", sum(successful$peak_max_24hr > 37.5, na.rm = TRUE)))
  print(paste("Days with peak > 100:", sum(successful$peak_max_24hr > 100, na.rm = TRUE)))
  
  print("\n--- PROLONGED EXPOSURE ---")
  print(paste("Days with 6+ hours > 37.5:", sum(successful$cells_6hrs_over_37.5 > 0, na.rm = TRUE)))
  print(paste("Days with 12+ hours > 37.5:", sum(successful$cells_12hrs_over_37.5 > 0, na.rm = TRUE)))
  print(paste("Days with 24+ hours > 37.5:", sum(successful$cells_24hrs_over_37.5 > 0, na.rm = TRUE)))
  
  print("\n--- SPATIAL EXTENT (Percentage) ---")
  print(paste("Mean % of BC over WHO AQG (15):", 
              round(mean(successful$pct_bc_over_15, na.rm = TRUE), 1), "%"))
  print(paste("Mean % of BC over WHO IT-3 (37.5):", 
              round(mean(successful$pct_bc_over_37.5, na.rm = TRUE), 1), "%"))
  
  print("\n--- SPATIAL EXTENT (Area, resolution-normalized) ---")
  print(paste("Mean area over 15 µg/m³:", 
              round(mean(successful$area_over_15_km2, na.rm = TRUE), 0), "km²"))
  print(paste("Max area over 15 µg/m³:", 
              round(max(successful$area_over_15_km2, na.rm = TRUE), 0), "km²"))
  print(paste("Mean area over 37.5 µg/m³:", 
              round(mean(successful$area_over_37.5_km2, na.rm = TRUE), 0), "km²"))
  print(paste("Max area over 37.5 µg/m³:", 
              round(max(successful$area_over_37.5_km2, na.rm = TRUE), 0), "km²"))
  
  # Find worst days
  if(sum(results_df$success) > 0) {
    successful <- results_df %>% filter(success == TRUE)
    
    worst_days <- successful %>%
      arrange(desc(peak_max_24hr)) %>%  # FIXED: was peak_max
      select(date, peak_max_24hr, mean_mean_24hr, consecutive_hrs_over_37.5, 
             area_over_37.5_km2, pct_bc_over_37.5) %>%
      head(10)
    
    print(worst_days)
    
    # Save worst days
    worst_days_file <- paste0("worst_days_", year, ".csv")
    write.csv(worst_days, worst_days_file, row.names = FALSE)
    print(paste("\nWorst days saved to:", worst_days_file))
    
  }
 
}

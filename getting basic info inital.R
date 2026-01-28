library(terra)
library(sf)
library(dplyr)

# PM2.5 EXTRACTION SCRIPT
# peaks, duration (at the peak/jump) & and temporal patterns

extract_pm25_comprehensive <- function(file_path, bc_boundary) {
  tryCatch({
    # Load raster (51 hourly layers)
    raster <- rast(file_path)
    bc_transformed <- st_transform(bc_boundary, crs = crs(raster))
    raster_bc_all <- mask(crop(raster, bc_transformed), bc_transformed)
    
    #  Peak exposure (max across all 51 hours)
    raster_max_51hr <- max(raster_bc_all, na.rm = TRUE)
    vals_max_51hr <- values(raster_max_51hr, na.rm = FALSE)
    vals_max_51hr_clean <- vals_max[!is.na(vals_max_51hr)]
    
    # Mean across 51 hours 
    raster_mean_51hr <- mean(raster_bc_all, na.rm = TRUE)
    vals_mean_51hr <- values(raster_mean_51hr, na.rm = FALSE)
    vals_mean_51hr_51hr_clean <- vals_mean_51hr[!is.na(vals_mean_51hr)]
    
    # 24-hour metrics
    raster_24hr <- raster_bc_all[[1:24]]
    raster_max_24hr <- max(raster_24hr, na.rm = TRUE)
    raster_mean_24hr <- mean(raster_24hr, na.rm = TRUE)
    
    vals_max_24hr <- values(raster_max_24hr, na.rm = FALSE)
    vals_max_24hr_clean <- vals_max_24hr[!is.na(vals_max_24hr)]
    
    vals_mean_24hr <- values(raster_mean_24hr, na.rm = FALSE)
    vals_mean_24hr_clean <- vals_mean_24hr[!is.na(vals_mean_24hr)]
    
    # For each cell, count how many hours exceeded thresholds
    hours_over_15 <- app(raster_bc_all, fun = function(x) sum(x > 15, na.rm = TRUE))
    hours_over_37.5 <- app(raster_bc_all, fun = function(x) sum(x > 37.5, na.rm = TRUE))
    hours_over_50 <- app(raster_bc_all, fun = function(x) sum(x > 50, na.rm = TRUE))
    
    # BC-wide max for each hour to see temporal evolution
    hourly_maxes <- global(raster_bc_all, "max", na.rm = TRUE)[,1]
    hourly_means <- global(raster_bc_all, "mean", na.rm = TRUE)[,1]
    
    # when peak occurred
    peak_hour <- which.max(hourly_maxes)
    
    # Check duration at different thresholds
    # Who AGQ for 24 hrs (15 ug/m^3)
    event_15 <- hourly_maxes > 15 
    consecutive_15 <- if(any(event_15)) {
      rle_result <- rle(event_15)
      max(rle_result$lengths[rle_result$values], na.rm = TRUE)
    } else {0}
    
    # WHO IT-3 (37.5 ug/m^3)
    event_37.5 <- hourly_maxes > 37.5 
    consecutive_37.5 <- if(any(event_37.5)) {
      rle_result <- rle(event_37.5)
      max(rle_result$lengths[rle_result$values], na.rm = TRUE)
    } else {0}
    
    # WHO IT-2 (50 ug/m^3)
    event_50 <- hourly_maxes > 50
    consecutive_50 <- if(any(event_50)) {
      rle_result <- rle(event_50)
      max(rle_result$lengths[rle_result$values], na.rm = TRUE)
    } else {0}
    
    # Spatial extent - PEAK values with WHO thresholds
    non_na_cells <- sum(!is.na(vals_max_51hr))
    cells_over_0 <- sum(vals_max_51hr_clean > 0, na.rm = TRUE)
    cells_over_10 <- sum(vals_max_51hr_clean > 10, na.rm = TRUE)
    cells_over_15 <- sum(vals_max_51hr_clean > 15, na.rm = TRUE)    # WHO AQG
    cells_over_25 <- sum(vals_max_51hr_clean > 25, na.rm = TRUE)    # WHO IT-4
    cells_over_37.5 <- sum(vals_max_51hr_clean > 37.5, na.rm = TRUE)    # WHO IT-3
    cells_over_50 <- sum(vals_max_51hr_clean > 50, na.rm = TRUE)    # WHO IT-2
    cells_over_75 <- sum(vals_max_51hr_clean > 75, na.rm = TRUE)    # WHO IT-1
    cells_over_100 <- sum(vals_max_51hr_clean > 100, na.rm = TRUE)  # Very Unhealthy
    cells_over_250 <- sum(vals_max_51hr_clean > 250, na.rm = TRUE)  # Hazardous
    
    # Spatial extent - MEAN values (24-hr average equivalent)
    cells_mean_51hr_over_15 <- sum(vals_mean_51hr_clean > 15, na.rm = TRUE)  
    cells_mean_51_over_37.5 <- sum(vals_mean_51hr_clean > 37.5, na.rm = TRUE)  
    cells_mean_over_15 <- sum(vals_mean_24hr_clean > 15, na.rm = TRUE)  
    cells_mean_over_37.5 <- sum(vals_mean_24hr_clean > 37.5, na.rm = TRUE)
    
    # Duration-based spatial metrics
    cells_6hrs_over_15 <- sum(values(hours_over_15) >= 6, na.rm = TRUE)
    cells_12hrs_over_15 <- sum(values(hours_over_15) >= 12, na.rm = TRUE)
    cells_24hrs_over_15 <- sum(values(hours_over_15) >= 24, na.rm = TRUE)
    
    cells_6hrs_over_37.5 <- sum(values(hours_over_37.5) >= 6, na.rm = TRUE)
    cells_12hrs_over_37.5 <- sum(values(hours_over_37.5) >= 12, na.rm = TRUE)
    cells_24hrs_over_37.5 <- sum(values(hours_over_37.5) >= 24, na.rm = TRUE)
    
    # Mean duration across affected cells
    mean_hours_over_15 <- mean(values(
      hours_over_15)[values(hours_over_15) > 0], na.rm = TRUE)
    
    mean_hours_over_37.5 <- mean(values(
      hours_over_37.5)[values(hours_over_37.5) > 0], na.rm = TRUE)
    
    # Determine resolution and calculate respective geo area
    # check if May-June 2024!!
    date_str <- gsub(".*_(\\d{8}_.*", "\\1", basename(file_path))
    file_date <- as.Date(date_str, format = "%Y%m%d")
    file_month <- month(file_date)
    
    #assign resolution
    resolution_km <-if(year %in% c(2021, 2022, 2023)) {
      12
    } else if(year == 2024 && file_month %in% c(5, 6)) {
      4
    } else {
      NA_real_
    }
    
    cell_area_km2 <- resolution_km^2
    
    # Data source label
    data_source <- if(year %in% c(2021, 2022, 2023)) {
      "12km (standard)"
    } else if (year == 2024 && file_month %in% c(5, 6)) {
      "12km (May-June gap-fill)"
    } else if (year %in% c(2024, 2025, 2026)) {
      "4km (standard)"
    } else {
      NA_character_
    }
    
    # Convert Cell Counts to Area (km^2) 
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
    area_mean_51hr_over_37.5_km2 <- cells_mean_51_over_37.5 * cell_area_km2
    area_mean_24hr_over_15_km2 <- cells_mean_24hr_over_15 * cell_area_km2
    area_mean_24hr_over_37.5_km2 <- cells_mean_24hr_over_37.5 * cell_area_km2
    
    area_6hrs_over_15_km2 <- cells_6hrs_over_15 * cell_area_km2
    area_12hrs_over_15_km2 <- cells_12hrs_over_15 * cell_area_km2
    area_24hrs_over_15_km2 <- cells_24hrs_over_15 * cell_area_km2
    area_6hrs_over_37.5_km2 <- cells_6hrs_over_37.5 * cell_area_km2
    area_12hrs_over_37.5_km2 <- cells_12hrs_over_37.5 * cell_area_km2
    area_24hrs_over_37.5_km2 <- cells_24hrs_over_37.5 * cell_area_km2
    
    
    #compile results 
    data.frame(
      filename = basename(file_path),
      n_hours = nlyr(raster),
      
      resolution_km = as.numeric(resolution_km),
      cell_area_km2 = as.numeric(cell_area_km2),
      data_source = as.character(data_source),
      
      # PEAK METRICS (max across hours)
      peak_max_51hr = as.numeric(max(vals_max_51hr_clean, na.rm = TRUE)),
      peak_mean_51hr = as.numeric(mean(vals_max_51hr_clean, na.rm = TRUE)),
      peak_median_51hr = as.numeric(median(vals_max_51hr_clean, na.rm = TRUE)),
      peak_p95_51hr = as.numeric(quantile(vals_max_51hr_clean, 0.95, na.rm = TRUE)),
      
      # MEAN METRICS (average across hours - proxy for mean)
      peak_max_51hr = as.numeric(max(vals_max_51hr_clean, na.rm = TRUE)),
      mean_mean_51hr = as.numeric(mean(vals_mean_51hr_clean, na.rm = TRUE)),
      mean_median_51hr = as.numeric(median(vals_mean_51hr_clean, na.rm = TRUE)),
      
      peak_max_24hr = as.numeric(max(vals_max_24hr_clean, na.rm = TRUE)),
      mean_mean_24hr = as.numeric(mean(vals_mean_24hr_clean, na.rm = TRUE)),
      mean_median_24hr = as.numeric(median(vals_mean_24hr_clean, na.rm = TRUE)),
      
      # TEMPORAL METRICS
      hour_of_peak = as.integer(peak_hour),
      consecutive_hrs_over_15 = as.integer(consecutive_15), 
      consecutive_hrs_over_37.5 = as.integer(consecutive_37.5),
      consecutive_hrs_over_50 = as.integer(consecutive_50),
      bc_wide_peak = as.numeric(max(hourly_maxes, na.rm = TRUE)),
      
      # SPATIAL EXTENT - Peak exposure thresholds (raw cells)
      total_cells = as.integer(non_na_cells),
      cells_any_pollution = as.integer(cells_over_0),
      cells_over_10 = as.integer(cells_over_10),
      cells_over_15_WHO_AQG = as.integer(cells_over_15),
      cells_over_25_WHO_IT4 = as.integer(cells_over_25),
      cells_over_37.5_WHO_IT3 = as.integer(cells_over_37.5),
      cells_over_50_WHO_IT2 = as.integer(cells_over_50),
      cells_over_75_WHO_IT1 = as.integer(cells_over_75),
      cells_over_50 = as.integer(cells_over_50),
      cells_over_250 = as.integer(cells_over_250),
      
      # SPATIAL EXTENT - AREA (km², resolution-normalized)
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
      
      # SPATIAL EXTENT - Mean exposure (raw cells)
      cells_mean_over_15 = as.integer(cells_mean_over_15),
      cells_mean_over_37.5 = as.integer(cells_mean_over_37.5),
      
      # === SPATIAL EXTENT - MEAN (area) ===
      area_mean_51hr_over_15_km2 = as.numeric(area_mean_51hr_over_15_km2),
      area_mean_51hr_over_37.5_km2 = as.numeric(area_mean_51hr_over_37.5_km2),
      area_mean_24hr_over_15_km2 = as.numeric(area_mean_24hr_over_15_km2),
      area_mean_24hr_over_37.5_km2 = as.numeric(area_mean_24hr_over_37.5_km2),
      
      # DURATION + SPATIAL (prolonged exposure)
      cells_6hrs_over_15 = as.integer(cells_6hrs_over_15),
      cells_12hrs_over_15 = as.integer(cells_12hrs_over_15),
      cells_24hrs_over_15 = as.integer(cells_24hrs_over_15),
      cells_6hrs_over_37.5 = as.integer(cells_6hrs_over_37.5),
      cells_12hrs_over_37.5 = as.integer(cells_12hrs_over_37.5),
      cells_24hrs_over_37.5 = as.integer(cells_24hrs_over_37.5),
      
      # === DURATION + SPATIAL (area) ===
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
      
      success = TRUE,
      error = NA,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    # Return NA dataframe with same structure
    data.frame(
      filename = basename(file_path),
      n_hours = NA_integer_,
      resolution_km = NA_real_,
      cell_area_km2 = NA_real_,
      data_source = NA_character_,
      peak_max = NA_real_, peak_mean = NA_real_, peak_median = NA_real_, peak_p95 = NA_real_,
      mean_max = NA_real_, mean_mean = NA_real_, mean_median_51hr = NA_real_,
      hour_of_peak = NA_integer_,
      consecutive_hrs_over_15 = NA_integer_, 
      consecutive_hrs_over_37.5 = NA_integer_, 
      consecutive_hrs_over_50 = NA_integer_,
      bc_wide_peak = NA_real_, 
      total_cells = NA_integer_, 
      cells_any_pollution = NA_integer_,
      cells_over_10 = NA_integer_, 
      cells_over_15_WHO_AQG = NA_integer_, 
      cells_over_25_WHO_IT4 = NA_integer_,
      cells_over_37.5_WHO_IT3 = NA_integer_, 
      cells_over_50_WHO_IT2 = NA_integer_, 
      cells_over_75_WHO_IT1 = NA_integer_, 
      cells_over_50 = NA_integer_, 
      cells_over_250 = NA_integer_,
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
      cells_mean_51hr_over_15 = NA_integer_, 
      cells_mean_51hr_over_37.5 = NA_integer_,
      cells_mean_24hr_over_15 = NA_integer_,
      cells_mean_24hr_over_37.5 = NA_integer_,
      area_mean_51hr_over_15_km2 = NA_real_,
      area_mean_51hr_over_37.5_km2 = NA_real_,
      area_mean_24hr_over_15_km2 = NA_real_,
      area_mean_24hr_over_37.5_km2 = NA_real_,
      cells_6hrs_over_15 = NA_integer_, 
      cells_12hrs_over_15 = NA_integer_, 
      cells_24hrs_over_15 = NA_integer_,
      cells_6hrs_over_37.5 = NA_integer_, 
      cells_12hrs_over_37.5 = NA_integer_, 
      cells_24hrs_over_37.5 = NA_integer_,
      area_6hrs_over_15_km2 = NA_real_,
      area_12hrs_over_15_km2 = NA_real_,
      area_24hrs_over_15_km2 = NA_real_,
      area_6hrs_over_37.5_km2 = NA_real_,
      area_12hrs_over_37.5_km2 = NA_real_,
      area_24hrs_over_37.5_km2 = NA_real_,
      mean_hrs_over_15 = NA_real_,
      mean_hrs_over_37.5 = NA_real_,
      pct_bc_over_15 = NA_real_, 
      pct_bc_over_37.5 = NA_real_, 
      pct_bc_over_50 = NA_real_,
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
year <- 2026

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
    print(paste("Processed", i, "of", length(files), "files"))
  }
  
  results_list[[i]] <- extract_pm25_comprehensive(files[i], bc_boundary)
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

# SUMMARY STATISTICS

print("COMPREHENSIVE SUMMARY")
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
  
  print("\n--- PEAK EXPOSURE ---")
  print(paste("Highest PM2.5 peak 51hr:", round(max(successful$peak_max_51hr, na.rm = TRUE), 2), "µg/m³"))
  print(paste("Mean of daily peaks 51hr:", round(mean(successful$peak_max_51hr, na.rm = TRUE), 2), "µg/m³"))
  print(paste("Median of daily peaks 51hr:", round(median(successful$peak_max_51hr, na.rm = TRUE), 2), "µg/m³"))
  
  print(paste("Highest PM2.5 peak 24hr:", round(max(successful$peak_max_24hr, na.rm = TRUE), 2), "µg/m³"))
  print(paste("Mean of daily peaks 24hr:", round(mean(successful$peak_max_24hr, na.rm = TRUE), 2), "µg/m³"))
  
  print("\n--- MEAN EXPOSURE (24-hr equivalent) ---")
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
  
  print(paste("Mean consecutive hours > 50:", 
              round(mean(successful$consecutive_hrs_over_50, na.rm = TRUE), 1), "hours"))
  print(paste("Max consecutive hours > 50:", 
              max(successful$consecutive_hrs_over_50, na.rm = TRUE), "hours"))
  
  print("\n--- EXPOSURE DAYS ---")
  print(paste("Days with peak > 15 (WHO AQG):", sum(successful$peak_max_24hr > 15, na.rm = TRUE)))
  print(paste("Days with peak > 37.5 (WHO IT-3):", sum(successful$peak_max_24hr > 37.5, na.rm = TRUE)))
  print(paste("Days with peak > 100:", sum(successful$peak_max_24hr > 100, na.rm = TRUE)))
  print(paste("Days with peak > 250:", sum(successful$peak_max_24hr > 250, na.rm = TRUE)))
  
  print("\n--- PROLONGED EXPOSURE ---")
  print(paste("Days with 6+ hours > 37.5:", sum(successful$cells_6hrs_over_37.5 > 0, na.rm = TRUE)))
  print(paste("Days with 12+ hours > 37.5:", sum(successful$cells_12hrs_over_37.5 > 0, na.rm = TRUE)))
  print(paste("Days with 24+ hours > 37.5:", sum(successful$cells_24hrs_over_37.5 > 0, na.rm = TRUE)))
  
  print("\n--- SPATIAL EXTENT ---")
  print(paste("Mean % of BC over WHO AQG (15):", 
              round(mean(successful$pct_bc_over_15, na.rm = TRUE), 1), "%"))
  print(paste("Mean % of BC over WHO IT-3 (37.5):", 
              round(mean(successful$pct_bc_over_37.5, na.rm = TRUE), 1), "%"))
  print(paste("Mean % of BC over 50:", 
              round(mean(successful$pct_bc_over_50, na.rm = TRUE), 1), "%"))
  
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
  worst_days <- successful %>%
    arrange(desc(peak_max)) %>%
    select(date, peak_max_24hr, mean_mean_24hr, consecutive_hrs_over_37.5, 
           cells_over_37.5_km2, pct_bc_over_37.5) %>%
    head(10)
  
  print("\n--- TOP 10 WORST DAYS ---")
  print(worst_days)
  
  worst_days_file <- paste0("worst_days_", year, ".csv")
  write.csv(worst_days, worst_days_file, row.names = FALSE)
  print(paste("\nWorst days saved to:", worst_days_file))
}




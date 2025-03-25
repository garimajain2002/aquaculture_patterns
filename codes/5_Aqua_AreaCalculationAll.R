#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Calculate Aquaculture and dry aquaculture areas
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(raster)
library(terra)
library(dplyr)
library(ggplot2)
library(sp)
library(sf)
library(exactextractr)


# Define state and district codes 
# # ALL
state_codes <- c("OD", "AP", "TN")
districts <- list(
  OD = 1:6,
  AP = 1:9,
  TN = 1:12
)
years <- 2013:2025


# Create output directory if it doesn't exist
dir.create("outputs", showWarnings = FALSE)

# Function to process a single district for a specific year
process_district <- function(state_code, dist_code, year) {
  # Construct file paths
  boundary_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  aqua_path <- sprintf("data/tif/Landsat8_classified_dtype_compress/%s/%s_%d/%s_%d_classified_%d.tif", 
                       state_code, state_code, dist_code, state_code, dist_code, year)
 
  # Output file paths
  output_prefix <- sprintf("outputs/%s_%d_%d", state_code, dist_code, year)
  projected_aqua_path <- sprintf("%s_classified_projected.tif", output_prefix)
  aqua_area_path <- sprintf("%s_aqua_area_by_zone.csv", output_prefix)
  
  # Log processing start
  cat(sprintf("Processing %s_%d for year %d\n", state_code, dist_code, year))
  
  # Skip if output already exists (optional, for resuming interrupted runs)
  if (file.exists(projected_aqua_path)) {
    cat(sprintf("Output for %s_%d_%d already exists, skipping.\n", state_code, dist_code, year))
    return(FALSE)
  }
  
  # Check if input files exist
  if (!file.exists(boundary_path) || !file.exists(aqua_path)) {
    cat(sprintf("One or more input files for %s_%d_%d do not exist, skipping.\n", 
                state_code, dist_code, year))
    return(FALSE)
  }
  
  # Load data with error handling
  tryCatch({
    # Read input data
    boundaries <- st_read(boundary_path, quiet = TRUE)
    aqua_image <- brick(aqua_path)
    
    
    # 1. Calculate aqua area
    boundary_crs <- st_crs(boundaries)
    projected_raster <- projectRaster(aqua_image, crs = boundary_crs$wkt)
    
    # Create binary raster for aquaculture pixels
    aqua_raster <- projected_raster == 2
    
    # Create binary raster for dry aquaculture pixels
    dryaqua_raster <- projected_raster == 1
    
    
    # Calculate zonal statistics for aqua
    zonal_stats_aqua <- exact_extract(aqua_raster, boundaries, fun = c('sum'))
    boundaries$aqua_pixel_count <- zonal_stats_aqua
    boundaries$aqua_area_m2 <- boundaries$aqua_pixel_count * 30 * 30
    
    # Calculate zonal statistics for dry aqua
    zonal_stats_dryaqua <- exact_extract(dryaqua_raster, boundaries, fun = c('sum'))
    boundaries$dryaqua_pixel_count <- zonal_stats_dryaqua
    boundaries$dryaqua_area_m2 <- boundaries$dryaqua_pixel_count * 30 * 30
    
    # Save the results
    zonal_data <- st_drop_geometry(boundaries) %>%
      select(UniqueID, aqua_pixel_count, aqua_area_m2, dryaqua_pixel_count, dryaqua_area_m2)
    
    write.csv(zonal_data, aqua_area_path, row.names = FALSE)
    
    cat(sprintf("Successfully processed %s_%d for year %d\n", state_code, dist_code, year))
    return(TRUE)
    
  }, error = function(e) {
    # Log the error and continue with the next iteration
    cat(sprintf("Error processing %s_%d_%d: %s\n", state_code, dist_code, year, e$message))
    return(FALSE)
  })
}

# Initialize counters
total_combinations <- 0
successful_combinations <- 0

# Create a log file
log_file_aqua <- "outputs/processing_log.txt"
cat("Processing started at", format(Sys.time()), "\n", file = log_file_aqua)

# Main loop to process all combinations
for (state_code in state_codes) {
  for (dist_code in districts[[state_code]]) {
    for (year in years) {
      total_combinations <- total_combinations + 1
      result <- process_district(state_code, dist_code, year)
      if (result) successful_combinations <- successful_combinations + 1
      
      # Update log file after each combination
      cat(sprintf("Processed %s_%d_%d: %s\n", 
                  state_code, dist_code, year, 
                  ifelse(result, "SUCCESS", "FAILED")), 
          file = log_file, append = TRUE)
      
      # Optional: Add a small delay to prevent system overload
      Sys.sleep(1)
    }
  }
}

# Final summary
summary_message <- sprintf(
  "Processing complete!\nSuccessfully processed %d out of %d combinations (%.1f%%)",
  successful_combinations, total_combinations, 
  (successful_combinations/total_combinations) * 100
)

cat(summary_message, "\n")
cat("\n", summary_message, "\n", file = log_file_aqua, append = TRUE)
cat("Processing finished at", format(Sys.time()), "\n", file = log_file_aqua, append = TRUE)

# Create a summary CSV with all results
results_files <- list.files("outputs", pattern = ".*_aqua_area_by_zone\\.csv$", full.names = TRUE)
if (length(results_files) > 0) {
  all_results <- data.frame()
  
  for (file in results_files) {
    # Extract state, district, and year from filename
    file_parts <- basename(file)
    file_parts <- gsub("_aqua_area_by_zone.csv", "", file_parts)
    parts <- strsplit(file_parts, "_")[[1]]
    
    state <- parts[1]
    district <- as.numeric(parts[2])
    year <- as.numeric(parts[3])
    
    # Read the data
    data <- read.csv(file)
    
    # Add state, district, and year columns
    data$State <- state
    data$District <- district
    data$Year <- year
    
    # Append to the combined dataframe
    all_results <- rbind(all_results, data)
  }
  
  # Save the combined results
  write.csv(all_results, "outputs/all_districts_aquaculture_summary.csv", row.names = FALSE)
  cat("Created summary file with all results: outputs/all_districts_aquaculture_summary.csv\n")
  cat("Created summary file with all results: outputs/all_districts_aquaculture_summary.csv\n", 
      file = log_file_aqua, append = TRUE)
}


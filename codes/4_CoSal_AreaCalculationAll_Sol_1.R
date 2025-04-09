#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Model Salinity for all villages in all years and calculate areas
#+ Parallelized version for a multicore environment
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(raster)
library(terra)
library(dplyr)
library(ggplot2)
library(sp)
library(sf)
library(exactextractr)
library(parallel)
library(future)


setwd()

getwd()

# Get available cores and set number to use
total_cores <- future::availableCores()
# Use 90% of available cores to leave some resources for system operations
cores_to_use <- max(1, floor(total_cores * 0.9))
cat(sprintf("System has %d cores, using %d cores for processing\n", total_cores, cores_to_use))

# Load the ensemble model once (outside the loop)
MultiStratEnsemble <- readRDS("codes/ensemble_geomedian.rds")
best_threshold = 0.475

# Define state and district codes 
# ALL
# state_codes <- c("OD", "AP", "TN")
# districts <- list(
#   OD = 1:6,
#   AP = 1:9,
#   TN = 1:12
# )
# years <- 2013:2025

# # TEST
state_codes <- c("TN")
districts <- list(
  TN = 1:12
)
years <- 1990:2012

# Create output directory if it doesn't exist
dir.create("outputs", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# Function to process a single district for a specific year
process_district <- function(params) {
  state_code <- params$state_code
  dist_code <- params$dist_code
  year <- params$year
  
  # Redirect log output to a file specific to this task
  log_path <- sprintf("logs/%s_%d_%d.log", state_code, dist_code, year)
  log_con <- file(log_path, "w")
  
  # Helper function to write to both console and log file
  log_message <- function(msg) {
    message <- sprintf("[%s][%s_%d_%d] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
                       state_code, dist_code, year, msg)
    cat(message, "\n")
    cat(message, "\n", file = log_con)
  }
  
  # Construct file paths
  boundary_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  aqua_path <- sprintf("data/tif/Landsat5_classified_dtype_compress/%s_%d_classified_%d.tif", 
                       state_code, dist_code, year)
  landsat_path <- sprintf("data/tif/Landsat5_classified_dtype_compress/%s/%s_%d_geomedian_%d_predicted_predicted.tif", 
                          state_code, state_code, dist_code, year)
  
  # Output file paths
  output_prefix <- sprintf("outputs/%s_%d_%d", state_code, dist_code, year)
  final_raster_path <- sprintf("%s_predicted_ECAqua_raster.tif", output_prefix)
  projected_raster_path <- sprintf("%s_predicted_ECAqua_raster_projected.tif", output_prefix)
  saline_area_path <- sprintf("%s_saline_area_by_zone.csv", output_prefix)
  ec_map_path <- sprintf("%s_predicted_EC_map.png", output_prefix)
  #aqua_map_path <- sprintf("%s_predicted_Aqua_map.png", output_prefix)
  ecaqua_map_path <- sprintf("%s_predicted_ECAqua_map.png", output_prefix)
  
  # Log processing start
  log_message(sprintf("Processing started"))
  
  # Skip if output already exists (optional, for resuming interrupted runs)
  if (file.exists(final_raster_path)) {
    log_message(sprintf("Output already exists, skipping"))
    close(log_con)
    return(list(state_code = state_code, dist_code = dist_code, year = year, status = "SKIPPED"))
  }
  
  # Check if input files exist
  if (!file.exists(boundary_path) || !file.exists(aqua_path) || !file.exists(landsat_path)) {
    log_message(sprintf("One or more input files do not exist, skipping"))
    close(log_con)
    return(list(state_code = state_code, dist_code = dist_code, year = year, status = "MISSING_FILES"))
  }
  
  # Load data with error handling
  tryCatch({
    # Read input data
    log_message("Reading input data")
    boundaries <- st_read(boundary_path, quiet = TRUE)
    landsat_image <- brick(landsat_path)
    aqua_image <- brick(aqua_path)
    
    # 1. Prepare landsat multiband image
    log_message("Preparing landsat data")
    landsat_df <- as.data.frame(landsat_image, na.rm = FALSE)
    colnames(landsat_df) <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")
    
    # run for harmonised images - they need to be scaled back to SR 
    landsat_df$Blue <- (landsat_df$Blue * 0.0000275) - 0.2
    landsat_df$Red <- (landsat_df$Red * 0.0000275) - 0.2
    landsat_df$Green <- (landsat_df$Green * 0.0000275) - 0.2
    landsat_df$NIR <- (landsat_df$NIR * 0.0000275) - 0.2
    landsat_df$SWIR1 <- (landsat_df$SWIR1 * 0.0000275) - 0.2
    landsat_df$SWIR2 <- (landsat_df$SWIR2 * 0.0000275) - 0.2
    
    
    landsat_df$ID <- seq_len(nrow(landsat_df))
    
    landsat_df <- landsat_df[, c("ID", "Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")] 
    
    # Drop masked cells or where value is 0 (surface water cells)
    landsat_df <- subset(landsat_df, Blue != 0)
    
    # Ensure no NA values in landsat_df
    landsat_df <- landsat_df[complete.cases(landsat_df), ]
    
    # Calculate additional indices
    log_message("Calculating indices")
    # Blue and red 
    landsat_df$NBR <- (landsat_df$Blue-landsat_df$Red)/(landsat_df$Blue+landsat_df$Red)
    # Blue and green 
    landsat_df$NBG <- (landsat_df$Blue-landsat_df$Green)/(landsat_df$Blue+landsat_df$Green)
    # Blue and NIR
    landsat_df$NBNIR <- (landsat_df$Blue-landsat_df$NIR)/(landsat_df$Blue+landsat_df$NIR)
    # Blue and SWIR1
    landsat_df$NBSWIR1 <- (landsat_df$Blue-landsat_df$SWIR1)/(landsat_df$Blue+landsat_df$SWIR1)
    # Blue and SWIR2 
    landsat_df$NBSWIR2 <- (landsat_df$Blue-landsat_df$SWIR2)/(landsat_df$Blue+landsat_df$SWIR2)
    # Red and Green (also NDVI) 
    landsat_df$NDVI <- (landsat_df$Red-landsat_df$Green)/(landsat_df$Red+landsat_df$Green)
    # Red and NIR (NDSI2 or Normalised Difference Salinity Index 2)
    landsat_df$NDSI2 <- (landsat_df$Red-landsat_df$NIR)/(landsat_df$Red+landsat_df$NIR)
    # Red and SWIR1
    landsat_df$NRSWIR1 <- (landsat_df$Red-landsat_df$SWIR1)/(landsat_df$Red+landsat_df$SWIR1)
    # Red and SWIR2
    landsat_df$NRSWIR2 <- (landsat_df$Red-landsat_df$SWIR2)/(landsat_df$Red+landsat_df$SWIR2)
    # Green and NIR (also NDWI) 
    landsat_df$NDWI <- (landsat_df$Green-landsat_df$NIR)/(landsat_df$Green+landsat_df$NIR)
    # Green and SWIR1
    landsat_df$NGSWIR1 <- (landsat_df$Green-landsat_df$SWIR1)/(landsat_df$Green+landsat_df$SWIR1)
    # Green and SWIR2
    landsat_df$NGSWIR2 <- (landsat_df$Green-landsat_df$SWIR2)/(landsat_df$Green+landsat_df$SWIR2)
    # NIR and SWIR1
    landsat_df$NNIRSWIR1 <- (landsat_df$NIR-landsat_df$SWIR1)/(landsat_df$NIR+landsat_df$SWIR1)
    # NIR and SWIR2 
    landsat_df$NNIRSWIR2 <- (landsat_df$NIR-landsat_df$SWIR2)/(landsat_df$NIR+landsat_df$SWIR2)
    # SWIR1 and SWIR2 (also NDSI)
    landsat_df$NDSI1 <- (landsat_df$SWIR1-landsat_df$SWIR2)/(landsat_df$SWIR1+landsat_df$SWIR2)
    # Salinity Index 1 = sqrt(green^2+red^2)
    landsat_df$SI1 <- sqrt((landsat_df$Green)^2 + (landsat_df$Red)^2) 
    # Salinity Index 2 = sqrt(green x red)
    landsat_df$SI2 <- sqrt(landsat_df$Green * landsat_df$Red)
    # Salinity Index 3 = sqrt(blue x red) 
    landsat_df$SI3 <- sqrt(landsat_df$Blue * landsat_df$Red)
    # salinity index 4 = red x NIR / green 
    landsat_df$SI4 <- (landsat_df$Red * landsat_df$NIR / landsat_df$Green)  
    # salinity index 5 = blue/red 
    landsat_df$SI5 <- (landsat_df$Blue / landsat_df$Red)   
    # Soil Adjusted Vegetation Index (SAVI)
    landsat_df$SAVI <- (1.5 * landsat_df$NIR) - (0.5 * landsat_df$Red) + landsat_df$NIR + landsat_df$Red
    # Vegetation Soil Salinity Index (VSSI)
    landsat_df$VSSI <- (2 * landsat_df$Green) - 5 * (landsat_df$Red + landsat_df$NIR)
    
    landsat_df <- landsat_df[, c("ID", "Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2",
                                 "NDVI", "NDWI", "NDSI1", "NDSI2", "SI1", "SI2", "SI3", "SI4", "SI5", "SAVI", "VSSI",
                                 "NBR", "NBG", "NBNIR", "NBSWIR1", "NBSWIR2", "NRSWIR1", "NRSWIR2", "NGSWIR1", "NGSWIR2",
                                 "NNIRSWIR1", "NNIRSWIR2")]
    
    # Handle any Inf values
    landsat_df <- as.data.frame(lapply(landsat_df, function(x) {
      if (is.numeric(x)) x[is.infinite(x)] <- NA
      return(x)
    }))
    
    # Remove rows with NA values
    landsat_df <- na.omit(landsat_df)
    
    # Check if there's sufficient data after cleaning
    if (nrow(landsat_df) == 0) {
      log_message("No valid data after cleaning, skipping")
      close(log_con)
      return(list(state_code = state_code, dist_code = dist_code, year = year, status = "NO_DATA"))
    }
    
    # 2. Apply salinity model
    log_message("Applying salinity model")
    predictions <- predict(MultiStratEnsemble, newdata = landsat_df)
    landsat_df$predicted_EC <- ifelse(predictions[, "X1"] > best_threshold, 1, 0)
    
    # 3. Map predictions back to raster
    log_message("Creating prediction raster")
    predicted_raster <- raster(landsat_image[[1]])
    values(predicted_raster)[landsat_df$ID] <- landsat_df$predicted_EC
    
    # Create and save EC plot
    log_message("Creating EC plot")
    predicted_df <- as.data.frame(predicted_raster, xy = TRUE)
    colnames(predicted_df) <- c("x", "y", "predicted_EC")
    
    predicted_EC <- ggplot(predicted_df, aes(x = x, y = y, fill = factor(predicted_EC))) +
      geom_raster() +
      scale_fill_manual(values = c("0" = "black", "1" = "white")) +
      theme_minimal() +
      labs(title = sprintf("Predicted Electrical Conductivity - %s_%d_%d", 
                           state_code, dist_code, year), 
           fill = "Salinity Class")
    
    ggsave(ec_map_path, plot = predicted_EC, width = 8, height = 6, dpi = 300)
    
    # 4. Mask Aquaculture ponds
    log_message("Processing aquaculture data")
    aqua_df <- as.data.frame(aqua_image, xy = TRUE)
    colnames(aqua_df) <- c("x", "y", "classification")
    
    # Create mask and apply to predictions
    aqua_mask_raster <- aqua_image == 2  # Binary mask of aquaculture areas
    predicted_raster[aqua_mask_raster] <- 2  # Mark aquaculture areas
    values(predicted_raster) <- as.numeric(values(predicted_raster))
    
    # 5. Plot and save final maps
    log_message("Creating final maps")
    final_df <- as.data.frame(predicted_raster, xy = TRUE)
    colnames(final_df) <- c("x", "y", "predicted_EC")
    
    predicted_ECAqua <- ggplot(final_df, aes(x = x, y = y, fill = factor(predicted_EC))) +
      geom_raster() +
      scale_fill_manual(values = c("0" = "black", "1" = "orange", "2" = "blue")) +
      theme_minimal() +
      labs(title = sprintf("Predicted EC in Aquaculture context - %s_%d_%d", 
                           state_code, dist_code, year), 
           fill = "Salinity Class")
    
    ggsave(ecaqua_map_path, plot = predicted_ECAqua, width = 8, height = 6, dpi = 300)
    writeRaster(predicted_raster, final_raster_path, format = "GTiff", overwrite = TRUE)
    
    # 6. Calculate salinity area
    log_message("Calculating salinity area")
    boundary_crs <- st_crs(boundaries)
    projected_raster <- projectRaster(predicted_raster, crs = boundary_crs$wkt)
    writeRaster(projected_raster, projected_raster_path, format = "GTiff", overwrite = TRUE)
    
    # Create binary raster for saline pixels
    saline_raster <- projected_raster == 1
    
    # Calculate zonal statistics
    zonal_stats <- exact_extract(saline_raster, boundaries, fun = c('sum'))
    boundaries$saline_pixel_count <- zonal_stats
    boundaries$saline_area_m2 <- boundaries$saline_pixel_count * 30 * 30
    
    # Save the results
    zonal_data <- st_drop_geometry(boundaries) %>%
      select(UniqueID, saline_pixel_count, saline_area_m2)
    
    write.csv(zonal_data, saline_area_path, row.names = FALSE)
    
    log_message("Processing completed successfully")
    close(log_con)
    return(list(state_code = state_code, dist_code = dist_code, year = year, status = "SUCCESS"))
    
  }, error = function(e) {
    # Log the error and continue with the next iteration
    log_message(sprintf("Error: %s", e$message))
    close(log_con)
    return(list(state_code = state_code, dist_code = dist_code, year = year, 
                status = "ERROR", error_message = e$message))
  })
}

# Create parameter list for all combinations
all_params <- list()
for (state_code in state_codes) {
  for (dist_code in districts[[state_code]]) {
    for (year in years) {
      all_params <- c(all_params, list(list(
        state_code = state_code,
        dist_code = dist_code,
        year = year
      )))
    }
  }
}

# Create a master log file
master_log_file <- "master_processing_log.txt"
cat("Processing started at", format(Sys.time()), "\n", 
    "Using", cores_to_use, "cores out of", total_cores, "available\n",
    "Total tasks:", length(all_params), "\n",
    file = master_log_file)

# Set memory limits for parallel processing
options(mc.cores = cores_to_use)
# Note: For very large rasters, you might need to adjust the memory limits
# options(future.globals.maxSize = 1000 * 1024^2)  # Set limit to 1GB

all_params <- all_params[1:5]

print(all_params)

# Process all combinations in parallel
cat("Starting parallel processing with", cores_to_use, "cores for", length(all_params), "tasks\n")
start_time <- Sys.time()

# Use mclapply for parallel processing on Unix/Linux/Mac
if (.Platform$OS.type == "unix") {
  results <- mclapply(all_params, process_district, mc.cores = cores_to_use)
} else {
  # For Windows, use parLapply instead
  cl <- makeCluster(cores_to_use)
  clusterExport(cl, c("MultiStratEnsemble", "best_threshold"))
  clusterEvalQ(cl, {
    library(raster)
    library(terra)
    library(dplyr)
    library(ggplot2)
    library(sp)
    library(sf)
    library(exactextractr)
  })
  
  results <- parLapply(cl, all_params, fun = process_district)
  stopCluster(cl)
}

end_time <- Sys.time()
processing_time <- difftime(end_time, start_time, units = "mins")

# Count successes and failures
statuses <- sapply(results, function(x) x$status)
successes <- sum(statuses == "SUCCESS")
skipped <- sum(statuses == "SKIPPED")
errors <- sum(statuses %in% c("ERROR", "MISSING_FILES", "NO_DATA"))

# Final summary
summary_message <- sprintf(
  "Processing complete in %.1f minutes!\n%d successful, %d skipped, %d failed out of %d total tasks (%.1f%% success rate)",
  as.numeric(processing_time),
  successes, skipped, errors, length(all_params),
  (successes/length(all_params)) * 100
)

cat(summary_message, "\n")
cat("\n", summary_message, "\n", file = master_log_file, append = TRUE)
cat("Processing finished at", format(Sys.time()), "\n\n", file = master_log_file, append = TRUE)

# Save detailed results
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    state_code = x$state_code,
    dist_code = x$dist_code,
    year = x$year,
    status = x$status,
    error_message = ifelse(is.null(x$error_message), NA, x$error_message),
    stringsAsFactors = FALSE
  )
}))

write.csv(results_df, "task_results_summary.csv", row.names = FALSE)
cat("Detailed processing results saved to: task_results_summary.csv\n")
cat("Detailed processing results saved to: task_results_summary.csv\n", 
    file = master_log_file, append = TRUE)

# Create a summary CSV with all results
results_files <- list.files("outputs", pattern = ".*_saline_area_by_zone\\.csv$", full.names = TRUE)
if (length(results_files) > 0) {
  cat("Aggregating result files...\n")
  all_results <- data.frame()
  
  for (file in results_files) {
    # Extract state, district, and year from filename
    file_parts <- basename(file)
    file_parts <- gsub("_saline_area_by_zone.csv", "", file_parts)
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
  write.csv(all_results, "outputs/all_districts_salinity_summary.csv", row.names = FALSE)
  cat("Created summary file with all results: outputs/all_districts_salinity_summary.csv\n")
  cat("Created summary file with all results: outputs/all_districts_salinity_summary.csv\n", 
      file = master_log_file, append = TRUE)
}

# Display memory usage information
gc_stats <- gc()
mem_used <- sum(gc_stats[,2]) * 1024^2  # in MB
cat(sprintf("Peak memory usage: %.2f MB\n", mem_used))
cat(sprintf("Peak memory usage: %.2f MB\n", mem_used), file = master_log_file, append = TRUE)

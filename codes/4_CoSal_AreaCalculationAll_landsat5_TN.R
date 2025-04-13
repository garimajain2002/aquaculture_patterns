#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+ Model Salinity for all villages in all years and calculate areas
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

library(caret)
library(caretEnsemble)
library(pROC)
library(data.table)
library(scales)
library(glmnet)
library(xgboost)
library(tidyverse)
library(gtsummary) # For summary tables
library(modelsummary)# For summary tables
library(mgcv) # GAM model fit 
library(randomForest) # to apply machine learning frameworks
library(datawizard) # for normalize()
library(nnet) # For ANN
library(neuralnet) # For more control on the architecture of ANN
library(glmnet) # For lasso regression 
library(MASS) # for stepwise regression
library(scales) # for scaling data from 0 to 1
library(kernlab)
library(naivebayes)



getwd()

# Load the ensemble model once (outside the loop)
MultiStratEnsemble <- readRDS("ensemble_geomedian.rds")
best_threshold = 0.475

# Define state and district codes 
# # ALL
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
  TN = 7:11
)
years <- 1990:2001

# Create output directory if it doesn't exist
dir.create("outputs", showWarnings = FALSE)

# Function to process a single district for a specific year
process_district <- function(state_code, dist_code, year) {
  # Construct file paths
  boundary_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  aqua_path <- sprintf("data/tif/Landsat5_classified_dtype_compress/%s_%d_classified_%d.tif", 
                       state_code, dist_code, year)
  landsat_path <- sprintf("G:/My Drive/GEE_Exports/Landsat5_harmonised_dtype_compress/%s/%s_%d_geomedian_%d_predicted_predicted.tif", 
                          state_code, state_code, dist_code, year)
  
  # Output file paths
  output_prefix <- sprintf("outputs/%s_%d_%d", state_code, dist_code, year)
  final_raster_path <- sprintf("%s_predicted_ECAqua_raster.tif", output_prefix)
  projected_raster_path <- sprintf("%s_predicted_ECAqua_raster_projected.tif", output_prefix)
  saline_area_path <- sprintf("%s_saline_area_by_zone.csv", output_prefix)
  ec_map_path <- sprintf("%s_predicted_EC_map.png", output_prefix)
  aqua_map_path <- sprintf("%s_predicted_Aqua_map.png", output_prefix)
  ecaqua_map_path <- sprintf("%s_predicted_ECAqua_map.png", output_prefix)
  
  # Log processing start
  cat(sprintf("Processing %s_%d for year %d\n", state_code, dist_code, year))
  
  # Skip if output already exists (optional, for resuming interrupted runs)
  if (file.exists(final_raster_path)) {
    cat(sprintf("Output for %s_%d_%d already exists, skipping.\n", state_code, dist_code, year))
    return(FALSE)
  }
  
  # Check if input files exist
  if (!file.exists(boundary_path) || !file.exists(aqua_path) || !file.exists(landsat_path)) {
    cat(sprintf("One or more input files for %s_%d_%d do not exist, skipping.\n", 
                state_code, dist_code, year))
    return(FALSE)
  }
  
  # Load data with error handling
  tryCatch({
    # Read input data
    boundaries <- st_read(boundary_path, quiet = TRUE)
    landsat_image <- brick(landsat_path)
    aqua_image <- brick(aqua_path)
    
    # 1. Prepare landsat multiband image
    landsat_df <- as.data.frame(landsat_image, na.rm = FALSE)
    colnames(landsat_df) <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")
    
    # run only for harmonised images - they need to be scaled back to SR 
    landsat_df$Blue <- (landsat_df$Blue * 0.0000275) - 0.2
    landsat_df$Red <- (landsat_df$Red * 0.0000275) - 0.2
    landsat_df$Green <- (landsat_df$Green * 0.0000275) - 0.2
    landsat_df$NIR <- (landsat_df$NIR * 0.0000275) - 0.2
    landsat_df$SWIR1 <- (landsat_df$SWIR1 * 0.0000275) - 0.2
    landsat_df$SWIR2 <- (landsat_df$SWIR2 * 0.0000275) - 0.2
    
    # # Check if the landsat image has the expected number of bands
    # expected_bands <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2", "NDVI", "NDWI", "NDSI")
    # if (ncol(landsat_df) != length(expected_bands)) {
    #   cat(sprintf("Warning: Landsat image for %s_%d_%d has %d bands, expected %d. Skipping.\n", 
    #               state_code, dist_code, year, ncol(landsat_df), length(expected_bands)))
    #   return(FALSE)
    # }
    # 
    # colnames(landsat_df) <- expected_bands
    
    landsat_df$ID <- seq_len(nrow(landsat_df))
    
    landsat_df <- landsat_df[, c("ID", "Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")] 
    
    # Drop masked cells or where value is 0 (surface water cells)
    landsat_df <- subset(landsat_df, Blue != 0)
    
    # Ensure no NA values in landsat_df
    landsat_df <- landsat_df[complete.cases(landsat_df), ]
    
    # Calculate additional indices
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
      cat(sprintf("No valid data after cleaning for %s_%d_%d, skipping.\n", 
                  state_code, dist_code, year))
      return(FALSE)
    }
    
    # 2. Apply salinity model
    predictions <- predict(MultiStratEnsemble, newdata = landsat_df)
    landsat_df$predicted_EC <- ifelse(predictions[, "X1"] > best_threshold, 1, 0)
    
    # 3. Map predictions back to raster
    predicted_raster <- raster(landsat_image[[1]])
    values(predicted_raster)[landsat_df$ID] <- landsat_df$predicted_EC
    
    # Create and save EC plot
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
    aqua_df <- as.data.frame(aqua_image, xy = TRUE)
    colnames(aqua_df) <- c("x", "y", "classification")
    
    # Plot aquaculture
    predicted_Aqua <- ggplot(aqua_df, aes(x = x, y = y, fill = factor(classification))) +
      geom_raster() +
      scale_fill_manual(values = c("0" = "black", "1" = "black", "2" = "blue")) +
      theme_minimal() +
      labs(title = sprintf("Classified Aquaculture ponds - %s_%d_%d", 
                           state_code, dist_code, year), 
           fill = "Aquaculture ponds")
    
    #ggsave(aqua_map_path, plot = predicted_Aqua, width = 8, height = 6, dpi = 300)
    
    # Create mask and apply to predictions
    aqua_mask_raster <- aqua_image == 2  # Binary mask of aquaculture areas
    predicted_raster[aqua_mask_raster] <- 2  # Mark aquaculture areas
    values(predicted_raster) <- as.numeric(values(predicted_raster))
    
    # 5. Plot and save final maps
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
      dplyr::select(UniqueID, saline_pixel_count, saline_area_m2)
    
    write.csv(zonal_data, saline_area_path, row.names = FALSE)
    
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
log_file <- "processing_log_TN.txt"
cat("Processing started at", format(Sys.time()), "\n", file = log_file)

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
cat("\n", summary_message, "\n", file = log_file, append = TRUE)
cat("Processing finished at", format(Sys.time()), "\n", file = log_file, append = TRUE)

# Create a summary CSV with all results
results_files <- list.files("outputs", pattern = ".*_saline_area_by_zone\\.csv$", full.names = TRUE)
if (length(results_files) > 0) {
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
      file = log_file, append = TRUE)
}


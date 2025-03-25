# Calculate % agriculture and % urban per village for a ref. year (2015) 

library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(purrr)

getwd()

# Define state codes and districts
state_codes <- c("OD", "AP", "TN")
districts <- list(
  OD = 1:6,
  AP = 1:9,
  TN = 1:12
)

# Function to calculate percentage of agriculture and urban area for each village
calculate_lulc_percentages <- function(state_code, dist_code) {
  # Construct file paths
  boundary_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  lulc_path <- sprintf("data/tif/LandCover_dtype_compress/%s_%d_LC.tif", state_code, dist_code)
  
  # Read shapefile and LULC raster
  cat(sprintf("Processing %s district %d...\n", state_code, dist_code))
  
  # Error handling for file reading
  tryCatch({
    villages <- sf::st_read(boundary_path, quiet = TRUE)
    lulc <- raster::raster(lulc_path)
    
    # Check and align projections
    cat(sprintf("Checking projections for %s district %d...\n", state_code, dist_code))
    villages_crs <- st_crs(villages)
    lulc_crs <- crs(lulc)
    
    # Print CRS information for debugging
    cat(sprintf("Village CRS: %s\n", as.character(villages_crs$wkt)))
    cat(sprintf("LULC CRS: %s\n", as.character(lulc_crs)))
    
    # Align projections if they don't match
    if (!compareCRS(villages_crs, lulc_crs)) {
      cat("Projections don't match. Reprojecting LULC to match village boundaries...\n")
      lulc <- projectRaster(lulc, crs = villages_crs$wkt)
      cat("Reprojection complete.\n")
    } else {
      cat("Projections match. No reprojection needed.\n")
    }
    
    # Create an empty dataframe to store results
    result <- data.frame(
      state_code = character(),
      district_code = integer(),
      UniqueID = character(),
      total_area_ha = numeric(),
      agriculture_area_ha = numeric(),
      urban_area_ha = numeric(),
      agriculture_percent = numeric(),
      urban_percent = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Process each village
    for (i in 1:nrow(villages)) {
      village <- villages[i, ]
      
      # Extract UniqueID
      if ("UniqueID" %in% names(village)) {
        unique_id <- as.character(village$UniqueID)
      } else if ("UNIQUE_ID" %in% names(village)) {
        unique_id <- as.character(village$UNIQUE_ID)
      } else if ("unique_id" %in% names(village)) {
        unique_id <- as.character(village$unique_id)
      } else {
        # If no specific UniqueID column found, generate one
        unique_id <- paste0(state_code, "_", dist_code, "_", i)
        cat(sprintf("Warning: No UniqueID found for village %d in %s district %d. Generated: %s\n", 
                    i, state_code, dist_code, unique_id))
      }
      
      # Extract village boundary
      village_boundary <- village
      
      # Calculate total area in hectares
      total_area_ha <- sf::st_area(village_boundary) / 10000
      
      # Check if village has valid geometry
      if (sf::st_is_valid(village_boundary)) {
        # Crop LULC raster to village boundary
        village_lulc <- tryCatch({
          # Get spatial extent of village with a small buffer
          ext <- extent(as(sf::st_buffer(village_boundary, 10), "Spatial"))
          raster::crop(lulc, ext)
        }, error = function(e) {
          cat(sprintf("Error cropping raster for village %s: %s\n", unique_id, e$message))
          return(NULL)
        })
        
        # Mask LULC raster with village boundary
        if (!is.null(village_lulc)) {
          village_lulc <- tryCatch({
            raster::mask(village_lulc, as(village_boundary, "Spatial"))
          }, error = function(e) {
            cat(sprintf("Error masking raster for village %s: %s\n", unique_id, e$message))
            return(NULL)
          })
        }
        
        # Calculate areas
        if (!is.null(village_lulc)) {
          # Count pixels by land use type
          lulc_counts <- raster::freq(village_lulc)
          
          # Extract agriculture (40) and urban (50) areas
          agriculture_pixels <- lulc_counts[lulc_counts[, 1] == 40, 2]
          urban_pixels <- lulc_counts[lulc_counts[, 1] == 50, 2]
          
          # Handle NA values
          agriculture_pixels <- ifelse(length(agriculture_pixels) == 0, 0, agriculture_pixels)
          urban_pixels <- ifelse(length(urban_pixels) == 0, 0, urban_pixels)
          
          # Calculate total pixels
          total_pixels <- sum(lulc_counts[, 2], na.rm = TRUE)
          
          # Fixed pixel size: 100m x 100m = 1 hectare
          pixel_area_ha <- 1  # Each pixel is exactly 1 hectare (100m x 100m)
          agriculture_area_ha <- agriculture_pixels * pixel_area_ha
          urban_area_ha <- urban_pixels * pixel_area_ha
          
          # Calculate percentages
          agriculture_percent <- (agriculture_area_ha / as.numeric(total_area_ha)) * 100
          urban_percent <- (urban_area_ha / as.numeric(total_area_ha)) * 100
        } else {
          agriculture_area_ha <- 0
          urban_area_ha <- 0
          agriculture_percent <- 0
          urban_percent <- 0
        }
      } else {
        # Invalid geometry
        cat(sprintf("Warning: Invalid geometry for village %s. Skipping...\n", unique_id))
        agriculture_area_ha <- 0
        urban_area_ha <- 0
        agriculture_percent <- 0
        urban_percent <- 0
      }
      
      # Add results to dataframe
      result <- bind_rows(result, data.frame(
        state_code = state_code,
        district_code = dist_code,
        UniqueID = unique_id,
        total_area_ha = as.numeric(total_area_ha),
        agriculture_area_ha = agriculture_area_ha,
        urban_area_ha = urban_area_ha,
        agriculture_percent = agriculture_percent,
        urban_percent = urban_percent,
        stringsAsFactors = FALSE
      ))
    }
    
    return(result)
  }, error = function(e) {
    cat(sprintf("Error processing %s district %d: %s\n", state_code, dist_code, e$message))
    return(NULL)
  })
}

# Process all districts for all states
results_list <- list()

for (state_code in state_codes) {
  for (dist_code in districts[[state_code]]) {
    result <- calculate_lulc_percentages(state_code, dist_code)
    if (!is.null(result)) {
      print(table(duplicated(result$UniqueID)))  # Check for duplicates per district
      results_list[[length(results_list) + 1]] <- result
    }
  }
}


# Combine all results
all_results <- do.call(rbind, results_list)


# Check and remove duplicates if any 
print("Data dimensions:")
print(dim(all_results)) # Should be 31088 
print(colnames(all_results))

table(duplicated(all_results$UniqueID))  # TRUE means duplicates exist
# ! check why
# select the first value only
all_results <- all_results %>% distinct(UniqueID, .keep_all = TRUE)

print("Data dimensions:")
print(dim(all_results)) # Should be 31088


# Save results to CSV
write.csv(all_results_unique, "outputs/village_lulc_percentages.csv", row.names = FALSE)

# Print summary statistics
cat("\nSummary of results:\n")
summary_stats <- all_results %>%
  group_by(state_code, district_code) %>%
  summarize(
    total_villages = n(),
    mean_agriculture_percent = mean(agriculture_percent, na.rm = TRUE),
    mean_urban_percent = mean(urban_percent, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

write.csv(summary_stats, "outputs/village_lulc_summarystats.csv", row.names = FALSE)

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

getwd()

# Define state and district codes 
# # ALL
# state_codes <- c("OD", "AP", "TN")
# districts <- list(
#   OD = 1:6,
#   AP = 1:9,
#   TN = 1:12
# )
# years <- 1990:2013

# Test
state_codes <- c("OD")
districts <- list(
  OD = 1:1
)
years <- 2013:2013

# 
# # Test start
# r <- raster("data/tif/Landsat5_classified_dtype_compress/OD_6_classified_2000.tif")
# unique_vals <- unique(values(r))
# print(unique_vals)
# 
# dataType(r)
# 
# # Load one shapefile and raster
# poly <- st_read("data/shp/OD_6_AllData.shp", quiet = TRUE)
# 
# # Reproject raster to match polygon (if needed)
# r <- projectRaster(r, crs = st_crs(poly)$wkt)
# 
# # Extract values for the first polygon
# vals <- exact_extract(r, poly[1, ], fun = function(vals, cov) round(vals))
# 
# # Look at frequency of values
# table(unlist(vals))
# 
# plot(r)
# plot(poly, add = TRUE)
# 
# # Function to count pixels in classes 0, 1, 2
# extract_class_counts <- function(values, coverage) {
#   values <- values[!is.na(values)]
#   tab <- table(factor(values, levels = 0:2))
#   # Return named values
#   out <- as.integer(tab)
#   names(out) <- c("Other", "DryAqua", "Aqua")
#   as.data.frame(as.list(out))
# }
# 
# # Run on the full shapefile
# result <- exact_extract(r, poly, fun = extract_class_counts)
# 
# # Combine with UniqueID
# final <- bind_cols(UniqueID = poly$UniqueID, result)
# 
# # Calculate area in m²
# final <- final %>%
#   mutate(DryAqua_area_m2 = DryAqua * 900,
#          Aqua_area_m2 = Aqua * 900)
# 
# print(head(final))
# 
# # Test End




# Create output directory if it doesn't exist
dir.create("outputs", showWarnings = FALSE)

# Function to process a single district for a specific year
process_district <- function(state_code, dist_code, year) {
  # Construct file paths
  boundary_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  
  # aqua_path <- sprintf("data/tif/Landsat5_classified_dtype_compress/%s_%d_classified_%d.tif", # change path for Landsat8 relevant folders for 2013-25 years and run again
  #                      state_code, dist_code, year)

  aqua_path <- sprintf("data/tif/Landsat8_classified_dtype_compress/%s/%s_%d/%s_%d_classified_%d.tif",
                       state_code, state_code, dist_code, state_code, dist_code, year)

  aqua_area_path <- sprintf("outputs/%s_%d_%d_aqua_area_by_zone.csv", 
                            state_code, dist_code, year)
  
  cat(sprintf("Processing %s_%d for year %d\n", state_code, dist_code, year))
  
  if (!file.exists(boundary_path) || !file.exists(aqua_path)) {
    cat(sprintf("Missing files for %s_%d_%d\n", state_code, dist_code, year))
    return(FALSE)
  }
  
  tryCatch({
    # Load
    poly <- st_read(boundary_path, quiet = TRUE)
    r <- raster(aqua_path)
    r <- projectRaster(r, crs = st_crs(poly)$wkt)
    
    # Check UniqueID
    if (!"UniqueID" %in% colnames(poly)) stop("UniqueID column missing")
    
    # Extraction
    extract_class_counts <- function(values, coverage) {
      values <- values[!is.na(values)]
      tab <- table(factor(values, levels = 0:2))
      out <- as.integer(tab)
      names(out) <- c("Other", "DryAqua", "Aqua")
      as.data.frame(as.list(out))
    }
    
    result <- exact_extract(r, poly, fun = extract_class_counts)
    
    # Combine with IDs
    out_df <- bind_cols(UniqueID = poly$UniqueID, result) %>%
      mutate(DryAqua_area_m2 = DryAqua * 900,
             Aqua_area_m2 = Aqua * 900)
    
    # Save
    write.csv(out_df, aqua_area_path, row.names = FALSE)
    cat(sprintf("Success: %s_%d_%d\n", state_code, dist_code, year))
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("Error processing %s_%d_%d: %s\n", state_code, dist_code, year, e$message))
    return(FALSE)
  })
}


# Initialize counters
total_combinations <- 0
successful_combinations <- 0

# Create a log file
log_file <- "processing_log.txt"
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run the following after all are done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  write.csv(all_results, "outputs/all_districts_aqua_summary.csv", row.names = FALSE)
  cat("Created summary file with all results: outputs/all_districts_aqua_summary.csv\n")
  cat("Created summary file with all results: outputs/all_districts_aqua_summary.csv\n",
      file = log_file, append = TRUE)
}

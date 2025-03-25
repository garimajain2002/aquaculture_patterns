# Convert wide form aquaculture data into long form

library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

getwd()

# Define state codes and districts
state_codes <- c("OD", "AP", "TN")
districts <- list(
  OD = 1:6,
  AP = 1:9,
  TN = 1:12
)

# Function to process each shapefile and convert to long format
process_shapefile_to_long <- function(state_code, dist_code) {
  # Construct file path
  boundary_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  
  # Read shapefile
  cat(sprintf("Processing %s district %d...\n", state_code, dist_code))
  
  # Error handling for file reading
  tryCatch({
    # Read the shapefile
    villages <- sf::st_read(boundary_path, quiet = TRUE)
    
    # Print all column names for debugging
    cat("All columns in shapefile:", paste(names(villages), collapse=", "), "\n")
    
    # Extract UniqueID from shapefile
    if ("UniqueID" %in% names(villages)) {
      villages$id <- villages$UniqueID
    } else if ("UNIQUE_ID" %in% names(villages)) {
      villages$id <- villages$UNIQUE_ID
    } else if ("unique_id" %in% names(villages)) {
      villages$id <- villages$unique_id
    } else {
      # Generate IDs if none found
      villages$id <- paste0(state_code, "_", dist_code, "_", 1:nrow(villages))
      cat(sprintf("Warning: No UniqueID found in %s district %d. Generated IDs.\n", 
                  state_code, dist_code))
    }
    
    # Convert villages to data frame (dropping geometry)
    villages_df <- st_drop_geometry(villages)
    
    # Based on the pattern observed in the column names:
    # - Aqua_YYYY: These are read correctly
    # - DryAqua_YYYY: These are truncated to DryAqua_YY where YY is a number starting from 20 and then gets renamed after
    
    # Identify Aqua columns (these should be correct)
    aqua_cols <- grep("^Aqua_[0-9]{4}$", names(villages_df), value = TRUE)
    
    # Extract years from Aqua columns
    aqua_years <- as.integer(sub("^Aqua_([0-9]{4})$", "\\1", aqua_cols))
    
    # For DryAqua columns, we need to map the truncated names to actual years
    dry_aqua_cols <- grep("^DryAqua_[0-9]{2}$", names(villages_df), value = TRUE)
    
    # Create a map from truncated DryAqua columns to years
    dry_aqua_year_map <- list()
    
    # If we have the same number of DryAqua and Aqua columns, we can map them directly
    if (length(dry_aqua_cols) == length(aqua_cols)) {
      # Sort both lists to ensure correct pairing
      dry_aqua_cols <- sort(dry_aqua_cols)
      aqua_cols <- sort(aqua_cols)
      aqua_years <- sort(aqua_years)
      
      for (i in 1:length(aqua_years)) {
        year <- aqua_years[i]
        dry_aqua_year_map[[as.character(year)]] <- dry_aqua_cols[i]
      }
    } else {
      # Try to match based on the expected pattern
      for (col in dry_aqua_cols) {
        # Extract the number from the column name
        num <- as.integer(sub("^DryAqua_([0-9]{2})$", "\\1", col))
        
        # Map to years starting from 2013 (if num is 20) through 2025 (if num is 32)
        if (num >= 20 && num <= 32) {
          year <- 2013 + (num - 20)
          dry_aqua_year_map[[as.character(year)]] <- col
        }
      }
    }
    
    # Print the mapping for verification
    cat("DryAqua column mapping:\n")
    for (year in names(dry_aqua_year_map)) {
      cat(sprintf("  Year %s -> Column %s\n", year, dry_aqua_year_map[[year]]))
    }
    
    # Create a list to store year datasets
    year_dfs <- list()
    
    # Process each year
    for (year in 2013:2025) {
      year_str <- as.character(year)
      aqua_col <- paste0("Aqua_", year)
      
      # Check if we have data for this year
      has_data <- FALSE
      
      # Create subset with id, year, and values
      year_df <- data.frame(
        UniqueID = villages_df$id,
        state_code = state_code,
        district_code = dist_code,
        year = year
      )
      
      # Add DryAqua data if we have a mapping for this year
      if (year_str %in% names(dry_aqua_year_map)) {
        dry_col <- dry_aqua_year_map[[year_str]]
        year_df$DryAqua <- villages_df[[dry_col]]
        has_data <- TRUE
        cat(sprintf("Added DryAqua data for year %d from column %s\n", year, dry_col))
      } else {
        year_df$DryAqua <- NA
      }
      
      # Add Aqua data if column exists
      if (aqua_col %in% names(villages_df)) {
        year_df$Aqua <- villages_df[[aqua_col]]
        has_data <- TRUE
        cat(sprintf("Added Aqua data for year %d from column %s\n", year, aqua_col))
      } else {
        year_df$Aqua <- NA
      }
      
      # Only add to the list if we have data for this year
      if (has_data) {
        year_dfs[[year_str]] <- year_df
      }
    }
    
    # Combine all years into one data frame
    if (length(year_dfs) > 0) {
      combined_df <- do.call(rbind, year_dfs)
      return(combined_df)
    } else {
      cat(sprintf("No valid year data found for %s district %d.\n", state_code, dist_code))
      return(NULL)
    }
    
  }, error = function(e) {
    cat(sprintf("Error processing %s district %d: %s\n", state_code, dist_code, e$message))
    return(NULL)
  })
}

# Process all districts for all states
results_list <- list()

for (state_code in state_codes) {
  for (dist_code in districts[[state_code]]) {
    result <- process_shapefile_to_long(state_code, dist_code)
    if (!is.null(result)) {
      results_list[[length(results_list) + 1]] <- result
    }
  }
}

# Combine all results
all_results <- do.call(rbind, results_list)

# Check and clean data
all_results <- all_results %>%
  # Set appropriate types
  mutate(
    UniqueID = as.character(UniqueID),
    state_code = as.character(state_code),
    district_code = as.integer(district_code),
    year = as.integer(year),
    DryAqua = as.numeric(DryAqua),
    Aqua = as.numeric(Aqua)
  ) %>%
  # Sort by UniqueID and year
  arrange(UniqueID, year)

# Save to CSV
write.csv(all_results, "outputs/aquaculture_by_year_long_format.csv", row.names = FALSE)

# Print summary
cat("\nData conversion complete. Total records:", nrow(all_results), "\n")
cat("Years included:", paste(sort(unique(all_results$year)), collapse=", "), "\n")
cat("Total unique villages:", length(unique(all_results$UniqueID)), "\n")

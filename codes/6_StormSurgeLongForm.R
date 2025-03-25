# Combine all surge data in one csv 

# Load required libraries
library(sf)
library(dplyr)
library(tidyr)
library(readr)

getwd() 

# Define state and district codes
state_codes <- c("OD", "AP", "TN")
district_codes <- list(
  OD = 1:6,
  AP = 1:9,
  TN = 1:12
)

# Define surge mapping
surge_mapping <- list(
  "Surge_A12_1964_Clipped" = c("TN_5", "TN_6", "TN_7", "TN_8", "TN_9"),
  "Surge_A26_1971_Clipped" = c("OD_1", "OD_2", "OD_3", "OD_4"),
  "Surge_A29_1974_Clipped" = c("OD_1", "OD_2", "OD_3"),
  "Surge_A33_1976_Clipped" = c("AP_1", "AP_2", "AP_3", "AP_4"),
  "Surge_A34_1977_Clipped" = c("AP_4", "AP_5", "AP_6", "AP_7", "AP_8"),
  "Surge_A37_1981_Clipped" = c("OD_1", "OD_2", "OD_3", "OD_4", "OD_5"),
  "Surge_A38_1982_Clipped" = c("OD_1", "OD_2", "OD_3", "OD_4"),
  "Surge_A39_1984_Clipped" = c("OD_1", "OD_2", "OD_3", "OD_4", "OD_5"),
  "Surge_A41_1986_Clipped" = c("AP_1", "AP_2", "AP_3", "AP_4"),
  "Surge_A43_1989_Clipped" = c("AP_6", "AP_7", "AP_8", "AP_9"),
  "Surge_A44_1990_Clipped" = c("AP_5", "AP_6", "AP_7", "AP_8"),
  "Surge_A46_1992_Clipped" = c("TN_9", "TN_10", "TN_11"),
  "Surge_A47_1993_Clipped" = c("TN_3", "TN_4", "TN_5", "TN_6", "TN_7"),
  "Surge_A48_1996_Clipped" = c("AP_3", "AP_4", "AP_5"),
  "Surge_A50_1999_Clipped" = c("OD_1", "OD_2", "OD_3", "OD_4", "OD_5"),
  "Surge_A55_2011_Clipped" = c("TN_2", "TN_3", "TN_4", "TN_5")
)

# Function to read shapefile
read_shapefile <- function(state_code, dist_code) {
  file_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  if (file.exists(file_path)) {
    print(paste("Reading", file_path))
    return(st_read(file_path))
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
}

# Initialize empty list to store all data frames
all_data_frames <- list()

# Read all shapefiles
for (state_code in state_codes) {
  for (dist_code in district_codes[[state_code]]) {
    shp_data <- read_shapefile(state_code, dist_code)
    if (!is.null(shp_data)) {
      # Add a source column to identify the origin (needed for processing but will be removed later)
      shp_data$source <- paste0(state_code, "_", dist_code)
      
      # Remove geometry column for merging
      shp_data_no_geom <- st_drop_geometry(shp_data)
      
      # Keep only UniqueID, source, and surge columns
      surge_cols <- grep("Surge_", colnames(shp_data_no_geom), value = TRUE)
      shp_data_filtered <- shp_data_no_geom %>%
        select(UniqueID, source, all_of(surge_cols))
      
      # Add to list
      all_data_frames[[paste0(state_code, "_", dist_code)]] <- shp_data_filtered
    }
  }
}

# Check if we have any data
if (length(all_data_frames) == 0) {
  stop("No shapefiles were loaded successfully. Check your file paths.")
}

# Create a set of all unique surge columns across all data frames
all_surge_cols <- unique(unlist(lapply(all_data_frames, function(df) {
  grep("Surge_", colnames(df), value = TRUE)
})))

# Initialize merged data with the first data frame
merged_data <- all_data_frames[[1]]

# Add missing surge columns to first data frame
for (col in all_surge_cols) {
  if (!col %in% colnames(merged_data)) {
    merged_data[[col]] <- NA
  }
}

# Merge all remaining data frames
for (i in 2:length(all_data_frames)) {
  df <- all_data_frames[[i]]
  
  # Add missing surge columns to this data frame
  for (col in all_surge_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA
    }
  }
  
  # Perform full join to merge with existing data
  merged_data <- full_join(merged_data, df, by = "UniqueID", suffix = c("", paste0(".", i)))
  
  # Handle the duplicate columns from the join
  # Keep the source column from the first dataset
  if (paste0("source.", i) %in% colnames(merged_data)) {
    merged_data$source <- ifelse(is.na(merged_data$source), 
                                 merged_data[[paste0("source.", i)]], 
                                 merged_data$source)
    merged_data <- select(merged_data, -paste0("source.", i))
  }
  
  # For surge columns, prefer non-NA values
  for (col in all_surge_cols) {
    col_i <- paste0(col, ".", i)
    if (col_i %in% colnames(merged_data)) {
      merged_data[[col]] <- ifelse(is.na(merged_data[[col]]), 
                                   merged_data[[col_i]], 
                                   merged_data[[col]])
      merged_data <- select(merged_data, -col_i)
    }
  }
}

# Fill missing values based on the surge mapping
for (surge_col in names(surge_mapping)) {
  col_name <- surge_col
  if (!col_name %in% colnames(merged_data)) {
    # Skip if the column doesn't exist in the merged data
    next
  }
  
  applicable_districts <- surge_mapping[[surge_col]]
  
  # For districts not in the mapping, set the value to 0 (not affected)
  merged_data <- merged_data %>%
    mutate(!!col_name := ifelse(is.na(!!sym(col_name)) & !source %in% applicable_districts, 0, !!sym(col_name)))
}

# Check for any remaining NA values and convert to 0
merged_data[is.na(merged_data)] <- 0

# Clean up the data frame - remove source column and keep only UniqueID and surge columns
merged_data <- merged_data %>%
  select(UniqueID, all_of(all_surge_cols)) %>%
  distinct(UniqueID, .keep_all = TRUE)

# Write to CSV
write_csv(merged_data, "outputs/master_surge_dataset.csv")

# Print summary
cat("Created master dataset with", nrow(merged_data), "rows and", ncol(merged_data), "columns\n")
cat("Surge columns:", paste(all_surge_cols, collapse = ", "), "\n")

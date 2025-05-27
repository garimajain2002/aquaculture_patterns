# Bind all villages census data (2001) 

# Load required libraries
library(sf)
library(dplyr)
library(readr)

getwd() 

# Define state and district codes
state_codes <- c("OD", "AP", "TN")
district_codes <- list(
  OD = 1:6,
  AP = 1:9,
  TN = 1:12
)

# Define the fields to keep
fields_to_keep <- c(
  "UniqueID", "SID", "DID", "TID", "C_CODE01", "Name", "TRU", 
  "No_HH", "TOT_P", "TOT_M", "TOT_F", "P_LIT", "M_LIT", "F_LIT", 
  "P_ILL", "M_ILL", "F_ILL", "TOT_WORK_P", "TOT_WORK_M", "TOT_WORK_F", 
  "MAINWORK_P", "MAINWORK_M", "MAINWORK_F", "MARGWORK_P", "MARGWORK_M", 
  "MARGWORK_F", "SC_P", "SC_M", "SC_F", "ST_P", "ST_M", "ST_F", 
  "TOT_IRR", "UN_IRR", "CULT_WASTE", "DEM_avg", "VHCF_NDMA", "Shape_Area"
)

# Function to read shapefile and keep only specified fields
read_shapefile <- function(state_code, dist_code) {
  file_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  if (file.exists(file_path)) {
    print(paste("Reading", file_path))
    
    # Read the shapefile
    shp_data <- st_read(file_path)
    
    # Check if all required fields exist
    missing_fields <- setdiff(fields_to_keep, colnames(shp_data))
    if (length(missing_fields) > 0) {
      warning(paste("File", file_path, "is missing these fields:", 
                    paste(missing_fields, collapse = ", ")))
    }
    
    # Keep only the fields that exist in this shapefile
    existing_fields <- intersect(fields_to_keep, colnames(shp_data))
    shp_data_filtered <- shp_data %>% select(all_of(existing_fields))
    
    # Add state and district info
    shp_data_filtered$state_code <- state_code
    shp_data_filtered$dist_code <- dist_code
    
    return(shp_data_filtered)
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
}

# Initialize empty list to store all data frames
all_village_boundaries <- list()

# Read all shapefiles
for (state_code in state_codes) {
  for (dist_code in district_codes[[state_code]]) {
    village_data <- read_shapefile(state_code, dist_code)
    if (!is.null(village_data)) {
      # Add to list
      all_village_boundaries[[paste0(state_code, "_", dist_code)]] <- village_data
    }
  }
}

# Check if we have any data
if (length(all_village_boundaries) == 0) {
  stop("No shapefiles were loaded successfully. Check file paths.")
}

# Bind all data frames together
# We use bind_rows to automatically handle columns that might be missing in some files
merged_villages <- bind_rows(all_village_boundaries)

# Check for duplicate UniqueIDs
duplicate_ids <- merged_villages %>%
  st_drop_geometry() %>%
  count(UniqueID) %>%
  filter(n > 1)

if (nrow(duplicate_ids) > 0) {
  warning(paste(nrow(duplicate_ids), "duplicate UniqueIDs found"))
  # Keep only the first occurrence of each UniqueID
  merged_villages <- merged_villages %>%
    group_by(UniqueID) %>%
    slice(1) %>%
    ungroup()
}

# Add a count of total fields to ensure we have all the expected fields
expected_field_count <- length(fields_to_keep)
actual_field_count <- sum(fields_to_keep %in% colnames(merged_villages))

cat("Expected", expected_field_count, "fields, found", actual_field_count, "fields\n")
if (actual_field_count < expected_field_count) {
  missing_fields <- setdiff(fields_to_keep, colnames(merged_villages))
  warning(paste("Missing fields in final dataset:", paste(missing_fields, collapse = ", ")))
}

# Reorder columns to match the requested order
available_fields <- intersect(fields_to_keep, colnames(merged_villages))
merged_villages <- merged_villages %>%
  select(all_of(available_fields), everything())

# Save as CSV without geometry for easier data analysis
merged_villages_df <- st_drop_geometry(merged_villages)
write_csv(merged_villages_df, "outputs/census_2001.csv")

# Print summary
cat("Created merged village boundaries dataset with", nrow(merged_villages), "villages\n")
cat("Dataset contains", ncol(merged_villages) - 1, "attributes \n")


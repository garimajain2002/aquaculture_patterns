# Combine all surge data in one csv 

# Load required libraries
library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

getwd() 

# Define state and district codes
state_codes <- c("OD", "AP", "TN")
district_codes <- list(
  OD = 1:6,
  AP = 1:9,
  TN = 1:12
)

# !!! Does not include tsunami affected villages yet 


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
  stop("No shapefiles were loaded successfully. Check file paths.")
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

surge <- merged_data 

# First rename the Surge fields to again include years (.shp truncates after 10 characters)
surge$Surge_A12_1964	<- surge$Surge_A12_
surge$Surge_A26_1971  <- surge$Surge_A26_	
surge$Surge_A29_1974	<- surge$Surge_A29_
surge$Surge_A33_1976	<- surge$Surge_A33_
surge$Surge_A34_1977	<- surge$Surge_A34_
surge$Surge_A37_1981	<- surge$Surge_A37_
surge$Surge_A38_1982	<- surge$Surge_A38_
surge$Surge_A39_1984	<- surge$Surge_A39_
surge$Surge_A41_1986	<- surge$Surge_A41_
surge$Surge_A43_1989	<- surge$Surge_A43_
surge$Surge_A44_1990	<- surge$Surge_A44_
surge$Surge_A46_1992  <- surge$Surge_A46_
surge$Surge_A47_1993	<- surge$Surge_A47_
surge$Surge_A48_1996	<- surge$Surge_A48_
surge$Surge_A50_1999	<- surge$Surge_A50_
surge$Surge_A55_2011	<- surge$Surge_A55_	

head(surge)

surge <- surge %>%
  select( UniqueID, Surge_A12_1964,Surge_A26_1971, Surge_A29_1974, Surge_A33_1976, Surge_A34_1977, Surge_A37_1981, Surge_A38_1982, Surge_A39_1984, 
         Surge_A41_1986, Surge_A43_1989, Surge_A44_1990, Surge_A46_1992, Surge_A47_1993, Surge_A48_1996, Surge_A50_1999, Surge_A55_2011 )

head(surge)

#----------------------
# Add Tsunami data
#----------------------

#----------------------
# Long form Conversion
#----------------------

villageYears <- read.csv("data/All_villageYears.csv")

# Extract the surge years
surge_long <- surge %>%
  # Select only UniqueID and surge columns
  select(UniqueID, starts_with("Surge_")) %>%
  # Pivot to long format
  pivot_longer(cols = starts_with("Surge_"), 
               names_to = "surge_column", 
               values_to = "affected") %>%
  # Keep only rows where affected = 1
  filter(affected == 1) %>%
  # Extract year from column name
  mutate(surge_year = as.numeric(gsub(".*_([0-9]{4})$", "\\1", surge_column))) %>%
  # For each UniqueID, find the earliest surge year
  group_by(UniqueID) %>%
  summarize(first_surge_year = min(surge_year))

# Join with villageYears and create postSurge
villageYears <- villageYears %>%
  left_join(surge_long, by = "UniqueID") %>%
  mutate(
    # If there's no surge year, or the current year is before the first surge,
    # postSurge is 0, otherwise 1
    postSurge = case_when(
      is.na(first_surge_year) ~ 0,
      year < first_surge_year ~ 0,
      TRUE ~ 1
    )
  )

head(villageYears)

table(villageYears$first_surge_year)

#----------------------
# Graph the trend
#----------------------

# Cumulative
# Filter rows with a known first surge year
surge_trend <- villageYears %>%
  filter(!is.na(first_surge_year)) %>%
  distinct(UniqueID, first_surge_year) %>%  # Ensure one entry per village
  group_by(first_surge_year) %>%
  summarize(new_affected_villages = n()) %>%
  arrange(first_surge_year) %>%
  mutate(cumulative_villages = cumsum(new_affected_villages))

# Plot the cumulative number of affected villages over time
ggplot(surge_trend, aes(x = first_surge_year, y = cumulative_villages)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "darkred") +
  scale_x_continuous(breaks = seq(1960, 2025, 5)) +
  labs(
    title = "Cumulative Number of Villages Affected by Storm Surges",
    x = "Year",
    y = "Cumulative Number of Villages"
  ) +
  theme_minimal()

ggsave("outputs/StormAffectedVillagesTrend_1964-2011.png", width = 14, height = 6, dpi = 300)


# By state
# Prepare cumulative surge trend by state
surge_trend_by_state <- villageYears %>%
  filter(!is.na(first_surge_year)) %>%
  distinct(UniqueID, state_code, first_surge_year) %>%
  group_by(state_code, first_surge_year) %>%
  summarize(new_affected_villages = n(), .groups = "drop") %>%
  arrange(state_code, first_surge_year) %>%
  group_by(state_code) %>%
  mutate(cumulative_villages = cumsum(new_affected_villages))

# Plot cumulative surge trend by state
ggplot(surge_trend_by_state, aes(x = first_surge_year, y = cumulative_villages, color = state_code)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1960, 2025, 5)) +
  labs(
    title = "Cumulative Number of Villages Affected by Storm Surges by State",
    x = "Year",
    y = "Cumulative Number of Villages",
    color = "State"
  ) +
  theme_minimal()
ggsave("outputs/StormAffectedVillagesTrend_ByState_1964-2011.png", width = 14, height = 6, dpi = 300)


#----------------------
# Export file
#----------------------

# Write to CSV
write_csv(villageYears, "outputs/master_surge_dataset.csv")

# Print summary
cat("Created master dataset with", nrow(villageYears), "rows and", ncol(villageYears), "columns\n")


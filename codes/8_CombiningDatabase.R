# Combine Aquaculture (Y), Salinity (X1) and Storm Data to create master database 

# First combine Aqua for multi years with salinity with multi-years. 
# Second, for surge years (all prior to 2011 i.e. pre-landsat 8 period), if a village is ever affected by any of the 16 surges, then mark storm yes. 
# Discuss with Dylan how to deal with LULC and Census (match to all or specific years only). 

# Convert the variables to usable values (absolute to %s etc.) 


library(dplyr)

getwd() 

#_________________________________________________________________________
# Read relevant dataframes
# ________________________________________________________________________

# Read relevant data (2013-2025 first, later repeat for 1990-2012 data)
aqua <- read.csv("outputs/all_districts_aquaculture_summary.csv")
salinity <- read.csv ("outputs/all_districts_salinity_summary.csv")
surge <- read.csv("outputs/master_surge_dataset.csv")
lulc <- read.csv("outputs/village_lulc_percentages.csv")
census <- read.csv("outputs/merged_village_boundaries.csv")
missing <- read.csv("logs/missingSatelliteData.csv") # villages and years either missing or affected by clouds in satellite data 
distance <- read.csv("outputs/all_SeaDistances.csv")


# Assess data before merging 
head(aqua)
print("Aqua data dimensions:")
print(dim(aqua))
print("Aqua data column names:")
print(colnames(aqua))

head(salinity)
print("Salinity data dimensions:")
print(dim(salinity))
print("Salinity data column names:")
print(colnames(salinity))

head(surge)
print("Surge data dimensions:")
print(dim(surge))
print("Surge data column names:")
print(colnames(surge))

head(lulc)
print("LULC data dimensions:")
print(dim(lulc))
print("LULC data column names:")
print(colnames(lulc))

head(census)
print("Census data dimensions:")
print(dim(census))
print("Census data column names:")
print(colnames(census))


#_________________________________________________________________________
# Remove missing/cloudy satellite areas - villages by years
# ________________________________________________________________________

# Make missing ID-year combinations NA in both salinity and aquaculture
# AQUACULTURE
# Function to replace erroneous aquaculture data with NA
update_aqua_data <- function(aqua, missing) {
  # Make a copy of the original dataframe to avoid modifying the input directly
  result_df <- aqua
  
  # Create a matching key for both dataframes
  aqua_keys <- paste(aqua$UniqueID, aqua$Year, sep = "_")
  missing_keys <- paste(missing$UniqueID, missing$Year, sep = "_")
  
  # Find rows in aqua that match with missing
  matches <- which(aqua_keys %in% missing_keys)
  
  # List of columns to set to NA
  aqua_cols <- c("aqua_pixel_count", "aqua_area_m2", 
                 "dryaqua_pixel_count", "dryaqua_area_m2")
  
  # Set the specified columns to NA for the matching rows
  if (length(matches) > 0) {
    for (col in aqua_cols) {
      if (col %in% names(result_df)) {
        result_df[matches, col] <- NA
      }
    }
    
    # Print number of rows changed
    cat("Modified", length(matches), "rows in the aqua dataframe\n")
  } else {
    cat("No matching rows found between the dataframes\n")
  }
  
  return(result_df)
}

# Direct usage with your existing dataframes
updated_aqua <- update_aqua_data(aqua, missing)

# Check how many values were changed for specific columns
for (col in c("aqua_pixel_count", "aqua_area_m2", "dryaqua_pixel_count", "dryaqua_area_m2")) {
  if (col %in% names(aqua)) {
    changed <- sum(is.na(updated_aqua[[col]]) & !is.na(aqua[[col]]))
    cat("Changed", changed, "values to NA in column", col, "\n")
  }
}

aqua <- updated_aqua

# SALINITY 
# Function to replace erroneous aquaculture data with NA
update_salinity_data <- function(salinity, missing) {
  # Make a copy of the original dataframe to avoid modifying the input directly
  result_df_salinity <- salinity
  
  # Create a matching key for both dataframes
  salinity_keys <- paste(salinity$UniqueID, salinity$Year, sep = "_")
  missing_keys <- paste(missing$UniqueID, missing$Year, sep = "_")
  
  # Find rows in salinity that match with missing
  matches <- which(salinity_keys %in% missing_keys)
  
  # List of columns to set to NA
  salinity_cols <- c("saline_pixel_count", "saline_area_m2")
  
  # Set the specified columns to NA for the matching rows
  if (length(matches) > 0) {
    for (col in salinity_cols) {
      if (col %in% names(result_df_salinity)) {
        result_df_salinity[matches, col] <- NA
      }
    }
    
    # Print number of rows changed
    cat("Modified", length(matches), "rows in the aqua dataframe\n")
  } else {
    cat("No matching rows found between the dataframes\n")
  }
  
  return(result_df_salinity)
}

# Direct usage with your existing dataframes
updated_salinity <- update_salinity_data(salinity, missing)

# Check how many values were changed for specific columns
for (col in c("salinity_pixel_count", "salinity_area_m2")) {
  if (col %in% names(salinity)) {
    changed <- sum(is.na(updated_salinity[[col]]) & !is.na(salinity[[col]]))
    cat("Changed", changed, "values to NA in column", col, "\n")
  }
}

salinity <- updated_salinity


#_________________________________________________________________________
# Merging Dataframes
# ________________________________________________________________________

# Check if the column names match 
# Check the structure of the join keys in both datasets
str(aqua[, c("UniqueID", "Row_no", "Year")])
str(salinity[, c("UniqueID", "Row_no", "Year")])

# Check for case differences in character columns
if(is.character(aqua$UniqueID)) {
  print("Sample aqua UniqueIDs:")
  print(head(aqua$UniqueID))
}

if(is.character(salinity$UniqueID)) {
  print("Sample salinity UniqueIDs:")
  print(head(salinity$UniqueID))
}

# Check for value ranges and data types
summary(aqua[, c("UniqueID", "Row_no", "Year")])
summary(salinity[, c("UniqueID", "Row_no", "Year")])


# Check for duplicates in aqua
aqua_duplicates <- aqua %>%
  group_by(UniqueID, Row_no, Year) %>%
  filter(n() > 1) %>%
  ungroup()

# Check for duplicates in salinity
salinity_duplicates <- salinity %>%
  group_by(UniqueID, Row_no, Year) %>%
  filter(n() > 1) %>%
  ungroup()

# View the results
print(paste("Duplicate rows in aqua:", nrow(aqua_duplicates)))
print(paste("Duplicate rows in salinity:", nrow(salinity_duplicates)))
# Confirmed - no duplicates 

# Perform the merge
aqua_salinity_merged <- aqua %>%
  left_join(salinity, by = c("UniqueID", "Row_no", "Year"), suffix = c("_aqua", "_salinity"))

print("Merged data dimensions:")
print(dim(aqua_salinity_merged))
print("Merged data column names:")
print(colnames(aqua_salinity_merged))

head(aqua_salinity_merged)


# Combine Surge dataset 
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

# For Landsat8 time period (post 2013), create a combined surge variable and mark all year instances as 1 if affected by any of these surges 
# More checks will be needed for Landsat 5 years (1990-2012) 

surge$postSurge <- ifelse(
  surge$Surge_A12_1964 == 1 | surge$Surge_A26_1971 == 1 | 
    surge$Surge_A29_1974 == 1 | surge$Surge_A33_1976 == 1 | 
    surge$Surge_A34_1977 == 1 | surge$Surge_A37_1981 == 1 | 
    surge$Surge_A38_1982 == 1 | surge$Surge_A39_1984 == 1 | 
    surge$Surge_A41_1986 == 1 | surge$Surge_A43_1989 == 1 | 
    surge$Surge_A44_1990 == 1 | surge$Surge_A46_1992 == 1 | 
    surge$Surge_A47_1993 == 1 | surge$Surge_A48_1996 == 1 | 
    surge$Surge_A50_1999 == 1 | surge$Surge_A55_2011 == 1, 
  1, 0
)

# Add surge data (postSurge only for now)
aqua_salinity_surge <- merge(aqua_salinity_merged, surge[, c("UniqueID", "postSurge")], by = "UniqueID", all.x = TRUE)
print("Merged data dimensions:")
print(dim(aqua_salinity_surge))

# Add total area from LULC calculations 
aqua_salinity_surge <- merge(aqua_salinity_surge, lulc[ , c("UniqueID", "total_area_ha")], by = "UniqueID", all.x = TRUE)
print("Merged data dimensions:")
print(dim(aqua_salinity_surge))

# Add distance from the sea (NEAR_DIST)
aqua_salinity_surge <- merge(aqua_salinity_surge, distance[ , c("UniqueID", "NEAR_DIST")], by = "UniqueID", all.x = TRUE)
print("Merged data dimensions:")
print(dim(aqua_salinity_surge))
# Create a categorical variable for distance 
# Define breaks based on the range of NEAR_DIST
breaks <- seq(min(aqua_salinity_surge$NEAR_DIST), max(aqua_salinity_surge$NEAR_DIST), length.out = 6)
# Create the categorical variable
aqua_salinity_surge$Sea_Dist <- cut(aqua_salinity_surge$NEAR_DIST, breaks = breaks, labels = c("Very Near", "Near", "Medium", "Far", "Very Far"), include.lowest = TRUE)
# Check distribution
table(aqua_salinity_surge$Sea_Dist)

head(aqua_salinity_surge)


# Add DEM_avg from census  
aqua_salinity_surge <- merge(aqua_salinity_surge, census[ , c("UniqueID", "DEM_avg")], by = "UniqueID", all.x = TRUE)
print("Merged data dimensions:")
print(dim(aqua_salinity_surge))

#_________________________________________________________________________
# Calculate relevant metrics for aqua and salinity
# ________________________________________________________________________

# Convert all areas to same unit to make them comparable (and check against total area) 
aqua_salinity_surge$DryAqua_ha <- aqua_salinity_surge$dryaqua_area_m2 / 10000
aqua_salinity_surge$Aqua_ha <- aqua_salinity_surge$aqua_area_m2 / 10000
aqua_salinity_surge$Saline_ha <- aqua_salinity_surge$saline_pixel_count*900 / 10000

aqua_salinity_surge$Aqua_perc <- aqua_salinity_surge$Aqua_ha / aqua_salinity_surge$total_area_ha * 100
aqua_salinity_surge$DryAqua_perc <- aqua_salinity_surge$DryAqua_ha / aqua_salinity_surge$total_area_ha * 100
aqua_salinity_surge$Saline_perc <- aqua_salinity_surge$Saline_ha / aqua_salinity_surge$total_area_ha * 100

head(aqua_salinity_surge)

summary(aqua_salinity_surge$Saline_perc)
summary(aqua_salinity_surge$Aqua_perc)


# Calculate % change in aquaculture, dry aquaculture and salinity from previous year 
aqua_salinity_surge <- aqua_salinity_surge %>%
  arrange(Row_no) %>%  # Ensure data is sorted by Row_no which is already sorted for UniqueID and Year
  group_by(UniqueID) %>%
  mutate(
    pct_change_aqua = (Aqua_ha - lag(Aqua_ha)) / lag(Aqua_ha) * 100,
    pct_change_dryaqua = (DryAqua_ha - lag(DryAqua_ha)) / lag(DryAqua_ha) * 100,
    pct_change_salinity = (Saline_ha - lag(Saline_ha)) / lag(Saline_ha) * 100
  ) %>%
  ungroup()

# Create Saline binary 
median_saline <- median(aqua_salinity_surge$Saline_perc, na.rm = TRUE)
aqua_salinity_surge$Saline <- ifelse(aqua_salinity_surge$Saline_perc >= median_saline, 1, 0)

# Categorize villages into four groups
aqua_salinity_surge <- aqua_salinity_surge %>%
  mutate(
    Saline_Storm_Category = factor(case_when(
      Saline == 1 & postSurge == 1 ~ "Flooded & High Saline",
      Saline == 0 & postSurge == 1 ~ "Flooded & Low Saline",
      Saline == 1 & postSurge == 0 ~ "Non-Flooded & High Saline",
      Saline == 0 & postSurge == 0 ~ "Non-Flooded & Low Saline"
    ), levels = c(
      "Flooded & High Saline",
      "Non-Flooded & High Saline",
      "Flooded & Low Saline",
      "Non-Flooded & Low Saline"
    ))
  )




head(aqua_salinity_surge)

# Create another variable - land unsuitable for agriculture as a sum of saline and aquaculture 
# This is because salne % will automatically reduce as aquaculture increases in a village (those pixels get masked) 
aqua_salinity_surge$NotforRice_ha <- aqua_salinity_surge$Aqua_ha + aqua_salinity_surge$Saline_ha
aqua_salinity_surge$NotforRice_perc <- aqua_salinity_surge$Aqua_perc + aqua_salinity_surge$Saline_perc


#_________________________________________________________________________
# Clean database and save
# ________________________________________________________________________

aqua_salinity_surge$State <- aqua_salinity_surge$State_salinity
aqua_salinity_surge$District <- aqua_salinity_surge$District_salinity

aqua_salinity_surge <- aqua_salinity_surge %>% select(Row_no, UniqueID, Year, State, District, total_area_ha, NEAR_DIST, Sea_Dist, DEM_avg, Aqua_ha, DryAqua_ha, Saline_ha, Aqua_perc, DryAqua_perc, Saline_perc, pct_change_aqua, pct_change_dryaqua, pct_change_salinity, postSurge, Saline, Saline_Storm_Category, aqua_area_m2, dryaqua_area_m2, saline_area_m2, aqua_pixel_count, dryaqua_pixel_count, saline_pixel_count, UniqueID_Row_Year)
head(aqua_salinity_surge)


write.csv(aqua_salinity_surge, "outputs/aqua_salinity_surge_2013-2025.csv")


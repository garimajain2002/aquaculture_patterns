# Step 1 - create long data form with all uniqueIDs and years
# Step 2 - Combine salinity (X1) data with this (already in long form) 
# Step 3 - Combine surge (X2) data with this (already converted to long form in StormSurgeLongForm)
# Step 4 - Combine aquaculture and dry aquaculture areas (Y) 
# Step 5 - Merge LULC, DEM, Sea distance and Census data to this (Controls) 

# Convert the variables to usable values (absolute to %s etc.) 


library(dplyr)
library(tidyr)
getwd() 


#_________________________________________________________________________________________
# Step 1 - Make a long form of all villages with all years as the base to bind other data
#_________________________________________________________________________________________

unique_ids <- read.csv("data/All_uniqueIDs.csv")

unique_ids <- unique_ids %>%
  select(UniqueID, state_code, district_code) %>%
  distinct()

head(unique_ids)

summary(unique_ids)


# Create a dataframe for the years 1990-2025
years <- 1990:2025


# Generate all combinations of UniqueID and historical years
villageYears <- unique_ids %>%
  crossing(year = years) %>%
  # Keep only state_code and district_code columns from unique_ids
  select(UniqueID, year, state_code, district_code)

head(villageYears)
summary(villageYears)

write.csv(villageYears, "data/All_villageYears.csv")


#_________________________________________________________________________
# Read relevant dataframes
# ________________________________________________________________________

# Read relevant data
aqua <- read.csv("outputs/all_districts_aqua_summary.csv")
salinity <- read.csv ("outputs/all_districts_salinity_summary.csv")
surge <- read.csv("outputs/master_surge_dataset.csv")
lulc <- read.csv("outputs/village_lulc_percentages.csv")
census <- read.csv("outputs/census_2001.csv")
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
# convert year to Year and X to Row_no. to remain comparable with other datasets 
surge$Year <- surge$year
surge$Row_no <- surge$X

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
str(aqua[, c("UniqueID", "Year")])
str(salinity[, c("UniqueID", "Year")])

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
summary(aqua[, c("UniqueID", "Year")])
summary(salinity[, c("UniqueID", "Year")])


# Check for duplicates in aqua
aqua_duplicates <- aqua %>%
  group_by(UniqueID, Year) %>%
  filter(n() > 1) %>%
  ungroup()

# Check for duplicates in salinity
salinity_duplicates <- salinity %>%
  group_by(UniqueID, Year) %>%
  filter(n() > 1) %>%
  ungroup()

# View the results
print(paste("Duplicate rows in aqua:", nrow(aqua_duplicates)))
print(paste("Duplicate rows in salinity:", nrow(salinity_duplicates)))

# Assess and Address Duplicates 
# Take the average of the aqua and dry aqua for each of the same IDs 
aqua_clean <- aqua %>%
  group_by(UniqueID, Year) %>%
  summarise(
    Aqua_area_m2 = mean(Aqua_area_m2, na.rm = TRUE),
    DryAqua_area_m2 = mean(DryAqua_area_m2, na.rm = TRUE),
    .groups = "drop"
  )

salinity_clean <- salinity %>%
  group_by(UniqueID, Year) %>%
  summarise(
    saline_pixel_count = mean(saline_pixel_count, na.rm = TRUE),
    saline_areas_m2 = mean(saline_area_m2, na.rm = TRUE),
    State = first(State),
    District = first(District),
    .groups = "drop"
  )

# Confirm - no duplicates 
# Check for duplicates in aqua
aqua_duplicates <- aqua_clean %>%
  group_by(UniqueID, Year) %>%
  filter(n() > 1) %>%
  ungroup()

# Check for duplicates in salinity
salinity_duplicates <- salinity_clean %>%
  group_by(UniqueID, Year) %>%
  filter(n() > 1) %>%
  ungroup()



# Perform the merge
# Start with surge since that has all IDs and years 
aqua_surge <- surge %>% 
  left_join(aqua_clean, by = c("UniqueID", "Year"), suffix = c("_surge", "_aqua"))

print("Merged data dimensions:")
print(dim(aqua_surge))
print("Merged data column names:")
print(colnames(aqua_surge))

# merge salinity
aqua_salinity_surge <- aqua_surge %>%
  left_join(salinity_clean, by = c("UniqueID", "Year"), suffix = c("_aqua", "_salinity"))

print("Merged data dimensions:")
print(dim(aqua_salinity_surge))
print("Merged data column names:")
print(colnames(aqua_salinity_surge))

head(aqua_salinity_surge)


# Add total area, agri are and urban area from LULC calculations 
aqua_salinity_surge <- merge(aqua_salinity_surge, lulc[ , c("UniqueID", "total_area_ha", "agriculture_area_ha", "urban_area_ha", "agriculture_percent", "urban_percent")], by = "UniqueID", all.x = TRUE)
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


# Add census data including DEM_avg 
# Note: We are taking only 2001 Census and applying to all years. 
# To only be used as a control and not a variable of interest. 
aqua_salinity_surge <- merge(aqua_salinity_surge, census[ , c("UniqueID", "Shape_Area", "TRU", "No_HH", "TOT_P", "TOT_M", "TOT_F", "P_LIT", "M_LIT", "F_LIT", "P_ILL", "M_ILL", "F_ILL", "TOT_WORK_P","TOT_WORK_M", "TOT_WORK_F", "MAINWORK_P", "MAINWORK_M", "MAINWORK_F", "MARGWORK_P", "MARGWORK_M",
                                                              "MARGWORK_F", "SC_P", "SC_M","SC_F", "ST_P", "ST_M", "ST_F", "TOT_IRR",   
                                                              "UN_IRR", "CULT_WASTE","DEM_avg","VHCF_NDMA")], by = "UniqueID", all.x = TRUE)
print("Merged data dimensions:")
print(dim(aqua_salinity_surge))

head(aqua_salinity_surge)

# Compare Shape_Area with total_area_ha to check how different they are



#_________________________________________________________________________
# Calculate relevant metrics for aqua and salinity
# ________________________________________________________________________

# ! Get total area from a different source than lulc calculations (which may have projection variations)
# Try "shape_Area" from census boundaries for area

# Convert all areas to same unit to make them comparable (and check against total area) 
aqua_salinity_surge$DryAqua_ha <- aqua_salinity_surge$DryAqua_area_m2 / 10000
aqua_salinity_surge$Aqua_ha <- aqua_salinity_surge$Aqua_area_m2 / 10000
aqua_salinity_surge$Saline_ha <- aqua_salinity_surge$saline_areas_m2 / 10000

aqua_salinity_surge$Aqua_perc <- aqua_salinity_surge$Aqua_ha / aqua_salinity_surge$total_area_ha * 100
aqua_salinity_surge$DryAqua_perc <- aqua_salinity_surge$DryAqua_ha / aqua_salinity_surge$total_area_ha * 100
aqua_salinity_surge$Saline_perc <- aqua_salinity_surge$Saline_ha / (aqua_salinity_surge$total_area_ha - aqua_salinity_surge$Aqua_ha) * 100
# Note: Salinity area is calculated after masking aquaculture areas. This implies salinity perc should be calculated on the area that is remaining after aquaculture lanbd conversions are accounted for 

head(aqua_salinity_surge)

summary(aqua_salinity_surge$Saline_perc)
summary(aqua_salinity_surge$Aqua_perc)

morethan100_Aqua <- subset(aqua_salinity_surge, Aqua_perc > 100)
summary(morethan100_Aqua$Year)
morethan100_Aqua$Aqua_pixels_calc <- morethan100_Aqua$Aqua_ha * 10000 / 900

morethan100_Saline <- subset(aqua_salinity_surge, Saline_perc > 100)
summary(morethan100_Saline$Year)

lessthan0_saline <- subset(subset(aqua_salinity_surge, Saline_perc < 0))
summary(lessthan0_saline$Year)

# 3910 observations have more than 100% aquaculture. This could be owing to differences in area calculation projections. 
# 105 observations also have less than 0 saline%, likely owing to aquaculture area dependence in their calculations. 
# # Option 1: Limit perc max to 100% 
# aqua_salinity_surge$Aqua_perc <- pmin(aqua_salinity_surge$Aqua_perc, 100)
# aqua_salinity_surge$DryAqua_perc <- pmin(aqua_salinity_surge$DryAqua_perc, 100)
# aqua_salinity_surge$Saline_perc <- pmin(aqua_salinity_surge$Saline_perc, 100)
# aqua_salinity_surge$agriculture_percent <- pmin(aqua_salinity_surge$agriculture_percent, 100)
# aqua_salinity_surge$urban_percent <- pmin(aqua_salinity_surge$urban_percent, 100)

# Option 2: Make them NA
aqua_salinity_surge$Aqua_perc <- ifelse(aqua_salinity_surge$Aqua_perc > 100, NA, aqua_salinity_surge$Aqua_perc)
aqua_salinity_surge$DryAqua_perc <- ifelse(aqua_salinity_surge$DryAqua_perc > 100, NA, aqua_salinity_surge$DryAqua_perc)
aqua_salinity_surge$Saline_perc <- ifelse(aqua_salinity_surge$Saline_perc > 100, NA, aqua_salinity_surge$Saline_perc)
aqua_salinity_surge$Saline_perc <- ifelse(aqua_salinity_surge$Saline_perc < 0, NA, aqua_salinity_surge$Saline_perc)
aqua_salinity_surge$agriculture_percent <- ifelse(aqua_salinity_surge$agriculture_percent > 100, NA, aqua_salinity_surge$agriculture_percent)
aqua_salinity_surge$urban_percent <- ifelse(aqua_salinity_surge$urban_percent > 100, NA, aqua_salinity_surge$urban_percent)


# ! Districts do not have images available for 2012, except 1 (TN_8). This district is skewing the average significantly.
# Best to make the Aquaculture NA for these TN_8 2012 observations 
aqua_salinity_surge_2012 <- aqua_salinity_surge[aqua_salinity_surge$Year == 2012 & !is.na(aqua_salinity_surge$Aqua_perc), ]
# Convert these 364 observations into NA
table(aqua_salinity_surge_2012$state_code, aqua_salinity_surge_2012$district_code)

aqua_salinity_surge$Aqua_perc[aqua_salinity_surge$Year == 2012] <- NA
aqua_salinity_surge$DryAqua_perc[aqua_salinity_surge$Year == 2012] <- NA
aqua_salinity_surge$Saline_perc[aqua_salinity_surge$Year == 2012] <- NA #Salinity is already NA for 2012 but just adding to be sure
aqua_salinity_surge$Aqua_ha[aqua_salinity_surge$Year == 2012] <- NA
aqua_salinity_surge$DryAqua_ha[aqua_salinity_surge$Year == 2012] <- NA
aqua_salinity_surge$Saline_ha[aqua_salinity_surge$Year == 2012] <- NA


# ! Similarly, for year 1998, only TN_8 plus some parts of OD have satellite images that are skewing the overall results 
# Best to make the Aquaculture NA for these 2012 observations (TN n=440; OD n = 5692)
aqua_salinity_surge_1998 <- aqua_salinity_surge[aqua_salinity_surge$Year == 1998 & !is.na(aqua_salinity_surge$Aqua_perc), ]
table(aqua_salinity_surge_1998$state_code, aqua_salinity_surge_1998$district_code)

sum(aqua_salinity_surge_1998$Year == 1998 & aqua_salinity_surge_1998$state_code == "TN" & aqua_salinity_surge_1998$district_code == 8)
sum(aqua_salinity_surge_1998$Year == 1998 & aqua_salinity_surge_1998$state_code == "OD")

aqua_salinity_surge$Aqua_perc[aqua_salinity_surge$Year == 1998] <- NA
aqua_salinity_surge$DryAqua_perc[aqua_salinity_surge$Year == 1998] <- NA
aqua_salinity_surge$Saline_perc[aqua_salinity_surge$Year == 1998] <- NA 
aqua_salinity_surge$Aqua_ha[aqua_salinity_surge$Year == 1998] <- NA
aqua_salinity_surge$DryAqua_ha[aqua_salinity_surge$Year == 1998] <- NA
aqua_salinity_surge$Saline_ha[aqua_salinity_surge$Year == 1998] <- NA


# Also test 2011 (aqua areas seem very high)
aqua_salinity_surge_2011 <- aqua_salinity_surge[aqua_salinity_surge$Year == 2011 & !is.na(aqua_salinity_surge$Aqua_perc), ]
# Number of observations (29200 or about 3%)
# check how many from each state and district 
table(aqua_salinity_surge_2011$state_code, aqua_salinity_surge_2011$district_code)
# Almost all districts are in, so leaving it in


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

# # Create another variable - land unsuitable for agriculture as a sum of saline and aquaculture 
# # This is because salne % will automatically reduce as aquaculture increases in a village (those pixels get masked) 
# aqua_salinity_surge$NotforRice_ha <- aqua_salinity_surge$Aqua_ha + aqua_salinity_surge$Saline_ha
# aqua_salinity_surge$NotforRice_perc <- aqua_salinity_surge$Aqua_perc + aqua_salinity_surge$Saline_perc

aqua_salinity_surge$State <- aqua_salinity_surge$state_code
aqua_salinity_surge$District <- aqua_salinity_surge$district_code

#_________________________________________________________________________
# Clean database and save
# ________________________________________________________________________

aqua_salinity_surge <- aqua_salinity_surge %>% select(Row_no, UniqueID, Year, State, District, total_area_ha, Shape_Area, NEAR_DIST, Sea_Dist, DEM_avg, Aqua_ha, DryAqua_ha, Saline_ha, Aqua_perc, DryAqua_perc, Saline_perc, pct_change_aqua, pct_change_dryaqua, pct_change_salinity, postSurge, Saline, Saline_Storm_Category, 
                                                      agriculture_area_ha, agriculture_percent, urban_area_ha, urban_percent, TRU, No_HH, TOT_P, TOT_M, TOT_F, P_LIT, M_LIT, F_LIT, P_ILL, M_ILL, F_ILL, TOT_WORK_P,TOT_WORK_M, TOT_WORK_F, MAINWORK_P, MAINWORK_M, MAINWORK_F, MARGWORK_P, MARGWORK_M,
                                                      MARGWORK_F, SC_P, SC_M,SC_F, ST_P, ST_M, ST_F, TOT_IRR,UN_IRR,CULT_WASTE,VHCF_NDMA)
head(aqua_salinity_surge)


write.csv(aqua_salinity_surge, "outputs/aqua_salinity_surge_1990-2025.csv")


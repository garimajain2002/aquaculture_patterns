# Step 1 - create long data form with all uniqueIDs and years
# Step 2 - Combine salinity (X1) data with this (already in long form) 
# Step 3 - Combine surge (X2) data with this (already converted to long form in StormSurgeLongForm)
# Step 4 - Combine aquaculture and dry aquaculture areas (Y) 
# Step 5 - Merge LULC, DEM, Sea distance and Census data to this (Controls) 

# Convert the variables to usable values (absolute to %s etc.) 


library(dplyr)
library(tidyr)
library(zoo) #for rollmean() for calculating smoother percentage changes over 3-5 years
library(slider)

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
# ! Many non-urban areas are likely NA. Match with TRU from Census and convert to 0 if Rural and Urban_pct is NA 

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
# Same just in different units

# ! Many non-urban areas are likely NA. Match with TRU from Census and convert to 0 if Rural and Urban_pct is NA 
aqua_salinity_surge$urban_area_ha <- ifelse(aqua_salinity_surge$TRU=="Rural" & is.na(aqua_salinity_surge$urban_area_ha), 0, aqua_salinity_surge$urban_area_ha)
aqua_salinity_surge$urban_percent <- aqua_salinity_surge$urban_area_ha / aqua_salinity_surge$total_area_ha *100

head(aqua_salinity_surge)





#_________________________________________________________________________
# Calculate relevant metrics for aqua and salinity
# ________________________________________________________________________

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
  mutate(
    Year = as.integer(Year),
    across(c(Aqua_ha, DryAqua_ha, Saline_ha), as.numeric)
  ) %>%
  arrange(UniqueID, Year) %>%
  group_by(UniqueID) %>%
  # Join the previous year's value manually for each metric
  mutate(
    prev_year = Year - 1,
    prev_aqua = Aqua_ha[match(prev_year, Year)],
    prev_dryaqua = DryAqua_ha[match(prev_year, Year)],
    prev_saline = Saline_ha[match(prev_year, Year)],
    
    pct_change_aqua = case_when(
      is.na(prev_aqua) ~ NA_real_,
      prev_aqua == 0 & Aqua_ha == 0 ~ 0,
      prev_aqua == 0 & Aqua_ha != 0 ~ NA_real_,
      TRUE ~ ((Aqua_ha - prev_aqua) / prev_aqua) * 100
    ),
    pct_change_dryaqua = case_when(
      is.na(prev_dryaqua) ~ NA_real_,
      prev_dryaqua == 0 & DryAqua_ha == 0 ~ 0,
      prev_dryaqua == 0 & DryAqua_ha != 0 ~ NA_real_,
      TRUE ~ ((DryAqua_ha - prev_dryaqua) / prev_dryaqua) * 100
    ),
    pct_change_salinity = case_when(
      is.na(prev_saline) ~ NA_real_,
      prev_saline == 0 & Saline_ha == 0 ~ 0,
      prev_saline == 0 & Saline_ha != 0 ~ NA_real_,
      TRUE ~ ((Saline_ha - prev_saline) / prev_saline) * 100
    ),
    
    # Smooth percentage change using a sliding average
    smooth_pct_change_aqua_3 = slide_dbl(pct_change_aqua, ~mean(.x, na.rm = TRUE), .before = 2, .complete = TRUE),
    smooth_pct_change_aqua_5 = slide_dbl(pct_change_aqua, ~mean(.x, na.rm = TRUE), .before = 4, .complete = TRUE),
  
    smooth_pct_change_saline_3 = slide_dbl(pct_change_salinity, ~mean(.x, na.rm = TRUE), .before = 2, .complete = TRUE),
    smooth_pct_change_saline_5 = slide_dbl(pct_change_salinity, ~mean(.x, na.rm = TRUE), .before = 4, .complete = TRUE)
    
    
    ) %>%
  ungroup()

table(aqua_salinity_surge$pct_change_aqua)
summary(aqua_salinity_surge$pct_change_aqua)

table(aqua_salinity_surge$pct_change_dryaqua)
summary(aqua_salinity_surge$pct_change_dryaqua)

table(aqua_salinity_surge$pct_change_salinity)
summary(aqua_salinity_surge$pct_change_salinity)


# Add a flag for villages showing erratic changes 
# Step 1: Flag rows with large jump
aqua_salinity_surge <- aqua_salinity_surge %>%
  mutate(flag_large_jump = abs(pct_change_aqua) > 100)

# Step 2: Subset rows with the flag, and keep year, state, district for inspection
village_flags <- aqua_salinity_surge %>%
  filter(flag_large_jump == TRUE) %>%
  select(UniqueID, Year, State, District, Aqua_ha, pct_change_aqua)

table(aqua_salinity_surge$flag_large_jump, aqua_salinity_surge$Year, aqua_salinity_surge$State)

# Frequency by district (if available)
table(village_flags$State, village_flags$District)


# ________________________________________________________________________________
# Assess and address the level of discontinuity in landsat 5 vs landsat 8 periods 
# ________________________________________________________________________________
# Salinity 
lm_discont <- lm(Saline_perc ~ Year + I(Year >= 2013), data = aqua_salinity_surge)
summary(lm_discont)
# This captures salinity percentage (Saline_perc) as a function of the year (to capture the trend) and a binary dummy (I(Year >= 2013)) that captures the effect of the Landsat sensor switch.
# There's a strong structural break in the salinity values starting in 2013, likely due to the switch from Landsat 5 to Landsat 8 despite the harmonisation. 
# Interpretation: Salinity declines by about 0.1 percentage point per year on average, all else equal.
# After 2013, salinity values drop by 4.32 percentage points on average, even after controlling for the time trend.
# The Landsat 5 and 8 models are not on the same scale.

# Normalize Salinity by Sensor Period
# This will remove level differences between periods, but will retain within-period variation. 
aqua_salinity_surge <- aqua_salinity_surge %>%
  mutate(period = ifelse(Year <= 2012, "L5", "L8")) %>%
  group_by(period) %>%
  mutate(Saline_perc_norm = scale(Saline_perc)) %>%
  ungroup()

summary(aqua_salinity_surge$Saline_perc_norm)
# Note: Values are also -ve here because we have used z-score standardization, which centers the data around 0 within each sensor period.
# Z = (x-mean)/SD
# This is a standardised measure centered at 0. 
# Regression interpretation would be: A one standard deviation increase in salinity (relative to the Landsat period's mean) is associated with a x percentage point increase/decrease in aquaculture land area percentage, holding storm and year constant.


# Create Saline binary 
median_saline <- median(aqua_salinity_surge$Saline_perc_norm, na.rm = TRUE)
aqua_salinity_surge$Saline <- ifelse(aqua_salinity_surge$Saline_perc_norm >= median_saline_norm, 1, 0)


# Aquaculture 
lm_discont <- lm(Aqua_perc ~ Year + I(Year >= 2013), data = aqua_salinity_surge)
summary(lm_discont)
# Interpretation: A statistically significant but very small upward trend over time in aquaculture land share (0.0034 percentage point per year).
# The coefficient for the sensor break (post-2012) is small and not statistically significant. That means: There is no evidence that aquaculture values jump discontinuously at 2013.
# R-squared: ~0.00006 (negligible): The model explains almost none of the variation in Aqua_perc. This is expected, because you're only modeling it with year and a sensor-break dummy — no real predictors yet (like salinity or storm).
# Aquaculture does not have a discontinuity at 2013. 


# TEST 
yearly_summary <- aqua_salinity_surge %>%
  group_by(Year) %>%
  summarize(
    mean_aquaculture = mean(Aqua_perc, na.rm = TRUE),
    mean_salinity = mean(Saline_perc_norm, na.rm = TRUE),
    storm_affected_villages = sum(postSurge, na.rm = TRUE),
    n_villages = n_distinct(UniqueID)
  )
print(yearly_summary)

ggplot(yearly_summary, aes(x = Year, y = mean_salinity)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Saline Area Percentage Points by Year",
       x = "Year", y = "Average Saline Area (%)") +
  scale_x_continuous(breaks = seq(min(yearly_summary$Year), max(yearly_summary$Year), by = 1)) +  # Ensure integer year labels
  
  theme_minimal()


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

aqua_salinity_surge$State <- aqua_salinity_surge$state_code
aqua_salinity_surge$District <- aqua_salinity_surge$district_code


#_________________________________________________________________________
# Clean database and save
# ________________________________________________________________________

aqua_salinity_surge <- aqua_salinity_surge %>% select(Row_no, UniqueID, Year, State, District, total_area_ha, NEAR_DIST, Sea_Dist, DEM_avg, Aqua_ha, DryAqua_ha, Saline_ha, Aqua_perc, DryAqua_perc, Saline_perc_norm, pct_change_aqua, pct_change_dryaqua, pct_change_salinity, smooth_pct_change_aqua_3, smooth_pct_change_aqua_5, smooth_pct_change_saline_3, smooth_pct_change_saline_5, flag_large_jump, postSurge, Saline, Saline_Storm_Category, 
                                                      agriculture_area_ha, agriculture_percent, urban_area_ha, urban_percent, TRU, No_HH, TOT_P, TOT_M, TOT_F, P_LIT, M_LIT, F_LIT, P_ILL, M_ILL, F_ILL, TOT_WORK_P,TOT_WORK_M, TOT_WORK_F, MAINWORK_P, MAINWORK_M, MAINWORK_F, MARGWORK_P, MARGWORK_M,
                                                      MARGWORK_F, SC_P, SC_M,SC_F, ST_P, ST_M, ST_F, TOT_IRR,UN_IRR,CULT_WASTE,VHCF_NDMA)
head(aqua_salinity_surge)


write.csv(aqua_salinity_surge, "outputs/aqua_salinity_surge_1990-2025.csv")


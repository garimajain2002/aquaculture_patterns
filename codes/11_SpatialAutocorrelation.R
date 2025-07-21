# Hypothesis - Aquaculture is spatially clustering over time 
# Implication - the regression analysis will have to account for spatial clustering to address spatial autocorrelation 

# Conduct a Moran's I test for all districts for all years and see how the trend changes 
# If positive in 2025, use a spatial lag or spatial error model  

getwd()

library(sf)
library(spdep)
library(ggplot2)
library(dplyr)

# 1. Load shapefile (for one district to start with - combine two sources of aqua area calculations Landsat 8 period and Landsat 5 period)
gdb_path <- "C:/Users/garim/Desktop/ArcGISPro/Ch1/Ch1_Patterns/Ch1_Patterns.gdb"
st_layers(gdb_path)

TN_12_Villages <- st_read(dsn = gdb_path, layer = "TN_12_Villages")
TN_12_AquaStats1 <- st_read (dsn = gdb_path, layer = "TN_12_AquaStats")
TN_12_AquaStats2 <- st_read (dsn = gdb_path, layer = "TN_12_AquaStats_Landsat5")

TN_12_AllData <- TN_12_Villages %>%
  left_join(TN_12_AquaStats1, by = "UniqueID") %>%
  left_join(TN_12_AquaStats2, by = "UniqueID")

st_write(TN_12_AllData, "data/shp/TN_12_AllData.shp", delete_layer = TRUE)

villages_sf <- TN_12_AllData

summary(villages_sf$Shape_Area)

# 2. Project to UTM (for accurate area & distance)
villages_sf <- st_transform(villages_sf, crs = 32644)  

# 3. Calculate area in square meters
villages_sf$Shape_Area <- st_area(villages_sf)

summary(villages_sf$Shape_Area)


# 4. Calculate Aqua_perc_<year> for each year
for (yr in 1990:2025) {
  varname_abs <- paste0("Aqua_", yr)
  varname_pct <- paste0("Aqua_perc_", yr)
  
  if (varname_abs %in% colnames(villages_sf) && !all(is.na(villages_sf[[varname_abs]]))) {
    villages_sf[[varname_pct]] <- as.numeric(villages_sf[[varname_abs]]) / as.numeric(villages_sf$Shape_Area)
  } else {
    villages_sf[[varname_pct]] <- NA_real_
  }
}

# Drop erroneous %
for (yr in 1990:2025) {
  varname <- paste0("Aqua_perc_", yr)
  
  if (varname %in% names(villages_sf)) {
    n_bad <- sum(villages_sf[[varname]] > 1, na.rm = TRUE)
    message(sprintf("Year %s: %d values > 100%% were set to NA", yr, n_bad))
    villages_sf[[varname]][villages_sf[[varname]] > 1] <- NA
  }
}

# 5. Define spatial neighbors: distance-based within 5 km
coords <- st_coordinates(st_centroid(villages_sf))
nb <- dnearneigh(coords, 0, 5000)  # 0 to 5 km
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)    # Row-standardized weights, allow empty neighbors

table(card(nb) == 0)  # card() = number of neighbors

# 6. Compute Moran’s I for each year and store results
moran_results <- data.frame(
  year = 1990:2025,
  moran_I = NA_real_,
  p_value = NA_real_
)

for (yr in 1990:2025) {
  varname <- paste0("Aqua_perc_", yr)
  values <- villages_sf[[varname]]
  
  # Create logical vector of non-NA values
  valid_logical <- !is.na(values)
  
  if (sum(valid_logical) > 2) {  # Must have at least 3 values to compute Moran's I
    # Subset the spatial data to valid rows
    villages_valid <- villages_sf[valid_logical, ]
    values_clean <- values[valid_logical]
    
    # Recompute coordinates and neighbors for this subset
    coords_valid <- st_coordinates(st_centroid(villages_valid))
    nb_valid <- dnearneigh(coords_valid, 0, 5000)  # or whatever distance
    lw_valid <- nb2listw(nb_valid, style = "W", zero.policy = TRUE)
    
    # Now Moran’s I will work
    moran_test <- moran.test(values_clean, lw_valid, zero.policy = TRUE)
    
    moran_results$moran_I[moran_results$year == yr] <- moran_test$estimate[["Moran I statistic"]]
    moran_results$p_value[moran_results$year == yr] <- moran_test$p.value
  }
}

# 7. View results
print(moran_results)

# 8. Plot Moran’s I over time
ggplot(moran_results, aes(x = year, y = moran_I)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(aes(color = p_value < 0.05), size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 1)+
  scale_color_manual(values = c("grey60", "red"), labels = c("Not Significant", "p < 0.05")) +
  labs(
    title = "Moran’s I for % Aquaculture Area (1990–2025) in TN_12",
    x = "Year", y = "Moran’s I", color = "Significance"
  ) +
  theme_minimal()

ggsave("outputs/MoransI/MoransI_TN_12.png", width = 8, height = 6, units = "in", dpi = 300)


# 9. Moran scatterplot for a specific year (e.g., 2025) - first account for NA
# Get values for 2025
values <- villages_sf$Aqua_perc_2025

# Identify non-missing observations
valid_logical <- !is.na(values)

# Subset both values and spatial weights
values_clean <- values[valid_logical]
lw_clean <- subset(lw, subset = valid_logical)

# Plot Moran scatterplot
moran.plot(values_clean, lw_clean, zero.policy = TRUE)



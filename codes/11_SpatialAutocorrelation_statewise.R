# Hypothesis - Aquaculture is spatially clustering over time 
# Implication - the regression analysis will have to account for spatial clustering to address spatial autocorrelation 

# Conduct a Moran's I test for all states for all years and see how the trend changes 
# If positive (in 2025), use a spatial lag or spatial error model  

getwd()

library(sf)
library(spdep)
library(ggplot2)
library(dplyr)

# 1. Load shapefiles (use a loop to read all districts and then bind for each state)
village_geometry <- st_read("data/shp/States_AllData_geometry.shp")
village_data <- read.csv(unz("data/aqua_salinity_surge_1990-2025.zip", 
                             "aqua_salinity_surge_1990-2025.csv"))

# Convert from long form to wide form first 
head(village_data)
village_data_wide <- village_data %>%
  select(
    UniqueID, Year, State, Shape_Area, Longitude, Latitude,
    Aqua_ha, DryAqua_ha, Saline_ha,
    Aqua_perc, DryAqua_perc, Saline_perc_norm, avg_salinity_5yr
  ) %>%
  pivot_wider(
    names_from = Year,
    values_from = c(
      Aqua_ha, DryAqua_ha, Saline_ha,
      Aqua_perc, DryAqua_perc, Saline_perc_norm, avg_salinity_5yr
    ),
    names_sep = "_"
  )

head(village_data_wide)


head(village_geometry) # shape area and lat long are already in village_data so drop those before joining 
village_geometry <- village_geometry %>%
  select(UniqueID, geometry)


# merge the two to get all data with geometry 
village_all <- left_join(village_data_wide, village_geometry, by = "UniqueID")
village_all <- st_as_sf(village_all)
print(st_geometry_type(village_all))
st_crs(village_all)

village_all <- st_transform(village_all, crs = 32644)  # UTM Zone 44N in m
summary(village_all)


# make state subsets 
village_AP <- village_all %>% filter(State == "AP")
village_TN <- village_all %>% filter(State == "TN")
village_OD <- village_all %>% filter(State == "OD")  

# Check if Ok 
summary(village_OD) 

# 1. Pick all, or per state for Moran's I 
#--------------------------------------
villages_sf <- village_all

head(villages_sf)
summary(villages_sf$Shape_Area)


# 2. Use centroids to define coordinates
coords <- st_coordinates(st_centroid(villages_sf))

# 3. Define distance-based neighbors (within 5 km)
nb <- dnearneigh(coords, d1 = 0, d2 = 5000)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Optional: Check how many villages have 0 neighbors
cat("Villages with no neighbors (within 5km):", sum(card(nb) == 0), "\n")

# 4. Compute Moran's I across years
moran_results <- data.frame(
  year = 1990:2025,
  moran_I = NA_real_,
  p_value = NA_real_
)

for (yr in 1990:2025) {
  varname <- paste0("Aqua_perc_", yr)
  
  if (varname %in% names(villages_sf)) {
    values <- villages_sf[[varname]]
    valid <- !is.na(values)
    
    if (sum(valid) > 2 && var(values[valid]) > 0) {
      lw_valid <- subset(lw, subset = valid, zero.policy = TRUE)
      moran_test <- moran.test(values[valid], lw_valid, zero.policy = TRUE)
      
      moran_results$moran_I[moran_results$year == yr] <- moran_test$estimate[["Moran I statistic"]]
      moran_results$p_value[moran_results$year == yr] <- moran_test$p.value
    } else {
      message(sprintf("Year %s: No valid data or zero variance", yr))
    }
  }
}

# 5. View result
print(moran_results)

summary(moran_results$p_value)
table(moran_results$p_value < 0.05, useNA = "ifany")

moran_results <- moran_results %>%
  mutate(significance = case_when(
    is.na(p_value) ~ "Missing",
    p_value < 0.05 ~ "p < 0.05",
    TRUE ~ "Not Significant"
  ))


# 7. Plot Moran’s I over time
ggplot(moran_results, aes(x = year, y = moran_I)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(aes(color = significance), size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 1) +
  scale_color_manual(
    values = c("p < 0.05" = "red", "Not Significant" = "grey60", "Missing" = "black"),
    drop = FALSE
  ) +
  labs(
    title = "Moran’s I for % Aquaculture Area (1990–2025) in Coastal India (OD, AP, TN)",
    x = "Year", y = "Moran’s I", color = "Significance"
  ) +
  theme_minimal()

ggsave("outputs/MoransI/MoransI_OD.png", width = 8, height = 6, units = "in", dpi = 300)


# 8. Moran scatterplot for a specific year (e.g., 2025) - first account for NA
# Get values for 2025
values <- villages_sf$Aqua_perc_2025

# Identify non-missing observations
valid_logical <- !is.na(values)

# Subset both values and spatial weights
values_clean <- values[valid_logical]
lw_clean <- subset(lw, subset = valid_logical)

# Plot Moran scatterplot
moran.plot(values_clean, lw_clean, zero.policy = TRUE)

 



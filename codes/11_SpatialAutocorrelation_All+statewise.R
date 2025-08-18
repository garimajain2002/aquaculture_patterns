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
compute_moran_over_time <- function(village_data, lw, region) {
  out <- data.frame(
    year = 1990:2025,
    moran_I = NA_real_,
    p_value = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    State = region
  )
  
  for (yr in 1990:2025) {
    varname <- paste0("Aqua_perc_", yr)
    if (varname %in% names(village_data)) {
      values <- village_data[[varname]]
      valid <- !is.na(values)
      
      if (sum(valid) > 2 && var(values[valid]) > 0) {
        lw_valid <- subset(lw, subset = valid, zero.policy = TRUE)
        moran_test <- moran.test(values[valid], lw_valid, zero.policy = TRUE)
        
        out$moran_I[out$year == yr] <- moran_test$estimate[["Moran I statistic"]]
        out$p_value[out$year == yr] <- moran_test$p.value
        
        out$lower[out$year == yr] <- out$moran_I[out$year == yr] - 0.02
        out$upper[out$year == yr] <- out$moran_I[out$year == yr] + 0.02
      }
    }
  }
  
  out <- out %>%
    mutate(significance = case_when(
      is.na(p_value) ~ "Missing",
      p_value < 0.05 ~ "p < 0.05",
      TRUE ~ "Not Significant"
    ))
  
  return(out)
}

# Buld listw for each state 
build_weights <- function(villages_sf, dmax = 5000) {
  coords <- st_coordinates(st_centroid(villages_sf))
  nb <- dnearneigh(coords, d1 = 0, d2 = dmax)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  return(lw)
}

# -------- Build weights and compute Moran’s I for each region --------
# ALL INDIA
lw_all <- build_weights(village_all)
moran_all <- compute_moran_over_time(village_all, lw_all, "India")

# AP
lw_AP <- build_weights(village_AP)
moran_ap <- compute_moran_over_time(village_AP, lw_AP, "AP")

# TN
lw_TN <- build_weights(village_TN)
moran_tn <- compute_moran_over_time(village_TN, lw_TN, "TN")

# OD
lw_OD <- build_weights(village_OD)
moran_od <- compute_moran_over_time(village_OD, lw_OD, "OD")

# -------- Combine all results --------
moran_all_states <- bind_rows(moran_all, moran_ap, moran_tn, moran_od)


# 5. Plot Moran’s I over time

moran_all_states$State <- factor(
  moran_all_states$State,
  levels = c("India", "AP", "TN", "OD")
)

ggplot(moran_all_states, aes(x = year, y = moran_I, color = State, linetype = State)) +
  # Confidence intervals
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = State), alpha = 0.2, color = NA) +
  
  # Lines
  geom_line(
    data = subset(moran_all_states, !is.na(moran_I)),
    aes(group = State), size = 1
  ) +
  
  # Points for significance
  geom_point(
    data = subset(moran_all_states, significance != "Missing"),
    aes(shape = significance),
    size = 2.5
  ) +
  
  # Custom color and fill
  scale_color_manual(values = c("India" = "black", "AP" = "darkgreen", "TN" = "darkorange", "OD" = "steelblue")) +
  scale_fill_manual(values = c("India" = "black", "AP" = "darkgreen", "TN" = "darkorange", "OD" = "steelblue")) +
  
  # Custom linetype
  scale_linetype_manual(values = c("India" = "solid", "AP" = "dashed", "TN" = "dashed", "OD" = "dashed")) +
  
  # Point shape for p-value significance
  scale_shape_manual(values = c("p < 0.05" = 16, "Not Significant" = 1, "Missing" = NA), drop = FALSE) +
  
  labs(
    title = "Moran’s I for % Aquaculture Area (1990–2025)",
    subtitle = "Solid line: India | Dashed lines: States | Filled points = p < 0.05",
    x = "Year", y = "Moran’s I",
    color = "Region", linetype = "Region", fill = "Region", shape = "Significance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggsave("outputs/MoransI/MoransI_All.png", width = 12, height = 8, units = "in", dpi = 300)



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

 



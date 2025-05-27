library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(purrr)
library(exactextractr)

# Define state codes and districts
state_codes <- c("OD", "AP", "TN")
districts <- list(
  OD = 1:6,
  AP = 1:9,
  TN = 1:12
)

# Define target CRS (WGS 84 / UTM zone for India, or EPSG:4326 as default)
target_crs <- sf::st_crs(4326)

# Function to calculate percentage of agriculture and urban area
calculate_lulc_percentages <- function(state_code, dist_code) {
  boundary_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  lulc_path <- sprintf("data/tif/LandCover_dtype_compress/%s_%d_LC.tif", state_code, dist_code)
  
  cat(sprintf("Processing %s district %d...\n", state_code, dist_code))
  
  tryCatch({
    villages <- sf::st_read(boundary_path, quiet = TRUE)
    lulc <- raster::raster(lulc_path)
    
    # Reproject both to a consistent CRS
    villages <- st_transform(villages, target_crs)
    
    # Reproject LULC using a template for alignment and nearest neighbor
    lulc_proj <- raster::projectRaster(
      lulc,
      crs = as.character(target_crs$wkt),
      method = "ngb"
    )
    
    # Get pixel resolution in meters (assumes projected CRS)
    res_m <- res(lulc_proj)
    pixel_area_m2 <- res_m[1] * res_m[2]  # m²
    pixel_area_ha <- pixel_area_m2 / 10000  # convert to hectares
    
    # Exact extraction with fractional area * pixel area to get hectares
    stats <- exact_extract(lulc_proj, villages, function(values, coverage_fractions) {
      df <- data.frame(class = values, frac = coverage_fractions) %>%
        filter(!is.na(class)) %>%
        group_by(class) %>%
        summarize(area_ha = sum(frac) * pixel_area_ha) %>%
        pivot_wider(names_from = class, values_from = area_ha, values_fill = 0)
    }, progress = FALSE)
    
    stats$UniqueID <- villages$UniqueID
    stats <- as.data.frame(stats)
    
    # Compute totals and percentages
    stats <- stats %>%
      mutate(
        total_area_ha = rowSums(select(., -UniqueID), na.rm = TRUE),
        agriculture_area_ha = ifelse("40" %in% names(.), .[["40"]], 0),
        urban_area_ha = ifelse("50" %in% names(.), .[["50"]], 0),
        agriculture_percent = (agriculture_area_ha / total_area_ha) * 100,
        urban_percent = (urban_area_ha / total_area_ha) * 100,
        state_code = state_code,
        district_code = dist_code
      ) %>%
      select(state_code, district_code, UniqueID, total_area_ha,
             agriculture_area_ha, urban_area_ha,
             agriculture_percent, urban_percent)
    
    return(stats)
  }, error = function(e) {
    cat(sprintf("Error processing %s district %d: %s\n", state_code, dist_code, e$message))
    return(NULL)
  })
}

# Run for all districts
results_list <- list()

for (state_code in state_codes) {
  for (dist_code in districts[[state_code]]) {
    result <- calculate_lulc_percentages(state_code, dist_code)
    if (!is.null(result)) {
      results_list[[length(results_list) + 1]] <- result
    }
  }
}

# Combine results and ensure uniqueness
all_results <- bind_rows(results_list)

# Drop duplicates if any
all_results <- all_results %>% distinct(UniqueID, .keep_all = TRUE)

# Summary statistics
cat("\nSummary statistics:\n")
summary_stats <- all_results %>%
  group_by(state_code, district_code) %>%
  summarize(
    total_villages = n(),
    mean_agriculture_percent = mean(agriculture_percent, na.rm = TRUE),
    mean_urban_percent = mean(urban_percent, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# Save output
write.csv(all_results, "outputs/village_lulc_percentages_exact.csv", row.names = FALSE)
write.csv(summary_stats, "outputs/village_lulc_summarystats_exact.csv", row.names = FALSE)

# Check extreme values
morethan100_exact <- subset(all_results, agriculture_percent > 100 | urban_percent > 100)

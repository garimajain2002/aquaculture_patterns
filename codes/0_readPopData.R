# Extract 1990-2025 population data for all villages 

library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)
library(zoo)

getwd()

shp <- st_read("data/shp/States_AllData_geometry.shp")
villages <- st_transform(shp, crs = 4326) # should be same as .nc data 

base_raster_dir <- "data/population_GHS/"  # Parent folder with all year subfolders

# Define years (5 year steps)
years <- seq(1990, 2025, by = 5)


# Initialize output list
pop_list <- list()


# Loop through years
for (year in years) {
  # List the three folders for this year
  year_pattern <- paste0("GHS_POP_E", year)
  year_folders <- list.files(base_raster_dir, pattern = year_pattern, full.names = TRUE)
  
  # Find TIFF files inside each folder
  tif_paths <- list.files(year_folders, pattern = "\\.tif$", full.names = TRUE)
  
  # Read and merge rasters
  raster_list <- lapply(tif_paths, rast)
  merged_raster <- do.call(merge, raster_list)
  
  # Reproject villages if needed
  if (!st_crs(villages) == crs(merged_raster)) {
    villages_proj <- st_transform(villages, crs(merged_raster))
  } else {
    villages_proj <- villages
  }
  
  # Extract total population per village
  pop_data <- exact_extract(merged_raster, villages_proj, 'sum')
  
  # Combine with village attributes
  villages_year <- villages_proj %>%
    st_drop_geometry() %>%
    mutate(Year = year,
           Population = pop_data)
  
  pop_list[[as.character(year)]] <- villages_year
  cat("Finished processing year", year, "\n")
}



# Combine all years
village_population <- bind_rows(pop_list)


# Linear interpolation of missing years (using tidyr::complete()) 
all_years <- data.frame(Year = 1990:2025)

village_population_full <- village_population %>%
  group_by(UniqueID) %>%
  complete(Year = 1990:2025) %>%
  arrange(UniqueID, Year) %>%
  # Linear interpolation
  mutate(Population = zoo::na.approx(Population, x = Year, na.rm = FALSE)) %>%
  ungroup()


village_population_full <- village_population_full %>% 
  select("UniqueID", "Year", "Population")


# Export
write.csv(village_population_full, "data/population_GHS/village_population_timeseries.csv", row.names = FALSE)



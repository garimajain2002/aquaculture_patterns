# Extract rainfall data for all years for all districts/villages 

library(ncdf4)
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(tidyr)


getwd()

shp <- st_read("data/shp/States_AllData_geometry.shp")
villages <- st_transform(shp, crs = 4326) # should be same as .nc data 

# Test reading one .nc file 
data_1990 <- nc_open("data/climate_data_IMD/Yearly_Gridded_Rainfall_25/RF25_ind1990_rfp25.nc")
print(data_1990)

# This contains daily RAINFALL (in mm) values Time = 365 
# There’s a single variable with dimensions [LONGITUDE, LATITUDE, TIME]


# List all .25° NetCDF files
nc_files <- list.files("data/climate_data_IMD/Yearly_Gridded_Rainfall_25/", 
                       pattern = "\\.nc$", full.names = TRUE)

basename(nc_files)

gsub(".*?(\\d{4})[^0-9]*\\.nc$", "\\1", "RF25_ind1990_rfp25.nc")


# Helper function to get year from filename
extract_year <- function(path) {
  as.integer(sub(".*?(\\d{4}).*", "\\1", basename(path)))
}

# test if this is working 
extract_year("RF25_ind1990_rfp25.nc")

rainfall_list <- list()


# Loop through NetCDFs
for (i in seq_along(nc_files)) {
  year <- extract_year(nc_files[i])
  cat("Processing year:", year, "\n")
  
  # Read daily rainfall stack
  r <- rast(nc_files[i])
  r[r == -999] <- NA  # clean fill values
  
  # Sum over 365 days to get total rainfall (in mm/year)
  r_total <- sum(r, na.rm = TRUE)
  
  # Extract mean rainfall per village
  rainfall_vals <- exact_extract(r_total, villages, 'mean')
  
  rainfall_list[[i]] <- data.frame(
    UniqueID = villages$UniqueID,
    Year = as.integer(year),
    Rainfall_mm = rainfall_vals
  )
}

# 2025 data is cy unavailable. Repeat 2024 data. 
rainfall_2024 <- rainfall_list[[which(sapply(rainfall_list, function(x) 2024 %in% x$Year))]]
rainfall_2025 <- rainfall_2024
rainfall_2025$Year <- 2025
rainfall_list[[length(rainfall_list) + 1]] <- rainfall_2025

# Combine all into one long dataframe
rainfall_long <- bind_rows(rainfall_list)

head(rainfall_long)

summary(rainfall_long$Rainfall_mm)

# Save to CSV if needed
write.csv(rainfall_long, "data/climate_data_IMD/village_yearly_rainfall_long.csv", row.names = FALSE)





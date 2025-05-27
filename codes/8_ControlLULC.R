# Calculate % agriculture and % urban per village for a ref. year (2015) 

library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(purrr)

getwd()

# Define state codes and districts
state_codes <- c("OD", "AP", "TN")
districts <- list(
  OD = 1:6,
  AP = 1:9,
  TN = 1:12
)

# Function to calculate percentage of agriculture and urban area for each village
calculate_lulc_percentages <- function(state_code, dist_code) {
  # Construct file paths
  boundary_path <- sprintf("data/shp/%s_%d_AllData.shp", state_code, dist_code)
  lulc_path <- sprintf("data/tif/LandCover_dtype_compress/%s_%d_LC.tif", state_code, dist_code)
  
  # Read shapefile and LULC raster
  cat(sprintf("Processing %s district %d...\n", state_code, dist_code))
  
  tryCatch({
    villages <- sf::st_read(boundary_path, quiet = TRUE)
    lulc <- raster::raster(lulc_path)
    
    # Check and align projections
    cat(sprintf("Checking projections for %s district %d...\n", state_code, dist_code))
    villages_crs <- st_crs(villages)
    lulc_crs <- crs(lulc)
    
    print(res(lulc)) 
    print(crs(lulc))  
    
    cat(sprintf("Village CRS: %s\n", as.character(villages_crs$wkt)))
    cat(sprintf("LULC CRS: %s\n", as.character(lulc_crs)))
    
    # Align projections if they don't match
    if (!compareCRS(villages_crs, lulc_crs)) {
      cat("Projections don't match. Reprojecting LULC to match village boundaries...\n")
      lulc <- projectRaster(lulc, crs = villages_crs$wkt, method = "ngb", res = res(lulc))
      cat("Reprojection complete.\n")
    } else {
      cat("Projections match. No reprojection needed.\n")
    }
    
    result <- data.frame(
      state_code = character(),
      district_code = integer(),
      UniqueID = character(),
      total_area_ha = numeric(),
      agriculture_area_ha = numeric(),
      urban_area_ha = numeric(),
      agriculture_percent = numeric(),
      urban_percent = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:nrow(villages)) {
      village <- villages[i, ]
      
      unique_id <- if ("UniqueID" %in% names(village)) as.character(village$UniqueID) 
      else if ("UNIQUE_ID" %in% names(village)) as.character(village$UNIQUE_ID) 
      else if ("unique_id" %in% names(village)) as.character(village$unique_id) 
      else paste0(state_code, "_", dist_code, "_", i)
      
      village_boundary <- village
      total_area_ha <- sf::st_area(village_boundary) / 10000
      
      if (sf::st_is_valid(village_boundary)) {
        village_lulc <- tryCatch({
          ext <- extent(as(village_boundary, "Spatial"))
          raster::crop(lulc, ext)
        }, error = function(e) {
          cat(sprintf("Error cropping raster for village %s: %s\n", unique_id, e$message))
          return(NULL)
        })
        
        if (!is.null(village_lulc)) {
          village_lulc <- tryCatch({
            raster::mask(village_lulc, as(village_boundary, "Spatial"))
          }, error = function(e) {
            cat(sprintf("Error masking raster for village %s: %s\n", unique_id, e$message))
            return(NULL)
          })
        }
        
        if (!is.null(village_lulc)) {
          lulc_counts <- raster::freq(village_lulc)
          
          agriculture_pixels <- lulc_counts[lulc_counts[, 1] == 40, 2]
          urban_pixels <- lulc_counts[lulc_counts[, 1] == 50, 2]
          
          agriculture_pixels <- ifelse(length(agriculture_pixels) == 0, 0, agriculture_pixels)
          urban_pixels <- ifelse(length(urban_pixels) == 0, 0, urban_pixels)
          
          total_pixels <- sum(lulc_counts[, 2], na.rm = TRUE)
          
          pixel_area_ha <- 1
          
          agriculture_area_ha <- agriculture_pixels * pixel_area_ha
          urban_area_ha <- urban_pixels * pixel_area_ha
          
          agriculture_percent <- (agriculture_area_ha / as.numeric(total_area_ha)) * 100
          urban_percent <- (urban_area_ha / as.numeric(total_area_ha)) * 100
        } else {
          agriculture_area_ha <- 0
          urban_area_ha <- 0
          agriculture_percent <- 0
          urban_percent <- 0
        }
      } else {
        cat(sprintf("Warning: Invalid geometry for village %s. Skipping...\n", unique_id))
        agriculture_area_ha <- 0
        urban_area_ha <- 0
        agriculture_percent <- 0
        urban_percent <- 0
      }
      
      result <- bind_rows(result, data.frame(
        state_code = state_code,
        district_code = dist_code,
        UniqueID = unique_id,
        total_area_ha = as.numeric(total_area_ha),
        agriculture_area_ha = agriculture_area_ha,
        urban_area_ha = urban_area_ha,
        agriculture_percent = agriculture_percent,
        urban_percent = urban_percent,
        stringsAsFactors = FALSE
      ))
    }
    
    return(result)
  }, error = function(e) {
    cat(sprintf("Error processing %s district %d: %s\n", state_code, dist_code, e$message))
    return(NULL)
  })
}

# Process all districts for all states
results_list <- list()

for (state_code in state_codes) {
  for (dist_code in districts[[state_code]]) {
    result <- calculate_lulc_percentages(state_code, dist_code)
    if (!is.null(result)) {
      print(table(duplicated(result$UniqueID)))
      results_list[[length(results_list) + 1]] <- result
    }
  }
}

all_results <- do.call(rbind, results_list)

print("Data dimensions:")
print(dim(all_results))
print(colnames(all_results))
table(duplicated(all_results$UniqueID))
all_results <- all_results %>% distinct(UniqueID, .keep_all = TRUE)

print("Data dimensions:")
print(dim(all_results))

head(all_results)

summary(all_results$agriculture_percent)
summary(all_results$urban_percent)
morethan100 <- subset(all_results, agriculture_percent > 100 | urban_percent > 100)
# About 2300 values are over 100% due to small edge condition errors of pixels crossing village bloundaries. 
# Tried using exactextractr (see 8_ControlLULC_exactextractr for the other version) but has way more errors (>14k values over 100)
# Can fix these 2k values by either checking these total areas against other sources of information or limiting the max perc to 100. 

write.csv(all_results, "outputs/village_lulc_percentages.csv", row.names = FALSE)


cat("\nSummary of results:\n")
summary_stats <- all_results %>%
  group_by(state_code, district_code) %>%
  summarize(
    total_villages = n(),
    mean_agriculture_percent = mean(agriculture_percent, na.rm = TRUE),
    mean_urban_percent = mean(urban_percent, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)
write.csv(summary_stats, "outputs/village_lulc_summarystats.csv", row.names = FALSE)

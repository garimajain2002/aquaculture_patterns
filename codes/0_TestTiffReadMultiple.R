library(raster)
library(terra)
library(dplyr)
library(ggplot2)

getwd()

# Method: 
# Step 1: Create a loop to read TIFFs for districts for different years
# Step 2: Extract Aquaculture mask 
# Step 3: Run Salinity model 
# Step 4: Make the final raster by applying the aquaculture mask to the salinity raster 
# Step 5: Save the final output for each district and year into a GeoTiff and csv (to be read later for the database)

#======================================
# Load appropriate models and assets 
#======================================
# Load Salinity Model
MultiStratEnsemble <- readRDS("codes/ensemble_2501.rds") #percentile
best_threshold <- 0.768

# Define district IDs and years 
# Do this in batches of 5 years and one state at a time
#districts <- c(paste0("OD_", 1:6), paste0("AP_", 1:9), paste0("TN_", 1:12))

districts <- c(paste0("OD_", 1:1))
years <- 2024:2025


# # For Testing
# landsat_image <- stack("data\\GEE_Outputs_Landsat_Composites_CoastalDistricts\\OD_4_Landsat_2024.tif")
# aqua_image <- stack("data/GEE_Outputs_Classified_Maps/OD_4_Classified_2024.tif")




# Loop through each district and year
for (district in districts) {
  for (year in years) {

    # Construct file paths
    landsat_path <- paste0("data/GEE_Outputs_Landsat_Composites_CoastalDistricts/", district, "_Landsat_", year, ".tif")
    aqua_path <- paste0("data/GEE_Outputs_Classified_Maps/", district, "_Classified_", year, ".tif")

    # Check if files exist before proceeding
    if (file.exists(landsat_path) & file.exists(aqua_path)) {

      # Read Landsat and Classified Images
      landsat_image <- stack(landsat_path)
      aqua_image <- stack(aqua_path)


# 1. Mask Aquaculture ponds
      aqua_df <- as.data.frame(aqua_image, xy = TRUE)

      # Check df for any Infinite, or NaN values
      colSums(is.na(aqua_df))   # Check for NA values
      sum(is.infinite(as.matrix(aqua_df)))  # Check for Inf values
      sum(is.nan(as.matrix(aqua_df)))  # Check for NaN values
      
      #plot aquaculture
      predicted_Aqua <- ggplot(aqua_df, aes(x = x, y = y, fill = factor(classification))) +
        geom_raster() +
        scale_fill_manual(values = c("0" = "white", "1" = "darkgrey", "2" = "blue"),
                          labels = c("Other uses", "Dry Aquaculture", "Active Aquaculture"),
                          name = "Aquaculture") +
            ggtitle(paste("Aquaculture -", district, year)) +
        theme_minimal()

      plot(predicted_Aqua)

      #Note: Aqua == 2, Dry Aqua == 1, and others == 0
      unique(values(aqua_image))
      aqua_mask_raster <- aqua_image == 2  # Create a binary mask of aquaculture areas
      

# 2. Prepare landsat dataframe
      landsat_df <- as.data.frame(landsat_image, xy = TRUE)
      landsat_df$ID <- 1:nrow(landsat_df) # Assign IDs before droppping cells
      
      colnames(landsat_df) <- c("x","y","BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2", "NDVI", "NDWI", "NDSI", "ID")
      landsat_df <- landsat_df[complete.cases(landsat_df), ] # Ensure no NA values in landsat_df
      landsat_df <- subset(landsat_df, BLUE != 0) # Drop masked cells or where value is 0 (surface water cells)

      # Calculate additional indices
      # Blue and red
      landsat_df$NBR <- (landsat_df$BLUE-landsat_df$RED)/(landsat_df$BLUE+landsat_df$RED)

      # Blue and green
      landsat_df$NBG <- (landsat_df$BLUE-landsat_df$GREEN)/(landsat_df$BLUE+landsat_df$GREEN)

      # Blue and NIR
      landsat_df$NBNIR <- (landsat_df$BLUE-landsat_df$NIR)/(landsat_df$BLUE+landsat_df$NIR)

      # Blue and SWIR1
      landsat_df$NBSWIR1 <- (landsat_df$BLUE-landsat_df$SWIR1)/(landsat_df$BLUE+landsat_df$SWIR1)

      # Blue and SWIR2
      landsat_df$NBSWIR2 <- (landsat_df$BLUE-landsat_df$SWIR2)/(landsat_df$BLUE+landsat_df$SWIR2)

      # Red and Green (also NDVI)
      landsat_df$NDVI <- (landsat_df$RED-landsat_df$GREEN)/(landsat_df$RED+landsat_df$GREEN)

      # Red and NIR (NDSI2 or Normalised Difference Salinity Index 2 as per as per Khan et al 2001 in Nguyen et al 2020)
      landsat_df$NDSI2 <- (landsat_df$RED-landsat_df$NIR)/(landsat_df$RED+landsat_df$NIR)

      # Red and SWIR1
      landsat_df$NRSWIR1 <- (landsat_df$RED-landsat_df$SWIR1)/(landsat_df$RED+landsat_df$SWIR1)

      # Red and SWIR2
      landsat_df$NRSWIR2 <- (landsat_df$RED-landsat_df$SWIR2)/(landsat_df$RED+landsat_df$SWIR2)

      # Green and NIR (also NDWI)
      landsat_df$NDWI <- (landsat_df$GREEN-landsat_df$NIR)/(landsat_df$GREEN+landsat_df$NIR)

      # Green and SWIR1
      landsat_df$NGSWIR1 <- (landsat_df$GREEN-landsat_df$SWIR1)/(landsat_df$GREEN+landsat_df$SWIR1)

      # Green and SWIR2
      landsat_df$NGSWIR2 <- (landsat_df$GREEN-landsat_df$SWIR2)/(landsat_df$GREEN+landsat_df$SWIR2)

      # NIR and SWIR1
      landsat_df$NNIRSWIR1 <- (landsat_df$NIR-landsat_df$SWIR1)/(landsat_df$NIR+landsat_df$SWIR1)

      # NIR and SWIR2
      landsat_df$NNIRSWIR2 <- (landsat_df$NIR-landsat_df$SWIR2)/(landsat_df$NIR+landsat_df$SWIR2)

      # SWIR1 and SWIR2 (also NDSI as per the Index Database: https://www.indexdatabase.de/db/is.php?sensor_id=168 )
      landsat_df$NDSI1 <- (landsat_df$SWIR1-landsat_df$SWIR2)/(landsat_df$SWIR1+landsat_df$SWIR2)

      # Salinity Index 1 = sqrt(green^2+red^2)
      landsat_df$SI1 <- sqrt((landsat_df$GREEN)^2 + (landsat_df$RED)^2)

      # Salinity Index 2 = sqrt(green x red)
      landsat_df$SI2 <- sqrt(landsat_df$GREEN * landsat_df$RED)

      # Salinity Index 3 = sqrt(blue x red)
      landsat_df$SI3 <- sqrt(landsat_df$BLUE * landsat_df$RED)

      # salinity index 4 = red x NIR / green
      landsat_df$SI4 <- (landsat_df$RED * landsat_df$NIR / landsat_df$GREEN)

      # salinity index 5 = blue/red
      landsat_df$SI5 <- (landsat_df$BLUE / landsat_df$RED)

      # Soil Adjusted Vegetation Index (SAVI) = ((1.5)x NIR) - (red/0.5) + NIR + Red
      landsat_df$SAVI <- (1.5 * landsat_df$NIR) - (0.5 * landsat_df$RED) + landsat_df$NIR + landsat_df$RED

      # Vegetation Soil Salinity Index (VSSI) = (2 x green) - 5 x (red + NIR)
      landsat_df$VSSI <- (2 * landsat_df$GREEN) - 5 * (landsat_df$RED + landsat_df$NIR)

      # Rename some bands to match .rds use ("Blue_R", "Green_R", "Red_R", "NIR_R", "SWIR1_R", "SWIR2_R",)
      landsat_df$Blue_R <- landsat_df$BLUE
      landsat_df$Red_R <- landsat_df$RED
      landsat_df$Green_R <- landsat_df$GREEN
      landsat_df$NIR_R <- landsat_df$NIR
      landsat_df$SWIR1_R <- landsat_df$SWIR1
      landsat_df$SWIR2_R <- landsat_df$SWIR2


      landsat_df <- landsat_df[, c("ID", "Blue_R", "Green_R", "Red_R", "NIR_R", "SWIR1_R", "SWIR2_R",
                                   "NDVI", "NDWI", "NDSI1", "NDSI2", "SI1", "SI2", "SI3", "SI4", "SI5", "SAVI", "VSSI",
                                   "NBR", "NBG", "NBNIR", "NBSWIR1", "NBSWIR2", "NRSWIR1", "NRSWIR2", "NGSWIR1", "NGSWIR2",
                                   "NNIRSWIR1", "NNIRSWIR2")]

    
      # Remove NA values if any after processing and before applying the ensemble model
      length(landsat_df$ID)
      landsat_df <- landsat_df[complete.cases(landsat_df), ]
      length(landsat_df$ID)

      # Check df for any Infinite, or NaN values
      colSums(is.na(landsat_df))   # Check for NA values
      sum(is.infinite(as.matrix(landsat_df)))  # Check for Inf values
      sum(is.nan(as.matrix(landsat_df)))  # Check for NaN values

      # Find columns containing Inf values
      inf_cols <- sapply(landsat_df, function(x) any(is.infinite(x)))
      names(inf_cols[inf_cols == TRUE])

      # If any Red band is 0, SI5 will become Infinite. Omit those rows
      # Convert Inf values to NA first
      landsat_df <- as.data.frame(lapply(landsat_df, function(x) {
        if (is.numeric(x)) x[is.infinite(x)] <- NA  # Only apply to numeric columns
        return(x)
      }))
      # Omit rows with NA values
      landsat_df <- na.omit(landsat_df)
      # Check if Infinite values are addressed
      sum(is.infinite(as.matrix(landsat_df)))  # Should return 0


# 3. Apply salinity model
      predictions <- predict(MultiStratEnsemble, newdata = landsat_df)
      landsat_df$predicted_EC <- ifelse(predictions[, "X1"] > best_threshold, 1, 0)

      print(dim(landsat_df))   # Number of rows & columns
      print(colnames(landsat_df))  # Column names
      
      print(dim(predicted_df))   # Number of rows & columns
      print(colnames(predicted_df))  # Column names

      
# 4. Map predictions back to raster
      predicted_raster <- raster(landsat_image[[1]])  # Create an empty raster with the original landsat
      predicted_raster[aqua_mask_raster] <- 2 # Apply Aqua mask i.e. Replace predicted EC values with 2 where aquaculture mask is TRUE
      unique(values(predicted_raster))
      
      predicted_df <- as.data.frame(predicted_raster, xy = TRUE)
      predicted_df$ID <- 1:nrow(predicted_df)
      colnames(predicted_df) <- c("x", "y", "aqua", "ID")
      table(predicted_df$aqua)
      
      length(landsat_df$ID)
      length(predicted_df$ID)
      length(aqua_df$classification)
      
      predicted_df <- merge(predicted_df, landsat_df[, c("ID", "predicted_EC")], by = "ID", all.x = TRUE)
      # make a separate column such that it takes predicted_EC value if aqua is NA, else takes aqua value (i.e. 2) 
      predicted_df <- predicted_df %>%
        mutate(LandType = ifelse(is.na(aqua), predicted_EC, aqua))
      table(predicted_df$LandType)
      
      predicted_raster[] <- predicted_df$predicted_EC
     
      


# 5. Plot the combined raster (salinity + aquaculture)
        raster_plot <- ggplot(predicted_df, aes(x = x, y = y, fill = as.factor(LandType))) +
        geom_raster() +
        scale_fill_manual(values = c("white", "red", "blue"),
                          labels = c("Low Salinity", "High Salinity", "Aquaculture"),
                          name = "Land Type", 
                          na.value = "grey") +
        ggtitle(paste("Salinity & Aquaculture -", district, year)) +
        theme_minimal()

      plot(raster_plot)

      ggsave(paste0("outputs/", district, "_SalinityMap_", year, ".png"), raster_plot, width = 6, height = 4, dpi = 300)

      

# 8. Check before saving

      print("Summary of landsat_df:")
      print(summary(landsat_df))

      print("First few predictions:")
      print(head(predictions))

      print("Summary of predicted_raster:")
      print(summary(values(predicted_raster)))

      print("First few rows of predicted_df before writing to CSV:")
      print(head(predicted_df))
      
      print("First few rows of raster_df before writing to CSV:")
      print(head(raster_df))

# 9. Save files

      output_tiff <- paste0("outputs/", district, "_predicted_ECAqua_", year, ".tif")
      output_csv <- paste0("outputs/", district, "_predicted_ECAqua_", year, ".csv")

      writeRaster(predicted_raster, output_tiff, format = "GTiff", overwrite = TRUE)
      write.csv(predicted_df, output_csv, row.names = FALSE)

      print(paste("Processed:", district, year))


    }
  }
}


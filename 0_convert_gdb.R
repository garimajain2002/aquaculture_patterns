# Load required libraries
library(sf)
library(arcgisbinding)  # For reading .gdb files
library(dplyr)

convert_gdb_to_shp <- function(gdb_path, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize ArcGIS binding
  arc.check_product()
  
  # List all feature classes in the geodatabase
  gdb_layers <- arc.open(gdb_path) %>%
    arc.features() %>%
    names()
  
  # Print initial status
  cat(sprintf("Found %d layers in the geodatabase\n", length(gdb_layers)))
  
  # Process each layer
  for (layer_name in gdb_layers) {
    tryCatch({
      # Construct the full path to the feature class
      feature_path <- file.path(gdb_path, layer_name)
      
      # Read the feature class
      cat(sprintf("Processing %s...\n", layer_name))
      layer_data <- arc.open(feature_path) %>%
        arc.select() %>%
        arc.data2sf()
      
      # Clean layer name for file naming
      clean_name <- gsub("[^[:alnum:]]", "_", layer_name)
      output_path <- file.path(output_dir, paste0(clean_name, ".shp"))
      
      # Save as shapefile
      st_write(layer_data, output_path, append = FALSE)
      
      cat(sprintf("Successfully saved %s\n", output_path))
    }, error = function(e) {
      cat(sprintf("Error processing %s: %s\n", layer_name, e$message))
    })
  }
  
  cat("\nConversion complete!\n")
}
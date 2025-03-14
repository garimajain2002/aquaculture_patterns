# Convert district .gdb layers to shape files 

# Load required libraries
library(sf)

gdb_path <- "C:\\Users\\garim\\Desktop\\ArcGISPro\\Ch1\\Ch1_Patterns\\Ch1_Patterns.gdb"
output_dir <- "C:\\Users\\garim\\Desktop\\ArcGISPro\\Ch1\\Ch1_Patterns\\Ch1_shapefiles"

layers <- st_layers(gdb_path)
print(layers)


convert_gdb_to_shp <- function(gdb_path, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # List all layers in the geodatabase using sf
  gdb_layers <- st_layers(gdb_path)
  layer_names <- gdb_layers$name
  
  # Print initial status
  cat(sprintf("Found %d layers in the geodatabase\n", length(layer_names)))
  
  # Process each layer
  for (layer_name in layer_names) {
    tryCatch({
      # Read the layer using sf
      cat(sprintf("Processing %s...\n", layer_name))
      layer_data <- st_read(dsn = gdb_path, layer = layer_name)
      
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


convert_gdb_to_shp(gdb_path, output_dir)

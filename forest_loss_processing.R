library(raster)
library(seraphim)
library(rworldmap)
library(sf)
library(terra)

# Define the folder containing raster tiles from https://glad.earthengine.app/view/global-forest-change
tiles_folder <- "forest_loss/loss_year//"  # Replace path

# List all .tif files in the folder
tile_files <- list.files(tiles_folder, pattern = "\\.tif$", full.names = TRUE)

# Function to process a single raster
process_raster <- function(raster) {
  
  reclassified_raster<-calc( raster , function(x) { x[ x > 0 ] <- 1; return(x) } )
  reclassified_raster<-terra::aggregate(reclassified_raster,100)
  
  return(reclassified_raster)
}


# Function to process a single raster
library(terra)

process_raster_temporal <- function(r) {
  
  r <- terra::rast(r)
  
  # Define five-year periods
  periods <- list(
    "2001_2005" = c(1, 5),
    "2006_2010" = c(6, 10),
    "2011_2015" = c(11, 15),
    "2016_2020" = c(16, 20),
    "2021_2025" = c(21, 25)
  )
  
  # Create one aggregated raster per period
  period_rasters <- lapply(periods, function(year_range) {
    
    # Binary raster: 1 = forest loss during period, 0 = no loss
    loss <- ifel(
      r >= year_range[1] & r <= year_range[2],
      1,
      0
    )
    
    # Proportion of cells with loss in the period
    terra::aggregate(
      loss,
      fact = 100,
      fun = "mean",
      na.rm = TRUE
    )
  })
  
  # Combine the five periods into a multi-layer SpatRaster
  result <- rast(period_rasters)
  
  names(result) <- names(periods)
  
  return(result)
}

# Initialize an empty list to store processed rasters
processed_rasters <- list()

# Loop through each tile, process it, and store in the list
for (file in tile_files) {
  # Load the raster tile
  raster_tile <- raster(file)
  
  # Process the raster tile
  processed_raster <- process_raster_temporal(raster_tile)
  
  # Append to the list
  processed_rasters <- append(processed_rasters, list(processed_raster))
}

# Merge all processed rasters into a single raster
merged_raster <- do.call(mosaic, c(processed_rasters, list(fun = "mean")))

# Save the merged raster to disk
output_path <- "Africa_merged_forest_loss_medium_res_temporal.tif"  # Replace with your desired path
writeRaster(merged_raster, output_path, overwrite = TRUE)
merged_raster<-rast("Africa_merged_forest_loss_medium_res_temporal.tif")



library(geodata)
library(terra)
library(dplyr)

dir.create("data/ndvi", recursive = TRUE)
dir.create("data/landcover", recursive = TRUE)

#downloaded the coefficient of variation from https://www.earthenv.org/texture
evi_cv <- rast("/Users/monikam/Desktop/ebola_bat_models/Hypsignathus_monstrosus/data/evi/cv_01_05_1km_uint16.tif")

# Confirm it looks right over Africa
africa_ext <- ext(-18, 52, -35, 38)
evi_cv_africa <- crop(evi_cv, africa_ext)
plot(evi_cv_africa)
global(evi_cv_africa, fun = "range", na.rm = TRUE)

# Apply scaling
evi_cv_africa <- evi_cv_africa * 0.0001

# Check rescaled values
global(evi_cv_africa, fun = "range", na.rm = TRUE)
quantile(values(evi_cv_africa), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)

# Plot
plot(evi_cv_africa)
plot(evi_cv_africa, range = c(0, 0.2))

# Save the scaled and cropped EVI CV layer for Africa
writeRaster(evi_cv_africa, 
            "data/ndvi/evi_cv_africa.tif", 
            overwrite = TRUE)

#ESA land cover dataset
# This downloads fractional cover for each land cover class
#One of "trees", "grassland", "shrubs", "cropland", "built", "bare", 
# "snow", "water", "wetland", "mangroves", "moss"
lc <- landcover(var = "trees", path = "data/landcover")
plot(lc)

lc1 <- landcover(var = "grassland", path = "data/landcover")
plot(lc1)

lc2 <- landcover(var = "shrubs", path = "data/landcover")
plot(lc2)

lc3 <- landcover(var = "cropland", path = "data/landcover")
plot(lc3)

lc4 <- landcover(var = "built", path = "data/landcover")
plot(lc4)

lc5 <- landcover(var = "bare", path = "data/landcover")
plot(lc5)

lc6 <- landcover(var = "snow", path = "data/landcover")
plot(lc6)

lc7 <- landcover(var = "water", path = "data/landcover")
plot(lc7)

lc8 <- landcover(var = "wetland", path = "data/landcover")
plot(lc8)

lc9 <- landcover(var = "mangroves", path = "data/landcover")
plot(lc9)

lc10 <- landcover(var = "moss", path = "data/landcover")
plot(lc10)

# Stack all land cover layers together
lc_stack <- c(lc, lc1, lc2, lc3, lc4, lc5, lc6, lc7, lc8, lc9, lc10)

# Check the stack
lc_stack
nlyr(lc_stack)  # Should be 11


# Load all 19 bioclim layers
bioclim_files <- list.files("/Users/monikam/Desktop/ebola_bat_models/Hypsignathus_monstrosus/data/bioclim/climate/wc2.1_5m", 
                            pattern = ".tif$", 
                            full.names = TRUE)
bioclim_stack <- rast(bioclim_files)
plot(bioclim_stack)
# Check
bioclim_stack
nlyr(bioclim_stack)  # Should be 19
names(bioclim_stack)
# Reorder to numeric sequence
bioclim_stack <- bioclim_stack[[paste0("wc2.1_5m_bio_", 1:19)]]
names(bioclim_stack)  # Confirm order is now 1-19


###align all of the rasters
# Step 1 — define Africa extent
africa_ext <- ext(-18, 52, -35, 38)

# Step 2 — crop WorldClim to Africa
bioclim_africa <- crop(bioclim_stack, africa_ext)
plot(bioclim_africa)

# Step 3 — crop and resample land cover to match WorldClim
lc_africa <- crop(lc_stack, africa_ext)
lc_africa <- resample(lc_africa, bioclim_africa, method = "bilinear")
plot(lc_africa)

# Step 4 — resample EVI CV to match WorldClim
# it's already cropped to Africa but resolution needs to match
evi_cv_africa <- resample(evi_cv_africa, bioclim_africa, method = "bilinear")
plot(evi_cv_africa)

# Step 5 — stack everything together
predictors <- c(bioclim_africa, evi_cv_africa, lc_africa)

# Check final stack
predictors
nlyr(predictors)  # Should be 31 (19 bioclim + 1 evi + 11 landcover)
names(predictors)
names(predictors)[20] <- "evi_cv"

# Save the full aligned predictor stack
writeRaster(predictors, 
            "/Users/monikam/Desktop/ebola_bat_models/Hypsignathus_monstrosus/data/predictors_africa_full.tif",
            overwrite = TRUE)
predictors <- rast("data/predictors_africa_full.tif")

#multicollinearity

install.packages("usdm")
library(usdm)

#extracting bioclim variables from Koch et al 2020 paper: 
# Step 1 — subset predictors to paper-selected bioclim variables
# plus EVI CV and land cover layers
selected_bioclim <- c("wc2.1_5m_bio_1",   # Annual mean temperature
                      "wc2.1_5m_bio_4",   # Temperature seasonality
                      "wc2.1_5m_bio_5",   # Max temp warmest month
                      "wc2.1_5m_bio_6",   # Min temp coldest month
                      "wc2.1_5m_bio_12",  # Annual precipitation
                      "wc2.1_5m_bio_13",  # Precipitation wettest month
                      "wc2.1_5m_bio_14",  # Precipitation driest month
                      "wc2.1_5m_bio_15")  # Precipitation seasonality

# Get names of EVI and land cover layers
other_vars <- c("evi_cv", "trees", "grassland", "shrubs", "cropland",
                "built", "bare", "water", "wetland", 
                "mangroves", "moss")

# Subset the full predictor stack
predictors_sub <- predictors[[c(selected_bioclim, other_vars)]]
nlyr(predictors_sub)  # Should be 19 removed snow
names(predictors_sub)

# Step 2 — sample points for collinearity testing
set.seed(42)
pred_sample <- spatSample(predictors_sub, size = 50000,
                          method = "random",
                          na.rm = TRUE,
                          as.df = TRUE)

# Step 3 — correlation analysis on the full subset
v_cor <- vifcor(pred_sample, th = 0.7)
v_cor

# Step 4 — VIF analysis on remaining variables
v_vif <- vifstep(pred_sample, th = 10)
v_vif

# Final retained variable set - also removed mangroves
final_vars <- c("wc2.1_5m_bio_5", "wc2.1_5m_bio_6",
                "wc2.1_5m_bio_13", "wc2.1_5m_bio_14", "wc2.1_5m_bio_15",
                "evi_cv", "trees", "grassland", "shrubs", "cropland",
                "built", "water", "wetland")

# Subset predictor stack to final variables
predictors_final <- predictors_sub[[final_vars]]
nlyr(predictors_final)  # Should be 13
names(predictors_final)
plot(predictors_final)

# Save this final predictor stack
writeRaster(predictors_final,
            "/Users/monikam/Desktop/ebola_bat_models/Hypsignathus_monstrosus/data/predictors_africa_final.tif",
            overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))

library(randomForest)
library(maxnet)
library(dismo)
library(doParallel)
library(terra)
library(dplyr)
library(biomod2)
library(rnaturalearth)
library(sf)
library(tidyterra)
library(ggplot2)


# Reload predictors
predictors_final <- rast("data/predictors_africa_final.tif")

# Reformat data with disk strategy
# Samples PAs randomly but outside a minimum distance from presences
myBiomodData <- BIOMOD_FormatingData(
  resp.name     = "Hypsignathus.monstrosus",    #species name
  resp.var      = rep(1, nrow(hm_thinned)),
  resp.xy       = hm_thinned,
  expl.var      = predictors_final,
  PA.nb.rep     = 3,
  PA.nb.absences = 162,     #calculate based on number of occurrence points - presence-absence ratio
  PA.strategy   = "disk",
  PA.dist.min   = 150000,   # Minimum 150km from any presence
  PA.dist.max   = 1000000,  # Maximum 1000km from any presence
  na.rm         = TRUE,
  filter.raster = TRUE,
  seed.val      = 42
)

myBiomodData
plot(myBiomodData)


myBiomodModels <- BIOMOD_Modeling(
  bm.format    = myBiomodData,
  models       = c("GAM", "GBM", "RF", "MAXNET"),
  OPT.strategy = "bigboss",
  
  CV.strategy  = "kfold",
  CV.nb.rep    = 3,
  CV.k         = 5,
  
  metric.eval  = c("TSS", "AUCroc"),
  var.import   = 3,
  scale.models = FALSE,
  nb.cpu       = 8,   #edit for number of cores to be used
  seed.val     = 42
)


# Get evaluation scores
eval_scores <- get_evaluations(myBiomodModels)
head(eval_scores)

# Plot mean evaluation scores per algorithm
bm_PlotEvalMean(bm.out      = myBiomodModels,
                metric.eval = c("TSS", "AUCroc"))

# Boxplot to see variation across PA replicates and CV folds
bm_PlotEvalBoxplot(bm.out      = myBiomodModels,
                   group.by    = c("algo", "algo"))

# Get a cleaner summary of mean scores per algorithm
eval_scores %>%
  group_by(algo, metric.eval) %>%
  summarise(
    mean_calibration = mean(calibration, na.rm = TRUE),
    mean_validation  = mean(validation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(metric.eval, desc(mean_validation))

# project individual models
gc()
myBiomodProj <- BIOMOD_Projection(
  bm.mod              = myBiomodModels,
  new.env             = predictors_final,
  proj.name           = "Africa_current",
  models.chosen       = "all",
  metric.binary       = "TSS",
  build.clamping.mask = TRUE,
  nb.cpu              = 8
)

myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod               = myBiomodModels,
  models.chosen        = "all",
  em.by                = "all",
  em.algo              = c("EMwmean", "EMcv"),
  metric.select        = "TSS",
  metric.select.thresh = 0.8,
  metric.eval          = c("TSS","AUCroc"),
  var.import           = 3,
  seed.val             = 42
)

# Step 2 — ensemble forecasting from saved projections
gc()
myBiomodEMProj <- BIOMOD_EnsembleForecasting(
  bm.em         = myBiomodEM,
  bm.proj       = myBiomodProj,
  proj.name     = "Africa_current_EM",
  models.chosen = "all",
  metric.binary = "TSS",
  nb.cpu        = 8
)

# Base path
proj_path <- "Hypsignathus.monstrosus/proj_Africa_current_EM"

# Use the main ensemble tif 
hm_suit <- rast(file.path(proj_path,
                          "proj_Africa_current_EM_Hypsignathus.monstrosus_ensemble.tif"))

# Check
global(hm_suit, fun = "range", na.rm = TRUE)
nlyr(hm_suit)   # Check how many layers
names(hm_suit)  # Check layer names
plot(hm_suit)

# Separate the two layers
hm_wmean <- hm_suit[[1]]  # Habitat suitability
hm_cv    <- hm_suit[[2]]  # Uncertainty

# Rescale suitability to 0-1 (biomod2 uses 0-1000 scale)
hm_wmean <- hm_wmean / 1000

# Rename for clarity
names(hm_wmean) <- "habitat_suitability"
names(hm_cv)    <- "uncertainty"

# Verify
global(hm_wmean, fun = "range", na.rm = TRUE)
global(hm_cv,    fun = "range", na.rm = TRUE)

# Get Africa shapefile
africa <- ne_countries(continent = "Africa", returnclass = "sf")

# Plot habitat suitability
ggplot() +
  geom_spatraster(data = hm_wmean) +
  scale_fill_gradientn(
    colours   = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    na.value  = "grey90",
    name      = "Habitat\nsuitability",
    limits    = c(0, 1)
  ) +
  geom_sf(data = africa, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_point(data = hm_coords,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "red", size = 1, alpha = 0.7) +
  coord_sf(xlim = c(-18, 52), ylim = c(-35, 38)) +
  theme_minimal() +
  labs(title = "Hypsignathus monstrosus",
       subtitle = "EMwmean ensemble",
       x = "", y = "")

MyModels_var_import <- get_variables_importance(myBiomodModels)
MyModels_var_import
dimnames(MyModels_var_import)

MyModels_var_import %>% 
  #filter(algo == "GBM") %>% 
  ggplot()+
  geom_col(aes(x=algo, y=var.imp, fill= expl.var ), position = position_dodge())+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  facet_wrap(~run)

# Get mean validation scores per algorithm per species
get_evaluations(myBiomodModels) %>%
  group_by(algo, metric.eval) %>%
  summarise(mean_validation = round(mean(validation, na.rm = TRUE), 3),
            .groups = "drop") %>%
  arrange(metric.eval, desc(mean_validation))


# Create outputs folder 
dir.create("outputs/rasters", recursive = TRUE, showWarnings = FALSE)

# Export suitability raster
writeRaster(hm_wmean,
            "outputs/rasters/Hm_suitability_Africa.tif",
            overwrite = TRUE)



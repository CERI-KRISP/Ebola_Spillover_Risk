library(randomForest)
library(maxnet)
library(dismo)
library(doParallel)
library(terra)
library(dplyr)
library(biomod2)
library(rnaturalearth)
library(sf)
library(ggplot2)
library(tidyterra)

# ======tidyterra# ============================================================
# PREPARE INPUTS FOR BIOMOD_FormatingData
# ============================================================

# 1. RESPONSE VARIABLE — numeric vector: 1 = presence, 0 = absence
response <- ifelse(biomod_data$type == "Presence", 1, 0)

# 2. COORDINATES — numeric matrix, presences FIRST then absences
# biomod2 expects this exact column naming
coords <- as.matrix(biomod_data[, c("long", "lat")])
colnames(coords) <- c("X", "Y")

# 3. PREDICTOR RASTER — crop/mask your full env stack to M region
# (you already have env_M from earlier in the script — use that)
#load predictors
predictors_final <- rast("predictors_africa_final.tif")
#predictors_final <- env_M

# Quick checks before passing to biomod
cat("Response length:", length(response), "\n")        # Should be 267 (89 + 178)
cat("Coords rows:", nrow(coords), "\n")                # Should be 267
cat("Presences:", sum(response == 1), "\n")            # Should be 89
cat("Absences:", sum(response == 0), "\n")             # Should be 178
cat("Predictor layers:", nlyr(predictors_final), "\n") # Your env layer count

# ============================================================
# FORMAT DATA FOR BIOMOD2
# ============================================================

myBiomodData <- BIOMOD_FormatingData(
  resp.name     = "Ebola",
  resp.var      = response,
  resp.xy       = coords,
  expl.var      = predictors_final,
  na.rm         = TRUE,
  filter.raster = TRUE,
  seed.val      = 42
)

print(myBiomodData)


myBiomodModels <- BIOMOD_Modeling(
  bm.format    = myBiomodData,
  models       = c("GLM", "GBM", "RF", "MAXNET"),
  OPT.strategy = "bigboss",
  
  CV.strategy  = "kfold",
  CV.nb.rep    = 3,
  CV.k         = 5,
  
  metric.eval  = c("TSS", "AUCroc"),
  var.import   = 3,
  scale.models = FALSE,
  nb.cpu       = 8,
  seed.val     = 42
)

# Get evaluation scores
eval_scores <- get_evaluations(myBiomodModels)

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


#ensemble
gc()
myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod               = myBiomodModels,
  models.chosen        = "all",
  em.by                = "all",
  em.algo              = c("EMwmean", "EMcv"),
  metric.select        = "TSS",
  metric.select.thresh = 0.8,
  metric.eval          = c("TSS", "AUCroc"),
  var.import           = 10,
  seed.val             = 42
)

kept_models <- get_kept_models(myBiomodEM)
print(kept_models)
cat("Number of models kept:", length(kept_models), "\n")

gc()
myBiomodProj <- BIOMOD_Projection(
  bm.mod              = myBiomodModels,
  new.env             = predictors_final,
  proj.name           = "Africa_current",
  models.chosen       = "all",
  metric.binary       = "TSS",
  do.stack            = TRUE,
  build.clamping.mask = TRUE,
  nb.cpu              = 8
)

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
proj_path <- "Ebola/proj_Africa_current_EM"


suit <- rast(file.path(proj_path,
                       "proj_Africa_current_EM_Ebola_ensemble.tif"))

# Inspect
cat("Number of layers:", nlyr(suit), "\n")
cat("Layer names:\n")
print(names(suit))
global(suit, fun = "range", na.rm = TRUE)

# EMwmean = suitability, EMcv = coefficient of variation (uncertainty)
wmean_raw <- suit[[grep("EMwmean", names(suit))]]
cv_raw    <- suit[[grep("EMcv",    names(suit))]]

# Rescale from biomod2's 0-1000 scale to 0-1
wmean <- wmean_raw / 1000
cv    <- cv_raw    / 1000

names(wmean) <- "habitat_suitability"
names(cv)    <- "uncertainty"

# Verify ranges
cat("Suitability range:\n")
print(global(wmean, fun = "range", na.rm = TRUE))
cat("Uncertainty (CV) range:\n")
print(global(cv, fun = "range", na.rm = TRUE))

# Get Africa shapefile
africa <- ne_countries(continent = "Africa", returnclass = "sf")

# Plot habitat suitability
wmean_df <- as.data.frame(wmean, xy = TRUE) %>%
  rename(suitability = habitat_suitability) %>%
  filter(!is.na(suitability))

p_suit <- ggplot() +
  geom_tile(data = wmean_df, aes(x = x, y = y, fill = suitability)) +
  scale_fill_gradientn(
    colors = c('grey', '#fff3b0', '#e09f3e', '#9e2a2b', '#540b0e'),
    values = c(0, 0.3, 0.6, 0.9, 1.0),
    name      = "Habitat\nsuitability",
    na.value = "transparent",
    limits    = c(0, 1)
  ) +
  geom_sf(data = africa, fill = NA, colour = "grey40", linewidth = 0.3) +
  #geom_point(data = ebola_coords,
  #           aes(x = Longitude, y = Latitude),
  #           colour = "black", size = 2, alpha = 0.7) +
  #coord_sf(xlim = c(-18, 52), ylim = c(-35, 38)) +
  theme_void() +
  theme(legend.position = c(0.3, 0.4))+
labs(title = "Ebola zoonotic niche 2026",
    x = "", y = "")

p_suit

# --- Uncertainty map ---
cv_df <- as.data.frame(cv, xy = TRUE) %>%
  rename(uncertainty = uncertainty) %>%
  filter(!is.na(uncertainty))

ggplot() +
  geom_tile(data = cv_df, aes(x = x, y = y, fill = uncertainty)) +
  geom_sf(data = africa, fill = NA, color = "grey40", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#f7f7f7", "#f1a340", "#7b3294"),
    name    = "Uncertainty (CV)"
  ) +
  coord_sf(xlim = c(-25, 55), ylim = c(-35, 38)) +
  theme_minimal() +
  labs(title    = "Model Uncertainty",
       subtitle = "Coefficient of variation across ensemble members",
       x = "Longitude", y = "Latitude")


MyModels_var_import <- get_variables_importance(myBiomodModels)
MyModels_var_import
dimnames(MyModels_var_import)


MyModels_var_import %>% 
  filter(algo == "RF") %>% 
  ggplot()+
  geom_col(aes(x=algo, y=var.imp, fill= expl.var ), position = position_dodge())+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  facet_wrap(~run)


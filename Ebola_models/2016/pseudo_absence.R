# ============================================================
# Ebola SDM - Pseudo-absence Sampling Script
# ============================================================

library(terra)
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
install.packages("lwgeom")
library(lwgeom)
# Turn off S2 spherical geometry globally to avoid topology errors
sf_use_s2(FALSE)

# ============================================================
# 1. LOAD & PREPARE OCCURRENCE DATA
# ============================================================

dat <- read.csv("ebola_cleaned_occur.csv")

occ_sf <- st_as_sf(dat, coords = c("Longitude", "Latitude"), crs = 4326)

# ============================================================
# 2. LOAD ECOREGIONS AND CLIP TO AFRICA ONLY
# ============================================================

ecoregions <- st_read("/Users/monikam/Desktop/Ebola_2026_DRC/ebola niche model/modelling/Ecoregions2017/Ecoregions2017.shp")
ecoregions <- st_make_valid(ecoregions)

# Build Africa polygon for clipping
africa_union <- ne_countries(continent = "Africa", returnclass = "sf") %>%
  st_make_valid() %>%
  st_union() %>%
  st_make_valid()

# Clip ecoregions strictly to Africa — this is the key fix
# that prevents M region from including S. America / SE Asia
ecoregions_africa <- st_intersection(ecoregions, africa_union)
ecoregions_africa <- st_make_valid(ecoregions_africa)

# ============================================================
# 3. BUILD M REGION 
# ============================================================

# Join occurrences to African ecoregions only
occ_ecoreg    <- st_join(occ_sf, ecoregions_africa)
target_ecoreg <- unique(na.omit(occ_ecoreg$ECO_NAME))
cat("Target ecoregions:", length(target_ecoreg), "\n")

# Core M: ecoregions with presences
M_core <- ecoregions_africa %>%
  filter(ECO_NAME %in% target_ecoreg)

# Expanded M: all African ecoregions sharing same biomes as core
target_biomes <- unique(na.omit(M_core$BIOME_NAME))

M_region <- ecoregions_africa %>%
  filter(BIOME_NAME %in% target_biomes) %>%
  st_union() %>%
  st_make_valid()

# Verify M region stays within Africa
cat("M region bbox:\n")
print(st_bbox(M_region))
# xmin/xmax should be within -25 to 55
# ymin/ymax should be within -35 to 38

# Convert to terra for raster operations
M_vect <- vect(M_region)

# ============================================================
# 4. MASK ENVIRONMENTAL RASTER TO M REGION
# ============================================================

env <- rast("environmental_layers/elevation_PhyloCov_srtm_to_2025_scale5000m.tif")

# reproject M_vect to match raster CRS before crop/mask

M_vect_proj <- project(M_vect, crs(env))

env_M <- crop(env, M_vect_proj)
env_M <- mask(env_M, M_vect_proj)

cat("env_M extent:\n")
print(ext(env_M))

# ============================================================
# 5. SAMPLE PSEUDO-ABSENCE CANDIDATES FROM MASKED RASTER
# ============================================================

set.seed(42)
n_pres  <- nrow(dat)   # 53
n_pa    <- n_pres * 2  # 106 target pseudo-absences

pa_candidates <- spatSample(
  env_M,
  size   = n_pres * 30,  # oversample then filter
  method = "random",
  na.rm  = TRUE,         # IMPORTANT: skip NA cells (outside mask)
  xy     = TRUE,
  as.points = FALSE      # return as dataframe — more reliable
)

cat("Raw candidates sampled:", nrow(pa_candidates), "\n")

# Convert to sf using the raster's CRS, then reproject to WGS84
pa_sf <- st_as_sf(pa_candidates, coords = c("x", "y"), crs = crs(env_M))
pa_sf <- st_transform(pa_sf, 4326)

# Sanity check — all points should be in African coordinate range
coords_check <- st_coordinates(pa_sf)
cat("Candidate coordinate ranges:\n")
print(summary(coords_check))

# ============================================================
# 6. FILTER: REMOVE POINTS WITHIN 200 KM OF ANY PRESENCE
# ============================================================

# Ensure both layers are in WGS84
occ_sf <- st_transform(occ_sf, 4326)

# st_distance returns distances in metres when CRS is geographic
dist_matrix <- st_distance(pa_sf, occ_sf)
min_dist    <- apply(dist_matrix, 1, min)

pa_filtered <- pa_sf[min_dist > 200000, ]  # 200 km = 200,000 m
cat("Pseudo-absences after 100 km filter:", nrow(pa_filtered), "\n")

if (nrow(pa_filtered) < n_pa) {
  stop(paste("Only", nrow(pa_filtered), "pseudo-absences after filtering,",
             "but", n_pa, "required. Increase oversample multiplier."))
}

# ============================================================
# 7. FINAL RANDOM SAMPLE OF PSEUDO-ABSENCES
# ============================================================

set.seed(42)
pa_final       <- pa_filtered[sample(seq_len(nrow(pa_filtered)), n_pa), ]
pa_final$type  <- "Absence"

cat("Final pseudo-absences:", nrow(pa_final), "\n")

# ============================================================
# 8. COMBINE PRESENCE AND ABSENCE INTO ONE DATAFRAME
# ============================================================

pres_coords <- as.data.frame(st_coordinates(occ_sf))
pres_coords$type <- "Presence"

pa_coords   <- as.data.frame(st_coordinates(pa_final))
pa_coords$type   <- "Absence"

biomod_data <- rbind(pres_coords, pa_coords)
colnames(biomod_data) <- c("long", "lat", "type")

biomod_sf <- st_as_sf(biomod_data, coords = c("long", "lat"), crs = 4326)

# ============================================================
# 9. VISUALISE
# ============================================================

africa_map <- ne_countries(continent = "Africa", returnclass = "sf")

ggplot() +
  geom_sf(data = africa_map, fill = "grey95", color = "grey60") +
  geom_sf(
    data  = biomod_sf[biomod_sf$type == "Presence", ],
    color = "red", size = 2.5, shape = 16
  ) +
  geom_sf(
    data  = biomod_sf[biomod_sf$type == "Absence", ],
    color = "blue", size = 1.5, alpha = 0.5, shape = 16
  ) +
  coord_sf(xlim = c(-25, 55), ylim = c(-35, 38)) +
  theme_minimal() +
  labs(
    title    = "Ebola SDM: Presence & Pseudo-absence Points",
    subtitle = paste0("n presence = ", n_pres, "  |  n absence = ", n_pa),
    x = "Longitude", y = "Latitude"
  )

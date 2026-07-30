library(ggplot2)
library(terra)
library(terra)
library(sf)
library(rnaturalearth)
library(ggplot2)

# -----------------------------
# Loading area boundaries
# -----------------------------

#Whole Africa
Africa_map1 <- rnaturalearth::ne_countries(type = "countries", continent = "Africa",
                                           scale = "medium", returnclass = "sf")
# 1. Get all country names in the continent
africa_countries <- ne_countries(continent = "Africa", returnclass = "sf")

# 2. Extract the country names 
country_list <- unique(africa_countries$admin)

# 3. Get states (admin level 1) for those countries

Africa_map2 <- rnaturalearth::ne_states(country = country_list,
                                        returnclass = "sf")

Africa_vect <- vect(Africa_map1)


#Central Africa - selected countries
CentralAfrica_map1 <- rnaturalearth::ne_countries(type = "countries", country = c("Democratic Republic of the Congo", "Republic of the Congo", 
                                                                                  "Uganda", "Rwanda", "South Sudan", "Burundi", 
                                                                                  "United Republic of Tanzania", "Angola", "Central African Republic",
                                                                                  "Zambia", "Kenya","Gabon","Ethiopia"),
                                                  scale = "medium", returnclass = "sf")


CentralAfrica_map2 <- rnaturalearth::ne_states(country = c("Democratic Republic of the Congo", "Republic of the Congo", 
                                                           "Uganda", "Rwanda", "South Sudan", "Burundi", 
                                                           "United Republic of Tanzania", "Angola", "Central African Republic",
                                                           "Zambia", "Kenya","Gabon","Ethiopia"),
                                               returnclass = "sf")

drc_map2 <- rnaturalearth::ne_states(country = c("Democratic Republic of the Congo"),
                                     returnclass = "sf")

drc_map3<- ne_download(
  scale = 10, 
  type = 'admin_2_counties', 
  category = 'cultural', 
  returnclass = 'sf'
)


CentralAfrica_vect <- vect(CentralAfrica_map1)



# -----------------------------
# 2. Figure 2 - Updated vs published dataset Ebola maps
# -----------------------------

ebola_niche_2016_data<-rast("data/Ebola_zoo_niche_2016.tif")
ebola_niche_2016_data[ebola_niche_2016_data < 0.5] <- 0

ebola_niche_2026_data<-rast("data/Ebola_zoo_niche_2026.tif")
ebola_niche_2026_data[ebola_niche_2026_data < 0.5] <- 0

#Outbreak data

outbreaks<-read.csv('data/cdc_ebola_outbreaks_table_location_expanded.csv')
outbreaks<-subset(outbreaks, species!='Orthoebolavirus restonense')
# Define colours for Orthoebolavirus species
species_cols <- c(
  "Orthoebolavirus zairense"       = "#fdf0d5",
  "Orthoebolavirus sudanense"       = "#0072B2",
  "Orthoebolavirus bundibugyoense"  = "#009E73"
 # "Orthoebolavirus restonense"  = "#CC79A7"
)

point_cols <- adjustcolor(
  species_cols[outbreaks$species],
  alpha.f = 0.9
)

# Scale point sizes by number of cases
point_size <- scales::rescale(
  sqrt(outbreaks$reported_cases_total),
  to = c(1, 4)
)




par(mfrow = c(1, 3)) # Create a 2 x 2 plotting matrix
par(mar = c(0, 0, 0, 0)) # Set the margin on all sides to 2


plot(ebola_niche_2016_data, main = "Ebola ecological niche (2016 Data)",
     col = MetBrewer::met.brewer("Demuth", n = 256, type = "continuous", direction = -1),
     breaks = c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
 
plg = list(
  
  cex = 1.5,      # text size
  title.cex = 1.8, # title size
  size = 1.8      # width of the color bar
)
)
plot(Africa_map1,add=T,border='black',col=NA,lwd=0.1)
plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


#plot(ebola_niche_2026_data, main = "Ebola ecological niche",
#     col = MetBrewer::met.brewer("Demuth", n = 256, type = "continuous", direction = -1)
#)


# Plot ecological niche
plot(
  ebola_niche_2026_data,
  main = "Ebola ecological niche (2026 Data) & recorded outbreaks",
  col = MetBrewer::met.brewer(
    "Demuth",
    n = 256,
    type = "continuous",
    direction = -1
  ),
  breaks = c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
  plg = list(
    cex = 1.5,
    title.cex = 1.8,
    size = 1.8
  )
)

# Add administrative boundaries
plot(
  Africa_map1,
  add = TRUE,
  border = "black",
  col = NA,
  lwd = 0.1
)

plot(
  CentralAfrica_map2,
  add = TRUE,
  border = "grey50",
  col = NA,
  lwd = 0.1
)

plot(
  CentralAfrica_map1,
  add = TRUE,
  border = "grey90",
  col = NA
)

# Add outbreak locations
points(
  outbreaks$longitude,
  outbreaks$latitude,
  pch = 1,
  col = point_cols,
 # bg = NA,
  lwd = 1,
  cex = point_size
)

legend(
  x = -25,
  y = -5,
  legend = names(species_cols),
  col = adjustcolor(species_cols, alpha.f = 1),
  pch = 21,
  #col = "black",
  pt.cex = 1.5,
  title = "Orthoebolavirus species",
  bty = "n"
)

legend(
  x = -25,
  y = -20,
  legend = c(10, 100, 1000),
  pt.cex = scales::rescale(
    sqrt(c(10, 100, 1000)),
    to = c(1, 4)
  ),
  pch = 21,
  pt.bg = "white",
  col = "black",
  title = "Number of cases",
  bty = "n"
)


plot(ebola_niche_2026_data-ebola_niche_2016_data, main = "Ebola ecological niche - difference",
     col = MetBrewer::met.brewer("Troy", n = 256, type = "continuous", direction = -1),
     plg = list(
       
       cex = 1.5,      # text size
       title.cex = 1.8, # title size
       size = 0.3      # width of the color bar
     )
)
plot(Africa_map1,add=T,border='black',col=NA,lwd=0.1)
plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

#quartz.save(file = "FIG2_updatedebolaniche_comparison_v3.pdf", type = "pdf")


# -----------------------------
# Loading covariates rasters for spillover model
# -----------------------------
built <- rast("data/GHS_BUILT_S_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif")    # built-up surface
mining <- rast("data/distance_to_nearest_mining_site_africa_0p05deg_osm_plus_non_osm.tif")
conflict_exposure<-rast("data/ConflictEvents_MeanAnnual_2000_2026.tif")
conflict_exposure2<-rast("data/Conflict_mean_annual_events_2000-2010.tif")
conflict_exposure3<-rast("data/Conflict_mean_annual_events_2011-2021.tif")
conflict_exposure4<-rast("data/Conflict_mean_annual_events_2022-2026.tif")

conflict_exposure1<-rast("data/ConflictExposure_MeanAnnual_2014_2026.tif")
conflict_exposure1a<-rast("data/ConflictExposure_2014_2021_MeanAnnual.tif")
conflict_exposure1b<-rast("data/ConflictExposure_2022_2026_MeanAnnual.tif")

ebola_niche<-rast("data/Ebola_zoo_niche_2026.tif")
ebola_niche_2016<-rast("data/Ebola_zoo_niche_2016.tif")

forest_loss<-rast("data/Africa_merged_forest_loss_medium_res.tif")
forest_loss_temporal<-rast("data/Africa_merged_forest_loss_medium_res_temporal.tif")

bushmeat<-rast("data/bushmeat_activity_Binomial_predictions.tif")

# -----------------------------
# 5. Crop covariates to bounding box 
# -----------------------------
built <- mask(crop(built, CentralAfrica_vect),CentralAfrica_vect)    # built-up surface
mining <- mask(crop(mining, CentralAfrica_vect),CentralAfrica_vect)
conflict_exposure<-mask(crop(conflict_exposure, CentralAfrica_vect),CentralAfrica_vect)
bushmeat<-mask(crop(bushmeat, CentralAfrica_vect),CentralAfrica_vect)
bushmeat <- project(bushmeat, crs(CentralAfrica_vect))
ebola_niche<-mask(crop(ebola_niche, CentralAfrica_vect),CentralAfrica_vect)
ebola_niche_2016<-mask(crop(ebola_niche_2016, CentralAfrica_vect),CentralAfrica_vect)
forest_loss<-mask(crop(forest_loss, CentralAfrica_vect),CentralAfrica_vect)
forest_loss_temporal<-mask(crop(forest_loss_temporal, CentralAfrica_vect),CentralAfrica_vect)

#turn distance to mining sites into mining proximity (0-1)
k <- 10
mining_risk <- exp(-mining / k)

forest_loss <- project(forest_loss, ebola_niche)
forest_loss_temporal <- project(forest_loss_temporal, ebola_niche)
conflict_exposure <- project(conflict_exposure, ebola_niche)
built <- project(built, ebola_niche)
mining_risk <- project(mining_risk, ebola_niche)
bushmeat <- project(bushmeat, ebola_niche)


# -----------------------------
# Compute spillover risk and evaluate performance of scenarios
# -----------------------------

#normalizing built and conflict layers - others already normalised

built<-log1p(built)

built_min <- global(built, "min", na.rm = TRUE)[1,1]
built_max <- global(built, "max", na.rm = TRUE)[1,1]

built <- (built - built_min) / (built_max - built_min)


conflict_min <- global(conflict_exposure, "min", na.rm = TRUE)[1,1]
conflict_max <- global(conflict_exposure, "max", na.rm = TRUE)[1,1]

conflict_exposure <- (conflict_exposure - conflict_min) / (conflict_max - conflict_min)

#spillover risk scenarios
spillover_risk_0 <- 
  #log1p(pop) * 
  built * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche


spillover_risk_1 <- 
  #log1p(pop) * 
  built * 
  mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  forest_loss * 
  bushmeat *
  ebola_niche

spillover_risk_2 <- 
  #log1p(pop) * 
  built * 
  mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche

spillover_risk_3 <- 
  #log1p(pop) * 
  built * 
  #mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche

spillover_risk_4 <- 
  #log1p(pop) * 
  built * 
  mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche

spillover_risk_5 <- 
  #log1p(pop) * 
  built * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  forest_loss * 
  bushmeat *
  ebola_niche

spillover_risk_6 <- 
  #log1p(pop) * 
  built * 
  #mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  forest_loss * 
  bushmeat *
  ebola_niche

spillover_risk_7 <- 
  #log1p(pop) * 
  built * 
  mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  forest_loss * 
  bushmeat *
  ebola_niche



spillover_risk_8 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  bushmeat *
  ebola_niche

spillover_risk_9 <- 
  #log1p(pop) * 
  #log1p(built) * 
  mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche

spillover_risk_10 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche



spillover_risk_11 <- 
  #log1p(pop) * 
  built * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche



spillover_risk_12 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche





spillover_risk_13 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche_2016


spillover_risk_14 <- 
  #log1p(pop) * 
  #log1p(built) * 
  mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche_2016

spillover_risk_15 <- 
  #log1p(pop) * 
  built * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche_2016


spillover_risk_16 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  bushmeat *
  ebola_niche_2016


spillover_risk_17 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche_2016



spillover_risk_18 <- 
  #log1p(pop) * 
  built * 
  mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche


spillover_risk_19 <- 
  #log1p(pop) * 
  built * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  bushmeat *
  ebola_niche


spillover_risk_20 <- 
  #log1p(pop) * 
  built * 
  #mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche


spillover_risk_21 <- 
  #log1p(pop) * 
  #log1p(built) * 
  mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche


spillover_risk_22 <- 
  #log1p(pop) * 
  #log1p(built) * 
  mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche


spillover_risk_23 <- 
  #log1p(pop) * 
  #log1p(built) * 
  mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  bushmeat *
  ebola_niche


spillover_risk_24 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  forest_loss * 
  bushmeat *
  ebola_niche


spillover_risk_25 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche

# normalize to 0–1 (probability-like)

risk_min0 <- global(spillover_risk_0, "min", na.rm = TRUE)[1,1]
risk_max0 <- global(spillover_risk_0, "max", na.rm = TRUE)[1,1]

risk_prob0 <- (spillover_risk_0 - risk_min0) / (risk_max0 - risk_min0)
names(risk_prob0)<-"risk_prob_base_built_forestloss"


risk_min1 <- global(spillover_risk_1, "min", na.rm = TRUE)[1,1]
risk_max1 <- global(spillover_risk_1, "max", na.rm = TRUE)[1,1]

risk_prob1 <- (spillover_risk_1 - risk_min1) / (risk_max1 - risk_min1)
names(risk_prob1)<-"risk_prob_mining_conflict_bushmeat"

risk_min2 <- global(spillover_risk_2, "min", na.rm = TRUE)[1,1]
risk_max2 <- global(spillover_risk_2, "max", na.rm = TRUE)[1,1]

risk_prob2 <- (spillover_risk_2 - risk_min2) / (risk_max2 - risk_min2)
names(risk_prob2)<-"risk_prob2_mining_conflict"



risk_min3 <- global(spillover_risk_3, "min", na.rm = TRUE)[1,1]
risk_max3 <- global(spillover_risk_3, "max", na.rm = TRUE)[1,1]

risk_prob3 <- (spillover_risk_3 - risk_min3) / (risk_max3 - risk_min3)
names(risk_prob3)<-"risk_prob3_conflict"


risk_min4 <- global(spillover_risk_4, "min", na.rm = TRUE)[1,1]
risk_max4 <- global(spillover_risk_4, "max", na.rm = TRUE)[1,1]

risk_prob4 <- (spillover_risk_4 - risk_min4) / (risk_max4 - risk_min4)
names(risk_prob4)<-"risk_prob4_mining"


risk_min5 <- global(spillover_risk_5, "min", na.rm = TRUE)[1,1]
risk_max5 <- global(spillover_risk_5, "max", na.rm = TRUE)[1,1]

risk_prob5 <- (spillover_risk_5 - risk_min5) / (risk_max5 - risk_min5)

names(risk_prob5)<-"risk_prob5_bushmeat"

risk_min6 <- global(spillover_risk_6, "min", na.rm = TRUE)[1,1]
risk_max6 <- global(spillover_risk_6, "max", na.rm = TRUE)[1,1]

risk_prob6 <- (spillover_risk_6 - risk_min6) / (risk_max6 - risk_min6)
names(risk_prob6)<-"risk_prob6_bushmeat_conflict"

risk_min7 <- global(spillover_risk_7, "min", na.rm = TRUE)[1,1]
risk_max7 <- global(spillover_risk_7, "max", na.rm = TRUE)[1,1]

risk_prob7 <- (spillover_risk_7 - risk_min7) / (risk_max7 - risk_min7)
names(risk_prob7)<-"risk_prob7_bushmeat_mining"

risk_min8 <- global(spillover_risk_8, "min", na.rm = TRUE)[1,1]
risk_max8 <- global(spillover_risk_8, "max", na.rm = TRUE)[1,1]

risk_prob8 <- (spillover_risk_8 - risk_min8) / (risk_max8 - risk_min8)
names(risk_prob8)<-"risk_prob8_habitat_bushmeat"


risk_min9 <- global(spillover_risk_9, "min", na.rm = TRUE)[1,1]
risk_max9 <- global(spillover_risk_9, "max", na.rm = TRUE)[1,1]

risk_prob9 <- (spillover_risk_9 - risk_min9) / (risk_max9 - risk_min9)
names(risk_prob9)<-"risk_prob9_habitat_mining"



risk_min10 <- global(spillover_risk_10, "min", na.rm = TRUE)[1,1]
risk_max10 <- global(spillover_risk_10, "max", na.rm = TRUE)[1,1]

risk_prob10 <- (spillover_risk_10 - risk_min10) / (risk_max10 - risk_min10)
names(risk_prob10)<-"risk_prob10_habitat_conflict"



risk_min11 <- global(spillover_risk_11, "min", na.rm = TRUE)[1,1]
risk_max11 <- global(spillover_risk_11, "max", na.rm = TRUE)[1,1]

risk_prob11 <- (spillover_risk_11 - risk_min11) / (risk_max11 - risk_min11)
names(risk_prob11)<-"risk_prob11_habitat_built"


risk_min12 <- global(spillover_risk_12, "min", na.rm = TRUE)[1,1]
risk_max12 <- global(spillover_risk_12, "max", na.rm = TRUE)[1,1]

risk_prob12 <- (spillover_risk_12 - risk_min12) / (risk_max12 - risk_min12)
names(risk_prob12)<-"risk_prob12_habitat_forestloss"



risk_min13 <- global(spillover_risk_13, "min", na.rm = TRUE)[1,1]
risk_max13 <- global(spillover_risk_13, "max", na.rm = TRUE)[1,1]

risk_prob13 <- (spillover_risk_13 - risk_min13) / (risk_max13 - risk_min13)
names(risk_prob13)<-"risk_prob13_habitat2016_forestloss"



risk_min14 <- global(spillover_risk_14, "min", na.rm = TRUE)[1,1]
risk_max14 <- global(spillover_risk_14, "max", na.rm = TRUE)[1,1]

risk_prob14 <- (spillover_risk_14 - risk_min14) / (risk_max14 - risk_min14)
names(risk_prob14)<-"risk_prob14_habitat2016_mining"



risk_min15 <- global(spillover_risk_15, "min", na.rm = TRUE)[1,1]
risk_max15 <- global(spillover_risk_15, "max", na.rm = TRUE)[1,1]

risk_prob15 <- (spillover_risk_15 - risk_min15) / (risk_max15 - risk_min15)
names(risk_prob15)<-"risk_prob15_habitat2016_built"




risk_min16 <- global(spillover_risk_16, "min", na.rm = TRUE)[1,1]
risk_max16 <- global(spillover_risk_16, "max", na.rm = TRUE)[1,1]

risk_prob16 <- (spillover_risk_16 - risk_min16) / (risk_max16 - risk_min16)
names(risk_prob16)<-"risk_prob16_habitat2016_bushmeat"


risk_min17 <- global(spillover_risk_17, "min", na.rm = TRUE)[1,1]
risk_max17 <- global(spillover_risk_17, "max", na.rm = TRUE)[1,1]

risk_prob17 <- (spillover_risk_17 - risk_min17) / (risk_max17 - risk_min17)
names(risk_prob17)<-"risk_prob17_habitat2017_conflict"



risk_min18 <- global(spillover_risk_18, "min", na.rm = TRUE)[1,1]
risk_max18 <- global(spillover_risk_18, "max", na.rm = TRUE)[1,1]

risk_prob18 <- (spillover_risk_18 - risk_min18) / (risk_max18 - risk_min18)
names(risk_prob18)<-"risk_prob18_habitat_built_mining"


risk_min19 <- global(spillover_risk_19, "min", na.rm = TRUE)[1,1]
risk_max19 <- global(spillover_risk_19, "max", na.rm = TRUE)[1,1]

risk_prob19 <- (spillover_risk_19 - risk_min19) / (risk_max19 - risk_min19)
names(risk_prob19)<-"risk_prob19_habitat_built_bushmeat"


risk_min20 <- global(spillover_risk_20, "min", na.rm = TRUE)[1,1]
risk_max20 <- global(spillover_risk_20, "max", na.rm = TRUE)[1,1]

risk_prob20 <- (spillover_risk_20 - risk_min20) / (risk_max20 - risk_min20)
names(risk_prob20)<-"risk_prob20_habitat_built_conflict"


risk_min21 <- global(spillover_risk_21, "min", na.rm = TRUE)[1,1]
risk_max21 <- global(spillover_risk_21, "max", na.rm = TRUE)[1,1]

risk_prob21 <- (spillover_risk_21 - risk_min21) / (risk_max21 - risk_min21)
names(risk_prob21)<-"risk_prob21_habitat_mining_forestloss"



risk_min22 <- global(spillover_risk_22, "min", na.rm = TRUE)[1,1]
risk_max22 <- global(spillover_risk_22, "max", na.rm = TRUE)[1,1]

risk_prob22 <- (spillover_risk_22 - risk_min22) / (risk_max22 - risk_min22)
names(risk_prob22)<-"risk_prob22_habitat_mining_conflict"



risk_min23 <- global(spillover_risk_23, "min", na.rm = TRUE)[1,1]
risk_max23 <- global(spillover_risk_23, "max", na.rm = TRUE)[1,1]

risk_prob23 <- (spillover_risk_23 - risk_min23) / (risk_max23 - risk_min23)
names(risk_prob23)<-"risk_prob23_habitat_mining_bushmeat"


risk_min24 <- global(spillover_risk_24, "min", na.rm = TRUE)[1,1]
risk_max24 <- global(spillover_risk_24, "max", na.rm = TRUE)[1,1]

risk_prob24 <- (spillover_risk_24 - risk_min24) / (risk_max24- risk_min24)
names(risk_prob24)<-"risk_prob24_habitat_forestloss_bushmeant"


risk_min25 <- global(spillover_risk_25, "min", na.rm = TRUE)[1,1]
risk_max25 <- global(spillover_risk_25, "max", na.rm = TRUE)[1,1]

risk_prob25 <- (spillover_risk_25 - risk_min25) / (risk_max25- risk_min25)
names(risk_prob25)<-"risk_prob25_habitat_forestloss_conflict"







names(ebola_niche_2016)<-"Habitat_suitability_2016"

risk_stack0 <- c(ebola_niche,ebola_niche_2016,risk_prob0,risk_prob1,risk_prob2,risk_prob3,risk_prob4,
                risk_prob5,risk_prob6,risk_prob7,
                risk_prob8, risk_prob9, risk_prob10)


risk_stack1 <- c(ebola_niche,ebola_niche_2016,
                 risk_prob8,
                 risk_prob9, risk_prob10,
                 risk_prob11, risk_prob12,
                 risk_prob18,
                 risk_prob19, risk_prob20,
                 risk_prob21, risk_prob22,
                 risk_prob23,
                 risk_prob24, risk_prob25)



risk_stack <- c(ebola_niche,risk_prob8,
                risk_prob9, risk_prob10,
                risk_prob11, risk_prob12)


risk_stack2 <- c(ebola_niche_2016,risk_prob13,
                risk_prob14, risk_prob15,
                risk_prob16, risk_prob17)

template<-ebola_niche

library(readxl)

library(terra)
library(ecospat)
library(dplyr)
library(purrr)

calc_boyce <- function(r, occ){
  
  occ_vals <- terra::extract(r, occ)[,2]
  occ_vals <- occ_vals[!is.na(occ_vals)]
  
  bg_vals <- values(r)
  bg_vals <- bg_vals[!is.na(bg_vals)]
  
  if(length(occ_vals) < 5)
    return(NA)
  
  if(length(unique(occ_vals)) < 2)
    return(NA)
  
  boyce <- try(
    ecospat::ecospat.boyce(
      fit = bg_vals,
      obs = occ_vals,
      PEplot = FALSE
    ),
    silent = TRUE
  )
  
  if(inherits(boyce, "try-error"))
    return(NA)
  
  if("Spearman.cor" %in% names(boyce))
    return(boyce$Spearman.cor)
  
  if("cor" %in% names(boyce))
    return(boyce$cor)
  
  return(NA)
}

calc_point_corr <- function(r, occ){
  
  occ_vals <- terra::extract(
    r,
    occ
  )[,2]
  
  occ_vals <- occ_vals[!is.na(occ_vals)]
  
  bg_vals <- sample(
    values(r),
    min(10000, sum(!is.na(values(r))))
  )
  
  bg_vals <- bg_vals[!is.na(bg_vals)]
  
  wilcox.test(
    occ_vals,
    bg_vals
  )$p.value
}

calc_enrichment <- function(r, occ){
  
  occ_vals <- terra::extract(
    r,
    occ
  )[,2]
  
  occ_vals <- occ_vals[!is.na(occ_vals)]
  
  mean(occ_vals) /
    global(r, "mean", na.rm=TRUE)[1,1]
}


spillover_events<-read_excel("data/elife_advances_human_index_merged_with_new_info.xls")
spillover_events$Year.Start<-spillover_events$`Year Start`
#spillover_events<-read.csv("elife_advances_human_index_merged_with_new_info_CM.csv")
spillover_events_before2022<-subset(spillover_events, Year.Start<2022)
spillover_events_after2022<-subset(spillover_events, Year.Start>=2022)
spillover_events_before2010<-subset(spillover_events, Year.Start<2010)
spillover_events_after2015<-subset(spillover_events, Year.Start>=2015)
spillover_events_2010to2021<-subset(subset(spillover_events, Year.Start>=2010),Year.Start<2022)

spillover_events$Lon<-as.numeric(spillover_events$Long)
spillover_events$Lat<-as.numeric(spillover_events$Lat)
spillover_events_before2022$Lon<-as.numeric(spillover_events_before2022$Long)
spillover_events_before2022$Lat<-as.numeric(spillover_events_before2022$Lat)
spillover_events_after2022$Lon<-as.numeric(spillover_events_after2022$Long)
spillover_events_after2022$Lat<-as.numeric(spillover_events_after2022$Lat)
spillover_events_before2010$Lon<-as.numeric(spillover_events_before2010$Long)
spillover_events_before2010$Lat<-as.numeric(spillover_events_before2010$Lat)
spillover_events_after2015$Lon<-as.numeric(spillover_events_after2015$Long)
spillover_events_after2015$Lat<-as.numeric(spillover_events_after2015$Lat)
spillover_events_2010to2021$Lon<-as.numeric(spillover_events_2010to2021$Long)
spillover_events_2010to2021$Lat<-as.numeric(spillover_events_2010to2021$Lat)

spillover_events_cases <- vect(
  spillover_events,
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326"
)
spillover_events_cases_before2022 <- vect(
  spillover_events_before2022,
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326"
)

spillover_events_cases_after2022 <- vect(
  spillover_events_after2022,
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326"
)

spillover_events_cases_before2010 <- vect(
  spillover_events_before2010,
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326"
)

spillover_events_cases_after2015 <- vect(
  spillover_events_after2015,
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326"
)

spillover_events_cases_2010to2021 <- vect(
  spillover_events_2010to2021,
  geom = c("Lon", "Lat"),
  crs = "EPSG:4326"
)



spillover_events_cases <- project(spillover_events_cases, template)
spillover_events_cases_before2022 <- project(spillover_events_cases_before2022, template)
spillover_events_cases_after2022 <- project(spillover_events_cases_after2022, template)
spillover_events_cases_before2010 <- project(spillover_events_cases_before2010, template)
spillover_events_cases_after2015 <- project(spillover_events_cases_after2015, template)
spillover_events_cases_2010to2021 <- project(spillover_events_cases_2010to2021, template)



evaluationbefore2010 <- map_dfr(
  1:nlyr(risk_stack1),
  function(i){
    
    r <- risk_stack1[[i]]
    
    tibble(
      Scenario = names(r),
      
      Boyce =
        calc_boyce(r, spillover_events_cases_before2010),
      
      Enrichment =
        calc_enrichment(r, spillover_events_cases_before2010),
      
    )
  }
)
evaluationbefore2010$Period<-"1[<2010]"


evaluationbefore2010 %>%
  arrange(
    desc(Boyce),
    desc(Enrichment)
  )


evaluationafter2022 <- map_dfr(
  1:nlyr(risk_stack1),
  function(i){
    
    r <- risk_stack1[[i]]
    
    tibble(
      Scenario = names(r),
      
      Boyce =
        calc_boyce(r, spillover_events_cases_after2022),
      
      Enrichment =
        calc_enrichment(r, spillover_events_cases_after2022),
      
    )
  }
)
evaluationafter2022$Period<-"3[>=2022]"


evaluationafter2022 %>%
  arrange(
    desc(Boyce),
    desc(Enrichment)
  )

evaluation2010to2021 <- map_dfr(
  1:nlyr(risk_stack1),
  function(i){
    
    r <- risk_stack1[[i]]
    
    tibble(
      Scenario = names(r),
      
      Boyce =
        calc_boyce(r, spillover_events_cases_2010to2021),
      
      Enrichment =
        calc_enrichment(r, spillover_events_cases_2010to2021),
      
    )
  }
)
evaluation2010to2021$Period<-"2[2010-2021]"


evaluation2010to2021 %>%
  arrange(
    desc(Boyce),
    desc(Enrichment)
  )



all_evalutation<-rbind(evaluationbefore2010,evaluation2010to2021,evaluationafter2022)


aesthetic15 <- c(
  "#264653",  # deep teal
  "#2A9D8F",  # teal
  "#8AB17D",  # sage
  "#A7C957",  # olive green
  "#E9C46A",  # sand
  "#F4A261",  # peach
  "#E76F51",  # terracotta
  "#D1495B",  # coral red
  "#9C6644",  # brown
  "#6D597A",  # mauve
  "#B56576",  # dusty rose
  "#EAAC8B",  # blush
  "#457B9D",  # steel blue
  "#5F0F40",  # burgundy
  "#7F5539"   # chestnut
)

p_boyce<-ggplot(all_evalutation,aes(Period,Boyce,colour=Scenario,group=Scenario))+
  theme_bw()+
  geom_line(size=1)+
  geom_point(size=2)+
  theme(legend.position = "bottom", legend.direction="vertical",legend.text = element_text(size=8,margin = margin(l = 5)),legend.title = element_text(size=10,margin = margin(l = 5,w=5)),   legend.key.spacing.x = unit(1, "pt"),
        legend.key.spacing.y = unit(1, "pt"),
        theme(legend.key.height = unit(0.4, "null")),
        legend.key.size = unit(0.4, "lines")
        
  )+
  
  guides(color = guide_legend(ncol = 3))+
  scale_color_manual(values = aesthetic15, 
                     name='Spillover risk scenario',
                     labels=c('Habitat Suitability (2026 Data)','Habitat Suitability (2016 Data)','Habitat x Conflict', 'Habitat x Built', 'Habitat x Forest Loss', 'Habitat x Built x Mining', 'Habitat x Buit x Bushmeat', 'Habitat x Buit x Conflict', 'Habitat x Mining x Forest Loss',
                              'Habitat x Mining x Conflict', 'Habitat x Mining x Bushmeat', 'Habitat x Forest Loss x Bushmeat', 
                              'Habitat x Forest Loss x Conflict', 'Habitat x Bushmeat', 'Habitat x mining')
 )
p_enrichment<-ggplot(all_evalutation,aes(Period,Enrichment,colour=Scenario,group=Scenario))+
  theme_bw()+
  geom_line(size=1)+
  geom_point(size=2)+
  theme(legend.position = 'none')+
  scale_color_manual(values = aesthetic15)
cowplot::plot_grid(p_boyce,p_enrichment,ncol=1, rel_heights = c(0.55,0.45))

#quartz.save(file = "spillover_risk_boyce_enrichment_final_habitat_plus_1_and_2_factors.pdf", type = "pdf")






evaluationOverall <- map_dfr(
  1:nlyr(risk_stack),
  function(i){
    
    r <- risk_stack[[i]]
    
    tibble(
      Scenario = names(r),
      
      Boyce =
        calc_boyce(r, spillover_events_cases),
      
      Enrichment =
        calc_enrichment(r, spillover_events_cases),
      
    )
  }
)


evaluationOverall %>%
  arrange(
    desc(Boyce),
    desc(Enrichment)
  )



evaluationOverall2016 <- map_dfr(
  1:nlyr(risk_stack2),
  function(i){
    
    r <- risk_stack2[[i]]
    
    tibble(
      Scenario = names(r),
      
      Boyce =
        calc_boyce(r, spillover_events_cases),
      
      Enrichment =
        calc_enrichment(r, spillover_events_cases),
      
    )
  }
)


evaluationOverall2016 %>%
  arrange(
    desc(Boyce),
    desc(Enrichment)
  )

evaluationOverall2016$SuitabilityVersion<-"2016"
evaluationOverall2016[evaluationOverall2016 == "risk_prob15_habitat2016_built"] <- "Habitat x Built"
evaluationOverall2016[evaluationOverall2016 == "risk_prob14_habitat2016_mining"] <- "Habitat x Mining"
evaluationOverall2016[evaluationOverall2016 == "risk_prob16_habitat2016_bushmeat"] <- "Habitat x Bushmeat"
evaluationOverall2016[evaluationOverall2016 == "risk_prob17_habitat2017_conflict"] <- "Habitat x Conflict"
evaluationOverall2016[evaluationOverall2016 == "Habitat_suitability_2016"] <- "Habitat"
evaluationOverall2016[evaluationOverall2016 == "risk_prob13_habitat2016_forestloss"] <- "Habitat x Forest Loss"

evaluationOverall$SuitabilityVersion<-"2026"
evaluationOverall[evaluationOverall == "risk_prob11_habitat_built"] <- "Habitat x Built"
evaluationOverall[evaluationOverall == "risk_prob9_habitat_mining"] <- "Habitat x Mining"
evaluationOverall[evaluationOverall == "risk_prob8_habitat_bushmeat"] <- "Habitat x Bushmeat"
evaluationOverall[evaluationOverall == "risk_prob10_habitat_conflict"] <- "Habitat x Conflict"
evaluationOverall[evaluationOverall == "habitat_suitability"] <- "Habitat"
evaluationOverall[evaluationOverall == "risk_prob12_habitat_forestloss"] <- "Habitat x Forest Loss"

all_evalutationOverall<-rbind(evaluationOverall,evaluationOverall2016)

ggplot(all_evalutationOverall, aes(x = Boyce, y = reorder(Scenario,Boyce), group=Scenario)) +
  theme_bw()+
  geom_line() +
  geom_point(aes(color = SuitabilityVersion, size = Enrichment),alpha=0.5) +
  theme(legend.position = "bottom")+
  scale_colour_manual(values=c('goldenrod3','steelblue4'))+
  xlab('Predictive power of spillover occurrences (Boyce index)')+
  ylab('Spillover risk scenario')


#quartz.save(file = "spillover_risk_boyce_2016vs2026data.pdf", type = "pdf")





ggplot(
  all_evalutationOverall,
  aes(
    reorder(Scenario, Boyce),
    Boyce
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "",
    y = "Boyce Index"
  )


threshold <- quantile(
  values(risk_prob3),
  0.9,
  na.rm = TRUE
)

occ_vals <- extract(
  risk_prob3,
  spillover_events_cases
)[,2]

mean(occ_vals > threshold, na.rm=TRUE)

# Plot Figure 4
par(mfrow = c(3, 3)) # Create a 2 x 2 plotting matrix
par(mar = c(0, 0, 0, 0)) # Set the margin on all sides to 2
plot(risk_prob11, main = "Habitat x Built",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)



plot(risk_prob9, main = "Habitat x Mining",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)






plot(risk_prob8, main = "Habitat x Bushmeat",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)




plot(risk_prob12, main = "Habitat x Forest Loss",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)



plot(risk_prob21, main = "Habitat x Mining x Forest Loss",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)



plot(risk_prob22, main = "Habitat x Mining x Conflict",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)






plot(risk_prob23, main = "Habitat x Mining x Bushmeat",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)





plot(risk_prob24, main = "Habitat x Forest Loss x Bushmeat",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)










plot(risk_prob20, main = "Habitat x Built x Conflict",
     col = hcl.colors(100, "YlOrBr")
)
#plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.5)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)
points(spillover_events_cases_before2010,add=T,col="cornsilk3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_2010to2021,add=T,col="dodgerblue3",pch=1,cex=2,lwd=1)
points(spillover_events_cases_after2022,add=T,col="firebrick2",pch=1,cex=2,lwd=1)


#quartz.save(file = "spillover_risk_composite_maps_final_selection1.pdf", type = "pdf")



# -----------------------------
# Temporal visualization of conflict exposure
# -----------------------------
par(mfrow = c(2, 4)) # Create a 2 x 2 plotting matrix
par(mar = c(0, 0, 0, 0)) # Set the margin on all sides to 2

plot(conflict_exposure, main = "Conflict Events (2000-2026 mean)",
     col = MetBrewer::met.brewer("Gauguin", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(conflict_exposure2, main = "Conflict Events (2000-2010 mean)",
     col = MetBrewer::met.brewer("Gauguin", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(conflict_exposure3, main = "Conflict Events (2011-2021 mean)",
     col = MetBrewer::met.brewer("Gauguin", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(conflict_exposure4, main = "Conflict Events (2022-2026 mean)",
     col = MetBrewer::met.brewer("Gauguin", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(conflict_exposure1, main = "Conflict Population Exposure (2014-2026 mean)",
     col = MetBrewer::met.brewer("Degas", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(conflict_exposure1a, main = "Conflict Population Exposure (2014-2021 mean)",
     col = MetBrewer::met.brewer("Degas", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(conflict_exposure1b, main = "Conflict Population Exposure (2022-2026 mean)",
     col = MetBrewer::met.brewer("Degas", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

#quartz.save(file = "conflicts_supp_fig.pdf", type = "pdf")



# -----------------------------
# Plot forest loss temporally every 5 years 
# -----------------------------

par(mfrow = c(3, 2)) # Create a 2 x 2 plotting matrix
par(mar = c(0, 0, 0, 0)) # Set the margin on all sides to 2

plot(forest_loss_temporal[[1]],   range = c(0, 1),main="Loss in 2001-2005",
     col = MetBrewer::met.brewer("Monet", n = 256, type = "continuous", direction = 1)
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

plot(forest_loss_temporal[[2]],   range = c(0, 1),main="Loss in 2006-2010",
     col = MetBrewer::met.brewer("Monet", n = 256, type = "continuous", direction = 1)
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

plot(forest_loss_temporal[[3]],   range = c(0, 1),main="Loss in 2011-2015",
     col = MetBrewer::met.brewer("Monet", n = 256, type = "continuous", direction = 1)
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

plot(forest_loss_temporal[[4]],   range = c(0, 1),main="Loss in 2016-2020",
     col = MetBrewer::met.brewer("Monet", n = 256, type = "continuous", direction = 1)
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(forest_loss_temporal[[5]],     range = c(0, 1),main="Loss in 2021-2025",
     col = MetBrewer::met.brewer("Monet", n = 256, type = "continuous", direction = 1)
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

#quartz.save(file = "forest_loss_temporal_every_5years.pdf", type = "pdf")


# Mean forest-loss intensity across the study area
forest_loss_summary <- global(
  forest_loss_temporal,
  fun = "mean",
  na.rm = TRUE
)

forest_loss_summary

forest_loss_summary2 <- global(
  forest_loss_temporal,
  fun = "mean",
  na.rm = TRUE
)

forest_loss_summary2$percentage <- 
  forest_loss_summary2$mean * 100



# Calculate area of each raster cell in km²
cell_area <- terra::cellSize(
  forest_loss_temporal[[1]],
  unit = "km"
)

forest_loss_area <- forest_loss_temporal * cell_area

forest_loss_area_summary <- terra::global(
  forest_loss_area,
  fun = "sum",
  na.rm = TRUE
)

forest_loss_area_summary

forest_loss_area_summary$period <- names(forest_loss_temporal)

forest_loss_area_summary

ebola_mask <- ifel(
  ebola_niche >= 0.7,
  1,
  NA
)

forest_loss_ebola <- mask(
  forest_loss_temporal,
  ebola_mask
)



forest_loss_summary_ebola <- terra::global(
  forest_loss_ebola,
  fun = "sum",
  na.rm = TRUE
)

forest_loss_summary_ebola$period <- names(forest_loss_ebola)


ggplot()+
  theme_bw()+
  geom_line(data=forest_loss_area_summary,aes(period,sum,colour='Study area'),group=1,size=1)+
  geom_point(data=forest_loss_area_summary,aes(period,sum,colour='Study area'),size=3)+
  geom_line(data=forest_loss_summary_ebola,aes(period,sum*100,colour='Ebola niche'),group=2,size=1)+
  geom_point(data=forest_loss_summary_ebola,aes(period,sum*100,colour='Ebola niche'),size=3)+
  scale_colour_manual(values=c('palegreen4','darkseagreen3'))+
  theme(legend.position = "bottom", legend.title = element_blank())+
  scale_y_continuous(labels=scales::comma,name="Area affected by deforestation (km²)",
                     sec.axis = sec_axis(~ . / 100,name='Ebola niche only (suitability >= 0.7)'))+
  xlab("Five-year period")
#quartz.save(file = "forest_loss_temporal_summary_area_every_5years.pdf", type = "pdf")

# -----------------------------
# plot built-up environment temporally each 5 years
# -----------------------------

built_2000 <- rast("data/GHS_BUILT_S_E2000_GLOBE_R2023A_4326_30ss_V1_0.tif")    # built-up surface
built_2000 <- mask(crop(built_2000, CentralAfrica_vect),CentralAfrica_vect)    # built-up surface

built_2005 <- rast("data/GHS_BUILT_S_E2005_GLOBE_R2023A_4326_30ss_V1_0.tif")    # built-up surface
built_2005 <- mask(crop(built_2005, CentralAfrica_vect),CentralAfrica_vect)    # built-up surface

built_2010 <- rast("data/GHS_BUILT_S_E2010_GLOBE_R2023A_4326_30ss_V1_0.tif")    # built-up surface
built_2010 <- mask(crop(built_2010, CentralAfrica_vect),CentralAfrica_vect)    # built-up surface


built_2015 <- rast("data/GHS_BUILT_S_E2015_GLOBE_R2023A_4326_30ss_V1_0.tif")    # built-up surface
built_2015 <- mask(crop(built_2015, CentralAfrica_vect),CentralAfrica_vect)    # built-up surface


built_2020 <- rast("data/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif")    # built-up surface
built_2020 <- mask(crop(built_2020, CentralAfrica_vect),CentralAfrica_vect)    # built-up surface

built_2025 <- rast("data/GHS_BUILT_S_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif")    # built-up surface
built_2025 <- mask(crop(built_2025, CentralAfrica_vect),CentralAfrica_vect)    # built-up surface

par(mfrow = c(3, 2)) # Create a 2 x 2 plotting matrix
par(mar = c(0, 0, 0, 0)) # Set the margin on all sides to 2


plot(log1p(built_2000), main = "Built-up intensity (log) (2000)",
     col = hcl.colors(100, "BuPu")
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(log1p(built_2005), main = "Built-up intensity (log) (2005)",
     col = hcl.colors(100, "BuPu")
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(log1p(built_2010), main = "Built-up intensity (log) (2010)",
     col = hcl.colors(100, "BuPu")
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(log1p(built_2015), main = "Built-up intensity (log) (2015)",
     col = hcl.colors(100, "BuPu")
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

plot(log1p(built_2020), main = "Built-up intensity (log) (2020)",
     col = hcl.colors(100, "BuPu")
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

plot(log1p(built_2025), main = "Built-up intensity (log) (2025)",
     col = hcl.colors(100, "BuPu")
)
plot(CentralAfrica_map2,add=T,border='grey70',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

#quartz.save(file = "built-up_summary_area_every_5years.pdf", type = "pdf")


built_up_change <- c(
  built_2005 - built_2000,  # 2000–2005
  built_2010 - built_2005,  # 2005–2010
  built_2015 - built_2010,  # 2010–2015
  built_2020 - built_2015,  # 2015–2020
  built_2025 - built_2020   # 2020–2025
)

names(built_up_change) <- c(
  "2000_2005",
  "2005_2010",
  "2010_2015",
  "2015_2020",
  "2020_2025"
)

built_up_change_km2 <- built_up_change / 1e6

built_up_area_summary <- terra::global(
  built_up_change_km2,
  fun = "sum",
  na.rm = TRUE
)

built_up_area_summary$period <- names(built_up_change_km2)

built_up_area_summary

built_up_change_km2 <- project(built_up_change_km2, ebola_niche,  method = "sum")

built_up_ebola <- mask(
  built_up_change_km2,
  ebola_mask
)

built_up_summary_ebola <- terra::global(
  built_up_ebola,
  fun = "sum",
  na.rm = TRUE
)

built_up_summary_ebola$period <- names(built_up_ebola)


ggplot()+
  theme_bw()+
  geom_line(data=built_up_area_summary,aes(period,sum,group=1,colour='Study area'),size=1)+
  geom_point(data=built_up_area_summary,aes(period,sum,colour='Study area'),size=3)+
  geom_line(data=built_up_summary_ebola,aes(period,sum*10,colour='Ebola niche'),group=2,size=1)+
  geom_point(data=built_up_summary_ebola,aes(period,sum*10,colour='Ebola niche'),size=3)+
  scale_colour_manual(values=c('purple4','plum4'))+
  theme(legend.position = "bottom", legend.title = element_blank())+
  scale_y_continuous(labels=scales::comma,name="New built-up area (km²)",
                     sec.axis = sec_axis(~ . / 10,name='Ebola niche only (suitability >= 0.7)'))+
  xlab("Five-year period")

#quartz.save(file = "built-up-change_summary_area_every_5years.pdf", type = "pdf")



# -----------------------------
# Figure 3 - plot covariates of spillover model
# -----------------------------
par(mfrow = c(2, 3)) # Create a 2 x 2 plotting matrix
par(mar = c(0, 0, 0, 0)) # Set the margin on all sides to 2


plot(ebola_niche, main = "Ebola ecological niche",
     col = MetBrewer::met.brewer("Demuth", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


plot(log1p(built), main = "Built-up intensity (log)",
     col = hcl.colors(100, "BuPu")
)
plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

plot(forest_loss, main = "Forest loss",
     col = MetBrewer::met.brewer("Monet", n = 256, type = "continuous", direction = 1)
)
plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)

plot(bushmeat, main = "Bushmeat activity",
     col = MetBrewer::met.brewer("Tam", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)



plot(mining_risk, main = "Mining proximity",
     col = MetBrewer::met.brewer("Paquin", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)



plot(conflict_exposure, main = "Conflict exposure",
     col = MetBrewer::met.brewer("Troy", n = 256, type = "continuous", direction = -1)
)
plot(CentralAfrica_map2,add=T,border='grey50',col=NA,lwd=0.1)
plot(CentralAfrica_map1,add=T,border='grey90',col=NA)


#quartz.save(file = "covariates_maps_FIG3_updatedebolaniche.pdf", type = "pdf")






# -----------------------------
# Figure 5 - Comput continental risk maps - reload rasters
# -----------------------------

################# Updated vs published dataset Ebola maps

ebola_niche_2016_data<-rast("data/Ebola_zoo_niche_2016.tif")

ebola_niche_2026_data<-rast("data/Ebola_zoo_niche_2026.tif")



pop <- rast("data/gpw_v4_population_density_rev11_2015_30_sec.tif")        # population density
built <- rast("data/GHS_BUILT_S_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif")    # built-up surface
mining <- rast("data/distance_to_nearest_mining_site_africa_0p05deg_osm_plus_non_osm.tif")
conflict_exposure<-rast("data/ConflictEvents_MeanAnnual_2000_2026.tif")
forest_loss<-rast("data/Africa_merged_forest_loss_medium_res.tif")
bushmeat<-rast("data/bushmeat_activity_Binomial_predictions.tif")

#crop to Africa area
pop <- mask(crop(pop, Africa_vect),Africa_vect)       # population density
pop <- project(pop, crs(Africa_vect))

built <- mask(crop(built, Africa_vect),Africa_vect)    # built-up surface
mining <- mask(crop(mining, Africa_vect),Africa_vect)
conflict_exposure<-mask(crop(conflict_exposure, Africa_vect),Africa_vect)

bushmeat<-mask(crop(bushmeat, Africa_vect),Africa_vect)
bushmeat <- project(bushmeat, crs(Africa_vect))

ebola_niche<-mask(crop(ebola_niche_2026_data, Africa_vect),Africa_vect)
ebola_niche_2016<-mask(crop(ebola_niche_2016_data, Africa_vect),Africa_vect)
forest_loss<-mask(crop(forest_loss, Africa_vect),Africa_vect)
pop <- project(pop, ebola_niche)
k <- 10
mining_risk <- exp(-mining / k)
forest_loss <- project(forest_loss, ebola_niche)
conflict_exposure <- project(conflict_exposure, ebola_niche)
built <- project(built, ebola_niche)
bushmeat <- project(bushmeat, ebola_niche)
mining_risk <- project(mining_risk, ebola_niche)



built<-log1p(built)

built_min <- global(built, "min", na.rm = TRUE)[1,1]
built_max <- global(built, "max", na.rm = TRUE)[1,1]

built <- (built - built_min) / (built_max - built_min)




conflict_min <- global(conflict_exposure, "min", na.rm = TRUE)[1,1]
conflict_max <- global(conflict_exposure, "max", na.rm = TRUE)[1,1]

conflict_exposure <- (conflict_exposure - conflict_min) / (conflict_max - conflict_min)


continental_spillover_risk_1 <- 
  #log1p(pop) * 
  log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche


continental_spillover_risk_2 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  bushmeat *
  ebola_niche


continental_spillover_risk_3 <- 
  #log1p(pop) * 
  #log1p(built) * 
  mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche

continental_spillover_risk_4 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  conflict_exposure  * 
  #forest_loss * 
  #bushmeat *
  ebola_niche


continental_spillover_risk_5 <- 
  #log1p(pop) * 
  #log1p(built) * 
  #mining_risk * 
  #conflict_displacements * 
  #conflict_exposure  * 
  forest_loss * 
  #bushmeat *
  ebola_niche



continent_risk_min1 <- global(continental_spillover_risk_1, "min", na.rm = TRUE)[1,1]
continent_risk_max1 <- global(continental_spillover_risk_1, "max", na.rm = TRUE)[1,1]

continent_risk_prob1 <- (continental_spillover_risk_1 - continent_risk_min1) / (continent_risk_max1 - continent_risk_min1)
names(continent_risk_prob1)<-"risk_habitat_built"



continent_risk_min2 <- global(continental_spillover_risk_2, "min", na.rm = TRUE)[1,1]
continent_risk_max2 <- global(continental_spillover_risk_2, "max", na.rm = TRUE)[1,1]

continent_risk_prob2 <- (continental_spillover_risk_2 - continent_risk_min2) / (continent_risk_max2 - continent_risk_min2)
names(continent_risk_prob2)<-"risk_habitat_bushmeat"



continent_risk_min3 <- global(continental_spillover_risk_3, "min", na.rm = TRUE)[1,1]
continent_risk_max3 <- global(continental_spillover_risk_3, "max", na.rm = TRUE)[1,1]

continent_risk_prob3 <- (continental_spillover_risk_3 - continent_risk_min3) / (continent_risk_max3 - continent_risk_min3)
names(continent_risk_prob3)<-"risk_habitat_mining"


continent_risk_min4 <- global(continental_spillover_risk_4, "min", na.rm = TRUE)[1,1]
continent_risk_max4 <- global(continental_spillover_risk_4, "max", na.rm = TRUE)[1,1]

continent_risk_prob4 <- (continental_spillover_risk_4 - continent_risk_min4) / (continent_risk_max4 - continent_risk_min4)
names(continent_risk_prob4)<-"risk_habitat_conflict"

continent_risk_min5 <- global(continental_spillover_risk_5, "min", na.rm = TRUE)[1,1]
continent_risk_max5 <- global(continental_spillover_risk_5, "max", na.rm = TRUE)[1,1]

continent_risk_prob5 <- (continental_spillover_risk_5 - continent_risk_min5) / (continent_risk_max5 - continent_risk_min5)
names(continent_risk_prob5)<-"risk_habitat_forest_loss"



par(mfrow = c(2, 2)) # Create a 2 x 2 plotting matrix
par(mar = c(0, 0, 0, 0)) # Set the margin on all sides to 2

plot(continent_risk_prob1, main = "Spillover Risk Scenario 1: Habitat x Built",
     col = MetBrewer::met.brewer("OKeeffe2", n = 256, type = "continuous", direction = 1)
)
plot(Africa_map1,add=T,border='black',col=NA,lwd=0.2)

plot(continent_risk_prob2, main = "Spillover Risk Scenario 2: Habitat x Bushmeat",
     col = MetBrewer::met.brewer("OKeeffe2", n = 256, type = "continuous", direction = 1)
)
plot(Africa_map1,add=T,border='black',col=NA,lwd=0.2)

plot(continent_risk_prob3, main = "Spillover Risk Scenario 3: Habitat x Mining",
     col = MetBrewer::met.brewer("OKeeffe2", n = 256, type = "continuous", direction = 1)
)
plot(Africa_map1,add=T,border='black',col=NA,lwd=0.2)

plot(continent_risk_prob4, main = "Spillover Risk Scenario 4: Habitat x Conflict",
     col = MetBrewer::met.brewer("OKeeffe2", n = 256, type = "continuous", direction = 1)
)
plot(Africa_map1,add=T,border='black',col=NA,lwd=0.2)


#quartz.save(file = "FIG5_spillover_risk_continental_maps_final_selection1.pdf", type = "pdf")











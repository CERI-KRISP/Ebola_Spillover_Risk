library(CoordinateCleaner)
library(spThin)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)

occ_filtered <- Hm_occurrences %>%
  dplyr::select(gbifID, scientificName, decimalLongitude, decimalLatitude, 
                countryCode, basisOfRecord) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) #remove exact duplicates

nrow(occ_filtered)

occ_clean <- occ_filtered %>%
  clean_coordinates(
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    countries = "countryCode",
    species = "scientificName",
    tests = c(
      "centroids",     # Within 1km of a country/province centroid — common GBIF artefact
      "duplicates",    # Duplicate lat/lon per species
      "equal",         # Lat and lon are equal — obvious error
      "gbif",          # Coordinates at GBIF headquarters in Copenhagen
      "institutions",  # Within 100m of a biodiversity institution (zoos, herbaria)
      "zeros",         # Coordinates at 0,0 or lat/lon = 0
      "seas"           # Records falling in the ocean
    ),
    value = "spatialvalid"  # Returns a flagged dataframe rather than removing rows
  )

occ_final <- occ_clean %>%
  filter(.summary == TRUE)

nrow(occ_final)

occ_final <- occ_final %>%
  filter(is.na(countryCode) | countryCode != "US")  # remove records from US

# Get Africa shapefile
africa <- ne_countries(continent = "Africa", returnclass = "sf")

# Plot cleaned occurrences
ggplot() +
  geom_sf(data = africa, fill = "grey90", colour = "grey60") +
  geom_point(data = occ_final, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "firebrick", alpha = 0.6, size = 1.5) +
  theme_minimal() +
  labs(title = "Cleaned occurrences", 
       subtitle = paste(nrow(occ_final), "records"))

occ_final <- occ_final %>%
  mutate(species = "Hypsignathus monstrosus")

#thinning to one point per 10km 
thinned <- thin(
  loc.data = occ_final,
  lat.col = "decimalLatitude",
  long.col = "decimalLongitude",
  spec.col = "species",
  thin.par = 10,        # Minimum distance between records in km — adjust based on your raster resolution
  reps = 100,           # Run 100 times and keep the rep that retains most records
  locs.thinned.list.return = TRUE,
  write.files = FALSE
)

# Extract the best thinning result (most records retained)
occ_thinned <- thinned[[1]]
nrow(occ_thinned)

#add species name and format df
hm_thinned <- occ_thinned %>%
  rename(decimalLongitude = Longitude,
         decimalLatitude = Latitude) %>%
  mutate(species = "Hypsignathus monstrosus")

write_csv(hm_thinned, "Hypsignathus_cleaned_occur.csv")



library(rnaturalearth)
library(ggplot2)
library(sf)
library(tidyterra)
library(biomod2)
library(rnaturalearthdata)
library(spThin)

df <- read_xlsx('elife_occurrence_2016.xlsx')

#plot occurrence points on map of Africa
occ_points <- st_as_sf(
  df,
  coords = c("long", "lat"),  # change if needed
  crs = 4326,
  remove = FALSE
)

africa <- ne_countries(
  scale = "medium",
  continent = "Africa",
  returnclass = "sf"
)


ggplot() +
  geom_sf(data = africa, fill = "white", color = "grey40", linewidth = 0.3) +
  geom_sf(aes(color = host), data = occ_points, size = 2, alpha = 0.7) +
  scale_color_manual(values = c(
    "animal" = "#26547c",
    "human" = "#ef476f"))+
  coord_sf(xlim = c(-20, 55), ylim = c(-40, 40), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(color = "Host"
  )


df$species <-"ebola"


#thinning to one point per 5km 
thinned <- thin(
  loc.data = df,
  lat.col = "lat",
  long.col = "long",
  spec.col = "species",
  thin.par = 5,        # Minimum distance between records in km — adjust based on your raster resolution
  reps = 100,           # Run 100 times and keep the rep that retains most records
  locs.thinned.list.return = TRUE,
  write.files = FALSE
)

# Extract the best thinning result (most records retained)
ebola_thinned <- thinned[[1]]
nrow(ebola_thinned)

#add species name and format df
ebola_final <- ebola_thinned %>%
  mutate(species = "Ebola")

write_csv(ebola_final, "ebola_cleaned_occur.csv")

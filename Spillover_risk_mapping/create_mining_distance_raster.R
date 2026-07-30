#!/usr/bin/env Rscript

# Create a 0.05-degree raster of distance to nearest mining site in Africa.
#
# This is the R equivalent of scripts/create_mining_distance_raster.py.
# It reads the same OSM and optional non-OSM mining-site CSVs, rasterizes the
# points to a regular EPSG:4326 grid, calculates distance to the nearest mining
# site with terra::distance(), masks the result to African land, and writes a
# GeoTIFF, ESRI ASCII grid, metadata JSON, and preview PNG.

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

ROOT <- "data"
OSM_CSV <- file.path(ROOT, "osm_mining", "processed", "africa_mining_sites.csv")
NON_OSM_CSV <- file.path(ROOT, "non_osm_asm", "processed", "non_osm_asm_sites_africa.csv")
AFRICA_BASEMAP <- file.path(ROOT, "osm_mining", "basemaps", "ne_110m_admin_0_countries.geojson")
DEFAULT_OUT_DIR <- file.path(ROOT, "rasters", "mining_distance_r")
NODATA <- -9999

ARTISANAL_TERMS <- c(
  "artisanal",
  "small-scale",
  "small_scale",
  "informal",
  "\"asm\"",
  " asm ",
  "orpaillage",
  "orpailleur",
  "aurifère artisanal",
  "aurifere artisanal"
)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    resolution = 0.05,
    include_non_osm = TRUE,
    artisanal_only = FALSE,
    min_lon = -20.5,
    min_lat = -36.0,
    max_lon = 55.5,
    max_lat = 38.5,
    out_dir = DEFAULT_OUT_DIR,
    no_land_mask = FALSE
  )

  for (arg in args) {
    if (arg == "--include-non-osm") opts$include_non_osm <- TRUE
    if (arg == "--osm-only") opts$include_non_osm <- FALSE
    if (arg == "--artisanal-only") opts$artisanal_only <- TRUE
    if (arg == "--no-land-mask") opts$no_land_mask <- TRUE
    if (grepl("^--resolution=", arg)) opts$resolution <- as.numeric(sub("^--resolution=", "", arg))
    if (grepl("^--min-lon=", arg)) opts$min_lon <- as.numeric(sub("^--min-lon=", "", arg))
    if (grepl("^--min-lat=", arg)) opts$min_lat <- as.numeric(sub("^--min-lat=", "", arg))
    if (grepl("^--max-lon=", arg)) opts$max_lon <- as.numeric(sub("^--max-lon=", "", arg))
    if (grepl("^--max-lat=", arg)) opts$max_lat <- as.numeric(sub("^--max-lat=", "", arg))
    if (grepl("^--out-dir=", arg)) opts$out_dir <- sub("^--out-dir=", "", arg)
  }

  opts
}

has_text <- function(x, terms) {
  x <- tolower(paste(x, collapse = " "))
  any(vapply(terms, function(term) grepl(term, x, fixed = TRUE), logical(1)))
}

is_osm_likely_artisanal <- function(df) {
  apply(
    df[, intersect(c("name", "commodity", "all_tags"), names(df)), drop = FALSE],
    1,
    has_text,
    terms = ARTISANAL_TERMS
  )
}

is_non_osm_artisanal <- function(df) {
  cols <- intersect(c("classification", "site_type", "status", "description", "all_properties"), names(df))
  apply(
    df[, cols, drop = FALSE],
    1,
    has_text,
    terms = c("artisanal", "small scale", "small-scale", "small_scale", "asm")
  )
}

read_points <- function(include_non_osm = TRUE, artisanal_only = FALSE) {
  osm <- read.csv(OSM_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  if (artisanal_only) {
    osm <- osm[is_osm_likely_artisanal(osm), , drop = FALSE]
  }
  osm_points <- data.frame(
    lon = as.numeric(osm$lon),
    lat = as.numeric(osm$lat),
    source = "osm",
    stringsAsFactors = FALSE
  )

  points <- osm_points
  if (include_non_osm && file.exists(NON_OSM_CSV)) {
    non_osm <- read.csv(NON_OSM_CSV, stringsAsFactors = FALSE, check.names = FALSE)
    if (artisanal_only) {
      non_osm <- non_osm[is_non_osm_artisanal(non_osm), , drop = FALSE]
    }
    non_osm_points <- data.frame(
      lon = as.numeric(non_osm$lon),
      lat = as.numeric(non_osm$lat),
      source = "non_osm",
      stringsAsFactors = FALSE
    )
    points <- rbind(points, non_osm_points)
  }

  points <- points[is.finite(points$lon) & is.finite(points$lat), , drop = FALSE]
  points
}

write_preview <- function(r, out_png, title) {
  png(out_png, width = 1500, height = 1200, res = 150)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  })

  vals <- values(r, mat = FALSE)
  vals <- vals[is.finite(vals) & vals != NODATA]
  zlim <- c(0, as.numeric(stats::quantile(vals, 0.98, na.rm = TRUE)))
  pal <- hcl.colors(64, "Spectral", rev = TRUE)

  plot(
    r,
    main = title,
    col = pal,
    range = zlim,
    plg = list(title = "km to nearest site"),
    axes = FALSE,
    box = FALSE,
    mar = c(2, 2, 3, 6)
  )
}

main <- function() {
  opts <- parse_args()
  dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

  points <- read_points(opts$include_non_osm, opts$artisanal_only)
  if (nrow(points) == 0) {
    stop("No source points matched the requested filters.", call. = FALSE)
  }

  template <- rast(
    ext(opts$min_lon, opts$max_lon, opts$min_lat, opts$max_lat),
    resolution = opts$resolution,
    crs = "EPSG:4326"
  )

  points$value <- 1
  points_v <- vect(points, geom = c("lon", "lat"), crs = "EPSG:4326", keepgeom = TRUE)
  site_raster <- rasterize(points_v, template, field = "value", background = NA)

  # terra::distance returns distance to the nearest non-NA raster cell. Mining
  # site cells are 1 and all other cells are NA; unit = "km" returns kilometres.
  distance_km <- distance(site_raster, unit = "km", method = "haversine")

  land_mask_path <- NULL
  if (!opts$no_land_mask && file.exists(AFRICA_BASEMAP)) {
    africa <- vect(AFRICA_BASEMAP)
    africa <- africa[africa$CONTINENT == "Africa", ]
    distance_km <- mask(distance_km, africa)
    land_mask_path <- AFRICA_BASEMAP
  }

  distance_km[is.na(distance_km)] <- NODATA

  source_suffix <- if (opts$include_non_osm) "osm_plus_non_osm" else "osm"
  suffix <- if (opts$artisanal_only) {
    paste0("artisanal_", source_suffix)
  } else {
    source_suffix
  }

  base_name <- paste0("distance_to_nearest_mining_site_africa_0p05deg_", suffix)
  tif_path <- file.path(opts$out_dir, paste0(base_name, ".tif"))
  asc_path <- file.path(opts$out_dir, paste0(base_name, ".asc"))
  png_path <- file.path(opts$out_dir, paste0(base_name, "_preview.png"))
  meta_path <- file.path(opts$out_dir, paste0(base_name, "_metadata.json"))

  writeRaster(distance_km, tif_path, overwrite = TRUE, NAflag = NODATA)
  writeRaster(distance_km, asc_path, overwrite = TRUE, NAflag = NODATA, filetype = "AAIGrid")

  title <- if (opts$artisanal_only) {
    "Distance to nearest artisanal mining site in Africa"
  } else {
    "Distance to nearest mining site in Africa"
  }
  write_preview(distance_km, png_path, title)

  metadata <- list(
    description = "Distance to nearest rasterized mining-site cell in Africa",
    site_filter = if (opts$artisanal_only) "likely artisanal / ASM" else "all mining-related sites",
    distance_units = "kilometers",
    crs = "EPSG:4326",
    resolution_degrees = opts$resolution,
    extent = list(
      min_lon = opts$min_lon,
      min_lat = opts$min_lat,
      max_lon = opts$max_lon,
      max_lat = opts$max_lat
    ),
    rows = nrow(distance_km),
    cols = ncol(distance_km),
    source_points_total = nrow(points),
    source_cells_with_site = as.integer(global(!is.na(site_raster), "sum", na.rm = TRUE)[1, 1]),
    source_layers = c(
      OSM_CSV,
      if (opts$include_non_osm && file.exists(NON_OSM_CSV)) NON_OSM_CSV else NULL
    ),
    land_mask = land_mask_path,
    nodata = NODATA,
    outputs = list(
      geotiff = tif_path,
      ascii_grid = asc_path,
      preview_png = png_path
    ),
    method_note = paste(
      "Mining-site points are rasterized to a 0.05-degree grid.",
      "Distance to the nearest site cell is calculated with terra::distance using the haversine method and returned in kilometres.",
      "Cells outside the Africa land mask are assigned NoData."
    )
  )

  write_json(metadata, meta_path, pretty = TRUE, auto_unbox = TRUE)
  cat(toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), "\n")
}

main()

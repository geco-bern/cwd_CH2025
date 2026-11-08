# Biogeographical regions of Switzerland (FOEN, revision 2020) on the CH2025 1 km grid
# for analyses of the PCWD results per region, e.g.
#   r_max <- rast(here("data", "pcwd_annmax", "grid_ch_pcwd_annmax_<model>_<period>.nc"),
#                 subds = "pcwd_annmax")
#   zones <- biogeoregio_grid(r_max)              # region of every grid cell
#   by_region(r_max, zones, "mean")               # mean annual maximum per region and year
# needs: library(terra), library(here)

# shapefile in LV95 (EPSG:2056) from data.geo.admin.ch, sha256 from the STAC item (file:checksum)
bgr_url    <- "https://data.geo.admin.ch/ch.bafu.biogeographische_regionen/biogeographische_regionen/biogeographische_regionen_2056.shp.zip"
bgr_sha256 <- "f3c0567d3ed460295981ba73028b1e283b7d9ab3e4b6c070c9fde38aa2acaaf0"

# English names of regions 1-6 (the dataset only has German, French and Italian names)
bgr_names <- c("Jura", "Central Plateau", "Northern Alps", "Western Central Alps",
               "Eastern Central Alps", "Southern Alps")

# Read the 6 regions 
# returns a SpatVector in LV95 with the columns id and name
read_biogeoregio <- function(dir = here("data-raw", "biogeoregio")) {
  shp <- file.path(dir, "N2020_Revision_BiogeoRegion.shp")

  if (!file.exists(shp)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    zip_file <- file.path(dir, basename(bgr_url))
    download.file(bgr_url, zip_file, mode = "wb", method = "libcurl", quiet = TRUE)
    stopifnot("biogeographical regions checksum mismatch" =
                digest::digest(zip_file, algo = "sha256", file = TRUE) == bgr_sha256)
    unzip(zip_file, exdir = dir)
    unlink(zip_file)
  }

  v <- aggregate(vect(shp), by = "RegionNumm")      # dissolves the subregions of each region
  v <- v[order(v$RegionNumm), ]
  v$id   <- as.integer(v$RegionNumm)
  v$name <- bgr_names
  v[, c("id", "name")]
}

# Regions on the CH2025 1 km grid 
# template: a raster on the CH2025 grid (e.g. a PCWD file); default: the empty CH2025 grid
# every cell takes the region that contains its centre; cells outside Switzerland are NA
# returns a categorical SpatRaster with one layer "region"
biogeoregio_grid <- function(template = rast(nrows = 240, ncols = 370, xmin = 2474000, xmax = 2844000,
                                             ymin = 1064000, ymax = 1304000, crs = "EPSG:2056")) {
  v <- read_biogeoregio()
  z <- rasterize(v, rast(template, nlyrs = 1), field = "id")
  levels(z) <- as.data.frame(v)
  setNames(z, "region")
}

# ---- Statistic per region ------------------------------------------------------------
# r: raster on the CH2025 grid (e.g. annual maximum PCWD, one layer per year)
# zones: output of biogeoregio_grid(); fun: a function name or function accepted by terra::zonal()
# returns a long data frame: one row per region and layer, with the number of cells used
# (region is a factor in the order of the region ids)
by_region <- function(r, zones, fun = "mean") {
  stat <- zonal(r, zones, fun = fun, na.rm = TRUE)
  n    <- zonal(!is.na(r[[1]]), zones, fun = "sum")
  i    <- order(match(stat[[1]], bgr_names))
  data.frame(region  = factor(stat[[1]][i], levels = bgr_names),
             n_cells = n[[2]][match(stat[[1]][i], n[[1]])],
             layer   = rep(seq_len(nlyr(r)), each = length(i)),
             value   = unlist(stat[i, -1], use.names = FALSE))
}

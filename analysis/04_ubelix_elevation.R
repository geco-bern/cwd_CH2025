#!/usr/bin/env Rscript
#SBATCH --job-name=elevation
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=elevation_%j.log

# Download DHM25 and prepare the elevation (LV03, 200 m) for matching with the CH2025 gridded climate files (LV95, 1 km, NetCDF).
# creates a .tif file with 1 km grid (LV95) and one layer "elv" (elevation in m)
# creates a mismatch map and table (climate data but no elevation) in the figures and tables folders

library(terra)

# Paths
zip_url  <- "https://data.geo.admin.ch/ch.swisstopo.digitales-hoehenmodell_25/data.zip"
dir_dhm  <- file.path("data-raw", "DHM25")
dhm_file <- file.path(dir_dhm, "DHM200.asc")
pr_file  <- file.path("data-raw", "CH2025", "DAILY_GRIDDED",
                      "ogd-climate-scenarios-ch2025-grid_ch_pr_clmcom-cclm4-cccma_gwl1.5.nc")
out_file <- file.path("data", "elevation", "elevation_1k_grid.tif")
fig_file <- file.path("figures", "mismatch_climate_vs_elevation.jpg")
csv_file <- file.path("tables", "mismatch_climate_vs_elevation.csv")

for (d in c(dir_dhm, dirname(c(out_file, fig_file, csv_file)))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Download and unzip DHM25
if (!file.exists(dhm_file)) {
  zip_file <- file.path(dir_dhm, "data.zip")
  download.file(zip_url, zip_file, mode = "wb", method = "libcurl", quiet = TRUE)
  sha <- sub(" .*", "", system2("sha256sum", shQuote(zip_file), stdout = TRUE))
  stopifnot("DHM25 checksum mismatch" = any(grepl(paste0("^x-amz-meta-sha256: *", sha),
                                                  curlGetHeaders(zip_url), ignore.case = TRUE)))
  unzip(zip_file, exdir = dir_dhm)
  unlink(zip_file)
}

# Read DHM25
dem_raw <- rast(dhm_file)
crs(dem_raw) <- "EPSG:21781"       # the ASCII grid carries no CRS: DHM25 is in LV03
dem_raw[dem_raw < -100] <- NA      # nodata: values below the lowest point in Switzerland

# Resample to the climate grid: CRS LV03 -> LV95, mean of the 200 m cells within each 1 km cell
r0 <- rast(nrows = 240, ncols = 370, xmin = 2474000, xmax = 2844000,
           ymin = 1064000, ymax = 1304000, crs = "EPSG:2056")
dem <- project(dem_raw, r0, method = "average")
names(dem) <- "elv"

print(global(dem, c("min", "max", "mean"), na.rm = TRUE))  # expect ~200-4500 m
writeRaster(dem, out_file, overwrite = TRUE)

# Quality check: cells with climate data (first layer of a climate file) but no elevation
mismatch <- !is.na(rast(pr_file, lyrs = 1)) & is.na(dem)
write.csv(data.frame(n_cells = global(mismatch, "sum")$sum), csv_file, row.names = FALSE)

# Save mismatch map
cols <- c("lightgreen", "firebrick")

jpeg(fig_file, width = 1800, height = 1300, res = 200, type = "cairo") 
plot(mismatch, main = "Grid cells with climate data but no elevation (CH2025, 1 km grid)",
     col = cols, legend = FALSE)
legend("bottomright", fill = cols, bty = "n",
       legend = c("Climate data and elevation available", "Climate data, no elevation"))
dev.off()

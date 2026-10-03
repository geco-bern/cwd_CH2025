# ------------------------------------------------------------------
# Prepare DHM25 elevation (LV03, 200 m) for matching with the CH2025
# gridded climate files (LV95, 1 km, NetCDF).
#
# Output: elevation_1k_grid.tif  (1 layer "elv", same grid as the CH2025 .nc files)
#         mismatch_climate_vs_elevation.jpg
# ------------------------------------------------------------------

library(terra)
library(here)

# ---- 1. Paths -----------------------------------------------------
# here() starts at the project root (folder with .git, an .Rproj or a .here file).
# Each folder is a separate argument. EDIT the folder names if your structure differs.
dhm_file <- here("data-raw", "DHM25", "DHM200.asc")
pr_file  <- here("data-raw", "CH2025", "DAILY_GRIDDED",
                 "ogd-climate-scenarios-ch2025-grid_ch_pr_clmcom-cclm4-cccma_gwl1.5.nc")
out_file <- here("data", "elevation", "elevation_1k_grid.tif")
fig_file <- here("figures", "mismatch_climate_vs_elevation.jpg")

stopifnot("DHM25 file not found" = file.exists(dhm_file),
          "climate .nc file not found" = file.exists(pr_file))

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(fig_file), showWarnings = FALSE, recursive = TRUE)

# ---- 2. Read DHM25 ------------------------------------------------
dem_raw <- rast(dhm_file)

# ASCII grids often carry no CRS: DHM25 is in LV03 (EPSG:21781)
crs(dem_raw) <- "EPSG:21781"

# make sure nodata is NA (values below the lowest point in Switzerland)
dem_raw[dem_raw < -100] <- NA
names(dem_raw) <- "elv"

print(dem_raw)   # extent should be in LV03 range (x ~480000-840000, y ~70000-300000)
plot(dem_raw, main = "DHM25 (LV03)")

# ---- 3. Climate grid as template (matches CH2025 daily gridded data) ----
r0 <- rast(nrows = 240, ncols = 370,
           xmin = 2474000, xmax = 2844000,
           ymin = 1064000, ymax = 1304000,
           crs  = "EPSG:2056")

# ---- 4. Resample DHM25 to the climate grid -----------------------
# project() changes CRS (LV03 -> LV95) and snaps to the template grid.
# method = "average": mean elevation of all 200 m cells within each 1 km cell
dem <- project(dem_raw, r0, method = "average")
names(dem) <- "elv"

# ---- 5. Checks ----------------------------------------------------
stopifnot(compareGeom(r0, dem, stopOnError = FALSE))
print(global(dem, c("min", "max", "mean"), na.rm = TRUE))  # expect ~200-4500 m
plot(dem, main = "Elevation on 1km grid (m)")

# first layer of a climate file, to see where it has data
pr1 <- rast(pr_file, lyrs = 1)

ok_dem <- !is.na(dem)
ok_pr  <- !is.na(pr1)

print(global(ok_dem, "sum"))      # cells with elevation
print(global(ok_pr,  "sum"))      # cells with climate data

mismatch <- ok_pr & !ok_dem       # climate data but no elevation
print(global(mismatch, "sum"))
plot(mismatch)                    # where are they?

# ---- 6. Save mismatch map ----------------------------------------
cols <- c("lightgreen", "firebrick")

jpeg(fig_file, width = 1800, height = 1300, res = 200)
plot(mismatch,
     main   = "Grid cells with climate data but no elevation (CH2025, 1 km grid)",
     col    = cols,
     legend = FALSE)
legend("bottomright",
       legend = c("Climate data and elevation available",
                  "Climate data, no elevation"),
       fill = cols, bty = "n")
dev.off()

# ---- 7. Save elevation --------------------------------------------
writeRaster(dem, out_file, overwrite = TRUE)
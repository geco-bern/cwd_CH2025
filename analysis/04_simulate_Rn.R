# run 03_model_Rn.R first to create data/radiation_models/rg_model_rf.rds and rn_model_rf.rds
# applies the two-step random forest on the CH2025 daily grid:
# Tmin, Tmax, Prec, Rsw-cc, elevation -> Rg (step 1) -> Rn (step 2), written to one NetCDF file

library(here)
library(terra)
library(ncdf4)
library(meteoland)
library(ranger)

# ---- 1. Settings ----------------------------------------------------------
dir_grid   <- here("data-raw", "CH2025", "DAILY_GRIDDED")
model      <- "clmcom-cclm4-cccma"
period     <- "gwl1.5"
block_rows <- 10                 # rows of the grid processed at once
n_threads  <- max(1, parallel::detectCores() - 1)
run_full   <- TRUE              # FALSE: only the test block, TRUE: the whole grid

out_dir  <- here("data", "netradiation_sim")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("grid_ch_rn_%s_%s.nc", model, period))
rsw_file <- here("data", "rsw-cc", "meteoland_grid_ch_rswcc.nc")

nc_path <- function(var)
  file.path(dir_grid, sprintf("ogd-climate-scenarios-ch2025-grid_ch_%s_%s_%s.nc",
                              var, model, period))

# ---- 1b. Layout of the output NetCDF files ---------------------------------
# regular 1 km grid in CH1903+ / LV95, coordinates = cell centres
res_m  <- 1000
E_vals <- seq(2474500, 2843500, by = res_m)       # easting,  370 cells, west to east
N_vals <- seq(1064500, 1303500, by = res_m)       # northing, 240 cells, south to north
n_days <- 30 * 365                                # 30 years, 365_day calendar (no Feb 29)
fill   <- -999.99                                 # value outside Switzerland

dim_E    <- ncdim_def("E", "meters", E_vals, longname = "swiss easting (lv95)")
dim_N    <- ncdim_def("N", "meters", N_vals, longname = "swiss northing (lv95)")
dim_time <- ncdim_def("time", "days since 1-1-1 00:00:00", 0:(n_days - 1),
                      unlim = TRUE, calendar = "365_day", longname = "time")
dim_doy  <- ncdim_def("doy", "day of year", 1:365,
                      longname = "day of year (365_day calendar, no Feb 29)")

# CF grid mapping of LV95 (EPSG:2056)
crs_atts <- list(grid_mapping_name                 = "oblique_mercator",
                 longitude_of_projection_origin    = 7.43958333,
                 latitude_of_projection_origin     = 46.9524056,
                 azimuth_of_central_line           = 90,
                 scale_factor_at_projection_origin = 1,
                 false_easting                     = 2600000,
                 false_northing                    = 1200000,
                 reference_ellipsoid_name          = "Bessel 1841",
                 semi_major_axis                   = 6377397.155,
                 inverse_flattening                = 299.1528128,
                 horizontal_datum_name             = "CH1903+",
                 prime_meridian_name               = "Greenwich",
                 crs_wkt                           = crs("EPSG:2056"),
                 spatial_ref                       = crs("EPSG:2056"))

# ---- 2. Inputs: climate grids, elevation, models --------------------------
r_tasmin <- rast(nc_path("tasmin"))
r_tasmax <- rast(nc_path("tasmax"))
r_pr     <- rast(nc_path("pr"))
elv      <- rast(here("data", "elevation", "elevation_1k_grid.tif"))
stopifnot(compareGeom(r_tasmin, r_tasmax), compareGeom(r_tasmin, r_pr),
          compareGeom(r_tasmin, elv))

# the input grid must have the layout defined in 1b
stopifnot(ncol(r_tasmin) == length(E_vals), nrow(r_tasmin) == length(N_vals),
          nlyr(r_tasmin) == n_days, all(res(r_tasmin) == res_m),
          all(as.vector(ext(r_tasmin)) == c(range(E_vals), range(N_vals)) + c(-1, 1, -1, 1) * res_m / 2))

dir_model <- here("data", "radiation_models")
rg_model  <- readRDS(file.path(dir_model, "rg_model_rf.rds"))
rn_model  <- readRDS(file.path(dir_model, "rn_model_rf.rds"))
print(rg_model$predictors); print(rn_model$predictors)

doy <- rep(1:365, length.out = n_days)            # 365_day calendar, no Feb 29

# ---- 3. Latitude and elevation of every cell ------------------------------
elv_vals <- values(elv)[, 1]                      # grid order, north to south
nc_ll    <- nc_open(nc_path("tasmin"))
lon      <- ncvar_get(nc_ll, "lon")               # E x N, N runs south to north
lat      <- ncvar_get(nc_ll, "lat")
nc_close(nc_ll)
lat_vals <- as.vector(lat[, ncol(lat):1])         # flip to terra order (north first)

# ---- 4. Rsw-cc: same function as R/calc_rsw-cc.R (flat surface, W m-2, 365 days) --
d_2025 <- seq(as.Date("2025-01-01"), by = "day", length.out = 365)
J      <- radiation_dateStringToJulianDays(format(d_2025, "%Y-%m-%d"))
delta  <- vapply(J, radiation_solarDeclination, numeric(1))
sc     <- vapply(J, radiation_solarConstant, numeric(1))
calc_rsw_cc <- function(lat_deg) {
  pot <- mapply(function(s, dl)
    radiation_potentialRadiation(s, lat_deg * pi / 180, 0, 0, dl), sc, delta)
  pot * 1e6 / 86400                               # MJ m-2 day-1 -> W m-2
}
print(round(calc_rsw_cc(46.81)[c(1, 172, 365)], 2))  # PAY: compare with meteoland_pay_rswcc.csv

# ---- 4a. Function: create an empty NetCDF file with the layout from 1b ---------
# dim3: dim_time or dim_doy; var_atts / glob_atts: extra attributes (named lists)
create_nc <- function(file, var, units, longname, dim3, chunks = NA,
                      var_atts = list(), glob_atts = list()) {
  var_lon <- ncvar_def("lon", "degrees_east",  list(dim_E, dim_N), prec = "float",
                       longname = "longitude coordinate", compression = 4)
  var_lat <- ncvar_def("lat", "degrees_north", list(dim_E, dim_N), prec = "float",
                       longname = "latitude coordinate", compression = 4)
  var_crs <- ncvar_def("swiss_lv95_coordinates", "", list(), prec = "integer")
  var_out <- ncvar_def(var, units, list(dim_E, dim_N, dim3), missval = fill,
                       prec = "float", longname = longname,
                       chunksizes = chunks, compression = 4)

  nc <- nc_create(file, list(var_lon, var_lat, var_crs, var_out), force_v4 = TRUE)
  ncvar_put(nc, "lon", lon)
  ncvar_put(nc, "lat", lat)

  ncatt_put(nc, "E", "standard_name", "projection_x_coordinate"); ncatt_put(nc, "E", "axis", "X")
  ncatt_put(nc, "N", "standard_name", "projection_y_coordinate"); ncatt_put(nc, "N", "axis", "Y")
  if (dim3$name == "time") {
    ncatt_put(nc, "time", "standard_name", "time"); ncatt_put(nc, "time", "axis", "T")
  }
  ncatt_put(nc, "lon", "standard_name", "longitude")
  ncatt_put(nc, "lat", "standard_name", "latitude")
  for (a in names(crs_atts)) ncatt_put(nc, "swiss_lv95_coordinates", a, crs_atts[[a]])

  ncatt_put(nc, var, "grid_mapping", "swiss_lv95_coordinates")
  ncatt_put(nc, var, "coordinates",  "lat lon")
  for (a in names(var_atts)) ncatt_put(nc, var, a, var_atts[[a]])

  ncatt_put(nc, 0, "Conventions",   "CF-1.8")
  ncatt_put(nc, 0, "spatial",       "regular 1km grid")
  ncatt_put(nc, 0, "source",        "analysis/04_simulate_Rn.R")
  ncatt_put(nc, 0, "creation_date", format(Sys.Date()))
  for (a in names(glob_atts)) ncatt_put(nc, 0, a, glob_atts[[a]])
  nc_close(nc)
}

# ---- 4b. Save gridded Rsw-cc (365 days, identical for every year) -------------------
if (!file.exists(rsw_file)) {
  land     <- which(!is.na(values(r_tasmin[[1]])[, 1]))       # cells with data (Switzerland)
  rsw_grid <- matrix(NA_real_, nrow = ncell(r_tasmin), ncol = 365)
  rsw_grid[land, ] <- t(vapply(lat_vals[land], calc_rsw_cc, numeric(365)))

  arr <- array(rsw_grid, dim = c(length(E_vals), length(N_vals), 365))  # E x rows (north first) x doy
  arr <- arr[, length(N_vals):1, ]                                       # the file runs south to north
  arr[is.na(arr)] <- fill

  create_nc(rsw_file, "rsw_cc", "W m-2",
            "Daily potential shortwave radiation on a flat surface (Rsw-cc)", dim_doy,
            var_atts  = list(method = "meteoland::radiation_potentialRadiation, slope 0, aspect 0, latitude of each cell (as R/calc_rsw-cc.R)"),
            glob_atts = list(title = "Daily potential shortwave radiation (Rsw-cc) on the CH2025 1 km grid"))
  nc <- nc_open(rsw_file, write = TRUE)
  ncvar_put(nc, "rsw_cc", arr)
  nc_close(nc)
}

# ---- 5. Two-step prediction for one block of rows ---------------------------
# returns a matrix: cells of the block x days (NA outside Switzerland)
process_block <- function(s) {
  n <- min(block_rows, nrow(r_tasmin) - s + 1)

  readStart(r_tasmin); readStart(r_tasmax); readStart(r_pr)
  tmin <- readValues(r_tasmin, row = s, nrows = n, mat = TRUE)
  tmax <- readValues(r_tasmax, row = s, nrows = n, mat = TRUE)
  pr   <- readValues(r_pr,     row = s, nrows = n, mat = TRUE)
  readStop(r_tasmin); readStop(r_tasmax); readStop(r_pr)

  rn_blk <- matrix(NA_real_, nrow = nrow(tmin), ncol = n_days)
  idx    <- which(!is.na(tmin[, 1]))              # cells with data
  if (length(idx) == 0) return(rn_blk)

  cell_id <- (s - 1) * ncol(r_tasmin) + idx       # cell numbers in the full grid
  rsw     <- t(vapply(lat_vals[cell_id], calc_rsw_cc, numeric(365)))  # cells x 365

  # one year (365 days) at a time to keep the data frame small
  for (y in seq(1, n_days, by = 365)) {
    days <- y:(y + 364)
    d <- data.frame(tasmin = as.vector(tmin[idx, days]),
                    tasmax = as.vector(tmax[idx, days]),
                    pr_t   = log1p(as.vector(pr[idx, days])),
                    rsw_cc = as.vector(rsw[, doy[days]]),
                    elev   = rep(elv_vals[cell_id], times = 365))
    ok <- complete.cases(d)
    d  <- d[ok, ]

    # step 1: Rg (clipped at 0, as in 03_model_Rn.R)
    d$rg_sim <- pmax(predict(rg_model$model, d, num.threads = n_threads)$predictions, 0)
    # step 2: Rn
    rn <- rep(NA_real_, length(ok))
    rn[ok] <- predict(rn_model$model, d, num.threads = n_threads)$predictions

    rn_blk[idx, days] <- matrix(rn, nrow = length(idx))
  }
  rn_blk
}

# ---- 6. Create the Rn NetCDF file ------------------------------------------------
create_nc(out_file, "rn", "W m-2", "Daily mean net radiation", dim_time,
          chunks   = c(length(E_vals), block_rows, 365),   # one chunk = one block of rows x one year
          var_atts = list(standard_name = "surface_net_downward_radiative_flux",
                          method = "Two-step random forest: Rg ~ Tmin + Tmax + Prec + Rsw-cc + elevation (CH2025 DAILY_LOCAL), Rn ~ Rg + Tmin + Tmax + Prec + Rsw-cc (SwissMetNet ALT, MAG, PAY, TAE)"),
          glob_atts = list(title                 = "Daily net radiation simulated from Climate CH2025 DAILY-GRIDDED (tasmin, tasmax, pr)",
                           GWL                   = toupper(period),
                           modelchain            = toupper(model),
                           frequency             = "day",
                           qm_calibration_period = "1991-2020",
                           CORDEX_domain         = "EUR11"))

# ---- 7. Write one block into the NetCDF file ----------------------------------
write_block <- function(file, rn_blk, s) {
  ncols <- length(E_vals)
  n     <- nrow(rn_blk) / ncols                    # rows in this block

  arr <- array(rn_blk, dim = c(ncols, n, n_days))  # E x rows (north first) x time
  arr <- arr[, n:1, , drop = FALSE]                # the file runs south to north
  arr[is.na(arr)] <- fill

  nc <- nc_open(file, write = TRUE)
  ncvar_put(nc, "rn", arr,
            start = c(1, length(N_vals) - (s + n - 1) + 1, 1),
            count = c(ncols, n, n_days))
  nc_close(nc)
}

# ---- 8. Test on the block with the most land cells -----------------------------
starts     <- seq(1, nrow(r_tasmin), by = block_rows)
hd         <- as.matrix(!is.na(r_tasmin[[1]]), wide = TRUE)
n_by_block <- sapply(starts, function(s) sum(hd[s:min(s + block_rows - 1, nrow(hd)), ]))
s_test     <- starts[which.max(n_by_block)]

t_test <- system.time({
  rn_blk <- process_block(s_test)
  write_block(out_file, rn_blk, s_test)
})
print(t_test)
cat("estimated full run (min):", round(t_test[["elapsed"]] * sum(n_by_block) / max(n_by_block) / 60), "\n")

land <- !is.na(rn_blk[, 1])
print(summary(as.vector(rn_blk[land, ])))                   # Rn values (W m-2)
print(round(rowMeans(matrix(colMeans(rn_blk[land, ]), nrow = 365)), 1)[c(1, 91, 182, 274)])  # mean Rn on doy 1, 91, 182, 274

r_out <- rast(out_file)
print(r_out)
print(compareGeom(r_out, r_tasmin))                         # same grid as the input?
chk <- values(r_out, row = s_test, nrows = block_rows, mat = TRUE)
print(all.equal(chk[land, ], rn_blk[land, ], tolerance = 1e-4, check.attributes = FALSE))

# ---- 9. Full run over all blocks ------------------------------------------------
if (run_full) {
  for (s in starts[n_by_block > 0]) {
    message("block starting at row ", s, " (", Sys.time(), ")")
    write_block(out_file, process_block(s), s)
  }
}

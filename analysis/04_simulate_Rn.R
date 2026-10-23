# run 03_model_Rn.R first to create data/rg_model_rf.rds and data/rn_model_rf.rds
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
run_full   <- FALSE              # FALSE: only the test block, TRUE: the whole grid

out_dir  <- here("data", "netradiation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("grid_ch_rn_%s_%s.nc", model, period))

nc_path <- function(var)
  file.path(dir_grid, sprintf("ogd-climate-scenarios-ch2025-grid_ch_%s_%s_%s.nc",
                              var, model, period))

# ---- 2. Inputs: climate grids, elevation, models --------------------------
r_tasmin <- rast(nc_path("tasmin"))
r_tasmax <- rast(nc_path("tasmax"))
r_pr     <- rast(nc_path("pr"))
elv      <- rast(here("data", "elevation", "elevation_1k_grid.tif"))
stopifnot(compareGeom(r_tasmin, r_tasmax), compareGeom(r_tasmin, r_pr),
          compareGeom(r_tasmin, elv))

rg_model <- readRDS(here("data", "rg_model_rf.rds"))
rn_model <- readRDS(here("data", "rn_model_rf.rds"))
print(rg_model$predictors); print(rn_model$predictors)

n_days <- nlyr(r_tasmin)                          # 10950 = 30 years x 365 days
doy    <- rep(1:365, length.out = n_days)         # 365_day calendar, no Feb 29

# ---- 3. Latitude and elevation of every cell (grid order, north to south) --
elv_vals <- values(elv)[, 1]
nc_lat   <- nc_open(nc_path("tasmin"))
lat      <- ncvar_get(nc_lat, "lat")              # E x N, N runs south to north
nc_close(nc_lat)
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

# ---- 4b. Save gridded Rsw-cc (365 days, identical for every year) -------------------
rsw_file <- here("data", "rsw-cc", "meteoland_grid_ch_rswcc.nc")
if (!file.exists(rsw_file)) {
  land     <- which(!is.na(values(r_tasmin[[1]])[, 1]))       # cells with data (Switzerland)
  rsw_grid <- matrix(NA_real_, nrow = ncell(r_tasmin), ncol = 365)
  rsw_grid[land, ] <- t(vapply(lat_vals[land], calc_rsw_cc, numeric(365)))

  arr <- array(rsw_grid, dim = c(ncol(r_tasmin), nrow(r_tasmin), 365))  # E x rows (north first) x doy
  arr <- arr[, nrow(r_tasmin):1, ]                                       # the file runs south to north
  arr[is.na(arr)] <- -999.99

  nc_in   <- nc_open(nc_path("tasmin"))
  dim_E   <- ncdim_def("E", "meters", nc_in$dim$E$vals)
  dim_N   <- ncdim_def("N", "meters", nc_in$dim$N$vals)
  dim_doy <- ncdim_def("doy", "day of year", 1:365,
                       longname = "day of year (365_day calendar, no Feb 29)")
  var_lon <- ncvar_def("lon", "degrees_east",  list(dim_E, dim_N), prec = "float",
                       longname = "longitude coordinate", compression = 4)
  var_lat <- ncvar_def("lat", "degrees_north", list(dim_E, dim_N), prec = "float",
                       longname = "latitude coordinate", compression = 4)
  var_crs <- ncvar_def("swiss_lv95_coordinates", "", list(), prec = "integer")
  var_rsw <- ncvar_def("rsw_cc", "W m-2", list(dim_E, dim_N, dim_doy),
                       missval = -999.99, prec = "float",
                       longname = "Daily potential shortwave radiation on a flat surface (Rsw-cc)",
                       compression = 4)

  nc_out <- nc_create(rsw_file, list(var_lon, var_lat, var_crs, var_rsw), force_v4 = TRUE)
  ncvar_put(nc_out, "lon", ncvar_get(nc_in, "lon"))
  ncvar_put(nc_out, "lat", ncvar_get(nc_in, "lat"))
  ncvar_put(nc_out, "rsw_cc", arr)
  for (v in c("lon", "lat", "swiss_lv95_coordinates", "E", "N")) {
    att <- ncatt_get(nc_in, v)
    for (a in setdiff(names(att), c("units", "_FillValue", "missing_value")))
      ncatt_put(nc_out, v, a, att[[a]])
  }
  ncatt_put(nc_out, "rsw_cc", "grid_mapping", "swiss_lv95_coordinates")
  ncatt_put(nc_out, "rsw_cc", "coordinates",  "lat lon")
  ncatt_put(nc_out, "rsw_cc", "method",
            "meteoland::radiation_potentialRadiation, slope 0, aspect 0, latitude of each cell (as R/calc_rsw-cc.R)")
  ncatt_put(nc_out, 0, "title", "Daily potential shortwave radiation (Rsw-cc) on the CH2025 1 km grid")
  ncatt_put(nc_out, 0, "source", "analysis/04_simulate_Rn.R")
  ncatt_put(nc_out, 0, "creation_date", format(Sys.Date()))
  nc_close(nc_out); nc_close(nc_in)
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

# ---- 6. Create the NetCDF file (structure of the CH2025 grid files) -----------
nc_in <- nc_open(nc_path("tasmin"))

dim_E    <- ncdim_def("E", "meters", nc_in$dim$E$vals)
dim_N    <- ncdim_def("N", "meters", nc_in$dim$N$vals)
dim_time <- ncdim_def("time", nc_in$dim$time$units, nc_in$dim$time$vals,
                      unlim = TRUE, calendar = ncatt_get(nc_in, "time", "calendar")$value)

fill <- -999.99                                   # same fill value as CH2025

var_lon <- ncvar_def("lon", "degrees_east",  list(dim_E, dim_N), prec = "float",
                     longname = "longitude coordinate", compression = 4)
var_lat <- ncvar_def("lat", "degrees_north", list(dim_E, dim_N), prec = "float",
                     longname = "latitude coordinate", compression = 4)
var_crs <- ncvar_def("swiss_lv95_coordinates", "", list(), prec = "integer")
var_rn  <- ncvar_def("rn", "W m-2", list(dim_E, dim_N, dim_time),
                     missval = fill, prec = "float",
                     longname = "Daily mean net radiation",
                     chunksizes = c(370, block_rows, 365),  # one chunk = one block of rows x one year
                     compression = 4)

nc_out <- nc_create(out_file, list(var_lon, var_lat, var_crs, var_rn), force_v4 = TRUE)

# coordinates and projection copied from the input
ncvar_put(nc_out, "lon", ncvar_get(nc_in, "lon"))
ncvar_put(nc_out, "lat", ncvar_get(nc_in, "lat"))
for (v in c("lon", "lat", "swiss_lv95_coordinates", "E", "N", "time")) {
  att <- ncatt_get(nc_in, v)
  for (a in setdiff(names(att), c("units", "_FillValue", "missing_value", "calendar")))
    ncatt_put(nc_out, v, a, att[[a]])
}

# attributes of the Rn variable
ncatt_put(nc_out, "rn", "standard_name", "surface_net_downward_radiative_flux")
ncatt_put(nc_out, "rn", "grid_mapping",  "swiss_lv95_coordinates")
ncatt_put(nc_out, "rn", "coordinates",   "lat lon")
ncatt_put(nc_out, "rn", "method",
          "Two-step random forest: Rg ~ Tmin + Tmax + Prec + Rsw-cc + elevation (CH2025 DAILY_LOCAL), Rn ~ Rg + Tmin + Tmax + Prec + Rsw-cc (SwissMetNet ALT, MAG, PAY, TAE)")

# global attributes copied from the input, title and source adjusted
glob <- ncatt_get(nc_in, 0)
for (a in c("Conventions", "GWL", "frequency", "qm_calibration_period", "spatial",
            "modelchain", "CORDEX_domain", "license"))
  ncatt_put(nc_out, 0, a, glob[[a]])
ncatt_put(nc_out, 0, "title", "Daily net radiation simulated from Climate CH2025 DAILY-GRIDDED")
ncatt_put(nc_out, 0, "source", "analysis/04_simulate_Rn.R, input: Climate CH2025 DAILY-GRIDDED (tasmin, tasmax, pr)")
ncatt_put(nc_out, 0, "creation_date", format(Sys.Date()))

nc_close(nc_out); nc_close(nc_in)

# ---- 7. Write one block into the NetCDF file ----------------------------------
write_block <- function(file, rn_blk, s) {
  ncols <- ncol(r_tasmin)
  n     <- nrow(rn_blk) / ncols                    # rows in this block

  arr <- array(rn_blk, dim = c(ncols, n, n_days))  # E x rows (north first) x time
  arr <- arr[, n:1, , drop = FALSE]                # the file runs south to north
  arr[is.na(arr)] <- fill

  nc <- nc_open(file, write = TRUE)
  ncvar_put(nc, "rn", arr,
            start = c(1, nrow(r_tasmin) - (s + n - 1) + 1, 1),
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

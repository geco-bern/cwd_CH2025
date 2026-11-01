# run 03_model_Rn.R first to create data/radiation_models/rg_model_rf.rds and rn_model_rf.rds
# run R/calc_rsw-cc.R first to create data/rsw-cc/meteoland_grid_ch_rswcc.nc
# applies the two-step random forest on the CH2025 daily grid:
# Tmin, Tmax, Prec, Rsw-cc, elevation -> Rg (step 1) -> Rn (step 2), written to one NetCDF file

library(here)
library(terra)
library(ncdf4)
library(ranger)
source(here("R", "nc_grid_ch.R"))                 # grid layout, create_nc(), write_block()

# ---- 1. Settings ----------------------------------------------------------
dir_grid   <- here("data-raw", "CH2025", "DAILY_GRIDDED")
model      <- "clmcom-cclm4-cccma"
period     <- "ref91-20"
block_rows <- 10                 # rows of the grid processed at once
n_threads  <- max(1, parallel::detectCores() - 1)

out_dir  <- here("data", "netradiation_sim")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("grid_ch_rn_%s_%s.nc", model, period))
rsw_file <- here("data", "rsw-cc", "meteoland_grid_ch_rswcc.nc")

nc_path <- function(var)
  file.path(dir_grid, sprintf("ogd-climate-scenarios-ch2025-grid_ch_%s_%s_%s.nc",
                              var, model, period))

# ---- 2. Inputs: climate grids, Rsw-cc, elevation, models --------------------
r_tasmin <- rast(nc_path("tasmin"))
r_tasmax <- rast(nc_path("tasmax"))
r_pr     <- rast(nc_path("pr"))
r_rsw    <- rast(rsw_file)                        # 365 days, written with the layout of R/nc_grid_ch.R
elv      <- rast(here("data", "elevation", "elevation_1k_grid.tif"))
stopifnot(compareGeom(r_tasmin, r_tasmax), compareGeom(r_tasmin, r_pr),
          compareGeom(r_tasmin, r_rsw, lyrs = FALSE), compareGeom(r_tasmin, elv),
          nlyr(r_tasmin) == n_days)

dir_model <- here("data", "radiation_models")
rg_model  <- readRDS(file.path(dir_model, "rg_model_rf.rds"))
rn_model  <- readRDS(file.path(dir_model, "rn_model_rf.rds"))

doy      <- rep(1:365, length.out = n_days)       # 365_day calendar, no Feb 29
elv_vals <- values(elv)[, 1]                      # grid order, north to south

nc_ll <- nc_open(nc_path("tasmin"))               # coordinates for the output file
lon   <- ncvar_get(nc_ll, "lon")
lat   <- ncvar_get(nc_ll, "lat")
nc_close(nc_ll)

# ---- 3. Two-step prediction for one block of rows ---------------------------
# returns a matrix: cells of the block x days (NA outside Switzerland)
process_block <- function(s) {
  n <- min(block_rows, nrow(r_tasmin) - s + 1)

  readStart(r_tasmin); readStart(r_tasmax); readStart(r_pr); readStart(r_rsw)
  tmin <- readValues(r_tasmin, row = s, nrows = n, mat = TRUE)
  tmax <- readValues(r_tasmax, row = s, nrows = n, mat = TRUE)
  pr   <- readValues(r_pr,     row = s, nrows = n, mat = TRUE)
  rsw  <- readValues(r_rsw,    row = s, nrows = n, mat = TRUE)   # cells x 365
  readStop(r_tasmin); readStop(r_tasmax); readStop(r_pr); readStop(r_rsw)

  rn_blk <- matrix(NA_real_, nrow = nrow(tmin), ncol = n_days)
  idx    <- which(!is.na(tmin[, 1]))              # cells with data
  if (length(idx) == 0) return(rn_blk)

  cell_id <- (s - 1) * ncol(r_tasmin) + idx       # cell numbers in the full grid

  # one year (365 days) at a time to keep the data frame small
  for (y in seq(1, n_days, by = 365)) {
    days <- y:(y + 364)
    d <- data.frame(tasmin = as.vector(tmin[idx, days]),
                    tasmax = as.vector(tmax[idx, days]),
                    pr_t   = log1p(as.vector(pr[idx, days])),
                    rsw_cc = as.vector(rsw[idx, doy[days]]),
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

# ---- 4. Create the Rn NetCDF file ------------------------------------------------
create_nc(out_file, "rn", "W m-2", "Daily mean net radiation", dim_time,
          lon = lon, lat = lat,
          chunks   = c(length(E_vals), block_rows, 365),   # one chunk = one block of rows x one year
          var_atts = list(standard_name = "surface_net_downward_radiative_flux",
                          method = "Two-step random forest: Rg ~ Tmin + Tmax + Prec + Rsw-cc + elevation (CH2025 DAILY_LOCAL), Rn ~ Rg + Tmin + Tmax + Prec + Rsw-cc (SwissMetNet ALT, MAG, PAY, TAE)"),
          glob_atts = list(title                 = "Daily net radiation simulated from Climate CH2025 DAILY-GRIDDED (tasmin, tasmax, pr)",
                           GWL                   = toupper(period),
                           modelchain            = toupper(model),
                           frequency             = "day",
                           qm_calibration_period = "1991-2020",
                           CORDEX_domain         = "EUR11"),
          source = "analysis/04_simulate_Rn.R")

# ---- 5. Run over all blocks with cells with data -------------------------------
starts     <- seq(1, nrow(r_tasmin), by = block_rows)
hd         <- as.matrix(!is.na(r_tasmin[[1]]), wide = TRUE)
n_by_block <- sapply(starts, function(s) sum(hd[s:min(s + block_rows - 1, nrow(hd)), ]))

for (s in starts[n_by_block > 0]) {
  message("block starting at row ", s, " (", Sys.time(), ")")
  write_block(out_file, "rn", process_block(s), s)
}
message("Done: ", out_file)

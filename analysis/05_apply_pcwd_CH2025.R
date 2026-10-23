# run 04_simulate_Rn.R first to create the gridded net radiation
# potential cumulative water deficit (PCWD) on the CH2025 daily grid, following
# geco-bern/cwd_global (ERA5Land workflow):
#   PET     = cwd::pet(Rn, tas, patm)          (Priestley-Taylor, patm from elevation)
#   snow    = pr on days with tas < temp_snow, cwd::simulate_snow() -> liquid_to_soil
#   wbal    = liquid_to_soil - PET
#   PCWD    = cwd::cwd(wbal), reset every year on doy_reset
# only cells with tas, pr, Rn and elevation are computed, all others stay NA

library(here)
library(terra)
library(ncdf4)
library(parallel)
source(here("R", "nc_grid_ch.R"))                 # grid layout, create_nc(), write_block()

# ---- 1. Settings ----------------------------------------------------------
dir_grid   <- here("data-raw", "CH2025", "DAILY_GRIDDED")
model      <- "clmcom-cclm4-cccma"
period     <- "gwl1.5"
temp_snow  <- 0                  # deg C: precipitation below this falls as snow (as cwd_global)
doy_reset  <- 1                  # deficit reset every year on 1 January
block_rows <- 10                 # rows of the grid processed at once
n_cores    <- max(1, detectCores() - 1)

rn_file  <- here("data", "netradiation_sim", sprintf("grid_ch_rn_%s_%s.nc", model, period))
out_dir  <- here("data", "pcwd")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("grid_ch_pcwd_%s_%s.nc", model, period))

nc_path <- function(var)
  file.path(dir_grid, sprintf("ogd-climate-scenarios-ch2025-grid_ch_%s_%s_%s.nc",
                              var, model, period))

# ---- 2. Inputs ---------------------------------------------------------------
r_tas <- rast(nc_path("tas"))
r_pr  <- rast(nc_path("pr"))
r_rn  <- rast(rn_file)
elv   <- rast(here("data", "elevation", "elevation_1k_grid.tif"))
stopifnot(compareGeom(r_tas, r_pr), compareGeom(r_tas, r_rn), compareGeom(r_tas, elv),
          nlyr(r_tas) == n_days, ncol(r_tas) == length(E_vals), nrow(r_tas) == length(N_vals))

elv_vals <- values(elv)[, 1]                      # grid order, north to south
nc_ll    <- nc_open(nc_path("tas"))
lon      <- ncvar_get(nc_ll, "lon")               # E x N, for the output file
lat      <- ncvar_get(nc_ll, "lat")
nc_close(nc_ll)

# ---- 3. Dates: CH2025 years 1-30 = assumed 2025-2054, no Feb 29 (as R/calc_rsw-cc.R) --
dates <- seq(as.Date("2025-01-01"), as.Date("2054-12-31"), by = "day")
dates <- dates[format(dates, "%m-%d") != "02-29"]
stopifnot(length(dates) == n_days)

# ---- 4. Air pressure from elevation (Pa), as in cwd_global R/calc_patm.R -----------
calc_patm <- function(elv)
  101325 * (1 - 0.0065 * elv / 298.15)^(9.80665 * 28.963e-3 / (8.31446262 * 0.0065))

# ---- 5. PCWD of one cell (vectors of 10950 days) -------------------------------
pcwd_cell <- function(tas, pr, pet) {
  d <- data.frame(date = dates, tas = tas,
                  rain = ifelse(tas < temp_snow, 0, pr),
                  snow = ifelse(tas < temp_snow, pr, 0))
  d <- cwd::simulate_snow(d, varnam_temp = "tas", varnam_prec = "rain", varnam_snow = "snow")
  d$wbal <- d$liquid_to_soil - pet
  cwd::cwd(d, varname_wbal = "wbal", varname_date = "date",
           thresh_drop = 0, doy_reset = doy_reset)$df$deficit
}

# workers for the cell-by-cell part (snow and cwd loop over days in R)
cl <- makeCluster(n_cores)
clusterExport(cl, c("pcwd_cell", "dates", "temp_snow", "doy_reset"))

# ---- 6. PCWD for one block of rows ------------------------------------------------
# returns a matrix: cells of the block x days (NA outside the domain)
process_block <- function(s) {
  n <- min(block_rows, nrow(r_tas) - s + 1)

  readStart(r_tas); readStart(r_pr); readStart(r_rn)
  tas <- readValues(r_tas, row = s, nrows = n, mat = TRUE)
  pr  <- readValues(r_pr,  row = s, nrows = n, mat = TRUE)
  rn  <- readValues(r_rn,  row = s, nrows = n, mat = TRUE)
  readStop(r_tas); readStop(r_pr); readStop(r_rn)

  pcwd_blk <- matrix(NA_real_, nrow = nrow(tas), ncol = n_days)
  cell_id  <- (s - 1) * ncol(r_tas) + seq_len(nrow(tas))
  # domain: all inputs present on every day
  idx <- which(rowSums(is.na(tas)) == 0 & rowSums(is.na(pr)) == 0 &
               rowSums(is.na(rn)) == 0 & !is.na(elv_vals[cell_id]))
  if (length(idx) == 0) return(pcwd_blk)

  # PET for all cells at once (cwd::pet returns mm per second -> x 86400 = mm per day)
  patm <- matrix(calc_patm(elv_vals[cell_id[idx]]), nrow = length(idx), ncol = n_days)
  pet  <- cwd::pet(netrad = rn[idx, ], tc = tas[idx, ], patm = patm) * 86400

  # snow and PCWD cell by cell, split over the workers
  part  <- ceiling(seq_along(idx) / ceiling(length(idx) / n_cores))
  parts <- lapply(split(seq_along(idx), part), function(k)
    list(tas = tas[idx[k], , drop = FALSE], pr = pr[idx[k], , drop = FALSE],
         pet = pet[k, , drop = FALSE]))
  res <- parLapply(cl, parts, function(p)
    t(vapply(seq_len(nrow(p$tas)), function(i) pcwd_cell(p$tas[i, ], p$pr[i, ], p$pet[i, ]),
             numeric(length(dates)))))

  pcwd_blk[idx, ] <- do.call(rbind, res)
  pcwd_blk
}

# ---- 7. Create the PCWD NetCDF file ---------------------------------------------
create_nc(out_file, "pcwd", "mm", "Daily potential cumulative water deficit", dim_time,
          lon = lon, lat = lat,
          chunks   = c(length(E_vals), block_rows, 365),   # one chunk = one block of rows x one year
          var_atts = list(method = sprintf("cwd::cwd(), reset on doy %d, thresh_drop 0; wbal = liquid_to_soil (cwd::simulate_snow, snow below %g deg C) - PET (cwd::pet, Priestley-Taylor with simulated Rn and patm from elevation)",
                                           doy_reset, temp_snow)),
          glob_atts = list(title      = "Daily potential cumulative water deficit from Climate CH2025 DAILY-GRIDDED (tas, pr) and simulated net radiation",
                           GWL        = toupper(period),
                           modelchain = toupper(model),
                           frequency  = "day"),
          source = "analysis/05_apply_pcwd_CH2025.R")

# ---- 8. Run over all blocks with cells with data -------------------------------
starts     <- seq(1, nrow(r_tas), by = block_rows)
hd         <- as.matrix(!is.na(r_tas[[1]]), wide = TRUE)
n_by_block <- sapply(starts, function(s) sum(hd[s:min(s + block_rows - 1, nrow(hd)), ]))

for (s in starts[n_by_block > 0]) {
  message("block starting at row ", s, " (", Sys.time(), ")")
  write_block(out_file, "pcwd", process_block(s), s)
}
message("Done: ", out_file)

stopCluster(cl)

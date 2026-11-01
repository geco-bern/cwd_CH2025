#!/usr/bin/env Rscript
#SBATCH --job-name=apply_pcwd
#SBATCH --array=1-128%8
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=apply_pcwd_%A_%a.log

# calculate potential cumulative water deficit (PCWD) on the grid using cwd::cwd() with simulated net radiation (Rn) and air pressure from elevation
#   PET  = cwd::pet(Rn, tas, patm)          (Priestley-Taylor, patm from elevation)
#   snow = pr on days with tas < temp_snow, cwd::simulate_snow() -> liquid_to_soil
#   wbal = liquid_to_soil - PET
#   PCWD = cwd::cwd(wbal), reset once the deficit drops thresh_reset mm below the event peak
# only cells with tas, pr, Rn and elevation are computed, all others stay NA
# one array task per model chain and period; finished files are skipped, so resubmit after a timeout

library(here)
library(terra)
library(ncdf4)
library(parallel)
source(here("R", "nc_grid_ch.R"))                 

# Settings
in_dir       <- here("data-raw", "CH2025", "DAILY_GRIDDED")
temp_snow    <- 0                # deg C: precipitation below this falls as snow (as cwd_global)
thresh_reset <- 50               # mm: deficit reset this far below the event peak
block_rows   <- 10               # rows of the grid processed at once
n_cores      <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# model chain and period of this array task, taken from the file names
parts  <- strsplit(sub("\\.nc$", "", list.files(in_dir, pattern = "_tas_.*\\.nc$")), "_")
task   <- parts[[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]]
model  <- task[4]
period <- task[5]

out_dir  <- here("data", "pcwd")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("grid_ch_pcwd_%s_%s.nc", model, period))
tmp_file <- paste0(out_file, ".tmp")              # renamed to out_file once complete

if (file.exists(out_file)) {
  message("Already done: ", out_file)
  quit(save = "no")
}

nc_path <- function(var)
  file.path(in_dir, sprintf("ogd-climate-scenarios-ch2025-grid_ch_%s_%s_%s.nc", var, model, period))

# Inputs
r_tas <- rast(nc_path("tas"))
r_pr  <- rast(nc_path("pr"))
r_rn  <- rast(here("data", "netradiation_sim", sprintf("grid_ch_rn_%s_%s.nc", model, period)))
elv   <- rast(here("data", "elevation", "elevation_1k_grid.tif"))
stopifnot(compareGeom(r_tas, r_pr), compareGeom(r_tas, r_rn), compareGeom(r_tas, elv),
          nlyr(r_tas) == n_days, ncol(r_tas) == length(E_vals), nrow(r_tas) == length(N_vals))

elv_vals <- values(elv)[, 1]                      # grid order, north to south
nc_ll    <- nc_open(nc_path("tas"))
lon      <- ncvar_get(nc_ll, "lon")               # E x N, for the output file
lat      <- ncvar_get(nc_ll, "lat")
nc_close(nc_ll)

# Dates: CH2025 years 1-30 = assumed 2025-2054 (365-day calendar)
dates <- seq(as.Date("2025-01-01"), as.Date("2054-12-31"), by = "day")
dates <- dates[format(dates, "%m-%d") != "02-29"]
stopifnot(length(dates) == n_days)

# Air pressure from elevation (Pa), as in cwd_global R/calc_patm.R
calc_patm <- function(elv)
  101325 * (1 - 0.0065 * elv / 298.15)^(9.80665 * 28.963e-3 / (8.31446262 * 0.0065))

# PCWD of one cell (vectors of 10950 days)
pcwd_cell <- function(tas, pr, pet) {
  d <- data.frame(date = dates, tas = tas,
                  rain = ifelse(tas < temp_snow, 0, pr),
                  snow = ifelse(tas < temp_snow, pr, 0))
  d <- cwd::simulate_snow(d, varnam_temp = "tas", varnam_prec = "rain", varnam_snow = "snow")
  d$wbal <- d$liquid_to_soil - pet
  cwd::cwd(d, varname_wbal = "wbal", varname_date = "date",
           thresh_drop = 0, thresh_reset_absolute = thresh_reset)$df$deficit
}

# PCWD for one block of rows
# returns a matrix: cells of the block x days (NA outside the domain)
process_block <- function(s) {
  n <- min(block_rows, nrow(r_tas) - s + 1)
  v <- lapply(list(tas = r_tas, pr = r_pr, rn = r_rn), function(r) {
    readStart(r)
    on.exit(readStop(r))
    readValues(r, row = s, nrows = n, mat = TRUE)
  })

  pcwd_blk <- matrix(NA_real_, nrow = nrow(v$tas), ncol = n_days)
  cell_id  <- (s - 1) * ncol(r_tas) + seq_len(nrow(v$tas))
  # domain: all inputs present on every day
  idx <- which(rowSums(is.na(v$tas)) == 0 & rowSums(is.na(v$pr)) == 0 &
               rowSums(is.na(v$rn)) == 0 & !is.na(elv_vals[cell_id]))
  if (length(idx) == 0) return(pcwd_blk)

  # PET for all cells at once (cwd::pet returns mm per second -> x 86400 = mm per day)
  patm <- matrix(calc_patm(elv_vals[cell_id[idx]]), nrow = length(idx), ncol = n_days)
  pet  <- cwd::pet(netrad = v$rn[idx, , drop = FALSE], tc = v$tas[idx, , drop = FALSE],
                   patm = patm) * 86400
  pet  <- matrix(pet, nrow = length(idx))

  # snow and PCWD cell by cell, split over the cores
  res <- mclapply(seq_along(idx), function(i)
    pcwd_cell(v$tas[idx[i], ], v$pr[idx[i], ], pet[i, ]), mc.cores = n_cores)
  stopifnot(vapply(res, is.numeric, NA))          # a failed worker returns an error object

  pcwd_blk[idx, ] <- do.call(rbind, res)
  pcwd_blk
}

# Create the PCWD NetCDF file
create_nc(tmp_file, "pcwd", "mm", "Daily potential cumulative water deficit", dim_time,
          lon = lon, lat = lat,
          chunks   = c(length(E_vals), block_rows, 365),   # one chunk = one block of rows x one year
          var_atts = list(method = sprintf("cwd::cwd(), thresh_reset_absolute %g mm, thresh_drop 0; wbal = liquid_to_soil (cwd::simulate_snow, snow below %g deg C) - PET (cwd::pet, Priestley-Taylor with simulated Rn and patm from elevation)",
                                           thresh_reset, temp_snow)),
          glob_atts = list(title      = "Daily potential cumulative water deficit from Climate CH2025 DAILY-GRIDDED (tas, pr) and simulated net radiation",
                           GWL        = toupper(period),
                           modelchain = toupper(model),
                           frequency  = "day"),
          source = "analysis/07_ubelix_apply_pcwd_CH2025.R")

# Run over all blocks (blocks without data stay NA)
for (s in seq(1, nrow(r_tas), by = block_rows)) {
  message("block starting at row ", s, " (", Sys.time(), ")")
  write_block(tmp_file, "pcwd", process_block(s), s)
}
file.rename(tmp_file, out_file)
message("Done: ", out_file)

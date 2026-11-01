# run 05_apply_pcwd_CH2025.R first to create the daily PCWD grid
# annual maximum PCWD of every cell on the CH2025 1 km grid, following geco-bern/cwd_global
# (analysis/ERA5Land/03_get_pcwd_annmax_ERA5Land_ubelix.R with R/get_cwd_annmax_byLON.R):
#   PCWD_annmax = max of the daily deficit within each year
# cwd_global loops over one RDS file per longitude, here the daily PCWD is one NetCDF file,
# so the grid is processed in blocks of rows (as in 04 and 05); cells outside the domain stay NA

library(here)
library(terra)
library(ncdf4)
source(here("R", "nc_grid_ch.R"))                 # grid layout, create_nc(), write_block()

# ---- 1. Settings ----------------------------------------------------------
model      <- "clmcom-cclm4-cccma"
period     <- "ref91-20"
block_rows <- 10                 # rows of the grid processed at once
n_years    <- n_days / 365       # 30 years, 365_day calendar (no Feb 29)

pcwd_file <- here("data", "pcwd", sprintf("grid_ch_pcwd_%s_%s.nc", model, period))
out_dir   <- here("data", "pcwd_annmax")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file  <- file.path(out_dir, sprintf("grid_ch_pcwd_annmax_%s_%s.nc", model, period))

# ---- 2. Input: daily PCWD ----------------------------------------------------
r_pcwd <- rast(pcwd_file, subds = "pcwd")
stopifnot(nlyr(r_pcwd) == n_days, ncol(r_pcwd) == length(E_vals), nrow(r_pcwd) == length(N_vals))

nc_in       <- nc_open(pcwd_file)                 # coordinates and method for the output file
lon         <- ncvar_get(nc_in, "lon")            # E x N
lat         <- ncvar_get(nc_in, "lat")
pcwd_method <- ncatt_get(nc_in, "pcwd", "method")$value
nc_close(nc_in)

# ---- 3. Years: CH2025 years 1-30 = assumed 2025-2054 (as R/calc_rsw-cc.R) ----------
year_id  <- rep(seq_len(n_years), each = 365)     # year of every day
# one time step per year, on its first day (same units and calendar as the daily files)
dim_year <- ncdim_def("time", "days since 1-1-1 00:00:00", (seq_len(n_years) - 1) * 365,
                      calendar = "365_day", longname = "time (first day of each year)")

# ---- 4. Annual maximum for one block of rows ---------------------------------------
# returns a matrix: cells of the block x years (NA outside the domain)
annmax_block <- function(s) {
  n <- min(block_rows, nrow(r_pcwd) - s + 1)

  readStart(r_pcwd)
  pcwd <- readValues(r_pcwd, row = s, nrows = n, mat = TRUE)   # cells x days
  readStop(r_pcwd)

  # cells in the domain have PCWD on every day (05_apply_pcwd_CH2025.R) and all other
  # cells are NA on every day, so max() without na.rm keeps the cells outside NA
  vapply(seq_len(n_years),
         function(y) apply(pcwd[, year_id == y, drop = FALSE], 1, max),
         numeric(nrow(pcwd)))
}

# ---- 5. Create the annual maximum NetCDF file -------------------------------------
create_nc(out_file, "pcwd_annmax", "mm", "Annual maximum potential cumulative water deficit", dim_year,
          lon = lon, lat = lat,
          chunks   = c(length(E_vals), block_rows, n_years),   # one chunk = one block of rows x all years
          var_atts = list(cell_methods = "time: maximum",
                          method       = "maximum of the daily PCWD within each year of 365 days",
                          pcwd_method  = pcwd_method),
          glob_atts = list(title      = "Annual maximum potential cumulative water deficit from Climate CH2025 DAILY-GRIDDED (tas, pr) and simulated net radiation",
                           GWL        = toupper(period),
                           modelchain = toupper(model),
                           frequency  = "year",
                           years      = "years 0001-0030 of the CH2025 365_day calendar, assumed 2025-2054"),
          source = "analysis/06_get_pcwd_annmax_CH2025.R")

# ---- 6. Run over all blocks with cells with data -------------------------------
starts     <- seq(1, nrow(r_pcwd), by = block_rows)
hd         <- as.matrix(!is.na(r_pcwd[[1]]), wide = TRUE)
n_by_block <- sapply(starts, function(s) sum(hd[s:min(s + block_rows - 1, nrow(hd)), ]))

for (s in starts[n_by_block > 0]) {
  message("block starting at row ", s, " (", Sys.time(), ")")
  write_block(out_file, "pcwd_annmax", annmax_block(s), s)
}
message("Done: ", out_file)

#!/usr/bin/env Rscript
#SBATCH --job-name=pcwd_maxlen
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker
#SBATCH --partition=cpu-invest
#SBATCH --array=1-128%32
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=pcwd_maxlen_%A_%a.log

# maximum length of a PCWD event within each year for every cell on the CH2025 1 km grid
# length = longest run of consecutive days with PCWD > 0, counted separately in each year,
# so an event that crosses 31 December is split at the end of the year
# the grid is processed in blocks of rows, cells outside the domain stay NA
# one array task per model chain and period; finished files are skipped, so resubmit after a timeout

library(here)
library(terra)
library(ncdf4)
source(here("R", "nc_grid_ch.R"))                 # grid layout, create_nc(), write_block()

# Settings
in_dir     <- here("data-raw", "CH2025", "DAILY_GRIDDED")
block_rows <- 10                 # rows of the grid processed at once
n_years    <- n_days / 365       # 30 years, 365_day calendar

# model chain and period of this array task, taken from the file names (same order as 07)
parts  <- strsplit(sub("\\.nc$", "", list.files(in_dir, pattern = "_tas_.*\\.nc$")), "_")
task   <- parts[[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]]
model  <- task[4]
period <- task[5]

pcwd_file <- here("data", "pcwd", sprintf("grid_ch_pcwd_%s_%s.nc", model, period))
out_dir   <- here("data", "pcwd_maxlen")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file  <- file.path(out_dir, sprintf("grid_ch_pcwd_maxlen_%s_%s.nc", model, period))
tmp_file  <- paste0(out_file, ".tmp")             # renamed to out_file once complete

if (file.exists(out_file)) {
  message("Already done: ", out_file)
  quit(save = "no")
}

# Input: daily PCWD
r_pcwd <- rast(pcwd_file, subds = "pcwd")
stopifnot(nlyr(r_pcwd) == n_days, ncol(r_pcwd) == length(E_vals), nrow(r_pcwd) == length(N_vals))

nc_in       <- nc_open(pcwd_file)                 # coordinates and method for the output file
lon         <- ncvar_get(nc_in, "lon")            # E x N
lat         <- ncvar_get(nc_in, "lat")
pcwd_method <- ncatt_get(nc_in, "pcwd", "method")$value
nc_close(nc_in)

# Years: CH2025 placeholder years 0001-0030 (365-day calendar, no Feb 29)
year_id  <- rep(seq_len(n_years), each = 365)     # year of every day
# one time step per year, on its first day (same units and calendar as the daily files)
dim_year <- ncdim_def("time", "days since 1-1-1 00:00:00", (seq_len(n_years) - 1) * 365,
                      calendar = "365_day", longname = "time (first day of each year)")

# Maximum event length for one block of rows
# returns a matrix: cells of the block x years (NA outside the domain)
# run counts the current consecutive days with PCWD > 0 and starts again on 1 January;
# cells outside the domain are NA on every day and stay NA
maxlen_block <- function(s) {
  readStart(r_pcwd)
  on.exit(readStop(r_pcwd))
  pcwd <- readValues(r_pcwd, row = s, nrows = min(block_rows, nrow(r_pcwd) - s + 1), mat = TRUE)
  vapply(seq_len(n_years), function(y) {
    run <- longest <- numeric(nrow(pcwd))
    for (j in which(year_id == y)) {
      run     <- (run + 1) * (pcwd[, j] > 0)
      longest <- pmax(longest, run)
    }
    longest
  }, numeric(nrow(pcwd)))
}

# Create the maximum event length NetCDF file
create_nc(tmp_file, "pcwd_maxlen", "days", "Maximum length of a potential cumulative water deficit event within the year", dim_year,
          lon = lon, lat = lat,
          chunks   = c(length(E_vals), block_rows, n_years),   # one chunk = one block of rows x all years
          var_atts = list(method      = "longest run of consecutive days with daily PCWD > 0 within each year of 365 days, events are split at the end of the year",
                          pcwd_method = pcwd_method),
          glob_atts = list(title      = "Maximum PCWD event length from Climate CH2025 DAILY-GRIDDED (tas, pr) and simulated net radiation",
                           GWL        = toupper(period),
                           modelchain = toupper(model),
                           frequency  = "year",
                           years      = "placeholder years 0001-0030 of the CH2025 365_day calendar, not calendar years"),
          source = "analysis/09_ubelix_get_pcwd_maxlen.R")

# Run over all blocks (blocks without data stay NA)
for (s in seq(1, nrow(r_pcwd), by = block_rows)) {
  message("block starting at row ", s, " (", Sys.time(), ")")
  write_block(tmp_file, "pcwd_maxlen", maxlen_block(s), s)
}
file.rename(tmp_file, out_file)
message("Done: ", out_file)

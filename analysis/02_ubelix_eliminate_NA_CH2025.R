#!/usr/bin/env Rscript
#SBATCH --job-name=eliminate_NA_CH2025
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker
#SBATCH --partition=cpu-invest
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=eliminate_NA_CH2025_%j.log

# Fill missing values in all DAILY_GRIDDED files with monthly means
# Cells that are NA on every day (outside Switzerland) are ignored.
# Original netCDF files are overwritten.
# Report stored under tables/qc_na_filled_gridded.csv

library(ncdf4)
library(parallel)

dir_inout  <- file.path("data-raw", "CH2025", "DAILY_GRIDDED")
dir_qc  <- file.path("tables", "qc_na_filled_gridded")     
out_csv <- file.path("tables", "qc_na_filled_gridded.csv")
dir.create(dir_qc, showWarnings = FALSE, recursive = TRUE)

qc_file <- function(path) file.path(dir_qc, paste0(basename(path), ".csv"))

# month of a day of the year (1..365, calendar "365_day")
month_end    <- cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
doy_to_month <- function(doy) findInterval(doy - 1, month_end) + 1

# n days starting at time index t0, as a matrix: cells x days
read_days <- function(nc, var, t0, n) {
  x <- ncvar_get(nc, var, start = c(1, 1, t0), count = c(-1, -1, n), collapse_degen = FALSE)
  matrix(x, ncol = n)
}

# check and fill one file, write its report (one row per year)
fill_file <- function(path) {
  var <- strsplit(basename(path), "_")[[1]][3]
  nc  <- nc_open(path)

  # time axis: "days since 1-1-1", calendar "365_day" -> year 1, 2, ...
  tvals <- ncvar_get(nc, "time")
  stopifnot(identical(ncatt_get(nc, "time", "calendar")$value, "365_day"), !is.unsorted(tvals))
  year  <- tvals %/% 365 + 1
  month <- doy_to_month(tvals %% 365 + 1)
  years <- unique(year)

  nE    <- nc$dim$E$len
  nN    <- nc$dim$N$len
  ncell <- nE * nN

  # pass 1: count missing values per cell and year, sums for the monthly climatology
  n_valid  <- numeric(ncell)
  n_na     <- matrix(0, ncell, length(years))
  clim_sum <- matrix(0, ncell, 12)
  clim_n   <- matrix(0, ncell, 12)

  for (k in seq_along(years)) {
    idx <- which(year == years[k])
    x   <- read_days(nc, var, idx[1], length(idx))
    ok  <- !is.na(x)

    n_valid   <- n_valid + rowSums(ok)
    n_na[, k] <- rowSums(!ok)
    for (m in 1:12) {
      j <- month[idx] == m
      clim_sum[, m] <- clim_sum[, m] + rowSums(x[, j, drop = FALSE], na.rm = TRUE)
      clim_n[, m]   <- clim_n[, m]   + rowSums(ok[, j, drop = FALSE])
    }
  }
  nc_close(nc)

  domain   <- n_valid > 0                      # cells with data (Switzerland)
  na_found <- colSums(n_na[domain, , drop = FALSE])
  replaced <- numeric(length(years))

  # pass 2: fill a copy, which replaces the original once it is complete
  if (sum(na_found) > 0) {
    clim <- clim_sum / clim_n                  # NaN where a month never has data
    clim[clim_n == 0] <- NA

    tmp <- paste0(path, ".tmp")
    stopifnot(file.copy(path, tmp, overwrite = TRUE))
    nc <- nc_open(tmp, write = TRUE)

    for (k in which(na_found > 0)) {
      idx <- which(year == years[k])
      x   <- read_days(nc, var, idx[1], length(idx))

      for (m in 1:12) {
        j     <- which(month[idx] == m)
        xm    <- x[, j, drop = FALSE]
        mmean <- rowMeans(xm, na.rm = TRUE)    # same month, same year
        mmean[is.nan(mmean)] <- clim[is.nan(mmean), m]

        gap <- is.na(xm) & domain              # recycled over the days of the month
        xm[gap] <- matrix(mmean, nrow = ncell, ncol = length(j))[gap]
        x[, j] <- xm
      }

      replaced[k] <- na_found[k] - sum(is.na(x[domain, ]))
      ncvar_put(nc, var, array(x, dim = c(nE, nN, length(idx))),
                start = c(1, 1, idx[1]), count = c(nE, nN, length(idx)))
    }
    nc_close(nc)
    file.rename(tmp, path)
  }

  write.csv(data.frame(file = basename(path), variable = var,
                       gwl = sub(".*_(.*)\\.nc$", "\\1", basename(path)),
                       year = years, n_na_found = na_found, n_replaced = replaced),
            qc_file(path), row.names = FALSE)
}

files <- list.files(dir_inout, pattern = "\\.nc$", full.names = TRUE)
todo  <- files[!file.exists(qc_file(files))]
message(length(todo), " of ", length(files), " files to process")

res <- mclapply(todo, fill_file, mc.preschedule = FALSE,
                mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
for (e in Filter(function(r) inherits(r, "try-error"), res)) message(e)

failed <- files[!file.exists(qc_file(files))]
if (length(failed) > 0) stop(length(failed), " failed, resubmit:\n", paste(failed, collapse = "\n"))

results <- do.call(rbind, lapply(qc_file(files), read.csv))
write.csv(results, out_csv, row.names = FALSE)
print(aggregate(cbind(n_na_found, n_replaced) ~ file, data = results, FUN = sum))
message("Done. Table saved to ", out_csv)

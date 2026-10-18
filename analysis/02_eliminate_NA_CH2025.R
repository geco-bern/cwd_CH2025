# ------------------------------------------------------------------
# Fill missing values (NA / empty cells) in the DAILY-GRIDDED CH2025 files
#
# For every selected file the script
#   1. reads the file year by year and counts the missing values per cell
#      and year (cells outside Switzerland, i.e. NA on every day, are ignored)
#   2. if there are missing values: replaces each one with the mean of the
#      same cell in the same month of the same year. If that whole month is
#      missing, the mean of that calendar month over all years is used.
#   3. writes the filled values back into the original .nc file (overwrite)
#   4. writes a table to tables/ with the number of replaced values per
#      file and year
#
# NOTE: the .nc files are changed in place. A second run finds nothing to
# fill and reports 0 for those files.
# ------------------------------------------------------------------

library(here)
library(ncdf4)

# ---- 1. Choices (EDIT) -----------------------------------------------
dir_grid  <- here("data-raw", "CH2025", "DAILY_GRIDDED")
model     <- "clmcom-cclm4-cccma"
parameter <- c("pr", "tas", "tasmin", "tasmax")   # variable(s) to check
gwl       <- "gwl1.5"             # period(s) to check

out_csv   <- here("tables", "qc_na_filled_gridded.csv")

# ---- 2. Files ---------------------------------------------------------
files <- expand.grid(variable = parameter, gwl = gwl, stringsAsFactors = FALSE)
files$file <- sprintf("ogd-climate-scenarios-ch2025-grid_ch_%s_%s_%s.nc",
                      files$variable, model, files$gwl)
files$path <- file.path(dir_grid, files$file)

missing <- !file.exists(files$path)
if (any(missing)) stop("File(s) not found:\n", paste(files$file[missing], collapse = "\n"))

# ---- 3. Helpers ---------------------------------------------------------
# month of a day of the year (1..365, calendar "365_day")
month_end    <- cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
doy_to_month <- function(doy) findInterval(doy - 1, month_end) + 1

# n days starting at time index t0, as a matrix: cells x days
read_days <- function(nc, var, t0, n) {
  x <- ncvar_get(nc, var, start = c(1, 1, t0), count = c(-1, -1, n), collapse_degen = FALSE)
  matrix(x, ncol = n)
}

# ---- 4. Function: check and fill one file -------------------------------
# returns a data frame: one row per year with the missing / replaced values
fill_file <- function(path, var) {
  nc <- nc_open(path, write = TRUE)
  on.exit(nc_close(nc))

  # time axis: "days since 1-1-1", calendar "365_day" -> year 1, 2, ...
  tvals <- ncvar_get(nc, "time")
  cal   <- ncatt_get(nc, "time", "calendar")$value
  stopifnot(identical(cal, "365_day"), !is.unsorted(tvals))
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

    n_valid  <- n_valid + rowSums(ok)
    n_na[, k] <- rowSums(!ok)
    for (m in 1:12) {
      j <- month[idx] == m
      clim_sum[, m] <- clim_sum[, m] + rowSums(x[, j, drop = FALSE], na.rm = TRUE)
      clim_n[, m]   <- clim_n[, m]   + rowSums(ok[, j, drop = FALSE])
    }
  }

  domain   <- n_valid > 0                      # cells with data (Switzerland)
  na_found <- colSums(n_na[domain, , drop = FALSE])
  replaced <- numeric(length(years))
  message(sprintf("  %d cells with data, %d missing values", sum(domain), sum(na_found)))

  # pass 2: fill the years with missing values and write them back
  if (sum(na_found) > 0) {
    clim <- clim_sum / clim_n                  # NaN where a month never has data
    clim[clim_n == 0] <- NA

    for (k in which(na_found > 0)) {
      idx <- which(year == years[k])
      x   <- read_days(nc, var, idx[1], length(idx))

      for (m in 1:12) {
        j    <- which(month[idx] == m)
        xm   <- x[, j, drop = FALSE]
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
    nc_sync(nc)

    if (sum(replaced) < sum(na_found))
      warning(sum(na_found) - sum(replaced), " missing values could not be filled in ",
              basename(path), " (no data for that month in any year)")
  }

  data.frame(year = years, n_na_found = na_found, n_replaced = replaced)
}

# ---- 5. Run -------------------------------------------------------------
results <- vector("list", nrow(files))

for (i in seq_len(nrow(files))) {
  message("[", i, "/", nrow(files), "] ", files$file[i])
  res <- fill_file(files$path[i], files$variable[i])
  results[[i]] <- data.frame(file     = files$file[i],
                             variable = files$variable[i],
                             gwl      = files$gwl[i],
                             res)
}
results <- do.call(rbind, results)

# ---- 6. Table -----------------------------------------------------------
# keep the rows of files that were not processed in this run
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
if (file.exists(out_csv)) {
  old     <- read.csv(out_csv, stringsAsFactors = FALSE)
  results <- rbind(old[!old$file %in% files$file, ], results)
}
write.csv(results, out_csv, row.names = FALSE)

print(aggregate(cbind(n_na_found, n_replaced) ~ file, data = results, FUN = sum))
message("Done. Table saved to ", out_csv)

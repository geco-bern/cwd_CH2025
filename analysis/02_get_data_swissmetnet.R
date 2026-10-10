# ------------------------------------------------------------------
# Download daily SwissMetNet (SMN) station data (historical files only) for all
# stations that record ALL of the parameters listed in `needed`:
#   gre000d0  global radiation, daily mean (W/m2)
#   tre200d0  air temperature 2 m, daily mean (deg C)
#   tre200dx  air temperature 2 m, daily maximum (deg C)
#   tre200dn  air temperature 2 m, daily minimum (deg C)
#   rre150d0  precipitation, daily total 6 UTC - 6 UTC following day (mm)
#
#   1. read the SMN metadata (data inventory, parameters, stations)
#   2. find the stations that record all needed parameters
#   3. download ogd-smn_[station]_d_historical.csv for each of them
#      (the file holds all daily parameters of the station, from the start of the
#       record to the end of last year)
#
# Files are saved in data-raw/swissmetnet/
# ------------------------------------------------------------------

library(here)
library(readr)
library(dplyr)

# ---- 1. Settings (EDIT) -----------------------------------------------
needed <- c("gre000d0", "tre200d0", "tre200dx", "tre200dn", "rre150d0")


out_dir  <- here("data-raw", "swissmetnet")
base_url <- "https://data.geo.admin.ch/ch.meteoschweiz.ogd-smn"

options(timeout = 600)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Metadata ---------------------------------------------------------
# uses the metadata files in out_dir; downloads them if they are not there yet
get_meta <- function(name) {
  fname <- paste0("ogd-smn_meta_", name, ".csv")
  path  <- file.path(out_dir, fname)
  if (!file.exists(path)) {
    message("Downloading ", fname)
    download.file(paste0(base_url, "/", fname), path, mode = "wb", quiet = TRUE)
  }
  path
}

# the inventory is UTF-8; the parameter and station tables use Windows-1252 (umlauts)
inventory  <- read_delim(get_meta("datainventory"), delim = ";",
                         locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
param_meta <- read_delim(get_meta("parameters"), delim = ";",
                         locale = locale(encoding = "windows-1252"), show_col_types = FALSE)
stations   <- read_delim(get_meta("stations"), delim = ";",
                         locale = locale(encoding = "windows-1252"), show_col_types = FALSE)

# what are the parameters?
p <- param_meta |> filter(parameter_shortname %in% needed)
stopifnot("a needed parameter is not in the parameter table" = nrow(p) == length(needed))
print(as.data.frame(p[c("parameter_shortname", "parameter_description_en", "parameter_unit")]),
      row.names = FALSE)

# ---- 3. Stations that record all needed parameters ----------------------------
station_ids <- inventory |>
  filter(parameter_shortname %in% needed) |>
  group_by(station_abbr) |>
  summarise(n_param = n_distinct(parameter_shortname), .groups = "drop") |>
  filter(n_param == length(needed)) |>
  pull(station_abbr)
stopifnot("some stations are missing in the station table" = all(station_ids %in% stations$station_abbr))

# station table for later use (location, elevation)
stations_sel <- stations |>
  filter(station_abbr %in% station_ids) |>
  select(station_abbr, station_name, station_canton, station_height_masl,
         station_coordinates_lv95_east, station_coordinates_lv95_north,
         station_coordinates_wgs84_lat, station_coordinates_wgs84_lon)

message(length(station_ids), " of ", nrow(stations), " stations record all ",
        length(needed), " parameters")

# ---- 4. Files to download ------------------------------------------------------
files <- data.frame(station = tolower(station_ids), stringsAsFactors = FALSE)
files$name <- paste0("ogd-smn_", files$station, "_d_historical.csv")
files$url  <- paste0(base_url, "/", files$station, "/", files$name)
files$dest <- file.path(out_dir, files$name)

# ---- 5. Download (files already in the folder are skipped) ---------------------
files$status <- NA_character_

for (i in seq_len(nrow(files))) {
  if (file.exists(files$dest[i]) && file.size(files$dest[i]) > 0) {
    files$status[i] <- "already there"
    next
  }

  message("[", i, "/", nrow(files), "] ", files$name[i])
  tryCatch(suppressWarnings(download.file(files$url[i], files$dest[i], mode = "wb", quiet = TRUE)),
           error = function(e) NULL)

  ok <- file.exists(files$dest[i]) && file.size(files$dest[i]) > 0
  if (!ok && file.exists(files$dest[i])) file.remove(files$dest[i])     # remove empty leftovers
  files$status[i] <- if (ok) "downloaded" else "FAILED"
}

# ---- 6. Report -----------------------------------------------------------------------
print(table(files$status))

failed <- files$name[files$status == "FAILED"]
if (length(failed) > 0) {
  message(length(failed), " files could not be downloaded:")
  message(paste0("  ", failed, collapse = "\n"))
}
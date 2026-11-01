#!/usr/bin/env Rscript
#SBATCH --job-name=get_data_CH2025
#SBATCH --time=96:00:00
#SBATCH --mem=2G
#SBATCH --output=get_CH2025_%j.log

# Download Climate CH2025 files via the STAC API into <working directory>/data-raw/CH2025/
# DAILY_GRIDDED/ all netCDF files
# DAILY-LOCAL/   all csv files of relevant stations to model radiation
# Files already on disk are skipped

library(jsonlite)
options(timeout = 7200)     # seconds allowed per download (large files)

stations <- tolower(c(
  "ABO", "AIG", "ALT", "BAS", "BER", "BUS", "CDF", "CGI", "CHA", "CHU",
  "CIM", "COM", "COV", "DAV", "DIS", "DOL", "ENG", "EVO", "FAH", "FRE",
  "GLA", "GRH", "GSB", "GUE", "GUT", "GVE", "INT", "KLO", "LUG", "LUZ",
  "MAG", "MLS", "MVE", "NAP", "NEU", "OTL", "PAY", "PIL", "PIO", "PLF",
  "PUY", "REH", "ROB", "ROE", "RUE", "SAE", "SAM", "SBO", "SCU", "SHA",
  "SIO", "SMA", "STG", "TAE", "ULR", "VAD", "VIS", "WAE", "WFJ", "WYN", "ZER"
))

# all file urls of a collection (the listing is split into pages)
list_urls <- function(collection) {
  url  <- paste0("https://data.geo.admin.ch/api/stac/v0.9/collections/ch.meteoschweiz.",
                 collection, "/items")
  urls <- character(0)
  while (!is.null(url)) {
    resp <- fromJSON(url, simplifyVector = FALSE)
    urls <- c(urls, unlist(lapply(resp$features, function(f) lapply(f$assets, `[[`, "href"))))
    nxt  <- Filter(function(l) l$rel == "next", resp$links)
    url  <- if (length(nxt) > 0) nxt[[1]]$href
  }
  unname(urls)
}

# download to a .part file, keep it only if the checksum matches; returns the failed urls
download_all <- function(urls, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  message(length(urls), " files for ", dir)
  failed <- character(0)
  for (url in urls) {
    dest <- file.path(dir, basename(url))
    tmp  <- paste0(dest, ".part")
    if (file.exists(dest)) next
    message("Downloading ", basename(url))
    ok <- tryCatch({
      download.file(url, tmp, mode = "wb", method = "libcurl", quiet = TRUE)
      sha <- sub(" .*", "", system2("sha256sum", shQuote(tmp), stdout = TRUE))
      any(grepl(paste0("^x-amz-meta-sha256: *", sha), curlGetHeaders(url), ignore.case = TRUE))
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok) file.rename(tmp, dest) else failed <- c(failed, url)
  }
  failed
}

gridded <- grep("\\.nc$", list_urls("ogd-climate-scenarios-ch2025-grid"), value = TRUE)
local   <- grep(paste0("ch2025_(", paste(stations, collapse = "|"), ")_.*\\.csv$"),
                list_urls("ogd-climate-scenarios-ch2025"), value = TRUE)

failed <- c(download_all(gridded, file.path("data-raw", "CH2025", "DAILY_GRIDDED")),
            download_all(local,   file.path("data-raw", "CH2025", "DAILY-LOCAL")))

if (length(failed) > 0) stop(length(failed), " failed, resubmit:\n", paste(failed, collapse = "\n"))
message("Done.")

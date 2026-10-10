#!/usr/bin/env Rscript

# Download daily surface radiation budget (SRB.d) NetCDF files from the
# CH archive item in the MeteoSwiss satellite-derived-grid STAC collection.

library(here)
library(jsonlite)

# ---- Settings ---------------------------------------------------------
collection_id <- "ch.meteoschweiz.ogd-satellite-derived-grid"
item_id <- "archive-ch"
file_prefix <- "ogd-satellite-derived-grid-archive.msg.srb.d"
archive_url <- paste0(
  "https://data.geo.admin.ch/api/stac/v1/collections/",
  collection_id, "/items/", item_id
)
output_dir <- here("data-raw", "msradiation")

# Leave TRUE to preview the matched files. Set to FALSE to download them.
dry_run <- FALSE

# ---- Read the CH archive item -----------------------------------------
archive <- jsonlite::fromJSON(archive_url, simplifyVector = FALSE)
if (!identical(archive$id, item_id)) {
  stop("Expected STAC item '", item_id, "', received '", archive$id, "'.")
}
if (is.null(archive$assets) || length(archive$assets) == 0) {
  stop("The CH archive STAC item has no assets.")
}

asset_names <- names(archive$assets)
matches <- startsWith(asset_names, file_prefix) &
  grepl("\\.nc$", asset_names, ignore.case = TRUE)
file_names <- asset_names[matches]

if (length(file_names) == 0) {
  stop("No NetCDF assets starting with '", file_prefix, "' were found.")
}

assets <- archive$assets[file_names]
urls <- vapply(file_names, function(file_name) {
  href <- assets[[file_name]]$href
  if (is.null(href) || !grepl("^https://data\\.geo\\.admin\\.ch/", href)) {
    stop("Unexpected or missing download URL for asset: ", file_name)
  }
  href
}, character(1))

if (any(basename(file_names) != file_names)) {
  stop("An archive asset name is not a plain file name; refusing to save it.")
}

destinations <- file.path(output_dir, file_names)

message("CH archive item: ", archive_url)
message("Matching SRB.d NetCDF files: ", length(file_names))
print(data.frame(
  file = file_names,
  status = ifelse(file.exists(destinations), "already exists", "to download")
), row.names = FALSE)

if (dry_run) {
  message(
    "Preview only. Set dry_run <- FALSE near the top of this script to download."
  )
} else {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Install the 'digest' package to verify downloads: install.packages('digest')")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (i in seq_along(file_names)) {
    file_name <- file_names[i]
    asset <- assets[[file_name]]
    destination <- destinations[i]
    checksum <- asset[["file:checksum"]]
    expected_sha256 <- if (
      is.character(checksum) &&
      length(checksum) == 1 &&
      grepl("^1220[[:xdigit:]]{64}$", checksum)
    ) {
      substring(checksum, 5)
    } else {
      NA_character_
    }

    if (file.exists(destination)) {
      if (!is.na(expected_sha256)) {
        actual_sha256 <- digest::digest(
          destination, algo = "sha256", file = TRUE
        )
        if (!identical(tolower(actual_sha256), tolower(expected_sha256))) {
          stop(
            "SHA-256 mismatch for existing file: ", destination,
            ". Move or remove this file before retrying."
          )
        }
      }
      message("Already downloaded: ", file_name)
      next
    }

    message("[", i, "/", length(file_names), "] Downloading ", file_name)
    status <- utils::download.file(
      urls[i], destination, mode = "wb", method = "libcurl"
    )
    if (!identical(status, 0L)) {
      stop("Download failed for: ", file_name)
    }

    if (!is.na(expected_sha256)) {
      actual_sha256 <- digest::digest(
        destination, algo = "sha256", file = TRUE
      )
      if (!identical(tolower(actual_sha256), tolower(expected_sha256))) {
        stop(
          "SHA-256 mismatch after download: ", file_name,
          ". The downloaded file may be incomplete or corrupted."
        )
      }
    } else {
      warning("No SHA-256 checksum available for ", file_name)
    }
  }

  message("SRB.d NetCDF files saved in: ", output_dir)
}


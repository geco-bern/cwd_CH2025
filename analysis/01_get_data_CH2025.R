# ------------------------------------------------------------------
# Download Climate CH2025 files (DAILY-LOCAL or DAILY-GRIDDED) via the STAC API
#
# R version of the official MeteoSwiss script download_ch2025_interactive.sh.
# All choices are made in section 1. The script then
#   1. lists all files of the collection (follows the paged results)
#   2. keeps the files that match the choices
#   3. shows which of them are already on disk and the size of the rest
#   4. asks for confirmation, downloads only the missing files and checks
#      the SHA256 checksum
#
# WARNING: gridded files are about 2 GB each.
# ------------------------------------------------------------------

library(here)
library(jsonlite)

# ---- 1. Choices (EDIT) -----------------------------------------------
mode      <- "gridded"      # "local" (station time series) or "gridded" (1 km grids)

# local only: station codes, e.g. "evo" or c("evo", "abo"), or "all"
station   <- "evo"

# variable(s): "pr", "tas", "tasmin", "tasmax"  (local also "rsds", "hurs", "sfcwind"), or "all"
parameter <- c("all") 

# gwl: "ref91-20" (historical reference), "gwl1.5", "gwl2.0",
# "gwl2.5", "gwl3.0", or "all"
gwl    <- "ref91-20"

# gridded only: model chain, e.g. "clmcom-cclm4-cccma", or "all"
model     <- "clmcom-cclm4-cccma"

# local only: ".csv" or ".zip" (gridded files are always ".nc")
format    <- ".csv"

show_size    <- TRUE        # ask the server for the file sizes (slow for many files)
check_sha256 <- TRUE        # compare the checksum after the download (needs `digest`)

# where the data are saved: <out_root>/DAILY_LOCAL/<station>/  or  <out_root>/DAILY_GRIDDED/
# NOTE: the project is in OneDrive; consider a folder outside OneDrive for large gridded data
out_root <- here("data-raw", "CH2025")

options(timeout = 7200)     # seconds allowed per download (large files)

# ---- 2. Prepare -------------------------------------------------------
stopifnot(mode %in% c("local", "gridded"))
station   <- tolower(station)
parameter <- tolower(parameter)
gwl    <- tolower(gwl)
model     <- tolower(model)
format    <- if (mode == "gridded") ".nc" else tolower(format)
if (mode == "local") stopifnot(format %in% c(".csv", ".zip"))

available_gwls <- c("ref91-20", "gwl1.5", "gwl2.0", "gwl2.5", "gwl3.0")
if (!identical(gwl, "all") && !all(gwl %in% available_gwls)) {
  stop(
    "Invalid GWL selection. Choose from: ",
    paste(c(available_gwls, "all"), collapse = ", ")
  )
}

collection_id <- if (mode == "local") "ogd-climate-scenarios-ch2025" else "ogd-climate-scenarios-ch2025-grid"
collection    <- paste0("ch.meteoschweiz.", collection_id)
api           <- "https://data.geo.admin.ch/api/stac/v0.9"     # version used by the official script

# target folder of a file: station (local) or variable (gridded), taken from the file name
target_dir <- function(fname) {
  tokens <- strsplit(sub(paste0("^", collection_id, "_"), "", fname), "_")[[1]]
  if (mode == "local") file.path(out_root, "DAILY_LOCAL", tokens[1])
  else                 file.path(out_root, "DAILY_GRIDDED")
}

# ---- 3. List all files of the collection ------------------------------
# all "href" entries anywhere in the response (like jq's  .. | .href?)
collect_hrefs <- function(x) {
  if (!is.list(x)) return(character(0))
  h <- x[["href"]]
  c(if (is.character(h)) h, unlist(lapply(x, collect_hrefs), use.names = FALSE))
}

hrefs    <- character(0)
next_url <- paste0(api, "/collections/", collection, "/items")

message("Querying STAC API for collection: ", collection)
while (!is.null(next_url)) {
  resp     <- fromJSON(next_url, simplifyVector = FALSE)
  hrefs    <- c(hrefs, collect_hrefs(resp))
  nxt      <- Filter(function(l) identical(l[["rel"]], "next"), resp$links)
  next_url <- if (length(nxt) > 0) nxt[[1]][["href"]] else NULL
}
hrefs <- unique(hrefs)

# ---- 4. Select the files that match the choices -----------------------
# A file matches if its name contains "_<value>_" or ends in "_<value>.nc".
matches <- function(fname, values, right = "_") {
  if (identical(values, "all")) return(TRUE)
  any(vapply(values, function(v) grepl(paste0("_", v, right), fname, fixed = TRUE), logical(1)))
}

keep <- vapply(basename(hrefs), function(f) {
  startsWith(f, paste0(collection_id, "_")) &&
    endsWith(f, format) &&
    matches(f, parameter) &&
  matches(f, gwl, right = format) &&
    (if (mode == "local") matches(f, station) else matches(f, model))
}, logical(1), USE.NAMES = FALSE)

urls <- hrefs[keep]
if (length(urls) == 0) stop("No matching files found. Check the choices in section 1.")

# ---- 5. What is on disk, what is missing ------------------------------
header_field <- function(url, field) {
  h <- curlGetHeaders(url)
  m <- grep(paste0("^", field, ":"), h, ignore.case = TRUE, value = TRUE)
  if (length(m) == 0) return(NA_character_)
  trimws(sub("^[^:]+:", "", tail(m, 1)))
}

size_bytes <- rep(NA_real_, length(urls))
if (show_size) {
  message("Reading file sizes from the server ...")
  size_bytes <- vapply(urls, function(u) as.numeric(header_field(u, "content-length")),
                       numeric(1), USE.NAMES = FALSE)
}

dest <- file.path(vapply(basename(urls), target_dir, character(1), USE.NAMES = FALSE),
                  basename(urls))

# a file counts as present if it exists and (when the server size is known) has that size
size_ok <- !is.na(file.size(dest)) & file.size(dest) == size_bytes
have    <- file.exists(dest) & (is.na(size_bytes) | size_ok)

status <- data.frame(file    = basename(urls),
                     on_disk = ifelse(have, "yes", "NO"),
                     size_GB = round(size_bytes / 1024^3, 2))
print(status, row.names = FALSE)

if (all(is.na(size_bytes))) {
  warning("File sizes could not be read from the server. Existing files are assumed to be complete.")
} else {
  cat(sprintf("Missing files: %d, size to download: %.2f GB\n", sum(!have),
              sum(size_bytes[!have], na.rm = TRUE) / 1024^3))
}

# ---- 6. Download the missing files -------------------------------------
if (sum(!have) > 0 &&
    !isTRUE(askYesNo(paste0("Download ", sum(!have), " missing files?")))) {
  stop("Download cancelled", call. = FALSE)
}

for (i in seq_along(urls)) {
  fname <- basename(urls[i])

  if (have[i]) {
    message("Already on disk: ", fname)
  } else {
    dir.create(dirname(dest[i]), showWarnings = FALSE, recursive = TRUE)
    message("[", i, "/", length(urls), "] Downloading ", fname)
    download.file(urls[i], dest[i], mode = "wb", method = "libcurl")
  }

  if (check_sha256) {
    remote_sha <- header_field(urls[i], "x-amz-meta-sha256")
    if (is.na(remote_sha)) {
      warning("No checksum available for ", fname)
    } else if (!identical(digest::digest(dest[i], algo = "sha256", file = TRUE), remote_sha)) {
      stop("SHA256 mismatch: ", fname, " is not intact, delete it and download again")
    }
  }
}
message("Done. Files saved below ", out_root)

# ------------------------------------------------------------------
# Step 1: fetch CH2025 DAILY-LOCAL data of one station from the STAC API
#
# Downloads the CSV files (all model chains in one file per variable and
# warming level) of one station into data-raw/CH2025/DAILY_LOCAL/<station>/
# ------------------------------------------------------------------

library(here)
library(jsonlite)

# ---- 1. Settings (EDIT) -------------------------------------------
station <- "evo"                       # station abbreviation (Evolene)
vars    <- c("pr", "tas", "rsds")      # variables to download, only the ones needed for cwd
gwls    <- c("1.5", "2.0", "2.5", "3.0")

out_dir <- here("data-raw", "CH2025", "DAILY_LOCAL", station)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Read the STAC item of the station -------------------------
item_url <- paste0("https://data.geo.admin.ch/api/stac/v1/collections/",
                   "ch.meteoschweiz.ogd-climate-scenarios-ch2025/items/", station)

item  <- fromJSON(item_url, simplifyVector = FALSE)
hrefs <- vapply(item$assets, function(a) a$href, character(1))
cat("Files (assets) listed in the STAC item:", length(hrefs), "\n")

# ---- 3. Select the files we need ----------------------------------
# file names look like ogd-climate-scenarios-ch2025_evo_pr_gwl1.5.csv
wanted <- paste0("ogd-climate-scenarios-ch2025_", station, "_",
                 rep(vars, each = length(gwls)), "_gwl", gwls, ".csv")

missing <- setdiff(wanted, basename(hrefs))
if (length(missing) > 0) stop("Not found in the STAC item:\n", paste(missing, collapse = "\n"))

sel <- hrefs[basename(hrefs) %in% wanted]

# ---- 4. Download (skips files that already exist) -----------------
options(timeout = 600)
for (u in sel) {
  dest <- file.path(out_dir, basename(u))
  if (!file.exists(dest)) {
    message("Downloading ", basename(u))
    download.file(u, dest, mode = "wb", quiet = TRUE)   # mode = "wb" is needed on Windows
  }
}

# ---- 5. Check what we have ----------------------------------------
n_files <- length(list.files(out_dir, pattern = "\\.csv$"))
stopifnot("station folder does not contain the expected number of csv files" =
            n_files == length(wanted))
message("OK: ", n_files, " files in ", out_dir)

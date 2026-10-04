# ------------------------------------------------------------------
# Step 1: fetch CH2025 DAILY-LOCAL data of one station from the STAC API
#
# Downloads the CSV files (all model chains in one file per variable and
# warming level) of one station into data-raw/CH2025/DAILY_LOCAL/<station>/
# ------------------------------------------------------------------

library(here)
library(jsonlite)

# ---- 1. Settings  -------------------------------------------
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

# ------------------------------------------------------------------
# Step 2: one data frame per station x model simulation x warming level
#
# Every folder in data-raw/CH2025/DAILY_LOCAL/ is treated as a station.
# For each station, 102 data frames (columns: date, wbal) are created:
#   CH2025_[STATION]_[MODEL]_GWL_[x]
#   e.g. CH2025_EVO_CLMCOM_CCLM4_CCCMA_GWL_1_5
# Each holds 10,950 days (0001-01-01 to 0030-12-31, 365-day calendar,
# no leap days). wbal (water balance) is empty and filled in later steps.
# ------------------------------------------------------------------

library(here)
library(dplyr)

# ---- 1. Model simulations and warming levels -----------------------
models <- c(
  "CLMCOM_CCLM4_CCCMA", "CLMCOM_CCLM4_ECEARTH", "CLMCOM_CCLM4_HADGEM",
  "CLMCOM_CCLM4_MIROC", "CLMCOM_CCLM4_MPIESM", "CNRM_ALADIN_CNRM",
  "CNRM_ALADIN_HADGEM", "CNRM_ALADIN_MPIESM", "CNRM_ALADIN_NORESM",
  "DMI_HIRHAM_ECEARTH", "DMI_HIRHAM_HADGEM", "GERICS_REMO_CNRM",
  "ICTP_REGCM_HADGEM", "IPSL_WRF_IPSL", "MOHC_HADREM_CNRM",
  "MOHC_HADREM_ECEARTH", "MOHC_HADREM_HADGEM", "MOHC_HADREM_MPIESM",
  "MOHC_HADREM_NORESM", "MPICSC_REMO1_MPIESM", "MPICSC_REMO2_MPIESM",
  "SMHI_RCA_HADGEM", "SMHI_RCA_ECEARTH", "SMHI_RCA_IPSL",
  "SMHI_RCA_MPIESM", "SMHI_RCA_NORESM"
)

# models available per warming level (GWL_3 has 24: two models are missing)
models_by_gwl <- list(
  GWL_1_5 = models,
  GWL_2   = models,
  GWL_2_5 = models,
  GWL_3   = setdiff(models, c("CLMCOM_CCLM4_MPIESM", "MPICSC_REMO1_MPIESM"))
)

# ---- 2. Stations = folders in DAILY_LOCAL ---------------------------
dir_local <- here("data-raw", "CH2025", "DAILY_LOCAL")
stations  <- list.dirs(dir_local, full.names = FALSE, recursive = FALSE)
stopifnot("no station folders found" = length(stations) > 0)
message("Stations found: ", paste(stations, collapse = ", "))

# ---- 3. Dates: 30 years x 365 days, no leap days ----------------------
dates <- seq(as.Date("0001-01-01"), as.Date("0030-12-31"), by = "day")
dates <- dates[format(dates, "%m-%d") != "02-29"]
stopifnot(length(dates) == 30 * 365)                      # 10,950

# ---- 4. Create the data frames ---------------------------------------
for (station in stations) {
  for (g in names(models_by_gwl)) {
    for (m in models_by_gwl[[g]]) {
      assign(paste0("CH2025_", toupper(station), "_", m, "_", g),
             tibble(date = dates, wbal = NA_real_),
             envir = globalenv())
    }
  }

  n <- length(ls(envir = globalenv(), pattern = paste0("^CH2025_", toupper(station), "_")))
  message(station, ": ", n, " data frames created")      # expect 102
}

View(df)
# list the data frames of one station:  ls(pattern = "^CH2025_EVO_")
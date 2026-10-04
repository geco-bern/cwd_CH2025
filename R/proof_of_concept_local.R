# ------------------------------------------------------------------
# Step 1: fetch CH2025 DAILY-LOCAL data of one station from the STAC API
#
# Downloads the CSV files (all model chains in one file per variable and
# warming level) of one station into data-raw/CH2025/DAILY_LOCAL/<station>/
# ------------------------------------------------------------------

library(here)
library(jsonlite)

# ---- Settings  -------------------------------------------
station <- "evo"                       # station abbreviation (Evolene)
vars    <- c("pr", "tas", "rsds")      # variables to download, only the ones needed for cwd
gwls    <- c("1.5", "2.0", "2.5", "3.0")

out_dir <- here("data-raw", "CH2025", "DAILY_LOCAL", station)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read the STAC item of the station -------------------------
item_url <- paste0("https://data.geo.admin.ch/api/stac/v1/collections/",
                   "ch.meteoschweiz.ogd-climate-scenarios-ch2025/items/", station)

item  <- fromJSON(item_url, simplifyVector = FALSE)
hrefs <- vapply(item$assets, function(a) a$href, character(1))
cat("Files (assets) listed in the STAC item:", length(hrefs), "\n")

# ---- Select the files we need ----------------------------------
# file names look like ogd-climate-scenarios-ch2025_evo_pr_gwl1.5.csv
wanted <- paste0("ogd-climate-scenarios-ch2025_", station, "_",
                 rep(vars, each = length(gwls)), "_gwl", gwls, ".csv")

missing <- setdiff(wanted, basename(hrefs))
if (length(missing) > 0) stop("Not found in the STAC item:\n", paste(missing, collapse = "\n"))

sel <- hrefs[basename(hrefs) %in% wanted]

# ---- Download (skips files that already exist) -----------------
options(timeout = 600)
for (u in sel) {
  dest <- file.path(out_dir, basename(u))
  if (!file.exists(dest)) {
    message("Downloading ", basename(u))
    download.file(u, dest, mode = "wb", quiet = TRUE)   # mode = "wb" is needed on Windows
  }
}

# ---- Check what we have ----------------------------------------
n_files <- length(list.files(out_dir, pattern = "\\.csv$"))
stopifnot("station folder does not contain the expected number of csv files" =
            n_files == length(wanted))
message("OK: ", n_files, " files in ", out_dir)

# ------------------------------------------------------------------
# Step 2: one data frame per station x model simulation x warming level
#
# Every folder in data-raw/CH2025/DAILY_LOCAL/ is treated as a station.
# For each station, 102 data frames are created:
#   CH2025_[STATION]_[MODEL]_GWL_[x]
#   e.g. CH2025_EVO_CLMCOM_CCLM4_CCCMA_GWL_1_5
# Each holds 10,950 days (0001-01-01 to 0030-12-31, 365-day calendar,
# no leap days) with the columns date, pr, tas, rsds, rain, snow, snow_pool,
# liquid_to_soil, pet, wbal.
# Apart from date, all columns are empty placeholders filled in later steps.
# ------------------------------------------------------------------

library(here)
library(dplyr)

# ---- Model simulations and warming levels -----------------------
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

# ---- Stations = folders in DAILY_LOCAL ---------------------------
dir_local <- here("data-raw", "CH2025", "DAILY_LOCAL")
stations  <- list.dirs(dir_local, full.names = FALSE, recursive = FALSE)
stopifnot("no station folders found" = length(stations) > 0)
message("Stations found: ", paste(stations, collapse = ", "))

# ---- Dates: 30 years x 365 days, no leap days ----------------------
dates <- seq(as.Date("0001-01-01"), as.Date("0030-12-31"), by = "day")
dates <- dates[format(dates, "%m-%d") != "02-29"]
stopifnot(length(dates) == 30 * 365)                      # 10,950

# ---- Create the data frames ---------------------------------------
for (station in stations) {
  for (g in names(models_by_gwl)) {
    for (m in models_by_gwl[[g]]) {
      assign(paste0("CH2025_", toupper(station), "_", m, "_", g),
             tibble(date           = dates,
                    pr             = NA_real_,   # precipitation (mm/day)
                    tas            = NA_real_,   # mean temperature (deg C)
                    rsds           = NA_real_,   # incoming shortwave radiation
                    rain           = NA_real_,   # precipitation falling as rain (mm/day)
                    snow           = NA_real_,   # precipitation falling as snow, water equivalent (mm/day)
                    snow_pool      = NA_real_,   # snow mass (mm water equivalent)
                    liquid_to_soil = NA_real_,   # rain + snowmelt (mm/day)
                    pet            = NA_real_,   # potential evapotranspiration (mm/day)
                    wbal           = NA_real_),  # liquid_to_soil - pet
             envir = globalenv())
    }
  }

  n <- length(ls(envir = globalenv(), pattern = paste0("^CH2025_", toupper(station), "_")))
  message(station, ": ", n, " data frames created")      # expect 102
}

# ------------------------------------------------------------------
# Step 3: fill a variable (pr, tas or rsds) into the data frames from step 2
#
# For every station folder in data-raw/CH2025/DAILY_LOCAL/ and every warming
# level, the file
#   ogd-climate-scenarios-ch2025_[station]_[var]_gwl[x].csv
# is read once. Each model column is copied into the column `var` of the
# matching data frame CH2025_[STATION]_[MODEL]_GWL_[x].
#
# ------------------------------------------------------------------

# ---- Settings ---------------------------------------------------
# warming level as written in the object names and in the file names
gwl_in_file <- c(GWL_1_5 = "1.5", GWL_2 = "2.0", GWL_2_5 = "2.5", GWL_3 = "3.0")

# ---- Reader for the CH2025 csv files ----------------------------
# The files start with metadata lines; the table starts at the line "DATE;..."
read_ch2025 <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  hdr   <- which(startsWith(lines, "DATE;"))[1]
  out   <- read.delim(text = lines[hdr:length(lines)], sep = ";",
                      check.names = FALSE, na.strings = c("", "NA", "NaN"))
  out$DATE <- as.Date(out$DATE)      # years 0001-0030 are placeholders (365_day calendar)
  attr(out, "unit") <- sub("^UNIT;", "", lines[startsWith(lines, "UNIT;")][1])
  out
}

# ---- Function: fill one variable into all data frames -------------
fill_variable <- function(var) {
  for (station in stations) {
    for (g in names(gwl_in_file)) {

      f <- file.path(dir_local, station,
                     paste0("ogd-climate-scenarios-ch2025_", station, "_", var,
                            "_gwl", gwl_in_file[[g]], ".csv"))
      if (!file.exists(f)) stop("File not found: ", f, " (run step 1 for this station)")

      x      <- read_ch2025(f)
      models <- setdiff(names(x), "DATE")      # model names in the file use "-"

      for (m in models) {
        nm <- paste0("CH2025_", toupper(station), "_", gsub("-", "_", m), "_", g)
        if (!exists(nm, envir = globalenv())) stop(nm, " does not exist, run step 2 first")

        d <- get(nm, envir = globalenv())
        if (length(d$date) != length(x$DATE) || !all(d$date == x$DATE))
          stop(nm, ": dates in the data frame and in the file differ")

        d[[var]] <- x[[m]]
        assign(nm, d, envir = globalenv())
      }

      message(station, " | ", var, " | GWL ", gwl_in_file[[g]], ": ",
              length(models), " models filled (unit: ", attr(x, "unit"), ")")

      # models with missing days in this file
      n_na <- colSums(is.na(x[models]))
      if (any(n_na > 0))
        message("   missing days: ",
                paste0(names(n_na)[n_na > 0], " (", n_na[n_na > 0], ")", collapse = ", "))
    }
  }
}

# ---- Run for pr, tas and rsds ----------------------------------------------------------
fill_variable("pr")
fill_variable("tas")
fill_variable("rsds")

# ------------------------------------------------------------------
# Step 4: rain/snow split and snow routine
#
# For every data frame CH2025_[STATION]_[MODEL]_GWL_[x] (all stations):
#   1. split pr into rain and snow with a temperature threshold (tas)
#   2. run cwd::simulate_snow() -> snow_pool and liquid_to_soil
# The results go into the existing columns rain, snow, snow_pool and
# liquid_to_soil of the same data frame (no new objects are created).
#
# The cwd package comes from GitHub:  remotes::install_github("geco-bern/cwd")
# ------------------------------------------------------------------

# ---- Setting ----------------------------------------------------
temp_snow <- 0      # deg C: precipitation below this temperature falls as snow

# ---- Loop over all data frames ----------------------------------
done    <- 0
skipped <- character(0)

for (nm in ls(pattern = "^CH2025_", envir = globalenv())) {
  d <- get(nm, envir = globalenv())

  if (all(is.na(d$tas)) || all(is.na(d$pr)))
    stop(nm, ": pr or tas not filled yet, run step 3 first")

  # cwd::simulate_snow() cannot handle missing values: skip these series
  if (anyNA(d$tas) || anyNA(d$pr)) {
    skipped <- c(skipped, nm)
    next
  }

  # rain and snow (mm/day)
  d$rain <- ifelse(d$tas < temp_snow, 0, d$pr)
  d$snow <- ifelse(d$tas < temp_snow, d$pr, 0)

  # snow routine: fills snow_pool and liquid_to_soil
  d <- cwd::simulate_snow(d, varnam_temp = "tas", varnam_prec = "rain", varnam_snow = "snow")

  assign(nm, d, envir = globalenv())
  done <- done + 1
}

# ---- Report -------------------------------------------------------
message(done, " data frames processed")
if (length(skipped) > 0) {
  message(length(skipped), " data frames skipped (missing values in pr or tas):")
  message(paste0("  ", skipped, collapse = "\n"))
}


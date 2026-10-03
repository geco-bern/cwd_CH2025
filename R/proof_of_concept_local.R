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
# Step 2: one data frame per warming level (GWL) x model chain (GCM-RCM)
#
# Each data frame has the columns `date` and `wbal` (water balance, still empty).
# Later steps add liquid_to_soil and pet and then fill wbal.
#
# Names: ch2025_[station]_[varname]_[model]_gwl[GWL]
#   e.g. ch2025_evo_wbal_CLMCOM_CCLM4_CCCMA_gwl1_5
# ("-" and "." are replaced by "_" so the names work without backticks)
# ------------------------------------------------------------------

library(dplyr)

# ---- 1. Settings ---------------------------------------------------
varname    <- "wbal"                 # part of the object names
n_expected <- c("1.5" = 26, "2.0" = 26, "2.5" = 26, "3.0" = 24)   # model chains per GWL

dir_station <- here("data-raw", "CH2025", "DAILY_LOCAL", station)

# path of one downloaded file, e.g. ogd-climate-scenarios-ch2025_evo_pr_gwl1.5.csv
file_for <- function(var, gwl) {
  file.path(dir_station,
            paste0("ogd-climate-scenarios-ch2025_", station, "_", var, "_gwl", gwl, ".csv"))
}

# ---- 2. Reader for the CH2025 csv files ----------------------------
# metadata lines first; the table starts at the line "DATE;..."
read_ch2025 <- function(path) {
  stopifnot("file not found, run step 1 first" = file.exists(path))
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  hdr   <- which(startsWith(lines, "DATE;"))[1]
  out   <- read.delim(text = lines[hdr:length(lines)], sep = ";",
                      check.names = FALSE, na.strings = c("", "NA", "NaN"))
  out$DATE <- as.Date(out$DATE)      # years 0001-0030 are placeholders (365_day calendar)
  out
}

# ---- 3. Naming -----------------------------------------------------
safe     <- function(x) gsub("[-.]", "_", x)
obj_name <- function(model, g) {
  paste0("ch2025_", station, "_", varname, "_", safe(model), "_gwl", safe(g))
}

# ---- 4. Create the data frames -------------------------------------
template_names <- character(0)       # names of all created data frames
models_by_gwl  <- list()

for (g in gwls) {
  pr     <- read_ch2025(file_for("pr", g))
  models <- setdiff(names(pr), "DATE")

  if (length(models) != n_expected[[g]])
    stop("GWL", g, ": found ", length(models), " model chains, expected ", n_expected[[g]])

  # tas and rsds are needed later: same models and same dates as pr?
  for (v in c("tas", "rsds")) {
    x <- read_ch2025(file_for(v, g))
    if (!setequal(names(x), names(pr)) || !identical(x$DATE, pr$DATE))
      stop("GWL", g, ": ", v, " has different models or dates than pr")
  }

  for (m in models) {
    nm <- obj_name(m, g)
    assign(nm, tibble::tibble(date = pr$DATE, wbal = NA_real_), envir = globalenv())
    template_names <- c(template_names, nm)
  }
  models_by_gwl[[g]] <- models
}

# ---- 5. Checks -----------------------------------------------------
stopifnot(length(template_names) == sum(n_expected),     # 102
          !anyDuplicated(template_names))

print(lengths(models_by_gwl))                            # expect 26, 26, 26, 24
cat("Models missing at GWL3.0 (compared with GWL2.5):",
    paste(setdiff(models_by_gwl[["2.5"]], models_by_gwl[["3.0"]]), collapse = ", "), "\n")
cat("Rows per data frame:",
    unique(vapply(template_names, function(n) nrow(get(n)), integer(1))), "\n")   # expect 10950

head(template_names)

# ---- How to use the data frames in the next steps -------------------
# one data frame:      ch2025_evo_wbal_CLMCOM_CCLM4_CCCMA_gwl1_5
# loop over all:       for (nm in template_names) { df <- get(nm); ...; assign(nm, df) }
# list them in R:      ls(pattern = "^ch2025_evo_wbal_")
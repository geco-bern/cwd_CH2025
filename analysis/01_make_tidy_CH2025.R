#!/usr/bin/env Rscript

# Convert the four Swiss CH2025 daily gridded NetCDF files into chunked,
# tidy time-series files. Output is written by variable to data/tidy/CH2025.

library(here)
library(map2tidy)

input_dir <- here("data-raw", "CH2025", "DAILY_GRIDDED")
output_dir <- here("data", "tidy", "CH2025", "GWL_1.5")

variables <- c("pr", "tas", "tasmax", "tasmin")
model <- "clmcom-cclm4-cccma"
warming_level <- "1.5"

nc_files <- setNames(
  file.path(
    input_dir,
    paste0(
      "ogd-climate-scenarios-ch2025-grid_ch_",
      variables,
      "_", model, "_gwl", warming_level, ".nc"
    )
  ),
  variables
)

missing_files <- nc_files[!file.exists(nc_files)]
if (length(missing_files) > 0) {
  stop(
    "Could not find the expected CH2025 NetCDF file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

for (variable in variables) {
  variable_output_dir <- file.path(output_dir, variable)
  dir.create(variable_output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Converting ", variable, " from ", basename(nc_files[[variable]]))

  map2tidy(
    nclist = nc_files[[variable]],
    varnam = variable,
    lonnam = "lon",
    latnam = "lat",
    timenam = "time",
    do_chunks = TRUE,
    outdir = variable_output_dir,
    fileprefix = paste0("CH2025_GWL_", warming_level, "_", variable),
    ncores = 1,
    overwrite = FALSE
  )
}

message("Finished converting the four CH2025 variables to tidy chunks in: ",
        output_dir)

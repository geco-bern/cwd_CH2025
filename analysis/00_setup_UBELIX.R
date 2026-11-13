# installs all R packages needed to run the complete analysis

use_pkgs <- c(
  "jsonlite",
  "ncdf4",
  "terra",
  "meteoland",
  "here",
  "dplyr",
  "ggplot2",
  "ranger",
  "remotes"
)

new_pkgs <- use_pkgs[!(use_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

# cwd R package for snow, PET and PCWD from GECO Bern
if (!requireNamespace("cwd", quietly = TRUE)) remotes::install_github("geco-bern/cwd")

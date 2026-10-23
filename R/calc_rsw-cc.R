## Calculate daily potential shortwave radiation (Rsw-cc) for the CH2025 stations (1) and grid (2)

library(here)
library(meteoland)
library(terra)
library(ncdf4)

# 1. Calculate daily potential shortwave radiation (Rsw-cc) for the CH2025 stations 

dir_in  <- here("data-raw", "CH2025", "DAILY_LOCAL")
dir_out <- here("data", "rsw-cc")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

stations <- c("abo", "aig", "alt", "bas", "ber", "bus", "cdf", "cgi", "cha", "chu",
              "cim", "com", "cov", "dav", "dis", "dol", "eng", "evo", "fah", "fre",
              "gla", "grh", "gsb", "gue", "gut", "gve", "int", "klo", "lug", "luz",
              "mag", "mls", "mve", "nap", "neu", "otl", "pay", "pil", "pio", "plf",
              "puy", "reh", "rob", "roe", "rue", "sae", "sam", "sbo", "scu", "sha",
              "sio", "sma", "stg", "tae", "ulr", "vad", "vis", "wae", "wfj", "wyn",
              "zer")

# latitude from the CH2025 station metadata
meta <- read.csv(file.path(dir_in, "ogd-climate-scenarios-ch2025_meta_stations.csv"),
                 sep = ";", fileEncoding = "UTF-8-BOM")
meta <- meta[match(toupper(stations), meta$station_abbr), ]

# daily potential shortwave radiation (rsw-cc) in W m-2 for 365 days, flat surface
d     <- seq(as.Date("2025-01-01"), by = "day", length.out = 365)
J     <- radiation_dateStringToJulianDays(format(d, "%Y-%m-%d"))
delta <- vapply(J, radiation_solarDeclination, numeric(1))
sc    <- vapply(J, radiation_solarConstant, numeric(1))

calc_rsw_cc <- function(lat_deg) {
  pot <- mapply(function(s, dl)
    radiation_potentialRadiation(s, lat_deg * pi / 180, 0, 0, dl), sc, delta)

  pot * 1e6 / 86400                                  # MJ m-2 day-1 -> W m-2
}

# dates: CH2025 placeholder years 0001-0030 = assumed 2025-2054, no Feb 29
d_assumed <- seq(as.Date("2025-01-01"), as.Date("2054-12-31"), by = "day")
d_assumed <- d_assumed[format(d_assumed, "%m-%d") != "02-29"]
d_ph      <- sprintf("%04d-%s", as.integer(format(d_assumed, "%Y")) - 2024,
                     format(d_assumed, "%m-%d"))
doy       <- rep(1:365, times = 30)

# one csv per station: meteoland_[station]_rswcc.csv
for (i in seq_along(stations)) {
  rsw <- calc_rsw_cc(meta$station_coordinates_wgs84_lat[i])

  write.csv(data.frame(date         = d_ph,
                       date_assumed = d_assumed,
                       doy          = doy,
                       rsw_cc       = round(rsw[doy], 2)),
            file.path(dir_out, paste0("meteoland_", stations[i], "_rswcc.csv")),
            row.names = FALSE)
}

# 2. rsw-cc on the CH2025 1 km grid: meteoland_grid_ch_rswcc.nc

source(here("R", "nc_grid_ch.R"))                    # grid layout, create_nc()

file_grid <- here("data-raw", "CH2025", "DAILY_GRIDDED",
                  "ogd-climate-scenarios-ch2025-grid_ch_tasmin_clmcom-cclm4-cccma_gwl1.5.nc")
nc   <- nc_open(file_grid)
lon  <- ncvar_get(nc, "lon")                         # E x N (N south to north)
lat  <- ncvar_get(nc, "lat")
with_data <- which(!is.na(ncvar_get(nc, "tasmin", start = c(1, 1, 1), count = c(-1, -1, 1))))  # cells with data
nc_close(nc)

rsw_grid <- matrix(fill, nrow = length(lat), ncol = 365)   # cells (file order) x doy
rsw_grid[with_data, ] <- t(vapply(lat[with_data], calc_rsw_cc, numeric(365)))

file_rsw <- file.path(dir_out, "meteoland_grid_ch_rswcc.nc")
create_nc(file_rsw, "rsw_cc", "W m-2",
          "Daily potential shortwave radiation on a flat surface (Rsw-cc)", dim_doy,
          lon = lon, lat = lat,
          var_atts  = list(method = "meteoland::radiation_potentialRadiation, slope 0, aspect 0, latitude of each cell"),
          glob_atts = list(title = "Daily potential shortwave radiation (Rsw-cc) on the CH2025 1 km grid"),
          source    = "R/calc_rsw-cc.R")
nc <- nc_open(file_rsw, write = TRUE)
ncvar_put(nc, "rsw_cc", array(rsw_grid, dim = c(dim(lat), 365)))
nc_close(nc)
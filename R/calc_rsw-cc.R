library(here)
library(meteoland)

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
calc_rsw_cc <- function(lat_deg) {
  d     <- seq(as.Date("2025-01-01"), by = "day", length.out = 365)
  J     <- radiation_dateStringToJulianDays(format(d, "%Y-%m-%d"))
  delta <- vapply(J, radiation_solarDeclination, numeric(1))
  sc    <- vapply(J, radiation_solarConstant, numeric(1))

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
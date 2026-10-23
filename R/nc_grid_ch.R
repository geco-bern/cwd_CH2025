# Layout and writer for NetCDF files on the CH2025 1 km grid
# (same layout as in analysis/04_simulate_Rn.R)
# needs: library(ncdf4), library(terra)

# ---- Layout -------------------------------------------------------------------
# regular 1 km grid in CH1903+ / LV95, coordinates = cell centres
res_m  <- 1000
E_vals <- seq(2474500, 2843500, by = res_m)       # easting,  370 cells, west to east
N_vals <- seq(1064500, 1303500, by = res_m)       # northing, 240 cells, south to north
n_days <- 30 * 365                                # 30 years, 365_day calendar (no Feb 29)
fill   <- -999.99                                 # value outside the domain

dim_E    <- ncdim_def("E", "meters", E_vals, longname = "swiss easting (lv95)")
dim_N    <- ncdim_def("N", "meters", N_vals, longname = "swiss northing (lv95)")
dim_time <- ncdim_def("time", "days since 1-1-1 00:00:00", 0:(n_days - 1),
                      unlim = TRUE, calendar = "365_day", longname = "time")
dim_doy  <- ncdim_def("doy", "day of year", 1:365,
                      longname = "day of year (365_day calendar, no Feb 29)")

# CF grid mapping of LV95 (EPSG:2056)
crs_atts <- list(grid_mapping_name                 = "oblique_mercator",
                 longitude_of_projection_origin    = 7.43958333,
                 latitude_of_projection_origin     = 46.9524056,
                 azimuth_of_central_line           = 90,
                 scale_factor_at_projection_origin = 1,
                 false_easting                     = 2600000,
                 false_northing                    = 1200000,
                 reference_ellipsoid_name          = "Bessel 1841",
                 semi_major_axis                   = 6377397.155,
                 inverse_flattening                = 299.1528128,
                 horizontal_datum_name             = "CH1903+",
                 prime_meridian_name               = "Greenwich",
                 crs_wkt                           = crs("EPSG:2056"),
                 spatial_ref                       = crs("EPSG:2056"))

# ---- Create an empty NetCDF file with this layout ---------------------------------
# lon, lat: E x N matrices (N south to north); dim3: dim_time or another 3rd dimension
create_nc <- function(file, var, units, longname, dim3, lon, lat, chunks = NA,
                      var_atts = list(), glob_atts = list(), source = NA) {
  var_lon <- ncvar_def("lon", "degrees_east",  list(dim_E, dim_N), prec = "float",
                       longname = "longitude coordinate", compression = 4)
  var_lat <- ncvar_def("lat", "degrees_north", list(dim_E, dim_N), prec = "float",
                       longname = "latitude coordinate", compression = 4)
  var_crs <- ncvar_def("swiss_lv95_coordinates", "", list(), prec = "integer")
  var_out <- ncvar_def(var, units, list(dim_E, dim_N, dim3), missval = fill,
                       prec = "float", longname = longname,
                       chunksizes = chunks, compression = 4)

  nc <- nc_create(file, list(var_lon, var_lat, var_crs, var_out), force_v4 = TRUE)
  ncvar_put(nc, "lon", lon)
  ncvar_put(nc, "lat", lat)

  ncatt_put(nc, "E", "standard_name", "projection_x_coordinate"); ncatt_put(nc, "E", "axis", "X")
  ncatt_put(nc, "N", "standard_name", "projection_y_coordinate"); ncatt_put(nc, "N", "axis", "Y")
  if (dim3$name == "time") {
    ncatt_put(nc, "time", "standard_name", "time"); ncatt_put(nc, "time", "axis", "T")
  }
  ncatt_put(nc, "lon", "standard_name", "longitude")
  ncatt_put(nc, "lat", "standard_name", "latitude")
  for (a in names(crs_atts)) ncatt_put(nc, "swiss_lv95_coordinates", a, crs_atts[[a]])

  ncatt_put(nc, var, "grid_mapping", "swiss_lv95_coordinates")
  ncatt_put(nc, var, "coordinates",  "lat lon")
  for (a in names(var_atts)) ncatt_put(nc, var, a, var_atts[[a]])

  ncatt_put(nc, 0, "Conventions",   "CF-1.8")
  ncatt_put(nc, 0, "spatial",       "regular 1km grid")
  if (!is.na(source)) ncatt_put(nc, 0, "source", source)
  ncatt_put(nc, 0, "creation_date", format(Sys.Date()))
  for (a in names(glob_atts)) ncatt_put(nc, 0, a, glob_atts[[a]])
  nc_close(nc)
}

# ---- Write one block of rows (cells x days, terra order: north first) ---------------
write_block <- function(file, var, blk, s) {
  ncols <- length(E_vals)
  n     <- nrow(blk) / ncols                       # rows in this block

  arr <- array(blk, dim = c(ncols, n, ncol(blk)))  # E x rows (north first) x time
  arr <- arr[, n:1, , drop = FALSE]                # the file runs south to north
  arr[is.na(arr)] <- fill

  nc <- nc_open(file, write = TRUE)
  ncvar_put(nc, var, arr,
            start = c(1, length(N_vals) - (s + n - 1) + 1, 1),
            count = c(ncols, n, ncol(blk)))
  nc_close(nc)
}

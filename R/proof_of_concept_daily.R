library(here)
library(terra)

# ---- Settings ----------------------------------------------------
dir_grid <- here("data-raw", "CH2025", "DAILY_GRIDDED")
model    <- "clmcom-cclm4-cccma"

# ---- Settings (period) ----------------------------------------
period <- "gwl1.5"

# helper: path of one file
nc_path <- function(var, period)
  file.path(dir_grid, sprintf("ogd-climate-scenarios-ch2025-grid_ch_%s_%s_%s.nc",
                              var, model, period))

stopifnot(file.exists(nc_path("tas", period)), file.exists(nc_path("pr", period)))

# ---- tas: unit and value range --------------------------------
r_tas <- rast(nc_path("tas", period))
print(r_tas)
print(unique(units(r_tas)))                                       
print(global(r_tas[[1]], c("min", "max"), na.rm = TRUE))          

# ---- Elevation grid (from prepare_elevation.R) ----------------
elv <- rast(here("data", "elevation", "elevation_1k_grid.tif"))
print(compareGeom(r_tas, elv))                                    

# ---- Cells with data (Switzerland) vs. cells without ----------
has_data <- !is.na(r_tas[[1]])
print(table(values(has_data)))                                    

# ---- Block design --------------------------------------------------
block_rows <- 10                                   # rows per block
starts     <- seq(1, nrow(r_tas), by = block_rows) # first row of each block

# number of cells with data per block
hd         <- as.matrix(has_data, wide = TRUE)     # 240 x 370 (TRUE = Switzerland)
n_by_block <- sapply(starts, function(s) sum(hd[s:min(s + block_rows - 1, nrow(hd)), ]))
print(n_by_block)

# ---- Test: read the block with the most data -----------------------
b <- which.max(n_by_block)
s <- starts[b]
n <- min(block_rows, nrow(r_tas) - s + 1)

r_pr <- rast(nc_path("pr", period))

r_pr <- rast(nc_path("pr", period))

readStart(r_tas); readStart(r_pr)          # open the files for reading

t_read <- system.time({
  tas_blk <- readValues(r_tas, row = s, nrows = n, mat = TRUE)   # cells x days
  pr_blk  <- readValues(r_pr,  row = s, nrows = n, mat = TRUE)
})

readStop(r_tas); readStop(r_pr)            # close them again
print(t_read)

print(dim(tas_blk))                      # expect: (n * 370) x 10950
keep <- !is.na(tas_blk[, 1])             # cells with data in this block
print(sum(keep))
print(c(NA_tas = sum(is.na(tas_blk[keep, ])), NA_pr = sum(is.na(pr_blk[keep, ]))))
print(range(pr_blk[keep, ], na.rm = TRUE))

# ---- Settings -------------------------------------------------------
temp_snow <- 1          # deg C: precipitation below this falls as snow

# ---- Test on the first 50 cells with data in the block -----------------
cells_test <- which(keep)[1:50]

liquid_test <- matrix(NA_real_, nrow = length(cells_test), ncol = ncol(tas_blk))

t_snow <- system.time({
  for (k in seq_along(cells_test)) {
    i     <- cells_test[k]
    tas_i <- tas_blk[i, ]
    pr_i  <- pr_blk[i, ]

    d <- data.frame(tas  = tas_i,
                    rain = ifelse(tas_i < temp_snow, 0, pr_i),
                    snow = ifelse(tas_i < temp_snow, pr_i, 0))

    d <- cwd::simulate_snow(d, varnam_temp = "tas", varnam_prec = "rain", varnam_snow = "snow")
    liquid_test[k, ] <- d$liquid_to_soil
  }
})
print(t_snow)

print(summary(rowMeans(pr_blk[cells_test, ]) * 365))   # mean annual precipitation (mm/yr)
print(summary(rowMeans(liquid_test) * 365))             # mean annual liquid_to_soil (mm/yr)


# ---- Air pressure from elevation (Pa) -------------------------------------
# same formula as the internal cwd:::calc_patm (standard atmosphere)
calc_patm <- function(elv)
  101325 * (1 - 0.0065 * elv / 298.15)^(9.80665 * 28.963e-3 / (8.31446262 * 0.0065))

# ---- Elevation of the test cells ----------------------------------------------
cell_id <- (s - 1) * ncol(r_tas) + cells_test      # cell numbers in the full grid
z_test  <- values(elv)[cell_id]
print(summary(z_test))
stopifnot(!anyNA(z_test))                           # every cell with data needs an elevation

# ---- PET (cwd::pet returns mm per SECOND -> x 86400 = mm per day) ----------------
tas_test  <- tas_blk[cells_test, ]                  # cells x days
patm_test <- matrix(calc_patm(z_test), nrow = length(z_test), ncol = ncol(tas_test))

pet_test  <- cwd::pet(netrad = 200, tc = tas_test, patm = patm_test) * 86400

# ---- Water balance (mm per day) ---------------------------------------------------
wbal_test <- liquid_test - pet_test

# ---- Sanity checks -------------------------------------------------------------------
print(cwd::pet(200, 15, 101325) * 86400)            # one value: 15 deg C at sea level (expect about 5.5 mm/day)
print(summary(rowMeans(pet_test)  * 365))           # mean annual PET (mm/yr)
print(summary(rowMeans(wbal_test) * 365))           # mean annual water balance (mm/yr)
print(cor(z_test, rowMeans(pet_test)))              # should be negative (higher = cooler = less PET)


library(ncdf4)

# ---- Output location ----------------------------------------------------
out_dir  <- here("data", "wbal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("wbal_%s_%s.nc", model, period))

# ---- Take axes and projection info from the input file -----------------------
nc_in <- nc_open(nc_path("tas", period))

cal <- ncatt_get(nc_in, "time", "calendar")           # calendar attribute, if present
dim_E    <- ncdim_def("E", "meters", nc_in$dim$E$vals)
dim_N    <- ncdim_def("N", "meters", nc_in$dim$N$vals)
dim_time <- ncdim_def("time", nc_in$dim$time$units, nc_in$dim$time$vals,
                      calendar = if (cal$hasatt) cal$value else NA)

crs_att  <- ncatt_get(nc_in, "swiss_lv95_coordinates")  # projection attributes
nc_close(nc_in)

# ---- Define the variables -----------------------------------------------------
fill <- -9999                                         # value for cells outside Switzerland

var_wbal <- ncvar_def("wbal", units = "mm day-1",
                      dim = list(dim_E, dim_N, dim_time),
                      missval = fill, prec = "float",
                      longname = "Daily water balance (liquid_to_soil - PET)",
                      chunksizes = c(370, block_rows, 365),   # one chunk = one block of rows
                      compression = 4)
var_crs  <- ncvar_def("swiss_lv95_coordinates", "", list(), prec = "integer")

# ---- Create the file and copy the projection attributes ---------------------------
nc_out <- nc_create(out_file, list(var_wbal, var_crs), force_v4 = TRUE)

for (a in names(crs_att)) ncatt_put(nc_out, "swiss_lv95_coordinates", a, crs_att[[a]])
ncatt_put(nc_out, "wbal", "grid_mapping", "swiss_lv95_coordinates")
ncatt_put(nc_out, "wbal", "comment",
          "PET: Priestley-Taylor with constant net radiation 200 W m-2 (proof of concept)")

nc_close(nc_out)
print(file.exists(out_file))

# ---- Settings ---------------------------------------------------------------
netrad <- 200                                     # W m-2, constant (proof of concept)
elv_vals <- values(elv)[, 1]                      # elevation of every cell, in grid order

# ---- Function 1: water balance for one block of rows ------------------------
# returns a matrix: cells of the block x days (NA outside Switzerland)
process_block <- function(s) {
  n <- min(block_rows, nrow(r_tas) - s + 1)

  readStart(r_tas); readStart(r_pr)
  tas_blk <- readValues(r_tas, row = s, nrows = n, mat = TRUE)
  pr_blk  <- readValues(r_pr,  row = s, nrows = n, mat = TRUE)
  readStop(r_tas); readStop(r_pr)

  wbal_blk <- matrix(NA_real_, nrow = nrow(tas_blk), ncol = ncol(tas_blk))
  idx <- which(!is.na(tas_blk[, 1]))              # cells with data
  if (length(idx) == 0) return(wbal_blk)

  # rain/snow split and snow routine, cell by cell
  liquid <- matrix(NA_real_, nrow = length(idx), ncol = ncol(tas_blk))
  for (k in seq_along(idx)) {
    tas_i <- tas_blk[idx[k], ]
    pr_i  <- pr_blk[idx[k], ]
    d <- data.frame(tas  = tas_i,
                    rain = ifelse(tas_i < temp_snow, 0, pr_i),
                    snow = ifelse(tas_i < temp_snow, pr_i, 0))
    d <- cwd::simulate_snow(d, varnam_temp = "tas", varnam_prec = "rain", varnam_snow = "snow")
    liquid[k, ] <- d$liquid_to_soil
  }

  # PET (mm per day) from tas, elevation and constant net radiation
  z     <- elv_vals[(s - 1) * ncol(r_tas) + idx]
  stopifnot(!anyNA(z))
  tas_k <- tas_blk[idx, , drop = FALSE]
  patm  <- matrix(calc_patm(z), nrow = length(z), ncol = ncol(tas_k))
  pet   <- cwd::pet(netrad = netrad, tc = tas_k, patm = patm) * 86400

  wbal_blk[idx, ] <- liquid - pet
  wbal_blk
}

# ---- Function 2: write one block into the NetCDF file ---------------------------
write_block <- function(file, wbal_blk, s) {
  ncols <- ncol(r_tas)
  n     <- nrow(wbal_blk) / ncols                  # rows in this block
  ntime <- ncol(wbal_blk)

  arr <- array(wbal_blk, dim = c(ncols, n, ntime)) # E x rows (north first) x time
  arr <- arr[, n:1, , drop = FALSE]                # the file runs south to north
  arr[is.na(arr)] <- fill

  nc <- nc_open(file, write = TRUE)
  ncvar_put(nc, "wbal", arr,
            start = c(1, nrow(r_tas) - (s + n - 1) + 1, 1),
            count = c(ncols, n, ntime))
  nc_close(nc)
}

# ---- Test on the biggest block ----------------------------------------------------
print(system.time({
  wbal_blk <- process_block(s)
  write_block(out_file, wbal_blk, s)
}))

# ---- Read back with terra and compare -------------------------------------------------
r_out <- rast(out_file)
print(r_out)
print(compareGeom(r_out, r_tas))                   # same grid as the input?
print(crs(r_out, describe = TRUE)$code)            # expect "2056"

chk  <- values(r_out, row = s, nrows = n, mat = TRUE)   # same rows, read from the file
land <- !is.na(wbal_blk[, 1])
print(all.equal(chk[land, ], wbal_blk[land, ], tolerance = 1e-4))   # file = what we computed?
print(all.equal(wbal_blk[cells_test, ], wbal_test, check.attributes = FALSE))  # = earlier 50-cell test?


#!/usr/bin/env Rscript

library(here)
library(terra)

target_date <- as.Date("2004-06-30")
input_file <- here(
	"data-raw", "msradiation",
	"ogd-satellite-derived-grid-archive.msg.srb.d_ch02.lonlat_20040101000000_20041231000000.nc"
)
output_file <- here("figures", "srb_2004-06-30.png")

if (!file.exists(input_file)) {
	stop("SRB file not found: ", input_file, call. = FALSE)
}

srb <- terra::rast(input_file)
dates <- as.Date(terra::time(srb))
layer_index <- match(target_date, dates)
if (is.na(layer_index)) {
	stop("No SRB layer is available for ", target_date, ".", call. = FALSE)
}
daily_srb <- srb[[layer_index]]

plot_srb <- function() {
	terra::plot(
		daily_srb,
		col = grDevices::hcl.colors(100, palette = "Blue-Red 3"),
		main = paste("Surface Radiative Budget (SRB) |", target_date),
		xlab = "Longitude",
		ylab = "Latitude",
		plg = list(title = "SRB (W/m2)")
	)
}

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
grDevices::png(output_file, width = 1200, height = 800, res = 120)
plot_srb()
grDevices::dev.off()

plot_srb()
message("Map saved to: ", output_file)

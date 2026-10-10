#!/usr/bin/env Rscript

# Compare CH2025 daily climate predictors with MeteoSwiss surface radiative budget.

library(here)
library(terra)

# ---- Settings ---------------------------------------------------------
ch2025_dir <- here("data-raw", "CH2025", "DAILY_GRIDDED")
msradiation_dir <- here("data-raw", "msradiation")

model_chain <- "clmcom-cclm4-cccma"
reference_period <- "ref91-20"
test_start_year <- 2017L
sampled_cells <- 300L
max_train_rows <- 50000L
max_test_rows <- 25000L
knn_k <- 10L
random_forest_trees <- 500L
random_seed <- 20261010L

predictor_names <- c("pr", "tas", "tasmax", "tasmin")
target_name <- "netradiation"
variable_files <- setNames(
	file.path(
		ch2025_dir,
		paste0(
			"ogd-climate-scenarios-ch2025-grid_ch_", predictor_names, "_",
			model_chain, "_", reference_period, ".nc"
		)
	),
	predictor_names
)

missing_files <- variable_files[!file.exists(variable_files)]
if (length(missing_files) > 0) {
	stop(
		"Required CH2025 files are not available yet:\n",
		paste(missing_files, collapse = "\n"),
		call. = FALSE
	)
}

if (!requireNamespace("randomForest", quietly = TRUE) ||
		!requireNamespace("FNN", quietly = TRUE)) {
	stop(
		"Install the modeling packages first: install.packages(c('randomForest', 'FNN'))",
		call. = FALSE
	)
}

as_dates <- function(x, label) {
	if (inherits(x, "Date")) return(x)
	if (inherits(x, "POSIXt")) return(as.Date(x, tz = "UTC"))
	dates <- suppressWarnings(as.Date(x))
	if (length(dates) == 0 || anyNA(dates)) {
		stop("Could not read valid daily timestamps for ", label, ".", call. = FALSE)
	}
	dates
}

get_raster_dates <- function(raster, label) {
	dates <- as_dates(terra::time(raster), label)
	if (length(dates) != terra::nlyr(raster) || anyDuplicated(dates)) {
		stop("Invalid or duplicate timestamps in ", label, ".", call. = FALSE)
	}
	dates
}

extract_daily_values <- function(raster, raster_dates, dates, points, label) {
	layer_index <- match(dates, raster_dates)
	if (anyNA(layer_index)) {
		stop("Some requested dates are missing from ", label, ".", call. = FALSE)
	}
	values <- terra::extract(raster[[layer_index]], points, ID = FALSE)
	values <- as.matrix(values)
	if (ncol(values) != length(dates)) {
		stop("Unexpected extracted layer count for ", label, ".", call. = FALSE)
	}
	as.vector(t(values))
}

message("Reading CH2025 reference grids...")
climate_rasters <- lapply(variable_files, terra::rast)
climate_dates <- Map(get_raster_dates, climate_rasters, names(climate_rasters))
common_climate_dates <- Reduce(intersect, climate_dates)
if (length(common_climate_dates) == 0) {
	stop("The CH2025 predictor files have no common dates.", call. = FALSE)
}

rad_files <- sort(list.files(
	msradiation_dir,
	pattern = "^ogd-satellite-derived-grid-archive\\.msg\\.srb\\.d_ch02\\.lonlat_.*\\.nc$",
	full.names = TRUE
))
if (length(rad_files) == 0) {
	stop("No MeteoSwiss SRB NetCDF files found in ", msradiation_dir, ".", call. = FALSE)
}

rad_rasters <- lapply(rad_files, terra::rast)
rad_dates <- Map(get_raster_dates, rad_rasters, basename(rad_files))
all_rad_dates <- do.call(c, rad_dates)
if (anyDuplicated(all_rad_dates)) {
	stop("The MeteoSwiss files contain overlapping dates.", call. = FALSE)
}

dates <- sort(intersect(common_climate_dates, all_rad_dates))
if (length(dates) == 0) {
	stop("No dates overlap between CH2025 and MeteoSwiss SRB.", call. = FALSE)
}
if (!any(as.integer(format(dates, "%Y")) >= test_start_year) ||
		!any(as.integer(format(dates, "%Y")) < test_start_year)) {
	stop("The selected dates do not provide both training and test years.", call. = FALSE)
}

message(
	"Overlapping dates: ", min(dates), " to ", max(dates),
	" (", length(dates), " days); test period starts ", test_start_year, "."
)

# Sample CH2025 cell centers that also fall inside the satellite grid extent.
set.seed(random_seed)
template <- climate_rasters[["tas"]]
rad_template <- rad_rasters[[1]]
candidate_points <- terra::spatSample(
	template[[1]],
	size = max(5000L, sampled_cells * 10L),
	method = "random",
	na.rm = TRUE,
	as.points = TRUE,
	values = FALSE
)
candidate_rad <- terra::project(candidate_points, terra::crs(rad_template))
candidate_xy <- terra::crds(candidate_rad)
rad_extent <- terra::ext(rad_template)
inside_rad_extent <- candidate_xy[, 1] >= terra::xmin(rad_extent) &
	candidate_xy[, 1] <= terra::xmax(rad_extent) &
	candidate_xy[, 2] >= terra::ymin(rad_extent) &
	candidate_xy[, 2] <= terra::ymax(rad_extent)
candidate_points <- candidate_points[which(inside_rad_extent), ]
if (nrow(candidate_points) < sampled_cells) {
	stop(
		"Only ", nrow(candidate_points), " candidate cells overlap the SRB extent. ",
		"Reduce sampled_cells or check the source grids.",
		call. = FALSE
	)
}
selected <- sample(seq_len(nrow(candidate_points)), sampled_cells)
sample_points <- candidate_points[selected, ]
cell_ids <- seq_len(nrow(sample_points))

message("Extracting predictors and SRB at ", length(cell_ids), " sampled cells...")
data_values <- list(
	date = rep(dates, each = length(cell_ids)),
	cell_id = rep(cell_ids, times = length(dates))
)
for (name in predictor_names) {
	data_values[[name]] <- extract_daily_values(
		climate_rasters[[name]], climate_dates[[name]], dates,
		sample_points, name
	)
}

target_values <- rep(NA_real_, length(dates) * length(cell_ids))
for (i in seq_along(rad_rasters)) {
	file_dates <- rad_dates[[i]]
	date_positions <- match(file_dates, dates)
	keep <- which(!is.na(date_positions))
	if (length(keep) == 0) next

	radiation_points <- terra::project(sample_points, terra::crs(rad_rasters[[i]]))
	values <- terra::extract(
		rad_rasters[[i]][[keep]], radiation_points, ID = FALSE
	)
	values <- as.matrix(values)
	for (j in seq_along(keep)) {
		row_index <- (date_positions[keep[j]] - 1L) * length(cell_ids) + cell_ids
		target_values[row_index] <- values[, j]
	}
}
data_values[[target_name]] <- target_values
model_data <- as.data.frame(data_values)
model_data <- model_data[stats::complete.cases(model_data), ]
model_data <- model_data[is.finite(model_data[[target_name]]), ]
if (nrow(model_data) < 1000L) {
	stop("Fewer than 1,000 complete cell-days remain after matching.", call. = FALSE)
}

# Hold out later years to evaluate temporal generalization.
test_rows <- as.integer(format(model_data$date, "%Y")) >= test_start_year
train_data <- model_data[!test_rows, c(predictor_names, target_name)]
test_data <- model_data[test_rows, c("date", "cell_id", predictor_names, target_name)]
if (nrow(train_data) > max_train_rows) {
	train_data <- train_data[sample.int(nrow(train_data), max_train_rows), ]
}
if (nrow(test_data) > max_test_rows) {
	test_data <- test_data[sample.int(nrow(test_data), max_test_rows), ]
}
if (nrow(train_data) < 100L || nrow(test_data) < 100L) {
	stop("Too few complete rows in the training or test period.", call. = FALSE)
}

correlations <- stats::cor(
	model_data[, c(predictor_names, target_name)],
	use = "pairwise.complete.obs"
)
correlation_file <- here("tables", "netradiation_predictor_correlations.csv")
dir.create(dirname(correlation_file), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(correlations, correlation_file, quote = TRUE)

observed <- test_data[[target_name]]
predictor_frame <- test_data[, predictor_names, drop = FALSE]
predictions <- list()

linear_model <- stats::lm(
	stats::reformulate(predictor_names, response = target_name),
	data = train_data
)
predictions[["Linear regression"]] <- as.numeric(
	stats::predict(linear_model, newdata = predictor_frame)
)

train_x <- as.matrix(train_data[, predictor_names, drop = FALSE])
test_x <- as.matrix(predictor_frame)
train_center <- colMeans(train_x)
train_scale <- apply(train_x, 2, stats::sd)
train_scale[!is.finite(train_scale) | train_scale == 0] <- 1
train_x <- sweep(sweep(train_x, 2, train_center, "-"), 2, train_scale, "/")
test_x <- sweep(sweep(test_x, 2, train_center, "-"), 2, train_scale, "/")
knn_fit <- FNN::knn.reg(
	train = train_x,
	test = test_x,
	y = train_data[[target_name]],
	k = knn_k
)
predictions[["KNN regression"]] <- as.numeric(knn_fit$pred)

set.seed(random_seed)
forest_model <- randomForest::randomForest(
	stats::reformulate(predictor_names, response = target_name),
	data = train_data,
	ntree = random_forest_trees,
	importance = TRUE
)
predictions[["Random forest"]] <- as.numeric(
	stats::predict(forest_model, newdata = predictor_frame)
)

score_predictions <- function(prediction, observed) {
	error <- prediction - observed
	c(
		MAE = mean(abs(error)),
		RMSE = sqrt(mean(error^2)),
		BIAS = mean(error),
		R2 = 1 - sum(error^2) / sum((observed - mean(observed))^2)
	)
}
metrics <- as.data.frame(do.call(rbind, lapply(predictions, score_predictions, observed)))
metrics$model <- rownames(metrics)
rownames(metrics) <- NULL
metrics <- metrics[, c("model", "MAE", "RMSE", "BIAS", "R2")]
metrics <- metrics[order(metrics$RMSE), ]

metrics_file <- here("tables", "netradiation_model_metrics.csv")
utils::write.csv(metrics, metrics_file, row.names = FALSE)
prediction_file <- here("tables", "netradiation_model_test_predictions.csv")
prediction_table <- data.frame(
	date = test_data$date,
	cell_id = test_data$cell_id,
	observed = observed,
	as.data.frame(predictions, check.names = FALSE)
)
utils::write.csv(prediction_table, prediction_file, row.names = FALSE)

plot_file <- here("figures", "netradiation_model_diagnostics.png")
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
grDevices::png(plot_file, width = 1400, height = 1100, res = 150)
graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot_limits <- range(c(observed, unlist(predictions)), finite = TRUE)
for (name in names(predictions)) {
	graphics::plot(
		observed, predictions[[name]],
		xlim = plot_limits, ylim = plot_limits,
		xlab = "Observed SRB (W/m2)", ylab = "Predicted SRB (W/m2)",
		main = name, pch = 16, cex = 0.35, col = grDevices::adjustcolor("#236B68", 0.35)
	)
	graphics::abline(0, 1, col = "#B44936", lwd = 2)
}
grDevices::dev.off()

message("Complete cell-days: ", nrow(model_data))
message("Training rows: ", nrow(train_data), "; test rows: ", nrow(test_data))
print(metrics, row.names = FALSE)
message("Metrics: ", metrics_file)
message("Correlations: ", correlation_file)
message("Test predictions: ", prediction_file)
message("Diagnostic plot: ", plot_file)

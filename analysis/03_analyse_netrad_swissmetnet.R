# ------------------------------------------------------------------
# Predict daily radiation from daily temperature and precipitation
#
# Target (EDIT `target`):
#   "net_radiation": gre000d0 - osr000d0 + oli000d0 - olo000d0   (W/m2, daily mean)
#                    -> only possible at the few stations that record all components
#   "gre000d0":      global radiation, daily mean (W/m2)
#                    -> available at 123 stations
# Predictors: tas (tre200d0), tasmax (tre200dx), tasmin (tre200dn), pr (rre150d0)
#   feature set "base"        : tas, tasmax, tasmin, log(1 + pr)
#   feature set "plus_season" : base + daily temperature range (tasmax - tasmin)
#                               + day of year
# Models: lm (reference), KNN, random forest, GAM (recommendation)
# Validation: leave-station-out cross-validation (grouped by station),
#             i.e. every model predicts stations it has never seen.
# Metrics: MAE, RMSE, Bias (mean of prediction - observation), R2 (1 - SSE/SST)
#
# ------------------------------------------------------------------

library(here)
library(readr)
library(dplyr)
library(ggplot2)

# ---- 1. Settings (EDIT) -----------------------------------------------
target    <- "net_radiation"            # "net_radiation" or "gre000d0"
smn_dir   <- here("data-raw", "swissmetnet")

knn_k     <- 15                         # number of neighbours
rf_trees  <- 300                        # trees of the random forest
n_folds   <- 10                         # max. number of station groups (<= number of stations)
max_train <- 200000                     # max. training rows per fold (speed, mainly for gre000d0)
seed      <- 42

feature_sets <- list(
  base        = c("tas", "tasmax", "tasmin", "pr_log"),
  plus_season = c("tas", "tasmax", "tasmin", "pr_log", "dtr", "doy_sin", "doy_cos")
)

set.seed(seed)

# SwissMetNet parameter names
v_pred <- c(tas = "tre200d0", tasmax = "tre200dx", tasmin = "tre200dn", pr = "rre150d0")
v_rad  <- c(gre = "gre000d0", osr = "osr000d0", oli = "oli000d0", olo = "olo000d0")
v_target <- if (target == "net_radiation") v_rad else v_rad["gre"]
v_needed <- unname(c(v_pred, v_target))

# ---- 2. Which stations can be used? -----------------------------------------
inventory <- read_delim(file.path(smn_dir, "ogd-smn_meta_datainventory.csv"), delim = ";",
                        locale = locale(encoding = "UTF-8"), show_col_types = FALSE)
stations_meta <- read_delim(file.path(smn_dir, "ogd-smn_meta_stations.csv"), delim = ";",
                            locale = locale(encoding = "windows-1252"), show_col_types = FALSE)

stations_use <- inventory |>
  filter(parameter_shortname %in% v_needed) |>
  group_by(station_abbr) |>
  summarise(n = n_distinct(parameter_shortname), .groups = "drop") |>
  filter(n == length(v_needed)) |>
  pull(station_abbr)

st_info <- stations_meta |>
  filter(station_abbr %in% stations_use) |>
  select(station_abbr, station_name, station_height_masl)
message(length(stations_use), " stations record all needed parameters (", target, ")")
if (length(stations_use) <= 10) print(as.data.frame(st_info), row.names = FALSE)
message("Elevation of the training stations: ",
        min(st_info$station_height_masl), " to ", max(st_info$station_height_masl), " m a.s.l.")
if (length(stations_use) < 10)
  warning("Only ", length(stations_use), " stations: the results cannot be transferred to ",
          "elevations or climates outside the range of these stations.")

# ---- 3. Read the daily files ------------------------------------------------------
read_smn <- function(station) {
  f <- file.path(smn_dir, paste0("ogd-smn_", tolower(station), "_d_historical.csv"))
  if (!file.exists(f)) stop("File not found: ", f, " (run download_swissmetnet.R)")

  d <- read_delim(f, delim = ";", na = c("", "NA", "-"),
                  col_types = cols(.default = col_character()), show_col_types = FALSE)

  miss <- setdiff(c("reference_timestamp", v_needed), names(d))
  if (length(miss) > 0) stop(station, ": columns missing in the file: ", paste(miss, collapse = ", "))

  # date, usually written as "01.01.1981 00:00"
  date <- as.Date(strptime(d$reference_timestamp, "%d.%m.%Y %H:%M", tz = "UTC"))
  if (all(is.na(date))) date <- as.Date(substr(d$reference_timestamp, 1, 10))
  if (all(is.na(date))) stop(station, ": cannot read the dates, first value: ", d$reference_timestamp[1])

  num <- function(v) if (v %in% names(d)) as.numeric(d[[v]]) else rep(NA_real_, nrow(d))

  tibble(station = station, date = date,
         tas = num(v_pred["tas"]), tasmax = num(v_pred["tasmax"]),
         tasmin = num(v_pred["tasmin"]), pr = num(v_pred["pr"]),
         gre = num(v_rad["gre"]), osr = num(v_rad["osr"]),
         oli = num(v_rad["oli"]), olo = num(v_rad["olo"]))
}

dat <- bind_rows(lapply(stations_use, read_smn))

# ---- 4. Target and predictors ---------------------------------------------------------
dat <- dat |>
  mutate(y       = if (target == "net_radiation") gre - osr + oli - olo else gre,
         dtr     = tasmax - tasmin,
         pr_log  = log1p(pr),
         doy     = as.integer(format(date, "%j")),
         doy_sin = sin(2 * pi * doy / 365.25),
         doy_cos = cos(2 * pi * doy / 365.25))

n_all <- nrow(dat)
dat <- dat |>
  filter(!is.na(y), !is.na(tas), !is.na(tasmax), !is.na(tasmin), !is.na(pr),
         dtr >= 0, pr >= 0) |>
  as.data.frame()
message(nrow(dat), " of ", n_all, " days are complete (target and all predictors)")

# ---- 5. Relationships: correlations with the target -----------------------------------
cat("\nPearson correlation with the target (", target, "):\n", sep = "")
print(round(cor(dat[c("tas", "tasmax", "tasmin", "pr_log", "dtr", "doy_sin", "doy_cos")], dat$y)[, 1], 2))

# ---- 6. Models ---------------------------------------------------------------------------
gam_formula <- function(f) {
  smooth <- setdiff(f, c("doy_sin", "doy_cos"))
  terms  <- paste0("s(", smooth, ", k = 8)")
  if (any(c("doy_sin", "doy_cos") %in% f)) terms <- c(terms, "s(doy, bs = 'cc', k = 8)")  # cyclic over the year
  as.formula(paste("y ~", paste(terms, collapse = " + ")))
}

models <- list(
  lm = function(train, test, f) {
    as.numeric(predict(lm(reformulate(f, "y"), data = train), newdata = test))
  },
  knn = function(train, test, f) {                        # predictors are standardised
    mu  <- colMeans(train[f]);  sdv <- apply(train[f], 2, sd)
    FNN::knn.reg(scale(as.matrix(train[f]), center = mu, scale = sdv),
                 scale(as.matrix(test[f]),  center = mu, scale = sdv),
                 y = train$y, k = knn_k)$pred
  },
  rf = function(train, test, f) {
    m <- ranger::ranger(reformulate(f, "y"), data = train, num.trees = rf_trees,
                        min.node.size = 5, seed = seed)
    predict(m, data = test)$predictions
  },
  gam = function(train, test, f) {                        # recommendation
    m <- mgcv::bam(gam_formula(f), data = train, method = "fREML", discrete = TRUE)
    as.numeric(predict(m, newdata = test))
  }
)

# ---- 7. Leave-station-out cross-validation ----------------------------------------------------
st      <- unique(dat$station)
n_folds <- min(n_folds, length(st))
fold_of <- setNames(sample(rep(seq_len(n_folds), length.out = length(st))), st)
dat$fold <- unname(fold_of[dat$station])
message("Cross-validation: ", n_folds, " folds, grouped by station")

res <- list()
for (fs in names(feature_sets)) {
  f <- feature_sets[[fs]]
  for (mod in names(models)) {
    message("  ", fs, " | ", mod)
    pred <- rep(NA_real_, nrow(dat))

    for (i in seq_len(n_folds)) {
      test_idx <- which(dat$fold == i)
      train    <- dat[dat$fold != i, ]
      if (nrow(train) > max_train) train <- train[sample(nrow(train), max_train), ]
      pred[test_idx] <- models[[mod]](train, dat[test_idx, ], f)
    }

    res[[paste(fs, mod)]] <- data.frame(station = dat$station, date = dat$date,
                                        obs = dat$y, pred = pred,
                                        model = mod, feature_set = fs)
  }
}
pred_all <- bind_rows(res)

# ---- 8. Metrics -------------------------------------------------------------------------------------
summarise_metrics <- function(x) {
  x |> summarise(n    = n(),
                 MAE  = mean(abs(pred - obs)),
                 RMSE = sqrt(mean((pred - obs)^2)),
                 Bias = mean(pred - obs),
                 R2   = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
                 .groups = "drop")
}

cat("\n=== Overall (out-of-station predictions; ", target, ", W/m2) ===\n", sep = "")
tab_overall <- pred_all |> group_by(feature_set, model) |> summarise_metrics()
print(as.data.frame(tab_overall |> mutate(across(c(MAE, RMSE, Bias, R2), ~ round(.x, 2)))), row.names = FALSE)

cat("\n=== By season ===\n")
tab_season <- pred_all |>
  mutate(month  = as.integer(format(date, "%m")),
         season = case_when(month %in% c(12, 1, 2) ~ "DJF", month %in% 3:5 ~ "MAM",
                            month %in% 6:8 ~ "JJA", TRUE ~ "SON")) |>
  group_by(feature_set, model, season) |>
  summarise_metrics()
print(as.data.frame(tab_season |> mutate(across(c(MAE, RMSE, Bias, R2), ~ round(.x, 2)))), row.names = FALSE)

cat("\n=== By station ===\n")
tab_station <- pred_all |> group_by(feature_set, model, station) |> summarise_metrics()
print(as.data.frame(tab_station |> mutate(across(c(MAE, RMSE, Bias, R2), ~ round(.x, 2)))), row.names = FALSE)

# ---- 9. Plot: observed vs predicted ---------------------------------------------------------------------
p <- ggplot(pred_all, aes(obs, pred)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, colour = "firebrick") +
  facet_grid(feature_set ~ model) +
  coord_equal() +
  labs(x = paste0("Observed ", target, " (W/m2)"), y = "Predicted (W/m2)") +
  theme_bw()
print(p)

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", paste0("model_comparison_", target, ".png")), p, width = 11, height = 6, dpi = 200)

# ---- 10. Relationships in the GAM (fitted on all data) --------------------------------------------------
gam_full <- mgcv::bam(gam_formula(feature_sets$plus_season), data = dat, method = "fREML", discrete = TRUE)
print(summary(gam_full)$r.sq)
plot(gam_full, pages = 1, scale = 0, shade = TRUE)    # effect of every predictor on the target
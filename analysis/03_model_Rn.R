# run calc_rsw-cc.R first to create the rsw-cc files in data/rsw-cc

# three step approach: first train a random forest on the CH2025 daily local data, then use the trained model to predict net radiation on four swissmetnet stations

# Step 1: Train a random forest on the CH2025 daily local data

library(here)
library(dplyr)
library(ggplot2)

# Settings
dir_in_local <- here("data-raw", "CH2025", "DAILY_LOCAL")
dir_in_rsw   <- here("data", "rsw-cc")
vars         <- c("pr", "tasmin", "tasmax", "rsds")
max_stations <- Inf          
n_per_group  <- 40
n_trees <- 100; seed <- 1
run_cv_rg    <- FALSE        # TRUE to rerun the step 1 CV (slow, only needed when step 1 changes)
n_threads    <- max(1, parallel::detectCores() - 1)
x_vars  <- c("pr_t", "tasmin", "tasmax", "rsw_cc", "elev")
fml     <- reformulate(x_vars, response = "rsds")
dir_model  <- here("data", "radiation_models")
file_train <- file.path(dir_model, "rg_training.rds")
for (d in c("data", "tables", "figures")) dir.create(here(d), showWarnings = FALSE)
dir.create(dir_model, showWarnings = FALSE)

# Reader for one CH2025 csv (table starts at the line "DATE;")
read_ch2025 <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  meta  <- function(key) sub(paste0("^", key, ";"), "", lines[startsWith(lines, paste0(key, ";"))][1])
  out   <- read.delim(text = lines[which(startsWith(lines, "DATE;"))[1]:length(lines)],
                      sep = ";", check.names = FALSE, na.strings = c("", "NA", "NaN"))
  out$DATE <- as.Date(out$DATE)
  attr(out, "elev") <- as.numeric(meta("STATION_HEIGHT"))
  attr(out, "lat")  <- as.numeric(meta("LATITUDE"))
  out
}

# File inventory
files <- list.files(dir_in_local, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
mm    <- regmatches(basename(files), regexec(
  "^ogd-climate-scenarios-ch2025_([A-Za-z0-9]+)_([a-z]+)_(.+)\\.csv$", basename(files)))
parts <- do.call(rbind, mm[lengths(mm) == 4])
inv   <- data.frame(path = files[lengths(mm) == 4], station = toupper(parts[, 2]),
                    var = parts[, 3], period = parts[, 4])
inv   <- inv[inv$var %in% vars, ]
inv   <- inv[inv$station %in% head(sort(unique(inv$station)), max_stations), ]

# Training subset: n_per_group random complete days per station x period x chain
build_training <- function(inv) {
  set.seed(seed)
  sets <- unique(inv[, c("station", "period")])
  tr <- bind_rows(lapply(seq_len(nrow(sets)), function(i) {
    sub <- inv[inv$station == sets$station[i] & inv$period == sets$period[i], ]
    x   <- lapply(setNames(vars, vars), function(v) read_ch2025(sub$path[sub$var == v]))
    chains <- Reduce(intersect, lapply(x, function(d) setdiff(names(d), "DATE")))
    bind_rows(lapply(chains, function(ch) {
      d   <- data.frame(pr = x$pr[[ch]], tasmin = x$tasmin[[ch]],
                        tasmax = x$tasmax[[ch]], rsds = x$rsds[[ch]])
      idx <- which(complete.cases(d))
      idx <- idx[sample.int(length(idx), min(n_per_group, length(idx)))]
      data.frame(station = sets$station[i], elev = attr(x$pr, "elev"),
                 lat = attr(x$pr, "lat"), period = sets$period[i], chain = ch,
                 date = x$pr$DATE[idx], d[idx, ])
    }))
  }))
  tr$pr_t <- log1p(tr$pr)
  tr
}

train <- if (file.exists(file_train)) readRDS(file_train)
if (is.null(train) || !all(unique(inv$station) %in% train$station)) {
  train <- build_training(inv); saveRDS(train, file_train)
}
train <- train[train$station %in% unique(inv$station), ]

# Join Rsw-cc
rsw <- bind_rows(lapply(unique(train$station), function(s) {
  x <- read.csv(file.path(dir_in_rsw, paste0("meteoland_", tolower(s), "_rswcc.csv")))
  data.frame(station = s, date = as.Date(x$date), rsw_cc = x$rsw_cc)
}))
train <- left_join(select(train, -any_of("rsw_cc")), rsw, by = c("station", "date"))
stopifnot(!anyNA(train$rsw_cc))

# Random forest
fit_rf <- function(d) {
  set.seed(seed)
  ranger::ranger(fml, data = d, num.trees = n_trees, num.threads = n_threads)
}
predict_rf <- function(fit, newdata) {
  pmax(predict(fit, newdata, num.threads = n_threads)$predictions, 0)
}
summ <- function(d) summarise(d, MAE = mean(abs(pred - rsds)), RMSE = sqrt(mean((pred - rsds)^2)),
  Bias = mean(pred - rsds), R2 = 1 - sum((pred - rsds)^2) / sum((rsds - mean(rsds))^2))

# Leave-one-station-out CV (only if run_cv_rg = TRUE)
if (run_cv_rg) {
  cv <- bind_rows(lapply(sort(unique(train$station)), function(st) {
    message("CV: ", st)
    trn <- train[train$station != st, ]; tst <- train[train$station == st, ]
    data.frame(tst[c("station", "elev", "rsds")], pred = predict_rf(fit_rf(trn), tst))
  }))
  cv$band <- cut(cv$elev, c(0, 500, 1000, 1500, 2000, Inf),
                 labels = c("<500", "500-1000", "1000-1500", "1500-2000", ">2000"))

  overall    <- cv |> summ()
  by_station <- cv |> group_by(station, elev) |> summ()
  by_elev    <- cv |> group_by(band) |> summ()
  print(overall); print(by_elev)
  write.csv(overall,    here("tables", "rg_cv_overall.csv"),      row.names = FALSE)
  write.csv(by_station, here("tables", "rg_cv_station.csv"),   row.names = FALSE)
  write.csv(by_elev,    here("tables", "rg_cv_elevation.csv"), row.names = FALSE)

  # Figures
  ggsave(here("figures", "rg_cv_obs_vs_pred.jpg"),
         ggplot(cv, aes(rsds, pred)) + geom_bin2d(bins = 60) + geom_abline(colour = "red") +
           scale_fill_viridis_c(trans = "log10") + coord_equal() +
           labs(x = "Observed rsds (W m-2)", y = "Predicted rsds (W m-2)", fill = "n") + theme_bw(),
         width = 6, height = 6, dpi = 200)
  ggsave(here("figures", "rg_cv_rmse_by_elevation.jpg"),
         ggplot(by_station, aes(elev, RMSE)) + geom_point() +
           labs(x = "Station elevation (m)", y = "RMSE (W m-2)") + theme_bw(),
         width = 6, height = 4, dpi = 200)
}

# Final random forest on all stations
saveRDS(list(model = fit_rf(train), predictors = x_vars),
        file.path(dir_model, "rg_model_rf.rds"))


# Step 2: Train a random forest for net radiation (Rn) on four SwissMetNet stations
# Rn ~ Rg + Tmin + Tmax + Prec + Rsw-cc, with Rg simulated by the step 1 model

# Step 2.1: Settings
dir_in_smn   <- here("data-raw", "swissmetnet")
smn_stations <- c("ALT", "MAG", "PAY", "TAE")   # only stations measuring all four components
x_vars_rn    <- c("rg_sim", "tasmin", "tasmax", "pr_t", "rsw_cc")
fml_rn       <- reformulate(x_vars_rn, response = "rn")

# Step 2.2: Read SwissMetNet daily data and derive observed Rn (all in W m-2, daily means)
# Rn = global radiation (gre000d0) - reflected shortwave (osr000d0)
#      + incoming longwave (oli000d0) - outgoing longwave (olo000d0)
smn_meta <- read.delim(file.path(dir_in_smn, "ogd-smn_meta_stations.csv"),
                       sep = ";", fileEncoding = "latin1")
read_smn <- function(s) {
  x <- read.delim(file.path(dir_in_smn, paste0("ogd-smn_", tolower(s), "_d_historical.csv")),
                  sep = ";", na.strings = c("", "NA"))
  data.frame(station = s, date = as.Date(x$reference_timestamp, "%d.%m.%Y"),
             elev   = smn_meta$station_height_masl[smn_meta$station_abbr == s],
             pr     = x$rre150d0, tasmin = x$tre200dn, tasmax = x$tre200dx,
             rg_obs = x$gre000d0,
             rn     = x$gre000d0 - x$osr000d0 + x$oli000d0 - x$olo000d0)
}
smn <- bind_rows(lapply(smn_stations, read_smn))
smn <- smn[complete.cases(smn), ]
smn$pr_t <- log1p(smn$pr)
print(smn |> group_by(station, elev) |> summarise(n = n(), from = min(date), to = max(date)))

# Step 2.3: Join Rsw-cc by day of year (rsw-cc files cover doy 1-365, so Feb 29 takes Feb 28)
rsw_doy <- bind_rows(lapply(smn_stations, function(s) {
  x <- read.csv(file.path(dir_in_rsw, paste0("meteoland_", tolower(s), "_rswcc.csv")))
  x |> group_by(doy) |> summarise(rsw_cc = mean(rsw_cc)) |> mutate(station = s)
}))
yr   <- as.integer(format(smn$date, "%Y"))
leap <- yr %% 4 == 0 & (yr %% 100 != 0 | yr %% 400 == 0)
smn$doy <- as.integer(format(smn$date, "%j"))
smn$doy <- ifelse(leap & smn$doy > 59, smn$doy - 1L, smn$doy)
smn <- left_join(smn, rsw_doy, by = c("station", "doy"))
stopifnot(!anyNA(smn$rsw_cc))

# Step 2.4: Simulate Rg with the step 1 model (and check it against observed Rg)
rg_model   <- readRDS(file.path(dir_model, "rg_model_rf.rds"))
smn$rg_sim <- predict_rf(rg_model$model, smn)
rg_check   <- smn |> transmute(station, rsds = rg_obs, pred = rg_sim) |> group_by(station) |> summ()
print(rg_check)
write.csv(rg_check, here("tables", "rg_check_smn.csv"), row.names = FALSE)

# Step 2.5: Random forest for Rn (no clipping at 0, Rn can be negative)
fit_rf_rn <- function(d) {
  set.seed(seed)
  ranger::ranger(fml_rn, data = d, num.trees = n_trees, num.threads = n_threads)
}
predict_rf_rn <- function(fit, newdata) {
  predict(fit, newdata, num.threads = n_threads)$predictions
}

# Step 2.6: Leave-one-station-out CV over the four stations
cv_rn <- bind_rows(lapply(smn_stations, function(st) {
  message("CV Rn: ", st)
  trn <- smn[smn$station != st, ]; tst <- smn[smn$station == st, ]
  data.frame(tst[c("station", "elev", "date", "rn")], pred = predict_rf_rn(fit_rf_rn(trn), tst))
}))
summ_rn <- function(d) summarise(d, MAE = mean(abs(pred - rn)), RMSE = sqrt(mean((pred - rn)^2)),
  Bias = mean(pred - rn), R2 = 1 - sum((pred - rn)^2) / sum((rn - mean(rn))^2))
overall_rn    <- cv_rn |> summ_rn()
by_station_rn <- cv_rn |> group_by(station, elev) |> summ_rn()
print(overall_rn); print(by_station_rn)
write.csv(overall_rn,    here("tables", "rn_cv_overall.csv"),    row.names = FALSE)
write.csv(by_station_rn, here("tables", "rn_cv_station.csv"), row.names = FALSE)

# Step 2.7: Figures
ggsave(here("figures", "rn_cv_obs_vs_pred.jpg"),
       ggplot(cv_rn, aes(rn, pred)) + geom_bin2d(bins = 60) + geom_abline(colour = "red") +
         scale_fill_viridis_c(trans = "log10") + facet_wrap(~station) + coord_equal() +
         labs(x = "Observed Rn (W m-2)", y = "Predicted Rn (W m-2)", fill = "n") + theme_bw(),
       width = 8, height = 8, dpi = 200)

# Step 2.8: Final Rn random forest on all four stations
saveRDS(list(model = fit_rf_rn(smn), predictors = x_vars_rn),
        file.path(dir_model, "rn_model_rf.rds"))


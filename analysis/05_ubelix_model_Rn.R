#!/usr/bin/env Rscript
#SBATCH --job-name=model_Rn
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --output=model_Rn_%j.log

# run analysis/03_ubelix_calc_rsw-cc.R first to create the rsw-cc files in data/rsw-cc

# two step approach: first train a random forest on the CH2025 daily local data then use the trained model to predict net radiation on four swissmetnet stations

# Step 1: Train a random forest on the CH2025 daily local data

library(here)
library(dplyr)
library(ggplot2)

# Settings
dir_in_local <- here("data-raw", "CH2025", "DAILY_LOCAL")
dir_in_rsw   <- here("data", "rsw-cc")
vars         <- c("pr", "tasmin", "tasmax", "rsds")
n_per_group  <- 40
n_trees <- 100; seed <- 1
n_threads    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
x_vars  <- c("pr_t", "tasmin", "tasmax", "rsw_cc", "elev")
fml     <- reformulate(x_vars, response = "rsds")
dir_model  <- here("data", "radiation_models")
file_train <- file.path(dir_model, "rg_training.rds")
for (d in c(dir_model, here("tables"), here("figures"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Reader for one CH2025 csv (table starts at the line "DATE;")
read_ch2025 <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  out   <- read.delim(text = lines[which(startsWith(lines, "DATE;"))[1]:length(lines)],
                      sep = ";", check.names = FALSE, na.strings = c("", "NA", "NaN"))
  out$DATE <- as.Date(out$DATE)
  elev <- lines[startsWith(lines, "STATION_HEIGHT;")][1]
  attr(out, "elev") <- as.numeric(sub("^STATION_HEIGHT;", "", elev))
  out
}

# Reader for the Rsw-cc csv of one station
read_rsw <- function(s) {
  x <- read.csv(file.path(dir_in_rsw, paste0("meteoland_", tolower(s), "_rswcc.csv")))
  data.frame(station = s, date = as.Date(x$date), doy = x$doy, rsw_cc = x$rsw_cc)
}

# File inventory
files <- list.files(dir_in_local, recursive = TRUE, full.names = TRUE,
                    pattern = "^ogd-climate-scenarios-ch2025_[A-Za-z0-9]+_[a-z]+_.+\\.csv$")
parts <- do.call(rbind, strsplit(sub("\\.csv$", "", basename(files)), "_"))
inv   <- data.frame(path = files, station = toupper(parts[, 2]),
                    var = parts[, 3], period = parts[, 4])
inv   <- inv[inv$var %in% vars, ]

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
                 period = sets$period[i], chain = ch, date = x$pr$DATE[idx], d[idx, ])
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
rsw   <- bind_rows(lapply(unique(train$station), read_rsw))
train <- left_join(train, rsw, by = c("station", "date"))
stopifnot(!anyNA(train$rsw_cc))

# Random forest; min = 0 clips Rg predictions (Rn can be negative)
fit_rf <- function(d, formula) {
  set.seed(seed)
  ranger::ranger(formula, data = d, num.trees = n_trees, num.threads = n_threads)
}
predict_rf <- function(fit, newdata, min = -Inf) {
  pmax(predict(fit, newdata, num.threads = n_threads)$predictions, min)
}
summ <- function(d, obs) {
  summarise(d, MAE = mean(abs(pred - {{ obs }})), RMSE = sqrt(mean((pred - {{ obs }})^2)),
            Bias = mean(pred - {{ obs }}),
            R2 = 1 - sum((pred - {{ obs }})^2) / sum(({{ obs }} - mean({{ obs }}))^2))
}

# Final random forest on all stations
rg_fit <- fit_rf(train, fml)
saveRDS(list(model = rg_fit, predictors = x_vars), file.path(dir_model, "rg_model_rf.rds"))

# Leave-one-station-out CV
cv <- bind_rows(lapply(sort(unique(train$station)), function(st) {
  message("CV: ", st)
  trn <- train[train$station != st, ]; tst <- train[train$station == st, ]
  data.frame(tst[c("station", "elev", "rsds")], pred = predict_rf(fit_rf(trn, fml), tst, min = 0))
}))
cv$band <- cut(cv$elev, c(0, 500, 1000, 1500, 2000, Inf),
               labels = c("<500", "500-1000", "1000-1500", "1500-2000", ">2000"))

by_station <- cv |> group_by(station, elev) |> summ(rsds)
write.csv(cv |> summ(rsds),                  here("tables", "rg_cv_overall.csv"),   row.names = FALSE)
write.csv(by_station,                        here("tables", "rg_cv_station.csv"),   row.names = FALSE)
write.csv(cv |> group_by(band) |> summ(rsds), here("tables", "rg_cv_elevation.csv"), row.names = FALSE)

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


# Step 2: Train a random forest for net radiation (Rn) on four SwissMetNet stations that contain necessary variables incoming longwave, outgoing longwave, reflected shortwave, and incoming shortwave.
# Rn ~ Rg + Tmin + Tmax + Prec + Rsw-cc, with Rg simulated by the step 1 model

# Step 2.1: Settings
dir_in_smn   <- here("data-raw", "swissmetnet")
smn_stations <- c("ALT", "MAG", "PAY", "TAE")   # only stations measuring all four components
x_vars_rn    <- c("rg_sim", "tasmin", "tasmax", "pr_t", "rsw_cc")
fml_rn       <- reformulate(x_vars_rn, response = "rn")

# Step 2.2: Download the historical daily csv of the four stations
smn_url  <- "https://data.geo.admin.ch/ch.meteoschweiz.ogd-smn"
smn_code <- tolower(smn_stations)
dir.create(dir_in_smn, showWarnings = FALSE, recursive = TRUE)

for (url in file.path(smn_url, smn_code, paste0("ogd-smn_", smn_code, "_d_historical.csv"))) {
  dest <- file.path(dir_in_smn, basename(url))
  if (!file.exists(dest)) {
    tmp <- paste0(dest, ".part")
    download.file(url, tmp, mode = "wb", method = "libcurl", quiet = TRUE)
    file.rename(tmp, dest)
  }
}

# Step 2.3: Read SwissMetNet daily data and derive observed Rn (all in W m-2, daily means)
# Rn = global radiation (gre000d0) - reflected shortwave (osr000d0)
#      + incoming longwave (oli000d0) - outgoing longwave (olo000d0)
# elevation comes from the training data, the four stations are all Rg model stations
stopifnot(smn_stations %in% train$station)
read_smn <- function(s) {
  x <- read.delim(file.path(dir_in_smn, paste0("ogd-smn_", tolower(s), "_d_historical.csv")),
                  sep = ";", na.strings = c("", "NA"))
  data.frame(station = s, date = as.Date(x$reference_timestamp, "%d.%m.%Y"),
             elev   = train$elev[match(s, train$station)],
             pr     = x$rre150d0, tasmin = x$tre200dn, tasmax = x$tre200dx,
             rg_obs = x$gre000d0,
             rn     = x$gre000d0 - x$osr000d0 + x$oli000d0 - x$olo000d0)
}
smn <- bind_rows(lapply(smn_stations, read_smn))
smn <- smn[complete.cases(smn), ]
smn$pr_t <- log1p(smn$pr)

# Step 2.4: Join Rsw-cc by day of year (rsw-cc files cover doy 1-365, so Feb 29 takes Feb 28)
rsw_doy <- bind_rows(lapply(smn_stations, read_rsw)) |>
  group_by(station, doy) |> summarise(rsw_cc = mean(rsw_cc), .groups = "drop")
doy  <- as.integer(format(smn$date, "%j"))
leap <- as.integer(format(as.Date(format(smn$date, "%Y-12-31")), "%j")) == 366
smn$doy <- doy - (leap & doy > 59)
smn <- left_join(smn, rsw_doy, by = c("station", "doy"))
stopifnot(!anyNA(smn$rsw_cc))

# Step 2.5: Simulate Rg with the step 1 model (and check it against observed Rg)
smn$rg_sim <- predict_rf(rg_fit, smn, min = 0)
write.csv(smn |> mutate(pred = rg_sim) |> group_by(station) |> summ(rg_obs),
          here("tables", "rg_check_smn.csv"), row.names = FALSE)

# Step 2.6: Leave-one-station-out CV over the four stations (no clipping, Rn can be negative)
cv_rn <- bind_rows(lapply(smn_stations, function(st) {
  message("CV Rn: ", st)
  trn <- smn[smn$station != st, ]; tst <- smn[smn$station == st, ]
  data.frame(tst[c("station", "elev", "date", "rn")], pred = predict_rf(fit_rf(trn, fml_rn), tst))
}))
write.csv(cv_rn |> summ(rn), here("tables", "rn_cv_overall.csv"), row.names = FALSE)
write.csv(cv_rn |> group_by(station, elev) |> summ(rn), here("tables", "rn_cv_station.csv"),
          row.names = FALSE)

# Step 2.7: Final Rn random forest on all four stations
saveRDS(list(model = fit_rf(smn, fml_rn), predictors = x_vars_rn),
        file.path(dir_model, "rn_model_rf.rds"))

# Step 2.8: Figures
ggsave(here("figures", "rn_cv_obs_vs_pred.jpg"),
       ggplot(cv_rn, aes(rn, pred)) + geom_bin2d(bins = 60) + geom_abline(colour = "red") +
         scale_fill_viridis_c(trans = "log10") + facet_wrap(~station) + coord_equal() +
         labs(x = "Observed Rn (W m-2)", y = "Predicted Rn (W m-2)", fill = "n") + theme_bw(),
       width = 8, height = 8, dpi = 200)

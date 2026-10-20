# run calc_rsw-cc.R first to create the rsw-cc files in data/rsw-cc

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
n_threads    <- max(1, parallel::detectCores() - 1)
x_vars  <- c("pr_t", "tasmin", "tasmax", "rsw_cc", "elev")
fml     <- reformulate(x_vars, response = "rsds")
file_train <- here("data", "rg_training.rds")
for (d in c("data", "tables", "figures")) dir.create(here(d), showWarnings = FALSE)

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

# Leave-one-station-out CV
cv <- bind_rows(lapply(sort(unique(train$station)), function(st) {
  message("CV: ", st)
  trn <- train[train$station != st, ]; tst <- train[train$station == st, ]
  data.frame(tst[c("station", "elev", "rsds")], pred = predict_rf(fit_rf(trn), tst))
}))
cv$band <- cut(cv$elev, c(0, 500, 1000, 1500, 2000, Inf),
               labels = c("<500", "500-1000", "1000-1500", "1500-2000", ">2000"))

summ <- function(d) summarise(d, MAE = mean(abs(pred - rsds)), RMSE = sqrt(mean((pred - rsds)^2)),
  Bias = mean(pred - rsds), R2 = 1 - sum((pred - rsds)^2) / sum((rsds - mean(rsds))^2))
overall    <- cv |> summ()
by_station <- cv |> group_by(station, elev) |> summ()
by_elev    <- cv |> group_by(band) |> summ()
print(overall); print(by_elev)
write.csv(overall,    here("tables", "rg_cv_metrics.csv"),      row.names = FALSE)
write.csv(by_station, here("tables", "rg_cv_by_station.csv"),   row.names = FALSE)
write.csv(by_elev,    here("tables", "rg_cv_by_elevation.csv"), row.names = FALSE)

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

# Final random forest on all stations
saveRDS(list(model = fit_rf(train), predictors = x_vars),
        here("data", "rg_model_rf.rds"))
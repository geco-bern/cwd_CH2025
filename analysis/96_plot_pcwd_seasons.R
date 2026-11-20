#!/usr/bin/env Rscript
#SBATCH --job-name=pcwd_seasons
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker
#SBATCH --partition=cpu-invest
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --output=ubelix_logs/pcwd_seasons_%j.log

# run 07_apply_pcwd_UBELIX.R first to create the daily PCWD grids of all model chains and periods
# seasonal mean PCWD over Switzerland per model chain and period:
#   1. per file (model chain x period): mean of the daily PCWD over all cells of the domain and all days
#      of the season in the 30 years (winter DJF, spring MAM, summer JJA, autumn SON)
#   2. per period and season: ensemble mean, interquartile range and range over the model chains
# one panel per season (2 x 2, same y axis), x axis = period, y axis = PCWD,
# bars = ensemble mean (reference grey, GWLs yellow to dark red)
# per-file results are kept in data/pcwd_seasons/, so a resubmitted job skips finished files;
# the combined table is saved to tables/, the figure to figures/

library(here)
library(terra)
library(ncdf4)
library(parallel)
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout (n_days)

# ---- 1. Settings ----------------------------------------------------------
block_rows <- 10                                  # rows of the grid read at once
n_cores    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

periods <- c("ref91-20", "gwl1.5", "gwl2.0", "gwl2.5", "gwl3.0")
labels  <- c("ref91-20" = "Ref.\n1991-2020",
             "gwl1.5"   = "GWL\n1.5 °C",
             "gwl2.0"   = "GWL\n2.0 °C",
             "gwl2.5"   = "GWL\n2.5 °C",
             "gwl3.0"   = "GWL\n3.0 °C")
seasons <- c(DJF = "Winter (DJF)", MAM = "Spring (MAM)", JJA = "Summer (JJA)", SON = "Autumn (SON)")
# reference grey, GWLs from yellow to dark red
cols    <- c("ref91-20" = "grey60",
             "gwl1.5"   = "#f6c342",
             "gwl2.0"   = "#ef8a2c",
             "gwl2.5"   = "#d4442a",
             "gwl3.0"   = "#8b1a1a")

dir_cache <- here("data", "pcwd_seasons")         # one small csv per model chain and period
out_csv   <- here("tables", "pcwd_seasons_ch.csv")
dir.create(dir_cache, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Inventory of the PCWD files: grid_ch_pcwd_<model chain>_<period>.nc ---------
files <- list.files(here("data", "pcwd"), pattern = "^grid_ch_pcwd_.*\\.nc$", full.names = TRUE)
parts <- do.call(rbind, strsplit(sub("\\.nc$", "", basename(files)), "_"))
inv   <- data.frame(path = files, chain = parts[, 4], period = parts[, 5])
inv   <- inv[inv$period %in% periods, ]
inv$cache <- file.path(dir_cache, sub("\\.nc$", ".csv", basename(inv$path)))
print(table(inv$period))                          # 26 model chains per period (24 for GWL 3.0)

# season of every day (365_day calendar, 30 years)
month_end <- cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
month     <- findInterval(rep(0:364, times = n_days / 365), month_end) + 1
season    <- c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")[month]

# ---- 3. Seasonal mean over the domain, one file at a time ---------------------------
# sum and number of valid cell-days per season, accumulated over the blocks of rows;
# cells outside the domain are NA and not counted, so every cell of the domain has the same weight
season_means <- function(i) {
  r <- rast(inv$path[i], subds = "pcwd")
  stopifnot(nlyr(r) == n_days)
  tot <- setNames(numeric(4), names(seasons))
  cnt <- setNames(numeric(4), names(seasons))

  readStart(r)
  on.exit(readStop(r))
  for (s in seq(1, nrow(r), by = block_rows)) {
    x <- readValues(r, row = s, nrows = min(block_rows, nrow(r) - s + 1), mat = TRUE)   # cells x days
    day_sum <- colSums(x, na.rm = TRUE)
    day_n   <- colSums(!is.na(x))
    tot <- tot + tapply(day_sum, season, sum)[names(seasons)]
    cnt <- cnt + tapply(day_n,   season, sum)[names(seasons)]
  }

  res <- data.frame(chain = inv$chain[i], period = inv$period[i], season = names(seasons),
                    pcwd = as.numeric(tot / cnt), n_cell_days = as.numeric(cnt))
  write.csv(res, inv$cache[i], row.names = FALSE)
  message("Done: ", basename(inv$path[i]), " (", Sys.time(), ")")
  TRUE
}

todo <- which(!file.exists(inv$cache))
message(length(todo), " of ", nrow(inv), " files to process")
res  <- mclapply(todo, function(i) try(season_means(i)), mc.cores = n_cores, mc.preschedule = FALSE)
for (e in Filter(function(r) inherits(r, "try-error"), res)) message(e)

failed <- inv$path[!file.exists(inv$cache)]
if (length(failed) > 0) stop(length(failed), " failed, resubmit:\n", paste(failed, collapse = "\n"))

chains <- do.call(rbind, lapply(inv$cache, read.csv))
dir.create(here("tables"), showWarnings = FALSE)
write.csv(chains, out_csv, row.names = FALSE)

# ---- 4. Ensemble statistics over the model chains, per period and season ------------
ens <- do.call(rbind, lapply(split(chains, list(chains$period, chains$season), drop = TRUE), function(x)
  data.frame(period = x$period[1], season = x$season[1], n_chains = nrow(x),
             mean = mean(x$pcwd), min = min(x$pcwd), max = max(x$pcwd),
             q25  = quantile(x$pcwd, 0.25, names = FALSE),
             q75  = quantile(x$pcwd, 0.75, names = FALSE))))
print(ens[order(ens$season, match(ens$period, periods)), ], row.names = FALSE)

# factors for the panel and x axis order, colour per period
for (v in c("chains", "ens")) {
  x <- get(v)
  x$season_lab <- factor(seasons[x$season], levels = seasons)        # Winter, Spring / Summer, Autumn
  x$period_lab <- factor(labels[x$period], levels = labels[periods])
  x$col        <- cols[x$period]
  assign(v, x)
}

# ---- 5. Plot: one panel per season, 2 x 2 ------------------------------------------
p <- ggplot(ens, aes(period_lab)) +
  geom_col(aes(y = mean, fill = col), width = 0.75) +                             # ensemble mean
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.2, linewidth = 0.4,
                colour = "grey15") +                                               # range
  geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 2, colour = "grey15") +  # interquartile range
  geom_point(data = chains, aes(y = pcwd), position = position_jitter(width = 0.15, height = 0, seed = 1),
             size = 0.8, alpha = 0.5, colour = "grey15") +                         # single model chains
  facet_wrap(~season_lab, ncol = 2) +                                              # same y axis in all panels
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title    = "Seasonal mean PCWD over Switzerland",
       subtitle = "Mean of the daily PCWD over all grid cells and all days of the season in the 30 years.\nBar: ensemble mean; thick line: interquartile range; whiskers: range; dots: single model chains.",
       x = NULL, y = "PCWD (mm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border       = element_rect(fill = NA, colour = "grey80"),
        strip.text         = element_text(hjust = 0, face = "bold"),
        plot.subtitle      = element_text(colour = "#52514e", size = 9))

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "pcwd_seasons_ch.jpg"), p, width = 9, height = 7.5, dpi = 300, bg = "white")
message("Saved figure")

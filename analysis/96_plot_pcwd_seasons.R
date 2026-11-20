# run 07_apply_pcwd_UBELIX.R first to create the daily PCWD grids of all model chains and periods
# seasonal mean PCWD over Switzerland per model chain and period:
#   1. per file (model chain x period): mean of the daily PCWD over all cells of the domain and all days
#      of the season in the 30 years (winter DJF, spring MAM, summer JJA, autumn SON)
#   2. per period and season: ensemble mean and spread over the model chains
# one panel per season (2 x 2, same y axis), in the style of the CH2025 figures:
# box = spread of the model chains, line = ensemble mean; reference grey, GWLs yellow to dark red
# per-file results are kept in data/pcwd_seasons/ and reused; the combined table is saved to tables/,
# the figure to figures/

library(here)
library(ncdf4)
library(terra)                                    # needed by R/nc_grid_ch.R
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout (n_days)

# ---- 1. Settings ----------------------------------------------------------
block_rows <- 10                                  # rows of the grid read at once
spread     <- c(0, 1)                             # box: quantiles of the model chains (0, 1 = full range;
                                                  #      0.25, 0.75 = interquartile range)

periods <- c("ref91-20", "gwl1.5", "gwl2.0", "gwl2.5", "gwl3.0")
labels  <- c("ref91-20" = "1991-2020",
             "gwl1.5"   = "GWL1.5",
             "gwl2.0"   = "GWL2.0",
             "gwl2.5"   = "GWL2.5",
             "gwl3.0"   = "GWL3.0")
cols    <- c("ref91-20" = "grey55",               # reference grey, GWLs yellow to dark red
             "gwl1.5"   = "#d4a537",
             "gwl2.0"   = "#e07b39",
             "gwl2.5"   = "#c8432b",
             "gwl3.0"   = "#7f1d1d")
seasons <- c(DJF = "winter (DJF)", MAM = "spring (MAM)", JJA = "summer (JJA)", SON = "autumn (SON)")

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
  nc <- nc_open(inv$path[i])
  on.exit(nc_close(nc))
  nN <- nc$dim$N$len
  stopifnot(nc$dim$time$len == n_days)
  tot <- setNames(numeric(4), names(seasons))
  cnt <- setNames(numeric(4), names(seasons))

  for (s in seq(1, nN, by = block_rows)) {                  # blocks match the chunks of 07 (10 rows)
    n <- min(block_rows, nN - s + 1)
    x <- ncvar_get(nc, "pcwd", start = c(1, s, 1), count = c(-1, n, -1), collapse_degen = FALSE)
    x <- matrix(x, ncol = n_days)                           # cells x days
    tot <- tot + tapply(colSums(x, na.rm = TRUE), season, sum)[names(seasons)]
    cnt <- cnt + tapply(colSums(!is.na(x)),       season, sum)[names(seasons)]
  }

  res <- data.frame(chain = inv$chain[i], period = inv$period[i], season = names(seasons),
                    pcwd = as.numeric(tot / cnt), n_cell_days = as.numeric(cnt))
  write.csv(res, inv$cache[i], row.names = FALSE)
  message("Done: ", basename(inv$path[i]))
}

todo <- which(!file.exists(inv$cache))
message(length(todo), " of ", nrow(inv), " files to process")
for (i in todo) season_means(i)

chains <- do.call(rbind, lapply(inv$cache, read.csv))
dir.create(here("tables"), showWarnings = FALSE)
write.csv(chains, out_csv, row.names = FALSE)

# ---- 4. Ensemble statistics over the model chains, per period and season ------------
ens <- do.call(rbind, lapply(split(chains, list(chains$period, chains$season), drop = TRUE), function(x)
  data.frame(period = x$period[1], season = x$season[1], n_chains = nrow(x),
             mean = mean(x$pcwd),
             lo   = quantile(x$pcwd, spread[1], names = FALSE),
             hi   = quantile(x$pcwd, spread[2], names = FALSE))))
print(ens[order(ens$season, match(ens$period, periods)), ], row.names = FALSE)

ens$x      <- match(ens$period, periods)                     # x position of each period
ens$col    <- cols[ens$period]
ens$season <- factor(seasons[ens$season], levels = seasons)  # winter, spring / summer, autumn

# ---- 5. Plot: one panel per season, 2 x 2 ------------------------------------------
w <- 0.32                                                    # half width of the boxes

p <- ggplot(ens) +
  geom_rect(aes(xmin = x - w, xmax = x + w, ymin = lo, ymax = hi, fill = col), alpha = 0.45) +  # spread
  geom_segment(aes(x = x - w, xend = x + w, y = mean, yend = mean, colour = col),
               linewidth = 1.3) +                                                              # mean
  # "free" only to draw the axes in every panel; the limits below keep the same scale in all panels
  facet_wrap(~season, ncol = 2, scales = "free") +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(breaks = seq_along(periods), labels = labels[periods],
                     limits = c(0.5, length(periods) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(ens$hi)), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "PCWD (mm)") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey65", linetype = "dashed", linewidth = 0.4),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title.y       = element_text(face = "bold"),
        axis.text          = element_text(colour = "black"),
        panel.spacing      = unit(1.5, "lines"))

print(p)

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "pcwd_seasons_ch.jpg"), p, width = 10, height = 7, dpi = 300, bg = "white")

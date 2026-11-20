# run 07_apply_pcwd_UBELIX.R first to create the daily PCWD grids of all model chains and periods
# maps of the mean PCWD per grid cell, ensemble mean over all model chains:
#   1. per file (model chain x period): mean of the daily PCWD per cell over all days of the season
#      in the 30 years (winter DJF, spring MAM, summer JJA, autumn SON) and over all days (annual)
#   2. per period: mean of these maps over the model chains
# one figure: rows = periods (reference on top, then GWL 1.5 to 3.0), columns = DJF, MAM, JJA, SON, annual,
# one shared binned colour scale; the figure is saved to figures/

library(here)
library(ncdf4)
library(terra)                                    # needed by R/nc_grid_ch.R
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout (E_vals, N_vals, n_days)

# ---- 1. Settings ----------------------------------------------------------
block_rows <- 10                                  # rows of the grid read at once (= chunks of 07)

periods <- c("ref91-20", "gwl1.5", "gwl2.0", "gwl2.5", "gwl3.0")
labels  <- c("ref91-20" = "Reference\n1991-2020",
             "gwl1.5"   = "GWL 1.5 °C",
             "gwl2.0"   = "GWL 2.0 °C",
             "gwl2.5"   = "GWL 2.5 °C",
             "gwl3.0"   = "GWL 3.0 °C")
cols_season <- c("DJF", "MAM", "JJA", "SON", "annual")

# class limits of the colour scale in mm (finer at low values, where most of the map lies);
# limits above the largest value are dropped, the largest value is added if it is above all limits
breaks_all <- c(0, 1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 400, 600)

# ---- 2. Inventory of the PCWD files: grid_ch_pcwd_<model chain>_<period>.nc ---------
files <- list.files(here("data", "pcwd"), pattern = "^grid_ch_pcwd_.*\\.nc$", full.names = TRUE)
if (length(files) == 0) stop("No PCWD files in ", here("data", "pcwd"))
parts <- do.call(rbind, strsplit(sub("\\.nc$", "", basename(files)), "_"))
inv   <- data.frame(path = files, chain = parts[, 4], period = parts[, 5])
inv   <- inv[inv$period %in% periods, ]
print(table(inv$period))                          # 26 model chains per period (24 for GWL 3.0)

# season of every day (365_day calendar, 30 years)
month_end <- cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
month     <- findInterval(rep(0:364, times = n_days / 365), month_end) + 1
season    <- c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")[month]

# ---- 3. Mean PCWD per cell and season, one file at a time -------------------------
# returns a matrix: cells (E fastest, N south to north, as in the file) x DJF, MAM, JJA, SON, annual;
# cells outside the domain stay NA
n_E <- length(E_vals)
n_N <- length(N_vals)

cell_means <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  stopifnot(nc$dim$time$len == n_days, nc$dim$E$len == n_E, nc$dim$N$len == n_N)
  out <- matrix(NA_real_, n_E * n_N, length(cols_season), dimnames = list(NULL, cols_season))

  for (s in seq(1, n_N, by = block_rows)) {
    n   <- min(block_rows, n_N - s + 1)
    x   <- ncvar_get(nc, "pcwd", start = c(1, s, 1), count = c(-1, n, -1), collapse_degen = FALSE)
    x   <- matrix(x, ncol = n_days)                     # cells of the block x days
    idx <- (s - 1) * n_E + seq_len(n_E * n)             # position of these cells in the grid
    for (k in cols_season[1:4]) out[idx, k] <- rowMeans(x[, season == k, drop = FALSE])
    out[idx, "annual"] <- rowMeans(x)
  }
  out
}

# ensemble mean per period: sum over the model chains / number of model chains
ens <- list()
for (per in intersect(periods, inv$period)) {
  paths <- inv$path[inv$period == per]
  tot   <- 0
  for (f in paths) {
    message(per, ": ", basename(f))
    tot <- tot + cell_means(f)
  }
  ens[[per]] <- tot / length(paths)
}

# ---- 4. Long table for the plot -------------------------------------------------
xy <- expand.grid(x = E_vals, y = N_vals)             # same cell order as in the files
d  <- do.call(rbind, lapply(names(ens), function(per) do.call(rbind, lapply(cols_season, function(k)
  data.frame(xy, period = per, season = k, pcwd = ens[[per]][, k])))))
d  <- d[!is.na(d$pcwd), ]

d$period <- factor(labels[d$period], levels = labels[names(ens)])
d$season <- factor(d$season, levels = cols_season)

# classes of the colour scale
v_max <- max(d$pcwd)
upper <- if (any(breaks_all >= v_max)) min(breaks_all[breaks_all >= v_max]) else ceiling(v_max)
brks  <- c(breaks_all[breaks_all < v_max], upper)     # upper limit at or above the largest value
d$bin <- cut(d$pcwd, brks, include.lowest = TRUE, dig.lab = 4)
bin_cols <- setNames(hcl.colors(nlevels(d$bin), "YlOrRd", rev = TRUE), levels(d$bin))
bin_labs <- paste(head(brks, -1), "-", brks[-1])      # legend labels, e.g. "10 - 20"

# panel letters (a), (b), ... row by row, as in the example
pan <- expand.grid(season = levels(d$season), period = levels(d$period))
pan$letter <- sprintf("(%s)", letters[seq_len(nrow(pan))])
pan$season <- factor(pan$season, levels = levels(d$season))
pan$period <- factor(pan$period, levels = levels(d$period))

# ---- 5. Plot ---------------------------------------------------------------
x_lim <- range(d$x) + c(-500, 500)                    # extent of the cells with data (cell edges)
y_lim <- range(d$y) + c(-500, 500)

p <- ggplot(d, aes(x, y, fill = bin)) +
  geom_raster(show.legend = TRUE) +                    # TRUE: also show classes without cells
  geom_text(data = pan, aes(x = x_lim[1], y = y_lim[2], label = letter),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5) +
  facet_grid(period ~ season, switch = "y") +          # periods as row labels on the left
  coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_manual(values = bin_cols, labels = bin_labs, drop = FALSE,
                    name = "PCWD (mm)") +
  guides(fill = guide_legend(reverse = TRUE, keyheight = unit(0.7, "cm"))) +
  theme_void(base_size = 12) +
  theme(strip.text.x      = element_text(size = 13, margin = margin(b = 4)),
        strip.text.y.left = element_text(size = 12, angle = 90, margin = margin(r = 4)),
        strip.placement   = "outside",
        panel.spacing     = unit(0.4, "lines"),
        legend.title      = element_text(size = 11),
        legend.text       = element_text(size = 10),
        plot.background   = element_rect(fill = "white", colour = NA))

print(p)

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "pcwd_mean_map_ch.jpg"), p, width = 16, height = 10, dpi = 300, bg = "white")

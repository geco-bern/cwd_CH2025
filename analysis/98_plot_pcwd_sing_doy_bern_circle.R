# run 07_apply_pcwd_UBELIX.R first to create the daily PCWD grids of all model chains and periods
# daily PCWD for the grid cell of Bern-Zollikofen as circle plots: one full turn = one year
# (1 January at the top, clockwise), distance from the centre = PCWD;
# one circle per period (reference period and the four GWLs) with all 30 years of all model chains
# as thin lines and the mean over all years and model chains per day of year as a thick line
# one figure per period and one figure with all periods stacked, all with the same radial scale;
# reference grey, GWLs yellow to dark red (as in 94_ and 96_)

library(here)
library(ncdf4)
library(terra)                                    # needed by R/nc_grid_ch.R
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout (E_vals, N_vals, n_days)

# ---- 1. Settings ----------------------------------------------------------
site    <- "Bern-Zollikofen"
site_E  <- 2601934                                # CH1903+ / LV95
site_N  <- 1204410

periods <- c("ref91-20", "gwl1.5", "gwl2.0", "gwl2.5", "gwl3.0")
labels  <- c("ref91-20" = "Reference 1991-2020",
             "gwl1.5"   = "GWL 1.5 °C",
             "gwl2.0"   = "GWL 2.0 °C",
             "gwl2.5"   = "GWL 2.5 °C",
             "gwl3.0"   = "GWL 3.0 °C")
cols    <- c("ref91-20" = "grey55",               # single years; the mean is drawn in black
             "gwl1.5"   = "#d4a537",
             "gwl2.0"   = "#e07b39",
             "gwl2.5"   = "#c8432b",
             "gwl3.0"   = "#7f1d1d")
alpha_years <- 0.12                               # transparency of the single years (many lines)

# ---- 2. Inventory of the PCWD files: grid_ch_pcwd_<model chain>_<period>.nc ---------
files <- list.files(here("data", "pcwd"), pattern = "^grid_ch_pcwd_.*\\.nc$", full.names = TRUE)
if (length(files) == 0) stop("No PCWD files in ", here("data", "pcwd"))
parts <- do.call(rbind, strsplit(sub("\\.nc$", "", basename(files)), "_"))
inv   <- data.frame(path = files, chain = parts[, 4], period = parts[, 5])
inv   <- inv[inv$period %in% periods, ]
print(table(inv$period))                          # 26 model chains per period (24 for GWL 3.0)

# ---- 3. Daily PCWD of the cell that contains the site, one line per year -------------
iE <- which.min(abs(E_vals - site_E))             # E_vals, N_vals are cell centres
iN <- which.min(abs(N_vals - site_N))
n_years <- n_days / 365

# day 366 = 1 January of the following year, so that each line closes the circle
# (not available after the last year, that line ends on 31 December)
read_years <- function(i) {
  nc   <- nc_open(inv$path[i])
  pcwd <- ncvar_get(nc, "pcwd", start = c(iE, iN, 1), count = c(1, 1, -1))
  nc_close(nc)
  stopifnot(length(pcwd) == n_days, !anyNA(pcwd))
  do.call(rbind, lapply(seq_len(n_years), function(y) {
    idx <- (y - 1) * 365 + 1:366
    idx <- idx[idx <= n_days]
    data.frame(chain = inv$chain[i], period = inv$period[i], year = y,
               doy = seq_along(idx), pcwd = pcwd[idx])
  }))
}
d <- do.call(rbind, lapply(seq_len(nrow(inv)), read_years))
d$line <- paste(d$chain, d$year)                  # one line per model chain and year

# mean over all years and model chains per day of year, day 366 = mean of 1 January
d_mean <- aggregate(pcwd ~ period + doy, d[d$doy <= 365, ], mean)
d_mean <- rbind(d_mean, transform(d_mean[d_mean$doy == 1, ], doy = 366))

# panel titles with the number of model chains and years
n_ch  <- tapply(inv$chain, inv$period, function(x) length(unique(x)))
strip <- setNames(sprintf("%s (%d model chains, %d years)", labels[names(n_ch)], n_ch, n_ch * n_years),
                  names(n_ch))
strip <- strip[intersect(periods, names(strip))]  # panel order: reference first
d$panel      <- factor(strip[d$period], levels = strip)
d_mean$panel <- factor(strip[d_mean$period], levels = strip)
r_max <- max(d$pcwd)                              # same radial scale in all figures

# ---- 4. Plot function ----------------------------------------------------------
month_start <- cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))   # doy of the 1st

plot_circle <- function(dd, dm) {
  ggplot(dd, aes(doy, pcwd)) +
    geom_line(aes(group = line, colour = cols[period]), linewidth = 0.25, alpha = alpha_years) +
    geom_line(data = dm, colour = "black", linewidth = 1) +
    coord_polar(theta = "x", start = 0, direction = 1) +        # 1 January at the top, clockwise
    scale_colour_identity() +
    scale_x_continuous(breaks = month_start, labels = month.abb, limits = c(1, 366)) +
    scale_y_continuous(limits = c(0, r_max * 1.02)) +
    labs(x = NULL, y = "PCWD (mm)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey88"),
          axis.text.x      = element_text(size = 10, colour = "grey30"),
          strip.text       = element_text(face = "bold", size = 12),
          plot.title       = element_text(face = "bold", size = 12, hjust = 0.5))
}

dir.create(here("figures"), showWarnings = FALSE)

# ---- 5. One figure per period ---------------------------------------------------
for (per in names(strip)) {
  p <- plot_circle(d[d$period == per, ], d_mean[d_mean$period == per, ]) + labs(title = strip[[per]])
  ggsave(here("figures", sprintf("pcwd_doy_bern_circle_%s.jpg", per)),
         p, width = 6.5, height = 6.5, dpi = 300, bg = "white")
  message("Saved figure for ", per)
}

# ---- 6. All periods stacked ---------------------------------------------------
p_all <- plot_circle(d, d_mean) + facet_wrap(~panel, ncol = 1)
print(p_all)
ggsave(here("figures", "pcwd_doy_bern_circle_all.jpg"),
       p_all, width = 6.5, height = 6 * length(strip) + 0.5, dpi = 300, bg = "white", limitsize = FALSE)

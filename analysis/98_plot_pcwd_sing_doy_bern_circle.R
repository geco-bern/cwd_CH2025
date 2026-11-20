# run 07_apply_pcwd_UBELIX.R first to create the daily PCWD grids
# daily PCWD for the grid cell of Bern-Zollikofen as a circle plot: one full turn = one year
# (1 January at the top, clockwise), distance from the centre = PCWD;
# reference period as 30-year mean per day of year (blue), GWL 2.0 as three single years (orange)
# the plot is shown and saved to figures/

library(here)
library(ncdf4)
library(terra)                                    # needed by R/nc_grid_ch.R
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout (E_vals, N_vals, n_days)

# ---- 1. Settings ----------------------------------------------------------
model    <- "clmcom-cclm4-cccma"
site     <- "Bern-Zollikofen"
site_E   <- 2601934                               # CH1903+ / LV95
site_N   <- 1204410

ref      <- "ref91-20"
gwl      <- "gwl2.0"
gwl_lab  <- "GWL 2.0 °C"
gwl_yrs  <- c(1, 15, 30)                          # single years of the GWL period to show

col_ref  <- "#0d366b"                             # reference mean
col_yrs  <- c("#f0a27a", "#d0612c", "#7a2a0c")    # GWL years, light to dark

pcwd_file <- function(period) here("data", "pcwd", sprintf("grid_ch_pcwd_%s_%s.nc", model, period))

# ---- 2. Daily PCWD of the cell that contains the site --------------------------
iE <- which.min(abs(E_vals - site_E))             # E_vals, N_vals are cell centres
iN <- which.min(abs(N_vals - site_N))
n_years <- n_days / 365

read_pcwd <- function(period) {
  nc   <- nc_open(pcwd_file(period))
  pcwd <- ncvar_get(nc, "pcwd", start = c(iE, iN, 1), count = c(1, 1, -1))
  nc_close(nc)
  stopifnot(length(pcwd) == n_days, !anyNA(pcwd))
  pcwd
}

# day 366 = 1 January of the next year, so that the line closes the circle
# reference: 30-year mean per day of year (365_day calendar), day 366 = mean of 1 January
p_ref  <- rowMeans(matrix(read_pcwd(ref), nrow = 365, ncol = n_years))
d_ref  <- data.frame(series = sprintf("Reference 1991-2020, 30-year mean"),
                     doy = 1:366, pcwd = c(p_ref, p_ref[1]))

# GWL: single years, day 366 = 1 January of the following year (not available after the last year)
p_gwl  <- read_pcwd(gwl)
d_gwl  <- do.call(rbind, lapply(gwl_yrs, function(y) {
  idx <- (y - 1) * 365 + 1:366
  idx <- idx[idx <= n_days]
  data.frame(series = sprintf("%s, year %d", gwl_lab, y), doy = seq_along(idx), pcwd = p_gwl[idx])
}))

d <- rbind(d_ref, d_gwl)
series_cols <- setNames(c(col_ref, col_yrs), c(d_ref$series[1], unique(d_gwl$series)))
d$series    <- factor(d$series, levels = names(series_cols))

# ---- 3. Plot ---------------------------------------------------------------
month_start <- cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))   # doy of the 1st

p <- ggplot(d, aes(doy, pcwd, colour = series)) +
  geom_line(aes(linewidth = series)) +
  coord_polar(theta = "x", start = 0, direction = 1) +        # 1 January at the top, clockwise
  scale_x_continuous(breaks = month_start, labels = month.abb, limits = c(1, 366)) +
  scale_y_continuous(limits = c(0, max(d$pcwd) * 1.05)) +
  scale_colour_manual(values = series_cols, name = NULL) +
  scale_linewidth_manual(values = setNames(c(1.2, 0.6, 0.6, 0.6), names(series_cols)), name = NULL) +
  guides(colour = guide_legend(ncol = 1), linewidth = guide_legend(ncol = 1)) +
  labs(title    = sprintf("Daily PCWD at %s, %s", site, toupper(model)),
       subtitle = sprintf("Reference 1991-2020 (30-year mean) and %s (single years %s)\nGrid cell E %d / N %d (LV95), distance from the centre = PCWD (mm)",
                          gwl_lab, paste(gwl_yrs, collapse = ", "), E_vals[iE], N_vals[iN]),
       x = NULL, y = "PCWD (mm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey88"),
        axis.text.x      = element_text(size = 10, colour = "grey30"),
        plot.subtitle    = element_text(colour = "#52514e", size = 9),
        legend.position  = "bottom",
        legend.text      = element_text(size = 9))

print(p)

# ---- 4. Save ---------------------------------------------------------------
dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", sprintf("pcwd_doy_bern_circle_%s_%s_%s.jpg", model, ref, gwl)),
       p, width = 7.5, height = 8.5, dpi = 300, bg = "white")

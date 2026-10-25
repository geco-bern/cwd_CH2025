# run 05_apply_pcwd_CH2025.R first (for every period) to create the daily PCWD grids
# daily PCWD against day of year for the grid cell of Bern-Zollikofen,
# one light line per year (all 30 years of each period in one plot) and a dark line for the
# 30-year mean per day of year; reference period in blue, GWL in orange;
# the plot is shown and saved to figures/

library(here)
library(ncdf4)
library(terra)                                    # needed by R/nc_grid_ch.R
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout (E_vals, N_vals, n_days)

# ---- 1. Settings ----------------------------------------------------------
model   <- "clmcom-cclm4-cccma"
site    <- "Bern-Zollikofen"
site_E  <- 2601934                                # CH1903+ / LV95
site_N  <- 1204410

# one row per period: label, colour of the single years (light) and of the mean (dark)
periods <- data.frame(
  period    = c("ref91-20", "gwl1.5"),
  label     = c("Reference 1991-2020", "GWL 1.5"),
  col_years = c("#b7d3f6", "#f5c3a9"),
  col_mean  = c("#0d366b", "#7a2a0c"))

pcwd_file <- function(period) here("data", "pcwd", sprintf("grid_ch_pcwd_%s_%s.nc", model, period))

# ---- 2. Daily PCWD of the cell that contains the site --------------------------
iE <- which.min(abs(E_vals - site_E))             # E_vals, N_vals are cell centres
iN <- which.min(abs(N_vals - site_N))
n_years <- n_days / 365

d <- do.call(rbind, lapply(seq_len(nrow(periods)), function(i) {
  nc   <- nc_open(pcwd_file(periods$period[i]))
  pcwd <- ncvar_get(nc, "pcwd", start = c(iE, iN, 1), count = c(1, 1, -1))
  nc_close(nc)
  stopifnot(length(pcwd) == n_days, !anyNA(pcwd))

  data.frame(period = periods$label[i],
             year   = rep(1:n_years, each = 365),   # 365_day calendar, no Feb 29
             doy    = rep(1:365, times = n_years),
             pcwd   = pcwd)
}))
d$line   <- paste(d$period, d$year)               # one line per period and year
d$series <- paste0(d$period, ", single years")

# 30-year mean per day of year of each period
d_mean <- aggregate(pcwd ~ period + doy, d, mean)
d_mean$series <- paste0(d_mean$period, ", 30-year mean")

# colours and legend order: single years and mean of each period
series_cols <- setNames(c(rbind(periods$col_years, periods$col_mean)),
                        c(rbind(paste0(periods$label, ", single years"),
                                paste0(periods$label, ", 30-year mean"))))

# ---- 3. Plot ---------------------------------------------------------------
month_start <- cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))   # doy of the 1st

# draw the periods in the order of `periods` (reference first), means on top of all years
d$line        <- factor(d$line, levels = unique(d$line))
d$series      <- factor(d$series, levels = names(series_cols))
d_mean$series <- factor(d_mean$series, levels = names(series_cols))

p <- ggplot(d, aes(doy, pcwd, colour = series)) +
  geom_line(aes(group = line), linewidth = 0.3) +
  geom_line(data = d_mean, aes(group = period), linewidth = 1) +
  scale_x_continuous(breaks = month_start, labels = month.abb,
                     limits = c(1, 365), expand = expansion(add = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = series_cols, breaks = names(series_cols), name = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = rep(c(0.6, 1.2), nrow(periods))))) +
  labs(title    = sprintf("Daily PCWD at %s, %s, %s", site, toupper(model),
                          paste(periods$label, collapse = " and ")),
       subtitle = sprintf("Grid cell E %d / N %d (LV95), all 30 years of each period and their mean per day of year, deficit reset on 1 January",
                          E_vals[iE], N_vals[iN]),
       x = NULL, y = "PCWD (mm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        panel.grid.major.y = element_line(colour = "grey92"),
        plot.subtitle      = element_text(colour = "#52514e"),
        legend.text        = element_text(size = 9))

print(p)

# ---- 4. Save ---------------------------------------------------------------
dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", sprintf("pcwd_doy_bern_%s_%s.jpg", model, paste(periods$period, collapse = "_"))),
       p, width = 10, height = 5.5, dpi = 300, bg = "white")

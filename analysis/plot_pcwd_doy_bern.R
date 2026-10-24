# run 05_apply_pcwd_CH2025.R first to create the daily PCWD grid
# daily PCWD against day of year for the grid cell of Bern-Zollikofen,
# one line per year (all 30 years in one plot); the plot is only shown, not saved

library(here)
library(ncdf4)
library(terra)                                    # needed by R/nc_grid_ch.R
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout (E_vals, N_vals, n_days)

# ---- 1. Settings ----------------------------------------------------------
model   <- "clmcom-cclm4-cccma"
period  <- "gwl1.5"
site    <- "Bern-Zollikofen"
site_E  <- 2601934                                # CH1903+ / LV95
site_N  <- 1204410

pcwd_file <- here("data", "pcwd", sprintf("grid_ch_pcwd_%s_%s.nc", model, period))

# ---- 2. Daily PCWD of the cell that contains the site --------------------------
iE <- which.min(abs(E_vals - site_E))             # E_vals, N_vals are cell centres
iN <- which.min(abs(N_vals - site_N))

nc   <- nc_open(pcwd_file)
pcwd <- ncvar_get(nc, "pcwd", start = c(iE, iN, 1), count = c(1, 1, -1))
nc_close(nc)
stopifnot(length(pcwd) == n_days, !anyNA(pcwd))

d <- data.frame(year = rep(1:(n_days / 365), each = 365),   # 365_day calendar, no Feb 29
                doy  = rep(1:365, times = n_days / 365),
                pcwd = pcwd)

# ---- 3. Plot ---------------------------------------------------------------
month_start <- cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))   # doy of the 1st

p <- ggplot(d, aes(doy, pcwd, colour = year, group = year)) +
  geom_line(linewidth = 0.4) +
  scale_x_continuous(breaks = month_start, labels = month.abb,
                     limits = c(1, 365), expand = expansion(add = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_gradient(low = "#86b6ef", high = "#0d366b",
                        breaks = c(1, 10, 20, 30),
                        name = "Year\n(1 = 2025)") +
  labs(title    = sprintf("Daily PCWD at %s, %s %s", site, toupper(model), toupper(period)),
       subtitle = sprintf("Grid cell E %d / N %d (LV95), all 30 years, deficit reset on 1 January",
                          E_vals[iE], N_vals[iN]),
       x = NULL, y = "PCWD (mm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        panel.grid.major.y = element_line(colour = "grey92"),
        plot.subtitle      = element_text(colour = "#52514e"),
        legend.title       = element_text(size = 9, colour = "#52514e"))

print(p)

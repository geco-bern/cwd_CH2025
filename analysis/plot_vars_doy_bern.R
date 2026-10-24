# run 04_simulate_Rn.R and 05_apply_pcwd_CH2025.R first
# PCWD, net radiation, mean temperature and precipitation against day of year for the
# grid cell of Bern-Zollikofen: one panel per variable, one line per year (all 30 years),
# black line = mean over the 30 years for each day of year; the plot is only shown, not saved

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

dir_grid <- here("data-raw", "CH2025", "DAILY_GRIDDED")
vars <- data.frame(
  var   = c("pcwd", "rn", "tas", "pr"),
  label = c("PCWD (mm)", "Net radiation (W m-2)", "Mean temperature (°C)", "Precipitation (mm day-1)"),
  file  = c(here("data", "pcwd", sprintf("grid_ch_pcwd_%s_%s.nc", model, period)),
            here("data", "netradiation_sim", sprintf("grid_ch_rn_%s_%s.nc", model, period)),
            file.path(dir_grid, sprintf("ogd-climate-scenarios-ch2025-grid_ch_tas_%s_%s.nc", model, period)),
            file.path(dir_grid, sprintf("ogd-climate-scenarios-ch2025-grid_ch_pr_%s_%s.nc", model, period))))

# ---- 2. Daily values of the cell that contains the site --------------------------
iE <- which.min(abs(E_vals - site_E))             # all files: E x N (south to north) x time
iN <- which.min(abs(N_vals - site_N))

read_cell <- function(file, var) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  x <- ncvar_get(nc, var, start = c(iE, iN, 1), count = c(1, 1, -1))
  stopifnot(length(x) == n_days, !anyNA(x))
  x
}

d <- do.call(rbind, lapply(seq_len(nrow(vars)), function(i)
  data.frame(variable = vars$label[i],
             year     = rep(1:(n_days / 365), each = 365),   # 365_day calendar, no Feb 29
             doy      = rep(1:365, times = n_days / 365),
             value    = read_cell(vars$file[i], vars$var[i]))))
d$variable <- factor(d$variable, levels = vars$label)          # panel order as in vars

d_mean <- aggregate(value ~ variable + doy, d, mean)           # 30-year mean per day of year

# ---- 3. Plot ---------------------------------------------------------------
month_start <- cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))   # doy of the 1st

p <- ggplot(d, aes(doy, value)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_line(aes(colour = year, group = year), linewidth = 0.3, alpha = 0.7) +
  geom_line(data = d_mean, colour = "#0b0b0b", linewidth = 0.6) +
  facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_x_continuous(breaks = month_start, labels = month.abb,
                     limits = c(1, 365), expand = expansion(add = 2)) +
  scale_colour_gradient(low = "#86b6ef", high = "#0d366b",
                        breaks = c(1, 10, 20, 30),
                        name = "Year\n(1 = 2025)") +
  labs(title    = sprintf("Daily values at %s, %s %s", site, toupper(model), toupper(period)),
       subtitle = sprintf("Grid cell E %d / N %d (LV95), one line per year, black = 30-year mean per day of year",
                          E_vals[iE], N_vals[iN]),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        strip.text       = element_text(hjust = 0, face = "bold"),
        plot.subtitle    = element_text(colour = "#52514e"),
        legend.title     = element_text(size = 9, colour = "#52514e"),
        legend.key.height = unit(1.6, "cm"))

print(p)

# run 04_simulate_Rn.R and 05_apply_pcwd_CH2025.R first (for every period)
# monthly means of PCWD, net radiation, mean temperature and precipitation for the
# grid cell of Bern-Zollikofen: one panel per variable, one light line per year
# (all 30 years of each period) and a dark line for the 30-year mean per month;
# reference period in blue, GWL in orange (as in plot_pcwd_doy_bern.R);
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

# one row per variable (panel order); units of the monthly means
vars <- data.frame(
  var   = c("pcwd", "rn", "tas", "pr"),
  label = c("PCWD (mm)", "Net radiation (W m-2)", "Mean temperature (°C)", "Precipitation (mm day-1)"))

dir_grid <- here("data-raw", "CH2025", "DAILY_GRIDDED")
var_file <- function(var, period) switch(var,
  pcwd = here("data", "pcwd", sprintf("grid_ch_pcwd_%s_%s.nc", model, period)),
  rn   = here("data", "netradiation_sim", sprintf("grid_ch_rn_%s_%s.nc", model, period)),
  file.path(dir_grid, sprintf("ogd-climate-scenarios-ch2025-grid_ch_%s_%s_%s.nc", var, model, period)))

# ---- 2. Monthly means of the cell that contains the site --------------------------
iE <- which.min(abs(E_vals - site_E))             # all files: E x N (south to north) x time
iN <- which.min(abs(N_vals - site_N))

n_years <- n_days / 365
month   <- rep(1:12, times = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))  # 365_day calendar

read_cell <- function(file, var) {
  nc <- nc_open(file)
  on.exit(nc_close(nc))
  x <- ncvar_get(nc, var, start = c(iE, iN, 1), count = c(1, 1, -1))
  stopifnot(length(x) == n_days, !anyNA(x))
  x
}

d <- do.call(rbind, lapply(seq_len(nrow(periods)), function(i)
  do.call(rbind, lapply(seq_len(nrow(vars)), function(j) {
    daily <- data.frame(year  = rep(1:n_years, each = 365),
                        month = rep(month, times = n_years),
                        value = read_cell(var_file(vars$var[j], periods$period[i]), vars$var[j]))
    m <- aggregate(value ~ year + month, daily, mean)       # mean of each month of each year
    cbind(period = periods$label[i], variable = vars$label[j], m)
  }))))

d_mean <- aggregate(value ~ period + variable + month, d, mean)   # 30-year mean per month

# panel order as in vars; colours and legend order: single years and mean of each period
series_cols <- setNames(c(rbind(periods$col_years, periods$col_mean)),
                        c(rbind(paste0(periods$label, ", single years"),
                                paste0(periods$label, ", 30-year mean"))))
d$variable      <- factor(d$variable, levels = vars$label)
d_mean$variable <- factor(d_mean$variable, levels = vars$label)
d$series        <- factor(paste0(d$period, ", single years"), levels = names(series_cols))
d_mean$series   <- factor(paste0(d_mean$period, ", 30-year mean"), levels = names(series_cols))

# draw the periods in the order of `periods` (reference first), means on top of all years
d$line <- factor(paste(d$period, d$year), levels = unique(paste(rep(periods$label, each = n_years),
                                                                 1:n_years)))

# ---- 3. Plot ---------------------------------------------------------------
p <- ggplot(d, aes(month, value, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_line(aes(group = line), linewidth = 0.3) +
  geom_line(data = d_mean, aes(group = period), linewidth = 1) +
  facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_x_continuous(breaks = 1:12, labels = month.abb, expand = expansion(add = 0.2)) +
  scale_colour_manual(values = series_cols, breaks = names(series_cols), name = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = rep(c(0.6, 1.2), nrow(periods))))) +
  labs(title    = sprintf("Monthly means at %s, %s, %s", site, toupper(model),
                          paste(periods$label, collapse = " and ")),
       subtitle = sprintf("Grid cell E %d / N %d (LV95), one line per year, thick line = 30-year mean per month",
                          E_vals[iE], N_vals[iN]),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        strip.text       = element_text(hjust = 0, face = "bold"),
        plot.subtitle    = element_text(colour = "#52514e"),
        legend.position  = "top",
        legend.text      = element_text(size = 9))

print(p)

# ---- 4. Save ---------------------------------------------------------------
dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", sprintf("vars_monthly_bern_%s_%s.jpg", model, paste(periods$period, collapse = "_"))),
       p, width = 9, height = 11, dpi = 300, bg = "white")

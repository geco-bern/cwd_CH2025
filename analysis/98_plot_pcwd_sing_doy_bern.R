# daily PCWD against day of year for the grid cell of Bern-Zollikofen, one GWL period of one model chain:
# one line per year, all years in grey except a few highlighted years in different colours,
# labelled at their annual maximum (needs the ggrepel package)
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

period  <- "gwl2.0"
label   <- "GWL 2.0 °C"

# highlighted years and their colours (Okabe-Ito palette, distinguishable with colour vision deficiency)
hl_years <- c(1, 5, 10, 15, 20, 25, 30)
hl_cols  <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000")
col_other <- "grey80"                             # all other years

pcwd_file <- here("data", "pcwd", sprintf("grid_ch_pcwd_%s_%s.nc", model, period))

# ---- 2. Daily PCWD of the cell that contains the site --------------------------
iE <- which.min(abs(E_vals - site_E))             # E_vals, N_vals are cell centres
iN <- which.min(abs(N_vals - site_N))
n_years <- n_days / 365

nc   <- nc_open(pcwd_file)
pcwd <- ncvar_get(nc, "pcwd", start = c(iE, iN, 1), count = c(1, 1, -1))
nc_close(nc)
stopifnot(length(pcwd) == n_days, !anyNA(pcwd))

d <- data.frame(year = rep(1:n_years, each = 365),   # 365_day calendar, no Feb 29
                doy  = rep(1:365, times = n_years),
                pcwd = pcwd)

# series: "Year 1", "Year 5", ... for the highlighted years, "Other years" for all others
series_cols <- c(setNames(hl_cols, paste("Year", hl_years)), "Other years" = col_other)
d$series <- ifelse(d$year %in% hl_years, paste("Year", d$year), "Other years")
d$series <- factor(d$series, levels = names(series_cols))

# labels of the highlighted years, placed at the annual maximum of each year
d_lab <- do.call(rbind, lapply(hl_years, function(y) {
  x <- d[d$year == y, ]
  x[which.max(x$pcwd), ]
}))

# ---- 3. Plot ---------------------------------------------------------------
month_start <- cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))   # doy of the 1st

p <- ggplot(d, aes(doy, pcwd, colour = series, group = year)) +
  geom_line(data = d[!d$year %in% hl_years, ], linewidth = 0.3) +        # grey years below
  geom_line(data = d[d$year %in% hl_years, ], linewidth = 0.7) +         # highlighted years on top
  geom_point(data = d_lab, size = 1.5) +
  ggrepel::geom_text_repel(data = d_lab, aes(label = series), size = 3.5, fontface = "bold",
                           nudge_y = 0.06 * max(d$pcwd), min.segment.length = 0,
                           segment.size = 0.3, box.padding = 0.4, seed = 1,
                           bg.colour = "white", bg.r = 0.15) +                  # white outline

  scale_x_continuous(breaks = month_start, labels = month.abb,
                     limits = c(1, 365), expand = expansion(add = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +            # room for the labels
  scale_colour_manual(values = series_cols) +
  labs(title = sprintf("Daily PCWD at %s, %s, %s", site, toupper(model), label),
       x = NULL, y = "PCWD (mm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        panel.grid.major.y = element_line(colour = "grey92"),
        legend.position    = "none")

print(p)

# ---- 4. Save ---------------------------------------------------------------
dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", sprintf("pcwd_doy_bern_%s_%s_years.jpg", model, period)),
       p, width = 10, height = 5.5, dpi = 300, bg = "white")

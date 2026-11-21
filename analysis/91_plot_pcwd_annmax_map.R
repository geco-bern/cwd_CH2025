# run 08_get_pcwd_annmax_UBELIX.R first to create the annual maximum PCWD grids
# maps of the mean annual maximum PCWD per grid cell, ensemble mean over all model chains
# (same layout as 93_plot_pcwd_mean_map.R, but annual only: one panel per period, one row)
#   1. per file (model chain x period): mean of the annual maximum PCWD over the 30 years, per cell
#   2. per period: mean over the model chains, per cell
# one shared binned colour scale, the same as in 93 (R/pcwd_colours.R); the figure is saved to figures/

library(here); library(ncdf4); library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout (E_vals, N_vals)
source(here("R", "pcwd_colours.R"))               # pcwd_bin(), scale_fill_pcwd()

periods <- c("ref91-20" = "Reference\n1991-2020", "gwl1.5" = "GWL 1.5 °C", "gwl2.0" = "GWL 2.0 °C",
             "gwl2.5" = "GWL 2.5 °C", "gwl3.0" = "GWL 3.0 °C")

# ---- 1. Ensemble mean maps (one value per cell and period) -------------------------
files <- list.files(here("data", "pcwd_annmax"), pattern = "^grid_ch_pcwd_annmax_.*\\.nc$", full.names = TRUE)
per   <- sapply(strsplit(sub("\\.nc$", "", basename(files)), "_"), `[`, 6)
files <- files[per %in% names(periods)]; per <- per[per %in% names(periods)]
if (length(files) == 0) stop("No annual maximum files in ", here("data", "pcwd_annmax"))
print(table(per))                                 # 26 model chains per period (24 for GWL 3.0)

# mean annual maximum over the 30 years per cell of one file (the files are small, so read at once);
# cells outside the domain are NA on every year and stay NA
cell_means <- function(f) {
  nc <- nc_open(f); on.exit(nc_close(nc))
  x  <- ncvar_get(nc, "pcwd_annmax", collapse_degen = FALSE)   # E x N x years
  rowMeans(matrix(x, ncol = dim(x)[3]))                        # cells (E fastest, N south to north)
}
res <- lapply(files, cell_means)
ens <- lapply(split(res, per), function(r) Reduce(`+`, r) / length(r))   # mean over the model chains

# ---- 2. Long table, colour classes, panel letters -------------------------------
xy <- expand.grid(x = E_vals, y = N_vals)           # cell order of the files (E fastest, N south to north)
d  <- do.call(rbind, lapply(intersect(names(periods), names(ens)), function(p)
  data.frame(xy, period = p, pcwd = ens[[p]])))
d  <- d[!is.na(d$pcwd), ]
d$period <- factor(periods[d$period], levels = periods[names(periods) %in% names(ens)])

d$bin <- pcwd_bin(d$pcwd)                          # colour classes of R/pcwd_colours.R

pan <- data.frame(period = factor(levels(d$period), levels = levels(d$period)))   # (a), (b), ... left to right
pan$letter <- sprintf("(%s)", letters[seq_len(nrow(pan))])

x_lim <- range(d$x) + c(-500, 500)                  # extent of the cells with data
y_lim <- range(d$y) + c(-500, 500)

# ---- 3. Plot ---------------------------------------------------------------
p <- ggplot(d, aes(x, y, fill = bin)) +
  geom_raster(show.legend = TRUE) +                  # TRUE: also show classes without cells
  geom_text(data = pan, aes(x_lim[1], y_lim[2], label = letter), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  facet_wrap(~period, nrow = 1) +
  coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_pcwd(d$bin, "Annual maximum PCWD (mm)") +
  # legend in one row below the maps (a vertical legend is taller than one row of maps)
  guides(fill = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5, label.position = "bottom",
                             keywidth = unit(1.5, "cm"), keyheight = unit(0.4, "cm"))) +
  theme_void(base_size = 12) +
  theme(strip.text      = element_text(size = 13, margin = margin(b = 4)),
        panel.spacing   = unit(0.4, "lines"),
        legend.position = "bottom",
        legend.text     = element_text(size = 9),
        plot.margin     = margin(10, 10, 10, 10),
        plot.background = element_rect(fill = "white", colour = NA))

print(p)

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "pcwd_annmax_map_ch.jpg"), p, width = 16, height = 4.8, dpi = 300, bg = "white")
message("Saved figure")

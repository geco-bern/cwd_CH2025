# run 09_get_pcwd_maxlen_UBELIX.R first to create the yearly maximum PCWD event lengths of all model chains and periods
# mean maximum PCWD event length over Switzerland per model chain and period:
#   1. per file (model chain x period): mean of the yearly maximum event length over all cells of the domain
#      and the 30 years
#   2. per period: ensemble mean and spread over the model chains
# one panel for the complete year, in the style of the CH2025 figures (same layout as 96_plot_pcwd_seasons.R):
# box = spread of the model chains, line = ensemble mean; reference grey, GWLs yellow to dark red
# second figure: violin per period over all years and model chains (same layout as 94_plot_pcwd_annual.R)
# third figure: the violins once per biogeographical region (R/biogeoregio_ch.R)
# all figures are saved to figures/

library(here)
library(ncdf4)
library(terra)                                    # needed by R/biogeoregio_ch.R
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout, needed by biogeoregio_nc()
source(here("R", "biogeoregio_ch.R"))             # biogeoregio_nc(), by_region()

# ---- 1. Settings ----------------------------------------------------------
spread  <- c(0, 1)                                # box: quantiles of the model chains (0, 1 = full range;
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

# ---- 2. Inventory of the files: grid_ch_pcwd_maxlen_<model chain>_<period>.nc --------
files <- list.files(here("data", "pcwd_maxlen"), pattern = "^grid_ch_pcwd_maxlen_.*\\.nc$", full.names = TRUE)
stopifnot(length(files) > 0)
parts <- do.call(rbind, strsplit(sub("\\.nc$", "", basename(files)), "_"))
inv   <- data.frame(path = files, chain = parts[, 5], period = parts[, 6])
inv   <- inv[inv$period %in% periods, ]
print(table(inv$period))                          # 26 model chains per period (24 for GWL 3.0)

# ---- 3. Mean over the domain per year, one file at a time --------------------------
# the files are small (one value per cell and year), so each is read at once;
# cells outside the domain are NA and not counted, so every cell of the domain has the same weight
read_maxlen <- function(i) {
  nc <- nc_open(inv$path[i])
  on.exit(nc_close(nc))
  x <- ncvar_get(nc, "pcwd_maxlen", collapse_degen = FALSE)   # E x N x years
  x <- matrix(x, ncol = dim(x)[3])                            # cells x years
  message("Done: ", basename(inv$path[i]))
  data.frame(chain = inv$chain[i], period = inv$period[i], year = seq_len(ncol(x)),
             maxlen = colMeans(x, na.rm = TRUE))
}
d <- do.call(rbind, lapply(seq_len(nrow(inv)), read_maxlen))

# mean over the 30 years per model chain (the domain is the same in every year,
# so this equals the mean over all cells and years)
chains <- aggregate(maxlen ~ chain + period, d, mean)

# ---- 4. Ensemble statistics over the model chains, per period ----------------------
ens <- do.call(rbind, lapply(split(chains, chains$period), function(x)
  data.frame(period = x$period[1], n_chains = nrow(x),
             mean = mean(x$maxlen),
             lo   = quantile(x$maxlen, spread[1], names = FALSE),
             hi   = quantile(x$maxlen, spread[2], names = FALSE))))
ens <- ens[order(match(ens$period, periods)), ]
print(ens, row.names = FALSE)

ens$x   <- match(ens$period, periods)             # x position of each period
ens$col <- cols[ens$period]

# ---- 5. Plot: one panel for the complete year --------------------------------------
w <- 0.32                                         # half width of the boxes

p <- ggplot(ens) +
  geom_rect(aes(xmin = x - w, xmax = x + w, ymin = lo, ymax = hi, fill = col), alpha = 0.45) +  # spread
  geom_segment(aes(x = x - w, xend = x + w, y = mean, yend = mean, colour = col),
               linewidth = 1.3) +                                                              # mean
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(breaks = seq_along(periods), labels = labels[periods],
                     limits = c(0.5, length(periods) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(150, ens$hi)), breaks = seq(0, 360, by = 30),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Max. PCWD event length (days)") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey65", linetype = "dashed", linewidth = 0.4),
        axis.title.y       = element_text(face = "bold"),
        axis.text          = element_text(colour = "black"))

print(p)

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "pcwd_maxlen_ch.jpg"), p, width = 6, height = 4.5, dpi = 300, bg = "white")

# ---- 6. Violin plot per period over all years and model chains (as in 94_plot_pcwd_annual.R) ----
# per period all 30 years of all model chains are pooled (the placeholder years 0001-0030 are
# different simulated years in every model chain, so they are not compared year by year)
# outline: distribution of all years and model chains (kernel density estimate), all violins
# with the same maximum width (scale = "width") and cut at the smallest and largest value (trim);
# inner box: median and interquartile range; dashed line: median of the reference period
labels_v <- c("ref91-20" = "Reference\n1991-2020",
              "gwl1.5"   = "GWL 1.5 °C",
              "gwl2.0"   = "GWL 2.0 °C",
              "gwl2.5"   = "GWL 2.5 °C",
              "gwl3.0"   = "GWL 3.0 °C")
d$period_lab <- factor(labels_v[d$period], levels = labels_v[periods])
n_val   <- aggregate(maxlen ~ period_lab, d, length)                # values per violin
ref_med <- median(d$maxlen[d$period == "ref91-20"])

p_violin <- ggplot(d, aes(period_lab, maxlen, fill = cols[period], colour = cols[period])) +
  geom_hline(yintercept = ref_med, colour = cols[["ref91-20"]], linetype = "dashed", linewidth = 0.6) +
  geom_violin(alpha = 0.45, linewidth = 0.6, width = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1, linewidth = 0.5, outlier.shape = NA) +
  geom_text(data = n_val, aes(period_lab, 0, label = paste0("n = ", maxlen)),
            inherit.aes = FALSE, vjust = -0.5, size = 3.2, colour = "grey30") +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 360, by = 30),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Max. PCWD event length (days)") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey65", linetype = "dashed", linewidth = 0.4),
        axis.title.y       = element_text(face = "bold"),
        axis.text          = element_text(colour = "black"))

print(p_violin)
ggsave(here("figures", "pcwd_maxlen_violin_ch.jpg"), p_violin, width = 7, height = 5, dpi = 300, bg = "white")

# ---- 7. Mean maximum event length per biogeographical region, per year --------------------
# same as section 3, but the mean is taken over the cells of each region instead of the whole domain;
# every cell takes the region that contains its centre (cells outside Switzerland are not counted)
read_maxlen_rast <- function(path) {
  r <- rast(path, subds = "pcwd_maxlen")
  crs(r) <- "EPSG:2056"                           # same CRS as the regions, also if it is not read from the file
  r
}
zones <- biogeoregio_nc()                         # region of every cell, from data/biogeoregio/grid_ch_biogeoregio.nc

read_maxlen_region <- function(i) {
  x <- by_region(read_maxlen_rast(inv$path[i]), zones, "mean")
  data.frame(chain = inv$chain[i], period = inv$period[i], region = x$region,
             n_cells = x$n_cells, year = x$layer, maxlen = x$value)
}
dr <- do.call(rbind, lapply(seq_len(nrow(inv)), read_maxlen_region))
dir.create(here("tables"), showWarnings = FALSE)
write.csv(dr, here("tables", "pcwd_maxlen_years_region_ch.csv"), row.names = FALSE)

# ---- 8. Violin plot per period, one panel per biogeographical region ---------------------
# same as section 6; all panels have the same y axis; dashed line: median of the reference period
# of each region
dr$period_lab <- factor(labels_v[dr$period], levels = labels_v[periods])
n_val_r   <- aggregate(maxlen ~ region + period_lab, dr, length)     # values per violin
ref_med_r <- aggregate(maxlen ~ region, dr[dr$period == "ref91-20", ], median)

p_violin_r <- ggplot(dr, aes(period_lab, maxlen, fill = cols[period], colour = cols[period])) +
  geom_hline(data = ref_med_r, aes(yintercept = maxlen), colour = cols[["ref91-20"]],
             linetype = "dashed", linewidth = 0.6) +
  geom_violin(alpha = 0.45, linewidth = 0.6, width = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1, linewidth = 0.5, outlier.shape = NA) +
  geom_text(data = n_val_r, aes(period_lab, 0, label = paste0("n = ", maxlen)),
            inherit.aes = FALSE, vjust = -0.5, size = 3.2, colour = "grey30") +
  # "free" only to draw the axes in every panel; the limits below keep the same scale in all panels
  facet_wrap(~region, ncol = 2, scales = "free") +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_y_continuous(limits = c(0, max(dr$maxlen)), breaks = seq(0, 360, by = 30),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Max. PCWD event length (days)") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey65", linetype = "dashed", linewidth = 0.4),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title.y       = element_text(face = "bold"),
        axis.text          = element_text(colour = "black"),
        panel.spacing      = unit(1.5, "lines"))

print(p_violin_r)
ggsave(here("figures", "pcwd_maxlen_violin_region_ch.jpg"), p_violin_r,
       width = 11, height = 12, dpi = 300, bg = "white")

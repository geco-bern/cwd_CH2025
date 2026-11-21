# run 08_get_pcwd_annmax_UBELIX.R first to create the annual maximum PCWD grids
# annual maximum PCWD over Switzerland, one boxplot and one violin per period (two figures):
#   1. per file (model chain x period): mean of the annual maximum PCWD over all cells of the domain,
#      for each of the 30 years
#   2. per period: all 30 years of all model chains pooled in one boxplot
#      (the placeholder years 0001-0030 are different simulated years in every model chain,
#      so they are not compared year by year)
# reference grey, GWLs yellow to dark red (as in 96_plot_pcwd_seasons.R);
# dashed line: median of the reference period
# third figure: the violins once per biogeographical region (R/biogeoregio_ch.R)

library(here)
library(ncdf4)
library(terra)                                    # needed by R/biogeoregio_ch.R
library(ggplot2)
source(here("R", "nc_grid_ch.R"))                 # grid layout, needed by biogeoregio_nc()
source(here("R", "biogeoregio_ch.R"))             # biogeoregio_nc(), by_region()

# ---- 1. Settings ----------------------------------------------------------
periods <- c("ref91-20", "gwl1.5", "gwl2.0", "gwl2.5", "gwl3.0")
labels  <- c("ref91-20" = "Reference\n1991-2020",
             "gwl1.5"   = "GWL 1.5 °C",
             "gwl2.0"   = "GWL 2.0 °C",
             "gwl2.5"   = "GWL 2.5 °C",
             "gwl3.0"   = "GWL 3.0 °C")
cols    <- c("ref91-20" = "grey55",               # reference grey, GWLs yellow to dark red
             "gwl1.5"   = "#d4a537",
             "gwl2.0"   = "#e07b39",
             "gwl2.5"   = "#c8432b",
             "gwl3.0"   = "#7f1d1d")

# ---- 2. Inventory of the annual maximum files: grid_ch_pcwd_annmax_<model chain>_<period>.nc ----
files <- list.files(here("data", "pcwd_annmax"), pattern = "^grid_ch_pcwd_annmax_.*\\.nc$", full.names = TRUE)
if (length(files) == 0) stop("No annual maximum files in ", here("data", "pcwd_annmax"),
                             "\nCheck that here() points to cwd_CH2025 and that 08_get_pcwd_annmax_UBELIX.R has run.")
parts <- do.call(rbind, strsplit(sub("\\.nc$", "", basename(files)), "_"))
inv   <- data.frame(path = files, chain = parts[, 5], period = parts[, 6])
inv   <- inv[inv$period %in% periods, ]
print(table(inv$period))                          # 26 model chains per period (24 for GWL 3.0)

# ---- 3. Mean annual maximum over the domain, per year -------------------------------
# cells outside the domain are NA and not counted, so every cell of the domain has the same weight
read_annmax <- function(i) {
  nc <- nc_open(inv$path[i])
  x  <- ncvar_get(nc, "pcwd_annmax", collapse_degen = FALSE)   # E x N x years
  nc_close(nc)
  x  <- matrix(x, ncol = dim(x)[3])                            # cells x years
  data.frame(chain = inv$chain[i], period = inv$period[i], year = seq_len(ncol(x)),
             pcwd = colMeans(x, na.rm = TRUE))
}
d <- do.call(rbind, lapply(seq_len(nrow(inv)), read_annmax))

dir.create(here("tables"), showWarnings = FALSE)
write.csv(d, here("tables", "pcwd_annmax_years_ch.csv"), row.names = FALSE)

# ---- 4. Boxplot per period over all years and model chains ---------------------------
d$period_lab <- factor(labels[d$period], levels = labels[periods])
n_val   <- aggregate(pcwd ~ period_lab, d, length)                  # values per box
ref_med <- median(d$pcwd[d$period == "ref91-20"])

p_box <- ggplot(d, aes(period_lab, pcwd, fill = cols[period], colour = cols[period])) +
  geom_hline(yintercept = ref_med, colour = cols[["ref91-20"]], linetype = "dashed", linewidth = 0.6) +
  geom_boxplot(alpha = 0.45, width = 0.6, linewidth = 0.6, outlier.size = 1, outlier.alpha = 0.6) +
  geom_text(data = n_val, aes(period_lab, 0, label = paste0("n = ", pcwd)),
            inherit.aes = FALSE, vjust = -0.5, size = 3.2, colour = "grey30") +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Annual maximum PCWD (mm)") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey65", linetype = "dashed", linewidth = 0.4),
        axis.title.y       = element_text(face = "bold"),
        axis.text          = element_text(colour = "black"))

print(p_box)

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "pcwd_annmax_box_ch.jpg"), p_box, width = 7, height = 5, dpi = 300, bg = "white")

# ---- 5. Violin plot per period over all years and model chains -------------------------
# outline: distribution of all years and model chains (kernel density estimate), all violins
# with the same maximum width (scale = "width") and cut at the smallest and largest value (trim);
# inner box: median and interquartile range
p_violin <- ggplot(d, aes(period_lab, pcwd, fill = cols[period], colour = cols[period])) +
  geom_hline(yintercept = ref_med, colour = cols[["ref91-20"]], linetype = "dashed", linewidth = 0.6) +
  geom_violin(alpha = 0.45, linewidth = 0.6, width = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1, linewidth = 0.5, outlier.shape = NA) +
  geom_text(data = n_val, aes(period_lab, 0, label = paste0("n = ", pcwd)),
            inherit.aes = FALSE, vjust = -0.5, size = 3.2, colour = "grey30") +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Annual maximum PCWD (mm)") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey65", linetype = "dashed", linewidth = 0.4),
        axis.title.y       = element_text(face = "bold"),
        axis.text          = element_text(colour = "black"))

print(p_violin)
ggsave(here("figures", "pcwd_annmax_violin_ch.jpg"), p_violin, width = 7, height = 5, dpi = 300, bg = "white")

# ---- 6. Mean annual maximum per biogeographical region, per year ----------------------
# same as section 3, but the mean is taken over the cells of each region instead of the whole domain;
# every cell takes the region that contains its centre (cells outside Switzerland are not counted)
read_annmax_rast <- function(path) {
  r <- rast(path, subds = "pcwd_annmax")
  crs(r) <- "EPSG:2056"                           # same CRS as the regions, also if it is not read from the file
  r
}
zones <- biogeoregio_nc()                         # region of every cell, from data/grid_ch_biogeoregio.nc

read_annmax_region <- function(i) {
  x <- by_region(read_annmax_rast(inv$path[i]), zones, "mean")
  data.frame(chain = inv$chain[i], period = inv$period[i], region = x$region,
             n_cells = x$n_cells, year = x$layer, pcwd = x$value)
}
dr <- do.call(rbind, lapply(seq_len(nrow(inv)), read_annmax_region))
write.csv(dr, here("tables", "pcwd_annmax_years_region_ch.csv"), row.names = FALSE)

# ---- 7. Violin plot per period, one panel per biogeographical region ---------------------
# same as section 5; all panels have the same y axis; dashed line: median of the reference period
# of each region
dr$period_lab <- factor(labels[dr$period], levels = labels[periods])
n_val_r   <- aggregate(pcwd ~ region + period_lab, dr, length)       # values per violin
ref_med_r <- aggregate(pcwd ~ region, dr[dr$period == "ref91-20", ], median)

p_violin_r <- ggplot(dr, aes(period_lab, pcwd, fill = cols[period], colour = cols[period])) +
  geom_hline(data = ref_med_r, aes(yintercept = pcwd), colour = cols[["ref91-20"]],
             linetype = "dashed", linewidth = 0.6) +
  geom_violin(alpha = 0.45, linewidth = 0.6, width = 0.8, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1, linewidth = 0.5, outlier.shape = NA) +
  geom_text(data = n_val_r, aes(period_lab, 0, label = paste0("n = ", pcwd)),
            inherit.aes = FALSE, vjust = -0.5, size = 3.2, colour = "grey30") +
  # "free" only to draw the axes in every panel; the limits below keep the same scale in all panels
  facet_wrap(~region, ncol = 2, scales = "free") +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_y_continuous(limits = c(0, max(dr$pcwd)), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Annual maximum PCWD (mm)") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey65", linetype = "dashed", linewidth = 0.4),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title.y       = element_text(face = "bold"),
        axis.text          = element_text(colour = "black"),
        panel.spacing      = unit(1.5, "lines"))

print(p_violin_r)
ggsave(here("figures", "pcwd_annmax_violin_region_ch.jpg"), p_violin_r,
       width = 11, height = 12, dpi = 300, bg = "white")

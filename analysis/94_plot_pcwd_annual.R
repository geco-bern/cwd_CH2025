# run 08_get_pcwd_annmax_UBELIX.R first to create the annual maximum PCWD grids
# annual maximum PCWD over Switzerland, one boxplot and one violin per period (two figures):
#   1. per file (model chain x period): mean of the annual maximum PCWD over all cells of the domain,
#      for each of the 30 years
#   2. per period: all 30 years of all model chains pooled in one boxplot
#      (the placeholder years 0001-0030 are different simulated years in every model chain,
#      so they are not compared year by year)
# reference grey, GWLs yellow to dark red (as in 96_plot_pcwd_seasons.R);
# dashed line: median of the reference period

library(here)
library(ncdf4)
library(ggplot2)

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

# run 09_get_pcwd_maxlen_UBELIX.R first to create the yearly maximum PCWD event lengths of all model chains and periods
# mean maximum PCWD event length over Switzerland per model chain and period:
#   1. per file (model chain x period): mean of the yearly maximum event length over all cells of the domain
#      and the 30 years
#   2. per period: ensemble mean and spread over the model chains
# one panel for the complete year, in the style of the CH2025 figures (same layout as 96_plot_pcwd_seasons.R):
# box = spread of the model chains, line = ensemble mean; reference grey, GWLs yellow to dark red
# the figure is saved to figures/

library(here)
library(ncdf4)
library(ggplot2)

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

# ---- 3. Mean over the domain and the 30 years, one file at a time ------------------
# the files are small (one value per cell and year), so each is read at once;
# cells outside the domain are NA and not counted, so every cell of the domain has the same weight
maxlen_mean <- function(i) {
  nc <- nc_open(inv$path[i])
  on.exit(nc_close(nc))
  x <- ncvar_get(nc, "pcwd_maxlen")               # E x N x years
  message("Done: ", basename(inv$path[i]))
  data.frame(chain = inv$chain[i], period = inv$period[i], maxlen = mean(x, na.rm = TRUE))
}

chains <- do.call(rbind, lapply(seq_len(nrow(inv)), maxlen_mean))

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
  scale_y_continuous(limits = c(0, max(ens$hi)), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Max. PCWD event length (days)") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey65", linetype = "dashed", linewidth = 0.4),
        axis.title.y       = element_text(face = "bold"),
        axis.text          = element_text(colour = "black"))

print(p)

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "pcwd_maxlen_ch.jpg"), p, width = 6, height = 4.5, dpi = 300, bg = "white")

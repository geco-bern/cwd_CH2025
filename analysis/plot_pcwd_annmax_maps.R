# run 06_get_pcwd_annmax_CH2025.R first (for every period) to create the annual maximum PCWD grids
# maps of the mean annual maximum PCWD over all 30 years, one panel per period
# (reference period on top), one shared colour scale; the plot is shown and saved to figures/

library(here)
library(terra)
library(ggplot2)

# ---- 1. Settings ----------------------------------------------------------
model   <- "clmcom-cclm4-cccma"
periods <- c("ref91-20", "gwl1.5")                # panel order, top to bottom
labels  <- c("ref91-20" = "Reference period 1991-2020",
             "gwl1.5"   = "GWL 1.5 °C")

annmax_file <- function(period)
  here("data", "pcwd_annmax", sprintf("grid_ch_pcwd_annmax_%s_%s.nc", model, period))

# ---- 2. Mean annual maximum PCWD over the 30 years of each period ------------------
r <- rast(lapply(periods, function(per) {
  ra <- rast(annmax_file(per), subds = "pcwd_annmax")
  stopifnot(nlyr(ra) == 30)
  mean(ra)                                        # NA outside the domain stays NA
}))
names(r) <- labels[periods]

d <- as.data.frame(r, xy = TRUE, na.rm = FALSE, wide = FALSE)   # x, y, layer, values
d <- d[!is.na(d$values), ]
d$layer <- factor(d$layer, levels = names(r))     # panels in the order of periods

# ---- 3. Plot ---------------------------------------------------------------
p <- ggplot(d, aes(x, y, fill = values)) +
  geom_raster() +
  facet_wrap(~layer, ncol = 1) +
  coord_equal(expand = FALSE) +
  scale_fill_gradientn(colours = c("#fde3d6", "#eb6834", "#7a2a0c"),
                       limits = c(0, max(d$values)),
                       name = "PCWD\n(mm)") +
  labs(title    = sprintf("Mean annual maximum PCWD over 30 years, %s", toupper(model)),
       subtitle = "Same colour scale in both panels, 1 km grid (LV95)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid      = element_blank(),
        axis.text       = element_blank(),
        strip.text      = element_text(hjust = 0, face = "bold"),
        plot.subtitle   = element_text(colour = "#52514e"),
        legend.title    = element_text(size = 9, colour = "#52514e"),
        legend.key.height = unit(1.6, "cm"))

print(p)

# ---- 4. Save ---------------------------------------------------------------
dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", sprintf("pcwd_annmax_mean_%s_%s.jpg", model, paste(periods, collapse = "_"))),
       p, width = 7, height = 7.5, dpi = 300, bg = "white")

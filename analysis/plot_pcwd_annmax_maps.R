# run 06_get_pcwd_annmax_CH2025.R first to create the annual maximum PCWD grid
# maps of the annual maximum PCWD in year 1, 15 and 30, one shared colour scale;
# the plot is only shown, not saved

library(here)
library(terra)
library(ggplot2)

# ---- 1. Settings ----------------------------------------------------------
model  <- "clmcom-cclm4-cccma"
period <- "gwl1.5"
years  <- c(1, 15, 30)                            # CH2025 years 1-30 = assumed 2025-2054

annmax_file <- here("data", "pcwd_annmax",
                    sprintf("grid_ch_pcwd_annmax_%s_%s.nc", model, period))

# ---- 2. Annual maximum PCWD of the selected years ------------------------------
r <- rast(annmax_file, subds = "pcwd_annmax")[[years]]
names(r) <- sprintf("Year %d (%d)", years, years + 2024)

d <- as.data.frame(r, xy = TRUE, na.rm = FALSE, wide = FALSE)   # x, y, layer, values
d <- d[!is.na(d$values), ]
d$layer <- factor(d$layer, levels = names(r))     # panels in time order

# ---- 3. Plot ---------------------------------------------------------------
p <- ggplot(d, aes(x, y, fill = values)) +
  geom_raster() +
  facet_wrap(~layer, ncol = 1) +
  coord_equal(expand = FALSE) +
  scale_fill_gradientn(colours = c("#fde3d6", "#eb6834", "#7a2a0c"),
                       limits = c(0, max(d$values)),
                       name = "PCWD\n(mm)") +
  labs(title    = sprintf("Annual maximum PCWD, %s %s", toupper(model), toupper(period)),
       subtitle = "Same colour scale in all panels, 1 km grid (LV95)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid      = element_blank(),
        axis.text       = element_blank(),
        strip.text      = element_text(hjust = 0, face = "bold"),
        plot.subtitle   = element_text(colour = "#52514e"),
        legend.title    = element_text(size = 9, colour = "#52514e"),
        legend.key.height = unit(1.6, "cm"))

print(p)

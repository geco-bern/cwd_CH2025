# run 02_eliminate_NA_CH2025_UBELIX.R first to create tables/qc_na_filled_gridded.csv
# heatmap of the number of replaced values (n_replaced) per model chain and year,
# one figure per variable with one panel per GWL; the plots are shown and saved to figures/
# white = no values replaced, grey = no file for this model chain and GWL

library(here)
library(ggplot2)

# ---- 1. Settings ----------------------------------------------------------
gwls   <- c("ref91-20", "gwl1.5", "gwl2.0", "gwl2.5", "gwl3.0")   # panel order, left to right
labels <- c("ref91-20" = "Reference 1991-2020",
            "gwl1.5"   = "GWL 1.5 °C",
            "gwl2.0"   = "GWL 2.0 °C",
            "gwl2.5"   = "GWL 2.5 °C",
            "gwl3.0"   = "GWL 3.0 °C")
vars   <- c("tas", "tasmin", "tasmax", "pr")                       # panel order, top to bottom

# ---- 2. Read the report ------------------------------------------------------
qc <- read.csv(here("tables", "qc_na_filled_gridded.csv"))

# file name: ogd-climate-scenarios-ch2025-grid_ch_<variable>_<model chain>_<gwl>.nc
qc$chain <- sapply(strsplit(qc$file, "_"), `[`, 4)

# every combination of model chain, variable, GWL and year, so missing files show up as NA
all <- expand.grid(chain = unique(qc$chain), variable = vars, gwl = gwls,
                   year = sort(unique(qc$year)), stringsAsFactors = FALSE)
d <- merge(all, qc[, c("chain", "variable", "gwl", "year", "n_replaced")], all.x = TRUE)

# model chains sorted by driving GCM, then RCM (GCM = last part of the name)
chains  <- unique(d$chain)
gcm     <- sub(".*-", "", chains)
d$chain <- factor(d$chain, levels = rev(chains[order(gcm, chains)]))   # first chain on top
d$variable <- factor(d$variable, levels = vars)
d$gwl      <- factor(d$gwl, levels = gwls, labels = labels[gwls])

# ---- 3. Summary in the console ------------------------------------------------
s <- aggregate(n_replaced ~ chain + variable, data = d[!is.na(d$n_replaced), ],
               FUN = function(x) c(total = sum(x), min_year = min(x), max_year = max(x)))
s <- do.call(data.frame, s)
print(s[s$n_replaced.total > 0, ], row.names = FALSE)

# ---- 4. Plot and save, one figure per variable ---------------------------------
dir.create(here("figures"), showWarnings = FALSE)

for (v in vars) {
  p <- ggplot(d[d$variable == v, ], aes(year, chain, fill = n_replaced)) +
    geom_tile() +
    facet_wrap(~gwl, nrow = 1) +
    scale_x_continuous(breaks = c(1, 10, 20, 30), expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "white", high = "#08306b", na.value = "grey80",
                        labels = scales::label_comma(), name = "Replaced\nvalues") +
    labs(title    = sprintf("Missing values replaced in the CH2025 DAILY-GRIDDED files: %s", v),
         subtitle = "Number of replaced cell-days per year; white = none replaced, grey = no file",
         x = "Year of the 30-year period", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid      = element_blank(),
          panel.border    = element_rect(fill = NA, colour = "grey60", linewidth = 0.3),
          panel.spacing   = unit(0.3, "lines"),
          axis.text.y     = element_text(size = 7),
          axis.text.x     = element_text(size = 7),
          strip.text      = element_text(face = "bold"),
          plot.subtitle   = element_text(colour = "#52514e"),
          legend.title    = element_text(size = 9, colour = "#52514e"),
          legend.key.height = unit(1.2, "cm"))

  print(p)
  ggsave(here("figures", sprintf("qc_na_filled_gridded_%s.jpg", v)),
         p, width = 12, height = 5.5, dpi = 300, bg = "white")
}

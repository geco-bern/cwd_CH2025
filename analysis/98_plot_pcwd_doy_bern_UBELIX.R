# daily PCWD against day of year for the grid cell of Bern-Zollikofen, model chain ensemble:
#   1. per model chain: 30-year mean PCWD per day of year
#   2. per period: ensemble mean over the model chains (line), interquartile range (dark band)
#      and range min-max (light band) of the model chain means
# one figure per period (reference period and the four GWLs), all with the same y axis;
# the GWL figures also show the ensemble mean of the reference period (dashed) for comparison;
# one more figure with all periods stacked (reference on top, then GWL 1.5 to 3.0)
# the 30-year means per model chain are saved to tables/, the figures to figures/

library(here)
library(ncdf4)
library(terra)                                    # needed by R/nc_grid_ch.R
library(ggplot2)
library(parallel)
source(here("R", "nc_grid_ch.R"))                 # grid layout (E_vals, N_vals, n_days)

# ---- 1. Settings ----------------------------------------------------------
site    <- "Bern-Zollikofen"
site_E  <- 2601934                                # CH1903+ / LV95
site_N  <- 1204410
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

periods <- c("ref91-20", "gwl1.5", "gwl2.0", "gwl2.5", "gwl3.0")
labels  <- c("ref91-20" = "Reference 1991-2020",
             "gwl1.5"   = "GWL 1.5 °C",
             "gwl2.0"   = "GWL 2.0 °C",
             "gwl2.5"   = "GWL 2.5 °C",
             "gwl3.0"   = "GWL 3.0 °C")
col_ref <- "#0d366b"                              # reference period
col_gwl <- "#7a2a0c"                              # GWLs

# ---- 2. Inventory of the PCWD files: grid_ch_pcwd_<model chain>_<period>.nc ---------
pcwd_dir <- here("data", "pcwd")
files    <- list.files(pcwd_dir, pattern = "^grid_ch_pcwd_.*\\.nc$", full.names = TRUE)
parts    <- do.call(rbind, strsplit(sub("\\.nc$", "", basename(files)), "_"))
inv      <- data.frame(path = files, chain = parts[, 4], period = parts[, 5])
inv      <- inv[inv$period %in% periods, ]
print(table(inv$period))                          # 26 model chains per period (24 for GWL 3.0)

# ---- 3. 30-year mean PCWD per day of year, per model chain and period ---------------
iE <- which.min(abs(E_vals - site_E))             # E_vals, N_vals are cell centres
iN <- which.min(abs(N_vals - site_N))
n_years <- n_days / 365

read_cell <- function(i) {
  nc   <- nc_open(inv$path[i])
  pcwd <- ncvar_get(nc, "pcwd", start = c(iE, iN, 1), count = c(1, 1, -1))
  nc_close(nc)
  stopifnot(length(pcwd) == n_days, !anyNA(pcwd))
  data.frame(chain  = inv$chain[i], period = inv$period[i], doy = 1:365,
             pcwd   = rowMeans(matrix(pcwd, nrow = 365, ncol = n_years)))   # 365_day calendar
}

res <- mclapply(seq_len(nrow(inv)), read_cell, mc.cores = n_cores)
stopifnot(vapply(res, is.data.frame, NA))         # a failed worker returns an error object
chains <- do.call(rbind, res)

dir.create(here("tables"), showWarnings = FALSE)
write.csv(chains, here("tables", "pcwd_doy_bern_chains.csv"), row.names = FALSE)

# ---- 4. Ensemble statistics over the model chains, per period and day of year -------
ens <- do.call(rbind, lapply(split(chains, list(chains$period, chains$doy), drop = TRUE), function(x)
  data.frame(period = x$period[1], doy = x$doy[1], n_chains = nrow(x),
             mean = mean(x$pcwd), min = min(x$pcwd), max = max(x$pcwd),
             q25  = quantile(x$pcwd, 0.25, names = FALSE),
             q75  = quantile(x$pcwd, 0.75, names = FALSE))))

y_max  <- max(ens$max)                            # same y axis in all figures
ens_ref <- ens[ens$period == "ref91-20", ]

# ---- 5. One figure per period -------------------------------------------------
month_start <- cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))   # doy of the 1st
dir.create(here("figures"), showWarnings = FALSE)

for (per in periods) {
  e   <- ens[ens$period == per, ]
  if (nrow(e) == 0) { message("No files for ", per); next }
  col <- if (per == "ref91-20") col_ref else col_gwl

  p <- ggplot(e, aes(doy)) +
    geom_ribbon(aes(ymin = min, ymax = max, fill = "Range of the model chains"), alpha = 0.25) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = "Interquartile range"), alpha = 0.45) +
    geom_line(aes(y = mean, colour = "Ensemble mean", linetype = "Ensemble mean"), linewidth = 1)

  if (per != "ref91-20") {
    p <- p + geom_line(data = ens_ref, aes(y = mean, colour = "Ensemble mean, reference 1991-2020",
                                           linetype = "Ensemble mean, reference 1991-2020"),
                       linewidth = 0.6)
  }

  p <- p +
    scale_fill_manual(values = c("Range of the model chains" = col, "Interquartile range" = col),
                      breaks = c("Interquartile range", "Range of the model chains"), name = NULL) +
    scale_colour_manual(values = c("Ensemble mean" = col, "Ensemble mean, reference 1991-2020" = col_ref),
                        name = NULL) +
    scale_linetype_manual(values = c("Ensemble mean" = "solid", "Ensemble mean, reference 1991-2020" = "dashed"),
                          name = NULL) +
    guides(fill = guide_legend(override.aes = list(alpha = c(0.45, 0.25)), order = 2),
           colour = guide_legend(order = 1), linetype = guide_legend(order = 1)) +
    scale_x_continuous(breaks = month_start, labels = month.abb,
                       limits = c(1, 365), expand = expansion(add = 2)) +
    scale_y_continuous(limits = c(0, y_max * 1.05), expand = c(0, 0)) +
    labs(title    = sprintf("Daily PCWD at %s, %s", site, labels[per]),
         subtitle = sprintf("30-year mean per day of year of %d model chains; grid cell E %d / N %d (LV95), deficit reset 50 mm below the event peak",
                            e$n_chains[1], E_vals[iE], N_vals[iN]),
         x = NULL, y = "PCWD (mm)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major   = element_line(colour = "grey92"),
          plot.subtitle      = element_text(colour = "#52514e", size = 9),
          legend.position    = "bottom",
          legend.text        = element_text(size = 9))

  ggsave(here("figures", sprintf("pcwd_doy_bern_ensemble_%s.jpg", per)),
         p, width = 10, height = 5.5, dpi = 300, bg = "white")
  message("Saved figure for ", per)
}

# ---- 6. All periods stacked in one figure ----------------------------------------
# panel labels with the number of model chains, panels in the order of `periods`
per_avail <- periods[periods %in% ens$period]
n_ch      <- sapply(per_avail, function(per) ens$n_chains[ens$period == per][1])
strip     <- setNames(sprintf("%s (%d model chains)", labels[per_avail], n_ch), per_avail)

ens$panel <- factor(strip[ens$period], levels = strip)
ens$col   <- ifelse(ens$period == "ref91-20", col_ref, col_gwl)

# reference ensemble mean repeated in every GWL panel
ref_lines <- do.call(rbind, lapply(setdiff(per_avail, "ref91-20"), function(per)
  data.frame(ens_ref[c("doy", "mean")], panel = strip[[per]])))
ref_lines$panel <- factor(ref_lines$panel, levels = strip)

p_all <- ggplot(ens, aes(doy)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = col, alpha = "Range of the model chains")) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = col, alpha = "Interquartile range")) +
  geom_line(aes(y = mean, colour = col, linetype = "Ensemble mean"), linewidth = 0.8) +
  geom_line(data = ref_lines, aes(y = mean, linetype = "Ensemble mean, reference 1991-2020"),
            colour = col_ref, linewidth = 0.5) +
  facet_wrap(~panel, ncol = 1) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_alpha_manual(values = c("Interquartile range" = 0.45, "Range of the model chains" = 0.25),
                     name = NULL) +
  scale_linetype_manual(values = c("Ensemble mean" = "solid", "Ensemble mean, reference 1991-2020" = "dashed"),
                        name = NULL) +
  guides(linetype = guide_legend(override.aes = list(colour = c("grey20", col_ref)), order = 1),
         alpha    = guide_legend(override.aes = list(fill = "grey30"), order = 2)) +
  scale_x_continuous(breaks = month_start, labels = month.abb,
                     limits = c(1, 365), expand = expansion(add = 2)) +
  scale_y_continuous(limits = c(0, y_max * 1.05), expand = c(0, 0)) +
  labs(title    = sprintf("Daily PCWD at %s, model chain ensemble", site),
       subtitle = sprintf("30-year mean per day of year; grid cell E %d / N %d (LV95), deficit reset 50 mm below the event peak",
                          E_vals[iE], N_vals[iN]),
       x = NULL, y = "PCWD (mm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major   = element_line(colour = "grey92"),
        panel.spacing      = unit(0.8, "lines"),
        strip.text         = element_text(hjust = 0, face = "bold"),
        plot.subtitle      = element_text(colour = "#52514e", size = 9),
        legend.position    = "bottom",
        legend.text        = element_text(size = 9))

ggsave(here("figures", "pcwd_doy_bern_ensemble_all.jpg"),
       p_all, width = 9, height = 2.3 * length(per_avail) + 1.5, dpi = 300, bg = "white")
message("Saved stacked figure")

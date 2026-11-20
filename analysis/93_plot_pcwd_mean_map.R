#!/usr/bin/env Rscript
#SBATCH --job-name=pcwd_mean_map
#SBATCH --account=invest
#SBATCH --qos=job_icpu-stocker
#SBATCH --partition=cpu-invest
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=ubelix_logs/pcwd_mean_map_%j.log

# maps of the mean daily PCWD per grid cell, ensemble mean over all model chains
# rows = periods, columns = DJF, MAM, JJA, SON, annual; one shared binned colour scale
# the ensemble mean maps are kept in data/pcwd_mean_maps.rds: delete it to recompute,
# keep it to only redraw the figure (fast, also in the RStudio console)

library(here); library(ncdf4); library(terra); library(ggplot2); library(parallel)
source(here("R", "nc_grid_ch.R"))                 # grid layout (E_vals, N_vals, n_days)

periods <- c("ref91-20" = "Reference\n1991-2020", "gwl1.5" = "GWL 1.5 °C", "gwl2.0" = "GWL 2.0 °C",
             "gwl2.5" = "GWL 2.5 °C", "gwl3.0" = "GWL 3.0 °C")
seasons <- c("DJF", "MAM", "JJA", "SON", "annual")
breaks_all <- c(0, 1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 400, 600)   # colour classes (mm)
file_rds   <- here("data", "pcwd_mean_maps.rds")

# ---- 1. Ensemble mean maps (cells x seasons per period) ---------------------------
if (!file.exists(file_rds)) {
  files <- list.files(here("data", "pcwd"), pattern = "^grid_ch_pcwd_.*\\.nc$", full.names = TRUE)
  per   <- sapply(strsplit(sub("\\.nc$", "", basename(files)), "_"), `[`, 5)
  files <- files[per %in% names(periods)]; per <- per[per %in% names(periods)]
  print(table(per))

  # weights days x seasons: 1 / number of days of the season, so that x %*% w = mean per season
  month  <- findInterval(rep(0:364, n_days / 365), cumsum(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))) + 1
  season <- c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")[month]
  w <- sapply(seasons, function(k) if (k == "annual") rep(1, n_days) else as.numeric(season == k))
  w <- sweep(w, 2, colSums(w), "/")

  # mean per cell and season of one file, read in blocks of 10 rows (= chunks of 07); NA outside the domain
  cell_means <- function(f) {
    nc <- nc_open(f); on.exit(nc_close(nc))
    do.call(rbind, lapply(seq(1, length(N_vals), by = 10), function(s) {
      n <- min(10, length(N_vals) - s + 1)
      matrix(ncvar_get(nc, "pcwd", start = c(1, s, 1), count = c(-1, n, -1), collapse_degen = FALSE),
             ncol = n_days) %*% w
    }))
  }
  res <- mclapply(files, cell_means, mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
                  mc.preschedule = FALSE)
  bad <- !sapply(res, is.matrix)                    # error or killed worker (e.g. out of memory)
  if (any(bad)) stop(sum(bad), " files failed, resubmit:\n", paste(files[bad], collapse = "\n"))

  ens <- lapply(split(res, per), function(r) Reduce(`+`, r) / length(r))   # mean over the model chains
  saveRDS(ens, file_rds)
}
ens <- readRDS(file_rds)

# ---- 2. Long table, colour classes, panel letters -------------------------------
xy <- expand.grid(x = E_vals, y = N_vals)           # cell order of the files (E fastest, N south to north)
d  <- do.call(rbind, lapply(intersect(names(periods), names(ens)), function(p)
  do.call(rbind, lapply(seasons, function(k) data.frame(xy, period = p, season = k, pcwd = ens[[p]][, k])))))
d  <- d[!is.na(d$pcwd), ]
d$period <- factor(periods[d$period], levels = periods[names(periods) %in% names(ens)])
d$season <- factor(d$season, levels = seasons)

v_max <- max(d$pcwd)
brks  <- c(breaks_all[breaks_all < v_max],
           if (any(breaks_all >= v_max)) min(breaks_all[breaks_all >= v_max]) else ceiling(v_max))
d$bin <- cut(d$pcwd, brks, include.lowest = TRUE)

pan <- expand.grid(season = levels(d$season), period = levels(d$period))   # (a), (b), ... row by row
pan$letter <- sprintf("(%s)", letters[seq_len(nrow(pan))])

x_lim <- range(d$x) + c(-500, 500)                  # extent of the cells with data
y_lim <- range(d$y) + c(-500, 500)

# ---- 3. Plot ---------------------------------------------------------------
p <- ggplot(d, aes(x, y, fill = bin)) +
  geom_raster(show.legend = TRUE) +                  # TRUE: also show classes without cells
  geom_text(data = pan, aes(x_lim[1], y_lim[2], label = letter), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.5) +
  facet_grid(period ~ season, switch = "y") +
  coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  scale_fill_manual(values = setNames(hcl.colors(nlevels(d$bin), "YlOrRd", rev = TRUE), levels(d$bin)),
                    labels = paste(head(brks, -1), "-", brks[-1]), drop = FALSE, name = "PCWD (mm)") +
  guides(fill = guide_legend(reverse = TRUE, keyheight = unit(0.7, "cm"))) +
  theme_void(base_size = 12) +
  theme(strip.text.x      = element_text(size = 13, margin = margin(b = 4)),
        strip.text.y.left = element_text(size = 12, angle = 90, margin = margin(r = 4)),
        strip.placement   = "outside",
        panel.spacing     = unit(0.4, "lines"),
        plot.background   = element_rect(fill = "white", colour = NA))

dir.create(here("figures"), showWarnings = FALSE)
ggsave(here("figures", "pcwd_mean_map_ch.jpg"), p, width = 16, height = 10, dpi = 300, bg = "white")
message("Saved figure")

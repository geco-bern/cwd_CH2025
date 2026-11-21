# Shared binned colour scale for PCWD maps (mm), so that the same value has the same colour in every figure
# (93_plot_pcwd_mean_map.R, 91_plot_pcwd_annmax_map.R)
#   d$bin <- pcwd_bin(d$pcwd)                     # colour class of every cell
#   ... + scale_fill_pcwd(d$bin, "PCWD (mm)")
# the legend shows the classes from 0 up to the class of the largest value
# needs: library(ggplot2)

pcwd_breaks <- c(0, 1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 400, 600, Inf)
pcwd_cols   <- c(hcl.colors(14, "YlOrRd", rev = TRUE),   # 0 - 600 mm, light yellow to dark red
                 "#3b0a14")                              # > 600 mm, almost black
pcwd_labels <- c(paste(head(pcwd_breaks, -2), "-", pcwd_breaks[2:(length(pcwd_breaks) - 1)]),
                 paste(">", pcwd_breaks[length(pcwd_breaks) - 1]))

# colour class of every value; the levels end with the class of the largest value
pcwd_bin <- function(v) {
  b <- cut(v, pcwd_breaks, include.lowest = TRUE)
  factor(b, levels = levels(b)[seq_len(max(as.integer(b), na.rm = TRUE))])
}

# fill scale for the output of pcwd_bin(); classes without cells are shown in the legend as well
scale_fill_pcwd <- function(bin, name = "PCWD (mm)") {
  k <- nlevels(bin)
  scale_fill_manual(values = setNames(pcwd_cols[seq_len(k)], levels(bin)),
                    labels = pcwd_labels[seq_len(k)], drop = FALSE, name = name)
}

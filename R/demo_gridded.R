if(!require(terra)) install.packages("terra")
if(!require(ncdf4)) install.packages("ncdf4")
library(terra)
library(ncdf4)

r <- rast("C:/Users/jucke/OneDrive - Universitaet Bern/BSc Geographie/05 HS 26/BT/CH2025/Gridded/Sample/ogd-climate-scenarios-ch2025-grid_ch_pr_clmcom-cclm4-cccma_gwl1.5.nc")
r                 # zeigt Anzahl Layer, Auflösung, Extent, CRS
names(r)          # Namen der einzelnen Layer (z. B. Zeitpunkte oder Perzentile)

plot(r[[1]], main = "Niederschlag, GWL 1.5°C (CLMcom-CCLM4 / CanESM2)")

plot(r)  # alle Layer als Panel-Plot


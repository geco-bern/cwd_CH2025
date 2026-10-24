# Drought projections for Switzerland (Work in Progress)

Evaluating the CH2025 climate projections for changes in drought across Switzerland, quantified using the potential cumulative water deficit (PCWD). PCWD is the cumulative difference between potential evapotranspiration (PET) and the liquid water reaching the soil (Helpap et al., 2026). PET is calculated from net radiation following Priestley & Taylor (1972) as implemented by Davis et al. (2017). Because CH2025 provides no radiation on its 1 km grid, net radiation is modelled from the gridded temperature and precipitation. The thesis will compare drought frequency and magnitude across the global warming levels (GWL) 1.5, 2.0, 2.5 and 3.0 °C and aggregate the results to the biogeographic regions of Switzerland.

This repository contains the code for a BSc thesis in Geography at the University of Bern. The thesis proposal is in `vignettes/Proposal.Rmd`.

## Project overview

Contains data and code for creating:

- Gap-filled CH2025 DAILY-GRIDDED climate files and an elevation grid that matches the CH2025 1 km grid. Outputs are stored in `data-raw/` (filled in place) and `data/elevation/`.
- Daily clear-sky shortwave radiation (Rsw-cc) for the CH2025 stations and the 1 km grid. Outputs are stored in `data/rsw-cc/`.
- A two-step random forest model for daily global radiation (Rg) and net radiation (Rn), with cross-validation tables and figures. Outputs are stored in `data/radiation_models/`, `tables/` and `figures/`.
- Daily simulated net radiation on the CH2025 grid. Outputs are stored in `data/netradiation_sim/`.
- Daily PCWD and annual maximum PCWD for every grid cell. Outputs are stored in `data/pcwd/` and `data/pcwd_annmax/`.

The full workflow is reproduced by running the scripts in this order:

1. `analysis/01_get_data_CH2025.R` downloads the CH2025 DAILY-LOCAL and DAILY-GRIDDED files from the MeteoSwiss STAC API.
2. `analysis/02_eliminate_NA_CH2025.R` fills missing values in the gridded files with the monthly mean of the same cell.
3. `R/elevation.R` resamples DHM25 to the CH2025 1 km grid.
4. `R/calc_rsw-cc.R` calculates Rsw-cc for the stations and the grid with meteoland.
5. `analysis/03_model_Rn.R` trains the Rg model on the CH2025 stations and the Rn model on four SwissMetNet stations, both evaluated with leave-one-site-out cross-validation.
6. `analysis/04_simulate_Rn.R` applies both models to the gridded data.
7. `analysis/05_apply_pcwd_CH2025.R` calculates daily PCWD on the grid with the cwd R package.
8. `analysis/06_get_pcwd_annmax_CH2025.R` takes the annual maximum PCWD of every grid cell.

The model chain and GWL are set in the settings section at the top of scripts 04 to 06. The gridded outputs are NetCDF files on the CH2025 grid in CH1903+ / LV95, written with the helpers in `R/nc_grid_ch.R`.

## Data

The following data were obtained from external sources and are stored in `data-raw/`. Because of their size they are not included in this repository. The CH2025 files can be downloaded with `analysis/01_get_data_CH2025.R`.

### Climate CH2025

- DAILY-GRIDDED: daily mean, minimum and maximum temperature (tas, tasmin, tasmax) and precipitation (pr) on a 1 km grid. There is one NetCDF file per variable, model chain and GWL, each covering 30 years in a 365-day calendar (MeteoSwiss, 2025).
- DAILY-LOCAL: the same variables plus daily global radiation (rsds) at stations. Rsds is used to train the Rg model (MeteoSwiss, 2025).

### SwissMetNet radiation

- Daily measurements of incoming and outgoing shortwave and longwave radiation at Altdorf, Magadino, Payerne and Aadorf, the four stations that measure all four components. Observed Rn is derived from these and used to train the Rn model (MeteoSwiss, n.d.).

### Elevation

- DHM25, the digital height model of swisstopo with a grid size of 200 m. It is resampled to the 1 km grid to calculate atmospheric pressure (swisstopo, n.d.).

## Project structure

```         
├── README.md                <- The top-level README includes instructions to use this repository
|
├── cwd_CH2025.Rproj         <- R project file
| 
├── .gitignore               <- file indicating which files should be ignored when pushing
|
├── analysis/                <- R scripts of the workflow, numbered in the order they are run
|
├── data-raw/                <- folder for data downloaded from the external sources (not on GitHub)
|
├── data/                    <- folder for data produced by the repository (only data/elevation/
|                               and data/rsw-cc/ are on GitHub)
│
├── figures/                 <- folder for figure files produced by the repository
│
├── tables/                  <- folder for tables produced by the repository (gap filling,
|                               cross-validation of the radiation models)
│
├── vignettes/               <- R markdown files
│   ├── Proposal.Rmd         <- Thesis proposal: context, objective, methods, workflow and timeline
│   └── references.bib       <- bibliography file for .Rmd compilation
|
└── R/                       <- R functions and helper scripts used in the project
```

## Dependencies

Install all required R libraries by:

``` r
use_pkgs <- c(
  "here",
  "jsonlite",
  "digest",
  "ncdf4",
  "terra",
  "meteoland",
  "dplyr",
  "ggplot2",
  "ranger",
  "readr",
  "stringr",
  "leaflet",
  "htmltools",
  "htmlwidgets",
  "tidyverse",
  "rmarkdown",
  "remotes"
  )

new_pkgs <- use_pkgs[!(use_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

# cwd R package for snow, PET and PCWD (Stocker, 2026)
remotes::install_github("geco-bern/cwd")
```

## License

This is published under a [CC BY-SA license](https://creativecommons.org/licenses/by-sa/4.0/).

## References

Davis, T. W., Prentice, I. C., Stocker, B. D., Thomas, R. T., Whitley, R. J., Wang, H., Evans, B. J., Gallego-Sala, A. V., Sykes, M. T., and Cramer, W.: Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of radiation, evapotranspiration and plant-available moisture, Geoscientific Model Development, 10, 689–708, <https://doi.org/10.5194/gmd-10-689-2017>, 2017.

De Cáceres, M., Martin-StPaul, N., Turco, M., Cabon, A., and Granda, V.: Estimating daily meteorological data and downscaling climate models over landscapes, Environmental Modelling & Software, 108, 186–196, <https://doi.org/10.1016/j.envsoft.2018.08.003>, 2018.

Helpap, P., Brönnimann, S., Hand, R., Franke, J., Raible, C. C., and Stocker, B. D.: Are modern droughts unprecedented?, AGU Advances, 7, e2026AV002536, <https://doi.org/10.1029/2026AV002536>, 2026.

MeteoSwiss: Climate CH2025 datasets, <https://www.meteoswiss.admin.ch/climate/climate-change/swiss-climate-scenarios/reports-data-and-graphs-of-climate-change-scenarios/climate-ch2025-datasets.html>, last access: 24 September 2026, 2025.

MeteoSwiss: Strahlungsbeobachtungen, <https://www.meteoschweiz.admin.ch/klima/klima-der-schweiz/strahlungsbeobachtungen.html>, last access: 10 October 2026, n.d.

MeteoSwiss and ETH Zurich: Climate CH2025 – Scientific Report, Federal Office of Meteorology and Climatology MeteoSwiss, Zurich, <https://doi.org/10.18751/climate/scenarios/ch2025/sr/1.0/>, 2025.

Priestley, C. H. B. and Taylor, R. J.: On the assessment of surface heat flux and evaporation using large-scale parameters, Monthly Weather Review, 100, 81–92, <https://doi.org/10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2>, 1972.

Stocker, B.: cwd v3.0.2: R package for cumulative water deficit calculation, Zenodo, <https://doi.org/10.5281/zenodo.21679200>, 2026.

swisstopo: DHM25, <https://www.swisstopo.admin.ch/en/height-model-dhm25>, last access: 3 October 2026, n.d.

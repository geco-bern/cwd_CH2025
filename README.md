# Drought projections for Switzerland

Evaluating the CH2025 climate projections for trends in drought, quantified using the potential cumulative water deficit (PCWD). PCWD can be calculated using net radiation and the potential evapotranspiration (PET). Here, we calculate PET based on net radiation and the approach by Priestley & Taylor (1972) as implemented by Davis et al. (2017).

use local dataset and not gridded. local includes radiation

## Project overview

Contains data and code for creating:

- Homogenised time series of atmospheric CO<sub>2</sub> and temperature from multiple sources, together covering the past (800 ka), until today, and the future (climate scenario to 2100). Outputs (homogenised time series) are stored in `data/`
- Visualisation of the parallel evolution of temperature and CO<sub>2</sub>. Outputs (figure files) are stored in `figures/`.

The full workflow can be reproduced by running `vignettes/past_to_future_CO2_temperature.Rmd`.

## Data

The following data were obtained from external sources and are contained in this repository, stored in `data-raw/`.

### CO<sub>2</sub>

- Bereiter et al., 2015
- Meinshausen et al., 2011

### Temperature (based on d<sup>18</sup>O)

- Jouzel et al., 2007
- Neukom et al., 2009

## Project structure

```         
├── README.md                <- The top-level README includes instructions to use this repository
|                               and the project proposal for the Proseminar
│
├── inaug_demo.Rproj.        <- R project file
| 
├── .gitignore               <- file indicating which files should be ignored when pushing
|
├── data-raw/                <- folder for data downloaded from the external sources, unprocessed
|
├── data/                    <- folder for data produced by the repository
│
├── figures/                 <- folder for figure files produced by the repository 
│
├── vignettes/               <- R markdown files
│   ├── past_to_future_co2_temperature.Rmd  <- Contains the full workflow of data read, processing, and visualisation.
│   └── references.bib       <- bibliography file for .Rmd compilation
|
└── R/                       <- R functions used in the project, contains one function per script
```

## Dependencies

Install all required R libraries by:

``` r
use_pkgs <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "lubridate",
  "readr",
  "ggplot2",
  "gganimate",
  "gifski",
  "here",
  "readr"
  )

new_pkgs <- use_pkgs[!(use_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
```

## License

This is published under a [CC BY-SA license](https://creativecommons.org/licenses/by-sa/4.0/).

## References

Bereiter, B., Eggleston, S., Schmitt, J., Nehrbass‐Ahles, C., Stocker, T. F., Fischer, H., Kipfstuhl, S., and Chappellaz, J.: Revision of the EPICA Dome C CO<sub>2</sub> record from 800 to 600 kyr before present, Geophysical Research Letters, 42, 542–549, <https://doi.org/10.1002/2014GL061957>, 2015.

Jouzel, J., Masson-Delmotte, V., Cattani, O., Dreyfus, G., Falourd, S., Hoffmann, G., Minster, B., Nouet, J., Barnola, J. M., Chappellaz, J., Fischer, H., Gallet, J. C., Johnsen, S., Leuenberger, M., Loulergue, L., Luethi, D., Oerter, H., Parrenin, F., Raisbeck, G., Raynaud, D., Schilt, A., Schwander, J., Selmo, E., Souchez, R., Spahni, R., Stauffer, B., Steffensen, J. P., Stenni, B., Stocker, T. F., Tison, J. L., Werner, M., and Wolff, E. W.: Orbital and Millennial Antarctic Climate Variability over the Past 800,000 Years, Science, 317, 793–796, <https://doi.org/10.1126/science.1141038>, 2007.

Meinshausen, M., S. J. Smith, K. V. Calvin, J. S. Daniel, M. L. T. Kainuma, J.-F. Lamarque, K. Matsumoto, S. A. Montzka, S. C. B. Raper, K. Riahi, A. M. Thomson, G. J. M. Velders and D. van Vuuren (2011). "The RCP Greenhouse Gas Concentrations and their Extension from 1765 to 2300." Climatic Change (Special Issue), DOI: 10.1007/s10584-011-0156-z, freely available online (PDF) (HTML)

Neukom, R., Barboza, L. A., Erb, M. P., Shi, F., Emile-Geay, J., Evans, M. N., Franke, J., Kaufman, D. S., Lücke, L., Rehfeld, K., Schurer, A., Zhu, F., Brönnimann, S., Hakim, G. J., Henley, B. J., Ljungqvist, F. C., McKay, N., Valler, V., von Gunten, L., and PAGES 2k Consortium: Consistent multidecadal variability in global temperature reconstructions and simulations over the Common Era, Nat. Geosci., 12, 643–649, <https://doi.org/10.1038/s41561-019-0400-0>, 2019.

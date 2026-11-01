# CH2025 stations: which variables are available at which station?
#
# read the station metadata
# extend the per-variable station lists to every station
# (one TRUE/FALSE column per variable)
# draw a dot matrix of all stations (figure for the thesis appendix)
# draw a map of all stations, highlighting the Rg model stations
# draw the elevation of the Rg model stations

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(here)

# Paths
meta_file <- here("data-raw", "CH2025", "DAILY_LOCAL", "ogd-climate-scenarios-ch2025_meta_stations.csv")
out_jpg   <- here("figures", "station_variable_availability.jpg")
dem_file  <- here("data", "elevation", "elevation_1k_grid.tif")
dir_gadm  <- here("data-raw", "gadm")                   # national border, downloaded once
out_map   <- here("figures", "station_map.jpg")
out_elev  <- here("figures", "station_elevation_rg.jpg")

# Read station metadata
# semicolon-separated, UTF-8 with BOM (readr removes the BOM)
stations <- read_delim(meta_file, delim = ";",
                       locale = locale(encoding = "UTF-8"),
                       show_col_types = FALSE)

# some codes carry trailing blanks in the file (e.g. "NABZUE ")
stations <- stations |>
  mutate(station_abbr = str_trim(station_abbr))

stopifnot(!anyDuplicated(stations$station_abbr))

# Variable lists (station codes per variable)
codes <- function(x) str_split(str_squish(x), " ")[[1]]

var_lists <- list(

  pr = codes("ABE ABG ABO AFI AFT AIE AIG AIR ALM ALS ALT ALV ALW AMW AND ANT APP APT ARB ARI
    ARO AST AUB AVA AVB AVE BAM BAS BAT BAU BAW BEP BER BEX BEY BEZ BIA BIE BIN BIO
    BIS BIV BIZ BLA BLS BLU BLZ BMU BNE BNU BOL BOS BOV BOZ BRA BRI BRL BRO BRP BRT
    BRW BRZ BSG BSP BUC BUD BUE BUF BUS CAC CAV CDF CDM CEV CGI CHA CHB CHD CHM
    CHU CHY CHZ CIM CMD CNZ COG COL COM COP COS COU COV COY CTA CTO CUE DAV DEH
    DEM DIB DIE DIS DIT DMA DOB DOL EBK ECH EFF EGL EGO EGR EHM EIN EIT EKO ELM ENG
    ENT EPT ERB ERI ESZ EUT EVO FAH FAI FIL FIO FIT FLI FLU FLW FRA FRC FRE FRF FRI FRU
    GAD GEA GEP GHS GIH GIN GLA GOA GOE GOS GRA GRC GRH GRI GRO GRY GSB GSG GSS
    GTT GUE GUG GUT GVE GWA HAI HAU HEB HEK HER HES HIW HLL HON HOY HTW HUT ILZ
    INF INN INT IST JAU JON JUS KAI KAR KAS KIB KIE KIS KLA KLO KLT KOP KRO KSE KUA KUB
    KUE LAB LAC LAF LAG LAN LAP LAT LEH LEI LEU LFB LGA LOB LOC LOH LON LOT LSN LTB
    LUG LUN LUT LUZ MAB MAC MAD MAG MAL MAR MAS MAT MBA MDO MER MES MEV MGB MGI
    MGL MLS MMO MOA MOB MOD MOE MOS MOU MSG MST MTE MTO MUB MUE MUL MUR MUS
    MUT MVE NAP NEU NIE OBD OBI OED OLI OPF ORS ORZ OTE OTL PAV PAY PDM PFA PIG PIL
    PIO PLF PON POT PSI PUD PUY RAG REG REH REM RIC RIH RIX ROB ROE ROM ROT RUE SAB
    SAE SAF SAG SAI SAM SAN SAP SAR SAX SBE SBO SCA SCD SCE SCH SCU SDO SED SEM
    SEV SGD SHA SIA SIE SIH SIM SIO SIS SKO SMA SMM SNS SOG SON SPA SPZ SRE SRN STB
    STE STG STP STU SUA SUS SVG SWA SWZ SZB TAE TAF TAM TAV TAW TDG TEU TFD THS
    THU TIC TRO TRU TST UBB ULR UNK UNS URB URN UST VAD VAE VAR VEL VEV VIG VIO VIS
    VIT VLS VRI VST VVI WAE WAG WAN WBR WCH WDO WEE WEF WET WFJ WHA WHF WID WIE
    WIN WIS WIT WIW WYN YVN ZER ZEV ZNZ ZOF ZUB ZWE ZWK"),
    
  tas = codes("ABO AIG ALT ANT ARO BAS BER BLA BRL BUF BUS CDF CGI CHA CHD CHM CHU CIM COM COV
    DAV DEM DIS DOL EBK EIN ELM ENG EVO FAH FRE GLA GRA GRC GRH GRO GSB GUE GUT
    GVE HAI HLL HOE INT JUN KLO KOP LAE LAG LUG LUZ MAG MER MLS MVE NAP NEU OTL PAY
    PIL PIO PLF PUY RAG REH ROB ROE RUE SAE SAM SBE SBO SCU SHA SIA SIO SMA SMM STG
    TAE ULR VAD VIS WAE WFJ WYN ZER NABZUE NABLAU"),
  
  tasmax = codes("ABO AIG ALT ARO BAS BER BLA BRL BUF BUS CDF CGI CHA CHD CHM CHU CIM COM COV DAV
    DEM DIS DOL EBK EIN ELM ENG EVO FAH FRE GLA GRA GRC GRH GRO GSB GUE GUT GVE
    HAI HLL HOE INT JUN KLO KOP LAE LAG LUG LUZ MAG MER MLS MVE NAP NEU OTL PAY PIL
    PIO PLF PUY RAG REH ROB ROE RUE SAE SAM SBE SBO SCU SHA SIA SIO SMA SMM STG TAE
    ULR VAD VIS WAE WFJ WYN ZER NABZUE NABLAU"),
  
  tasmin = codes("ABO AIG ALT ANT ARO BAS BER BLA BRL BUF BUS CDF CGI CHA CHD CHM CHU CIM COM COV
    DAV DEM DIS DOL EBK EIN ELM ENG EVO FAH FRE GLA GRA GRC GRH GRO GSB GUE GUT
    GVE HAI HLL HOE INT JUN KLO KOP LAE LAG LUG LUZ MAG MER MLS MVE NAP NEU OTL PAY
    PIL PIO PLF PUY RAG REH ROB ROE RUE SAE SAM SBE SBO SCU SHA SIA SIO SMA SMM STG
    TAE ULR VAD VIS WAE WFJ WYN ZER NABZUE NABLAU"),
  
  rsds = codes("ABO AIG ALT BAS BER BUS CDF CGI CHA CHU CIM COM COV DAV DIS DOL ENG EVO FAH FRE
    GLA GRH GSB GUE GUT GVE HOE INT JUN KLO LAE LUG LUZ MAG MLS MVE NAP NEU OTL
    PAY PIL PIO PLF PUY REH ROB ROE RUE SAE SAM SBE SBO SCU SHA SIO SMA STG TAE ULR
    VAD VIS WAE WFJ WYN ZER NABZUE NABLAU"),
  
  hurs = codes("ABO AIG ALT ANT ARO BAS BER BRL BUF BUS CDF CGI CHA CHD CHM CHU CIM COM COV
    DAV DEM DIS DOL EBK EIN ELM ENG EVO FAH FRE GLA GRA GRC GRH GRO GSB GUE GUT
    GVE HAI HLL HOE INT JUN KLO KOP LAE LAG LUG LUZ MAG MER MLS MVE NAP NEU OTL PAY
    PIL PIO PLF PUY RAG REH ROB ROE RUE SAE SAM SBE SBO SCU SHA SIA SIO SMA SMM STG
    TAE ULR VAD VIS WAE WFJ WYN ZER NABZUE NABLAU"),
  
  sfcWind = codes("ABO AIG ALT ARO BAS BER BRL BUF BUS CDF CGI CHA CHD CHM CHU CIM COM COV DAV
    DEM DIS DOL EBK EIN ELM ENG EVO FAH FRE GLA GOE GRA GRH GSB GUE GUT GVE HAI HLL
    HOE INT JUN KLO KOP LAE LAG LEI LUG LUZ MAG MER MLS MUB MVE NAP NEU OTL PAY PIL
    PIO PLF PUY RAG REH ROB ROE RUE SAE SAM SBE SBO SCU SHA SIO SMA SMM STG TAE
    ULR VAD VIS WAE WFJ WYN ZER NABZUE NABLAU")
)

var_labels <- c(
  pr      = "precipitation",
  tas     = "surface temperature",
  tasmax  = "surface maximum temperature",
  tasmin  = "surface minimum temperature",
  rsds    = "radiation",
  hurs    = "relative humidity",
  sfcWind = "surface wind"
)
stopifnot(identical(names(var_lists), names(var_labels)))

# Extend the lists to all stations
# one logical column per variable: TRUE = station has this variable
for (v in names(var_lists)) {
  stations[[v]] <- stations$station_abbr %in% var_lists[[v]]
}

# codes in the lists that do not exist in the metadata (should be none)
unknown <- setdiff(unlist(var_lists), stations$station_abbr)
if (length(unknown) > 0) {
  warning("Codes not found in station metadata: ", paste(unknown, collapse = ", "))
}

# Dot matrix
n_panels <- 5                                   # columns of stations on the page
vars     <- names(var_labels)
n_rows   <- ceiling(nrow(stations) / n_panels)

# stations used for the Rg model (shown in bold)
rg_stations <- toupper(c(

  "abo", "aig", "alt", "bas", "ber", "bus", "cdf", "cgi", "cha", "chu",
  "cim", "com", "cov", "dav", "dis", "dol", "eng", "evo", "fah", "fre",
  "gla", "grh", "gsb", "gue", "gut", "gve", "int", "klo", "lug", "luz",
  "mag", "mls", "mve", "nap", "neu", "otl", "pay", "pil", "pio", "plf",
  "puy", "reh", "rob", "roe", "rue", "sae", "sam", "sbe", "sbo", "scu", "sha",
  "sio", "sma", "stg", "tae", "ulr", "vad", "vis", "wae", "wfj", "wyn",
  "zer"
))
rg_rows <- stations$station_abbr %in% rg_stations
stopifnot(sum(rg_rows) == length(rg_stations),
          all(stations[rg_rows, c("rsds", "tasmax", "tasmin", "pr")]))

bold_codes <- function(x) {
  parse(text = ifelse(x %in% rg_stations, sprintf('bold("%s")', x), sprintf('"%s"', x)))
}

dots <- stations |>
  arrange(station_abbr) |>
  mutate(row   = row_number(),
         panel = ceiling(row / n_rows),
         code  = factor(station_abbr, levels = rev(station_abbr))) |>
  select(row, panel, code, all_of(vars)) |>
  pivot_longer(all_of(vars), names_to = "variable", values_to = "available") |>
  mutate(x     = match(variable, vars),
         state = ifelse(available, "Available", "Not available"),
         band  = (row %% n_rows) %% 2 == 0)    # light band on every other row of a panel

var_names <- paste0(vars, " (", colSums(as.matrix(stations[vars])), ")")
caption   <- paste0("Stations in bold contain the daily global radiation (rsds), temperature maximum ",
                    "(tasmax), temperature minimum (tasmin) and precipitation (pr) necessary for ",
                    "the Rg model. Number of stations per variable in brackets. ",
                    paste0(vars, " = ", var_labels, collapse = ", "), ".")

p <- ggplot(dots, aes(x, code)) +
  geom_tile(data = filter(dots, x == 1), aes(x = 4, fill = band), width = 7, height = 1) +
  geom_point(aes(colour = state, size = state)) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq_along(vars), labels = var_names, position = "top",
                     expand = expansion(add = 0.5)) +
  scale_y_discrete(labels = bold_codes) +
  scale_fill_manual(values = c("TRUE" = "#EEF0F3", "FALSE" = "white"), guide = "none") +
  scale_colour_manual(values = c("Available" = "#1F2A37", "Not available" = "#B8BEC7")) +
  scale_size_manual(values = c("Available" = 1.5, "Not available" = 0.35)) +
  labs(x = NULL, y = NULL, colour = NULL, size = NULL,
       caption = str_wrap(caption, 118)) +
  guides(colour = guide_legend(override.aes = list(size = c(2.4, 0.8)))) +
  theme_minimal(base_size = 8) +
  theme(panel.grid       = element_blank(),
        strip.text       = element_blank(),
        panel.spacing.x  = unit(5, "mm"),
        axis.text        = element_text(colour = "#1F2A37"),
        axis.text.x.top  = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6.5),
        axis.text.y      = element_text(size = 6),
        legend.position  = "bottom",
        legend.text      = element_text(size = 8),
        plot.caption     = element_text(hjust = 0, size = 7, colour = "#4B5563"),
        plot.caption.position = "plot",
        plot.margin      = margin(4, 4, 4, 4))

# Save
dir.create(dirname(out_jpg), showWarnings = FALSE, recursive = TRUE)
ggsave(out_jpg, p, width = 17, height = 24, units = "cm", dpi = 300, quality = 95, bg = "white")
message("Figure saved to: ", normalizePath(out_jpg))

# Map of all stations
# base map: national border and light relief from the 1 km elevation grid
dir.create(dir_gadm, showWarnings = FALSE, recursive = TRUE)
border <- terra::project(geodata::gadm("CHE", level = 0, path = dir_gadm), "EPSG:2056")

dem  <- terra::rast(dem_file)
hill <- terra::shade(terra::terrain(dem, "slope",  unit = "radians"),
                     terra::terrain(dem, "aspect", unit = "radians"), angle = 45, direction = 315)
hill <- as.data.frame(terra::mask(hill, border), xy = TRUE)
names(hill)[3] <- "shade"
border_lines <- terra::geom(terra::as.lines(border), df = TRUE)

group_rg <- sprintf("Rg model stations (%d)", sum(rg_rows))
group_ot <- sprintf("Other stations (%d)", sum(!rg_rows))

map_points <- stations |>
  mutate(east  = station_coordinates_lv95_east,
         north = station_coordinates_lv95_north,
         group = ifelse(rg_rows, group_rg, group_ot)) |>
  arrange(rg_rows)                               # Rg model stations on top

m <- ggplot() +
  geom_raster(data = hill, aes(x, y, fill = shade)) +
  geom_path(data = border_lines, aes(x, y, group = interaction(geom, part)),
            colour = "#6B7280", linewidth = 0.3) +
  geom_point(data = filter(map_points, group == group_rg), aes(east, north),
             colour = "white", size = 2.9) +                # white ring around the blue dots
  geom_point(data = map_points, aes(east, north, colour = group, size = group)) +
  ggrepel::geom_text_repel(data = filter(map_points, group == group_rg),
                           aes(east, north, label = station_abbr), size = 2, fontface = "bold",
                           colour = "#1F2A37", bg.color = "white", bg.r = 0.15,
                           point.size = 2.9, box.padding = 0.1, point.padding = 0,
                           min.segment.length = 0.3, segment.size = 0.2,
                           segment.colour = "#6B7280", max.overlaps = Inf, seed = 1) +
  annotate("segment", x = 2490000, xend = 2540000, y = 1078000, yend = 1078000,
           colour = "#1F2A37", linewidth = 0.6) +
  annotate("text", x = 2515000, y = 1083000, label = "50 km", size = 2.6, vjust = 0,
           colour = "#1F2A37") +
  scale_fill_gradient(low = "#D3D7DD", high = "#F8F9FA", guide = "none") +
  scale_colour_manual(values = setNames(c("#1F6FB2", "#8A9099"), c(group_rg, group_ot)),
                      breaks = c(group_rg, group_ot)) +
  scale_size_manual(values = setNames(c(2.1, 0.6), c(group_rg, group_ot)),
                    breaks = c(group_rg, group_ot)) +
  labs(colour = NULL, size = NULL) +
  coord_equal(expand = FALSE, xlim = c(2474000, 2844000), ylim = c(1064000, 1304000)) +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8, colour = "#1F2A37"))

ggsave(out_map, m, width = 17, height = 12, units = "cm", dpi = 300, quality = 95, bg = "white")
message("Map saved to: ", normalizePath(out_map))

# Elevation of the Rg model stations, highest on top
elev <- stations[rg_rows, ] |>
  arrange(station_height_masl) |>
  mutate(label = paste0(station_name, " (", station_abbr, ")"),
         label = factor(label, levels = label))
stopifnot(!anyNA(elev$station_height_masl))

e <- ggplot(elev, aes(station_height_masl, label)) +
  geom_segment(aes(x = 0, xend = station_height_masl, yend = label),
               colour = "#B8BEC7", linewidth = 0.3) +
  geom_point(colour = "#1F6FB2", size = 1.7) +
  scale_x_continuous(breaks = seq(0, 4000, by = 500), limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.02)), position = "top") +
  labs(x = "Elevation (m a.s.l.)", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "#E5E7EB", linewidth = 0.3),
        axis.text          = element_text(colour = "#1F2A37"),
        axis.text.y        = element_text(size = 6.5),
        plot.margin        = margin(4, 8, 4, 4))

ggsave(out_elev, e, width = 17, height = 20, units = "cm", dpi = 300, quality = 95, bg = "white")
message("Elevation plot saved to: ", normalizePath(out_elev))

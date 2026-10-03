# ------------------------------------------------------------------
# CH2025 stations: which variables are available at which station?
#
# 1. read the station metadata
# 2. extend the per-variable station lists to every station
#    (one TRUE/FALSE column per variable)
# 3. draw an interactive map (saved as HTML when run on its own;
#    the map object `m` is used directly when sourced from an .Rmd)
#
# Packages: install.packages(c("readr", "dplyr", "stringr", "leaflet",
#                              "htmltools", "htmlwidgets"))
# ------------------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(leaflet)
library(htmltools)
library(htmlwidgets)
library(here)

# TRUE when the script is sourced from an R Markdown document being knitted
knitting <- isTRUE(getOption("knitr.in.progress"))

# ---- 1. Paths -----------------------------------------------------
meta_file <- here("data-raw", "CH2025", "DAILY_LOCAL", "ogd-climate-scenarios-ch2025_meta_stations.csv")
out_html  <- here("figures", "station_variable_availability_map.html")

# ---- 2. Read station metadata ------------------------------------
# semicolon-separated, UTF-8 with BOM (readr removes the BOM)
stations <- read_delim(meta_file, delim = ";",
                       locale = locale(encoding = "UTF-8"),
                       show_col_types = FALSE)

# some codes carry trailing blanks in the file (e.g. "NABZUE ")
stations <- stations |>
  mutate(station_abbr = str_trim(station_abbr))

stopifnot(!anyDuplicated(stations$station_abbr))

# ---- 3. Variable lists (station codes per variable) --------------
codes <- function(x) str_split(str_squish(x), " ")[[1]]

var_lists <- list(
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
  tas     = "surface temperature",
  tasmax  = "surface maximum temperature",
  tasmin  = "surface minimum temperature",
  pr      = "precipitation",
  rsds    = "radiation",
  hurs    = "relative humidity",
  sfcWind = "surface wind"
)
stopifnot(identical(names(var_lists), names(var_labels)))

# ---- 4. Extend the lists to all stations -------------------------
# one logical column per variable: TRUE = station has this variable
for (v in names(var_lists)) {
  stations[[v]] <- stations$station_abbr %in% var_lists[[v]]
}
stations$n_vars <- rowSums(as.matrix(stations[names(var_lists)]))

# codes in the lists that do not exist in the metadata (should be none)
unknown <- setdiff(unlist(var_lists), stations$station_abbr)
if (length(unknown) > 0) {
  warning("Codes not found in station metadata: ", paste(unknown, collapse = ", "))
}


# ---- 5. Popup and label text -------------------------------------
tick  <- '<span style="color:#2a9d3f">&#10004;</span>'
cross <- '<span style="color:#bbbbbb">&#10008;</span>'

status <- sapply(names(var_labels), function(v) ifelse(stations[[v]], tick, cross))
status <- matrix(status, nrow = nrow(stations))   # one row per station

var_lines <- apply(status, 1, function(r) {
  paste0(r, " ", names(var_labels), " (", unname(var_labels), ")", collapse = "<br>")
})

stations$popup <- paste0(
  "<b>", htmlEscape(stations$station_name), "</b> (", stations$station_abbr, ")<br>",
  htmlEscape(stations$station_canton), ", ", stations$station_height_masl, " m a.s.l.<br>",
  htmlEscape(stations$station_type_en),
  "<hr style='margin:4px 0'>", var_lines
)
stations$label <- paste0(stations$station_name, " (", stations$station_abbr, "): ",
                         stations$n_vars, " of ", length(var_labels), " variables")

# ---- 6. Map -------------------------------------------------------
lv  <- as.character(0:length(var_labels))
pal <- colorFactor("viridis", levels = lv)

var_cols <- c(tas = "#e6550d", tasmax = "#d62728", tasmin = "#1f77b4", pr = "#08519c",
              rsds = "#e6ab02", hurs = "#17becf", sfcWind = "#7f7f7f")

grp_all  <- "All stations"
grp_core <- "pr + tas + rsds all available"
grp_var  <- paste0(names(var_labels), " (", unname(var_labels), ")")
names(grp_var) <- names(var_labels)

m <- leaflet(stations, options = leafletOptions(minZoom = 7)) |>
  addProviderTiles(providers$Esri.WorldGrayCanvas) |>
  addCircleMarkers(
    lng = ~station_coordinates_wgs84_lon, lat = ~station_coordinates_wgs84_lat,
    radius = 5, weight = 1, color = "white", fillColor = ~pal(as.character(n_vars)),
    fillOpacity = 0.9, popup = ~popup, label = ~label, group = grp_all
  )

# a coloured ring around every station that has the variable
for (v in names(var_labels)) {
  d <- stations[stations[[v]], ]
  m <- m |>
    addCircleMarkers(
      data = d, lng = ~station_coordinates_wgs84_lon, lat = ~station_coordinates_wgs84_lat,
      radius = 9, weight = 2.5, color = unname(var_cols[v]), fill = FALSE,
      popup = ~popup, label = ~label, group = unname(grp_var[v])
    )
}

core <- stations[stations$pr & stations$tas & stations$rsds, ]
m <- m |>
  addCircleMarkers(
    data = core, lng = ~station_coordinates_wgs84_lon, lat = ~station_coordinates_wgs84_lat,
    radius = 9, weight = 2.5, color = "black", fill = FALSE,
    popup = ~popup, label = ~label, group = grp_core
  ) |>
  addLayersControl(
    overlayGroups = c(grp_all, grp_core, unname(grp_var)),
    options       = layersControlOptions(collapsed = FALSE)
  ) |>
  hideGroup(c(grp_core, unname(grp_var))) |>
  addLegend(pal = pal, values = lv, position = "bottomright",
            title = paste0("Variables available<br>(of ", length(var_labels), ")")) |>
  addControl("<b>CH2025 stations: variable availability</b><br>
              Click a station for details; tick the ring layers to filter by variable",
             position = "topright") |>
  addScaleBar(position = "bottomleft")

# ---- 7. Save and show --------------------------------------------
# When knitting, the map object `m` is used directly in the document,
# so the HTML file is only written when the script is run on its own.
if (!knitting) {
  dir.create(dirname(out_html), showWarnings = FALSE, recursive = TRUE)
  # selfcontained = TRUE needs pandoc (included with RStudio)
  saveWidget(m, out_html, selfcontained = TRUE,
             title = "CH2025 station variable availability")
  message("Map saved to: ", normalizePath(out_html))
  print(m)   # show in the RStudio viewer
}
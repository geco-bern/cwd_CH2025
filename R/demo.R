library(rstac)

catalog <- stac("https://data.geo.admin.ch/api/stac/v1/")
query <- stac_search(catalog,
                     collections = c("ch.meteoschweiz.ogd-climate-scenarios-ch2025-grid"),
                     bbox = c(7.532035947766065, 46.71592755296092, 7.7209811162807345, 46.79608084001828))

result <- get_request(query)

entries <- result[["features"]]
if (!is.null(entries) && length(entries) > 0) {
  for (entry in entries) {
    if (!is.null(entry$id)) {
      cat(entry$id, "\n")
    }
  }
}
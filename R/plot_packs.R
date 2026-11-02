#' @title Plot Packs Interactively
#' @description Interactive map of the packs: kernel territories (95% isopleth, outlined) and core areas
#'   (50% isopleth, filled), and the samples coloured by pack. Individuals without a pack are drawn in white.
#'   Clicking a sample shows the individual, its pack and, when there are such columns, its status and date;
#'   clicking a territory, its pack, level, area and members.
#' @param obs An `sf` object of points (one row per sample) with `Individual` and the pack column.
#' @param pack Name of the pack column (default `"Pack"`; a bare name works too, as in older scripts).
#'   `NA` or `"Lone Individual"` = no pack.
#' @param ud A `wolfpack_ud` object from [space_use()]; computed from `obs` if `NULL`.
#' @param h Bandwidth passed to [space_use()] when `ud` is `NULL`.
#' @param levels Isopleth levels drawn (default: 95% territory and 50% core area).
#' @param tiles Background map(s): names of `leaflet.providers` (e.g. `"Esri.WorldTopoMap"`, `"CartoDB.Positron"`),
#'   `"OpenStreetMap"`, or tile URL templates (`"https://.../{z}/{x}/{y}.png"`). With several, the first one is shown
#'   and a selector lets you switch; the names of the vector are the labels of the selector.
#'   `"GeoportailFrance.plan"` and `"GeoportailFrance.orthos"` (IGN map and aerial photos, France) use the IGN
#'   Geoplateforme whatever the version of `leaflet.providers` (up to 2.0.0 it points to a closed server).
#'   OpenStreetMap tiles are refused (403) when the map is saved with `htmlwidgets::saveWidget()` and opened from the
#'   disk: their servers require a Referer, which a local file doesn't send. Default: Esri topographic map.
#' @param exclusive Draw the territories without overlap: where the isopleths of several packs overlap, the space goes
#'   to the pack using it most (see [pack_territories()]). Default: `FALSE`, the plain isopleths, which spill over
#'   into the neighbours' territories by about the bandwidth h.
#' @return An interactive leaflet map.
#' @importFrom sf st_is_longlat st_coordinates st_transform
#' @importFrom leaflet leaflet colorFactor addTiles addProviderTiles addLayersControl tileOptions addPolygons fitBounds addCircleMarkers addLegend
#' @importFrom dplyr %>%
#' @importFrom grDevices colorRampPalette
#' @export
#' @examples
#' data("samples_Carpathians"); data("relate_Carpathians")
#' res <- wolf_packs(samples_Carpathians, relate_Carpathians, "QG", 0.25)
#' plot_packs(res$samples, "Pack", ud = res$ud)
#' # several backgrounds, with a selector
#' plot_packs(res$samples, "Pack", ud = res$ud,
#'            tiles = c("Topo" = "Esri.WorldTopoMap", "Aerial" = "Esri.WorldImagery"))
#' # territories without overlap between neighbours
#' plot_packs(res$samples, "Pack", ud = res$ud, exclusive = TRUE)
plot_packs <- function(obs, pack = "Pack", ud = NULL, h = NULL, levels = c(0.95, 0.5),
                       tiles = c("Topo (Esri)" = "Esri.WorldTopoMap"), exclusive = FALSE) {
  pack <- .col_name(substitute(pack), function() pack, names(obs))

  # Checks
  if (!inherits(obs, "sf")) stop("The 'obs' object must be of class 'sf'.")
  if (nrow(obs) == 0) stop("The 'obs' object is empty.")
  if (!pack %in% names(obs)) stop("Column '", pack, "' doesn't exist.")
  if (!is.character(tiles) || !length(tiles)) stop("'tiles' must be one or several provider names or tile URLs.")
  lab <- as.character(obs[[pack]])
  lab[lab %in% "Lone Individual"] <- NA
  packs <- sort(unique(stats::na.omit(lab)))
  if (is.null(ud)) ud <- space_use(obs, h = h)

  # One colour per pack ("Set1", interpolated beyond 9 packs)
  n_packs <- length(packs)
  pal_cols <- if (n_packs <= 9) RColorBrewer::brewer.pal(max(3, n_packs), "Set1")[seq_len(max(1, n_packs))] else
    colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_packs)
  colors <- colorFactor(palette = pal_cols, domain = if (n_packs) packs else "none")

  # Map in WGS84, on the chosen background(s)
  pts <- if (isTRUE(st_is_longlat(obs))) obs else st_transform(obs, 4326)
  coords <- st_coordinates(pts)
  l <- .add_tiles(leaflet(), tiles)

  # Territories first, under the points: first level (95%) outlined with a light fill, the others (50%) filled
  if (n_packs) {
    members <- lapply(stats::setNames(packs, packs), function(p) unique(as.character(obs$Individual[lab %in% p])))
    terr <- st_transform(pack_territories(ud, members, levels, exclusive), 4326)
    terr$Individuals <- vapply(members[terr$Pack], paste, character(1), collapse = ", ")
    for (k in seq_len(nrow(terr))) {
      p <- terr$Pack[k]
      l <- l %>% addPolygons(data = terr[k, ], color = colors(p), fillColor = colors(p), weight = 2,
                             fillOpacity = if (terr$level[k] == levels[1]) 0.1 else 0.35, opacity = 0.8,
                             popup = ~paste0("Pack: ", Pack, " (", round(100 * level), "%, ", round(area_km2), " km2)",
                                             "<br>Individuals: ", Individuals))
    }
  }

  # Samples coloured by pack; individuals without a pack in white with a black border
  df <- data.frame(Individual = as.character(obs$Individual), Pack = lab, lon = coords[, 1], lat = coords[, 2],
                   stringsAsFactors = FALSE)
  df$info <- paste0("Individual: ", df$Individual, "<br>Pack: ", ifelse(is.na(df$Pack), "none", df$Pack))
  if ("Status" %in% names(obs)) df$info <- paste0(df$info, "<br>Status: ", obs$Status)
  if ("Date" %in% names(obs)) df$info <- paste0(df$info, "<br>Date: ", obs$Date)
  in_pack <- !is.na(df$Pack)
  if (any(in_pack)) {
    df_in <- df[in_pack, ]
    df_in$Color <- colors(df_in$Pack)
    l <- l %>% addCircleMarkers(data = df_in, lng = ~lon, lat = ~lat, radius = 4, color = ~Color,
                                stroke = FALSE, fillOpacity = 0.9, popup = ~info)
  }
  if (any(!in_pack)) {
    l <- l %>% addCircleMarkers(data = df[!in_pack, ], lng = ~lon, lat = ~lat, radius = 4, color = "black",
                                weight = 1, stroke = TRUE, fillColor = "white", fillOpacity = 0.9, popup = ~info)
  }

  # Frame and legend
  l <- l %>% fitBounds(lng1 = min(coords[, 1]), lat1 = min(coords[, 2]), lng2 = max(coords[, 1]), lat2 = max(coords[, 2]))
  leg_lab <- c(packs, if (any(!in_pack)) "No pack")
  leg_col <- c(if (n_packs) colors(packs), if (any(!in_pack)) "white")
  l %>% addLegend(position = "bottomright", colors = leg_col, labels = leg_lab, title = "Packs")
}

# IGN layers given by their URL on the Geoplateforme (open data, no key), so they don't depend on leaflet.providers
.ign_tiles <- c(
  GeoportailFrance.plan = paste0("https://data.geopf.fr/wmts?SERVICE=WMTS&REQUEST=GetTile&VERSION=1.0.0",
                                 "&LAYER=GEOGRAPHICALGRIDSYSTEMS.PLANIGNV2&STYLE=normal&FORMAT=image/png",
                                 "&TILEMATRIXSET=PM&TILEMATRIX={z}&TILEROW={y}&TILECOL={x}"),
  GeoportailFrance.orthos = paste0("https://data.geopf.fr/wmts?SERVICE=WMTS&REQUEST=GetTile&VERSION=1.0.0",
                                   "&LAYER=ORTHOIMAGERY.ORTHOPHOTOS&STYLE=normal&FORMAT=image/jpeg",
                                   "&TILEMATRIXSET=PM&TILEMATRIX={z}&TILEROW={y}&TILECOL={x}")
)

# Background maps of plot_packs(): IGN by URL, URL templates, plain OSM or leaflet.providers; a selector if several
.add_tiles <- function(l, tiles) {
  url <- grepl("^https?://", tiles)
  groups <- names(tiles)
  if (is.null(groups)) groups <- rep("", length(tiles))
  auto <- ifelse(url, paste("Map", seq_along(tiles)), tiles)          # no name: the provider name, or "Map 2" for a URL
  groups[is.na(groups) | groups == ""] <- auto[is.na(groups) | groups == ""]
  for (k in seq_along(tiles)) {
    t <- unname(tiles[k])
    l <- if (t %in% names(.ign_tiles)) {
      addTiles(l, .ign_tiles[[t]], group = groups[k], options = tileOptions(maxZoom = 19),
               attribution = '<a target="_blank" href="https://www.geoportail.gouv.fr/">IGN - Geoportail</a>')
    } else if (url[k]) {
      addTiles(l, t, group = groups[k])
    } else if (t == "OpenStreetMap") {
      addTiles(l, group = groups[k])                                  # fine in RStudio, refused in a saved html file
    } else addProviderTiles(l, t, group = groups[k])
  }
  if (length(tiles) > 1) l <- addLayersControl(l, baseGroups = groups, position = "topright")
  l
}

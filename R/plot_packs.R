#' @title Plot Packs Interactively
#' @description Interactive map of the packs: kernel territories (95% isopleth, outlined) and core areas
#'   (50% isopleth, filled), and the samples coloured by pack. Individuals without a pack are drawn in white.
#'   Clicking a sample shows the individual, its pack and, when there are such columns, its status and date.
#' @param obs An `sf` object of points (one row per sample) with `Individual` and the pack column.
#' @param pack Name of the pack column (default `"Pack"`; a bare name works too, as in older scripts).
#'   `NA` or `"Lone Individual"` = no pack.
#' @param ud A `wolfpack_ud` object from [space_use()]; computed from `obs` if `NULL`.
#' @param h Bandwidth passed to [space_use()] when `ud` is `NULL`.
#' @param levels Isopleth levels drawn (default: 95% territory and 50% core area).
#' @return An interactive leaflet map.
#' @importFrom sf st_is_longlat st_coordinates st_transform
#' @importFrom leaflet leaflet colorFactor addTiles addPolygons fitBounds addCircleMarkers addLegend
#' @importFrom dplyr %>%
#' @importFrom grDevices colorRampPalette
#' @export
#' @examples
#' data("samples_Carpathians"); data("relate_Carpathians")
#' res <- wolf_packs(samples_Carpathians, relate_Carpathians, "QG", 0.25)
#' plot_packs(res$samples, "Pack", ud = res$ud)
plot_packs <- function(obs, pack = "Pack", ud = NULL, h = NULL, levels = c(0.95, 0.5)) {
  pack <- .col_name(substitute(pack), function() pack, names(obs))

  # Checks
  if (!inherits(obs, "sf")) stop("The 'obs' object must be of class 'sf'.")
  if (nrow(obs) == 0) stop("The 'obs' object is empty.")
  if (!pack %in% names(obs)) stop("Column '", pack, "' doesn't exist.")
  lab <- as.character(obs[[pack]])
  lab[lab %in% "Lone Individual"] <- NA
  packs <- sort(unique(stats::na.omit(lab)))
  if (is.null(ud)) ud <- space_use(obs, h = h)

  # One colour per pack ("Set1", interpolated beyond 9 packs)
  n_packs <- length(packs)
  pal_cols <- if (n_packs <= 9) RColorBrewer::brewer.pal(max(3, n_packs), "Set1")[seq_len(max(1, n_packs))] else
    colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_packs)
  colors <- colorFactor(palette = pal_cols, domain = if (n_packs) packs else "none")

  # Map in WGS84
  pts <- if (isTRUE(st_is_longlat(obs))) obs else st_transform(obs, 4326)
  coords <- st_coordinates(pts)
  l <- leaflet() %>% addTiles()

  # Territories first, under the points: 95% outlined with a light fill, 50% filled
  for (p in packs) {
    members <- unique(as.character(obs$Individual[lab %in% p]))
    iso <- st_transform(isopleth(ud, members, levels), 4326)
    iso$Pack <- p
    iso$Individuals <- paste(members, collapse = ", ")
    for (k in seq_len(nrow(iso))) {
      l <- l %>% addPolygons(data = iso[k, ], color = colors(p), fillColor = colors(p), weight = 2,
                             fillOpacity = if (k == 1) 0.1 else 0.35, opacity = 0.8,
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

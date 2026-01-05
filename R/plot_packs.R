#' @title Plot Packs Interactively
#' @description Creates an interactive map visualizing the packs territories.
#' @param obs An `sf` object with individuals and their genetic groups.
#' @param ud_df A data.frame containing individuals, their genetic group, and their assigned pack.
#' @return An interactive leaflet map.
#' @importFrom sf st_is_longlat st_coordinates st_polygon st_sf st_sfc
#' @importFrom leaflet colorFactor addTiles addPolygons fitBounds addCircleMarkers addLegend
#' @export
#' @examples
#' map <- plot_groups(obs, ud_df = ud_df)
#' map
plot_packs <- function(obs, ud_df) {
  suppressPackageStartupMessages({
    library(leaflet)
    library(sf)
    library(RColorBrewer)
  })

  # Vérification des données d'entrée
  if (!inherits(obs, "sf")) {
    stop("The 'obs' object must be of class 'sf'.")
  }

  if (nrow(obs) == 0) {
    stop("The 'obs' object is empty.")
  }

  # Create a color palette for packs
  packs <- unique(ud_df$Pack[ud_df$Pack != "Lone Individual"])
  colors <- colorFactor(palette = "Set1", domain = packs)

  # Prepare data for leaflet
  obs$Pack <- ud_df$Pack[match(obs$Individual, ud_df$Individual)]
  obs$Color <- ifelse(obs$Pack == "Lone Individual", "red", colors(obs$Pack))

  # Ensure coordinates are in WGS84 (EPSG:4326)
  if (!st_is_longlat(obs)) {
    obs <- tryCatch(st_transform(obs, 4326), error = function(e) {
      message("Failed to transform coordinates to WGS84: ", e$message)
      return(obs)
    })
  }

  # Create the leaflet map
  l <- leaflet() %>%
    addTiles()

  # MCP for each pack (added first, under the points)
  for (pack in packs) {
    pack_obs <- obs[obs$Pack == pack, ]
    if (nrow(pack_obs) >= 3) {
      coords_pack <- st_coordinates(pack_obs)
      if (nrow(coords_pack) >= 3) {
        hull <- tryCatch({
          chull(coords_pack)
        }, error = function(e) {
          NULL
        })

        if (!is.null(hull) && length(hull) >= 3) {
          hull_coords <- rbind(coords_pack[hull, ], coords_pack[hull[1], ])
          polygon <- st_polygon(list(hull_coords))
          mcp_i <- st_sf(geometry = st_sfc(polygon), Pack = pack, Color = colors(pack))

          l <- l %>% addPolygons(
            data = mcp_i,
            fillColor = colors(pack),
            color = colors(pack),
            fillOpacity = 0.3,
            weight = 2,
            popup = pack
          )
        }
      }
    }
  }

  # Fit the map to the bounds of the observations
  coords <- st_coordinates(obs)
  l <- l %>% fitBounds(lng1 = min(coords[, 1]), lat1 = min(coords[, 2]),
                       lng2 = max(coords[, 1]), lat2 = max(coords[, 2]))

  # Add points to the map, colored by pack
  obs_df <- data.frame(
    Individual = obs$Individual,
    Pack = obs$Pack,
    Color = obs$Color,
    lon = coords[, 1],
    lat = coords[, 2],
    stringsAsFactors = FALSE
  )

  if (nrow(obs_df) > 0) {
    # Add regular circle markers for pack individuals
    l <- l %>% addCircleMarkers(
      lng = ~lon,
      lat = ~lat,
      radius = 3,
      color = ~Color,
      stroke = FALSE,
      fillOpacity = 0.8,
      data = obs_df[obs_df$Pack != "Lone Individual", ],
      popup = ~paste("Individual:", Individual, "<br>Pack:", Pack)
    )

    # Add circle markers with border for Lone Individuals
    lone_individuals_df <- obs_df[obs_df$Pack == "Lone Individual", ]
    if (nrow(lone_individuals_df) > 0) {
      l <- l %>% addCircleMarkers(
        lng = ~lon,
        lat = ~lat,
        radius = 3,
        color = "black",
        stroke = TRUE,
        fillOpacity = 0.8,
        weight = 1,
        fillColor = "white",
        data = lone_individuals_df,
        popup = ~paste("Individual:", Individual, "<br>Pack: Lone Individual")
      )
    }
  } else {
    warning("No points to plot.")
  }

  # Add legend
  pack_names <- c(packs, "Lone Individual")
  pack_colors <- c(colors(packs), "white")
  l <- l %>% addLegend(
    position = "bottomright",
    colors = pack_colors,
    labels = pack_names,
    title = "Packs"
  )

  return(l)

}

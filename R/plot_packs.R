#' @title Plot Packs Interactively
#' @description Creates an interactive map visualizing the packs territories.
#' @param obs An `sf` object with individuals and their packs.
#' @param pack The name of the column in `obs` that contains the pack information.
#' @return An interactive leaflet map.
#' @importFrom sf st_is_longlat st_coordinates st_polygon st_sf st_sfc st_as_sf
#' @importFrom leaflet colorFactor addTiles addPolylines addPolygons fitBounds addCircleMarkers addLegend
#' @export
#' @examples
#' data(samples)
#' samples$pack <- NA
#' samples[samples$Individual == "W1",]$pack <- "1"
#' samples[samples$Individual == "W2",]$pack <- "1"
#' samples[samples$Individual == "W3",]$pack <- "1"
#' samples[samples$Individual == "W4",]$pack <- "Lone Individual"
#' samples[samples$Individual == "W5",]$pack <- "2"
#' samples[samples$Individual == "W6",]$pack <- "3"
#' samples[samples$Individual == "W7",]$pack <- "2"
#' samples[samples$Individual == "W8",]$pack <- "2"
#' samples[samples$Individual == "W9",]$pack <- "3"
#' map <- plot_packs(samples, pack = pack)
#' map
plot_packs <- function(obs, pack) {
  pack_name <- deparse(substitute(pack))

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
  packs <- unique(obs[[pack_name]][obs[[pack_name]] != "Lone Individual"])
  colors <- colorFactor(palette = "Set1", domain = packs)

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
    pack_obs <- obs[obs[[pack_name]] == pack, ]
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
          mcp_sf <- st_sf(geometry = st_sfc(polygon), Pack = pack, Individuals = paste(unique(pack_obs$Individual), collapse = ", "))

          # Use the sf object directly in addPolygons
          l <- l %>% addPolygons(
            data = mcp_sf,
            fillColor = colors(pack),
            color = colors(pack),
            fillOpacity = 0.3,
            weight = 2,
            popup = ~paste("Pack:", Pack, "<br>Individuals:", Individuals)
          )
        }
      }
    }
  }

  # Fit the map to the bounds of the observations
  coords <- st_coordinates(obs)
  l <- l %>% fitBounds(lng1 = min(coords[, 1]), lat1 = min(coords[, 2]),
                       lng2 = max(coords[, 1]), lat2 = max(coords[, 2]))

  # Prepare data for leaflet points
  obs_df <- data.frame(
    Individual = obs$Individual,
    Pack = obs[[pack_name]],
    lon = coords[, 1],
    lat = coords[, 2],
    stringsAsFactors = FALSE
  )

  # Assign colors based on pack
  obs_df$Color <- ifelse(obs_df$Pack == "Lone Individual", "red", colors(obs_df$Pack))

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

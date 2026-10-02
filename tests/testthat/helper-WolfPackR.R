# Shared helpers for the tests: toy data built from the datasets shipped with the package

# Load a dataset of the package and hand it back (no LazyData, so no WolfPackR::samples)
pkg_data <- function(name) {
  env <- new.env()
  utils::data(list = name, package = "WolfPackR", envir = env)
  env[[name]]
}

# Samples + genetic groups, i.e. the usual starting point of the workflow
toy_obs <- function() {
  grDevices::pdf(NULL)                 # genetic_groups() draws the graph, keep it off-screen
  on.exit(grDevices::dev.off())
  gg <- genetic_groups(pkg_data("relate"), threshold = 0.4, estimator = "wang")
  merge(pkg_data("samples"), gg, by = "Individual")
}

# n made-up samples of one individual, scattered around (x, y) in lon/lat,
# with the same columns and geometry column name as `template`
fake_ind <- function(template, id, group, n, x, y, d = 0.005) {
  set.seed(42)
  pts <- sf::st_as_sf(data.frame(X = x + stats::runif(n, -d, d), Y = y + stats::runif(n, -d, d)),
                      coords = c("X", "Y"), crs = sf::st_crs(template))
  new <- template[rep(1, n), ]          # copy the structure of the template
  new$Individual <- id
  new$group <- group
  sf::st_geometry(new) <- sf::st_geometry(pts)
  new
}

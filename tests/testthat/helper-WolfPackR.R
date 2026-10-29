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
  gg <- genetic_groups(pkg_data("relate_example"), threshold = 0.4, estimator = "wang")
  merge(pkg_data("samples_example"), gg, by = "Individual")
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

# Calls a deprecated function without its deprecation warning (other warnings still go through)
quiet_deprecated <- function(expr) withCallingHandlers(expr, deprecatedWarning = function(w) invokeRestart("muffleWarning"))

# Simulated study area in metres (UTM 32N): family A around (0, 0), family B around (20 km, 0), plus special cases
#   D  related to A1 but sampled in B's area      -> putative disperser
#   U  related to nobody, sampled in A's area      -> unrelated resident
#   T  related to nobody, scattered far away       -> transient
#   L  related to A2, a single sample in B's area  -> too little data
#   M  related to A3, in A in P1 then in B in P2   -> disperser (dated)
sim_samples <- function() {
  set.seed(11)
  pt <- function(id, n, x, y, sd = 1500, per = "P1") data.frame(Individual = id, X = x + stats::rnorm(n, 0, sd),
                                                               Y = y + stats::rnorm(n, 0, sd), Campaign = per)
  d <- rbind(
    pt("A1", 4, 0, 0), pt("A1", 4, 0, 0, per = "P2"), pt("A2", 8, 0, 0), pt("A3", 6, 0, 0, per = "P2"),
    pt("B1", 4, 20000, 0), pt("B1", 4, 20000, 0, per = "P2"), pt("B2", 8, 20000, 0), pt("B3", 6, 20000, 0, per = "P2"),
    pt("D", 6, 20000, 0), pt("U", 5, 0, 0), pt("T", 5, 10000, 35000, sd = 7000), pt("L", 1, 20000, 0),
    pt("M", 3, 0, 0, per = "P1"), pt("M", 5, 20000, 0, per = "P2"))
  d$Sex <- ifelse(d$Individual %in% c("A1", "B1", "D", "T", "M"), "M", "F")
  sf::st_as_sf(d, coords = c("X", "Y"), crs = 32632)
}
sim_relate <- function() {
  ids <- c("A1", "A2", "A3", "B1", "B2", "B3", "D", "U", "T", "L", "M")
  p <- t(utils::combn(ids, 2))
  r <- data.frame(ind1 = p[, 1], ind2 = p[, 2], r = 0)
  kin <- rbind(c("A1", "A2"), c("A1", "A3"), c("A2", "A3"), c("B1", "B2"), c("B1", "B3"), c("B2", "B3"),
               c("A1", "D"), c("A2", "L"), c("A3", "M"))
  for (k in seq_len(nrow(kin))) r$r[(r$ind1 == kin[k, 1] & r$ind2 == kin[k, 2]) | (r$ind1 == kin[k, 2] & r$ind2 == kin[k, 1])] <- 0.5
  r
}

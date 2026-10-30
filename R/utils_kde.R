# Internal helpers shared by the KDE functions: projection, kernel, grid cells, spread of an individual, seeds.
# Nothing here is exported.

# Projects the samples to a metric CRS: UTM zone of the centre of the data when they're in lon/lat,
# left as they are when already projected (or without CRS, assumed to be in metres)
.as_metric <- function(obs) {
  crs <- sf::st_crs(obs)
  if (is.na(crs) || !isTRUE(sf::st_is_longlat(obs))) return(obs)
  ctr <- colMeans(sf::st_coordinates(obs)[, 1:2, drop = FALSE])           # rough centre is plenty to pick the zone
  zone <- floor((ctr[1] + 180) / 6) + 1
  sf::st_transform(obs, if (ctr[2] >= 0) 32600 + zone else 32700 + zone)
}

# Gaussian kernel density of points `xy` on the cell centres `grid`, normalised to sum to 1 over the grid
.kde <- function(xy, grid, h) {
  u <- numeric(nrow(grid))
  for (k in seq_len(nrow(xy))) {                                         # one point at a time: no cells x points matrix in memory
    u <- u + exp(-((grid[, 1] - xy[k, 1])^2 + (grid[, 2] - xy[k, 2])^2) / (2 * h^2))
  }
  u / sum(u)
}

# Index of the grid cell holding each point (points off the grid go to the nearest edge cell)
.cell <- function(ud, xy) {
  ix <- pmin(pmax(round((xy[, 1] - ud$x0) / ud$res), 0), ud$nx - 1)
  iy <- pmin(pmax(round((xy[, 2] - ud$y0) / ud$res), 0), ud$ny - 1)
  iy * ud$nx + ix + 1
}

# Radius of gyration (m): typical distance of the points to their centre; 0 with fewer than 2 points.
# With trim > 0, the farthest share of the points (rounded up, but keeping at least 3) is left out first
.rg <- function(xy, trim = 0) {
  if (nrow(xy) < 2) return(0)
  if (trim > 0 && nrow(xy) > 3) {
    k <- min(ceiling(trim * nrow(xy)), nrow(xy) - 3)                       # how many points to leave out
    d <- rowSums(sweep(xy, 2, colMeans(xy))^2)
    xy <- xy[order(d)[seq_len(nrow(xy) - k)], , drop = FALSE]
  }
  sqrt(mean(rowSums(sweep(xy, 2, colMeans(xy))^2)))
}

# Wide-ranging individuals: radius of gyration above `wide` x h, leaving out the farthest 10% of the samples
# (at least one, keeping at least 3), so that a single foray into a neighbouring territory doesn't make a resident
# look wide-ranging, while an animal spread over several places still is. With a period column, the spread is
# computed within each period having 3+ samples, so an animal that moved from one place to another between periods
# is not flagged (and 2 samples far apart in a season don't make a wide-ranger either)
.wide_ranging <- function(ud, wide, period = NULL, trim = 0.1) {
  s <- ud$samples
  rg <- sapply(ud$ids, function(i) {
    w <- which(s$Individual == i)
    pooled <- .rg(as.matrix(s[w, c("x", "y"), drop = FALSE]), trim)
    if (is.null(period)) return(pooled)
    per <- split(w, s[[period]][w])
    per <- per[lengths(per) >= 3]
    if (!length(per)) return(pooled)                                       # no well-sampled period: pooled spread
    max(vapply(per, function(r) .rg(as.matrix(s[r, c("x", "y"), drop = FALSE]), trim), numeric(1)))
  })
  names(rg) <- ud$ids
  rg > wide * ud$h
}

# Utilisation distribution of a set of individuals: their UDs mixed, each weighted by its number of samples
.mix_ud <- function(ud, members) {
  members <- intersect(members, ud$ids)
  if (!length(members)) return(rep(0, nrow(ud$grid)))
  w <- ud$n[members]
  as.vector(ud$U[, members, drop = FALSE] %*% w) / sum(w)
}

# Majority rule: the label holding more than min_share of the samples (a tie at the top never wins)
.majority <- function(labels, min_share = 0.5) {
  s <- prop.table(table(labels))
  top <- s[s == max(s)]
  if (length(top) == 1 && top >= min_share) names(top) else NA_character_
}

# Text summary of where the samples went, e.g. "3: 5, 4: 2, outside: 1"
.split_text <- function(labels) {
  t <- sort(table(labels), decreasing = TRUE)
  paste(names(t), t, sep = ": ", collapse = ", ")
}

# Runs `expr` with a given seed and puts the user's random number stream back afterwards
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  expr
}

# Column name given as a string, as a bare column name (older scripts did that) or as a variable holding the name
.col_name <- function(sub, value, cols) {
  if (is.character(sub)) return(sub)
  if (is.symbol(sub) && deparse(sub) %in% cols) return(deparse(sub))
  value()                                                                 # evaluated only now, e.g. p <- "Pack"
}

# Kin links between individuals: logical matrix, TRUE when relatedness > threshold
.kin_matrix <- function(relate, estimator, threshold, ids) {
  K <- matrix(FALSE, length(ids), length(ids), dimnames = list(ids, ids))
  r <- relate[as.character(relate$ind1) %in% ids & as.character(relate$ind2) %in% ids, ]
  r <- r[!is.na(r[[estimator]]) & r[[estimator]] > threshold & as.character(r$ind1) != as.character(r$ind2), ]
  K[cbind(as.character(r$ind1), as.character(r$ind2))] <- TRUE
  K | t(K)
}

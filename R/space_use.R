#' @title Space use of each individual (kernel density)
#' @description Estimates the utilisation distribution (UD) of every individual with a Gaussian kernel, on a grid
#'   shared by all individuals and with one bandwidth `h` for everybody. Because the bandwidth is common,
#'   an individual with a single sample still gets a UD (a bump of width `h` around its sample), so individuals
#'   with few samples are handled like the others, without special cases.
#'
#'   Everything is done with `sf` and base R. The result is the engine of [ud_overlap()], [isopleth()],
#'   [attribute_samples()], [assign_individuals()] and [spatial_groups()].
#' @param obs An `sf` object of points (one row per sample) with an `Individual` column. Data in lon/lat are
#'   projected to the UTM zone of the study area; projected data are used as they are (units: metres).
#' @param h Kernel bandwidth in metres. `NULL` (default): median of the reference bandwidths (href) of the
#'   individuals with at least `min_href` samples. You can also set it from the expected territory size:
#'   it is the scale at which samples are considered to come from the same area.
#' @param res Grid resolution in metres (default: `h / 4`).
#' @param buffer Margin around the samples, in multiples of `h` (default: 3).
#' @param min_href Minimum number of samples of an individual to use its href when `h = NULL` (default: 10).
#' @return An object of class `wolfpack_ud`: a list with `U` (matrix, one column per individual, each column
#'   sums to 1 over the grid), `grid` (cell centres), `h`, `res`, `crs`, `n` (samples per individual), `ids`
#'   and `samples` (the samples with their projected coordinates and grid cell).
#' @references Worton, B. J. (1989). Kernel methods for estimating the utilization distribution in home-range
#'   studies. Ecology, 70, 164-168.
#' @export
#' @examples
#' data("samples_Carpathians")
#' ud <- space_use(samples_Carpathians)
#' ud
space_use <- function(obs, h = NULL, res = NULL, buffer = 3, min_href = 10) {

  # Checks and projection
  if (!inherits(obs, "sf")) stop("'obs' must be an sf object of points.")
  if (!"Individual" %in% names(obs)) stop("'obs' needs an 'Individual' column.")
  if (nrow(obs) == 0) stop("'obs' is empty.")
  proj <- .as_metric(obs)
  xy <- sf::st_coordinates(proj)[, 1:2, drop = FALSE]
  ind <- as.character(proj$Individual)
  ids <- sort(unique(ind))
  n <- table(factor(ind, levels = ids))
  n <- stats::setNames(as.integer(n), ids)

  # Bandwidth: median href of the well-sampled individuals, unless given
  if (is.null(h)) {
    good <- ids[n >= min_href]
    if (!length(good)) good <- ids[n >= 3]
    if (!length(good)) stop("Not enough samples per individual to choose 'h': please give it.")
    href <- vapply(good, function(i) {
      p <- xy[ind == i, , drop = FALSE]
      sqrt(0.5 * (stats::var(p[, 1]) + stats::var(p[, 2]))) * nrow(p)^(-1 / 6)
    }, numeric(1))
    h <- stats::median(href[href > 0])
    message(sprintf("Bandwidth h = %.0f m (median href of %d individuals with >= %d samples).", h, length(href), min_href))
  }
  if (is.null(res)) res <- h / 4

  # Common grid: cell centres covering the samples plus a margin
  x0 <- min(xy[, 1]) - buffer * h; y0 <- min(xy[, 2]) - buffer * h
  nx <- ceiling((max(xy[, 1]) + buffer * h - x0) / res) + 1
  ny <- ceiling((max(xy[, 2]) + buffer * h - y0) / res) + 1
  grid <- cbind(x = x0 + res * rep(seq_len(nx) - 1, ny), y = y0 + res * rep(seq_len(ny) - 1, each = nx))

  # One UD per individual, as the columns of a matrix
  U <- vapply(ids, function(i) .kde(xy[ind == i, , drop = FALSE], grid, h), numeric(nrow(grid)))
  if (is.null(dim(U))) U <- matrix(U, ncol = 1, dimnames = list(NULL, ids))

  samples <- sf::st_drop_geometry(proj)
  samples$Individual <- ind
  samples$x <- xy[, 1]; samples$y <- xy[, 2]
  samples$row <- seq_len(nrow(samples))                  # position in obs, to join back
  ud <- structure(list(U = U, grid = grid, x0 = x0, y0 = y0, nx = nx, ny = ny, res = res, h = h,
                       crs = sf::st_crs(proj), n = n, ids = ids, samples = samples), class = "wolfpack_ud")
  ud$samples$cell <- .cell(ud, xy)
  ud
}

#' @export
print.wolfpack_ud <- function(x, ...) {
  cat(sprintf("Space use (kernel UD) of %d individuals from %d samples\n", length(x$ids), nrow(x$samples)))
  cat(sprintf("Bandwidth h = %.0f m, grid of %d x %d cells of %.0f m\n", x$h, x$nx, x$ny, x$res))
  invisible(x)
}

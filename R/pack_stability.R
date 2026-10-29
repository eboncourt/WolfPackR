#' @title How robust are the packs?
#' @description Runs [wolf_packs()] many times and measures how much the packs change:
#'   * **sensitivity**: one run per combination of the settings given in `vary` (e.g. several bandwidths and
#'     relatedness thresholds), everything else at its base value;
#'   * **bootstrap**: `n_boot` runs where the samples of each individual are resampled with replacement
#'     (same number of samples per individual), to see the effect of sampling chance.
#'
#'   For each individual, `stability` is the mean, over all runs, of the Jaccard index between its pack mates in
#'   the base run and in that run (1 = always with the same companions). `coassignment` gives, for every pair,
#'   the share of runs where both are in the same pack.
#' @param obs,relate,estimator,threshold As in [wolf_packs()].
#' @param vary Named list of settings to vary, e.g. `list(h = c(1500, 2000, 2500), threshold = c(0.2, 0.25, 0.3))`.
#'   Any argument of [wolf_packs()] works (`h`, `threshold`, `min_share`, `merge`, `wide`, `method`...).
#' @param n_boot Number of bootstrap runs (default: 0).
#' @param reference Optional reference partition (see [compare_partitions()]) to score every run.
#' @param seed Random seed (community detection and bootstrap).
#' @param ... Other arguments of [wolf_packs()], fixed for all runs.
#' @return A list with `runs` (one row per run: settings, number of packs, individuals with a pack, ARI with the
#'   base run and, with `reference`, ARI with it), `individuals` (base pack, `stability`, share of runs with a pack)
#'   and `coassignment` (matrix).
#' @export
#' @examples
#' \donttest{
#' data("samples_Carpathians"); data("relate_Carpathians"); data("reference_Carpathians")
#' st <- pack_stability(samples_Carpathians, relate_Carpathians, "QG", 0.25,
#'                      vary = list(threshold = c(0.2, 0.25, 0.3)), n_boot = 10, reference = reference_Carpathians)
#' st$runs
#' }
pack_stability <- function(obs, relate, estimator, threshold, vary = list(), n_boot = 0, reference = NULL, seed = 1, ...) {
  base_args <- list(obs = obs, relate = relate, estimator = estimator, threshold = threshold, seed = seed, ...)
  base <- suppressMessages(do.call(wolf_packs, base_args))
  ids <- base$individuals$Individual
  base_p <- stats::setNames(base$individuals$Pack, ids)

  # Every run: its settings and its packs
  runs <- list()
  if (length(vary)) {
    grid <- expand.grid(vary, stringsAsFactors = FALSE)
    for (k in seq_len(nrow(grid))) {
      a <- base_args; a[names(grid)] <- as.list(grid[k, , drop = FALSE])
      r <- suppressMessages(do.call(wolf_packs, a))
      runs[[length(runs) + 1]] <- list(type = "sensitivity", setting = paste(names(grid), unlist(grid[k, ]), sep = " = ", collapse = ", "),
                                       packs = stats::setNames(r$individuals$Pack, r$individuals$Individual))
    }
  }
  if (n_boot > 0) {
    boot <- .with_seed(seed, lapply(seq_len(n_boot), function(b) {
      rows <- unlist(lapply(split(seq_len(nrow(obs)), as.character(obs$Individual)), function(r) r[sample.int(length(r), length(r), replace = TRUE)]))
      a <- base_args; a$obs <- obs[rows, ]; a$h <- base$ud$h                      # same bandwidth, only the samples change
      r <- suppressMessages(do.call(wolf_packs, a))
      stats::setNames(r$individuals$Pack, r$individuals$Individual)
    }))
    for (b in seq_len(n_boot)) runs[[length(runs) + 1]] <- list(type = "bootstrap", setting = paste("replicate", b), packs = boot[[b]])
  }
  if (!length(runs)) stop("Nothing to compare: give 'vary' and/or 'n_boot'.")

  # Scores of each run
  ari <- function(p, q) compare_partitions(p, q)$summary$ARI
  run_tab <- do.call(rbind, lapply(runs, function(r) data.frame(
    type = r$type, setting = r$setting, n_packs = length(unique(stats::na.omit(r$packs))),
    n_with_pack = sum(!is.na(r$packs)), ARI_base = ari(r$packs, base_p[!is.na(base_p)]),
    ARI_reference = if (is.null(reference)) NA_real_ else compare_partitions(r$packs, reference)$summary$ARI,
    stringsAsFactors = FALSE)))
  run_tab <- rbind(data.frame(type = "base", setting = "base", n_packs = length(unique(stats::na.omit(base_p))),
    n_with_pack = sum(!is.na(base_p)), ARI_base = 1,
    ARI_reference = if (is.null(reference)) NA_real_ else compare_partitions(base_p, reference)$summary$ARI), run_tab)

  # Who stays with whom
  mates <- function(p, i) if (is.na(p[i])) i else names(p)[p %in% p[i]]
  C <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  stab <- numeric(length(ids)); names(stab) <- ids
  with_pack <- numeric(length(ids)); names(with_pack) <- ids
  for (r in runs) {
    p <- r$packs[ids]; names(p) <- ids
    same <- outer(p, p, "==") & !is.na(outer(p, p, "=="))
    C <- C + same
    for (i in ids) {
      a <- mates(base_p, i); b <- mates(p, i)
      stab[i] <- stab[i] + length(intersect(a, b)) / length(union(a, b))
    }
    with_pack <- with_pack + !is.na(p)
  }
  n <- length(runs)
  individuals <- data.frame(Individual = ids, Pack = unname(base_p), stability = round(unname(stab) / n, 2),
                            prop_with_pack = round(unname(with_pack) / n, 2), stringsAsFactors = FALSE)
  list(runs = run_tab, individuals = individuals, coassignment = C / n)
}

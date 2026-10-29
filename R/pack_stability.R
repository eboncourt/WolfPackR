#' @title How robust are the packs?
#' @description Runs [wolf_packs()] many times and measures how much the packs change:
#'   * **sensitivity**: one run per combination of the settings given in `vary` (e.g. several bandwidths and
#'     relatedness thresholds), everything else at its base value;
#'   * **bootstrap**: `n_boot` runs where the samples of each individual are resampled with replacement
#'     (same number of samples per individual), to see the effect of sampling chance;
#'   * **subsampling** (rarefaction): for each share in `fraction`, `n_sub` runs keeping only that share of the
#'     samples, drawn at random. It mimics a lighter sampling effort: individuals with few samples drop out,
#'     territories are drawn from fewer points. Use it to check whether your sampling is dense enough: if the packs
#'     barely change when a quarter of the samples is removed, more sampling would probably not change them much
#'     either; if they fall apart, the sampling is close to the limit.
#'
#'   Bootstrap and subsampling runs keep the bandwidth of the base run (only the samples change).
#'   For each individual, `stability` is the mean, over the runs where it was sampled, of the Jaccard index between
#'   its pack mates in the base run and in that run (1 = always with the same companions). `coassignment` gives,
#'   for every pair, the share of runs where both are in the same pack.
#' @param obs,relate,estimator,threshold As in [wolf_packs()].
#' @param vary Named list of settings to vary, e.g. `list(h = c(1500, 2000, 2500), threshold = c(0.2, 0.25, 0.3))`.
#'   Any argument of [wolf_packs()] works (`h`, `threshold`, `min_share`, `merge`, `wide`, `method`...).
#' @param n_boot Number of bootstrap runs (default: 0).
#' @param fraction Shares of the samples to keep for the subsampling runs, e.g. `c(0.75, 0.5)` (default: none).
#' @param n_sub Number of runs per share in `fraction` (default: 20).
#' @param reference Optional reference partition (see [compare_partitions()]) to score every run.
#' @param seed Random seed (community detection, bootstrap and subsampling).
#' @param ... Other arguments of [wolf_packs()], fixed for all runs.
#' @return A list with:
#'   * `runs`: one row per run, with its type, settings, number of samples and of individuals sampled, number of
#'     packs, individuals with a pack, ARI with the base run and, with `reference`, ARI with it;
#'   * `subsampling` (with `fraction`): one row per share, with the median and the 2.5% and 97.5% quantiles of
#'     these scores over its runs;
#'   * `individuals`: base pack, `stability`, share of runs with a pack, number of runs where it was sampled;
#'   * `coassignment`: matrix, share of runs where two individuals are in the same pack.
#' @export
#' @examples
#' \donttest{
#' data("samples_Carpathians"); data("relate_Carpathians"); data("reference_Carpathians")
#' st <- pack_stability(samples_Carpathians, relate_Carpathians, "QG", 0.25,
#'                      vary = list(threshold = c(0.2, 0.25, 0.3)), n_boot = 10,
#'                      fraction = c(0.75, 0.5), n_sub = 10, reference = reference_Carpathians)
#' st$runs
#' st$subsampling
#' }
pack_stability <- function(obs, relate, estimator, threshold, vary = list(), n_boot = 0, fraction = NULL, n_sub = 20,
                           reference = NULL, seed = 1, ...) {
  if (!is.null(fraction) && any(fraction <= 0 | fraction > 1)) stop("'fraction' must be between 0 and 1.")
  base_args <- list(obs = obs, relate = relate, estimator = estimator, threshold = threshold, seed = seed, ...)
  base <- suppressMessages(do.call(wolf_packs, base_args))
  ids <- base$individuals$Individual
  base_p <- stats::setNames(base$individuals$Pack, ids)

  # Every run: its settings, its samples and its packs
  runs <- list()
  add_run <- function(type, setting, r, n_obs) runs[[length(runs) + 1]] <<- list(
    type = type, setting = setting, n_samples = n_obs, sampled = unique(as.character(r$samples$Individual)),
    packs = stats::setNames(r$individuals$Pack, r$individuals$Individual))
  if (length(vary)) {
    grid <- expand.grid(vary, stringsAsFactors = FALSE)
    for (k in seq_len(nrow(grid))) {
      a <- base_args; a[names(grid)] <- as.list(grid[k, , drop = FALSE])
      add_run("sensitivity", paste(names(grid), unlist(grid[k, ]), sep = " = ", collapse = ", "),
              suppressMessages(do.call(wolf_packs, a)), nrow(obs))
    }
  }
  # Same bandwidth as the base run for bootstrap and subsampling: only the samples change
  resampled <- function(rows) { a <- base_args; a$obs <- obs[rows, ]; a$h <- base$ud$h; suppressMessages(do.call(wolf_packs, a)) }
  if (n_boot > 0) {
    boot <- .with_seed(seed, lapply(seq_len(n_boot), function(b)
      resampled(unlist(lapply(split(seq_len(nrow(obs)), as.character(obs$Individual)), function(r) r[sample.int(length(r), length(r), replace = TRUE)])))))
    for (b in seq_len(n_boot)) add_run("bootstrap", paste("replicate", b), boot[[b]], nrow(obs))
  }
  for (f in fraction) {
    n_keep <- round(f * nrow(obs))
    sub <- .with_seed(seed, lapply(seq_len(n_sub), function(b) resampled(sort(sample.int(nrow(obs), n_keep)))))
    for (b in seq_len(n_sub)) add_run("subsampling", paste0("fraction = ", f), sub[[b]], n_keep)
  }
  if (!length(runs)) stop("Nothing to compare: give 'vary', 'n_boot' and/or 'fraction'.")

  # Scores of each run
  ari <- function(p, q) compare_partitions(p, q)$summary$ARI
  score <- function(type, setting, n_obs, n_ind, p) data.frame(
    type = type, setting = setting, n_samples = n_obs, n_individuals = n_ind, n_packs = length(unique(stats::na.omit(p))),
    n_with_pack = sum(!is.na(p)), ARI_base = if (type == "base") 1 else ari(p, base_p[!is.na(base_p)]),
    ARI_reference = if (is.null(reference)) NA_real_ else compare_partitions(p, reference)$summary$ARI,
    stringsAsFactors = FALSE)
  run_tab <- rbind(score("base", "base", nrow(obs), length(unique(as.character(obs$Individual))), base_p),
                   do.call(rbind, lapply(runs, function(r) score(r$type, r$setting, r$n_samples, length(r$sampled), r$packs))))

  # Subsampling: median and 95% range of the scores for each share of the samples
  subsampling <- NULL
  if (length(fraction)) {
    sr <- run_tab[run_tab$type == "subsampling", ]
    q <- function(x) stats::quantile(x, c(0.5, 0.025, 0.975), na.rm = TRUE, names = FALSE)
    subsampling <- do.call(rbind, lapply(fraction, function(f) {
      x <- sr[sr$setting == paste0("fraction = ", f), ]
      v <- c(n_samples = x$n_samples[1], n_individuals = stats::median(x$n_individuals), n_packs = stats::median(x$n_packs),
             stats::setNames(q(x$ARI_base), c("ARI_base", "ARI_base_low", "ARI_base_high")),
             if (!is.null(reference)) stats::setNames(q(x$ARI_reference), c("ARI_reference", "ARI_reference_low", "ARI_reference_high")))
      data.frame(fraction = f, t(round(v, 2)))
    }))
  }

  # Who stays with whom (an individual absent from a subsample is left out of that run)
  mates <- function(p, i) if (is.na(p[i])) i else names(p)[p %in% p[i]]
  C <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  N <- C                                                                     # runs where both were sampled
  stab <- with_pack <- n_runs <- stats::setNames(numeric(length(ids)), ids)
  for (r in runs) {
    p <- r$packs[ids]; names(p) <- ids
    here <- ids %in% r$sampled
    same <- outer(p, p, "==") & !is.na(outer(p, p, "=="))
    C <- C + same; N <- N + outer(here, here)
    for (i in ids[here]) {
      a <- mates(base_p, i); b <- mates(p, i)
      stab[i] <- stab[i] + length(intersect(a, b)) / length(union(a, b))
    }
    with_pack <- with_pack + (!is.na(p) & here); n_runs <- n_runs + here
  }
  individuals <- data.frame(Individual = ids, Pack = unname(base_p), stability = round(unname(stab / n_runs), 2),
                            prop_with_pack = round(unname(with_pack / n_runs), 2), n_runs = unname(n_runs),
                            stringsAsFactors = FALSE)
  out <- list(runs = run_tab, subsampling = subsampling, individuals = individuals, coassignment = ifelse(N > 0, C / N, NA))
  if (is.null(subsampling)) out$subsampling <- NULL
  out
}

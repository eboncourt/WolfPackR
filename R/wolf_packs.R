#' @title Whole workflow in one call
#' @description Runs the typical workflow: genetic groups ([genetic_groups()]) -> space use ([space_use()]) ->
#'   assignment of each individual ([assign_individuals()]) -> spatially coherent packs ([spatial_groups()]).
#'   Every step can also be run on its own with other settings.
#' @param obs An `sf` object of points (one row per sample) with an `Individual` column (and `Sex`, `Date`, period if any).
#' @param relate,estimator,threshold Pairwise relatedness table (e.g. from [relatedness()]), estimator column and threshold.
#' @param method Community detection for the genetic groups and the spatial split (default: `"leiden"`).
#' @param h Kernel bandwidth in metres (default: chosen by [space_use()]).
#' @param period Optional period column of `obs` (e.g. `"Campaign"`): spread measured within each period, dated dispersal.
#' @param min_share,min_samples,wide,wide_trim,outside,max_iter See [assign_individuals()].
#' @param merge,min_affinity,min_group_size See [spatial_groups()].
#' @param seed Random seed for the community detection.
#' @param ud Optional `wolfpack_ud` object, to skip the kernel step.
#' @return An object of class `wolfpack`, a list with:
#'   * `individuals`: one row per individual, the output of [assign_individuals()] plus `Pack` and `Core`;
#'   * `samples`: `obs` with `Genetic_Group`, `Group`, `Pack` and `Status` added, ready for [summarise_packs()],
#'     [plot_packs()] and [pack_timeline()];
#'   * `genetic`, `assignment`, `packs`: the outputs of each step;
#'   * `ud`: the space use object;
#'   * `settings`: the settings used (threshold, bandwidth and every other one), to report them or run again.
#' @export
#' @examples
#' data("samples_Carpathians"); data("relate_Carpathians")
#' res <- wolf_packs(samples_Carpathians, relate_Carpathians, "QG", 0.25, period = "Campaign")
#' res
#' head(res$individuals)
wolf_packs <- function(obs, relate, estimator, threshold, method = c("leiden", "louvain", "walktrap", "components"),
                       h = NULL, period = NULL, min_share = 0.5, min_samples = 3, wide = 3,
                       merge = 0.5, min_affinity = 0.1, min_group_size = 2, seed = 1, ud = NULL,
                       outside = 1, wide_trim = 0.1, max_iter = 20) {
  method <- match.arg(method)

  # 1. Genetic groups (individuals missing from the relatedness table have none)
  gg <- genetic_groups(relate, estimator, threshold, method = method, seed = seed)
  obs$Genetic_Group <- as.character(gg$group[match(as.character(obs$Individual), gg$Individual)])

  # 2. Space use, 3. assignment, 4. packs
  if (is.null(ud)) ud <- space_use(obs, h = h)
  ai <- assign_individuals(obs, "Genetic_Group", ud = ud, relate = relate, estimator = estimator, threshold = threshold,
                           min_share = min_share, min_samples = min_samples, wide = wide, period = period,
                           max_iter = max_iter, outside = outside, wide_trim = wide_trim)
  obs$Group <- ai$Group[match(as.character(obs$Individual), ai$Individual)]
  sg <- spatial_groups(obs, "Group", ud = ud, method = method, min_affinity = min_affinity, merge = merge,
                       min_samples = min_samples, min_group_size = min_group_size, wide = wide, period = period,
                       min_share = min_share, seed = seed, outside = outside, wide_trim = wide_trim)

  # Everything per individual, and back on the samples
  ind <- ai
  ind$Pack <- sg$Pack[match(ind$Individual, sg$Individual)]
  ind$Core <- sg$Core[match(ind$Individual, sg$Individual)]
  obs$Pack <- ind$Pack[match(as.character(obs$Individual), ind$Individual)]
  obs$Status <- ind$Status[match(as.character(obs$Individual), ind$Individual)]

  settings <- list(estimator = estimator, threshold = threshold, method = method, h = ud$h, period = period,
                   min_share = min_share, min_samples = min_samples, wide = wide, wide_trim = wide_trim,
                   outside = outside, max_iter = max_iter, merge = merge, min_affinity = min_affinity,
                   min_group_size = min_group_size, seed = seed)
  structure(list(individuals = ind, samples = obs, genetic = gg, assignment = ai, packs = sg, ud = ud,
                 settings = settings), class = "wolfpack")
}

#' @export
print.wolfpack <- function(x, ...) {
  ind <- x$individuals
  sizes <- sort(table(ind$Pack), decreasing = TRUE)
  cat(sprintf("WolfPackR: %d individuals, %d packs (%s), %d without a pack\n", nrow(ind), length(sizes),
              paste(names(sizes), sizes, sep = ": ", collapse = ", "), sum(is.na(ind$Pack))))
  cat("Statuses:", .split_text(ind$Status), "\n")
  cat(sprintf("Kernel bandwidth h = %.0f m\n", x$ud$h))
  invisible(x)
}

#' @title Packs through time
#' @description Follows individuals and territories from one sampling period to the next (e.g. biological years),
#'   with territories estimated on all periods together (more samples, steadier outlines).
#'
#'   * **Individuals**: each sample is attributed to a territory ([attribute_samples()]), then each individual gets
#'     its majority territory in each period (`mixed` when no territory holds the majority). Leaving the mixed
#'     periods aside, the sequence tells whether it stayed (`stable`), moved for good to another pack
#'     (`changed pack`, with `change_period`), left the territories (`left`), went back and forth (`back and forth`),
#'     was never in a territory, or was seen in a single period (`single period`). Its spread (radius of gyration) in each
#'     period separates residents patrolling a big territory from individuals moving between places.
#'   * **Territories**: for each pack and period, the members detected in the territory (a member is detected in a
#'     period when one of its samples of that period falls in its pack's territory), the new ones and the ones
#'     missed (detected before and after, but not in that period). When the members have at least `min_territory`
#'     samples in the period (leaving out those then living in another pack's territory, e.g. before dispersal),
#'     the territory of that period is compared with the territory over all periods (Bhattacharyya's affinity and
#'     distance between the two centres) and with the territories of the other packs in the same period (share of
#'     its `level` isopleth inside theirs). Low affinities or large shifts point to a territory that moved or
#'     split: the analysis is then worth running period by period.
#'
#'   All this describes what was sampled: a pack or a member not detected in a period may simply not have been
#'   sampled there, and territories drawn from few samples are smaller, so they overlap less. Read the results
#'   in the light of the sampling effort of each period.
#' @param obs An `sf` object of points (one row per sample) with `Individual`, the pack column and the period column.
#' @param pack Name of the pack column (default `"Pack"`, `NA` = no pack).
#' @param period Name of the period column (default `"Campaign"`).
#' @param ud A `wolfpack_ud` object from [space_use()]; computed from `obs` if `NULL`.
#' @param h Bandwidth passed to [space_use()] when `ud` is `NULL`.
#' @param min_share Share of an individual's samples in a period needed to give it a territory for that period (default: 0.5).
#' @param min_territory Minimum number of samples of a pack in a period to draw its territory for that period (default: 5).
#' @param level Isopleth level for the overlap between packs in a period (default: 0.95).
#' @return A list with:
#'   * `individuals`: one row per individual, with its territory in each period, `pattern`, `change_period`
#'     and its radius of gyration in each period (km);
#'   * `samples`: the territory of each sample;
#'   * `territories`: one row per pack and period, with `n_samples` (samples of the members, leaving out the members
#'     that were then in another pack's territory), `n_detected`,
#'     `n_new` (`NA` in the first period), `n_missed`, `affinity_with_all_periods`, `shift_km` (distance between
#'     the centre of the territory in that period and over all periods) and `shared_pct` (share of the territory
#'     inside the territory of another pack in the same period; `NA` when no other pack has enough samples).
#' @references Kanishka, Koley, S., Khan, S., Hussain, Z., Sathyakumar, S., & Habib, B. (2026). The missing metric
#'   in ecological studies: Land tenure system of territorial carnivores. Journal of Animal Ecology.
#'   \doi{10.1111/1365-2656.70373}
#' @seealso [summarise_packs()] with `period`, which sums this up for each pack.
#' @export
#' @examples
#' data("samples_Carpathians"); data("relate_Carpathians")
#' res <- wolf_packs(samples_Carpathians, relate_Carpathians, "QG", 0.25, period = "Campaign")
#' tl <- pack_timeline(res$samples, ud = res$ud)
#' table(tl$individuals$pattern)
#' tl$territories
pack_timeline <- function(obs, pack = "Pack", period = "Campaign", ud = NULL, h = NULL,
                          min_share = 0.5, min_territory = 5, level = 0.95) {
  pack <- .col_name(substitute(pack), function() pack, names(obs))
  for (col in c(pack, period)) if (!col %in% names(obs)) stop("Column '", col, "' not found in 'obs'.")
  if (is.null(ud)) ud <- space_use(obs, h = h)
  ids <- ud$ids
  pk <- tapply(as.character(obs[[pack]]), as.character(obs$Individual), function(x) x[1])[ids]
  pk[pk %in% "Lone Individual"] <- NA
  terr <- split(ids[!is.na(pk)], pk[!is.na(pk)])

  # Territory of every sample, then of every individual in every period
  s <- attribute_samples(ud, terr)
  per <- sort(unique(stats::na.omit(s[[period]])))
  seq_mat <- t(vapply(ids, function(i) vapply(per, function(p) {
    t <- s$Territory[s$Individual == i & s[[period]] %in% p]
    if (!length(t)) return(NA_character_)
    m <- .majority(t, min_share); if (is.na(m)) "mixed" else m
  }, character(1)), character(length(per))))
  if (length(per) == 1) seq_mat <- matrix(seq_mat, ncol = 1)
  dimnames(seq_mat) <- list(ids, per)

  # What the sequence says ("mixed" periods, without a majority, are left out of the reading)
  classify <- function(L) {
    L <- L[!is.na(L)]
    if (length(L) <= 1) return(c("single period", NA))
    L <- L[L != "mixed"]
    if (!length(L)) return(c("unclear", NA))
    if (all(L == L[1])) return(c(if (L[1] %in% names(terr)) "stable" else "never in a territory", NA))
    first <- L[1]; last <- L[length(L)]
    back <- which(L == first)
    changed <- max(back) < length(L) && all(L[(max(back) + 1):length(L)] != first)
    if (!changed) return(c("back and forth", NA))
    change_at <- names(L)[max(back) + 1]
    if (last %in% names(terr) && all(L[(max(back) + 1):length(L)] == last)) c("changed pack", change_at) else c("left", change_at)
  }
  cls <- t(vapply(ids, function(i) { L <- seq_mat[i, ]; names(L) <- per; classify(L) }, character(2)))

  # Spread of each individual in each period (km)
  rg <- t(vapply(ids, function(i) vapply(per, function(p) {
    w <- s$Individual == i & s[[period]] %in% p
    if (sum(w) < 2) NA_real_ else round(.rg(as.matrix(s[w, c("x", "y")])) / 1000, 1)
  }, numeric(1)), numeric(length(per))))
  if (length(per) == 1) rg <- matrix(rg, ncol = 1)

  individuals <- data.frame(Individual = ids, Pack = unname(pk), seq_mat, pattern = cls[, 1], change_period = cls[, 2],
                            stringsAsFactors = FALSE, check.names = FALSE)
  rg_df <- as.data.frame(rg); names(rg_df) <- paste0("rg_km_", per)
  individuals <- cbind(individuals, rg_df)
  rownames(individuals) <- NULL

  # Territory of each pack in each period: who was detected in it, how much it moved, how much it was shared
  # (a member is detected in a period when one of its samples of that period falls in its pack's territory)
  s$Pack <- unname(pk[as.character(s$Individual)])
  home <- s[!is.na(s$Pack) & s$Territory == s$Pack, ]
  # Samples drawing the territory of pack k in period p: its members', except those then living in another pack's territory
  in_period <- function(k, p) {
    away <- ids[seq_mat[, p] %in% setdiff(names(terr), k)]
    s$Individual %in% setdiff(terr[[k]], away) & s[[period]] %in% p
  }
  per_ud <- lapply(stats::setNames(names(terr), names(terr)), function(k) lapply(per, function(p) {
    w <- in_period(k, p)
    if (sum(w) >= min_territory) .kde(as.matrix(s[w, c("x", "y")]), ud$grid, ud$h) else NULL
  }))
  centre <- function(u) colSums(ud$grid * u)
  territories <- do.call(rbind, lapply(names(terr), function(k) do.call(rbind, lapply(seq_along(per), function(b) {
    p <- per[b]
    seen <- unique(home$Individual[home$Pack == k & home[[period]] %in% p])
    before <- unique(home$Individual[home$Pack == k & home[[period]] %in% per[seq_len(b - 1)]])
    after <- unique(home$Individual[home$Pack == k & home[[period]] %in% per[-seq_len(b)]])
    u <- per_ud[[k]][[b]]; all_u <- .mix_ud(ud, terr[[k]])
    others <- Filter(Negate(is.null), lapply(per_ud[setdiff(names(terr), k)], `[[`, b))   # neighbours with enough samples that period
    shared <- if (is.null(u) || !length(others)) NA_real_ else {
      mine <- .hr_cells(u, level)
      100 * mean(Reduce(`|`, lapply(others, .hr_cells, level))[mine])
    }
    data.frame(Pack = k, period = p, n_samples = sum(in_period(k, p)),
               n_detected = length(seen),
               n_new = if (b == 1) NA_integer_ else length(setdiff(seen, before)),     # first detected in the territory
               n_missed = length(setdiff(intersect(before, after), seen)),             # seen before and after, not now
               affinity_with_all_periods = if (is.null(u)) NA_real_ else round(sum(sqrt(u * all_u)), 2),
               shift_km = if (is.null(u)) NA_real_ else round(sqrt(sum((centre(u) - centre(all_u))^2)) / 1000, 1),
               shared_pct = round(shared, 1), stringsAsFactors = FALSE)
  }))))
  rownames(territories) <- NULL

  list(individuals = individuals, samples = s, territories = territories)
}

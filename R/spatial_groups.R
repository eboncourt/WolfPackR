#' @title Split groups into spatially coherent packs
#' @description Within each group (e.g. the groups from [assign_individuals()]), finds the sets of individuals that
#'   share the same space, so that a family spread over two neighbouring territories gives two packs.
#'
#'   1. **Core individuals** (at least `min_samples` samples and not wide-ranging) are linked by the overlap of their
#'      UDs (Bhattacharyya's affinity, [ud_overlap()]); links below `min_affinity` are dropped.
#'   2. **Communities** are detected on this weighted graph (Leiden by default), which is much less prone than
#'      connected components to chain neighbouring territories through a single individual.
#'   3. **Merge**: two subgroups whose UDs overlap more than `merge` (affinity) are merged back, so a territory
#'      is not cut in two just because of its internal structure.
#'   4. **The other members** (few samples, wide-ranging) join the subgroup where most of their samples fall
#'      ([attribute_samples()]), or none.
#'   5. Subgroups with fewer than `min_group_size` individuals are dropped (a pair is enough by default).
#'
#'   Individuals without a group (`NA`, e.g. Transient or Undetermined) are pooled and go through the same steps:
#'   that's how a territory held by unrelated individuals (e.g. a new pair whose pups were not sampled) can show up.
#' @param obs An `sf` object of points (one row per sample) with an `Individual` column.
#' @param group Name of the column giving the group of each individual (`NA` = no group). `NULL`: everybody in a
#'   single group, i.e. a purely spatial analysis.
#' @param ud A `wolfpack_ud` object from [space_use()]; computed from `obs` if `NULL`.
#' @param h Bandwidth passed to [space_use()] when `ud` is `NULL`.
#' @param method Community detection: `"leiden"` (default), `"louvain"`, `"walktrap"` or `"components"`.
#' @param min_affinity Overlaps below this value are not links (default: 0.1).
#' @param merge Subgroups overlapping more than this are merged (default: 0.5).
#' @param min_samples Minimum number of samples of a core individual (default: 3).
#' @param min_group_size Minimum number of individuals in a pack (default: 2).
#' @param wide,period Wide-ranging individuals (radius of gyration above `wide` x h, measured within each `period`
#'   if given) are not core individuals (default `wide`: 3).
#' @param min_share Share of the samples needed for a non-core individual to join a pack (default: 0.5).
#' @param seed Random seed for Leiden and Louvain (your own random numbers are left alone).
#' @return A data.frame with one row per individual: `Individual`, `Group`, `Pack` (`NA` = no pack) and `Core`.
#' @export
#' @examples
#' data("samples_Carpathians")
#' # purely spatial analysis, no genetics
#' sg <- spatial_groups(samples_Carpathians)
#' table(sg$Pack, useNA = "ifany")
spatial_groups <- function(obs, group = NULL, ud = NULL, h = NULL,
                           method = c("leiden", "louvain", "walktrap", "components"),
                           min_affinity = 0.1, merge = 0.5, min_samples = 3, min_group_size = 2,
                           wide = 3, period = NULL, min_share = 0.5, seed = 1) {
  method <- match.arg(method)
  if (!is.null(group)) {
    group <- .col_name(substitute(group), function() group, names(obs))
    if (!group %in% names(obs)) stop("Column '", group, "' not found in 'obs'.")
  }
  if (is.null(ud)) ud <- space_use(obs, h = h)
  ids <- ud$ids
  grp <- if (is.null(group)) stats::setNames(rep("", length(ids)), ids) else
    tapply(as.character(obs[[group]]), as.character(obs$Individual), function(x) x[1])[ids]
  grp_key <- ifelse(is.na(grp), "new", grp)                                # individuals without a group are pooled
  is_wide <- .wide_ranging(ud, wide, period)
  core_ok <- ud$n[ids] >= min_samples & !is_wide[ids]
  BA <- ud_overlap(ud, "BA")

  # UD affinity between two sets of individuals
  set_affinity <- function(a, b) sum(sqrt(.mix_ud(ud, a) * .mix_ud(ud, b)))

  pack <- stats::setNames(rep(NA_character_, length(ids)), ids)
  for (g in unique(grp_key)) {
    members <- ids[grp_key == g]
    core <- members[core_ok[members]]

    # Communities of core individuals on the overlap graph
    subs <- if (length(core) >= 3) {
      W <- BA[core, core]; W[W < min_affinity] <- 0; diag(W) <- 0
      gr <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)
      memb <- .with_seed(seed, .communities(gr, method))
      split(core, memb)
    } else if (length(core) == 2 && BA[core[1], core[2]] < min_affinity) as.list(core)   # two cores that don't overlap stay apart
    else if (length(core)) list(core) else list()

    # Merge back the subgroups that overlap a lot, the most overlapping pair first
    while (length(subs) >= 2) {
      pr <- utils::combn(length(subs), 2)
      aff <- apply(pr, 2, function(p) set_affinity(subs[[p[1]]], subs[[p[2]]]))
      if (max(aff) <= merge) break
      p <- pr[, which.max(aff)]
      subs[[p[1]]] <- c(subs[[p[1]]], subs[[p[2]]]); subs[[p[2]]] <- NULL
    }
    if (!length(subs)) next
    names(subs) <- paste0("s", seq_along(subs))

    # The other members join the subgroup where most of their samples fall
    others <- setdiff(members, unlist(subs))
    if (length(others)) {
      att <- attribute_samples(ud, subs)
      for (i in others) {
        m <- .majority(att$Territory[att$Individual == i], min_share)
        if (!is.na(m) && m != "outside") subs[[m]] <- c(subs[[m]], i)
      }
    }

    # Keep the subgroups big enough, and name them
    subs <- subs[lengths(subs) >= min_group_size]
    if (!length(subs)) next
    lab <- if (length(subs) == 1) g else paste0(g, if (g != "") ".", seq_along(subs))
    if (g == "") lab <- if (length(subs) == 1) "1" else as.character(seq_along(subs))
    for (k in seq_along(subs)) pack[subs[[k]]] <- lab[k]
  }

  out <- data.frame(Individual = ids, Group = unname(if (is.null(group)) NA_character_ else grp),
                    Pack = unname(pack), Core = unname(core_ok[ids]), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' @title Assign each individual to the territory where its samples fall
#' @description Compares where each individual's samples fall with the territories of the genetic groups,
#'   and gives every individual a group and a status. Replaces `ugly_ducklings()`.
#'
#'   1. **Reference territories**: the UD of each genetic group is the mix of its members' UDs. Wide-ranging
#'      individuals (radius of gyration above `wide` x h, i.e. wider than a resident covering a whole territory)
#'      don't help draw them, and each individual is judged against territories drawn without itself. The radius
#'      of gyration leaves out the farthest 10% of the samples (`wide_trim`; at least one, keeping at least three),
#'      so that a single foray into a neighbouring territory doesn't make a resident look wide-ranging.
#'   2. **Sample attribution** ([attribute_samples()]): each sample goes to the territory where it is most
#'      likely, or "outside" when no territory beats a uniform use of the area.
#'   3. **Majority rule**: the individual goes where more than `min_share` of its samples are (a tie never wins).
#'   4. Groups are updated (an individual can move to another group, or be left out) and steps 1-3 are
#'      repeated until nobody moves. Individuals that keep going back and forth between two groups (the
#'      iterations cycle) are set Undetermined and listed in `attr(, "flipping")`.
#'
#'   Individuals without any relative above the threshold (genetic singletons) start without a group and are
#'   assigned like the others, by where their samples fall.
#'
#' @section Statuses:
#' * `Resident`: assigned to a group where it has at least one relative;
#' * `Putative disperser`: assigned to a group where it has no relative, while it has relatives elsewhere
#'   (its natal group is probably another one);
#' * `Disperser`: a putative disperser whose samples, with `period`, move for good from the territory of its
#'   relatives to its new one (`change_period` gives when);
#' * `Unrelated resident`: assigned to a group, no relative detected anywhere (immigrant, breeder from
#'   outside the study area, or relatives not sampled);
#' * `Transient`: most samples outside every territory (floater, disperser on the move);
#' * `Undetermined`: no majority, or fewer than `min_samples` samples for any status other than Resident.
#'
#' Relatives are pairs above `threshold` in `relate` when it is given; otherwise members of the same genetic group.
#' @param obs An `sf` object of points (one row per sample) with `Individual` and the genetic group column.
#' @param group Name of the genetic group column (e.g. from [genetic_groups()]). `NA` = no group.
#' @param ud A `wolfpack_ud` object from [space_use()]; computed from `obs` if `NULL`.
#' @param h Bandwidth passed to [space_use()] when `ud` is `NULL`.
#' @param relate,estimator,threshold Pairwise relatedness table, estimator column and threshold, to define
#'   relatives pair by pair (recommended). If `NULL`, relatives are the members of the same genetic group.
#' @param min_share Share of the samples needed to assign an individual (default: 0.5).
#' @param min_samples Below this number of samples, any status other than Resident becomes Undetermined (default: 3).
#' @param wide Individuals whose radius of gyration (leaving out the farthest `wide_trim` of their samples) exceeds
#'   `wide` x h don't help draw the territories (default: 3, about 5.6 km in the Carpathians, where a resident using
#'   a whole 180 km2 territory has a radius of gyration around 5.4 km).
#' @param wide_trim Share of the samples farthest from the centre left out when measuring the radius of gyration
#'   (rounded up, keeping at least 3; default: 0.1). 0 gives the plain radius of gyration (versions before 1.1.0).
#' @param period Optional column of `obs` with the sampling period (e.g. `"Campaign"`): the spread of each individual
#'   is then measured within each period, and dispersal can be dated.
#' @param max_iter Maximum number of iterations (default: 20).
#' @param outside A sample is "outside" every territory when its best territory density is below `outside` times a
#'   uniform use of the study area (default: 1, see [attribute_samples()]). Higher values send more samples
#'   outside, hence more Transient individuals.
#' @return A data.frame with one row per individual: `Individual`, `Genetic_Group`, `Group` (assigned group, `NA`
#'   for Transient and Undetermined), `Status`, `n_samples`, `share` (share of its samples in `Group`),
#'   `split` (where its samples went), `wide_ranging`, `low_data` and `change_period`. The attribution of every
#'   sample is in `attr(, "samples")`.
#' @export
#' @examples
#' data("samples_Carpathians"); data("relate_Carpathians")
#' gg <- genetic_groups(relate_Carpathians, "QG", 0.25, method = "leiden")
#' obs <- merge(samples_Carpathians, gg, by = "Individual")
#' ai <- assign_individuals(obs, "group", relate = relate_Carpathians, estimator = "QG", threshold = 0.25,
#'                          period = "Campaign")
#' table(ai$Status)
assign_individuals <- function(obs, group, ud = NULL, h = NULL, relate = NULL, estimator = NULL, threshold = NULL,
                               min_share = 0.5, min_samples = 3, wide = 3, period = NULL, max_iter = 20,
                               outside = 1, wide_trim = 0.1) {

  # Inputs: one genetic group per individual, the UDs, who is wide-ranging
  group <- .col_name(substitute(group), function() group, names(obs))
  if (!group %in% names(obs)) stop("Column '", group, "' not found in 'obs'.")
  if (!is.null(period) && !period %in% names(obs)) stop("Column '", period, "' not found in 'obs'.")
  if (is.null(ud)) ud <- space_use(obs, h = h)
  ids <- ud$ids
  gen <- tapply(as.character(obs[[group]]), as.character(obs$Individual), function(x) x[1])[ids]
  gsize <- table(gen)
  gen_start <- ifelse(!is.na(gen) & gsize[gen] >= 2, gen, NA)             # singletons start without a group
  names(gen_start) <- ids
  is_wide <- .wide_ranging(ud, wide, period, wide_trim)

  # Iterate: territories from the current groups, attribute samples, move individuals, until stable
  cur <- gen_start
  history <- list()
  flipping <- character(0)                                                # individuals going back and forth between groups
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    ok <- !is.na(cur) & !cur %in% c("outside", "undetermined")
    terr <- split(ids[ok & !is_wide], cur[ok & !is_wide])
    att <- attribute_samples(ud, terr, outside = outside)
    new <- vapply(ids, function(i) {
      m <- .majority(att$Territory[att$Individual == i], min_share)
      if (is.na(m)) "undetermined" else m
    }, character(1))
    new[flipping] <- "undetermined"
    if (identical(unname(new), unname(cur))) { converged <- TRUE; break }

    # Going round in circles (back to an earlier state): whoever keeps moving is left undetermined
    seen <- which(vapply(history, function(h) identical(unname(h), unname(new)), logical(1)))
    if (length(seen)) {
      cycle <- do.call(cbind, c(history[min(seen):length(history)], list(cur, new)))
      flipping <- union(flipping, ids[apply(cycle, 1, function(x) length(unique(x)) > 1)])
      new[flipping] <- "undetermined"
    }
    history[[length(history) + 1]] <- cur
    cur <- new
  }
  names(cur) <- ids
  if (!converged) warning("No convergence after ", max_iter, " iterations.")

  # Relatives: pair by pair when we have the relatedness table, otherwise same genetic group
  K <- if (!is.null(relate)) .kin_matrix(relate, estimator, threshold, ids) else
    outer(gen, gen, function(a, b) !is.na(a) & !is.na(b) & a == b & gsize[a] >= 2) & !diag(length(ids))
  dimnames(K) <- list(ids, ids)
  assigned <- stats::setNames(!cur %in% c("outside", "undetermined"), ids)

  status <- vapply(ids, function(i) {
    if (cur[i] == "outside") return("Transient")
    if (cur[i] == "undetermined") return("Undetermined")
    here <- setdiff(ids[assigned & cur == cur[i]], i)
    if (any(K[i, here])) return("Resident")
    if (any(K[i, setdiff(ids, c(here, i))])) return("Putative disperser")
    "Unrelated resident"
  }, character(1))

  # Dated dispersal: relatives' territory first, then the new territory for good
  samples <- att
  change <- stats::setNames(rep(NA_character_, length(ids)), ids)
  if (!is.null(period)) {
    for (i in ids[status == "Putative disperser"]) {
      natal <- unique(cur[assigned & K[i, ]])                               # territories of its relatives
      si <- samples[samples$Individual == i, ]
      per <- sort(unique(si[[period]]))
      maj <- vapply(per, function(p) { m <- .majority(si$Territory[si[[period]] == p], min_share); if (is.na(m)) "mixed" else m }, character(1))
      last_natal <- max(c(0, which(maj %in% natal)))
      after <- if (last_natal < length(maj)) maj[(last_natal + 1):length(maj)] else character(0)
      if (last_natal > 0 && length(after) && all(after == cur[i])) {
        status[i] <- "Disperser"
        change[i] <- per[last_natal + 1]
      }
    }
  }

  # Not enough samples to say anything but Resident
  n <- ud$n[ids]
  low <- n < min_samples
  status[low & status != "Resident"] <- "Undetermined"

  share <- vapply(ids, function(i) { t <- samples$Territory[samples$Individual == i]; if (assigned[i]) mean(t == cur[i]) else NA_real_ }, numeric(1))
  out <- data.frame(
    Individual = ids, Genetic_Group = unname(gen), Group = unname(ifelse(assigned & status != "Undetermined", cur, NA)),
    Status = unname(status), n_samples = unname(n), share = round(unname(share), 2),
    split = vapply(ids, function(i) .split_text(samples$Territory[samples$Individual == i]), character(1)),
    wide_ranging = unname(is_wide), low_data = unname(low), change_period = unname(change),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "samples") <- samples
  attr(out, "iterations") <- it
  attr(out, "flipping") <- flipping
  out
}

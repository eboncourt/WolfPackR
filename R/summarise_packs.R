#' @title Summarise packs
#' @description One line per pack, along the lines of the land tenure framework of Kanishka et al. (2026):
#'   * **Composition**: number of individuals (by sex), of samples, and statuses;
#'   * **Space use**: territory and core area (kernel isopleths), share of the territory also inside the territory
#'     of another pack, and nearest pack (distance between the centres of their UDs);
#'   * **Through time** (with `period` and/or `date`): when the pack was detected in its territory, how many members
#'     were detected there in more than one period, how much the territory changed from one period to the next
#'     (summed up from [pack_timeline()]) and the time from the first to the last detection in the territory;
#'   * **Kinship** (with `relate`): mean relatedness within the pack, tested against random packs of the same sizes
#'     (permutations of the pack labels among the individuals in a pack).
#'
#'   A member is detected in a period when one of its samples of that period falls in its pack's territory.
#'   The time columns describe what was sampled: they are minimum values (a pack was there at least from its first to
#'   its last detection), a pack or a member not detected in a period may simply not have been sampled there, and
#'   a study as long as the tenures can't measure them. Read them in the light of the sampling effort.
#'
#'   Earlier versions also gave "putative dominant" individuals (the most sampled male and female). This is gone:
#'   on the Carpathian data, the most sampled individual of each sex was the known breeder in only half of the
#'   cases. The number of samples per individual is still there, as a plain description, in `attr(, "individuals")`.
#' @param obs An `sf` object of points (one row per sample) with `Individual` and the pack column; `Sex` and
#'   `Status` columns are used when present.
#' @param pack Name of the pack column (default `"Pack"`). `NA` or `"Lone Individual"` = no pack.
#' @param ud A `wolfpack_ud` object from [space_use()]; computed from `obs` if `NULL`.
#' @param h Bandwidth passed to [space_use()] when `ud` is `NULL`.
#' @param levels Isopleth levels for the areas (default: 95% territory and 50% core area). The largest one is used
#'   for the space shared between packs.
#' @param relate,estimator Optional pairwise relatedness table and estimator column, for the within-pack relatedness.
#' @param n_perm Number of permutations for the relatedness test (default: 999).
#' @param period Optional name of the column of `obs` with the sampling period (e.g. `"Campaign"`).
#' @param date Optional name of the column of `obs` with the sampling date (class `Date`).
#' @param min_territory Minimum number of samples of a pack in a period to draw its territory for that period
#'   (default: 5), see [pack_timeline()].
#' @param sex_column,male_pattern,female_pattern Column holding the sex and regex patterns for males and females.
#' @param seed Random seed for the permutations (your own random numbers are left alone).
#' @return A data.frame with one row per pack:
#'   * `Pack`, `n_individuals`, `n_males`, `n_females`, `n_samples`, one `area<level>_km2` column per level and
#'     `statuses` (if `obs` has a `Status` column);
#'   * `shared_pct` (% of the territory inside the territory of at least one other pack), `nearest_pack`, `nearest_km`;
#'   * with `period`: `first_period`, `last_period` and `n_periods` (periods in which the pack was detected in its
#'     territory), `n_persistent` (members detected there in 2 periods or more), `min_affinity` (affinity between the
#'     territory of the least similar period and the territory over all periods) and `max_shift_km` (largest
#'     distance between the centre of the territory in a period and over all periods);
#'   * with `date`: `tenure_years`, time from the first to the last detection of a member in the territory;
#'   * with `relate`: `mean_r`, `r_z` and `r_p` (permutation p-value, one-sided: more related than random packs).
#'
#'   Per-individual details are in `attr(, "individuals")` (with `period` and `date`: first and last period in
#'   the territory, number of periods, pattern from [pack_timeline()], tenure in years), and, with `period`, the
#'   pack x period table of [pack_timeline()] is in `attr(, "periods")`.
#' @references Kanishka, Koley, S., Khan, S., Hussain, Z., Sathyakumar, S., & Habib, B. (2026). The missing metric
#'   in ecological studies: Land tenure system of territorial carnivores. Journal of Animal Ecology.
#'   \doi{10.1111/1365-2656.70373}
#' @export
#' @examples
#' data("samples_Carpathians"); data("relate_Carpathians")
#' res <- wolf_packs(samples_Carpathians, relate_Carpathians, "QG", 0.25, period = "Campaign")
#' sm <- summarise_packs(res$samples, "Pack", ud = res$ud, relate = relate_Carpathians, estimator = "QG",
#'                       n_perm = 199, period = "Campaign", date = "Date")
#' sm
#' attr(sm, "periods")
summarise_packs <- function(obs, pack = "Pack", ud = NULL, h = NULL, levels = c(0.95, 0.5),
                            relate = NULL, estimator = NULL, n_perm = 999,
                            period = NULL, date = NULL, min_territory = 5,
                            sex_column = "Sex", male_pattern = "^M", female_pattern = "^F", seed = 1) {
  pack <- .col_name(substitute(pack), function() pack, names(obs))
  if (!pack %in% names(obs)) stop("Column '", pack, "' doesn't exist.")
  for (col in c(period, date)) if (!col %in% names(obs)) stop("Column '", col, "' doesn't exist.")
  if (!is.null(date) && !inherits(obs[[date]], "Date")) stop("Column '", date, "' must be of class Date (see as.Date()).")
  if (is.null(ud)) ud <- space_use(obs, h = h)

  # One row per individual
  d <- sf::st_drop_geometry(obs)
  d$Individual <- as.character(d$Individual)
  first <- function(col) tapply(as.character(d[[col]]), d$Individual, function(x) x[1])
  ids <- sort(unique(d$Individual))
  ind <- data.frame(Individual = ids, Pack = first(pack)[ids], n_samples = as.integer(table(d$Individual)[ids]),
                    stringsAsFactors = FALSE)
  ind$Pack[ind$Pack %in% "Lone Individual"] <- NA
  if (sex_column %in% names(d)) ind$Sex <- first(sex_column)[ids]
  if ("Status" %in% names(d)) ind$Status <- first("Status")[ids]
  rownames(ind) <- NULL
  packs <- sort(unique(stats::na.omit(ind$Pack)))
  terr <- split(ind$Individual[!is.na(ind$Pack)], ind$Pack[!is.na(ind$Pack)])[packs]

  # Composition and areas
  out <- do.call(rbind, lapply(packs, function(p) {
    m <- ind[ind$Pack %in% p, ]
    row <- data.frame(Pack = p, n_individuals = nrow(m),
                      n_males = if (!is.null(m$Sex)) sum(grepl(male_pattern, m$Sex)) else NA_integer_,
                      n_females = if (!is.null(m$Sex)) sum(grepl(female_pattern, m$Sex)) else NA_integer_,
                      n_samples = sum(m$n_samples), stringsAsFactors = FALSE)
    iso <- isopleth(ud, m$Individual, levels)
    for (k in seq_along(levels)) row[[paste0("area", round(100 * levels[k]), "_km2")]] <- round(iso$area_km2[k], 1)
    if (!is.null(m$Status)) row$statuses <- .split_text(m$Status)
    row
  }))
  if (is.null(out)) out <- data.frame(Pack = character(0), stringsAsFactors = FALSE)       # no pack at all

  # Space shared with the other packs (largest isopleth level) and nearest pack (distance between UD centres)
  if (length(packs)) {
    H <- matrix(vapply(packs, function(p) .hr_cells(.mix_ud(ud, terr[[p]]), max(levels)), logical(nrow(ud$grid))),
                ncol = length(packs), dimnames = list(NULL, packs))         # cells of each territory
    out$shared_pct <- round(vapply(packs, function(p)
      if (length(packs) < 2) 0 else 100 * mean(rowSums(H[H[, p], setdiff(packs, p), drop = FALSE]) > 0), numeric(1)), 1)
    ctr <- t(vapply(packs, function(p) colSums(ud$grid * .mix_ud(ud, terr[[p]])), numeric(2)))
    D <- as.matrix(stats::dist(ctr)) / 1000; diag(D) <- Inf
    out$nearest_pack <- if (length(packs) < 2) NA_character_ else packs[apply(D, 1, which.min)]
    out$nearest_km <- if (length(packs) < 2) NA_real_ else round(apply(D, 1, min), 1)
  }

  # Through time: when each pack and each member was detected in the pack's territory
  tl <- NULL
  if (length(packs) && (!is.null(period) || !is.null(date))) {
    if (!is.null(period)) {
      tl <- pack_timeline(obs, pack, period, ud = ud, min_territory = min_territory, level = max(levels))
      s <- tl$samples
    } else s <- attribute_samples(ud, terr)
    s$Individual <- as.character(s$Individual)
    s$Pack <- ind$Pack[match(s$Individual, ind$Individual)]
    home <- s[!is.na(s$Pack) & s$Territory == s$Pack, ]                     # members' samples inside their own territory
    at_home <- function(i, col) home[[col]][home$Individual %in% i]

    if (!is.null(period)) {
      per <- sort(unique(stats::na.omit(s[[period]])))
      span <- function(x) {                                                 # first period, last period, number of periods
        x <- match(x, per)
        if (!length(x)) c(NA, NA, 0) else c(per[min(x)], per[max(x)], length(unique(x)))
      }
      sp <- t(vapply(ind$Individual, function(i) as.character(span(at_home(i, period))), character(3)))
      ind$first_period <- sp[, 1]; ind$last_period <- sp[, 2]; ind$n_periods <- as.integer(sp[, 3])
      ind$n_periods[is.na(ind$Pack)] <- NA
      ind$pattern <- tl$individuals$pattern[match(ind$Individual, tl$individuals$Individual)]

      sp <- t(vapply(packs, function(p) as.character(span(at_home(terr[[p]], period))), character(3)))
      out$first_period <- sp[, 1]; out$last_period <- sp[, 2]; out$n_periods <- as.integer(sp[, 3])
      out$n_persistent <- vapply(packs, function(p) sum(ind$n_periods[ind$Pack %in% p] >= 2), integer(1))
      tt <- tl$territories                                                  # pack x period, from pack_timeline()
      worst <- function(col, f) vapply(packs, function(p) {
        v <- tt[[col]][tt$Pack == p & !is.na(tt[[col]])]
        if (length(v)) f(v) else NA_real_
      }, numeric(1))
      out$min_affinity <- worst("affinity_with_all_periods", min)
      out$max_shift_km <- worst("shift_km", max)
    }
    if (!is.null(date)) {
      years <- function(x) { x <- x[!is.na(x)]; if (length(x) < 2) 0 else round(as.numeric(diff(range(x))) / 365.25, 1) }
      ind$tenure_years <- vapply(ind$Individual, function(i) years(at_home(i, date)), numeric(1))
      ind$tenure_years[is.na(ind$Pack)] <- NA
      out$tenure_years <- vapply(packs, function(p) years(at_home(terr[[p]], date)), numeric(1))
    }
  }

  # Within-pack relatedness vs random packs of the same sizes
  if (!is.null(relate) && length(packs)) {
    a <- ind$Individual[!is.na(ind$Pack)]
    R <- matrix(NA_real_, length(a), length(a), dimnames = list(a, a))
    r <- relate[as.character(relate$ind1) %in% a & as.character(relate$ind2) %in% a, ]
    R[cbind(as.character(r$ind1), as.character(r$ind2))] <- r[[estimator]]
    R[cbind(as.character(r$ind2), as.character(r$ind1))] <- r[[estimator]]
    diag(R) <- NA
    Rz <- R; Rz[is.na(Rz)] <- 0; Cn <- (!is.na(R)) * 1
    lab <- ind$Pack[match(a, ind$Individual)]
    within_mean <- function(l) {                                 # mean relatedness within each pack, all packs at once
      P <- outer(l, packs, "==") * 1
      diag(crossprod(P, Rz %*% P)) / diag(crossprod(P, Cn %*% P))
    }
    obs_r <- within_mean(lab)
    perm <- .with_seed(seed, replicate(n_perm, within_mean(sample(lab))))
    perm <- matrix(perm, nrow = length(packs))
    out$mean_r <- round(obs_r, 3)
    out$r_z <- round((obs_r - rowMeans(perm)) / apply(perm, 1, stats::sd), 2)
    out$r_p <- round((1 + rowSums(perm >= obs_r)) / (n_perm + 1), 3)
  }

  rownames(out) <- NULL
  ind[] <- lapply(ind, unname); out[] <- lapply(out, unname)
  attr(out, "individuals") <- ind
  if (!is.null(tl)) attr(out, "periods") <- tl$territories
  out
}

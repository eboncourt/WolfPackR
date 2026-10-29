#' @title Compare a pack assignment with a reference
#' @description Agreement between two ways of putting individuals into groups, e.g. WolfPackR packs vs packs
#'   known from pedigrees or experts, or two runs with different settings. Only individuals present in the
#'   reference are used.
#'   * `ARI`: adjusted Rand index (1 = same partition, about 0 = no better than chance);
#'   * `NMI`: normalised mutual information (1 = same partition, 0 = independent);
#'   * `n_correct`: individuals whose group corresponds to their reference group (each group found is matched
#'     to the reference group most of its members belong to).
#'   Individuals without a group in `x` count as alone (each one its own group) for the ARI and NMI.
#' @param x,reference Either named vectors (names = individuals) or data.frames with `Individual` and a `Pack`
#'   column (or, without `Pack`, the second column). `NA` and `"Lone Individual"` mean no group.
#' @return A list with `summary` (one-row data.frame) and `table` (reference groups x groups found).
#' @export
#' @examples
#' data("samples_Carpathians"); data("relate_Carpathians"); data("reference_Carpathians")
#' res <- wolf_packs(samples_Carpathians, relate_Carpathians, "QG", 0.25)
#' compare_partitions(res$individuals, reference_Carpathians)
compare_partitions <- function(x, reference) {
  x <- .as_partition(x); reference <- .as_partition(reference)
  reference <- reference[!is.na(reference)]
  ids <- names(reference)
  found <- x[ids]                                                 # NA if missing from x
  assigned <- !is.na(found)

  # Each group found stands for the reference group most of its members come from
  maj <- tapply(reference[assigned], found[assigned], function(r) names(which.max(table(r))))
  correct <- assigned & maj[found] %in% reference & maj[found] == reference
  alone <- ifelse(assigned, found, paste0(".alone", seq_along(found)))
  ri <- as.integer(factor(reference)); fi <- as.integer(factor(alone))

  summary <- data.frame(
    n_reference = length(ids), n_assigned = sum(assigned), n_correct = sum(correct, na.rm = TRUE),
    prop_correct = round(mean(correct %in% TRUE), 2),
    ARI = round(igraph::compare(ri, fi, method = "adjusted.rand"), 2),
    NMI = round(igraph::compare(ri, fi, method = "nmi"), 2),
    reference_groups = length(unique(reference)), reference_groups_found = length(unique(maj))
  )
  list(summary = summary, table = table(reference = reference, found = ifelse(assigned, found, "none")))
}

# Named vector of groups from a vector or a data.frame (Individual + Pack or second column)
.as_partition <- function(p) {
  if (is.data.frame(p)) {
    col <- if ("Pack" %in% names(p)) "Pack" else names(p)[2]
    p <- stats::setNames(as.character(p[[col]]), as.character(p$Individual))
  }
  if (is.null(names(p))) stop("Give a named vector or a data.frame with an 'Individual' column.")
  p <- stats::setNames(as.character(p), names(p))
  p[p %in% "Lone Individual"] <- NA
  p
}

#' @title Ugly ducklings (deprecated)
#' @description Deprecated in WolfPackR 1.0.0 and to be removed in 1.1: use [assign_individuals()], which works on
#'   kernel densities and gives every individual an explicit status. This wrapper calls [assign_individuals()]
#'   with `min_share = min_overlap`; `buffer` and `mcp.percent` no longer mean anything and are ignored.
#'   For older scripts, the output also has `Pack` and `new_group` columns (`new_group` = `Lone_1`, `Lone_2`, ...
#'   for the individuals left without a group).
#' @param obs An `sf` object with individuals and their groups.
#' @param group Name of the group column (string, or bare name as in older scripts).
#' @param min_overlap Passed to `min_share` of [assign_individuals()] (default: 0.5).
#' @param buffer,mcp.percent Ignored.
#' @return The output of [assign_individuals()] plus `Pack` and `new_group`.
#' @export
#' @keywords internal
ugly_ducklings <- function(obs, group, min_overlap = 0.5, buffer = 0, mcp.percent = 100) {
  .Deprecated("assign_individuals", package = "WolfPackR",
              msg = "ugly_ducklings() is deprecated and will be removed in WolfPackR 1.1: use assign_individuals() instead.")
  group <- .col_name(substitute(group), function() group, names(obs))
  out <- assign_individuals(obs, group, min_share = min_overlap)
  out$Pack <- out$Group
  lone <- is.na(out$Group)
  out$new_group <- out$Group
  out$new_group[lone] <- paste0("Lone_", seq_len(sum(lone)))
  out
}

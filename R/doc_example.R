#' samples_example
#'
#' A small made-up dataset to try the functions: 47 samples of 9 fictive wolves (W1 to W9), in two families
#' (W1-W3 and W4-W9, see [relate_example]). Formerly called `samples`.
#' @format An `sf` object (POINT, WGS84) with 47 rows and 2 variables:
#' \describe{
#'   \item{Individual}{character, individual ID}
#'   \item{Sex}{character, "Male" or "Female"}
#'   \item{geometry}{sfc_POINT}
#' }
#' @source Simulated.
#' @seealso [samples_Carpathians] for real data.
"samples_example"

#' relate_example
#'
#' Made-up pairwise relatedness between the 9 individuals of [samples_example], in the format expected by
#' [genetic_groups()] and [wolf_packs()]: one row per pair, the two individuals and one column per estimator
#' (here a single one, `wang`, as an estimator of Wang, 2002, would give). Formerly called `relate`.
#' @format A data.frame with 81 rows (all ordered pairs, including each individual with itself) and 3 variables:
#' \describe{
#'   \item{ind1, ind2}{character, the two individuals}
#'   \item{wang}{numeric, relatedness}
#' }
#' @source Simulated.
#' @seealso [relatedness()] to compute relatedness from genotypes, [relate_Carpathians] for real data.
"relate_example"

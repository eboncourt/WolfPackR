#' relate_Carpathians
#'
#' Pairwise relatedness between the 47 individuals of [samples_Carpathians], computed with
#' `relatedness(genotypes_Carpathians, freq_ids = <all but the hybrid EF.15J8>)`. One row per pair.
#' @format A data.frame with 1081 rows and 6 variables:
#' \describe{
#'   \item{ind1, ind2}{character, the two individuals}
#'   \item{QG}{numeric, Queller & Goodnight (1989) relatedness}
#'   \item{LR}{numeric, Lynch & Ritland (1999) relatedness}
#'   \item{LOD}{numeric, log10 likelihood ratio first-degree relatives vs unrelated}
#'   \item{n_loci}{integer, number of loci typed in both individuals}
#' }
#' @source Computed from [genotypes_Carpathians].
"relate_Carpathians"

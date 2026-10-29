#' genotypes_Carpathians
#'
#' Reference genotypes of the 47 individuals of [samples_Carpathians] at 22 microsatellite loci,
#' two columns per locus (allele sizes in bp). Ready to use with [relatedness()].
#' Individual EF.15J8 is a wolf x dog F1 hybrid: leave it out of allele frequencies with `freq_ids`.
#' Built by `data-raw/carpathians.R`.
#' @format A data.frame with 47 rows and 45 variables:
#' \describe{
#'   \item{Individual}{character, individual ID}
#'   \item{CPH9_1, CPH9_2, ..., AHTh171_1, AHTh171_2}{integer, the two alleles at each locus (`NA` if missing)}
#' }
#' @source Foundation Conservation Carpathia (2024). Wolf population size and composition in one of Europe's
#'   strongholds, the Romanian Carpathians (Dataset). Zenodo. \doi{10.5281/zenodo.14544266}
"genotypes_Carpathians"

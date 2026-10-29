#' samples_Carpathians
#'
#' Non-invasive genetic samples of wolves from the Romanian Carpathians (Piatra Craiului area), 2017-2020:
#' the 263 samples genotyped as wolf (262) or wolf x dog F1 hybrid (1, individual EF.15J8), i.e. 47 individuals.
#' Samples of dogs (including CC01CK, put in the dog cluster by Structure although three of its samples are
#' labelled "Wolf" in the source file), foxes and samples that could not be genotyped are left out: see [effort_Carpathians] for all of them.
#' Built by `data-raw/carpathians.R` from the Zenodo dataset of Foundation Conservation Carpathia (2024),
#' analysed by Iosif et al. (2025).
#' @format An `sf` object (POINT, WGS84) with 263 rows and 8 variables:
#' \describe{
#'   \item{Individual}{character, individual ID (genotype reference)}
#'   \item{Sex}{character, genetic sex ("M" or "F")}
#'   \item{Date}{Date, sampling date (missing for 2 samples)}
#'   \item{Campaign}{character, sampling season from July to June ("2017-18", "2018-19" or "2019-20"), from the date, or from the campaign label of the source file when the date is missing}
#'   \item{Sample}{character, sample ID}
#'   \item{Type}{character, sample type: Scat, Urine, Hair, Tissue, Blood or Blood-noninvasive}
#'   \item{pWolf}{numeric, wolf ancestry of the individual (Structure, 2 clusters): 0.51 for the hybrid EF.15J8, >= 0.9 for the others}
#'   \item{geometry}{sfc_POINT}
#' }
#' @source Foundation Conservation Carpathia (2024). Wolf population size and composition in one of Europe's
#'   strongholds, the Romanian Carpathians (Dataset). Zenodo. \doi{10.5281/zenodo.14544266}
#' @references Iosif, R., Skrbinšek, T., Erős, N., Konec, M., Boljte, B., Jan, M., & Promberger-Fürpass, B. (2025).
#'   Wolf population size and composition in one of Europe's strongholds, the Romanian Carpathians.
#'   Ecology and Evolution, 15(4), e71200. \doi{10.1002/ece3.71200}
"samples_Carpathians"

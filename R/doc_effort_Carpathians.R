#' effort_Carpathians
#'
#' All 505 samples collected for the Carpathian wolf survey (2017-2020), genotyped or not: handy to map the
#' sampling effort and the genotyping success, e.g. to see whether an area was little sampled or sampled
#' with poor success. The genotyped wolf samples are the ones in [samples_Carpathians].
#' Built by `data-raw/carpathians.R`.
#' @format An `sf` object (POINT, WGS84) with 505 rows and 8 variables:
#' \describe{
#'   \item{Sample}{character, sample ID}
#'   \item{Date}{Date, sampling date (missing for 2 samples)}
#'   \item{Campaign}{character, sampling season from July to June ("2017-18", "2018-19" or "2019-20")}
#'   \item{Type}{character, sample type: Scat, Urine, Hair, Tissue, Blood or Blood-noninvasive}
#'   \item{Verdict}{character, genotyping outcome: "Genotyped", "Poor" (failed), "Mixed" (several animals) or "NonTarget" (dog, fox)}
#'   \item{Species}{character, "Wolf", "Hybrid", "Dog", "Fox", or "Unknown" when not genotyped}
#'   \item{Individual}{character, individual ID for genotyped samples, `NA` otherwise}
#'   \item{geometry}{sfc_POINT}
#' }
#' @source Foundation Conservation Carpathia (2024). Wolf population size and composition in one of Europe's
#'   strongholds, the Romanian Carpathians (Dataset). Zenodo. \doi{10.5281/zenodo.14544266}
"effort_Carpathians"

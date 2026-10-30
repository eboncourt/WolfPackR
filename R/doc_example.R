#' samples_example
#'
#' Non-invasive genetic samples of a simulated wolf population followed for four biological years (July to June,
#' 2020-21 to 2023-24): 343 samples of 56 individuals. Since the population is simulated, the truth is known
#' ([truth_example]): handy to see what WolfPackR gets right and wrong, to try settings, or to teach.
#'
#' The story: 8 packs at the start (A to H in [truth_example]), 15 km apart, some of them related through their breeders;
#' a daughter of pack A founds pack I next door with a male from outside (one family, two territories); both breeders
#' of pack C die and one of their sons becomes the breeding male of pack F (a dated dispersal), while a daughter
#' becomes a floater; the breeding male of B and the breeding female of G are replaced by animals from outside;
#' a transient crosses the area between the territories and a few others leave a single sample. Adults leave 3% of
#' their samples in a neighbouring territory, sampling is lighter in 2022-23 and pack E, on the edge, is poorly sampled.
#' Coordinates are placed in the French Alps only to get a background map: they mean nothing.
#'
#' With the default settings, `wolf_packs(samples_example, relate_example, "QG", 0.25, period = "Campaign")` puts
#' 48 of the 51 pack members with their pack mates (adjusted Rand index 0.94 with the true packs, the median over
#' 40 populations simulated the same way). Built by `data-raw/simulate_example.R`; `simulate_population()`, in
#' `data-raw/simulate_population.R` on GitHub, simulates more populations with the same story.
#' Formerly a made-up set of 47 samples of 9 wolves called `samples`.
#' @format An `sf` object (POINT, WGS84) with 343 rows and 7 variables:
#' \describe{
#'   \item{Individual}{character, individual ID ("W001" to "W056", in order of first sampling)}
#'   \item{Sex}{character, "M" or "F"}
#'   \item{Date}{Date, sampling date}
#'   \item{Campaign}{character, sampling season from July to June ("2020-21" to "2023-24"), from the date}
#'   \item{Sample}{character, sample ID}
#'   \item{Type}{character, sample type: Scat, Urine, Hair, Saliva or Tissue (wolves found dead)}
#'   \item{geometry}{sfc_POINT}
#' }
#' @source Simulated (`data-raw/simulate_example.R`).
#' @seealso [genotypes_example], [relate_example], [truth_example]; [samples_Carpathians] for real data.
#' @examples
#' data("samples_example"); data("relate_example"); data("truth_example")
#' res <- wolf_packs(samples_example, relate_example, "QG", 0.25, period = "Campaign")
#' compare_partitions(res$individuals, truth_example)
"samples_example"

#' genotypes_example
#'
#' Simulated genotypes of the 56 individuals of [samples_example] at 22 microsatellite loci (made-up loci and allele
#' sizes, 5 to 10 alleles per locus, Mendelian inheritance along the pedigree of [truth_example]), two columns per
#' locus, with 4% of the loci not typed and rare genotyping errors. Ready to use with [relatedness()].
#' @format A data.frame with 56 rows and 45 variables:
#' \describe{
#'   \item{Individual}{character, individual ID}
#'   \item{MS01_1, MS01_2, ..., MS22_1, MS22_2}{integer, the two alleles at each locus (`NA` if missing)}
#' }
#' @source Simulated (`data-raw/simulate_example.R`).
"genotypes_example"

#' relate_example
#'
#' Pairwise relatedness between the 56 individuals of [samples_example], computed with
#' `relatedness(genotypes_example)`. One row per pair, in the format expected by [genetic_groups()] and [wolf_packs()].
#' Formerly a made-up table of 9 wolves with a single `wang` column.
#' @format A data.frame with 1540 rows and 6 variables:
#' \describe{
#'   \item{ind1, ind2}{character, the two individuals}
#'   \item{QG}{numeric, Queller & Goodnight (1989) relatedness}
#'   \item{LR}{numeric, Lynch & Ritland (1999) relatedness}
#'   \item{LOD}{numeric, log10 likelihood ratio first-degree relatives vs unrelated}
#'   \item{n_loci}{integer, number of loci typed in both individuals}
#' }
#' @source Computed from [genotypes_example].
"relate_example"

#' truth_example
#'
#' What really happened in the simulated population of [samples_example], one row per individual. `Pack` can be
#' given straight to [compare_partitions()] as the reference; the yearly columns can be set against [pack_timeline()].
#' @format A data.frame with 56 rows and 14 variables:
#' \describe{
#'   \item{Individual}{character, individual ID}
#'   \item{Sex}{character, "M" or "F"}
#'   \item{Pack}{character, pack ("A" to "I") where most of its samples come from, i.e. where it lived when sampled;
#'     `NA` for floaters and transients}
#'   \item{Role}{character, "Breeder" (breeder at some point of the study), "Offspring" (young living in its natal
#'     pack) or "Floater" (no territory)}
#'   \item{Origin}{character, "Founder" (breeder already there at the start), "Born" (born in a pack of the study
#'     area) or "Immigrant" (came from outside the study area)}
#'   \item{Natal_Pack}{character, pack where it was born (`NA` if born outside the study area or before it)}
#'   \item{Mother, Father}{character, its parents when they were sampled, `NA` otherwise}
#'   \item{Birth_Year}{integer, year of birth (May) for the individuals born in the study area}
#'   \item{Dispersal}{Date, when it left its natal pack for good, if during the study}
#'   \item{Pack_2020_21, Pack_2021_22, Pack_2022_23, Pack_2023_24}{character, where it lived for the longest time
#'     each year: a pack, "floater", or `NA` (not born yet, dead, or out of the study area)}
#' }
#' @source Simulated (`data-raw/simulate_example.R`).
"truth_example"

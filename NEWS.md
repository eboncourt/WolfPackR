# WolfPackR 1.1.1

## Bug fixes

* `isopleth()` (and so `plot_packs()` and the maps drawn from it): with some bandwidths, rounding left hairline gaps between the columns of grid cells, so a territory came out as many vertical strips, each with its own outline. Squares are now merged in grid units, then put back in metres: one outline per territory. Areas are unchanged.
* `spatial_groups()`: when a group had exactly two core individuals, they made one pack even if their kernels didn't overlap (with three or more, the overlap graph and `min_affinity` already decided). They now need an affinity of at least `min_affinity`. This mostly concerns the pool of individuals without a group, where two unrelated wanderers far apart could make a "new" pack. No change on the Carpathian data or the simulated population.

# WolfPackR 1.1.0

## Changes

* Wide-ranging individuals (`wide` in `assign_individuals()`, `spatial_groups()` and `wolf_packs()`): the radius of gyration now leaves out the farthest 10% of the samples (at least one, keeping at least three). A single foray into a neighbouring territory no longer makes a resident look wide-ranging, which kept it from drawing its territory and could make a small pack disappear; an animal spread over several places is still flagged. On the Carpathian data, one more wolf joins its pack (CC04FK, formerly transient): ARI with Iosif et al. (2025) 0.80 instead of 0.77, and the same result with or without `period` (0.68 without it before). Over 40 simulated populations, residents flagged as wide-ranging dropped from 6.6 to 3.0 on average. `pack_timeline()` still reports the plain radius of gyration.
* Every threshold can be set from `wolf_packs()`: new arguments `outside` (when a sample is outside every territory; until now only an argument of `attribute_samples()`) and `wide_trim` (share of the farthest samples left out of the radius of gyration; 0 = as before) in `assign_individuals()`, `spatial_groups()` and `wolf_packs()`, `outside` in `pack_timeline()` and `summarise_packs()`, and `max_iter` in `wolf_packs()`. `pack_stability(vary = )` can test them. `wolf_packs()` returns the settings used (`$settings`). The vignette lists all the settings.

## Data

* New simulated example population, with the truth: `samples_example` (343 samples of 56 wolves over four biological years, with dates, seasons and sample types), `genotypes_example` (22 microsatellites), `relate_example` (now computed by `relatedness()`: `QG`, `LR` and `LOD` columns instead of a made-up `wang` one) and `truth_example` (true packs, roles, origins, parents, dispersal dates and pack each year). It replaces the 9-wolf toy dataset and covers every status, families spread over two territories, a pack that disappears, a pack founded during the study, a dated dispersal and breeder turnover. Built by `data-raw/simulate_example.R`; `simulate_population()` (`data-raw/simulate_population.R`) simulates more populations with the same story.
* The vignette compares WolfPackR results with the truth of the simulated population (packs, statuses, land tenure, lighter sampling).

# WolfPackR 1.0.0

Territories are now estimated by kernel density instead of minimum convex polygons, and the workflow gives every individual an explicit status.

## Breaking changes

* `genetic_groups()` uses community detection (`method = "leiden"`) by default; `method = "components"` gives the former behaviour.
* `spatial_groups()` works on kernel densities. New arguments: `ud`, `h`, `method`, `min_affinity`, `merge`, `min_samples`, `min_group_size`, `wide`, `period`, `min_share`, `seed`. Gone: `percentile`, `buffer_radius`, `max_iterations`, `min_mcp_overlap`. It returns `Individual`, `Group`, `Pack` (`NA` = no pack) and `Core` instead of `Subgroup` / "Lone Individual". Pairs are valid packs (`min_group_size = 2`).
* `summarise_packs()` returns a data.frame with one row per pack (was a list). Putative dominant individuals are gone: on the Carpathian data, the most sampled individual of each sex was the known breeder in only half of the cases. Areas are kernel isopleths (95% and 50%). New columns on the space shared with other packs (`shared_pct`) and the nearest pack (`nearest_pack`, `nearest_km`); with the new `period` and `date` arguments, land tenure metrics (Kanishka et al., 2026): first and last period in the territory, members detected in several periods, stability and shift of the territory, time from first to last detection.
* `plot_packs()` draws kernel territories (95% outline, 50% core area); individuals without a pack are white.
* Column names are given as text (e.g. `"Pack"`); bare names still work.

## New functions

* `space_use()`: kernel UD of every individual on a common grid, with `sf` and base R; bandwidth from the well-sampled individuals or set by the user.
* `ud_overlap()`: Bhattacharyya affinity, volume of intersection, PHR and HR (Fieberg & Kochanny 2005).
* `isopleth()`: territory and core-area polygons of an individual or a group.
* `attribute_samples()`: the territory each sample falls in, or "outside".
* `assign_individuals()` replaces `ugly_ducklings()`: majority rule on the samples, iterated until stable, with statuses Resident, Putative disperser, Disperser (dated with `period`), Unrelated resident, Transient and Undetermined; relatives defined pair by pair; spread measured within each period; a low-data rule.
* `pack_timeline()`: individuals and territories from one sampling period to the next (members detected, new and missed, stability, shift and overlap of each territory in each period).
* `compare_partitions()`: ARI, NMI and share of individuals correctly placed against a reference.
* `pack_stability()`: sensitivity to the settings, bootstrap of the samples and subsampling (`fraction`: keep only part of the samples, to check whether the sampling is dense enough).
* `wolf_packs()`: the whole workflow in one call.
* Null models: `genetic_groups(n_null = )` tests the modularity against random graphs (z-score); `summarise_packs(relate = )` tests the within-pack relatedness by permutation.

## Data

* The example datasets `samples` and `relate` are renamed `samples_example` and `relate_example` (`relate` was also the name of an argument, and close to the `related` package).
* New `reference_Carpathians`: packs of 41 wolves read on Fig. 3 of Iosif et al. (2025).

## Deprecated

* `ugly_ducklings()` (now a wrapper of `assign_individuals()`), `mcp_sf()` and `calculate_overlap()` still work with a warning and will be removed in 1.1.

# WolfPackR 0.9.5

* New `relatedness()`: pairwise relatedness from microsatellite genotypes in base R (Queller & Goodnight, Lynch & Ritland, and a LOD first-degree vs unrelated), no COLONY / COANCESTRY / `related` needed.
* `genetic_groups()` gains `method` (`"components"`, the old behaviour and still the default, or community detection with `"leiden"`, `"louvain"`, `"walktrap"`), `seed` and `plot`. The graph is no longer drawn by default, and having no pair above the threshold no longer fails. Pairs kept with a relatedness <= 0 (threshold < 0) now raise an error, as igraph versions handle such links differently.
* Carpathian datasets rebuilt from the Zenodo file (`data-raw/carpathians.R`):
  * `samples_Carpathians` now holds only the 263 samples genotyped as wolf or hybrid (47 individuals). The samples of 11 dogs (including CC01CK, labelled "Wolf" in the source file but put in the dog cluster by Structure) and 1 sample absent from the genotyped list are gone. New columns: `Date`, `Campaign` (July-June season), `Sample`, `Type`, `pWolf`.
  * New `effort_Carpathians`: the 502 collected samples that have a position (3 recorded at 0, 0 are left out), genotyped or not, to map sampling effort and genotyping success.
  * New `genotypes_Carpathians` (22 microsatellites) and `relate_Carpathians` (QG, LR, LOD).
  * `Bruvo_Carpathians` restricted to the 47 wolves, kept for comparison with the first version of the paper.

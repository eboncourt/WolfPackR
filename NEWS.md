# WolfPackR 0.9.5

* New `relatedness()`: pairwise relatedness from microsatellite genotypes in base R (Queller & Goodnight, Lynch & Ritland, and a LOD first-degree vs unrelated), no COLONY / COANCESTRY / `related` needed.
* `genetic_groups()` gains `method` (`"components"`, the old behaviour and still the default, or community detection with `"leiden"`, `"louvain"`, `"walktrap"`), `seed` and `plot`. The graph is no longer drawn by default, and having no pair above the threshold no longer fails. Pairs kept with a relatedness <= 0 (threshold < 0) now raise an error, as igraph versions handle such links differently.
* Carpathian datasets rebuilt from the Zenodo file (`data-raw/carpathians.R`):
  * `samples_Carpathians` now holds only the 263 samples genotyped as wolf or hybrid (47 individuals). The samples of 11 dogs (including CC01CK, labelled "Wolf" in the source file but put in the dog cluster by Structure) and 1 sample absent from the genotyped list are gone. New columns: `Date`, `Campaign` (July-June season), `Sample`, `Type`, `pWolf`.
  * New `effort_Carpathians`: all 505 collected samples, genotyped or not, to map sampling effort and genotyping success.
  * New `genotypes_Carpathians` (22 microsatellites) and `relate_Carpathians` (QG, LR, LOD).
  * `Bruvo_Carpathians` restricted to the 47 wolves, kept for comparison with the first version of the paper.

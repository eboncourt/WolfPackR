# Builds the simulated example datasets of the package: samples_example, genotypes_example, relate_example, truth_example
# (the population simulated by simulate_population() with seed 22, see simulate_population.R for the story).
# Run from the package root: source("data-raw/simulate_example.R"). Fixed seed, so same data every time.
library(sf)
pkgload::load_all()                                     # for relatedness(), and wolf_packs() in the check below
source("data-raw/simulate_population.R")

# ---- The datasets of the package (seed 22)
ex <- simulate_population(seed = 22)
list2env(ex, environment())
cat(nrow(samples_example), "samples of", nrow(truth_example), "individuals\n")
for (nm in names(ex)) save(list = nm, file = file.path("data", paste0(nm, ".rda")), compress = "xz", version = 2)

# ---- How does WolfPackR do over many populations simulated the same way? (40 populations, a few minutes)
# Seed 22 was picked because every status shows up and its agreement with the truth is typical of the 40 runs.
if (FALSE) {
  reps <- do.call(rbind, lapply(1:40, function(s) {
    x <- simulate_population(seed = s)
    res <- suppressMessages(wolf_packs(x$samples_example, x$relate_example, "QG", 0.25, period = "Campaign"))
    cbind(seed = s, n_individuals = nrow(x$truth_example), n_samples = nrow(x$samples_example),
          compare_partitions(res$individuals, x$truth_example)$summary)
  }))
  summary(reps$ARI)
}

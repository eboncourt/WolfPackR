# Comparison with the submitted version (0.9.4, minimum convex polygons): ARI = 0.52, 24/41 correctly placed
# Version 0.9.4 is installed in a separate library and runs in a separate R session (callr), leaving the current version alone

lib094 <- file.path(tempdir(), "lib_WolfPackR_094")
dir.create(lib094, showWarnings = FALSE)
remotes::install_github("eboncourt/WolfPackR", ref = "53861872bbc67148e48ae919e300f6f462051e0f",   # last commit of 0.9.4 on master
                        lib = lib094, upgrade = "never", dependencies = FALSE)

# The former workflow, with the settings of the first version of the article
old <- callr::r(function() {
  library(WolfPackR); library(sf)
  grDevices::pdf(NULL)                                            # the old genetic_groups() draws the graph
  data("Bruvo_Carpathians"); data("samples_Carpathians")
  gg <- genetic_groups(Bruvo_Carpathians, threshold = 0.7, estimator = "1-Bruvo")
  ud0 <- suppressWarnings(ugly_ducklings(merge(samples_Carpathians, gg, by = "Individual"), group = group,
                                         min_overlap = 0.5, buffer = 0, mcp.percent = 100))
  sg0 <- suppressWarnings(spatial_groups(merge(samples_Carpathians, ud0, by = "Individual"), group = "new_group",
                                         percentile = 100, buffer_radius = 1e-7, max_iterations = 20, min_mcp_overlap = 0.5))
  c(packageVersion("WolfPackR") == "0.9.4", setNames(sg0$Subgroup, sg0$Individual))
}, libpath = c(lib094, .libPaths()))
stopifnot("The separate session did not load version 0.9.4" = old[1] == "TRUE")
old <- old[-1]
old[old %in% "Lone Individual"] <- NA

# And the comparison, with the current version
library(WolfPackR)
data("reference_Carpathians")
compare_partitions(old, reference_Carpathians)

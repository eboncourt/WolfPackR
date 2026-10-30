# Comparaison avec la version soumise (0.9.4, MCP) : ARI = 0.52, 24/41 bien placés
# La 0.9.4 est installée dans une bibliothèque à part et tourne dans une session R séparée (callr), pour ne pas toucher à la 1.0

lib094 <- file.path(tempdir(), "lib_WolfPackR_094")
dir.create(lib094, showWarnings = FALSE)
remotes::install_github("eboncourt/WolfPackR", ref = "53861872bbc67148e48ae919e300f6f462051e0f",   # dernier commit de la 0.9.4 sur master
                        lib = lib094, upgrade = "never", dependencies = FALSE)

# L'ancienne chaîne, avec les réglages de la première version de l'article
old <- callr::r(function() {
  library(WolfPackR); library(sf)
  grDevices::pdf(NULL)                                            # l'ancienne genetic_groups() dessine le graphe
  data("Bruvo_Carpathians"); data("samples_Carpathians")
  gg <- genetic_groups(Bruvo_Carpathians, threshold = 0.7, estimator = "1-Bruvo")
  ud0 <- suppressWarnings(ugly_ducklings(merge(samples_Carpathians, gg, by = "Individual"), group = group,
                                         min_overlap = 0.5, buffer = 0, mcp.percent = 100))
  sg0 <- suppressWarnings(spatial_groups(merge(samples_Carpathians, ud0, by = "Individual"), group = "new_group",
                                         percentile = 100, buffer_radius = 1e-7, max_iterations = 20, min_mcp_overlap = 0.5))
  c(packageVersion("WolfPackR") == "0.9.4", setNames(sg0$Subgroup, sg0$Individual))
}, libpath = c(lib094, .libPaths()))
stopifnot("La session séparée n'a pas chargé la 0.9.4" = old[1] == "TRUE")
old <- old[-1]
old[old %in% "Lone Individual"] <- NA

# Et la comparaison, avec la 1.0
library(WolfPackR)
data("reference_Carpathians")
compare_partitions(old, reference_Carpathians)

# Fig. 1 de l'article : le workflow de WolfPackR 1.0 (remplace l'ancienne Fig. 2)
# Juste ggplot2 : chaque boîte a ses coordonnées (x, y = centre ; w, h = largeur, hauteur), faciles à déplacer

library(ggplot2)
out <- "article/output"                                  # dossier de sortie, relatif à la racine du dépôt
dir.create(out, showWarnings = FALSE, recursive = TRUE)
largeur <- 17                                            # cm, pleine largeur de page

# ---- Fig. 1 : le workflow ------------------------------------------------------------------------------------
# Des boîtes et des flèches en ggplot2 : coordonnées à la main, faciles à déplacer
boite <- function(x, y, w, h, titre, texte, type) data.frame(x = x, y = y, w = w, h = h, titre = titre, texte = texte, type = type)
b1 <- rbind(
  boite(2.3, 9.45, 3.6, 0.8, "Genotypes", "microsatellites, one row per individual", "Input"),
  boite(7.0, 9.45, 3.6, 0.8, "Sample locations", "sf points: individual, sex, date, period", "Input"),
  boite(2.3, 8.2, 3.6, 0.9, "relatedness()", "QG, LR and LOD estimators\n(or any pairwise relatedness table)", "Genetics"),
  boite(2.3, 6.55, 3.6, 1.15, "1. genetic_groups()", "families: community detection on the\nrelatedness graph (Leiden by default),\nmodularity tested against a null model", "Main"),
  boite(7.0, 6.55, 3.6, 1.15, "2. space_use()", "kernel density of every individual\non a common grid, with a common\nbandwidth h", "Main"),
  boite(4.65, 4.95, 3.9, 1.15, "3. assign_individuals()", "each individual to the territory where\nmost of its samples fall; a status\nfrom where its relatives are", "Main"),
  boite(4.65, 3.4, 3.9, 0.95, "4. spatial_groups()", "families split into spatially coherent\npacks (UD overlap, communities)", "Main"),
  boite(1.35, 1.35, 2.5, 1.25, "summarise_packs()", "composition, areas,\nshared space, tenure,\nkinship test", "Output"),
  boite(3.85, 1.35, 2.3, 1.25, "plot_packs()", "interactive map,\n95% and 50%\nisopleths", "Output"),
  boite(6.3, 1.35, 2.4, 1.25, "pack_timeline()", "individuals and\nterritories, season\nby season", "Output"),
  boite(9.6, 1.35, 3.6, 1.25, "compare_partitions()\npack_stability()", "agreement with a reference;\nsensitivity, bootstrap, subsampling", "Assessment"),
  boite(10.6, 4.95, 2.6, 1.95, "Statuses", "Resident\nPutative disperser\nDisperser (dated)\nUnrelated resident\nTransient\nUndetermined", "Status"),
  boite(10.6, 8.2, 2.6, 1.25, "Building blocks", "ud_overlap()\nisopleth()\nattribute_samples()", "Auxiliary")
)
fl <- data.frame(x = c(2.3, 2.3, 7.0, 2.3, 7.0, 6.6, 4.65, 4.65, 4.65, 4.65, 4.65),
                 y = c(9.05, 7.75, 9.05, 5.975, 5.975, 4.95, 4.375, 2.925, 2.925, 2.925, 2.925),
                 xend = c(2.3, 2.3, 7.0, 3.6, 5.7, 9.3, 4.65, 1.35, 3.85, 6.3, 9.6),
                 yend = c(8.65, 7.125, 7.125, 5.525, 5.525, 4.95, 3.875, 1.975, 1.975, 1.975, 1.975))
coul <- c(Input = "#B8E3B0", Genetics = "#D9EAD3", Main = "#FFF2B3", Output = "#CFE2F3", Assessment = "#E4D4F0",
          Status = "#FFFFFF", Auxiliary = "#F9CB9C")
fig1 <- ggplot() +
  annotate("rect", xmin = 0.25, xmax = 9.05, ymin = 2.7, ymax = 7.4, fill = NA, colour = "grey45", linetype = "dashed") +
  annotate("text", x = 8.95, y = 2.85, label = "wolf_packs(): steps 1 to 4 in one call", hjust = 1, size = 2.6, colour = "grey30", fontface = "italic") +
  geom_segment(data = fl, aes(x, y, xend = xend, yend = yend), arrow = arrow(length = unit(1.6, "mm"), type = "closed"), colour = "grey35", linewidth = 0.4) +
  geom_rect(data = b1, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = type), colour = "grey40", linewidth = 0.3) +
  geom_text(data = b1, aes(x, y + h / 2 - 0.12, label = titre), vjust = 1, size = 2.8, fontface = "bold", family = "mono", lineheight = 0.9) +
  geom_text(data = b1, aes(x, y + h / 2 - ifelse(grepl("\n", titre), 0.62, 0.38), label = texte), vjust = 1, size = 2.35, lineheight = 0.95) +
  scale_fill_manual(values = coul, breaks = c("Input", "Main", "Output", "Assessment", "Auxiliary"),
                    labels = c("Input data", "Workflow", "Description", "Assessment", "Building blocks"), name = NULL) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0.6, 9.95), expand = FALSE) +
  theme_void(base_size = 8) + theme(legend.position = "bottom", plot.margin = margin(4, 4, 4, 4))
ggsave(file.path(out, "Fig1_workflow.png"), fig1, width = largeur, height = 13.5, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(out, "Fig1_workflow.pdf"), fig1, width = largeur, height = 13.5, units = "cm")

# Fig. 1 of the article on its own: the WolfPackR workflow, with the package logo in the top right corner
# Just ggplot2: each box has its coordinates (x, y = centre; w, h = width, height), easy to move. Run from the root of the repository.
# Saved as PNG and TIFF (600 dpi), PDF and SVG (vector; the SVG opens in Inkscape or Illustrator, text still editable)

library(ggplot2)
library(patchwork)                                       # to put the logo in the corner
out <- "article/output"                                  # output folder, relative to the root of the repository
dir.create(out, showWarnings = FALSE, recursive = TRUE)
width_cm <- 17; height_cm <- 13.5                        # full page width
logo <- png::readPNG("man/figures/logo.png")
sans <- "Liberation Sans"; mono <- "Liberation Mono"     # same widths as Arial and Courier New, the fonts written in the SVG

# ---- Boxes and arrows ------------------------------------------------------------------------------------------
box <- function(x, y, w, h, title, body, type) data.frame(x = x, y = y, w = w, h = h, title = title, body = body, type = type)
boxes <- rbind(
  box(2.3, 9.45, 3.6, 0.8, "Genotypes", "microsatellites, one row per individual", "Input"),
  box(7.0, 9.45, 3.6, 0.8, "Sample locations", "sf points: individual, sex, date, period", "Input"),
  box(2.3, 8.2, 3.6, 0.9, "relatedness()", "QG, LR and LOD estimators\n(or any pairwise relatedness table)", "Genetics"),
  box(2.3, 6.55, 3.6, 1.15, "1. genetic_groups()", "families: community detection on the\nrelatedness graph (Leiden by default),\nmodularity tested against a null model", "Main"),
  box(7.0, 6.55, 3.6, 1.15, "2. space_use()", "kernel density of every individual\non a common grid, with a common\nbandwidth h", "Main"),
  box(4.65, 4.95, 3.9, 1.15, "3. assign_individuals()", "each individual to the territory where\nmost of its samples fall; a status\nfrom where its relatives are", "Main"),
  box(4.65, 3.4, 3.9, 0.95, "4. spatial_groups()", "families split into spatially coherent\npacks (UD overlap, communities)", "Main"),
  box(1.35, 1.35, 2.5, 1.25, "summarise_packs()", "composition, areas,\nshared space, tenure,\nkinship test", "Output"),
  box(3.85, 1.35, 2.3, 1.25, "plot_packs()", "interactive map,\nterritories with or\nwithout overlap", "Output"),
  box(6.3, 1.35, 2.4, 1.25, "pack_timeline()", "individuals and\nterritories, season\nby season", "Output"),
  box(9.6, 1.35, 3.6, 1.25, "compare_partitions()\npack_stability()", "agreement with a reference;\nsensitivity, bootstrap, subsampling", "Assessment"),
  box(10.6, 4.95, 2.6, 1.95, "Statuses", "Resident\nPutative disperser\nDisperser (dated)\nUnrelated resident\nTransient\nUndetermined", "Status"),
  box(10.6, 7.1, 2.6, 1.55, "Building blocks", "ud_overlap()\nisopleth()\npack_territories()\nattribute_samples()", "Auxiliary")   # lowered to leave the corner to the logo
)
arrows <- data.frame(x = c(2.3, 2.3, 7.0, 2.3, 7.0, 6.6, 4.65, 4.65, 4.65, 4.65, 4.65),
                     y = c(9.05, 7.75, 9.05, 5.975, 5.975, 4.95, 4.375, 2.925, 2.925, 2.925, 2.925),
                     xend = c(2.3, 2.3, 7.0, 3.6, 5.7, 9.3, 4.65, 1.35, 3.85, 6.3, 9.6),
                     yend = c(8.65, 7.125, 7.125, 5.525, 5.525, 4.95, 3.875, 1.975, 1.975, 1.975, 1.975))
fills <- c(Input = "#B8E3B0", Genetics = "#D9EAD3", Main = "#FFF2B3", Output = "#CFE2F3", Assessment = "#E4D4F0",
           Status = "#FFFFFF", Auxiliary = "#F9CB9C")

# ---- Fig. 1 -----------------------------------------------------------------------------------------------------
fig1 <- ggplot() +
  annotate("rect", xmin = 0.25, xmax = 9.05, ymin = 2.7, ymax = 7.4, fill = NA, colour = "grey45", linetype = "dashed") +
  annotate("text", x = 8.95, y = 2.85, label = "wolf_packs(): steps 1 to 4 in one call", hjust = 1, size = 2.6, colour = "grey30",
           fontface = "italic", family = sans) +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), arrow = arrow(length = unit(1.6, "mm"), type = "closed"), colour = "grey35", linewidth = 0.4) +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = type), colour = "grey40", linewidth = 0.3) +
  geom_text(data = boxes, aes(x, y + h / 2 - 0.12, label = title), vjust = 1, size = 2.8, fontface = "bold", family = mono, lineheight = 0.9) +
  geom_text(data = boxes, aes(x, y + h / 2 - ifelse(grepl("\n", title), 0.62, 0.38), label = body), vjust = 1, size = 2.35, lineheight = 0.95, family = sans) +
  scale_fill_manual(values = fills, breaks = c("Input", "Main", "Output", "Assessment", "Auxiliary"),
                    labels = c("Input data", "Workflow", "Description", "Assessment", "Building blocks"), name = NULL) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0.6, 9.95), expand = FALSE) +
  theme_void(base_size = 8, base_family = sans) + theme(legend.position = "bottom", plot.margin = margin(4, 4, 4, 4))

# Logo in the top right corner, 2.3 cm high (its proportions are kept)
logo_h <- 2.3; logo_w <- logo_h * ncol(logo) / nrow(logo)
fig1_logo <- fig1 + inset_element(grid::rasterGrob(logo, interpolate = TRUE),
                                  left = 1 - (logo_w + 0.15) / width_cm, right = 1 - 0.15 / width_cm,
                                  bottom = 1 - (logo_h + 0.12) / height_cm, top = 1 - 0.12 / height_cm,
                                  align_to = "full", clip = FALSE, on_top = TRUE)

# ---- Files ------------------------------------------------------------------------------------------------------
f <- file.path(out, "Fig1_workflow")
ggsave(paste0(f, ".png"), fig1_logo, width = width_cm, height = height_cm, units = "cm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(paste0(f, ".tiff"), fig1_logo, width = width_cm, height = height_cm, units = "cm", dpi = 600, device = ragg::agg_tiff,
       compression = "lzw", bg = "white")                                       # what journals usually ask for
ggsave(paste0(f, ".pdf"), fig1_logo, width = width_cm, height = height_cm, units = "cm", device = cairo_pdf, bg = "white")   # fonts embedded
ggsave(paste0(f, ".svg"), fig1_logo, width = width_cm, height = height_cm, units = "cm", device = svglite::svglite, bg = "white",
       fix_text_size = FALSE)                                                    # plain text, easy to edit
svg <- readLines(paste0(f, ".svg"), warn = FALSE)                               # in the SVG: Arial and Courier New first, so it looks the same on any computer
svg <- gsub('font-family: *"Liberation Sans"', 'font-family: Arial, "Liberation Sans", Helvetica, sans-serif', svg)
svg <- gsub('font-family: *"Liberation Mono"', 'font-family: "Courier New", "Liberation Mono", monospace', svg)
writeLines(svg, paste0(f, ".svg"))

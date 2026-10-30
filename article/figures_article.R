# Figures, tables and numbers of the WolfPackR article (ESE revision): this is Appendix S1
# Needs WolfPackR 1.1.0 (radius of gyration without the farthest 10% of the samples). Everything is recomputed from the
# data shipped with the package; about 2 to 3 minutes in all. Run from the root of the repository.
# Fig. 1: workflow; Fig. 2: maps compared with Iosif et al. (2025); Fig. 3: robustness;
# Table 2: the packs; Tables S1-S2: appendix; numbers_in_text.txt: the numbers cited in the text

library(WolfPackR)
library(sf)
library(ggplot2)
library(patchwork)

out <- "article/output"                                  # output folder, relative to the root of the repository
dir.create(out, showWarnings = FALSE, recursive = TRUE)
width_cm <- 17                                           # full page width

data("samples_Carpathians"); data("effort_Carpathians"); data("relate_Carpathians"); data("reference_Carpathians")
rel <- relate_Carpathians
ref <- reference_Carpathians
ref_p <- setNames(ref$Pack, ref$Individual)
numbers <- character(0)
note <- function(...) { txt <- paste0(...); numbers <<- c(numbers, txt); message(txt) }   # keeps track of the numbers cited

# Colours of the packs of Iosif et al. (Okabe-Ito, colour-blind friendly); other groups get extra colours
col_ref <- c(a = "#E69F00", b = "#56B4E9", c = "#009E73", d = "#CC79A7", e = "#D55E00", f = "#0072B2")
col_extra <- c("#8C6D31", "#6A3D9A", "#B15928", "#33A02C", "#E31A1C")
legend_inside <- function(x, y) if (packageVersion("ggplot2") >= "3.5.0")    # legend inside the panel, whatever the ggplot2 version
  theme(legend.position = "inside", legend.position.inside = c(x, y)) else theme(legend.position = c(x, y))


# ---- 1. The data --------------------------------------------------------------------------------------------
note("Samples collected with a position: ", nrow(effort_Carpathians), "; genotyped as wolf: ", sum(effort_Carpathians$Verdict == "Genotyped"))
note("Wolf samples: ", nrow(samples_Carpathians), "; individuals: ", length(unique(samples_Carpathians$Individual)))
note("By season: ", paste(names(table(samples_Carpathians$Campaign)), table(samples_Carpathians$Campaign), collapse = ", "))
note("By type: ", paste(names(table(samples_Carpathians$Type)), table(samples_Carpathians$Type), collapse = ", "))
n_ind <- table(samples_Carpathians$Individual)
note("Samples per individual: median ", median(n_ind), ", 1 sample: ", sum(n_ind == 1), ", 2: ", sum(n_ind == 2), ", >= 10: ", sum(n_ind >= 10))
lat <- st_coordinates(effort_Carpathians)[, 2]
note("Genotyping success in the south (< 45.45°N): ", round(100 * mean(effort_Carpathians$Verdict[lat < 45.45] == "Genotyped")),
     "%; in the north: ", round(100 * mean(effort_Carpathians$Verdict[lat >= 45.45] == "Genotyped")), "%")


# ---- 2. The workflow step by step, then in one call ---------------------------------------------------------
gen <- genetic_groups(rel, "QG", 0.25, n_null = 199)
mod <- attr(gen, "modularity")
note("Genetic groups (Leiden, QG > 0.25): sizes ", paste(sort(table(gen$group), decreasing = TRUE), collapse = " "),
     "; modularity Q = ", round(mod["Q"], 2), ", z = ", round(mod["z"], 1))
ud <- space_use(samples_Carpathians)
note("Bandwidth h = ", round(ud$h), " m")

res <- wolf_packs(samples_Carpathians, rel, "QG", 0.25, period = "Campaign", ud = ud)
ind <- res$individuals
note("Statuses: ", paste(names(table(ind$Status)), table(ind$Status), collapse = ", "))
note("Packs: ", paste(names(table(ind$Pack)), table(ind$Pack), collapse = ", "), "; without a pack: ", sum(is.na(ind$Pack)))

# Three partitions to compare with Iosif et al.: genetics only, space only, complete workflow
fam <- setNames(as.character(gen$group), gen$Individual)
fam[table(fam)[fam] < 2] <- NA                                            # singletons don't make a group
spa <- spatial_groups(samples_Carpathians, ud = ud)
spa <- setNames(spa$Pack, spa$Individual)
full <- setNames(ind$Pack, ind$Individual)
cmp <- lapply(list(genetics = fam, space = spa, workflow = full), compare_partitions, reference = ref)
tab_cmp <- do.call(rbind, lapply(names(cmp), function(k) cbind(partition = k, cmp[[k]]$summary)))
print(tab_cmp)
for (k in names(cmp)) note("Agreement with Iosif, ", k, ": ARI = ", cmp[[k]]$summary$ARI, ", NMI = ", cmp[[k]]$summary$NMI,
                           ", correctly placed = ", cmp[[k]]$summary$n_correct, "/41, packs recovered = ", cmp[[k]]$summary$reference_groups_found)
print(cmp$workflow$table)

# The genetic settings of the first version of the article (Bruvo similarity, connected components) with the new engine
data("Bruvo_Carpathians")
old <- wolf_packs(samples_Carpathians, Bruvo_Carpathians, "1-Bruvo", 0.7, method = "components", period = "Campaign", ud = ud)
cmp_old <- compare_partitions(old$individuals, ref)$summary
note("Genetic settings of the first version + kernel engine: ARI = ", cmp_old$ARI, ", correctly placed = ", cmp_old$n_correct, "/41")


# ---- 3. The packs (Table 2) ----------------------------------------------------------------------------------
smry <- summarise_packs(res$samples, "Pack", ud = ud, relate = rel, estimator = "QG", n_perm = 999,
                        period = "Campaign", date = "Date")
maj <- apply(cmp$workflow$table[, smry$Pack, drop = FALSE], 2, function(x) names(which.max(x)))
tab2 <- data.frame(Pack = smry$Pack, Iosif = paste0(maj, " (", ref$PackName[match(maj, ref$Pack)], ")"),
                   Individuals = paste0(smry$n_individuals, " (", smry$n_males, " M, ", smry$n_females, " F)"),
                   Samples = smry$n_samples, Area95 = smry$area95_km2, Area50 = smry$area50_km2, Shared = smry$shared_pct,
                   Seasons = paste0(smry$first_period, " to ", smry$last_period), Persistent = smry$n_persistent,
                   Affinity = smry$min_affinity, Shift = smry$max_shift_km, Tenure = smry$tenure_years,
                   r = smry$mean_r, p = smry$r_p)
print(tab2)
write.csv(tab2, file.path(out, "Table2_packs.csv"), row.names = FALSE)
tl <- pack_timeline(res$samples, "Pack", period = "Campaign", ud = ud)
note("Individuals changing pack: ", paste(tl$individuals$Individual[tl$individuals$pattern == "changed pack"], collapse = ", "))
write.csv(attr(smry, "periods"), file.path(out, "TableS2_pack_periods.csv"), row.names = FALSE)


# ---- 4. Robustness: grid of settings, bootstrap, subsampling --------------------------------------------------
# The grid is 144 runs (about 40 s)
grid_runs <- pack_stability(samples_Carpathians, rel, "QG", 0.25, period = "Campaign", reference = ref,
                            vary = list(threshold = c(0.2, 0.25, 0.3), h = c(1200, 1500, round(ud$h), 2500),
                                        min_samples = c(1, 2, 3), method = c("leiden", "louvain", "walktrap", "components")))
g <- grid_runs$runs[grid_runs$runs$type == "sensitivity", ]
g$method <- sub(".*method = ", "", g$setting)
note("Grid (", nrow(g), " settings): median ARI with Iosif ", median(g$ARI_reference), " (", min(g$ARI_reference), "-", max(g$ARI_reference), ")")
for (m in unique(g$method)) note("  ", m, ": median ARI ", median(g$ARI_reference[g$method == m]),
                                  " (", min(g$ARI_reference[g$method == m]), "-", max(g$ARI_reference[g$method == m]), ")")
e_ids <- ref$Individual[ref$Pack == "e"]
co <- grid_runs$coassignment[e_ids, e_ids]
note("Pack e: share of runs where its members are together = ", round(mean(co[upper.tri(co)]), 2))
write.csv(g, file.path(out, "TableS1_sensitivity_grid.csv"), row.names = FALSE)

# Bootstrap (50 runs) and subsampling (20 runs per share kept), same seed
boot <- pack_stability(samples_Carpathians, rel, "QG", 0.25, period = "Campaign", reference = ref, n_boot = 50,
                       fraction = c(0.9, 0.75, 0.5, 0.25), n_sub = 20)
b <- boot$runs[boot$runs$type == "bootstrap", ]
note("Bootstrap: median ARI with Iosif ", median(b$ARI_reference), ", 95% interval ", paste(round(quantile(b$ARI_reference, c(0.025, 0.975)), 2), collapse = "-"))
print(boot$subsampling)
for (k in seq_len(nrow(boot$subsampling))) with(boot$subsampling[k, ], note(
  "Subsampling ", fraction * 100, "% (", n_samples, " samples): ARI with Iosif ", ARI_reference, " (", ARI_reference_low, "-",
  ARI_reference_high, "), individuals ", n_individuals, ", packs ", n_packs))
stab <- boot$individuals
stab$n <- as.integer(n_ind[stab$Individual])
stab$n_class <- cut(stab$n, c(0, 1, 2, 4, 9, Inf), labels = c("1", "2", "3-4", "5-9", "10+"))
st_cl <- aggregate(stability ~ n_class, stab[!is.na(stab$Pack), ], mean)
note("Mean stability (bootstrap + subsampling) by number of samples: ", paste(st_cl$n_class, round(st_cl$stability, 2), sep = " -> ", collapse = "; "))


# ---- 5. Fig. 1: the workflow ----------------------------------------------------------------------------------
# Boxes and arrows in ggplot2: coordinates set by hand, easy to move (x, y = centre; w, h = width, height)
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
  box(3.85, 1.35, 2.3, 1.25, "plot_packs()", "interactive map,\n95% and 50%\nisopleths", "Output"),
  box(6.3, 1.35, 2.4, 1.25, "pack_timeline()", "individuals and\nterritories, season\nby season", "Output"),
  box(9.6, 1.35, 3.6, 1.25, "compare_partitions()\npack_stability()", "agreement with a reference;\nsensitivity, bootstrap, subsampling", "Assessment"),
  box(10.6, 4.95, 2.6, 1.95, "Statuses", "Resident\nPutative disperser\nDisperser (dated)\nUnrelated resident\nTransient\nUndetermined", "Status"),
  box(10.6, 8.2, 2.6, 1.25, "Building blocks", "ud_overlap()\nisopleth()\nattribute_samples()", "Auxiliary")
)
arrows <- data.frame(x = c(2.3, 2.3, 7.0, 2.3, 7.0, 6.6, 4.65, 4.65, 4.65, 4.65, 4.65),
                     y = c(9.05, 7.75, 9.05, 5.975, 5.975, 4.95, 4.375, 2.925, 2.925, 2.925, 2.925),
                     xend = c(2.3, 2.3, 7.0, 3.6, 5.7, 9.3, 4.65, 1.35, 3.85, 6.3, 9.6),
                     yend = c(8.65, 7.125, 7.125, 5.525, 5.525, 4.95, 3.875, 1.975, 1.975, 1.975, 1.975))
fills <- c(Input = "#B8E3B0", Genetics = "#D9EAD3", Main = "#FFF2B3", Output = "#CFE2F3", Assessment = "#E4D4F0",
           Status = "#FFFFFF", Auxiliary = "#F9CB9C")
fig1 <- ggplot() +
  annotate("rect", xmin = 0.25, xmax = 9.05, ymin = 2.7, ymax = 7.4, fill = NA, colour = "grey45", linetype = "dashed") +
  annotate("text", x = 8.95, y = 2.85, label = "wolf_packs(): steps 1 to 4 in one call", hjust = 1, size = 2.6, colour = "grey30", fontface = "italic") +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), arrow = arrow(length = unit(1.6, "mm"), type = "closed"), colour = "grey35", linewidth = 0.4) +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = type), colour = "grey40", linewidth = 0.3) +
  geom_text(data = boxes, aes(x, y + h / 2 - 0.12, label = title), vjust = 1, size = 2.8, fontface = "bold", family = "mono", lineheight = 0.9) +
  geom_text(data = boxes, aes(x, y + h / 2 - ifelse(grepl("\n", title), 0.62, 0.38), label = body), vjust = 1, size = 2.35, lineheight = 0.95) +
  scale_fill_manual(values = fills, breaks = c("Input", "Main", "Output", "Assessment", "Auxiliary"),
                    labels = c("Input data", "Workflow", "Description", "Assessment", "Building blocks"), name = NULL) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0.6, 9.95), expand = FALSE) +
  theme_void(base_size = 8) + theme(legend.position = "bottom", plot.margin = margin(4, 4, 4, 4))
ggsave(file.path(out, "Fig1_workflow.png"), fig1, width = width_cm, height = 13.5, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(out, "Fig1_workflow.pdf"), fig1, width = width_cm, height = 13.5, units = "cm")


# ---- 6. Fig. 2: maps ------------------------------------------------------------------------------------------
pts <- res$samples                                                       # samples with Pack, Status...
pts$ref <- ref_p[as.character(pts$Individual)]
group_colours <- function(part) {                                        # each group takes the colour of its majority pack in Iosif et al.
  tb <- table(part[names(ref_p)], ref_p)
  tb <- tb[order(-apply(tb, 1, max)), , drop = FALSE]
  cols <- c(); free <- col_ref; k <- 0
  for (g in rownames(tb)) {
    best <- colnames(tb)[order(-tb[g, ])]; best <- best[best %in% names(free) & tb[g, best] > 0]
    if (length(best)) { cols[g] <- free[best[1]]; free <- free[names(free) != best[1]] } else { k <- k + 1; cols[g] <- col_extra[k] }
  }
  for (g in setdiff(unique(na.omit(part)), names(cols))) { k <- k + 1; cols[g] <- col_extra[k] }
  cols
}
outlines <- function(groups, level) {                                    # isopleths of the groups, in WGS84
  do.call(rbind, lapply(names(groups), function(g) { p <- isopleth(ud, groups[[g]], level); p$group <- g; p })) |> st_transform(4326)
}
map_theme <- list(theme_bw(base_size = 7), coord_sf(expand = FALSE, xlim = c(24.94, 25.36), ylim = c(45.21, 45.78)),
                  scale_x_continuous(breaks = c(25, 25.15, 25.3)), labs(x = NULL, y = NULL),
                  guides(colour = guide_legend(ncol = 2, title.position = "top"), shape = guide_legend(ncol = 2, title.position = "top")),
                  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 5.5),
                        legend.title = element_text(size = 6), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0),
                        axis.text = element_text(size = 5), plot.title = element_text(size = 7.5), panel.grid = element_line(linewidth = 0.2)))
pack_map <- function(part, title, legend_title, iso = TRUE, core = FALSE) {
  cols <- group_colours(part)
  d <- pts; d$g <- part[as.character(d$Individual)]; d$g[is.na(d$g)] <- "none"
  vals <- c(cols, none = "grey75")
  p <- ggplot()
  if (iso) {
    grp <- split(names(part)[!is.na(part)], part[!is.na(part)])
    p <- p + geom_sf(data = outlines(grp, 0.95), aes(colour = group), fill = NA, linewidth = 0.45, show.legend = FALSE)
    if (core) p <- p + geom_sf(data = outlines(grp, 0.5), aes(fill = group), colour = NA, alpha = 0.35, show.legend = FALSE)
  }
  p + geom_sf(data = d[d$g == "none", ], colour = "grey60", size = 0.5) +
    geom_sf(data = d[d$g != "none", ], aes(colour = g), size = 0.7) +
    scale_colour_manual(values = vals, breaks = sort(names(cols)), name = legend_title) +
    scale_fill_manual(values = vals, guide = "none") + ggtitle(title) + map_theme
}
eff <- effort_Carpathians
eff$Genotyped <- ifelse(eff$Verdict == "Genotyped", "wolf genotyped", "other")
pa <- ggplot() + geom_sf(data = eff[eff$Genotyped == "other", ], aes(colour = Genotyped), size = 0.5) +
  geom_sf(data = eff[eff$Genotyped != "other", ], aes(colour = Genotyped), size = 0.5) +
  scale_colour_manual(values = c("wolf genotyped" = "black", other = "grey70"), name = "Samples") + ggtitle("a. Sampling") + map_theme
pb <- suppressMessages(pack_map(ref_p, "b. Iosif et al. (2025)", "Pack", core = TRUE) +     # pack names are in the figure caption
  scale_colour_manual(values = c(col_ref, none = "grey75"), breaks = names(col_ref), name = "Pack (Iosif et al.)") +
  guides(colour = guide_legend(ncol = 3, title.position = "top")))
pc <- pack_map(fam, "c. Genetics only", "Family", iso = FALSE)
pd <- pack_map(spa, "d. Space only", "Pack", core = TRUE)
pe <- pack_map(full, "e. Complete workflow", "Pack", core = TRUE)

# f. agreement with Iosif et al.: each individual of the reference at the centre of its samples
maj_ref <- apply(cmp$workflow$table, 2, function(x) names(which.max(x)))
verdict <- ifelse(is.na(full[names(ref_p)]), "no pack", ifelse(maj_ref[full[names(ref_p)]] == ref_p, "same pack", "other pack"))
cen <- aggregate(st_coordinates(pts)[, 1:2], list(Individual = as.character(pts$Individual)), mean)
cen <- st_as_sf(cen[cen$Individual %in% names(ref_p), ], coords = c("X", "Y"), crs = 4326)
cen$verdict <- factor(verdict[cen$Individual], levels = c("same pack", "other pack", "no pack"))
cen$ref <- ref_p[cen$Individual]
note("Reference individuals misplaced or without a pack: ", sum(verdict != "same pack"), ", of which ", sum(verdict != "same pack" & n_ind[names(verdict)] <= 2),
     " with 1 or 2 samples (", paste(names(verdict)[verdict != "same pack"], collapse = ", "), ")")
grp <- split(names(full)[!is.na(full)], full[!is.na(full)])
pf <- ggplot() + geom_sf(data = outlines(grp, 0.95), colour = "grey55", fill = NA, linewidth = 0.35) +
  geom_sf(data = cen, aes(shape = verdict, colour = verdict), size = 1.3, stroke = 0.6) +
  geom_sf_text(data = cen[cen$verdict != "same pack", ], aes(label = ref), size = 2, nudge_x = 0.012, colour = "grey20") +
  scale_shape_manual(values = c("same pack" = 16, "other pack" = 4, "no pack" = 1), name = "vs Iosif") +
  scale_colour_manual(values = c("same pack" = "#009E73", "other pack" = "#D55E00", "no pack" = "black"), name = "vs Iosif") +
  ggtitle("f. Agreement with Iosif et al.") + map_theme
# Shared legend of the isopleths: outline = territory (95%), shading = core area (50%)
leg_iso <- ggplot() +
  annotate("rect", xmin = 0, xmax = 0.7, ymin = 0, ymax = 1, fill = NA, colour = "grey35", linewidth = 0.45) +
  annotate("text", x = 0.9, y = 0.5, label = "95% kernel isopleth (territory)", hjust = 0, size = 2.4) +
  annotate("rect", xmin = 6.5, xmax = 7.2, ymin = 0, ymax = 1, fill = "grey35", alpha = 0.35, colour = NA) +
  annotate("text", x = 7.4, y = 0.5, label = "50% kernel isopleth (core area)", hjust = 0, size = 2.4) +
  coord_cartesian(xlim = c(-2.5, 14), ylim = c(-0.4, 1.4), expand = FALSE) + theme_void()
fig2 <- (pa | pb | pc) / (pd | pe | pf) / leg_iso + plot_layout(heights = c(1, 1, 0.045))
ggsave(file.path(out, "Fig2_maps.png"), fig2, width = width_cm, height = 22, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(out, "Fig2_maps.pdf"), fig2, width = width_cm, height = 22, units = "cm")


# ---- 7. Fig. 3: robustness ------------------------------------------------------------------------------------
set.seed(1)                                                              # so that jittered points always land in the same place
base_ari <- cmp$workflow$summary$ARI
g$method <- factor(g$method, levels = c("leiden", "louvain", "walktrap", "components"))
p3a <- ggplot(g, aes(method, ARI_reference)) + geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey92", linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, size = 0.6, alpha = 0.6) + geom_hline(yintercept = base_ari, linetype = "dashed", colour = "#D55E00") +
  labs(x = "Community detection", y = "ARI with Iosif et al. (2025)", title = "a. Settings") + ylim(0, 1) + theme_bw(base_size = 7)
ss <- boot$subsampling
ss <- rbind(data.frame(fraction = 1, n_samples = nrow(samples_Carpathians), ARI_reference = base_ari, ARI_reference_low = base_ari,
                       ARI_reference_high = base_ari, ARI_base = 1, ARI_base_low = 1, ARI_base_high = 1),
            ss[, c("fraction", "n_samples", "ARI_reference", "ARI_reference_low", "ARI_reference_high", "ARI_base", "ARI_base_low", "ARI_base_high")])
ssl <- rbind(data.frame(fraction = ss$fraction, med = ss$ARI_reference, lo = ss$ARI_reference_low, hi = ss$ARI_reference_high, vs = "Iosif et al. (2025)"),
             data.frame(fraction = ss$fraction, med = ss$ARI_base, lo = ss$ARI_base_low, hi = ss$ARI_base_high, vs = "all samples"))
p3b <- ggplot(ssl, aes(100 * fraction, med, colour = vs)) + geom_ribbon(aes(ymin = lo, ymax = hi, fill = vs), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.5) + geom_point(size = 1) + scale_colour_manual(values = c("#0072B2", "grey30"), name = "ARI with") +
  scale_fill_manual(values = c("#0072B2", "grey30"), name = "ARI with") +
  labs(x = "Samples kept (%)", y = "ARI", title = "b. Subsampling") + ylim(0, 1) + theme_bw(base_size = 7) +
  legend_inside(0.7, 0.2) + theme(legend.key.size = unit(3, "mm"), legend.background = element_blank())
nb <- data.frame(Individual = names(verdict), verdict = factor(verdict, levels = c("same pack", "other pack", "no pack")))
nb$n_class <- cut(as.integer(n_ind[nb$Individual]), c(0, 1, 2, 4, 9, Inf), labels = c("1", "2", "3-4", "5-9", "10+"))
p3c <- ggplot(nb, aes(n_class, fill = verdict)) + geom_bar(width = 0.65, colour = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = c("same pack" = "#009E73", "other pack" = "#D55E00", "no pack" = "grey85"), name = "vs Iosif") +
  labs(x = "Samples per individual", y = "Individuals of the reference", title = "c. Agreement by sampling") +
  theme_bw(base_size = 7) + legend_inside(0.27, 0.8) + theme(legend.key.size = unit(3, "mm"), legend.background = element_blank())
fig3 <- p3a | p3b | p3c
ggsave(file.path(out, "Fig3_robustness.png"), fig3, width = width_cm, height = 6.5, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(out, "Fig3_robustness.pdf"), fig3, width = width_cm, height = 6.5, units = "cm")

writeLines(numbers, file.path(out, "numbers_in_text.txt"))
sessionInfo()

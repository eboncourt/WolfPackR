# Figures, tableaux et chiffres de l'article WolfPackR (révision ESE) : c'est l'Appendix S1
# Il faut WolfPackR 1.1.0 (rayon de giration sans les 10 % d'échantillons les plus éloignés). Tout est recalculé à partir des données du package ; compter 2 à 3 minutes en tout.
# Fig. 1 : le workflow ; Fig. 2 : cartes comparées à Iosif et al. (2025) ; Fig. 3 : robustesse ;
# Tableau 2 : description des meutes ; Tableaux S1-S2 : annexes ; chiffres_article.txt : les nombres cités dans le texte

library(WolfPackR)
library(sf)
library(ggplot2)
library(patchwork)

out <- "article/output"                                  # dossier de sortie, relatif à la racine du dépôt
dir.create(out, showWarnings = FALSE, recursive = TRUE)
largeur <- 17                                            # cm, pleine largeur de page

data("samples_Carpathians"); data("effort_Carpathians"); data("relate_Carpathians"); data("reference_Carpathians")
rel <- relate_Carpathians
ref <- reference_Carpathians
ref_p <- setNames(ref$Pack, ref$Individual)
chiffres <- character(0)
note <- function(...) { txt <- paste0(...); chiffres <<- c(chiffres, txt); message(txt) }   # garde une trace des nombres cités

# Couleurs des meutes de Iosif (Okabe-Ito, lisibles par les daltoniens) ; les autres groupes prennent des couleurs en plus
col_ref <- c(a = "#E69F00", b = "#56B4E9", c = "#009E73", d = "#CC79A7", e = "#D55E00", f = "#0072B2")
col_extra <- c("#8C6D31", "#6A3D9A", "#B15928", "#33A02C", "#E31A1C")
dans <- function(x, y) if (packageVersion("ggplot2") >= "3.5.0")             # légende dans le panneau, selon la version de ggplot2
  theme(legend.position = "inside", legend.position.inside = c(x, y)) else theme(legend.position = c(x, y))


# ---- 1. Les données ----------------------------------------------------------------------------------------
note("Echantillons collectés avec position : ", nrow(effort_Carpathians), " ; génotypés loup : ", sum(effort_Carpathians$Verdict == "Genotyped"))
note("Echantillons de loups : ", nrow(samples_Carpathians), " ; individus : ", length(unique(samples_Carpathians$Individual)))
note("Par saison : ", paste(names(table(samples_Carpathians$Campaign)), table(samples_Carpathians$Campaign), collapse = ", "))
note("Par type : ", paste(names(table(samples_Carpathians$Type)), table(samples_Carpathians$Type), collapse = ", "))
n_ind <- table(samples_Carpathians$Individual)
note("Echantillons par individu : médiane ", median(n_ind), ", 1 échantillon : ", sum(n_ind == 1), ", 2 : ", sum(n_ind == 2), ", >= 10 : ", sum(n_ind >= 10))
lat <- st_coordinates(effort_Carpathians)[, 2]
note("Succès de génotypage au sud (< 45.45°N) : ", round(100 * mean(effort_Carpathians$Verdict[lat < 45.45] == "Genotyped")),
     " % ; au nord : ", round(100 * mean(effort_Carpathians$Verdict[lat >= 45.45] == "Genotyped")), " %")


# ---- 2. Le workflow pas à pas, puis en un appel ------------------------------------------------------------
gen <- genetic_groups(rel, "QG", 0.25, n_null = 199)
mod <- attr(gen, "modularity")
note("Groupes génétiques (Leiden, QG > 0.25) : tailles ", paste(sort(table(gen$group), decreasing = TRUE), collapse = " "),
     " ; modularité Q = ", round(mod["Q"], 2), ", z = ", round(mod["z"], 1))
ud <- space_use(samples_Carpathians)
note("Bande passante h = ", round(ud$h), " m")

res <- wolf_packs(samples_Carpathians, rel, "QG", 0.25, period = "Campaign", ud = ud)
ind <- res$individuals
note("Statuts : ", paste(names(table(ind$Status)), table(ind$Status), collapse = ", "))
note("Meutes : ", paste(names(table(ind$Pack)), table(ind$Pack), collapse = ", "), " ; sans meute : ", sum(is.na(ind$Pack)))

# Trois partitions à comparer à Iosif : génétique seule, espace seul, workflow complet
fam <- setNames(as.character(gen$group), gen$Individual)
fam[table(fam)[fam] < 2] <- NA                                            # les singletons ne forment pas de groupe
spa <- spatial_groups(samples_Carpathians, ud = ud)
spa <- setNames(spa$Pack, spa$Individual)
full <- setNames(ind$Pack, ind$Individual)
cmp <- lapply(list(genetics = fam, space = spa, workflow = full), compare_partitions, reference = ref)
tab_cmp <- do.call(rbind, lapply(names(cmp), function(k) cbind(partition = k, cmp[[k]]$summary)))
print(tab_cmp)
for (k in names(cmp)) note("Comparaison à Iosif, ", k, " : ARI = ", cmp[[k]]$summary$ARI, ", NMI = ", cmp[[k]]$summary$NMI,
                           ", bien placés = ", cmp[[k]]$summary$n_correct, "/41, meutes retrouvées = ", cmp[[k]]$summary$reference_groups_found)
print(cmp$workflow$table)

# Les réglages génétiques de la première version de l'article (similarité de Bruvo, composantes connexes) avec le nouveau moteur
data("Bruvo_Carpathians")
old <- wolf_packs(samples_Carpathians, Bruvo_Carpathians, "1-Bruvo", 0.7, method = "components", period = "Campaign", ud = ud)
cmp_old <- compare_partitions(old$individuals, ref)$summary
note("Réglages génétiques de la 1re version + moteur KDE : ARI = ", cmp_old$ARI, ", bien placés = ", cmp_old$n_correct, "/41")


# ---- 3. Description des meutes (Tableau 2) -------------------------------------------------------------------
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
note("Individus qui changent de meute : ", paste(tl$individuals$Individual[tl$individuals$pattern == "changed pack"], collapse = ", "))
write.csv(attr(smry, "periods"), file.path(out, "TableS2_pack_periods.csv"), row.names = FALSE)


# ---- 4. Robustesse : grille de réglages, bootstrap, sous-échantillonnage -------------------------------------
# La grille fait 144 runs (environ 40 s)
grid_runs <- pack_stability(samples_Carpathians, rel, "QG", 0.25, period = "Campaign", reference = ref,
                            vary = list(threshold = c(0.2, 0.25, 0.3), h = c(1200, 1500, round(ud$h), 2500),
                                        min_samples = c(1, 2, 3), method = c("leiden", "louvain", "walktrap", "components")))
g <- grid_runs$runs[grid_runs$runs$type == "sensitivity", ]
g$method <- sub(".*method = ", "", g$setting)
note("Grille (", nrow(g), " réglages) : ARI avec Iosif médian ", median(g$ARI_reference), " (", min(g$ARI_reference), "-", max(g$ARI_reference), ")")
for (m in unique(g$method)) note("  ", m, " : ARI médian ", median(g$ARI_reference[g$method == m]),
                                  " (", min(g$ARI_reference[g$method == m]), "-", max(g$ARI_reference[g$method == m]), ")")
e_ids <- ref$Individual[ref$Pack == "e"]
co <- grid_runs$coassignment[e_ids, e_ids]
note("Meute e : part des runs où ses membres sont ensemble = ", round(mean(co[upper.tri(co)]), 2))
write.csv(g, file.path(out, "TableS1_sensitivity_grid.csv"), row.names = FALSE)

# Bootstrap (50 runs) et sous-échantillonnage (20 runs par part gardée), avec la même graine
boot <- pack_stability(samples_Carpathians, rel, "QG", 0.25, period = "Campaign", reference = ref, n_boot = 50,
                       fraction = c(0.9, 0.75, 0.5, 0.25), n_sub = 20)
b <- boot$runs[boot$runs$type == "bootstrap", ]
note("Bootstrap : ARI avec Iosif médian ", median(b$ARI_reference), ", IC 95 % ", paste(round(quantile(b$ARI_reference, c(0.025, 0.975)), 2), collapse = "-"))
print(boot$subsampling)
for (k in seq_len(nrow(boot$subsampling))) with(boot$subsampling[k, ], note(
  "Sous-échantillonnage ", fraction * 100, " % (", n_samples, " échantillons) : ARI Iosif ", ARI_reference, " (", ARI_reference_low, "-",
  ARI_reference_high, "), individus ", n_individuals, ", meutes ", n_packs))
stab <- boot$individuals
stab$n <- as.integer(n_ind[stab$Individual])
stab$classe <- cut(stab$n, c(0, 1, 2, 4, 9, Inf), labels = c("1", "2", "3-4", "5-9", "10+"))
st_cl <- aggregate(stability ~ classe, stab[!is.na(stab$Pack), ], mean)
note("Stabilité moyenne (bootstrap + sous-éch.) par nombre d'échantillons : ", paste(st_cl$classe, round(st_cl$stability, 2), sep = " -> ", collapse = " ; "))


# ---- 5. Fig. 1 : le workflow ------------------------------------------------------------------------------------
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


# ---- 6. Fig. 2 : cartes ------------------------------------------------------------------------------------------
pts <- res$samples                                                       # échantillons avec Pack, Status...
pts$ref <- ref_p[as.character(pts$Individual)]
couleurs <- function(part) {                                             # chaque groupe prend la couleur de la meute de Iosif majoritaire
  tb <- table(part[names(ref_p)], ref_p)
  tb <- tb[order(-apply(tb, 1, max)), , drop = FALSE]
  cols <- c(); libres <- col_ref; k <- 0
  for (g in rownames(tb)) {
    best <- colnames(tb)[order(-tb[g, ])]; best <- best[best %in% names(libres) & tb[g, best] > 0]
    if (length(best)) { cols[g] <- libres[best[1]]; libres <- libres[names(libres) != best[1]] } else { k <- k + 1; cols[g] <- col_extra[k] }
  }
  for (g in setdiff(unique(na.omit(part)), names(cols))) { k <- k + 1; cols[g] <- col_extra[k] }
  cols
}
contours <- function(groups, level) {                                   # isoplèthes des groupes, en WGS84
  do.call(rbind, lapply(names(groups), function(g) { p <- isopleth(ud, groups[[g]], level); p$group <- g; p })) |> st_transform(4326)
}
fond <- list(theme_bw(base_size = 7), coord_sf(expand = FALSE, xlim = c(24.94, 25.36), ylim = c(45.21, 45.78)),
             scale_x_continuous(breaks = c(25, 25.15, 25.3)), labs(x = NULL, y = NULL),
             guides(colour = guide_legend(ncol = 2, title.position = "top"), shape = guide_legend(ncol = 2, title.position = "top")),
             theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 5.5),
                   legend.title = element_text(size = 6), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-4, 0, 0, 0),
                   axis.text = element_text(size = 5), plot.title = element_text(size = 7.5), panel.grid = element_line(linewidth = 0.2)))
carte <- function(part, titre, legende, iso = TRUE, core = FALSE) {
  cols <- couleurs(part)
  d <- pts; d$g <- part[as.character(d$Individual)]; d$g[is.na(d$g)] <- "none"
  vals <- c(cols, none = "grey75")
  p <- ggplot()
  if (iso) {
    grp <- split(names(part)[!is.na(part)], part[!is.na(part)])
    p <- p + geom_sf(data = contours(grp, 0.95), aes(colour = group), fill = NA, linewidth = 0.45, show.legend = FALSE)
    if (core) p <- p + geom_sf(data = contours(grp, 0.5), aes(fill = group), colour = NA, alpha = 0.35, show.legend = FALSE)
  }
  p + geom_sf(data = d[d$g == "none", ], colour = "grey60", size = 0.5) +
    geom_sf(data = d[d$g != "none", ], aes(colour = g), size = 0.7) +
    scale_colour_manual(values = vals, breaks = sort(names(cols)), name = legende) +
    scale_fill_manual(values = vals, guide = "none") + ggtitle(titre) + fond
}
eff <- effort_Carpathians
eff$Genotyped <- ifelse(eff$Verdict == "Genotyped", "wolf genotyped", "other")
pa <- ggplot() + geom_sf(data = eff[eff$Genotyped == "other", ], aes(colour = Genotyped), size = 0.5) +
  geom_sf(data = eff[eff$Genotyped != "other", ], aes(colour = Genotyped), size = 0.5) +
  scale_colour_manual(values = c("wolf genotyped" = "black", other = "grey70"), name = "Samples") + ggtitle("a. Sampling") + fond
pb <- suppressMessages(carte(ref_p, "b. Iosif et al. (2025)", "Pack", core = TRUE) +              # noms des meutes dans la légende de la figure
  scale_colour_manual(values = c(col_ref, none = "grey75"), breaks = names(col_ref), name = "Pack (Iosif et al.)") +
  guides(colour = guide_legend(ncol = 3, title.position = "top")))
pc <- carte(fam, "c. Genetics only", "Family", iso = FALSE)
pd <- carte(spa, "d. Space only", "Pack", core = TRUE)
pe <- carte(full, "e. Complete workflow", "Pack", core = TRUE)

# f. accord avec Iosif : chaque individu de la référence au centre de ses échantillons
okmaj <- apply(cmp$workflow$table, 2, function(x) names(which.max(x)))
verdict <- ifelse(is.na(full[names(ref_p)]), "no pack", ifelse(okmaj[full[names(ref_p)]] == ref_p, "same pack", "other pack"))
cen <- aggregate(st_coordinates(pts)[, 1:2], list(Individual = as.character(pts$Individual)), mean)
cen <- st_as_sf(cen[cen$Individual %in% names(ref_p), ], coords = c("X", "Y"), crs = 4326)
cen$verdict <- factor(verdict[cen$Individual], levels = c("same pack", "other pack", "no pack"))
cen$ref <- ref_p[cen$Individual]
note("Individus de la référence mal placés ou sans meute : ", sum(verdict != "same pack"), " dont ", sum(verdict != "same pack" & n_ind[names(verdict)] <= 2),
     " avec 1 ou 2 échantillons (", paste(names(verdict)[verdict != "same pack"], collapse = ", "), ")")
grp <- split(names(full)[!is.na(full)], full[!is.na(full)])
pf <- ggplot() + geom_sf(data = contours(grp, 0.95), colour = "grey55", fill = NA, linewidth = 0.35) +
  geom_sf(data = cen, aes(shape = verdict, colour = verdict), size = 1.3, stroke = 0.6) +
  geom_sf_text(data = cen[cen$verdict != "same pack", ], aes(label = ref), size = 2, nudge_x = 0.012, colour = "grey20") +
  scale_shape_manual(values = c("same pack" = 16, "other pack" = 4, "no pack" = 1), name = "vs Iosif") +
  scale_colour_manual(values = c("same pack" = "#009E73", "other pack" = "#D55E00", "no pack" = "black"), name = "vs Iosif") +
  ggtitle("f. Agreement with Iosif et al.") + fond
# Légende commune des isoplèthes : contour = territoire (95 %), aplat = cœur de territoire (50 %)
leg_iso <- ggplot() +
  annotate("rect", xmin = 0, xmax = 0.7, ymin = 0, ymax = 1, fill = NA, colour = "grey35", linewidth = 0.45) +
  annotate("text", x = 0.9, y = 0.5, label = "95% kernel isopleth (territory)", hjust = 0, size = 2.4) +
  annotate("rect", xmin = 6.5, xmax = 7.2, ymin = 0, ymax = 1, fill = "grey35", alpha = 0.35, colour = NA) +
  annotate("text", x = 7.4, y = 0.5, label = "50% kernel isopleth (core area)", hjust = 0, size = 2.4) +
  coord_cartesian(xlim = c(-2.5, 14), ylim = c(-0.4, 1.4), expand = FALSE) + theme_void()
fig2 <- (pa | pb | pc) / (pd | pe | pf) / leg_iso + plot_layout(heights = c(1, 1, 0.045))
ggsave(file.path(out, "Fig2_maps.png"), fig2, width = largeur, height = 22, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(out, "Fig2_maps.pdf"), fig2, width = largeur, height = 22, units = "cm")


# ---- 7. Fig. 3 : robustesse --------------------------------------------------------------------------------------
set.seed(1)                                                              # pour que les points décalés (jitter) soient toujours au même endroit
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
  dans(0.7, 0.2) + theme(legend.key.size = unit(3, "mm"), legend.background = element_blank())
nb <- data.frame(Individual = names(verdict), verdict = factor(verdict, levels = c("same pack", "other pack", "no pack")))
nb$classe <- cut(as.integer(n_ind[nb$Individual]), c(0, 1, 2, 4, 9, Inf), labels = c("1", "2", "3-4", "5-9", "10+"))
p3c <- ggplot(nb, aes(classe, fill = verdict)) + geom_bar(width = 0.65, colour = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = c("same pack" = "#009E73", "other pack" = "#D55E00", "no pack" = "grey85"), name = "vs Iosif") +
  labs(x = "Samples per individual", y = "Individuals of the reference", title = "c. Agreement by sampling") +
  theme_bw(base_size = 7) + dans(0.27, 0.8) + theme(legend.key.size = unit(3, "mm"), legend.background = element_blank())
fig3 <- p3a | p3b | p3c
ggsave(file.path(out, "Fig3_robustness.png"), fig3, width = largeur, height = 6.5, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(out, "Fig3_robustness.pdf"), fig3, width = largeur, height = 6.5, units = "cm")

writeLines(chiffres, file.path(out, "chiffres_article.txt"))
sessionInfo()

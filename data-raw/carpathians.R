# Builds the Carpathian datasets of the package from the Zenodo file of Foundation Conservation Carpathia
# (https://doi.org/10.5281/zenodo.14544266). Put WolfSampleDataFCCSep2021Nandor.xlsx in data-raw/ and run from the package root.
# The xlsx itself is not shipped: everything useful from it ends up in the five datasets below.
library(sf)
library(readxl)
pkgload::load_all()                                     # for relatedness()

xlsx <- "data-raw/WolfSampleDataFCCSep2021Nandor.xlsx"
rout <- as.data.frame(read_excel(xlsx, sheet = "ROut"))                    # one row per collected sample (505)
rout <- rout[!(rout$X == 0 & rout$Y == 0), ]                               # 3 samples with no position (stored as 0, 0): a dog and 2 failed urines
meta <- as.data.frame(read_excel(xlsx, sheet = "SampleData"))              # campaign label
strc <- as.data.frame(read_excel(xlsx, sheet = "StructureRes2021_2"))      # wolf ancestry of each genotyped animal
ref  <- as.data.frame(read_excel(xlsx, sheet = "ReferenceGenotypes"))      # reference genotype of each animal

# --- Who is a dog? Animals put in the dog cluster by Structure (p2Wolf < 0.5).
# Catches CC01CK: its reference sample is a dog but three of its samples are labelled "Wolf" in ROut
p_wolf <- setNames(strc$p2Wolf, strc$Sample)
dogs <- names(p_wolf)[p_wolf < 0.5]
rout$Verdict[rout$AnimalRef %in% dogs] <- "NonTarget"
rout$Species[rout$AnimalRef %in% dogs] <- "Dog"

# --- Sampling season, July to June (2017-18, 2018-19, 2019-20), from the date; from the Zenodo campaign when there's no date
season_of <- function(date, sample) {
  d <- as.Date(date)
  y <- as.integer(format(d, "%Y")) - (as.integer(format(d, "%m")) < 7)
  out <- ifelse(is.na(y), NA, paste0(y, "-", substr(y + 1, 3, 4)))
  study <- meta$Study[match(sample, meta$Sample)]
  out[is.na(out) & grepl("19_20", study)] <- "2019-20"
  out[is.na(out) & grepl("18_19", study)] <- "2018-19"
  out[is.na(out) & grepl("17_18", study)] <- "2017-18"
  out
}

# --- Every collected sample, genotyped or not: sampling effort and genotyping success
effort_Carpathians <- data.frame(
  Sample = rout$Sample,
  Date = as.Date(rout$Date),
  Campaign = season_of(rout$Date, rout$Sample),
  Type = rout$SType,
  Verdict = rout$Verdict,                               # Genotyped, Poor, Mixed or NonTarget
  Species = rout$Species,                               # Wolf, Hybrid, Dog, Fox, or Unknown when not genotyped
  Individual = ifelse(rout$Verdict == "Genotyped", rout$AnimalRef, NA),
  lon = rout$Y, lat = rout$X,                           # careful: in the Zenodo file X is the latitude
  stringsAsFactors = FALSE
)
effort_Carpathians <- st_as_sf(effort_Carpathians, coords = c("lon", "lat"), crs = 4326)
effort_Carpathians <- effort_Carpathians[order(effort_Carpathians$Date, effort_Carpathians$Sample), ]
rownames(effort_Carpathians) <- NULL

# --- Samples of wolves (and the F1 hybrid EF.15J8): the ones the workflow runs on
g <- rout[rout$Verdict %in% "Genotyped", ]
samples_Carpathians <- data.frame(
  Individual = g$AnimalRef,                             # keep it first: some functions expect the ID in column 1
  Sex = g$GeneticSex,
  Date = as.Date(g$Date),                               # 2 samples have no date (season known though)
  Campaign = season_of(g$Date, g$Sample),
  Sample = g$Sample,
  Type = g$SType,                                       # scat, urine, hair, tissue, blood
  pWolf = unname(p_wolf[g$AnimalRef]),                  # wolf ancestry (Structure): 0.51 for the hybrid, 0.90 for CC02H8
  lon = g$Y, lat = g$X,
  stringsAsFactors = FALSE
)
samples_Carpathians <- st_as_sf(samples_Carpathians, coords = c("lon", "lat"), crs = 4326)
samples_Carpathians <- samples_Carpathians[order(samples_Carpathians$Individual, samples_Carpathians$Date), ]
rownames(samples_Carpathians) <- NULL
wolves <- sort(unique(samples_Carpathians$Individual))
stopifnot(!anyNA(samples_Carpathians$pWolf), !any(wolves %in% dogs))

# --- Genotypes: reference genotype of each individual, 22 microsatellites, two columns per locus
ref <- ref[ref$AnimalRef %in% wolves & ref$IsAnimalReference %in% TRUE, ]
loci <- names(ref)[which(names(ref) == "CPH9"):which(names(ref) == "AHTh171")]
split_alleles <- function(v) {                           # "144 150" -> c(144, 150); "267271" (missing space) -> c(267, 271)
  v <- trimws(as.character(v))
  if (is.na(v) || v == "") return(c(NA, NA))
  if (grepl(" ", v)) return(as.integer(strsplit(v, " +")[[1]][1:2]))
  if (nchar(v) == 6) return(as.integer(c(substr(v, 1, 3), substr(v, 4, 6))))
  c(NA, NA)
}
genotypes_Carpathians <- data.frame(Individual = ref$AnimalRef, stringsAsFactors = FALSE)
for (l in loci) {
  a <- t(vapply(ref[[l]], split_alleles, integer(2)))
  genotypes_Carpathians[[paste0(l, "_1")]] <- a[, 1]
  genotypes_Carpathians[[paste0(l, "_2")]] <- a[, 2]
}
genotypes_Carpathians <- genotypes_Carpathians[order(genotypes_Carpathians$Individual), ]
rownames(genotypes_Carpathians) <- NULL
stopifnot(setequal(genotypes_Carpathians$Individual, wolves))

# --- Pairwise relatedness (QG, LR, LOD); allele frequencies from wolves only, so without the hybrid EF.15J8
relate_Carpathians <- relatedness(genotypes_Carpathians, freq_ids = setdiff(wolves, "EF.15J8"))

# --- 1 - Bruvo similarity used in the first version of the paper, now restricted to the wolves (dogs removed)
load("data/Bruvo_Carpathians.rda")
Bruvo_Carpathians <- Bruvo_Carpathians[Bruvo_Carpathians$ind1 %in% wolves & Bruvo_Carpathians$ind2 %in% wolves, ]
rownames(Bruvo_Carpathians) <- NULL

# --- Reference packs of Iosif et al. (2025), read by hand on their Fig. 3 (parentage with COLONY + movements)
# Only the 41 individuals involved in a parent-offspring link are on the figure; a few labels are doubtful (Note)
reference_Carpathians <- data.frame(
  Individual = c("CC01ET", "EK.18AU", "CC01EC", "CC01EE", "CC01EH",
                 "CC00T7", "CC00T6", "CC01F1", "CC0655", "CC0630", "CC064U", "CC04FK",
                 "CC00X0", "CC00TC", "CC01EU", "CC01CY", "CC0666", "CC04K4", "CC04AU", "CC04EK", "CC01CJ", "CC01FE", "CC01F2",
                 "CC025E", "CC024X", "CC01AY", "CC01EP", "CC03E3", "CC00JK", "CC00JM", "CC05X7",
                 "CC063M", "CC063J", "EF.15J8",
                 "CC04C0", "CC04JC", "CC03KH", "CC0637", "CC04H8", "CC066H", "CC05X0"),
  Pack = rep(c("a", "b", "c", "d", "e", "f"), times = c(5, 7, 11, 8, 3, 7)),
  stringsAsFactors = FALSE
)
reference_Carpathians$PackName <- c(a = "Sercaita", b = "Dambovita-Barsa", c = "Piatra Craiului",
                                    d = "Dambovita-Raul Targului", e = "Stoenesti hybrids", f = "Stoenesti")[reference_Carpathians$Pack]
reference_Carpathians$Breeder <- c("M", "", "", "", "",  "F", "F", "M", "", "", "", "",  "F", "M", "M", rep("", 8),
                                   "M", "F", rep("", 6),  "F", "", "",  "F", "M", rep("", 5))
reference_Carpathians$Note <- ""
reference_Carpathians$Note[reference_Carpathians$Individual == "CC01EU"] <- "labelled both as breeder and as offspring"
reference_Carpathians$Note[reference_Carpathians$Individual %in% c("CC0655", "CC01F2")] <- "on the edge of two territories"
reference_Carpathians$Note[reference_Carpathians$Individual == "CC05X0"] <- "outside the ellipses, between e and f"
stopifnot(all(reference_Carpathians$Individual %in% wolves))

# --- Save (version 2 so older R versions can still read them; no LazyData in this package)
for (nm in c("samples_Carpathians", "effort_Carpathians", "genotypes_Carpathians", "relate_Carpathians", "Bruvo_Carpathians",
             "reference_Carpathians"))
  save(list = nm, file = file.path("data", paste0(nm, ".rda")), compress = "xz", version = 2)

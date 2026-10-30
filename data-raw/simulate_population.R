# simulate_population(): a made-up wolf population followed for four biological years (July to June, 2020-21 to
# 2023-24), where we know everything: who belongs to which pack, who is whose parent, who dispersed where and when.
# Handy to see what WolfPackR gets right and wrong. The example datasets of the package are simulate_population(seed = 22)
# (see simulate_example.R); other seeds give other populations with the same story, e.g. to test settings or sampling.
# `foray` is the share of the samples of adults left in a neighbouring territory.
# Needs sf and WolfPackR (for relatedness()) loaded:
#   library(sf); library(WolfPackR); source("data-raw/simulate_population.R"); x <- simulate_population(seed = 5)
# It returns a list with samples_example, genotypes_example, relate_example and truth_example.
#
# The story (pack letters are those of truth_example; WolfPackR numbers its packs its own way):
#   * 8 packs (A to H) at the start, on a hexagonal lattice 15 km apart, each with a breeding pair and offspring;
#   * the breeding male of B and the breeding female of D are full sibs, so are the breeding female of A and the
#     breeding male of G: neighbouring packs are related, as in real populations;
#   * a daughter of A leaves in December 2021 and founds pack I next door with a male from outside (a family over two territories);
#   * both breeders of C die (Oct 2021, Feb 2022): the young stay a while, then one son becomes the breeding male of F
#     in summer 2023 (no pup of his sampled: he has no relative in F) and one daughter becomes a floater;
#   * the breeding male of B dies in Aug 2022 and is replaced by a male from outside, who sires the 2023 litter;
#   * the breeding female of G dies in Nov 2022 and is replaced in March 2023 by a female from outside, too late to breed;
#   * a male from outside crosses the area between the territories in winter 2021-22, a few others leave a single sample;
#   * other young disperse out of the study area at 1 to 3 years old, a few stay as helpers;
#   * adults leave 3% of their samples in a neighbouring territory (forays);
#   * sampling is lighter in 2022-23 (half the effort) and pack E, on the edge, is poorly sampled.
# Coordinates are placed in the French Alps (Maurienne / Belledonne) only to get a background map: they mean nothing.
simulate_population <- function(seed = 22, foray = 0.03) {
  set.seed(seed)
  campaigns <- c("2020-21", "2021-22", "2022-23", "2023-24")
  start <- as.Date("2020-07-01"); end <- as.Date("2024-06-30"); END <- as.Date("2024-12-31")
  effort <- c("2020-21" = 1, "2021-22" = 1, "2022-23" = 0.5, "2023-24" = 1)     # sampling effort of each year
  campaign_of <- function(d) { y <- as.integer(format(d, "%Y")) - (as.integer(format(d, "%m")) < 7); paste0(y, "-", substr(y + 1, 3, 4)) }
  D <- function(x) as.Date(x)

  # ---- 1. Territories: centre, elongated shape and two core sites (den, rendezvous site), in km in a local frame
  terr <- data.frame(pack = c("A", "B", "C", "F", "D", "E", "G", "H", "I"),
                     x = c(0, 15, 30, 7.5, 22.5, 37.5, 0, 15, 7.5), y = c(26, 26, 26, 13, 13, 13, 0, 0, 39))
  n_t <- nrow(terr)
  terr <- within(terr, {
    x <- x + rnorm(n_t, 0, 1); y <- y + rnorm(n_t, 0, 1)                         # not a perfect lattice
    theta <- runif(n_t, 0, pi); sd1 <- runif(n_t, 3, 3.6); sd2 <- runif(n_t, 2.2, 2.7)
    cx1 <- x + rnorm(n_t, 0, 1.5); cy1 <- y + rnorm(n_t, 0, 1.5)
    cx2 <- x + rnorm(n_t, 0, 1.5); cy2 <- y + rnorm(n_t, 0, 1.5)
  })
  rownames(terr) <- terr$pack
  sample_factor <- c(A = 1, B = 1, C = 1, D = 1, E = 0.35, F = 1, G = 1, H = 1, I = 1)  # E is on the edge, poorly sampled
  neighbours <- lapply(stats::setNames(terr$pack, terr$pack), function(p)        # packs next door (centres < 20 km apart)
    setdiff(terr$pack[sqrt((terr$x - terr[p, "x"])^2 + (terr$y - terr[p, "y"])^2) < 20], p))

  terr_xy <- function(p, n) {                            # n locations in territory p: 35% around the core sites, the rest spread out
    t <- terr[p, ]
    u1 <- rnorm(n, 0, t$sd1); u2 <- rnorm(n, 0, t$sd2)
    xy <- cbind(t$x + u1 * cos(t$theta) - u2 * sin(t$theta), t$y + u1 * sin(t$theta) + u2 * cos(t$theta))
    core <- runif(n) < 0.35; k <- sample(1:2, n, TRUE)
    xy[core, ] <- cbind(c(t$cx1, t$cx2)[k], c(t$cy1, t$cy2)[k])[core, , drop = FALSE] + rnorm(2 * sum(core), 0, 0.9)
    xy
  }
  float_xy <- function() {                               # floaters: mostly between and around the territories
    repeat {
      p <- c(runif(1, -12, 47), runif(1, -10, 49))
      if (min(sqrt((terr$x - p[1])^2 + (terr$y - p[2])^2)) > 8 || runif(1) < 0.1) return(p)
    }
  }
  path <- rbind(c(-11, 19), c(7.5, 21.7), c(15, 17.3), c(22.5, 21.7), c(30, 17.3), c(47, 21))  # through the gaps of the lattice
  path_xy <- function(f) {                               # point at fraction f (0-1) of the transient's journey
    L <- c(0, cumsum(sqrt(rowSums(diff(path)^2)))); s <- f * max(L)
    i <- min(findInterval(s, L), nrow(path) - 1); w <- (s - L[i]) / (L[i + 1] - L[i])
    path[i, ] * (1 - w) + path[i + 1, ] * w + rnorm(2, 0, 1)
  }

  # ---- 2. Genetics: 22 microsatellite loci, 5 to 10 alleles each, made-up allele sizes and frequencies
  n_loci <- 22
  loci <- sprintf("MS%02d", seq_len(n_loci))
  alleles <- lapply(seq_len(n_loci), function(l) { k <- sample(5:10, 1); sample(seq(90, 260, 2), 1) + sample(c(2, 4), 1) * (0:(k - 1)) })
  freqs <- lapply(alleles, function(a) { g <- rgamma(length(a), 0.8); g / sum(g) })
  founder_geno <- function() sapply(seq_len(n_loci), function(l) sample(alleles[[l]], 2, TRUE, freqs[[l]]))   # 2 x loci
  child_geno <- function(m, f) sapply(seq_len(n_loci), function(l) c(m[sample(2, 1), l], f[sample(2, 1), l]))

  # ---- 3. Individuals, their stays (where they live, from when to when) and the breeders
  ind <- data.frame(id = character(), sex = character(), mother = character(), father = character(), born = as.Date(character()),
                    origin = character(), natal = character(), dispersal = as.Date(character()), stringsAsFactors = FALSE)
  geno <- list()
  stays <- data.frame(id = character(), place = character(), from = as.Date(character()), to = as.Date(character()))
  breeders <- data.frame(pack = character(), id = character(), from = as.Date(character()), to = as.Date(character()))
  new_ind <- function(id, sex, origin, mother = NA, father = NA, born = NA, natal = NA) {
    ind[nrow(ind) + 1, ] <<- list(id, sex, mother, father, as.Date(born), origin, natal, as.Date(NA))
    geno[[id]] <<- if (is.na(mother)) founder_geno() else child_geno(geno[[mother]], geno[[father]])
  }
  stay <- function(id, place, from, to) stays[nrow(stays) + 1, ] <<- list(id, place, as.Date(from), as.Date(to))
  breeder <- function(pack, id, from, to) { breeders[nrow(breeders) + 1, ] <<- list(pack, id, as.Date(from), as.Date(to)); stay(id, pack, from, to) }

  # Breeding pairs already there at the start (since 2017); two pairs of full sibs come from never-sampled grandparents
  for (g in c("gBD_f", "gBD_m", "gAG_f", "gAG_m")) new_ind(g, toupper(substr(g, 5, 5)), "Unsampled")
  sibs <- c(A_f = "AG", G_m = "AG", B_m = "BD", D_f = "BD")
  for (p in c("A", "B", "C", "D", "E", "F", "G", "H")) for (sx in c("f", "m")) {
    id <- paste0(p, "_", sx); s <- sibs[id]
    if (is.na(s)) new_ind(id, toupper(sx), "Founder") else
      new_ind(id, toupper(sx), "Founder", paste0("g", s, "_f"), paste0("g", s, "_m"), "2015-05-01")
  }
  deaths <- c(B_m = "2022-08-20", C_f = "2021-10-05", C_m = "2022-02-12", F_m = "2023-06-15", G_f = "2022-11-20")
  for (id in ind$id[ind$origin == "Founder"]) breeder(substr(id, 1, 1), id, "2017-01-01", if (id %in% names(deaths)) deaths[[id]] else END)

  # Newcomers from outside the study area
  new_ind("I_m", "M", "Immigrant"); breeder("I", "I_m", "2021-11-20", END)          # founds I with the daughter of A
  new_ind("B_m2", "M", "Immigrant"); breeder("B", "B_m2", "2022-11-01", END)        # replaces the dead breeder of B
  new_ind("G_f2", "F", "Immigrant"); breeder("G", "G_f2", "2023-03-01", END)        # replaces the dead breeder of G, too late to breed
  new_ind("T1", "M", "Immigrant"); stay("T1", "pass", "2021-11-05", "2022-02-25")   # crosses the area
  for (v in 1:3) {                                                                   # leave a single sample on their way
    new_ind(paste0("V", v), sample(c("M", "F"), 1), "Immigrant")
    d <- start + sample(0:(end - start), 1); stay(paste0("V", v), "float", d, d)
  }

  # The scripted young: the daughter of A who founds I, and two young of C (stays written now, so that I breeds in 2022)
  stay("dA", "A", "2019-05-01", "2021-12-14"); breeder("I", "dA", "2021-12-15", END)
  stay("dC1", "C", "2020-05-01", "2023-07-15"); breeder("F", "dC1", "2023-08-01", END)
  stay("dC2", "C", "2020-05-01", "2022-09-30"); stay("dC2", "float", "2022-10-01", END)

  # Litters: born in May when a breeding pair is together at mating time (mid-February)
  pair_on <- function(p, d) {
    b <- breeders[breeders$pack == p & breeders$from <= d & breeders$to >= d & breeders$id %in% ind$id, ]
    s <- ind$sex[match(b$id, ind$id)]
    if (any(s == "F") && any(s == "M")) c(f = b$id[s == "F"][1], m = b$id[s == "M"][1]) else NULL
  }
  scripted <- list("A.2019" = c(dA = "F"), "C.2020" = c(dC1 = "M", dC2 = "F"))    # pups whose fate is part of the story
  for (y in 2018:2023) for (p in terr$pack) {
    pr <- pair_on(p, D(paste0(y, "-02-15"))); if (is.null(pr)) next
    if (breeders$to[breeders$id == pr["f"]] < D(paste0(y, "-05-01"))) next        # mother dead before birth
    sc <- scripted[[paste(p, y, sep = ".")]]
    n <- max(1, rpois(1, 4.5)); born <- D(paste0(y, "-05-01")) + sample(-10:10, 1)
    ids <- c(names(sc), sprintf("%s%d_%d", p, y, seq_len(n)))
    for (k in seq_along(ids)) {
      id <- ids[k]
      new_ind(id, if (k <= length(sc)) sc[[k]] else sample(c("M", "F"), 1), "Born", pr["f"], pr["m"], born, p)
      if (k <= length(sc)) next                                                      # their stays are written above
      if (runif(1) > 0.6) next                                                       # dead before autumn: never sampled
      death <- D(paste0(y, "-10-01")) + round(rexp(1, -log(0.72)) * 365)             # 72% yearly survival
      age_disp <- sample(c(1, 2, 3, NA), 1, prob = c(0.2, 0.45, 0.25, 0.1))         # dispersal at 1-3 years, 10% stay as helpers
      disp <- if (is.na(age_disp)) END else D(paste0(y + age_disp - 1, "-11-01")) + sample(0:150, 1)
      stay(id, p, born, min(death, disp))
      if (disp < death && disp <= end) {
        ind$dispersal[ind$id == id] <- disp
        if (runif(1) < 0.15) stay(id, "float", disp, disp + 20)                        # sometimes one sample on its way out
      }
    }
  }
  ind$dispersal[match(c("dA", "dC1", "dC2"), ind$id)] <- D(c("2021-12-15", "2023-07-15", "2022-10-01"))

  # ---- 4. Sampling: detection and number of samples per individual and year, dates mostly in winter, locations from the stays
  month_w <- c(.14, .14, .11, .07, .04, .04, .04, .04, .06, .08, .11, .13)        # Jan to Dec: snow tracking in winter
  is_breeder <- function(id, d) any(breeders$id == id & breeders$from <= d & breeders$to >= d)
  forced <- list(dC1 = c("2020-21" = 0, "2021-22" = 2, "2022-23" = 1, "2023-24" = 6),  # few samples at home, many as a breeder in F
                 dC2 = c("2020-21" = 1, "2021-22" = 1, "2022-23" = 2, "2023-24" = 4),
                 G_f2 = c("2022-23" = 2, "2023-24" = 4), T1 = c("2021-22" = 6))
  dead_found <- c(C_f = "2021-10-05", B_m = "2022-08-20", G_f = "2022-11-20")       # found dead: a tissue sample
  smp <- data.frame(id = character(), date = as.Date(character()), place = character(), x = numeric(), y = numeric(), type = character())

  for (id in ind$id[ind$origin != "Unsampled"]) {
    st <- stays[stays$id == id, ]; if (!nrow(st)) next
    born <- ind$born[ind$id == id]
    place_of <- function(d) { i <- which(st$from <= d & st$to >= d); st$place[i[length(i)]] }
    for (cp in campaigns) {
      c0 <- D(paste0(substr(cp, 1, 4), "-07-01")); c1 <- D(paste0(as.integer(substr(cp, 1, 4)) + 1, "-06-30"))
      days <- do.call(c, lapply(seq_len(nrow(st)), function(i) { a <- max(st$from[i], c0); b <- min(st$to[i], c1)
        if (a <= b) seq(a, b, by = "day") else as.Date(character()) }))
      if (!is.na(born)) days <- days[days >= born + 150]                             # pups are sampled from October
      if (!length(days)) next
      mid <- days[ceiling(length(days) / 2)]; pl <- place_of(mid)
      role <- if (is_breeder(id, mid)) "breeder" else if (pl %in% c("float", "pass")) "floater" else
        if (!is.na(born) && mid - born < 365) "pup" else "young"
      fac <- effort[[cp]] * (if (pl %in% names(sample_factor)) sample_factor[[pl]] else 1) * length(days) / 365
      if (!is.null(forced[[id]])) n <- if (cp %in% names(forced[[id]])) forced[[id]][[cp]] else 0 else {
        p0 <- c(breeder = 0.95, young = 0.75, pup = 0.55, floater = 0.8)[[role]]
        n <- if (runif(1) < 1 - (1 - p0)^fac) 1 + rnbinom(1, size = 1.5, mu = c(breeder = 4, young = 1.8, pup = 1, floater = 2)[[role]] * fac) else 0
        if (st$place[1] == "float" && st$from[1] == st$to[1]) n <- 1                 # visitors: one sample
      }
      if (n == 0) next
      dd <- sort(days[sample(length(days), n, replace = TRUE, prob = month_w[as.integer(format(days, "%m"))])])
      for (d in as.list(dd)) {
        pl <- place_of(d)
        out <- pl %in% terr$pack && (is.na(born) || d - born > 365) && runif(1) < foray   # adults sometimes leave a sample next door
        xy <- if (pl == "pass") path_xy(as.numeric(d - D("2021-11-05")) / as.numeric(D("2022-02-25") - D("2021-11-05"))) else
          if (pl == "float") float_xy() else if (out) terr_xy(sample(neighbours[[pl]], 1), 1) else terr_xy(pl, 1)
        m <- as.integer(format(d, "%m"))
        ty <- sample(c("Scat", "Urine", "Hair", "Saliva"), 1, prob = if (m %in% c(12, 1, 2, 3)) c(.66, .2, .07, .07) else c(.84, 0, .08, .08))
        smp[nrow(smp) + 1, ] <- list(id, d, pl, xy[1], xy[2], ty)
      }
    }
  }
  for (id in names(dead_found)) {
    pl <- stays$place[stays$id == id][1]; xy <- terr_xy(pl, 1)
    smp[nrow(smp) + 1, ] <- list(id, D(dead_found[[id]]), pl, xy[1], xy[2], "Tissue")
  }

  # ---- 5. Genotypes of the sampled individuals: a few missing loci, rare errors
  first <- tapply(smp$date, smp$id, min)
  sampled <- names(first)[order(first, names(first))]
  new_id <- stats::setNames(sprintf("W%03d", seq_along(sampled)), sampled)         # IDs in order of first sampling
  G <- t(sapply(sampled, function(id) {
    g <- geno[[id]]
    for (k in which(runif(2 * n_loci) < 0.003)) g[k] <- sample(alleles[[ceiling(k / 2)]], 1)   # rare genotyping errors
    g[, runif(n_loci) < 0.03] <- NA                                                  # 3% of the loci not typed
    as.vector(g)                                                                     # MS01_1, MS01_2, MS02_1...
  }))
  genotypes <- data.frame(Individual = unname(new_id[sampled]), G, row.names = NULL)
  names(genotypes)[-1] <- paste0(rep(loci, each = 2), c("_1", "_2"))
  genotypes[-1] <- lapply(genotypes[-1], as.integer)

  # ---- 6. Samples as sf points (WGS84), relatedness, and the truth
  origin <- st_coordinates(st_transform(st_sfc(st_point(c(6.3, 45.15)), crs = 4326), 32632))
  smp <- smp[order(smp$date, smp$id), ]
  pts <- st_transform(st_as_sf(data.frame(X = origin[1] + smp$x * 1000, Y = origin[2] + smp$y * 1000), coords = c("X", "Y"), crs = 32632), 4326)
  geom <- st_sfc(lapply(st_geometry(pts), function(p) st_point(round(unclass(p), 5))), crs = 4326)
  samples <- st_sf(Individual = unname(new_id[smp$id]), Sex = ind$sex[match(smp$id, ind$id)], Date = smp$date,
                   Campaign = campaign_of(smp$date), Sample = sprintf("S%04d", seq_len(nrow(smp))), Type = smp$type,
                   geometry = geom)

  # Truth, one row per sampled individual: its pack (where most of its samples come from; NA for floaters), role,
  # origin, natal pack, parents (when sampled), dispersal date, and its pack in each year (where it lived longest)
  tr <- ind[match(sampled, ind$id), ]
  main <- tapply(smp$place, smp$id, function(p) { t <- table(p); t <- t[t == max(t)]; names(t)[length(t)] })[sampled]
  yearly <- sapply(campaigns, function(cp) {
    c0 <- D(paste0(substr(cp, 1, 4), "-07-01")); c1 <- D(paste0(as.integer(substr(cp, 1, 4)) + 1, "-06-30"))
    sapply(sampled, function(id) {
      st <- stays[stays$id == id, ]
      len <- pmax(0, as.numeric(pmin(st$to, c1) - pmax(st$from, c0)) + 1)
      if (!any(len > 0)) NA_character_ else { pl <- st$place[which.max(len)]; if (pl %in% c("float", "pass")) "floater" else pl }
    })
  })
  truth <- data.frame(
    Individual = unname(new_id[sampled]), Sex = tr$sex,
    Pack = unname(ifelse(main %in% c("float", "pass"), NA, main)),
    Role = ifelse(sampled %in% breeders$id, "Breeder", ifelse(main %in% c("float", "pass"), "Floater", "Offspring")),
    Origin = tr$origin, Natal_Pack = tr$natal,
    Mother = unname(new_id[tr$mother]), Father = unname(new_id[tr$father]),
    Birth_Year = ifelse(tr$origin == "Born", as.integer(format(tr$born, "%Y")), NA_integer_),
    Dispersal = tr$dispersal,
    stats::setNames(as.data.frame(yearly, row.names = NULL), paste0("Pack_", sub("-", "_", campaigns))),
    stringsAsFactors = FALSE, row.names = NULL)
  truth <- truth[order(truth$Individual), ]; rownames(truth) <- NULL

  list(samples_example = samples, genotypes_example = genotypes, relate_example = relatedness(genotypes), truth_example = truth)
}

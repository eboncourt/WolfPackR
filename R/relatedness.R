#' @title Pairwise relatedness from microsatellite genotypes
#' @description Computes pairwise relatedness between all individuals of a genotype table, in base R
#'   (no call to COLONY, COANCESTRY or the `related` package). Three estimators are available:
#'   * `QG`: Queller & Goodnight (1989) moment estimator;
#'   * `LR`: Lynch & Ritland (1999) moment estimator;
#'   * `LOD`: log10 likelihood ratio of a first-degree relationship (parent-offspring or full sibs,
#'     whichever fits best) against unrelated, with a small genotyping error rate.
#'
#'   QG and LR are averaged over both directions of the pair (x as reference, then y), so they are symmetric.
#'   They are expected to be about 0.5 for first-degree relatives, 0.25 for second-degree relatives and 0
#'   for unrelated individuals. LOD > 3 means the first-degree hypothesis is 1000 times more likely than "unrelated".
#' @param genotypes A data.frame with one row per individual: an ID column, then the alleles.
#'   Either two consecutive columns per locus (e.g. `CPH9_1`, `CPH9_2`) when `sep = NULL`,
#'   or one column per locus holding both alleles (e.g. `"144 150"`) when `sep` is given.
#' @param id Name or position of the ID column (default: the first column).
#' @param sep Separator between the two alleles when each locus is stored in a single column
#'   (e.g. `" "` or `"/"`). `NULL` (default) means two columns per locus.
#' @param estimators Which estimators to compute, any of `"QG"`, `"LR"`, `"LOD"` (default: all three).
#' @param freq_ids IDs of the individuals used to estimate allele frequencies (default: all).
#'   Handy to leave out hybrids or individuals from another population.
#' @param missing Values treated as missing alleles (default: `NA`, `0`, `"0"`, `""` and `"NA"`).
#' @param error_rate Genotyping error rate used by the LOD estimator (default: 0.01).
#' @return A data.frame with one row per pair: `ind1`, `ind2`, one column per estimator and `n_loci`,
#'   the number of loci typed in both individuals. It can be passed straight to [genetic_groups()].
#' @references
#' Lynch, M. & Ritland, K. (1999). Estimation of pairwise relatedness with molecular markers. Genetics, 152, 1753-1766.
#'
#' Queller, D. C. & Goodnight, K. F. (1989). Estimating relatedness using genetic markers. Evolution, 43, 258-275.
#' @export
#' @examples
#' geno <- data.frame(Individual = c("A", "B", "C", "D"),
#'                    L1_1 = c(100, 100, 102, 104), L1_2 = c(102, 104, 104, 106),
#'                    L2_1 = c(150, 150, 152, 150), L2_2 = c(152, 154, 154, 156))
#' relatedness(geno)
#' # same thing with one column per locus
#' geno2 <- data.frame(Individual = geno$Individual,
#'                     L1 = paste(geno$L1_1, geno$L1_2), L2 = paste(geno$L2_1, geno$L2_2))
#' relatedness(geno2, sep = " ")
relatedness <- function(genotypes, id = 1, sep = NULL, estimators = c("QG", "LR", "LOD"),
                        freq_ids = NULL, missing = c(NA, 0, "0", "", "NA"), error_rate = 0.01) {

  # Read the table: IDs + two allele matrices (individuals x loci), alleles kept as character
  genotypes <- as.data.frame(genotypes, stringsAsFactors = FALSE)
  estimators <- match.arg(estimators, c("QG", "LR", "LOD"), several.ok = TRUE)
  ids <- as.character(genotypes[[id]])
  if (anyDuplicated(ids)) stop("Duplicated IDs in 'genotypes'.")
  geno <- genotypes[, -match(if (is.numeric(id)) names(genotypes)[id] else id, names(genotypes)), drop = FALSE]

  if (is.null(sep)) {                                      # two columns per locus
    if (ncol(geno) %% 2 != 0) stop("With sep = NULL, 'genotypes' needs two columns per locus (odd number of allele columns).")
    first <- seq(1, ncol(geno), by = 2)
    A1 <- as.matrix(geno[, first, drop = FALSE]); A2 <- as.matrix(geno[, first + 1, drop = FALSE])
    loci <- sub("[._-]?[12abAB]$", "", names(geno)[first])  # CPH9_1 -> CPH9
  } else {                                                 # one column per locus: "144 150"
    parts <- lapply(geno, function(x) do.call(rbind, lapply(strsplit(trimws(as.character(x)), sep, fixed = TRUE),
                                                            function(a) { a <- a[a != ""]; if (length(a) == 2) a else c(NA, NA) })))
    A1 <- sapply(parts, function(p) p[, 1]); A2 <- sapply(parts, function(p) p[, 2])
    loci <- names(geno)
  }
  A1 <- matrix(trimws(as.character(A1)), nrow = length(ids)); A2 <- matrix(trimws(as.character(A2)), nrow = length(ids))
  A1[A1 %in% as.character(missing) | is.na(A1)] <- NA; A2[A2 %in% as.character(missing) | is.na(A2)] <- NA
  A1[is.na(A2)] <- NA; A2[is.na(A1)] <- NA                 # a locus with one allele missing is missing
  dimnames(A1) <- dimnames(A2) <- list(ids, loci)

  # Allele frequencies, from freq_ids only; an allele never seen there gets a tiny frequency instead of 0
  if (is.null(freq_ids)) freq_ids <- ids
  if (!all(freq_ids %in% ids)) stop("Some 'freq_ids' are not in 'genotypes'.")
  freq <- function(l) { x <- c(A1[freq_ids, l], A2[freq_ids, l]); x <- x[!is.na(x)]; table(x) / length(x) }
  fq <- lapply(loci, freq)
  pfun <- function(l, a) { f <- as.numeric(fq[[l]][a]); f[is.na(f) & !is.na(a)] <- 1 / (2 * length(freq_ids) + 1); f }

  # All pairs, then we add up each locus' contribution for every pair at once
  pr <- utils::combn(length(ids), 2)
  i <- pr[1, ]; j <- pr[2, ]; np <- length(i)
  qg_num <- qg_den <- lr_num <- lr_den <- matrix(0, np, 2)   # column 1: x = ind1, column 2: x = ind2
  ll <- matrix(0, np, 3, dimnames = list(NULL, c("U", "PO", "FS")))
  n_loci <- integer(np)

  for (l in seq_along(loci)) {
    a1 <- A1[, l]; a2 <- A2[, l]; p1 <- pfun(l, a1); p2 <- pfun(l, a2)
    ok <- !is.na(a1[i]) & !is.na(a1[j])                    # locus typed in both individuals
    n_loci <- n_loci + ok

    for (s in 1:2) {                                       # x = reference individual, y = the other one
      x <- if (s == 1) i else j; y <- if (s == 1) j else i
      a <- a1[x]; b <- a2[x]; c <- a1[y]; d <- a2[y]; pa <- p1[x]; pb <- p2[x]
      Sab <- a == b
      ac <- a == c; ad <- a == d; bc <- b == c; bd <- b == d

      # Queller & Goodnight: shared alleles minus what you'd expect by chance
      qn <- 0.5 * (ac + ad + bc + bd) - pa - pb
      qd <- 1 + Sab - pa - pb
      qg_num[ok, s] <- qg_num[ok, s] + qn[ok]
      qg_den[ok, s] <- qg_den[ok, s] + qd[ok]

      # Lynch & Ritland: per-locus estimate weighted by its inverse variance
      den <- (1 + Sab) * (pa + pb) - 4 * pa * pb
      use <- ok & den > 1e-12                              # degenerate loci (e.g. two alleles at 0.5) carry no info
      w <- den / (2 * pa * pb)
      r_l <- (pa * (bc + bd) + pb * (ac + ad) - 4 * pa * pb) / den
      lr_num[use, s] <- lr_num[use, s] + (w * r_l)[use]
      lr_den[use, s] <- lr_den[use, s] + w[use]
    }

    # Likelihoods (up to P(Gx), which cancels in the ratio): y given x, under 0, 1 or 2 alleles identical by descent
    if ("LOD" %in% estimators) {
      a <- a1[i]; b <- a2[i]; c <- a1[j]; d <- a2[j]; pc <- p1[j]; pd <- p2[j]
      P0 <- ifelse(c == d, pc^2, 2 * pc * pd)                                    # y drawn from the population
      other <- function(s) ifelse(c == s, pd, ifelse(d == s, pc, 0))            # freq of y's non-shared allele
      P1 <- 0.5 * (other(a) + other(b))                                         # one of x's alleles passed to y
      P2 <- as.numeric((a == c & b == d) | (a == d & b == c))                   # same genotype
      P1 <- (1 - error_rate) * P1 + error_rate * P0                             # genotyping error: some chance
      P2 <- (1 - error_rate) * P2 + error_rate * P0                             # of looking like a random genotype
      ll[ok, "U"] <- ll[ok, "U"] + log(P0[ok])
      ll[ok, "PO"] <- ll[ok, "PO"] + log(P1[ok])
      ll[ok, "FS"] <- ll[ok, "FS"] + log(0.25 * P0[ok] + 0.5 * P1[ok] + 0.25 * P2[ok])
    }
  }

  # Put it all together, one row per pair
  out <- data.frame(ind1 = ids[i], ind2 = ids[j], stringsAsFactors = FALSE)
  if ("QG" %in% estimators) out$QG <- rowMeans(qg_num / qg_den)
  if ("LR" %in% estimators) out$LR <- rowMeans(lr_num / lr_den)
  if ("LOD" %in% estimators) out$LOD <- (pmax(ll[, "PO"], ll[, "FS"]) - ll[, "U"]) / log(10)
  out$n_loci <- n_loci
  out[out$n_loci == 0, intersect(c("QG", "LR", "LOD"), names(out))] <- NA   # nothing in common: no estimate
  out
}

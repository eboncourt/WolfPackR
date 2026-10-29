# Simulated genotypes: random allele frequencies, unrelated founders, and kids drawn from two parents
sim_freqs <- function(n_loci, n_alleles = 8) {
  set.seed(7)
  lapply(seq_len(n_loci), function(l) { p <- stats::rgamma(n_alleles, 2); stats::setNames(p / sum(p), 100 + 2 * seq_len(n_alleles)) })
}
draw_ind <- function(freqs) unlist(lapply(freqs, function(p) as.integer(sample(names(p), 2, replace = TRUE, prob = p))))
draw_kid <- function(mum, dad) {                          # one allele from each parent, locus by locus
  n <- length(mum) / 2
  unlist(lapply(seq_len(n), function(l) c(mum[2 * l - sample(0:1, 1)], dad[2 * l - sample(0:1, 1)])))
}
as_geno <- function(rows, ids) {
  g <- as.data.frame(do.call(rbind, rows))
  names(g) <- paste0("L", rep(seq_len(ncol(g) / 2), each = 2), "_", 1:2)
  cbind(Individual = ids, g)
}

test_that("relatedness returns one row per pair with the requested estimators", {
  geno <- pkg_data("genotypes_Carpathians")[1:5, ]
  res <- relatedness(geno)
  expect_equal(nrow(res), choose(5, 2))
  expect_named(res, c("ind1", "ind2", "QG", "LR", "LOD", "n_loci"))
  expect_named(relatedness(geno, estimators = "QG"), c("ind1", "ind2", "QG", "n_loci"))
})

test_that("an individual and its copy have QG = LR = 1 and a high LOD", {
  geno <- pkg_data("genotypes_Carpathians")[1:4, ]
  geno <- rbind(geno, transform(geno[1, ], Individual = "copy"))
  res <- relatedness(geno)
  self <- res[res$ind1 == geno$Individual[1] & res$ind2 == "copy", ]
  expect_equal(self$QG, 1)
  expect_equal(self$LR, 1)
  expect_gt(self$LOD, 3)
})

test_that("relatedness is about 0.5 for parent-offspring and 0 for unrelated pairs", {
  freqs <- sim_freqs(200)                                  # lots of loci so the estimates are tight
  set.seed(1)
  founders <- replicate(12, draw_ind(freqs), simplify = FALSE)
  kids <- lapply(1:6, function(k) draw_kid(founders[[2 * k - 1]], founders[[2 * k]]))
  geno <- as_geno(c(founders, kids), c(paste0("F", 1:12), paste0("K", 1:6)))
  res <- relatedness(geno)
  fk <- res[res$ind1 %in% paste0("F", 1:12) & res$ind2 %in% paste0("K", 1:6), ]   # founder x kid pairs
  po <- fk[ceiling(as.integer(sub("F", "", fk$ind1)) / 2) == as.integer(sub("K", "", fk$ind2)), ]  # K1 = F1 x F2, etc.
  un <- res[res$ind1 %in% paste0("F", 1:12) & res$ind2 %in% paste0("F", 1:12), ]
  expect_equal(nrow(po), 12)
  for (e in c("QG", "LR")) {
    expect_lt(abs(mean(po[[e]]) - 0.5), 0.1)
    expect_lt(abs(mean(un[[e]])), 0.1)
  }
  expect_true(all(po$LOD > 3))
  expect_true(all(un$LOD < 0))
})

test_that("one column per locus with a separator gives the same result as two columns", {
  geno <- pkg_data("genotypes_Carpathians")[1:6, 1:9]     # 4 loci are enough
  one_col <- data.frame(Individual = geno$Individual)
  for (k in seq(2, ncol(geno), by = 2)) one_col[[sub("_1$", "", names(geno)[k])]] <- paste(geno[[k]], geno[[k + 1]], sep = "/")
  one_col[2, 2] <- "NA/NA"                                 # missing values written as text
  geno[2, 2:3] <- NA
  expect_equal(relatedness(one_col, sep = "/"), relatedness(geno))
})

test_that("missing loci are skipped, whatever the code used for missing", {
  geno <- pkg_data("genotypes_Carpathians")[1:3, ]
  n_full <- relatedness(geno)$n_loci[1]
  geno[1, 2:3] <- 0                                        # first locus missing for the first individual
  res <- relatedness(geno)
  expect_equal(res$n_loci[1], n_full - 1)
  expect_false(anyNA(res$QG))
})

test_that("relatedness checks its inputs", {
  geno <- pkg_data("genotypes_Carpathians")[1:3, ]
  expect_error(relatedness(geno[, -2]), "two columns per locus")
  expect_error(relatedness(geno, freq_ids = "nobody"), "freq_ids")
  expect_error(relatedness(rbind(geno, geno[1, ])), "Duplicated")
})

test_that("relate_Carpathians is what relatedness() gives on genotypes_Carpathians", {
  geno <- pkg_data("genotypes_Carpathians")
  rel <- relatedness(geno, freq_ids = setdiff(geno$Individual, "EF.15J8"))
  expect_equal(rel, pkg_data("relate_Carpathians"))
})

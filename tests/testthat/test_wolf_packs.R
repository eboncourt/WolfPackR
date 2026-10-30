test_that("wolf_packs runs the whole workflow on the Carpathians and matches Iosif et al. reasonably", {
  res <- suppressMessages(wolf_packs(pkg_data("samples_Carpathians"), pkg_data("relate_Carpathians"), "QG", 0.25,
                                     period = "Campaign"))
  expect_s3_class(res, "wolfpack")
  expect_equal(nrow(res$individuals), 47)
  expect_true(all(c("Pack", "Status", "Genetic_Group", "Group") %in% names(res$samples)))
  cmp <- compare_partitions(res$individuals, pkg_data("reference_Carpathians"))$summary
  expect_gte(cmp$ARI, 0.7)                                             # 0.80 when written
  expect_gte(cmp$reference_groups_found, 5)
  expect_output(print(res), "packs")
})

test_that("wolf_packs finds the true packs of the simulated population", {
  res <- suppressMessages(wolf_packs(pkg_data("samples_example"), pkg_data("relate_example"), "QG", 0.25,
                                     period = "Campaign"))
  expect_equal(nrow(res$individuals), 56)
  cmp <- compare_partitions(res$individuals, pkg_data("truth_example"))$summary
  expect_gte(cmp$ARI, 0.85)                                            # 0.94 when written
  expect_equal(cmp$reference_groups_found, 9)

  # the son of C who became the breeding male of F in 2023-24: a dated disperser
  tr <- pkg_data("truth_example")
  son <- tr$Individual[tr$Natal_Pack %in% "C" & tr$Pack %in% "F"]
  expect_equal(res$individuals$Status[res$individuals$Individual == son], "Disperser")
  expect_equal(res$individuals$change_period[res$individuals$Individual == son], "2023-24")
})

test_that("every threshold can be set from wolf_packs(), and the settings are kept", {
  obs <- sim_samples(); rel <- sim_relate()
  res <- suppressMessages(wolf_packs(obs, rel, "r", 0.25, h = 2000, period = "Campaign", method = "components"))
  expect_equal(res$settings$outside, 1)
  expect_equal(res$settings$wide_trim, 0.1)
  expect_equal(res$settings$h, 2000)
  # a much stricter "outside" rule sends more samples outside every territory
  strict <- suppressMessages(wolf_packs(obs, rel, "r", 0.25, h = 2000, period = "Campaign", method = "components",
                                        outside = 50))
  n_out <- function(r) sum(attr(r$assignment, "samples")$Territory == "outside")
  expect_gt(n_out(strict), n_out(res))
  expect_equal(strict$settings$outside, 50)
})

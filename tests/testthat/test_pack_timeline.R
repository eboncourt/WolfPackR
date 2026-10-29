test_that("pack_timeline follows individuals and territories from one period to the next", {
  obs <- sim_samples()
  obs$Pack <- ifelse(obs$Individual %in% c("A1", "A2", "A3", "U"), "A",
                     ifelse(obs$Individual %in% c("B1", "B2", "B3", "D", "M"), "B", NA))
  tl <- pack_timeline(obs, "Pack", "Campaign", h = 2000)
  ind <- tl$individuals
  expect_equal(ind$pattern[ind$Individual == "M"], "changed pack")
  expect_equal(ind$change_period[ind$Individual == "M"], "P2")
  expect_equal(ind$pattern[ind$Individual == "A1"], "stable")
  expect_equal(ind$pattern[ind$Individual == "A2"], "single period")
  expect_true(all(c("P1", "P2", "rg_km_P1", "rg_km_P2") %in% names(ind)))
  expect_equal(nrow(tl$samples), nrow(obs))
  expect_named(tl$territories, c("Pack", "period", "n_samples", "n_detected", "n_new", "n_missed",
                                 "affinity_with_all_periods", "shift_km", "shared_pct"))
  expect_true(all(tl$territories$affinity_with_all_periods > 0.5, na.rm = TRUE))
})

test_that("pack_timeline counts members per period and leaves out those living elsewhere at the time", {
  obs <- sim_samples()
  obs$Pack <- ifelse(obs$Individual %in% c("A1", "A2", "A3", "U"), "A",
                     ifelse(obs$Individual %in% c("B1", "B2", "B3", "D", "M"), "B", NA))
  tt <- pack_timeline(obs, "Pack", "Campaign", h = 2000)$territories
  B1 <- tt[tt$Pack == "B" & tt$period == "P1", ]; B2 <- tt[tt$Pack == "B" & tt$period == "P2", ]
  expect_equal(B1$n_detected, 3)                                     # B1, B2, D (M was still in A's area)
  expect_equal(B2$n_detected, 3)                                     # B1, B3 and M, just arrived
  expect_equal(B2$n_new, 2)                                          # B3 and M
  expect_true(is.na(B1$n_new))                                       # nothing to compare with in the first period
  expect_equal(B1$n_samples, 18)                                     # M's samples in A's area don't draw B's territory
  expect_equal(tt$shared_pct, rep(0, 4))                             # 20 km apart: no shared space
  expect_true(all(tt$shift_km < 2))
})

sim_with_packs <- function() {
  obs <- sim_samples()
  obs$Pack <- ifelse(obs$Individual %in% c("A1", "A2", "A3", "U"), "A", ifelse(obs$Individual %in% c("B1", "B2", "B3", "D"), "B", NA))
  obs
}

test_that("summarise_packs gives one row per pack, with areas and composition", {
  sm <- summarise_packs(sim_with_packs(), "Pack", h = 2000)
  expect_equal(sm$Pack, c("A", "B"))
  expect_equal(sm$n_individuals, c(4, 4))
  expect_true(all(sm$area95_km2 > sm$area50_km2))
  expect_equal(sm$n_males + sm$n_females, sm$n_individuals)
  expect_false(any(grepl("Dominant|dominant", names(sm))))           # no more putative dominants
  ind <- attr(sm, "individuals")
  expect_equal(nrow(ind), 11)
  expect_true(all(c("Individual", "Pack", "n_samples", "Sex") %in% names(ind)))
})

test_that("within-pack relatedness is tested against random packs", {
  sm <- summarise_packs(sim_with_packs(), "Pack", h = 2000, relate = sim_relate(), estimator = "r", n_perm = 99)
  expect_true(all(c("mean_r", "r_z", "r_p") %in% names(sm)))
  expect_true(all(sm$r_p > 0 & sm$r_p <= 1))
  expect_equal(sm$mean_r, c(0.25, 0.25))                              # 3 related pairs (0.5) out of 6 in each pack
})

test_that("summarise_packs accepts a bare column name and checks it", {
  obs <- sim_with_packs()
  expect_equal(summarise_packs(obs, Pack, h = 2000)$Pack, c("A", "B"))
  expect_error(summarise_packs(obs, "nope", h = 2000), "doesn't exist")
})

test_that("summarise_packs gives the space shared with the other packs and the nearest one", {
  sm <- summarise_packs(sim_with_packs(), "Pack", h = 2000)
  expect_equal(sm$shared_pct, c(0, 0))                                # the two areas are 20 km apart
  expect_equal(sm$nearest_pack, c("B", "A"))
  expect_true(all(abs(sm$nearest_km - 20) < 2))
  expect_false("first_period" %in% names(sm))                         # nothing about time without period or date
  expect_null(attr(sm, "periods"))
})

test_that("summarise_packs follows packs and members through time", {
  obs <- sim_with_packs()
  obs$Date <- as.Date(ifelse(obs$Campaign == "P1", "2020-01-15", "2021-01-15")) + seq_len(nrow(obs)) %% 30
  sm <- summarise_packs(obs, "Pack", h = 2000, period = "Campaign", date = "Date")
  expect_equal(sm$first_period, c("P1", "P1"))
  expect_equal(sm$last_period, c("P2", "P2"))
  expect_equal(sm$n_periods, c(2, 2))
  expect_equal(sm$n_persistent, c(1, 1))                              # only A1 and B1 were sampled in both periods
  expect_true(all(sm$min_affinity > 0.9 & sm$max_shift_km < 2))
  expect_true(all(sm$tenure_years >= 1 & sm$tenure_years < 1.2))
  ind <- attr(sm, "individuals")
  expect_equal(ind$n_periods[ind$Individual %in% c("A1", "A2", "A3")], c(2, 1, 1))
  expect_equal(ind$tenure_years[ind$Individual == "A2"], 0)
  expect_true(all(is.na(ind$tenure_years[is.na(ind$Pack)])))          # no pack, no tenure
  expect_equal(nrow(attr(sm, "periods")), 4)                          # 2 packs x 2 periods
})

test_that("summarise_packs works with dates alone, and checks the time columns", {
  obs <- sim_with_packs()
  obs$Date <- as.Date("2020-01-15") + seq_len(nrow(obs))
  sm <- summarise_packs(obs, "Pack", h = 2000, date = "Date")
  expect_true("tenure_years" %in% names(sm))
  expect_false("first_period" %in% names(sm))
  obs$Date <- as.character(obs$Date)
  expect_error(summarise_packs(obs, "Pack", h = 2000, date = "Date"), "class Date")
  expect_error(summarise_packs(obs, "Pack", h = 2000, period = "nope"), "doesn't exist")
})

test_that("summarise_packs copes with no pack, or a single one", {
  obs <- sim_with_packs()
  obs$Pack <- NA
  expect_equal(nrow(summarise_packs(obs, "Pack", h = 2000, period = "Campaign")), 0)
  obs <- sim_with_packs()
  obs$Pack[obs$Pack %in% "B"] <- NA
  sm <- summarise_packs(obs, "Pack", h = 2000, period = "Campaign")
  expect_equal(sm$shared_pct, 0)
  expect_true(is.na(sm$nearest_pack))
})

sim_assign <- function(...) {
  obs <- sim_samples()
  gg <- genetic_groups(sim_relate(), "r", 0.25, method = "components")
  obs <- merge(obs, gg, by = "Individual")
  assign_individuals(obs, "group", h = 2000, ...)
}

test_that("each special case gets its status", {
  ai <- sim_assign(relate = sim_relate(), estimator = "r", threshold = 0.25)
  st <- stats::setNames(ai$Status, ai$Individual)
  expect_true(all(st[c("A1", "A2", "A3", "B1", "B2", "B3")] == "Resident"))
  expect_equal(st[["D"]], "Putative disperser")
  expect_equal(st[["U"]], "Unrelated resident")
  expect_equal(st[["T"]], "Transient")
  expect_equal(st[["L"]], "Undetermined")                              # one sample: low data
  expect_true(ai$low_data[ai$Individual == "L"])
  g <- stats::setNames(ai$Group, ai$Individual)
  expect_equal(g[["D"]], g[["B1"]])                                    # moved to the group it lives with
  expect_equal(g[["U"]], g[["A1"]])
  expect_true(is.na(g[["T"]]))
})

test_that("with periods, a dated move from the relatives' territory makes a Disperser", {
  ai <- sim_assign(relate = sim_relate(), estimator = "r", threshold = 0.25, period = "Campaign")
  m <- ai[ai$Individual == "M", ]
  expect_equal(m$Status, "Disperser")
  expect_equal(m$change_period, "P2")
  expect_equal(ai$Status[ai$Individual == "D"], "Putative disperser")  # no earlier period in A: stays putative
})

test_that("without the relatedness table, relatives are the genetic group", {
  ai <- sim_assign()
  st <- stats::setNames(ai$Status, ai$Individual)
  expect_equal(st[["U"]], "Unrelated resident")
  expect_equal(st[["T"]], "Transient")
  # coarser rule: M, of the same genetic group as D, also lives in B, so D counts as having a relative there
  expect_equal(st[["D"]], "Resident")
})

test_that("the output has one row per individual and keeps the sample attribution", {
  ai <- sim_assign(relate = sim_relate(), estimator = "r", threshold = 0.25)
  expect_equal(nrow(ai), 11)
  expect_named(ai, c("Individual", "Genetic_Group", "Group", "Status", "n_samples", "share", "split",
                     "wide_ranging", "low_data", "change_period"))
  expect_equal(nrow(attr(ai, "samples")), nrow(sim_samples()))
  expect_error(assign_individuals(sim_samples(), "nope", h = 2000), "not found")
})

# A resident with one sample in a neighbouring territory is not wide-ranging; an animal split between two places is
test_that("one foray doesn't make a resident wide-ranging, two places do", {
  set.seed(3)
  pt <- function(id, n, x, sd = 1500) data.frame(Individual = id, X = x + stats::rnorm(n, 0, sd), Y = stats::rnorm(n, 0, sd))
  d <- rbind(pt("R", 12, 0), pt("R", 1, 25000, sd = 0),                   # resident + one sample 25 km away
             pt("M", 6, 0), pt("M", 6, 25000),                             # half here, half there
             pt("F", 2, 0), pt("F", 1, 25000, sd = 0))                    # 3 samples, one far: nothing left out
  ud <- space_use(sf::st_as_sf(d, coords = c("X", "Y"), crs = 32632), h = 2000)
  w <- .wide_ranging(ud, wide = 3)
  expect_false(w[["R"]])
  expect_true(w[["M"]])
  expect_true(w[["F"]])
  expect_true(.wide_ranging(ud, wide = 3, trim = 0)[["R"]])               # the plain radius of gyration flags it
})

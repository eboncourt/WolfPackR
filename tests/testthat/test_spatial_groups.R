test_that("a purely spatial analysis finds the two areas", {
  sg <- spatial_groups(sim_samples(), h = 2000)
  p <- stats::setNames(sg$Pack, sg$Individual)
  expect_equal(length(unique(stats::na.omit(p))), 2)
  expect_equal(length(unique(p[c("A1", "A2", "A3", "U")])), 1)
  expect_equal(length(unique(p[c("B1", "B2", "B3", "D")])), 1)
  expect_false(p[["A1"]] == p[["B1"]])
  expect_true(is.na(p[["T"]]))                                        # scattered far away: no pack
})

test_that("a group spread over two areas is split, and the packs are named after the group", {
  obs <- sim_samples()
  obs$g <- ifelse(obs$Individual == "T", NA, "fam")
  sg <- spatial_groups(obs, "g", h = 2000)
  expect_setequal(stats::na.omit(unique(sg$Pack)), c("fam.1", "fam.2"))
})

test_that("two subgroups overlapping a lot are merged back", {
  obs <- sim_samples()
  obs <- obs[obs$Individual %in% c("A1", "A2", "A3", "U"), ]
  obs$g <- "A"
  expect_equal(unique(spatial_groups(obs, "g", h = 2000)$Pack), "A")
  expect_true(length(unique(spatial_groups(obs, "g", h = 2000, merge = 1, min_affinity = 0.99)$Pack)) >= 1)
})

test_that("individuals without a group are pooled and can make a pack of their own", {
  obs <- sim_samples()
  obs$g <- ifelse(obs$Individual %in% c("A1", "A2", "A3", "U"), "A", NA)
  sg <- spatial_groups(obs, "g", h = 2000)
  p <- stats::setNames(sg$Pack, sg$Individual)
  expect_equal(p[["B1"]], "new")
  expect_equal(p[["A1"]], "A")
})

test_that("packs smaller than min_group_size are dropped", {
  obs <- sim_samples()
  obs$g <- ifelse(obs$Individual == "A1", "solo", "rest")
  sg <- spatial_groups(obs, "g", h = 2000)
  expect_true(is.na(sg$Pack[sg$Individual == "A1"]))
  expect_equal(sg$Pack[sg$Individual == "A1"], NA_character_)
  sg2 <- spatial_groups(obs, "g", h = 2000, min_group_size = 1)
  expect_equal(sg2$Pack[sg2$Individual == "A1"], "solo")
})

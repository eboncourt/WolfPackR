test_that("space_use gives one UD per individual, each summing to 1", {
  ud <- space_use(sim_samples(), h = 2000)
  expect_s3_class(ud, "wolfpack_ud")
  expect_equal(ncol(ud$U), 11)
  expect_equal(unname(colSums(ud$U)), rep(1, 11))
  expect_equal(ud$h, 2000)
  expect_equal(ud$res, 500)
  expect_equal(unname(ud$n["L"]), 1)                       # a single sample still gives a UD
})

test_that("space_use picks a bandwidth from the well-sampled individuals when h is missing", {
  expect_message(ud <- space_use(sim_samples()), "Bandwidth")
  expect_gt(ud$h, 500); expect_lt(ud$h, 5000)
})

test_that("lon/lat data are projected to UTM", {
  ud <- suppressMessages(space_use(pkg_data("samples_Carpathians")))
  expect_false(isTRUE(ud$crs$IsGeographic))
  expect_equal(ud$crs$epsg, 32635)                          # Carpathians: UTM 35N
})

test_that("space_use checks its input", {
  expect_error(space_use(data.frame(x = 1)), "sf")
  s <- sim_samples(); names(s)[1] <- "ID"
  expect_error(space_use(s, h = 2000), "Individual")
})

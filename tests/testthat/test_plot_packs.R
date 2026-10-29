obs_plot <- function() {
  obs <- sim_samples()
  obs$Pack <- ifelse(obs$Individual %in% c("A1", "A2", "A3"), "A", ifelse(obs$Individual %in% c("B1", "B2", "B3"), "B", NA))
  obs
}
legend_colors <- function(map) {
  calls <- map$x$calls
  calls[[which(sapply(calls, `[[`, "method") == "addLegend")]]$args[[1]]$colors
}

# Regression test: plot_packs used to call library(leaflet / sf / RColorBrewer),
# which attached those packages to the user's session (kept first in the file)
test_that("plot_packs leaves the search path alone", {
  attached <- c("package:leaflet", "package:RColorBrewer")
  skip_if(any(attached %in% search()), "leaflet or RColorBrewer already attached")
  plot_packs(obs_plot(), h = 2000)
  expect_false(any(attached %in% search()))
})

test_that("plot_packs returns a leaflet map, with a pack column given as text or bare name", {
  expect_s3_class(plot_packs(obs_plot(), "Pack", h = 2000), "leaflet")
  expect_s3_class(plot_packs(obs_plot(), Pack, h = 2000), "leaflet")
})

test_that("plot_packs checks its input", {
  expect_error(plot_packs(obs_plot()[0, ], h = 2000), "empty")
  expect_error(plot_packs(data.frame(), h = 2000), "sf")
  expect_error(plot_packs(obs_plot(), "nope", h = 2000), "doesn't exist")
})

test_that("plot_packs draws a territory and a core area per pack, and a legend entry for individuals without a pack", {
  map <- plot_packs(obs_plot(), h = 2000)
  methods <- sapply(map$x$calls, `[[`, "method")
  expect_equal(sum(methods == "addPolygons"), 4)                        # 2 packs x (95% + 50%)
  expect_length(legend_colors(map), 3)                                  # A, B, no pack
})

test_that("plot_packs gives a colour to every pack, even with more than 9 packs", {
  obs <- obs_plot()
  obs$Pack <- paste0("P", as.integer(factor(obs$Individual)))           # 11 packs of one individual
  map <- plot_packs(obs, h = 2000, levels = 0.95)
  expect_false(any(is.na(legend_colors(map))))
  expect_length(unique(legend_colors(map)), 11)
})

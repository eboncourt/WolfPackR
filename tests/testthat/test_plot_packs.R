# Prepare data
obs <- toy_obs()

legend_colors <- function(map) {
  calls <- map$x$calls
  calls[[which(sapply(calls, `[[`, "method") == "addLegend")]]$args[[1]]$colors
}

# Regression test: plot_packs used to call library(leaflet / sf / RColorBrewer),
# which attached those packages to the user's session
# (kept first in the file, before any other call to plot_packs)
test_that("plot_packs leaves the search path alone", {
  attached <- c("package:leaflet", "package:RColorBrewer")
  skip_if(any(attached %in% search()), "leaflet or RColorBrewer already attached")
  plot_packs(obs, group)
  expect_false(any(attached %in% search()))
})

test_that("plot_packs returns a leaflet object", {
  expect_s3_class(plot_packs(obs, group), "leaflet")
})

test_that("plot_packs returns an error if the input is empty", {
  expect_error(plot_packs(obs[0, ], group), "empty")
})

test_that("plot_packs returns an error if the input is not an sf object", {
  expect_error(plot_packs(data.frame(), group), "sf")
})

test_that("plot_packs can deal with lone individuals", {
  obs_lone <- obs
  obs_lone$group[obs_lone$group == 1] <- "Lone Individual"
  expect_s3_class(suppressWarnings(plot_packs(obs_lone, group)), "leaflet")
})

test_that("plot_packs gives a colour to every pack, even with more than 9 packs", {
  obs_many <- obs
  obs_many$group <- paste0("P", seq_len(nrow(obs_many)) %% 12)  # 12 packs
  map <- plot_packs(obs_many, group)
  expect_false(any(is.na(legend_colors(map))))
  expect_length(unique(legend_colors(map)), 13)                 # 12 packs + lone individuals
})

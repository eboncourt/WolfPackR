test_that("the old MCP functions and ugly_ducklings still work, with a deprecation warning", {
  pts <- sf::st_as_sf(data.frame(X = c(0, 1, 1, 0), Y = c(0, 0, 1, 1)), coords = c("X", "Y"))
  expect_warning(m <- mcp_sf(pts), class = "deprecatedWarning")
  expect_warning(calculate_overlap(m, m), class = "deprecatedWarning")
  obs <- merge(sim_samples(), genetic_groups(sim_relate(), "r", 0.25, method = "components"), by = "Individual")
  expect_warning(ud <- suppressMessages(ugly_ducklings(obs, group)), class = "deprecatedWarning")
  expect_true(all(c("new_group", "Pack", "Status") %in% names(ud)))
  expect_true(all(grepl("^Lone_", ud$new_group[is.na(ud$Group)])))
})

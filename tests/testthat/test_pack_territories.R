test_that("pack_territories gives the isopleths of each pack, the same as isopleth()", {
  ud <- space_use(sim_samples(), h = 2000)
  packs <- list(A = c("A1", "A2", "A3"), B = c("B1", "B2", "B3"))
  ter <- pack_territories(ud, packs)
  expect_s3_class(ter, "sf")
  expect_equal(ter$Pack, c("A", "A", "B", "B"))
  expect_equal(ter$level, c(0.95, 0.5, 0.95, 0.5))
  expect_equal(ter$area_km2[1:2], isopleth(ud, packs$A)$area_km2)
  expect_equal(sf::st_crs(ter), ud$crs)
  expect_error(pack_territories(ud, list(c("A1", "A2"))), "named list")
  expect_error(pack_territories(ud, list(A = "nobody")), "No member of A")
})

test_that("exclusive territories don't overlap, and each core stays inside its own territory", {
  ud <- space_use(sim_samples(), h = 2000)
  packs <- list(A = c("A1", "A2"), U = c("A3", "U"), B = c("B1", "B2", "B3"))   # A and U share the same area
  plain <- pack_territories(ud, packs)
  excl <- pack_territories(ud, packs, exclusive = TRUE)
  t95 <- function(x, p) x[x$Pack == p & x$level == 0.95, ]
  overlap <- function(x) as.numeric(sf::st_area(sf::st_intersection(sf::st_geometry(t95(x, "A")), sf::st_geometry(t95(x, "U")))))
  expect_gt(sum(overlap(plain)), 10e6)                                   # plain isopleths: tens of km2 in common
  expect_lt(sum(overlap(excl)), 1)                                       # exclusive: only a shared boundary
  m <- merge(sf::st_drop_geometry(plain), sf::st_drop_geometry(excl), by = c("Pack", "level"))
  expect_true(all(m$area_km2.y <= m$area_km2.x + 1e-9))                  # a pack can only lose cells
  expect_equal(t95(excl, "B")$area_km2, t95(plain, "B")$area_km2)        # B is far from the others: unchanged
  for (p in unique(excl$Pack)) {
    core <- excl[excl$Pack == p & excl$level == 0.5, ]
    if (nrow(core)) expect_true(sf::st_covered_by(core, sf::st_buffer(t95(excl, p), 1), sparse = FALSE)[1, 1])
  }
})

test_that("summarise_packs and plot_packs take exclusive", {
  obs <- sim_samples()
  obs$Pack <- ifelse(obs$Individual %in% c("A1", "A2"), "A", ifelse(obs$Individual %in% c("A3", "U"), "U",
                                                                    ifelse(obs$Individual %in% c("B1", "B2", "B3"), "B", NA)))
  ud <- space_use(obs, h = 2000)
  s0 <- summarise_packs(obs, "Pack", ud = ud)
  s1 <- summarise_packs(obs, "Pack", ud = ud, exclusive = TRUE)
  expect_true(all(s1$area95_km2 <= s0$area95_km2))
  expect_lt(s1$area95_km2[s1$Pack == "A"], s0$area95_km2[s0$Pack == "A"])
  expect_equal(s1$shared_pct, s0$shared_pct)                             # contested share, on the plain isopleths
  expect_equal(s0[, setdiff(names(s0), c("area95_km2", "area50_km2"))], s1[, setdiff(names(s1), c("area95_km2", "area50_km2"))])
  m <- plot_packs(obs, "Pack", ud = ud, exclusive = TRUE)
  expect_s3_class(m, "leaflet")
})

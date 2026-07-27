test_that("post-smoothing CI widths satisfy every accumulated floor", {
  # Mimic edge smoothing that narrowed three intervals below either the
  # minimum-relative-width floor or the stronger local-MAD floor.
  lo <- c(-0.10, -0.20, -0.30, -0.40, -0.50)
  hi <- c( 0.10,  0.20,  0.30,  0.40,  0.50)
  min_relative_floor <- rep(0.50, length(lo))
  local_mad_floor <- c(0.30, 0.70, 0.60, 1.10, 0.80)
  floor_width <- pmax(min_relative_floor, local_mad_floor)
  center_before <- (lo + hi) / 2

  out <- geneSCOPE:::.reapply_ci_width_floor(lo, hi, floor_width)

  expect_true(all((out$hi - out$lo) >= floor_width))
  expect_equal((out$lo + out$hi) / 2, center_before, tolerance = 0)
  expect_true(out$reapplied)
  expect_identical(out$reapplied_count, 3L)
  expect_identical(out$active_count, 5L)
  expect_equal(out$max_shortfall_before, 0.30, tolerance = 1e-15)
})

test_that("CI floor reapplication is an exact no-op when smoothing is disabled or safe", {
  # This is the no-smoothing path: the pre-smoothed intervals already satisfy
  # their floors and must not move, widen, or shift centre.
  lo <- c(-1.0, -0.5, 0.0)
  hi <- c( 1.0,  0.5, 2.0)
  floor_width <- c(1.0, 1.0, 0.0)

  out <- geneSCOPE:::.reapply_ci_width_floor(lo, hi, floor_width)

  expect_identical(out$lo, lo)
  expect_identical(out$hi, hi)
  expect_false(out$reapplied)
  expect_identical(out$reapplied_count, 0L)
  expect_equal(out$max_shortfall_before, 0, tolerance = 0)
})

test_that("CI floor helper rejects mismatched pointwise constraints", {
  expect_error(
    geneSCOPE:::.reapply_ci_width_floor(c(0, 1), c(1, 2), 1),
    "identical lengths"
  )
})

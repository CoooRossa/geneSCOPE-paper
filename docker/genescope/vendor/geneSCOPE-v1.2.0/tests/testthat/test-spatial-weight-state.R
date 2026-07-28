test_that("disabled weight outputs remove stale layer state", {
  obj <- make_lee_scope()
  obj@grid$grid1$listw <- structure(list(stale = TRUE), class = "listw")

  out <- geneSCOPE:::.compute_weights(
    obj, grid_name = "grid1", style = "B", topology = "queen",
    store_mat = FALSE, store_listw = FALSE, ncores = 1L,
    verbose = FALSE
  )

  expect_false("W" %in% names(out@grid$grid1))
  expect_false("listw" %in% names(out@grid$grid1))
  expect_false(out@grid$grid1$weight_provenance$store_mat)
  expect_false(out@grid$grid1$weight_provenance$store_listw)
  expect_identical(out@grid$grid1$weight_provenance$verification_source, "none")
})

test_that("binary and row-standardized weights retain matching semantics", {
  binary <- geneSCOPE:::.compute_weights(
    make_lee_scope(), grid_name = "grid1", style = "B", topology = "queen",
    store_mat = TRUE, store_listw = FALSE, ncores = 1L,
    verbose = FALSE
  )
  W_binary <- binary@grid$grid1$W
  p_binary <- binary@grid$grid1$weight_provenance
  expect_identical(attr(W_binary, "weight_style", exact = TRUE), "B")
  expect_true(all(W_binary@x == 1))
  expect_false(p_binary$row_normalized)
  expect_true(p_binary$binary_verified)
  expect_false("listw" %in% names(binary@grid$grid1))

  standardized <- geneSCOPE:::.compute_weights(
    binary, grid_name = "grid1", style = "W", topology = "queen",
    store_mat = TRUE, store_listw = FALSE, ncores = 1L,
    verbose = FALSE
  )
  W_standardized <- standardized@grid$grid1$W
  p_standardized <- standardized@grid$grid1$weight_provenance
  expect_identical(attr(W_standardized, "weight_style", exact = TRUE), "W")
  expect_true(p_standardized$row_normalized)
  expect_true(p_standardized$row_normalized_verified)
  expect_true(grepl("listw_B_omp", p_standardized$matrix_backend, fixed = TRUE))
  expect_equal(as.numeric(Matrix::rowSums(W_standardized)), rep(1, 4), tolerance = 1e-12)
  expect_false("listw" %in% names(standardized@grid$grid1))
  expect_identical(
    attr(W_standardized, "weight_provenance", exact = TRUE),
    p_standardized
  )
})

test_that("kernel weights are explicitly W and honor output storage flags", {
  kernel <- geneSCOPE:::.compute_weights(
    make_lee_scope(), grid_name = "grid1", style = "kernel_gaussian",
    topology = "queen", kernel_radius = 2L, kernel_sigma = 1,
    store_mat = TRUE, store_listw = FALSE, ncores = 1L,
    verbose = FALSE
  )
  W_kernel <- kernel@grid$grid1$W
  p_kernel <- kernel@grid$grid1$weight_provenance
  expect_identical(attr(W_kernel, "weight_style", exact = TRUE), "W")
  expect_identical(p_kernel$requested_style, "kernel_gaussian")
  expect_identical(p_kernel$weight_style, "W")
  expect_true(p_kernel$row_normalized_verified)
  expect_true(grepl("grid_weights_kernel_rect_omp", p_kernel$backend, fixed = TRUE))
  expect_equal(as.numeric(Matrix::rowSums(W_kernel)), rep(1, 4), tolerance = 1e-12)

  cleared <- geneSCOPE:::.compute_weights(
    kernel, grid_name = "grid1", style = "kernel_flat", topology = "queen",
    kernel_radius = 2L, store_mat = FALSE, store_listw = FALSE,
    ncores = 1L, verbose = FALSE
  )
  expect_false("W" %in% names(cleared@grid$grid1))
  expect_false("listw" %in% names(cleared@grid$grid1))
  expect_false(cleared@grid$grid1$weight_provenance$store_mat)
  expect_false(cleared@grid$grid1$weight_provenance$store_listw)
  expect_identical(cleared@grid$grid1$weight_provenance$weight_style, "W")
  expect_true(cleared@grid$grid1$weight_provenance$row_normalized_verified)
})

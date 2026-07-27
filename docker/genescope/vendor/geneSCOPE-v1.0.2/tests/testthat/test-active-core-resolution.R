test_that("computeL passes the safely resolved explicit request to the Lee backend", {
  detector_state <- new.env(parent = emptyenv())
  detector_state$value <- NA_integer_
  backend_calls <- integer()
  original_lee_l_cache <- geneSCOPE:::.lee_l_cache

  testthat::local_mocked_bindings(
    detectCores = function(logical = TRUE) detector_state$value,
    .lee_l_cache = function(Xz, W, n_threads = 1L) {
      backend_calls <<- c(backend_calls, as.integer(n_threads))
      original_lee_l_cache(Xz, W, n_threads = 1L)
    },
    .package = "geneSCOPE"
  )

  out_na <- computeL(
    make_lee_scope(), grid_name = "grid1", perms = 0L, ncores = 16L,
    cache_inputs = FALSE, verbose = FALSE
  )
  detector_state$value <- 8L
  out_eight <- computeL(
    make_lee_scope(), grid_name = "grid1", perms = 0L, ncores = 16L,
    cache_inputs = FALSE, verbose = FALSE
  )

  expect_identical(backend_calls, c(16L, 8L))
  expect_identical(out_na@stats$grid1$LeeStats_Xz$meta$ncores, 16L)
  expect_identical(out_eight@stats$grid1$LeeStats_Xz$meta$ncores, 8L)
})

test_that("computeCorrelation passes the safely resolved request to its native backend", {
  detector_state <- new.env(parent = emptyenv())
  detector_state$value <- NA_integer_
  backend_calls <- integer()

  obj <- make_lee_scope()
  counts <- matrix(
    c(1, 2, 4, 3, 4, 1, 3, 2, 2, 5, 1, 4),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), paste0("c", 1:4))
  )
  obj@cells$logCPM <- methods::as(counts, "dgCMatrix")

  testthat::local_mocked_bindings(
    detectCores = function(logical = TRUE) detector_state$value,
    .pearson_cor = function(X, bs, n_threads) {
      backend_calls <<- c(backend_calls, as.integer(n_threads))
      stats::cor(X)
    },
    .package = "geneSCOPE"
  )

  out_na <- computeCorrelation(
    obj, level = "cell", layer = "logCPM", ncores = 16L,
    compute_fdr = FALSE, use_bigmemory = FALSE, verbose = FALSE
  )
  detector_state$value <- 8L
  out_eight <- computeCorrelation(
    obj, level = "cell", layer = "logCPM", ncores = 16L,
    compute_fdr = FALSE, use_bigmemory = FALSE, verbose = FALSE
  )
  detector_state$value <- NA_integer_
  out_invalid <- computeCorrelation(
    obj, level = "cell", layer = "logCPM", ncores = NA_integer_,
    compute_fdr = FALSE, use_bigmemory = FALSE, verbose = FALSE
  )

  expect_identical(backend_calls, c(16L, 8L, 1L))
  expect_equal(
    out_na@cells$.pearson_cor,
    out_eight@cells$.pearson_cor,
    tolerance = 0
  )
  expect_equal(
    out_na@cells$.pearson_cor,
    out_invalid@cells$.pearson_cor,
    tolerance = 0
  )
})

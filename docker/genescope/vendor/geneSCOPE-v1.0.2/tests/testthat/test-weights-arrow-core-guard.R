test_that("computeWeights resolves explicit and default core requests safely", {
  detector_state <- new.env(parent = emptyenv())
  detector_state$mode <- "na"
  backend_calls <- integer()
  original_builder <- geneSCOPE:::.spw_build_neighbourhood

  testthat::local_mocked_bindings(
    detectCores = function(logical = TRUE) {
      if (identical(detector_state$mode, "error")) {
        stop("simulated detector failure")
      }
      if (identical(detector_state$mode, "na")) NA_integer_ else 8L
    },
    .spw_build_neighbourhood = function(meta,
                                        topo_info,
                                        ncores = 1L,
                                        verbose = FALSE,
                                        parent = "computeWeights",
                                        step = "S03") {
      backend_calls <<- c(backend_calls, as.integer(ncores))
      # Preserve the production builder and its W mathematics, but avoid
      # launching many processes for this four-grid fixture.
      original_builder(
        meta = meta,
        topo_info = topo_info,
        ncores = 1L,
        verbose = verbose,
        parent = parent,
        step = step
      )
    },
    .package = "geneSCOPE"
  )

  out_na <- computeWeights(
    make_lee_scope(), grid_name = "grid1", topology = "queen",
    store_listw = FALSE, ncores = 16L, verbose = FALSE
  )
  detector_state$mode <- "error"
  out_error <- computeWeights(
    make_lee_scope(), grid_name = "grid1", topology = "queen",
    store_listw = FALSE, ncores = 16L, verbose = FALSE
  )
  detector_state$mode <- "eight"
  out_eight <- computeWeights(
    make_lee_scope(), grid_name = "grid1", topology = "queen",
    store_listw = FALSE, ncores = 16L, verbose = FALSE
  )
  detector_state$mode <- "error"
  out_default_error <- computeWeights(
    make_lee_scope(), grid_name = "grid1", topology = "queen",
    store_listw = FALSE, verbose = FALSE
  )

  expect_identical(backend_calls, c(16L, 16L, 8L, 1L))
  expect_identical(out_na@grid$grid1$W, out_error@grid$grid1$W)
  expect_identical(out_na@grid$grid1$W, out_eight@grid$grid1$W)
  expect_identical(out_na@grid$grid1$W, out_default_error@grid$grid1$W)
  expect_null(formals(computeWeights)$ncores)
  expect_null(formals(geneSCOPE:::.compute_weights)$ncores)
})

test_that("spatial-weight guard never performs an unsafe independent probe", {
  testthat::local_mocked_bindings(
    detectCores = function(logical = TRUE) stop("simulated detector failure"),
    .package = "geneSCOPE"
  )

  guard <- NULL
  expect_silent({
    guard <- geneSCOPE:::.spw_thread_guard_v2(
      requested_ncores = 16L,
      builder_threads = 1L
    )
  })
  on.exit(if (!is.null(guard)) guard$restore(), add = TRUE)
  expect_identical(guard$config$openmp_threads, 1L)
  expect_identical(guard$config$arrow_threads, 1L)
})

test_that("Arrow core resolution is detector-safe and never expands one thread", {
  detector_na <- function(logical = TRUE) NA_integer_
  detector_error <- function(logical = TRUE) stop("simulated detector failure")
  detector_eight <- function(logical = TRUE) 8L

  expect_identical(
    geneSCOPE:::.deterministic_arrow_threads(16L, detector = detector_na),
    16L
  )
  expect_identical(
    geneSCOPE:::.deterministic_arrow_threads(16L, detector = detector_error),
    16L
  )
  expect_identical(
    geneSCOPE:::.deterministic_arrow_threads(16L, detector = detector_eight),
    8L
  )
  expect_identical(
    geneSCOPE:::.deterministic_arrow_threads(1L, detector = detector_eight),
    1L
  )
  expect_identical(
    geneSCOPE:::.deterministic_arrow_threads(NULL, detector = detector_error),
    1L
  )
})

test_that("Arrow setters receive the resolved count without a second detector", {
  thread_vars <- c(
    "OMP_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS", "GOTO_NUM_THREADS", "ATLAS_NUM_THREADS",
    "LAPACK_NUM_THREADS"
  )
  old_thread_env <- Sys.getenv(thread_vars, unset = NA_character_)
  on.exit({
    present <- !is.na(old_thread_env)
    if (any(present)) {
      do.call(Sys.setenv, as.list(old_thread_env[present]))
    }
    if (any(!present)) Sys.unsetenv(names(old_thread_env)[!present])
  }, add = TRUE)
  do.call(
    Sys.setenv,
    as.list(stats::setNames(rep("1", length(thread_vars)), thread_vars))
  )

  calls <- new.env(parent = emptyenv())
  calls$cpu <- integer()
  calls$io <- integer()
  fake_arrow <- new.env(parent = emptyenv())
  fake_arrow$set_cpu_count <- function(value) {
    calls$cpu <- c(calls$cpu, as.integer(value))
  }
  fake_arrow$set_io_thread_count <- function(value) {
    calls$io <- c(calls$io, as.integer(value))
  }
  detector_error <- function(logical = TRUE) stop("must not be called")

  applied_16 <- geneSCOPE:::.apply_thread_config(
    list(
      openmp_threads = 16L,
      r_threads = 1L,
      blas_threads = 1L,
      arrow_threads = 16L
    ),
    arrow_namespace = fake_arrow,
    detector = detector_error
  )
  applied_1 <- geneSCOPE:::.apply_thread_config(
    list(
      openmp_threads = 1L,
      r_threads = 1L,
      blas_threads = 1L,
      arrow_threads = 1L
    ),
    arrow_namespace = fake_arrow,
    detector = detector_error
  )

  expect_identical(applied_16, 16L)
  expect_identical(applied_1, 1L)
  expect_identical(calls$cpu, c(16L, 1L))
  expect_identical(calls$io, c(16L, 1L))

  # A failing CPU setter must not prevent the independent I/O setter.
  calls$io <- integer()
  fake_arrow$set_cpu_count <- function(value) stop("simulated setter failure")
  expect_silent(
    geneSCOPE:::.apply_thread_config(
      list(openmp_threads = 1L, r_threads = 1L, blas_threads = 1L),
      arrow_namespace = fake_arrow,
      detector = detector_error
    )
  )
  expect_identical(calls$io, 1L)
})

test_that("fixed-r native Delta uses the observed Pearson reference", {
  set.seed(11)
  X <- scale(matrix(rnorm(24), 6, 4))
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6),
    j = c(2, 1, 3, 2, 4, 3, 5, 4, 6, 5),
    x = 1, dims = c(6, 6)
  )
  idx <- cbind(c(1, 0, 2, 4, 3, 5), c(2, 1, 0, 5, 4, 3))
  pairs <- matrix(c(0, 1, 2, 3), ncol = 2, byrow = TRUE)
  pearson_ref <- c(-0.40, 0.25)
  S2 <- sum(Matrix::rowSums(W)^2)
  lee_pair <- function(ii, pair) {
    Xp <- X[ii + 1L, , drop = FALSE]
    WX <- W %*% Xp
    nrow(Xp) / S2 * sum(WX[, pair[1] + 1L] * WX[, pair[2] + 1L]) /
      sqrt(sum(Xp[, pair[1] + 1L]^2) * sum(Xp[, pair[2] + 1L]^2))
  }
  observed_L <- vapply(seq_len(nrow(pairs)), function(k) lee_pair(0:5, pairs[k, ]), numeric(1))
  delta_ref <- observed_L - pearson_ref
  oracle <- vapply(seq_len(nrow(pairs)), function(k) {
    sum(vapply(seq_len(ncol(idx)), function(b) {
      abs(lee_pair(idx[, b], pairs[k, ]) - pearson_ref[k]) >= abs(delta_ref[k])
    }, logical(1)))
  }, numeric(1))

  Wr <- methods::as(W, "RsparseMatrix")
  got <- geneSCOPE:::delta_l_fixed_r_perm_csr(
    X, as.integer(Wr@j), Wr@x, as.integer(Wr@p), idx,
    pairs, delta_ref, pearson_ref, 1L
  )
  expect_equal(as.numeric(got), oracle, tolerance = 0)
})

test_that("observed native Lee matches an independent canonical S2 oracle", {
  set.seed(20260728)
  X <- scale(matrix(stats::rnorm(35), nrow = 7L, ncol = 5L))
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7),
    j = c(2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6),
    x = c(1, 1, 0.5, 0.5, 1.5, 1.5, 0.75, 0.75, 1.25, 1.25, 1, 1),
    dims = c(7L, 7L)
  )

  WX <- as.matrix(W %*% X)
  S2 <- sum(as.numeric(Matrix::rowSums(W))^2)
  norms <- sqrt(colSums(X^2))
  oracle <- nrow(X) / S2 * crossprod(WX) / outer(norms, norms)
  native <- geneSCOPE:::lee_L(X, W, 1L)

  expect_lt(max(abs(native - oracle)), 1e-12)
})

test_that("legacy Delta native ABIs keep their arities and reject invalid permutations", {
  X <- scale(matrix(seq_len(8), 4, 2))
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1, dims = c(4, 4)
  )
  bad_idx <- matrix(c(0L, 0L, 2L, 3L), ncol = 1)
  pairs <- matrix(c(0L, 1L), nrow = 1)

  expect_error(
    geneSCOPE:::delta_lr_perm(
      X, W, bad_idx, pairs, 0, 1L, 1000L
    ),
    "must be a bijection"
  )
  expect_error(
    geneSCOPE:::delta_lr_perm(
      X, W, matrix(0:3, ncol = 1), pairs, 0, 0L, 1000L
    ),
    "n_threads must be at least 1"
  )

  Wr <- methods::as(W, "RsparseMatrix")
  expect_error(
    geneSCOPE:::delta_lr_perm_csr(
      X, as.integer(Wr@j), Wr@x, as.integer(Wr@p), bad_idx,
      pairs, 0, 1L, FALSE
    ),
    "must be a bijection"
  )
})

test_that("Lee cache fingerprint changes with every determining input", {
  X <- matrix(1:12, 4, 3, dimnames = list(paste0("s", 1:4), paste0("g", 1:3)))
  W <- Matrix::Diagonal(4)
  dimnames(W) <- list(rownames(X), rownames(X))
  attr(W, "weight_style") <- "B"
  gi <- data.frame(grid_id = rownames(X), gx = 1:4, gy = 1L)
  fp <- geneSCOPE:::.lee_input_fingerprint(X, W, gi, "Xz", 8L)

  X2 <- X; X2[1, 1] <- X2[1, 1] + 1
  expect_false(identical(fp, geneSCOPE:::.lee_input_fingerprint(X2, W, gi, "Xz", 8L)))
  W2 <- W; W2[1, 1] <- 2
  expect_false(identical(fp, geneSCOPE:::.lee_input_fingerprint(X, W2, gi, "Xz", 8L)))
  gi2 <- gi; gi2$gx[1] <- 9
  expect_false(identical(fp, geneSCOPE:::.lee_input_fingerprint(X, W, gi2, "Xz", 8L)))
  W3 <- W; attr(W3, "weight_style") <- "W"
  expect_false(identical(fp, geneSCOPE:::.lee_input_fingerprint(X, W3, gi, "Xz", 8L)))
})

test_that("perms=0 stores observed L without entering inference", {
  obj <- make_lee_scope()
  out <- computeL(
    obj, grid_name = "grid1", perms = 0, ncores = 1,
    use_bigmemory = FALSE, cache_inputs = FALSE, verbose = FALSE
  )
  ans <- out@stats$grid1$LeeStats_Xz
  expect_true(is.matrix(ans$L))
  expect_null(ans$P)
  expect_null(ans$Z)
  expect_null(ans$FDR)
  expect_identical(ans$meta$FDR_main_method, "unavailable (perms=0)")
})

test_that("unsupported legacy and rectangular permutation modes fail early", {
  obj <- make_lee_scope()
  expect_error(
    computeL(obj, grid_name = "grid1", legacy_formula = TRUE, verbose = FALSE),
    "legacy_formula=TRUE is no longer supported"
  )
  expect_error(
    computeL(obj, grid_name = "grid1", genes = c("g1", "g2"), within = FALSE,
      perms = 1, verbose = FALSE),
    "Permutation inference for within=FALSE"
  )
})

test_that("large-result pseudo-streaming fails closed before allocation", {
  obj <- make_lee_scope()
  expect_error(
    geneSCOPE:::.compute_lee_l(
      obj, grid_name = "grid1", norm_layer = "Xz",
      use_bigmemory = TRUE, mem_limit_GB = 0,
      cache_inputs = FALSE, ncores = 1
    ),
    "disabled.*RAM-backed"
  )
  expect_error(
    geneSCOPE:::.compute_lee_l(
      obj, grid_name = "grid1", norm_layer = "Xz",
      use_bigmemory = FALSE, mem_limit_GB = 0,
      cache_inputs = FALSE, ncores = 1
    ),
    "exceeds mem_limit_GB"
  )
})

test_that("safe in-memory Lee mode is the public and internal default", {
  expect_identical(formals(computeL)$use_bigmemory, FALSE)
  expect_identical(formals(geneSCOPE:::.compute_l)$use_bigmemory, FALSE)
  expect_identical(formals(geneSCOPE:::.compute_lee_l)$use_bigmemory, FALSE)
})

test_that("within-gene selection is computed on the selected native matrix", {
  obj <- make_lee_scope()
  full <- geneSCOPE:::.compute_lee_l(
    obj, grid_name = "grid1", norm_layer = "Xz", genes = NULL,
    within = TRUE, use_bigmemory = FALSE, cache_inputs = FALSE, ncores = 1
  )
  sub <- geneSCOPE:::.compute_lee_l(
    obj, grid_name = "grid1", norm_layer = "Xz", genes = c("g3", "g1"),
    within = TRUE, use_bigmemory = FALSE, cache_inputs = FALSE, ncores = 1
  )
  expect_identical(dim(sub$Lmat), c(2L, 2L))
  expect_identical(colnames(sub$X_used), c("g3", "g1"))
  expect_equal(
    unname(sub$Lmat),
    unname(full$Lmat[c("g3", "g1"), c("g3", "g1"), drop = FALSE]),
    tolerance = 1e-15
  )
})

test_that("observed native Lee fails closed for non-positive S2", {
  X <- scale(matrix(seq_len(12), nrow = 4L, ncol = 3L))
  W0 <- Matrix::sparseMatrix(
    i = integer(), j = integer(), x = numeric(), dims = c(4L, 4L)
  )
  expect_error(geneSCOPE:::lee_L(X, W0, 1L), "non-positive Lee S2")
  expect_error(geneSCOPE:::lee_L_cache(X, W0, 1L), "non-positive Lee S2")
  expect_error(geneSCOPE:::lee_L_cols(X, W0, 0L, 1L), "non-positive Lee S2")
})

test_that("native Lee entry points validate dimensions, threads, and selected columns", {
  X <- scale(matrix(seq_len(12), nrow = 4L, ncol = 3L))
  W <- Matrix::Diagonal(4L)
  expect_error(geneSCOPE:::lee_L(X, W, 0L), "n_threads must be at least 1")
  expect_error(geneSCOPE:::lee_L(X, Matrix::Diagonal(3L), 1L), "W must be")
  expect_error(geneSCOPE:::lee_L_cols(X, W, integer(), 1L), "at least one")
  expect_error(geneSCOPE:::lee_L_cols(X, W, 3L, 1L), "out-of-range")

  L <- geneSCOPE:::lee_L_cache(X, W, 1L)
  idx <- matrix(0:3, ncol = 1L)
  expect_error(
    geneSCOPE:::lee_perm_block(X, W, idx, rep(1L, 4L), L[-1, -1], 1L),
    "L_ref must be"
  )
})

test_that("Lee and fixed-r permutation counts are thread-count invariant", {
  set.seed(20260722)
  X <- scale(matrix(rnorm(30), nrow = 6L, ncol = 5L))
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6),
    j = c(2, 1, 3, 2, 4, 3, 5, 4, 6, 5),
    x = 1, dims = c(6, 6)
  )
  idx <- cbind(
    c(1, 0, 2, 4, 3, 5),
    c(2, 1, 0, 5, 4, 3),
    c(0, 2, 1, 3, 5, 4)
  )
  L <- geneSCOPE:::lee_L_cache(X, W, 1L)
  one <- geneSCOPE:::lee_perm_block(X, W, idx, rep(1L, 6L), L, 1L)
  two <- geneSCOPE:::lee_perm_block(X, W, idx, rep(1L, 6L), L, 2L)
  expect_equal(two, one, tolerance = 0)

  pairs <- matrix(c(0L, 1L, 2L, 4L), ncol = 2L, byrow = TRUE)
  pearson_ref <- c(stats::cor(X[, 1], X[, 2]), stats::cor(X[, 3], X[, 5]))
  delta_ref <- c(L[1, 2], L[3, 5]) - pearson_ref
  fixed_one <- geneSCOPE:::delta_l_fixed_r_perm(
    X, W, idx, pairs, delta_ref, pearson_ref, 1L
  )
  fixed_two <- geneSCOPE:::delta_l_fixed_r_perm(
    X, W, idx, pairs, delta_ref, pearson_ref, 2L
  )
  expect_equal(fixed_two, fixed_one, tolerance = 0)
})

test_that("Lee permutation treats zero-norm null statistics as zero and rejects bad references", {
  X <- cbind(
    g1 = as.numeric(scale(c(1, 3, 2, 4))),
    zero = 0,
    g3 = as.numeric(scale(c(4, 1, 3, 2)))
  )
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1, dims = c(4L, 4L)
  )
  idx <- cbind(0:3, c(1L, 0L, 3L, 2L))
  L <- geneSCOPE:::lee_L_cache(X, W, 1L)
  expect_equal(L[2L, ], rep(0, 3L), tolerance = 0)
  expect_equal(L[, 2L], rep(0, 3L), tolerance = 0)

  counts <- geneSCOPE:::lee_perm_block(X, W, idx, rep(1L, 4L), L, 2L)
  expect_equal(counts[2L, ], rep(ncol(idx), 3L), tolerance = 0)
  expect_equal(counts[, 2L], rep(ncol(idx), 3L), tolerance = 0)

  bad_ref <- L
  bad_ref[1L, 2L] <- NA_real_
  expect_error(
    geneSCOPE:::lee_perm_block(X, W, idx, rep(1L, 4L), bad_ref, 1L),
    "L_ref contains a non-finite"
  )
})

test_that("normalization restores empty grids in authoritative order", {
  obj <- make_lee_scope(include_empty_grid = TRUE)
  out <- normalizeMoleculesInGrid(obj, grid_name = "grid1", keep_zero_grids = FALSE, verbose = FALSE)
  expect_identical(rownames(out@grid$grid1$Xz), out@grid$grid1$grid_info$grid_id)
  expect_equal(unname(out@grid$grid1$Xz["s5", ]), rep(0, 3), tolerance = 0)
})

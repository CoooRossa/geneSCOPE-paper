test_that("L-vs-r source fingerprints detect values and gene ordering", {
  common <- c("g2", "g1", "g3")
  L <- matrix(c(
    0, 0.3, 0.4,
    0.3, 0, 0.5,
    0.4, 0.5, 0
  ), 3, 3, dimnames = list(common, common))
  r <- matrix(c(
    1, 0.1, -0.2,
    0.1, 1, 0.25,
    -0.2, 0.25, 1
  ), 3, 3, dimnames = list(common, common))
  upstream <- list(schema = "lee_input_fingerprint_v2", data = list(hash = "abc"))

  first <- geneSCOPE:::.lvsr_curve_source_fingerprint(
    L, r, common, "cell", "LeeStats_Xz", upstream
  )
  changed_r <- r
  changed_r[1, 2] <- changed_r[2, 1] <- 0.11
  changed <- geneSCOPE:::.lvsr_curve_source_fingerprint(
    L, changed_r, common, "cell", "LeeStats_Xz", upstream
  )
  reordered_genes <- rev(common)
  reordered <- geneSCOPE:::.lvsr_curve_source_fingerprint(
    L[reordered_genes, reordered_genes],
    r[reordered_genes, reordered_genes],
    reordered_genes, "cell", "LeeStats_Xz", upstream
  )

  expect_identical(first$schema, "lvsr_curve_source_fingerprint_v1")
  expect_identical(first$pearson$level, "cell")
  expect_identical(first$common_genes$count, 3L)
  expect_match(first$pearson$aligned_matrix_sha256, "^[0-9a-f]{64}$")
  expect_false(identical(
    first$pearson$aligned_matrix_sha256,
    changed$pearson$aligned_matrix_sha256
  ))
  expect_false(identical(
    first$common_genes$order_sha256,
    reordered$common_genes$order_sha256
  ))
  expect_identical(
    first$common_genes$set_sha256,
    reordered$common_genes$set_sha256
  )
})

test_that("L-vs-r fingerprints ignore dimname component labels only", {
  common <- c("g2", "g1", "g3")
  L <- matrix(c(
    0, 0.3, 0.4,
    0.3, 0, 0.5,
    0.4, 0.5, 0
  ), 3, 3, dimnames = list(common, common))
  r <- matrix(c(
    1, 0.1, -0.2,
    0.1, 1, 0.25,
    -0.2, 0.25, 1
  ), 3, 3, dimnames = list(common, common))
  upstream <- list(schema = "lee_input_fingerprint_v2", data = list(hash = "abc"))

  canonical <- geneSCOPE:::.lvsr_curve_source_fingerprint(
    L, r, common, "cell", "LeeStats_Xz", upstream
  )

  L_labeled <- L
  dimnames(L_labeled) <- list(row = common, col = common)
  expect_identical(names(dimnames(L_labeled)), c("row", "col"))
  expect_null(names(dimnames(r)))

  mixed_labels <- geneSCOPE:::.lvsr_curve_source_fingerprint(
    L_labeled, r, common, "cell", "LeeStats_Xz", upstream
  )
  expect_identical(
    mixed_labels$lee$aligned_matrix_sha256,
    canonical$lee$aligned_matrix_sha256
  )
  expect_identical(
    mixed_labels$pearson$aligned_matrix_sha256,
    canonical$pearson$aligned_matrix_sha256
  )
  expect_identical(
    mixed_labels$unordered_pair_universe$aligned_L_and_Pearson_sha256,
    canonical$unordered_pair_universe$aligned_L_and_Pearson_sha256
  )

  L_bad_row <- L_labeled
  rownames(L_bad_row) <- rev(common)
  expect_error(
    geneSCOPE:::.lvsr_curve_source_fingerprint(
      L_bad_row, r, common, "cell", "LeeStats_Xz", upstream
    ),
    "dimnames must match common in order",
    fixed = TRUE
  )

  r_bad_col <- r
  colnames(r_bad_col) <- c(common[-1L], "not-a-common-gene")
  expect_error(
    geneSCOPE:::.lvsr_curve_source_fingerprint(
      L, r_bad_col, common, "cell", "LeeStats_Xz", upstream
    ),
    "dimnames must match common in order",
    fixed = TRUE
  )
})

test_that("stored L-vs-r metadata records complete curve and RNG provenance", {
  genes <- paste0("g", seq_len(8L))
  ut <- upper.tri(matrix(0, 8L, 8L), diag = FALSE)
  pear_values <- seq(-0.85, 0.85, length.out = sum(ut))
  L_values <- 0.18 + 0.32 * pear_values + 0.07 * sin(seq_along(pear_values))

  pearson <- matrix(0, 8L, 8L, dimnames = list(genes, genes))
  pearson[ut] <- pear_values
  pearson <- pearson + t(pearson)
  diag(pearson) <- 1
  lees_l <- matrix(0, 8L, 8L, dimnames = list(genes, genes))
  lees_l[ut] <- L_values
  lees_l <- lees_l + t(lees_l)

  lee_fp <- list(
    schema = "lee_input_fingerprint_v2",
    norm_layer = "Xz",
    data = list(schema = "lee_observed_data_v1", hash = "fixture")
  )
  obj <- methods::new(
    "scope_object",
    grid = list(grid1 = list(fixture = TRUE)),
    stats = list(grid1 = list(
      LeeStats_Xz = list(
        L = lees_l,
        meta = list(
          formula_id = "Lee2009_S2_v1",
          input_fingerprint = lee_fp
        )
      ),
      .pearson_cor = pearson
    )),
    coord = list(), cells = list(), density = list(),
    meta.data = data.frame(row.names = genes)
  )

  set.seed(20260723)
  out <- geneSCOPE:::.compute_l_vs_r_curve(
    obj,
    grid_name = "grid1",
    level = "grid",
    B = 20L,
    span = 0.55,
    deg = 1L,
    ncores = 2L,
    length_out = 25L,
    downsample = 0.8,
    n_strata = 6L,
    k_max = 20L,
    jitter_eps = 0,
    ci_method = "percentile",
    ci_adjust = "none",
    min_rel_width = 0.15,
    widen_span = 0.5,
    curve_name = "provenance_fixture",
    verbose = FALSE
  )
  curve <- out@stats$grid1$LeeStats_Xz$provenance_fixture
  meta <- curve$meta[[1L]]

  expect_identical(meta$schema, "lvsr_curve_meta_v2")
  expect_identical(meta$level, "grid")
  expect_identical(meta$pearson_level, "grid")
  expect_identical(meta$B, 20L)
  expect_identical(meta$length_out, 25L)
  expect_identical(meta$n_strata, 6L)
  expect_true(meta$n_strata_effective <= meta$n_strata)
  expect_true(meta$n_strata_realized <= meta$n_strata_effective)
  expect_identical(meta$confidence_level, 0.95)

  rng <- meta$bootstrap_rng
  expect_identical(rng$schema, "lvsr_cpp_bootstrap_rng_v1")
  expect_identical(rng$engine, "std::mt19937_64")
  expect_identical(rng$base_seed_hex, "0x00000000B5297A4D")
  expect_identical(rng$bootstrap_iteration_origin, 0L)
  expect_true(rng$thread_invariant)
  expect_false(rng$uses_R_rng)
  expect_match(rng$seed_derivation, "SplitMix64")

  expect_identical(meta$common_genes_fingerprint$count, 8L)
  expect_identical(meta$pearson_fingerprint$level, "grid")
  expect_match(meta$pearson_fingerprint$aligned_matrix_sha256, "^[0-9a-f]{64}$")
  expect_identical(
    meta$source_fingerprint$lee$upstream_input_fingerprint,
    lee_fp
  )
  expect_true(meta$preprocessing_provenance$R_rng$used)
  expect_identical(
    meta$preprocessing_provenance$downsample$mode,
    "random_fraction_without_replacement"
  )
  expect_match(
    meta$preprocessing_provenance$processed_inputs$pair_vectors_sha256,
    "^[0-9a-f]{64}$"
  )
  expect_identical(meta$provenance$bootstrap$rng, meta$bootstrap_rng)
  expect_true(all(vapply(curve$meta, identical, logical(1), meta)))
})

test_that("historical repaired OpenMP ABI remains registered", {
  registered <- names(getDLLRegisteredRoutines("geneSCOPE")$.Call)
  expect_true("_geneSCOPERebuild_native_openmp_info" %in% registered)
  expect_true("_geneSCOPERebuild_native_openmp_set_threads" %in% registered)
  info <- .Call("_geneSCOPERebuild_native_openmp_info", PACKAGE = "geneSCOPE")
  expect_type(info, "list")
})

test_that("Lee fingerprints separate observed data from permutation configuration", {
  obj <- make_lee_scope()
  layer <- obj@grid$grid1
  fp8 <- geneSCOPE:::.lee_input_fingerprint(
    layer$Xz, layer$W, layer$grid_info, "Xz", 8L, use_blocks = TRUE
  )
  fp2 <- geneSCOPE:::.lee_input_fingerprint(
    layer$Xz, layer$W, layer$grid_info, "Xz", 2L, use_blocks = TRUE
  )

  expect_identical(fp8$schema, "lee_input_fingerprint_v2")
  expect_identical(fp8$data, fp2$data)
  expect_false(identical(fp8$permutation, fp2$permutation))
  expect_true(geneSCOPE:::.lee_same_observed_data(fp8, fp2))
})

test_that("computeL stores complete canonical Lee's L/S2 input provenance", {
  obj <- make_lee_scope()
  native <- geneSCOPE:::.compute_lee_l(
    obj, grid_name = "grid1", norm_layer = "Xz", ncores = 1L,
    block_side = 2L, cache_inputs = FALSE
  )
  expect_null(native$input_cache)
  expect_identical(native$input_fingerprint$schema, "lee_input_fingerprint_v2")
  expect_identical(native$weight_style, "B")
  expect_equal(native$S2, sum(Matrix::rowSums(native$W)^2), tolerance = 0)

  out <- computeL(
    obj, grid_name = "grid1", perms = 0L, ncores = 1L,
    block_side = 2L, cache_inputs = FALSE, verbose = FALSE
  )
  meta <- out@stats$grid1$LeeStats_Xz$meta
  expect_identical(meta$formula_id, "Lee2009_S2_v1")
  expect_identical(meta$norm_layer, "Xz")
  expect_identical(meta$weight_style, "B")
  expect_equal(meta$S2, native$S2, tolerance = 0)
  expect_identical(meta$data_fingerprint, meta$input_fingerprint$data)
  expect_identical(meta$permutation_fingerprint, meta$input_fingerprint$permutation)
  expect_identical(meta$permutation_backend, "skipped")
})

test_that("custom Lee layers must be numeric, finite, and centred but need not have unit variance", {
  obj <- make_lee_scope()
  obj@grid$grid1$CustomZ <- sweep(
    obj@grid$grid1$Xz, 2L, c(2, 5, 0.25), `*`
  )
  expect_silent(geneSCOPE:::.validate_lee_norm_layer(
    obj@grid$grid1$CustomZ, "CustomZ"
  ))
  expect_silent(geneSCOPE:::.compute_lee_l(
    obj, grid_name = "grid1", norm_layer = "CustomZ", ncores = 1L,
    cache_inputs = FALSE
  ))

  shifted <- obj@grid$grid1$CustomZ
  shifted[, 1L] <- shifted[, 1L] + 1
  expect_error(
    geneSCOPE:::.validate_lee_norm_layer(shifted, "CustomZ"),
    "column-centred"
  )
  nonfinite <- obj@grid$grid1$CustomZ
  nonfinite[1L, 1L] <- Inf
  expect_error(
    geneSCOPE:::.validate_lee_norm_layer(nonfinite, "CustomZ"),
    "non-finite"
  )
  expect_error(
    geneSCOPE:::.validate_lee_norm_layer(matrix("x", 2L, 2L), "CustomZ"),
    "numeric matrix"
  )
})

test_that("the default Xz layer also fails closed when it is not centred", {
  shifted_xz <- make_lee_scope()@grid$grid1$Xz
  shifted_xz[, 1L] <- shifted_xz[, 1L] + 0.5
  expect_error(
    geneSCOPE:::.validate_lee_norm_layer(shifted_xz, "Xz"),
    "column-centred"
  )
})

test_that("computeL and getTopLvsR reject invalid scalar controls before native work", {
  obj <- make_lee_scope()
  expect_error(computeL(obj, grid_name = "grid1", ncores = 0L), "ncores.*positive integer")
  expect_error(computeL(obj, grid_name = "grid1", block_side = 0L), "block_side.*positive integer")
  expect_error(computeL(obj, grid_name = "grid1", block_size = 0L), "block_size.*positive integer")
  expect_error(computeL(obj, grid_name = "grid1", within = 1), "within.*logical")
  expect_error(computeL(obj, grid_name = "grid1", use_blocks = 1), "use_blocks.*logical")
  expect_error(computeL(obj, grid_name = "grid1", mem_limit_GB = Inf), "mem_limit_GB.*positive finite")

  expect_error(getTopLvsR(obj, grid_name = "grid1", ncores = 0L), "ncores.*positive integer")
  expect_error(getTopLvsR(obj, grid_name = "grid1", block_side = 0L), "block_side.*positive integer")
  expect_error(getTopLvsR(obj, grid_name = "grid1", do_perm = 1), "do_perm.*logical")
  expect_error(getTopLvsR(obj, grid_name = "grid1", mem_limit_GB = NaN), "mem_limit_GB.*positive finite")
})

test_that("computeL defaults to an all-grid joint permutation", {
  expect_identical(formals(computeL)$use_blocks, FALSE)
  expect_identical(formals(geneSCOPE:::.compute_l)$use_blocks, FALSE)
  expect_identical(formals(geneSCOPE:::.compute_lee_l)$use_blocks, FALSE)
  expect_identical(formals(geneSCOPE:::.lee_input_fingerprint)$use_blocks, FALSE)

  obj <- make_lee_scope()
  global_native <- geneSCOPE:::.compute_lee_l(
    obj, grid_name = "grid1", ncores = 1L, block_side = 2L,
    use_blocks = FALSE, cache_inputs = FALSE
  )
  block_native <- geneSCOPE:::.compute_lee_l(
    obj, grid_name = "grid1", ncores = 1L, block_side = 2L,
    use_blocks = TRUE, cache_inputs = FALSE
  )
  expect_identical(unique(global_native$block_id), 1L)
  expect_gt(length(unique(block_native$block_id)), 1L)

  set.seed(29)
  global <- computeL(
    obj, grid_name = "grid1", perms = 19L, block_size = 7L,
    block_side = 2L, use_blocks = FALSE, ncores = 1L,
    cache_inputs = FALSE, verbose = FALSE
  )
  set.seed(29)
  one_block <- computeL(
    obj, grid_name = "grid1", perms = 19L, block_size = 7L,
    block_side = 100L, use_blocks = TRUE, ncores = 1L,
    cache_inputs = FALSE, verbose = FALSE
  )
  set.seed(29)
  spatial_blocks <- computeL(
    obj, grid_name = "grid1", perms = 19L, block_size = 7L,
    block_side = 2L, use_blocks = TRUE, ncores = 1L,
    cache_inputs = FALSE, verbose = FALSE
  )

  global_stats <- global@stats$grid1$LeeStats_Xz
  one_block_stats <- one_block@stats$grid1$LeeStats_Xz
  expect_identical(global_stats$P, one_block_stats$P)
  expect_identical(global_stats$FDR, one_block_stats$FDR)
  expect_false(identical(global_stats$P, spatial_blocks@stats$grid1$LeeStats_Xz$P))
  expect_identical(global_stats$meta$permutation_scheme, "global_joint_shuffle")
  expect_identical(one_block_stats$meta$permutation_scheme, "spatial_block_joint")
  expect_match(global_stats$meta$permutation_backend, "lee_perm \\(all-grid")
  expect_match(one_block_stats$meta$permutation_backend, "lee_perm_block")

  cached_global <- computeL(
    obj, grid_name = "grid1", perms = 0L, block_side = 2L,
    use_blocks = FALSE, ncores = 1L, cache_inputs = TRUE, verbose = FALSE
  )
  cached_block <- computeL(
    cached_global, grid_name = "grid1", perms = 0L, block_side = 2L,
    use_blocks = TRUE, ncores = 1L, cache_inputs = TRUE, verbose = FALSE
  )
  global_cache <- cached_global@stats$grid1[["_lee_input_cache"]]
  block_cache <- cached_block@stats$grid1[["_lee_input_cache"]]
  expect_false(global_cache$use_blocks)
  expect_true(block_cache$use_blocks)
  expect_gt(length(unique(block_cache$block_id)), 1L)
})

test_that("getTopLvsR accepts a new block size but rejects stale observed inputs", {
  obj <- computeL(
    make_lee_scope(), grid_name = "grid1", perms = 0L, ncores = 1L,
    block_side = 8L, cache_inputs = FALSE, verbose = FALSE
  )
  obj@stats$grid1$.pearson_cor <- stats::cor(obj@grid$grid1$Xz)

  expect_silent(getTopLvsR(
    obj, grid_name = "grid1", pear_level = "grid", top_n = 1L,
    do_perm = TRUE, perms = 2L, block_side = 2L, ncores = 1L,
    CI_rule = "none", verbose = FALSE
  ))

  stale_formula <- obj
  stale_formula@stats$grid1$LeeStats_Xz$meta$formula_id <- NULL
  expect_error(
    getTopLvsR(stale_formula, grid_name = "grid1", pear_level = "grid",
      do_perm = FALSE, verbose = FALSE),
    "Rerun computeL"
  )

  stale_data <- obj
  stale_data@grid$grid1$Xz[1L, 1L] <- stale_data@grid$grid1$Xz[1L, 1L] + 0.25
  stale_data@grid$grid1$Xz[2L, 1L] <- stale_data@grid$grid1$Xz[2L, 1L] - 0.25
  expect_error(
    getTopLvsR(stale_data, grid_name = "grid1", pear_level = "grid",
      do_perm = FALSE, verbose = FALSE),
    "do not match.*Rerun computeL"
  )
})

test_that("getTopLvsR uses the recorded custom norm layer without Xz fallback", {
  obj <- make_lee_scope()
  obj@grid$grid1$CustomZ <- sweep(obj@grid$grid1$Xz, 2L, c(2, 3, 4), `*`)
  obj <- computeL(
    obj, grid_name = "grid1", norm_layer = "CustomZ", perms = 0L,
    ncores = 1L, cache_inputs = FALSE, verbose = FALSE
  )
  obj@stats$grid1$.pearson_cor <- stats::cor(obj@grid$grid1$CustomZ)

  expect_error(
    getTopLvsR(obj, grid_name = "grid1", pear_level = "grid",
      lee_stats_layer = "LeeStats_CustomZ", expr_layer = "Xz",
      do_perm = FALSE, verbose = FALSE),
    "does not match the norm_layer"
  )

  no_custom <- obj
  no_custom@grid$grid1$CustomZ <- NULL
  expect_error(
    getTopLvsR(no_custom, grid_name = "grid1", pear_level = "grid",
      lee_stats_layer = "LeeStats_CustomZ", do_perm = FALSE, verbose = FALSE),
    "fallback to Xz is prohibited"
  )
})

test_that("Lee clustering inputs require corrected formula provenance only for the L slot", {
  obj <- computeL(
    make_lee_scope(), grid_name = "grid1", perms = 0L, ncores = 1L,
    cache_inputs = FALSE, verbose = FALSE
  )
  expect_silent(geneSCOPE:::.extract_scope_layers(
    obj, grid_name = "grid1", stats_layer = "LeeStats_Xz",
    similarity_slot = "L", use_significance = FALSE, verbose = FALSE
  ))

  stale <- obj
  stale@stats$grid1$LeeStats_Xz$meta$formula_id <- NULL
  expect_error(
    geneSCOPE:::.extract_scope_layers(
      stale, grid_name = "grid1", stats_layer = "LeeStats_Xz",
      similarity_slot = "L", use_significance = FALSE, verbose = FALSE
    ),
    "Rerun computeL.*before clusterGenes"
  )

  custom <- stale
  custom@stats$grid1$LeeStats_Xz$custom_similarity <- custom@stats$grid1$LeeStats_Xz$L
  expect_silent(geneSCOPE:::.extract_scope_layers(
    custom, grid_name = "grid1", stats_layer = "LeeStats_Xz",
    similarity_slot = "custom_similarity", use_significance = FALSE, verbose = FALSE
  ))
})

test_that("L-vs-r curve construction rejects legacy Lee matrices", {
  obj <- computeL(
    make_lee_scope(), grid_name = "grid1", perms = 0L, ncores = 1L,
    cache_inputs = FALSE, verbose = FALSE
  )
  obj@stats$grid1$LeeStats_Xz$meta$formula_id <- NULL
  expect_error(
    geneSCOPE:::.compute_l_vs_r_curve(
      obj, grid_name = "grid1", level = "grid", B = 20L,
      ncores = 1L, verbose = FALSE
    ),
    "Rerun computeL.*before computeLvsRCurve"
  )
})

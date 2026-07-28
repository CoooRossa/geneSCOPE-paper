test_that("Stage-2 early return graph contains only filtered assigned edges", {
  genes <- c("A", "B", "C", "D")
  filtered <- matrix(0, 4, 4, dimnames = list(genes, genes))
  filtered["A", "B"] <- filtered["B", "A"] <- 0.8
  filtered["C", "D"] <- filtered["D", "C"] <- 0.6
  raw <- filtered
  raw["A", "D"] <- raw["D", "A"] <- 0.95
  filtered_sp <- Matrix::Matrix(filtered, sparse = TRUE)
  stage1 <- list(
    kept_genes = genes,
    genes_all = genes,
    stage1_membership = setNames(c(1L, 1L, 2L, 2L), genes),
    stage1_consensus = filtered_sp,
    g_stage1 = igraph::graph_from_adjacency_matrix(filtered_sp, mode = "undirected", weighted = TRUE),
    mode_selected = "safe",
    stage1_membership_matrix = matrix(c(1L, 1L, 2L, 2L), ncol = 1L, dimnames = list(genes, NULL)),
    stage1_backend = "igraph",
    stage1_algo = "louvain",
    stage1_algo_per_run = "louvain",
    similarity_sub = filtered_sp,
    A_sub = filtered_sp
  )
  config <- list(
    use_mh_weight = FALSE, CI95_filter = FALSE,
    stage2_algo = NULL, algo = "louvain", prefer_fast = FALSE,
    min_cluster_size = 1L, use_log1p_weight = TRUE
  )
  out <- geneSCOPE:::.stage2_refine_workflow_v2(
    stage1, raw, config, verbose = FALSE
  )
  ed <- igraph::as_data_frame(out$consensus_graph, what = "edges")
  keys <- sort(paste(pmin(ed$from, ed$to), pmax(ed$from, ed$to), sep = "|"))
  expect_identical(keys, c("A|B", "C|D"))
  expect_equal(sort(ed$weight), sort(log1p(c(0.8, 0.6))), tolerance = 1e-15)
  expect_setequal(igraph::V(out$consensus_graph)$name, genes)
})

test_that("non-early MH backbone can recover only Stage-1 filtered assigned edges", {
  genes <- c("A", "B", "C", "D")
  raw <- matrix(0, 4, 4, dimnames = list(genes, genes))
  raw["A", "B"] <- raw["B", "A"] <- 0.40
  raw["A", "C"] <- raw["C", "A"] <- 0.99
  raw["B", "C"] <- raw["C", "B"] <- 0.70
  raw["A", "D"] <- raw["D", "A"] <- 0.80
  stage1_filtered <- matrix(0, 4, 4, dimnames = list(genes, genes))
  stage1_filtered["A", "B"] <- stage1_filtered["B", "A"] <- 0.40
  stage1_filtered["A", "D"] <- stage1_filtered["D", "A"] <- 0.80
  membership <- setNames(c(1L, 1L, 1L, NA_integer_), genes)
  mh_zero <- matrix(0, 4, 4, dimnames = list(genes, genes))
  config <- list(
    min_cutoff = 0, use_significance = FALSE, significance_max = 0.05,
    pct_min = "q0", CI95_filter = FALSE, CI_rule = "remove_within",
    use_mh_weight = TRUE, keep_stage1_backbone = TRUE,
    backbone_floor_q = 0.02, use_log1p_weight = FALSE,
    post_smooth = FALSE, post_smooth_quant = c(0.05, 0.95),
    post_smooth_power = 1
  )
  spec <- list(kept_genes = genes)
  expect_error(
    geneSCOPE:::.stage2_native_inputs_materialize(
      spec, Matrix::Matrix(raw, sparse = TRUE), FDR = NULL,
      aux_stats = NULL, pearson_matrix = NULL, mh_object = mh_zero,
      stage1_membership_labels = membership, config = config, verbose = FALSE
    ),
    "requires stage1\\$similarity_sub/A_sub"
  )
  out <- geneSCOPE:::.stage2_native_inputs_materialize(
    spec, Matrix::Matrix(raw, sparse = TRUE), FDR = NULL,
    aux_stats = NULL, pearson_matrix = NULL, mh_object = mh_zero,
    stage1_membership_labels = membership, config = config, verbose = FALSE,
    stage1_similarity = Matrix::Matrix(stage1_filtered, sparse = TRUE)
  )
  keys <- paste(
    pmin(out$edges_corr$from, out$edges_corr$to),
    pmax(out$edges_corr$from, out$edges_corr$to),
    sep = "|"
  )
  expect_identical(keys, "A|B")
  expect_lte(out$edges_corr$L_corr, stage1_filtered["A", "B"])
  expect_false(any(grepl("D", keys, fixed = TRUE)))
})

test_that("weight metadata validation rejects ambiguity", {
  layer <- list(
    grid_info = data.frame(grid_id = c("s1", "s2"), gx = c(1, 1), gy = c(1, 1)),
    xbins_eff = 2L, ybins_eff = 2L
  )
  expect_error(
    geneSCOPE:::.spw_validate_grid_metadata(layer, "grid1"),
    "duplicated"
  )
  expect_error(
    computeWeights(NULL, style = "C"),
    "Only B"
  )
  expect_error(
    computeWeights(NULL, style = "B", topology = "fuzzy_queen"),
    "fuzzy_queen/fuzzy_hex weights are disabled"
  )
})

test_that("pairwise BH adjusts only unique unordered off-diagonal pairs", {
  P <- matrix(c(
    NA, 0.01, 0.02,
    0.01, NA, 0.20,
    0.02, 0.20, NA
  ), 3, 3, byrow = TRUE)
  got <- geneSCOPE:::.lee_pairwise_p_adjust(P, "BH")
  expected <- p.adjust(c(0.01, 0.02, 0.20), "BH")
  expect_equal(got[upper.tri(got)], expected, tolerance = 0)
  expect_equal(got, t(got), tolerance = 0)
  expect_true(all(is.na(diag(got))))
})

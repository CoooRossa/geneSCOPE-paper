graph_test_keys <- function(graph) {
  edges <- igraph::as_data_frame(graph, what = "edges")
  sort(paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to), sep = "|"))
}

test_that("legacy L backbone and membership are unchanged while graph roles are explicit", {
  genes <- c("A", "B", "E", "C", "D")
  filtered <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  filtered["A", "E"] <- filtered["E", "A"] <- 0.8
  filtered["B", "E"] <- filtered["E", "B"] <- 0.7
  filtered["C", "D"] <- filtered["D", "C"] <- 0.6
  filtered["A", "C"] <- filtered["C", "A"] <- 0.5
  filtered_sp <- Matrix::Matrix(filtered, sparse = TRUE)

  restart_memberships <- cbind(
    c(1L, 1L, 1L, 2L, 2L),
    c(1L, 1L, 1L, 2L, 2L),
    c(1L, 1L, 1L, 2L, 3L),
    c(1L, 1L, 1L, 1L, 1L)
  )
  rownames(restart_memberships) <- genes
  consensus <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  consensus[cbind(c(1, 1, 2, 4), c(2, 3, 3, 5))] <- c(1, 1, 1, 0.75)
  consensus <- consensus + t(consensus)
  consensus_sp <- Matrix::Matrix(consensus, sparse = TRUE)

  membership_before <- setNames(c(1L, 1L, 1L, 2L, 2L), genes)
  stage1 <- list(
    kept_genes = genes,
    genes_all = genes,
    stage1_membership = membership_before,
    stage1_consensus = consensus_sp,
    g_stage1 = igraph::graph_from_adjacency_matrix(filtered_sp,
      mode = "undirected", weighted = TRUE, diag = FALSE
    ),
    mode_selected = "safe",
    stage1_membership_matrix = restart_memberships,
    stage1_backend = "igraph",
    stage1_algo = "louvain",
    stage1_algo_per_run = "louvain",
    similarity_sub = filtered_sp,
    A_sub = filtered_sp
  )
  config <- list(
    use_mh_weight = FALSE,
    CI95_filter = FALSE,
    stage2_algo = NULL,
    algo = "louvain",
    prefer_fast = FALSE,
    min_cluster_size = 1L,
    use_log1p_weight = TRUE,
    consensus_thr = 0.75
  )

  out <- geneSCOPE:::.stage2_refine_workflow_v2(
    stage1, filtered, config, verbose = FALSE
  )

  expect_identical(out$membership, membership_before)
  expect_identical(
    graph_test_keys(out$consensus_graph),
    sort(c("A|C", "A|E", "B|E", "C|D"))
  )

  backbone_edges <- igraph::as_data_frame(out$consensus_graph, what = "edges")
  backbone_key <- paste(
    pmin(backbone_edges$from, backbone_edges$to),
    pmax(backbone_edges$from, backbone_edges$to), sep = "|"
  )
  expected_raw <- c("A|C" = 0.5, "A|E" = 0.8, "B|E" = 0.7, "C|D" = 0.6)
  expected_frequency <- c("A|C" = 0.25, "A|E" = 1, "B|E" = 1, "C|D" = 0.75)
  expect_equal(backbone_edges$weight, unname(log1p(expected_raw[backbone_key])), tolerance = 0)
  expect_equal(backbone_edges$raw_L, unname(expected_raw[backbone_key]), tolerance = 0)
  expect_identical(backbone_edges$weight_transform, rep("log1p", nrow(backbone_edges)))
  expect_equal(
    backbone_edges$consensus_frequency,
    unname(expected_frequency[backbone_key]), tolerance = 0
  )
  expect_identical(
    backbone_edges$is_consensus_edge,
    unname(expected_frequency[backbone_key] >= 0.75)
  )

  expect_identical(
    graph_test_keys(out$consensus_frequency_graph),
    sort(c("A|B", "A|E", "B|E", "C|D"))
  )
  expect_identical(
    graph_test_keys(out$paper_consensus_graph),
    sort(c("A|E", "B|E", "C|D"))
  )
  expect_identical(
    igraph::graph_attr(out$consensus_graph, "graph_role"),
    "thresholded_L_backbone"
  )
  expect_identical(
    igraph::graph_attr(out$consensus_frequency_graph, "graph_role"),
    "consensus_frequency_graph"
  )
  expect_identical(
    igraph::graph_attr(out$paper_consensus_graph, "graph_role"),
    "paper_consensus_network"
  )
  expect_equal(
    igraph::graph_attr(out$paper_consensus_graph, "consensus_threshold"),
    0.75, tolerance = 0
  )

  layers <- geneSCOPE:::.store_cluster_graph_views(
    list(untouched = "sentinel"), "g_consensus", out
  )
  expect_identical(layers$untouched, "sentinel")
  expect_identical(layers$g_consensus, out$consensus_graph)
  expect_identical(
    layers$g_consensus__consensus_frequency,
    out$consensus_frequency_graph
  )
  expect_identical(
    layers$g_consensus__paper_consensus,
    out$paper_consensus_graph
  )
})

test_that("identity backbone weights receive matching transform metadata", {
  vertices <- c("A", "B")
  graph <- igraph::graph_from_data_frame(
    data.frame(from = "A", to = "B", weight = 0.4, sign = "pos"),
    directed = FALSE, vertices = vertices
  )
  raw <- matrix(c(0, 0.4, 0.4, 0), 2, 2,
    dimnames = list(vertices, vertices)
  )
  membership <- matrix(c(1L, 1L), ncol = 1L,
    dimnames = list(vertices, "run1")
  )
  consensus <- Matrix::Matrix(raw > 0, sparse = TRUE)
  dimnames(consensus) <- list(vertices, vertices)

  views <- geneSCOPE:::.build_consensus_graph_views(
    graph, consensus, membership, raw,
    consensus_thr = 0.95, use_log1p_weight = FALSE
  )
  edges <- igraph::as_data_frame(views$backbone_graph, what = "edges")
  expect_equal(edges$weight, 0.4, tolerance = 0)
  expect_equal(edges$raw_L, 0.4, tolerance = 0)
  expect_identical(edges$weight_transform, "identity")
  expect_identical(
    igraph::graph_attr(views$backbone_graph, "weight_semantics"),
    "raw_L"
  )
})

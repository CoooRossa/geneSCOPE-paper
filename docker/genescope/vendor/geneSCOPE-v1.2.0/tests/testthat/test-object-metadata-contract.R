make_legacy_object_metadata_scope <- function(
    platform = "10x Xenium",
    dataset = "P5") {
  obj <- make_lee_scope()
  genes <- rownames(obj@meta.data)

  true_gene_metadata <- data.frame(
    feature_type = rep("Gene Expression", length(genes)),
    final_cluster = setNames(c(1L, 1L, 2L), NULL),
    stringsAsFactors = FALSE,
    row.names = genes
  )
  obj@meta.data <- cbind(
    true_gene_metadata,
    platform = rep(NA_character_, length(genes)),
    dataset = rep(NA_character_, length(genes))
  )
  obj@meta.data <- rbind(
    obj@meta.data,
    data.frame(
      feature_type = NA_character_,
      final_cluster = NA_integer_,
      platform = platform,
      dataset = dataset,
      stringsAsFactors = FALSE,
      row.names = "__scope_platform__"
    )
  )

  L <- matrix(
    c(1, 0.8, 0.2, 0.8, 1, 0.4, 0.2, 0.4, 1),
    nrow = 3L,
    dimnames = list(genes, genes)
  )
  FDR <- matrix(
    c(NA, 0.01, 0.20, 0.01, NA, 0.03, 0.20, 0.03, NA),
    nrow = 3L,
    dimnames = list(genes, genes)
  )
  membership <- setNames(c(1L, 1L, 2L), genes)
  graph <- igraph::graph_from_adjacency_matrix(
    L * (L >= 0.4) * (1 - diag(3L)),
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  obj@stats$grid1 <- list(
    LeeStats_Xz = list(L = L, FDR = FDR),
    cluster_membership = membership,
    g_consensus = graph
  )

  list(
    object = obj,
    genes = genes,
    true_gene_metadata = true_gene_metadata
  )
}

test_that("new object metadata is kept out of gene metadata", {
  before <- make_lee_scope()
  state <- list(
    objects = list(scope_obj = before),
    params = list(data_type = "xenium", verbose = FALSE)
  )

  finalized <- geneSCOPE:::.finalize_output_object(state)$objects$scope_obj
  finalized <- geneSCOPE:::.set_scope_object_metadata(
    finalized,
    dataset = "P5"
  )

  expect_identical(rownames(finalized@meta.data), rownames(before@meta.data))
  expect_identical(finalized@meta.data, before@meta.data)
  expect_false(any(c("__scope_platform__", "__scope_metadata__") %in%
    rownames(finalized@meta.data)))
  expect_false(any(c("platform", "dataset") %in% names(finalized@meta.data)))
  expect_identical(finalized@stats$object_metadata$schema_version, "1.0")
  expect_identical(finalized@stats$object_metadata$platform, "Xenium")
  expect_identical(finalized@stats$object_metadata$dataset, "P5")

  # The scalar mirrors remain readable by v1.0.0/v1.0.1-era code without
  # introducing object-level values into the gene table.
  expect_identical(finalized@stats$platform, "Xenium")
  expect_identical(finalized@stats$dataset, "P5")
})

test_that("object metadata reader follows canonical, legacy stats, legacy gene-table order", {
  obj <- make_lee_scope()
  obj@meta.data$platform <- c("CosMx", "CosMx", NA_character_)
  obj@meta.data$dataset <- c("meta-dataset", "meta-dataset", NA_character_)
  obj@stats$platform <- "Visium"
  obj@stats$dataset <- "stats-dataset"
  obj@stats$object_metadata <- list(
    schema_version = "1.0",
    platform = "Xenium",
    dataset = "canonical-dataset"
  )

  expect_identical(
    geneSCOPE:::.get_scope_object_metadata(obj, "platform"),
    "Xenium"
  )
  expect_identical(
    geneSCOPE:::.get_scope_object_metadata(obj, "dataset"),
    "canonical-dataset"
  )

  obj@stats$object_metadata <- NULL
  expect_identical(
    geneSCOPE:::.get_scope_object_metadata(obj, "platform"),
    "Visium"
  )
  expect_identical(
    geneSCOPE:::.get_scope_object_metadata(obj, "dataset"),
    "stats-dataset"
  )

  obj@stats$platform <- NULL
  obj@stats$dataset <- NULL
  expect_identical(
    geneSCOPE:::.get_scope_object_metadata(obj, "platform"),
    "CosMx"
  )
  expect_identical(
    geneSCOPE:::.get_scope_object_metadata(obj, "dataset"),
    "meta-dataset"
  )
  expect_identical(
    geneSCOPE:::.get_scope_object_metadata(obj, "missing", default = "fallback"),
    "fallback"
  )
})

test_that("explicit legacy metadata migration is idempotent and analysis-invariant", {
  fixture <- make_legacy_object_metadata_scope()
  legacy <- fixture$object

  protected_before <- list(
    coord = legacy@coord,
    grid = legacy@grid,
    cells = legacy@cells,
    density = legacy@density,
    analysis = legacy@stats$grid1
  )

  migrated <- migrateScopeMetadata(legacy)
  migrated_twice <- migrateScopeMetadata(migrated)

  expect_identical(migrated_twice, migrated)
  expect_identical(rownames(migrated@meta.data), fixture$genes)
  expect_identical(
    migrated@meta.data[, names(fixture$true_gene_metadata), drop = FALSE],
    fixture$true_gene_metadata
  )
  expect_false("__scope_platform__" %in% rownames(migrated@meta.data))
  expect_false(any(c("platform", "dataset") %in% names(migrated@meta.data)))
  expect_identical(migrated@stats$object_metadata$platform, "Xenium")
  expect_identical(migrated@stats$object_metadata$dataset, "P5")

  expect_identical(migrated@coord, protected_before$coord)
  expect_identical(migrated@grid, protected_before$grid)
  expect_identical(migrated@cells, protected_before$cells)
  expect_identical(migrated@density, protected_before$density)
  expect_identical(migrated@stats$grid1, protected_before$analysis)

  # Name the scientific invariants explicitly so a later metadata refactor
  # cannot silently perturb the Lee kernel, W, clustering, or graph objects.
  expect_identical(
    migrated@stats$grid1$LeeStats_Xz$L,
    legacy@stats$grid1$LeeStats_Xz$L
  )
  expect_identical(
    migrated@stats$grid1$LeeStats_Xz$FDR,
    legacy@stats$grid1$LeeStats_Xz$FDR
  )
  expect_identical(migrated@grid$grid1$W, legacy@grid$grid1$W)
  expect_identical(
    migrated@stats$grid1$cluster_membership,
    legacy@stats$grid1$cluster_membership
  )
  expect_identical(
    migrated@stats$grid1$g_consensus,
    legacy@stats$grid1$g_consensus
  )
})

test_that("ambiguous legacy geneSCOPE labels require an explicit platform override", {
  fixture <- make_legacy_object_metadata_scope(
    platform = "geneSCOPE",
    dataset = "legacy-dataset"
  )
  legacy <- fixture$object
  legacy@stats$platform <- "geneSCOPE"
  legacy@stats$dataset <- "legacy-dataset"

  expect_warning(
    unresolved <- migrateScopeMetadata(legacy),
    "geneSCOPE"
  )
  expect_identical(
    unresolved@stats$object_metadata$platform,
    "geneSCOPE"
  )
  expect_false(unresolved@stats$object_metadata$platform %in%
    c("Xenium", "CosMx", "Visium"))

  corrected <- expect_silent(migrateScopeMetadata(
    legacy,
    platform = "xenium",
    dataset = "P5"
  ))
  corrected_twice <- expect_silent(migrateScopeMetadata(
    corrected,
    platform = "Xenium",
    dataset = "P5"
  ))

  expect_identical(corrected_twice, corrected)
  expect_identical(corrected@stats$object_metadata$platform, "Xenium")
  expect_identical(corrected@stats$object_metadata$dataset, "P5")
  expect_identical(corrected@stats$platform, "Xenium")
  expect_identical(corrected@stats$dataset, "P5")
  expect_identical(rownames(corrected@meta.data), fixture$genes)
  expect_false(any(c("platform", "dataset") %in% names(corrected@meta.data)))

  expect_identical(corrected@grid, legacy@grid)
  expect_identical(
    corrected@stats$grid1$LeeStats_Xz,
    legacy@stats$grid1$LeeStats_Xz
  )
  expect_identical(
    corrected@stats$grid1$cluster_membership,
    legacy@stats$grid1$cluster_membership
  )
  expect_identical(
    corrected@stats$grid1$g_consensus,
    legacy@stats$grid1$g_consensus
  )
})

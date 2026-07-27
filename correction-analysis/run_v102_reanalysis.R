#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: run_v102_reanalysis.R SAMPLE_ID XENIUM_OUTS RESULT_ROOT [NCORES]")
}

script_dir <- local({
  x <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", x[[1L]])))
})
source(file.path(script_dir, "..", "main-text-scripts", "freeze_helpers.R"))

suppressPackageStartupMessages({
  library(geneSCOPE)
  library(jsonlite)
})
configure_freeze_runtime()

sample_id <- toupper(args[[1L]])
input_dir <- normalizePath(args[[2L]], mustWork = TRUE)
result_root <- normalizePath(args[[3L]], mustWork = FALSE)
ncores <- if (length(args) >= 4L) as.integer(args[[4L]]) else 8L
if (!is.finite(ncores) || ncores < 1L) stop("NCORES must be a positive integer.")

config <- utils::read.delim(file.path(script_dir, "samples.tsv"), stringsAsFactors = FALSE)
cfg <- config[config$sample_id == sample_id, , drop = FALSE]
if (nrow(cfg) != 1L) stop("Unknown SAMPLE_ID: ", sample_id)
cfg <- cfg[1L, ]
roi_file <- normalizePath(file.path(script_dir, cfg$roi_file), mustWork = TRUE)
out_dir <- file.path(result_root, sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!identical(as.character(utils::packageVersion("geneSCOPE")), "1.0.2")) {
  stop("The correction workflow requires geneSCOPE 1.0.2.")
}
if (!identical(formals(geneSCOPE::computeL)$use_blocks, FALSE) ||
    !identical(formals(geneSCOPE::getTopLvsR)$use_blocks, FALSE) ||
    !identical(eval(formals(geneSCOPE::getTopLvsR)$p_adj_mode)[[1L]], "BH")) {
  stop("Installed geneSCOPE API defaults do not match the frozen v1.0.2 contract.")
}
freeze_source <- require_freeze_source_metadata()
paper_commit <- require_paper_commit_metadata()
gate_max_abs_L_diff <- canonical_lee_s2_gate()

seed <- 1L
perms <- 1000L
grid_name <- paste0("grid", cfg$grid_um)
cluster_col <- paste0(
  "correction_", cfg$cluster_pct, "_res", cfg$resolution, "_",
  grid_name, "_freq", cfg$consensus_threshold
)

scope_obj <- createSCOPE(
  data_dir = input_dir,
  grid_length = as.numeric(cfg$grid_um),
  seg_type = "cell",
  coord_file = roi_file,
  ncores = ncores
)
scope_obj <- addSingleCells(scope_obj = scope_obj, xenium_dir = input_dir)
scope_obj <- normalizeSingleCells(
  scope_obj = scope_obj,
  input_layer = "counts",
  output_layer = "logCPM",
  scale_factor = 1e4
)
scope_obj <- normalizeMoleculesInGrid(scope_obj = scope_obj, grid_name = grid_name)
scope_obj <- computeWeights(
  scope_obj = scope_obj,
  grid_name = grid_name,
  style = "B",
  topology = "auto",
  store_mat = TRUE,
  store_listw = TRUE,
  ncores = ncores
)

reset_freeze_rng(seed)
scope_obj <- computeL(
  scope_obj = scope_obj,
  grid_name = grid_name,
  ncores = ncores,
  perms = perms,
  use_blocks = FALSE,
  norm_layer = "Xz",
  use_bigmemory = FALSE
)
lee <- scope_obj@stats[[grid_name]][["LeeStats_Xz"]]
if (is.null(lee$meta) ||
    !lee$meta$formula_id %in% c("Lee2009_S2_v1", "Lee_S2_v1") ||
    !identical(lee$meta$use_blocks, FALSE) ||
    !identical(lee$meta$permutation_scheme, "global_joint_shuffle") ||
    !identical(as.integer(lee$meta$perms), perms)) {
  stop("Observed Lee layer failed the formula/permutation provenance gate.")
}

scope_obj <- computeCorrelation(
  scope_obj = scope_obj,
  level = "cell",
  layer = "logCPM",
  method = "pearson",
  blocksize = 2000,
  ncores = ncores
)
curve_name <- paste0("LR_curve_", cfg$grid_um)
reset_freeze_rng(seed)
scope_obj <- computeLvsRCurve(
  scope_obj = scope_obj,
  level = "cell",
  grid_name = grid_name,
  B = 1000,
  ncores = ncores,
  downsample = 0.05,
  k_max = 2000,
  n_strata = 1000,
  min_rel_width = 0.15,
  widen_span = 0.1,
  curve_name = curve_name
)

reset_freeze_rng(seed)
top_all <- getTopLvsR(
  scope_obj = scope_obj,
  grid_name = grid_name,
  pear_level = "cell",
  L_range = c(0, 1),
  top_n = as.integer(cfg$top_n),
  direction = "largest",
  do_perm = TRUE,
  perms = perms,
  use_blocks = FALSE,
  ncores = ncores,
  p_adj_mode = "BH_universe",
  pval_mode = "uniform",
  curve_layer = curve_name,
  CI_rule = "remove_within"
)
delta_provenance <- attr(top_all, "permutation_provenance")
if (is.null(delta_provenance) ||
    !identical(delta_provenance$p_adj_mode, "BH_universe") ||
    !identical(delta_provenance$use_blocks, FALSE) ||
    delta_provenance$selected_pairs != delta_provenance$total_universe) {
  stop(
    "Delta inference did not cover the complete eligible universe (eligible=",
    delta_provenance$total_universe, ", selected=", delta_provenance$selected_pairs, ")."
  )
}
top_display <- filter_display_pairs(top_all)
top6_gate <- assert_reference_top6(script_dir, sample_id, top_display)

reset_freeze_rng(seed)
scope_obj <- clusterGenes(
  scope_obj = scope_obj,
  grid_name = grid_name,
  L_min = 0,
  algo = "leiden",
  resolution = as.numeric(cfg$resolution),
  pct_min = cfg$cluster_pct,
  cluster_name = cluster_col,
  graph_slot_name = cluster_col,
  use_log1p_weight = TRUE,
  use_consensus = TRUE,
  consensus_thr = as.numeric(cfg$consensus_threshold),
  n_restart = as.integer(cfg$n_restart),
  ncores = ncores
)

membership <- scope_obj@meta.data[[cluster_col]]
membership_chr <- as.character(membership)
membership_gate <- assert_reference_membership(
  script_dir, sample_id, rownames(scope_obj@meta.data), membership_chr
)
display_mapping <- read_display_mapping(script_dir, sample_id)
display_membership <- apply_display_mapping(membership_chr, display_mapping)
assigned <- !is.na(membership_chr) & nzchar(membership_chr) & membership_chr != "-1"
n_modules <- length(unique(membership_chr[assigned]))
n_assigned <- sum(assigned)
if (n_modules != cfg$expected_modules || n_assigned != cfg$expected_assigned) {
  stop(
    sample_id, " anchor mismatch: observed ", n_modules, "/", n_assigned,
    "; expected ", cfg$expected_modules, "/", cfg$expected_assigned
  )
}

dendro_audit <- NULL
if (identical(sample_id, "P5")) {
  reset_freeze_rng(seed)
  dnet_obj <- plotDendroNetwork(
    scope_obj = scope_obj,
    lee_stats_layer = "LeeStats_Xz",
    grid_name = grid_name,
    use_consensus_graph = TRUE,
    graph_slot_name = cluster_col,
    cluster_vec = cluster_col,
    IDelta_col_name = NULL,
    node_size = 4,
    edge_width = 3,
    label_size = 4,
    seed = seed,
    max.overlaps = 10,
    title = " ",
    tree_mode = "radial"
  )
  dendro_graph <- dnet_obj$graph
  if (is.null(dendro_graph) || !inherits(dendro_graph, "igraph")) {
    stop("P5 dendrogram audit did not return an igraph tree.")
  }
  gene_module <- stats::setNames(membership_chr, rownames(scope_obj@meta.data))
  endpoints <- igraph::as_data_frame(dendro_graph, what = "edges")
  endpoints$module_from <- unname(gene_module[endpoints$from])
  endpoints$module_to <- unname(gene_module[endpoints$to])
  cross <- endpoints[
    !is.na(endpoints$module_from) & !is.na(endpoints$module_to) &
      endpoints$module_from != endpoints$module_to,
    c("module_from", "module_to"), drop = FALSE
  ]
  cross$key <- apply(cross, 1L, function(z) paste(sort(z), collapse = "--"))
  cross <- cross[!duplicated(cross$key), c("module_from", "module_to"), drop = FALSE]
  module_graph <- igraph::graph_from_data_frame(cross, directed = FALSE)
  start_module <- unname(gene_module[["C3"]])
  end_module <- unname(gene_module[["GPX2"]])
  module_path <- names(igraph::shortest_paths(
    module_graph, from = start_module, to = end_module
  )$vpath[[1L]])

  path_result <- getDendroWalkPaths(
    dnet_obj, gene = "C3", gene2 = "GPX2", cutoff = 20L,
    max_paths = 50000L, verbose = FALSE
  )
  endpoint_paths <- Filter(function(p) {
    ends <- c(p[[1L]], p[[length(p)]])
    identical(ends, c("C3", "GPX2")) || identical(ends, c("GPX2", "C3"))
  }, path_result$paths)
  if (!length(endpoint_paths)) stop("No C3--GPX2 dendrogram path was returned.")
  path_lengths <- vapply(endpoint_paths, length, integer(1))
  gene_path <- endpoint_paths[[which.min(path_lengths)]]
  if (!identical(gene_path[[1L]], "C3")) gene_path <- rev(gene_path)

  display_lookup <- stats::setNames(
    as.character(display_mapping$display_module),
    as.character(display_mapping$current_module)
  )
  display_path <- unname(display_lookup[module_path])
  raw_path_string <- paste(module_path, collapse = "->")
  display_path_string <- paste(display_path, collapse = "->")
  if (!identical(raw_path_string, "1->2->10->6") ||
      !identical(display_path_string, "1->3->2->5")) {
    stop("P5 dendrogram path anchor changed: raw=", raw_path_string,
         ", display=", display_path_string)
  }
  dendro_audit <- data.frame(
    gene_query = "C3--GPX2",
    gene_path = paste(gene_path, collapse = "->"),
    raw_module_path = raw_path_string,
    display_module_path = display_path_string,
    stem_positioned_between = TRUE,
    direct_stem_to_broad_adjacency = FALSE,
    stringsAsFactors = FALSE
  )
  utils::write.table(dendro_audit, file.path(out_dir, "P5_dendro_path_audit.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
}

membership_out <- data.frame(
  gene = rownames(scope_obj@meta.data),
  raw_module_id = membership_chr,
  display_module_id = as.character(display_membership),
  stringsAsFactors = FALSE
)
membership_path <- file.path(out_dir, paste0(sample_id, "_module_membership.tsv"))
top_all_path <- file.path(out_dir, paste0(sample_id, "_top_pairs_all.tsv"))
top_display_path <- file.path(out_dir, paste0(sample_id, "_top_pairs_display_filter.tsv"))
top6_path <- file.path(out_dir, paste0(sample_id, "_Top6.tsv"))
scope_path <- file.path(out_dir, paste0(sample_id, "_scope_v102.rds"))
utils::write.table(membership_out, membership_path,
                   sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(as.data.frame(top_all), top_all_path,
                   sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(top_display, top_display_path,
                   sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(utils::head(top_display, 6L), top6_path,
                   sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(scope_obj, scope_path)

raw_candidates <- c(
  roi = roi_file,
  cell_feature_matrix = file.path(input_dir, "cell_feature_matrix.h5"),
  cells_parquet = file.path(input_dir, "cells.parquet"),
  transcripts_parquet = file.path(input_dir, "transcripts.parquet")
)
raw_candidates <- raw_candidates[file.exists(raw_candidates)]
input_sha256 <- as.list(vapply(raw_candidates, sha256_file, character(1)))
output_sha256 <- as.list(vapply(
  c(membership = membership_path, top_all = top_all_path,
    top_display = top_display_path, top6 = top6_path, scope = scope_path),
  sha256_file, character(1)
))

manifest <- list(
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  sample_id = sample_id,
  package_version = as.character(utils::packageVersion("geneSCOPE")),
  gate_max_abs_L_diff = gate_max_abs_L_diff,
  formula_id = lee$meta$formula_id,
  formula_gate = "canonical Lee S2 provenance accepted",
  permutation = list(
    scheme = "global_joint_shuffle",
    permutations = perms,
    seed = seed,
    rng = "L'Ecuyer-CMRG",
    use_blocks = FALSE
  ),
  delta_permutation = delta_provenance,
  delta_adjustment = paste0(
    "eligible-universe-scaled BH for selected Top-N; frozen gate requires selected_pairs == total_universe"
  ),
  display_filter = list(q_Delta_max_exclusive = 0.05, L_min_exclusive = 0,
                        r_max_exclusive = 0.05, pct1_min_exclusive = 20,
                        pct2_min_exclusive = 20, ranking = "full-precision L-r"),
  clustering = list(pct_min = cfg$cluster_pct, resolution = cfg$resolution,
                    consensus_threshold = cfg$consensus_threshold,
                    n_restart = cfg$n_restart, modules = n_modules, assigned = n_assigned),
  weight_provenance = scope_obj@grid[[grid_name]]$weight_provenance,
  display_mapping_sha256 = sha256_file(file.path(
    script_dir, "display-mappings", paste0(sample_id, "_display_mapping.tsv")
  )),
  membership_freeze_gate = membership_gate,
  top6_freeze_gate = as.list(top6_gate),
  package_source_commit = freeze_source$source_commit,
  package_vendor_tree_sha256 = freeze_source$vendor_tree_sha256,
  paper_workflow_commit = paper_commit,
  p5_dendro_path = if (!is.null(dendro_audit)) as.list(dendro_audit[1L, ]) else NULL,
  inputs = list(xenium_outs = input_dir, roi_file = roi_file, sha256 = input_sha256),
  outputs = list(top_pairs = nrow(top_all), display_pairs = nrow(top_display),
                 sha256 = output_sha256)
)
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 17)
message(sample_id, " completed: modules=", n_modules, ", assigned=", n_assigned,
        ", display_pairs=", nrow(top_display))

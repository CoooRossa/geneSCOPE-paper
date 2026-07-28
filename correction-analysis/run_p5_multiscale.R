#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: run_p5_multiscale.R P5_XENIUM_OUTS RESULT_ROOT [NCORES]")
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

input_dir <- normalizePath(args[[1L]], mustWork = TRUE)
result_root <- normalizePath(args[[2L]], mustWork = FALSE)
ncores <- if (length(args) >= 3L) as.integer(args[[3L]]) else 8L
if (!is.finite(ncores) || ncores < 1L) stop("NCORES must be a positive integer.")
if (!identical(as.character(utils::packageVersion("geneSCOPE")), "1.2.0")) {
  stop("The P5 multiscale workflow requires geneSCOPE 1.2.0.")
}
if (!identical(formals(geneSCOPE::computeL)$use_blocks, FALSE)) {
  stop("computeL() global-shuffle default is not frozen.")
}
freeze_source <- require_freeze_source_metadata()
paper_commit <- require_paper_commit_metadata()
gate_max_abs_L_diff <- canonical_lee_s2_gate()

roi_file <- normalizePath(file.path(script_dir, "..", "ROI-coordinate-files", "P5_roi.csv"),
                          mustWork = TRUE)
anchors <- utils::read.delim(file.path(script_dir, "reference_p5_multiscale.tsv"),
                             stringsAsFactors = FALSE)
seed <- 1L
perms <- 1000L
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(anchors))) {
  grid_um <- as.integer(anchors$grid_um[[i]])
  grid_name <- paste0("grid", grid_um)
  out_dir <- file.path(result_root, grid_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  scope_obj <- createSCOPE(
    data_dir = input_dir,
    grid_length = grid_um,
    seg_type = "cell",
    coord_file = roi_file,
    ncores = ncores
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
    norm_layer = "Xz",
    perms = perms,
    use_blocks = FALSE,
    use_bigmemory = FALSE,
    ncores = ncores
  )
  lee_meta <- scope_obj@stats[[grid_name]][["LeeStats_Xz"]]$meta
  if (is.null(lee_meta) ||
      !lee_meta$formula_id %in% c("Lee2009_S2_v1", "Lee_S2_v1") ||
      !identical(lee_meta$use_blocks, FALSE) ||
      !identical(lee_meta$permutation_scheme, "global_joint_shuffle") ||
      !identical(as.integer(lee_meta$perms), perms)) {
    stop(grid_name, " failed Lee formula/permutation provenance.")
  }

  cluster_col <- paste0("correction_q95_res0.1_", grid_name, "_freq0.95")
  reset_freeze_rng(seed)
  scope_obj <- clusterGenes(
    scope_obj = scope_obj,
    grid_name = grid_name,
    L_min = 0,
    algo = "leiden",
    resolution = 0.1,
    pct_min = "q95",
    cluster_name = cluster_col,
    graph_slot_name = cluster_col,
    use_log1p_weight = TRUE,
    use_consensus = TRUE,
    consensus_thr = 0.95,
    n_restart = 1000,
    ncores = ncores
  )

  membership <- as.character(scope_obj@meta.data[[cluster_col]])
  assigned_mask <- !is.na(membership) & nzchar(membership) & membership != "-1"
  n_modules <- length(unique(membership[assigned_mask]))
  n_assigned <- sum(assigned_mask)
  expected_modules <- as.integer(anchors$modules[[i]])
  expected_assigned <- as.integer(anchors$assigned[[i]])
  if (n_modules != expected_modules || n_assigned != expected_assigned) {
    stop(grid_name, " anchor mismatch: observed ", n_modules, "/", n_assigned,
         "; expected ", expected_modules, "/", expected_assigned)
  }

  membership_out <- data.frame(
    gene = rownames(scope_obj@meta.data), module_id = membership,
    stringsAsFactors = FALSE
  )
  membership_path <- file.path(out_dir, paste0(grid_name, "_membership.tsv"))
  scope_path <- file.path(out_dir, paste0(grid_name, "_scope_v102.rds"))
  utils::write.table(
    membership_out, membership_path,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  saveRDS(scope_obj, scope_path)
  raw_candidates <- c(
    roi = roi_file,
    cell_feature_matrix = file.path(input_dir, "cell_feature_matrix.h5"),
    cells_parquet = file.path(input_dir, "cells.parquet"),
    transcripts_parquet = file.path(input_dir, "transcripts.parquet")
  )
  raw_candidates <- raw_candidates[file.exists(raw_candidates)]
  jsonlite::write_json(list(
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    sample_id = "P5",
    grid_um = grid_um,
    artifact_series = "v102",
    package_version = as.character(utils::packageVersion("geneSCOPE")),
    package_source_commit = freeze_source$source_commit,
    package_vendor_tree_sha256 = freeze_source$vendor_tree_sha256,
    reference_provenance = list(
      artifact_series = "v102",
      candidate_package_version = "1.0.2",
      note = "The v102 token identifies the frozen historical candidate-result series."
    ),
    paper_workflow_commit = paper_commit,
    gate_max_abs_L_diff = gate_max_abs_L_diff,
    formula_id = lee_meta$formula_id,
    permutation = list(scheme = "global_joint_shuffle", use_blocks = FALSE,
                       permutations = perms, seed = seed, rng = "L'Ecuyer-CMRG"),
    weight_provenance = scope_obj@grid[[grid_name]]$weight_provenance,
    clustering = list(pct_min = "q95", resolution = 0.1,
                      consensus_threshold = 0.95, n_restart = 1000,
                      modules = n_modules, assigned = n_assigned),
    inputs = list(xenium_outs = input_dir, roi_file = roi_file,
                  sha256 = as.list(vapply(raw_candidates, sha256_file, character(1)))),
    outputs = list(sha256 = list(
      membership = sha256_file(membership_path),
      scope = sha256_file(scope_path)
    ))
  ), file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = 17)
  message(grid_name, " completed: modules=", n_modules, ", assigned=", n_assigned)
}

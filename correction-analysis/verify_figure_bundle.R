#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: verify_figure_bundle.R OUTPUT_DIR P5|LN")
}

script_dir <- local({
  x <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", x[[1L]])))
})
paper_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(paper_root, "main-text-scripts", "freeze_helpers.R"))

output_root <- normalizePath(args[[1L]], mustWork = TRUE)
sample_id <- toupper(args[[2L]])
if (!sample_id %in% c("P5", "LN")) stop("Sample must be P5 or LN.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.")

manifest_path <- file.path(output_root, paste0(sample_id, "_figure_manifest.json"))
if (!file.exists(manifest_path)) stop("Missing figure manifest: ", manifest_path)
manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
scalar <- function(x) as.character(unlist(x, use.names = FALSE)[[1L]])

if (!identical(scalar(manifest$sample_id), sample_id) ||
    !identical(scalar(manifest$package_version), "1.0.2")) {
  stop("Figure manifest sample or package version mismatch.")
}

git_status <- system2(
  "git", c("-C", shQuote(paper_root), "status", "--porcelain", "--untracked-files=all"),
  stdout = TRUE, stderr = TRUE
)
if (!is.null(attr(git_status, "status")) || length(git_status)) {
  stop("Refusing to verify a figure bundle from a dirty paper worktree.")
}
paper_commit <- system2(
  "git", c("-C", shQuote(paper_root), "rev-parse", "HEAD"),
  stdout = TRUE, stderr = TRUE
)
if (!is.null(attr(paper_commit, "status")) || length(paper_commit) != 1L ||
    !grepl("^[0-9a-f]{40}$", paper_commit)) {
  stop("Could not resolve the clean paper commit.")
}
if (!identical(scalar(manifest$paper_workflow_commit), paper_commit[[1L]])) {
  stop("Figure manifest paper commit does not match the current commit.")
}

vendor_dir <- file.path(paper_root, "docker", "genescope", "vendor", "geneSCOPE-v1.0.2")
source_commit <- trimws(readLines(file.path(vendor_dir, ".freeze-source-commit"), warn = FALSE))
vendor_tree <- trimws(readLines(file.path(vendor_dir, ".freeze-tree-sha256"), warn = FALSE))
if (!identical(scalar(manifest$package_source_commit), source_commit) ||
    !identical(scalar(manifest$package_vendor_tree_sha256), vendor_tree)) {
  stop("Figure manifest package provenance does not match the frozen vendor.")
}

workflow_path <- scalar(manifest$workflow$path)
mapping_path <- scalar(manifest$inputs$display_mapping$path)
if (!file.exists(workflow_path) ||
    !identical(sha256_file(workflow_path), scalar(manifest$workflow$sha256)) ||
    !file.exists(mapping_path) ||
    !identical(sha256_file(mapping_path),
               scalar(manifest$inputs$display_mapping$sha256))) {
  stop("Figure manifest workflow or display-mapping hash mismatch.")
}

input_root <- scalar(manifest$inputs$xenium_outs)
input_paths <- c(
  roi = scalar(manifest$inputs$roi_file),
  cell_feature_matrix = file.path(input_root, "cell_feature_matrix.h5"),
  cells_parquet = file.path(input_root, "cells.parquet"),
  transcripts_parquet = file.path(input_root, "transcripts.parquet")
)
recorded_inputs <- unlist(manifest$inputs$sha256, use.names = TRUE)
input_paths <- input_paths[names(input_paths) %in% names(recorded_inputs)]
if (!length(input_paths) || !all(file.exists(input_paths))) {
  stop("Figure manifest raw inputs are missing.")
}
observed_inputs <- vapply(input_paths, sha256_file, character(1L))
if (!identical(unname(observed_inputs), unname(recorded_inputs[names(input_paths)]))) {
  stop("Figure manifest raw-input hash mismatch.")
}

files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
files <- files[file.exists(files) & !dir.exists(files)]
files <- files[grepl("\\.(png|pdf|tsv|csv|rds)$", files, ignore.case = TRUE)]
files <- setdiff(files, manifest_path)
files <- sort(files, method = "radix")
relative <- substring(files, nchar(output_root) + 2L)
recorded_outputs <- unlist(manifest$outputs_sha256, use.names = TRUE)
if (!setequal(relative, names(recorded_outputs))) {
  stop("Figure manifest output inventory does not match the bundle.")
}
observed_outputs <- stats::setNames(vapply(files, sha256_file, character(1L)), relative)
if (!identical(unname(observed_outputs[names(recorded_outputs)]),
               unname(recorded_outputs))) {
  stop("Figure manifest output hash mismatch.")
}

if (identical(sample_id, "P5")) {
  cluster_col <- "q95_res0.1_grid30_log1p_freq0.95"
  required <- c(
    "LvsR/LvsR_grid30.png",
    file.path("grid30", "network", paste0("network_", cluster_col, ".png")),
    file.path("grid30", "network", paste0("dendro_network_", cluster_col, ".png")),
    "grid30/grid30_boundary.png",
    file.path("grid30", "heatmap", paste0("LeeL_heatmap_grid30_", cluster_col, ".png")),
    file.path("grid30", "idelta", paste0(cluster_col, "_idelta_by_cluster.png")),
    "P5_top_pairs_all.tsv", "P5_top_pairs_display_filter.tsv", "P5_Top6.tsv",
    "P5_dendro_path_audit.tsv"
  )
  assert_required_figure_outputs(output_root, sample_id, required, 205L)
  dendro <- utils::read.delim(file.path(output_root, "P5_dendro_path_audit.tsv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
  reference <- utils::read.delim(file.path(script_dir, "reference_p5_dendro_path.tsv"),
                                 stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(dendro) != 1L || nrow(reference) != 1L ||
      !identical(dendro$raw_module_path[[1L]], reference$raw_module_path[[1L]]) ||
      !identical(dendro$display_module_path[[1L]], reference$display_module_path[[1L]]) ||
      !isTRUE(dendro$stem_positioned_between[[1L]]) ||
      isTRUE(dendro$direct_stem_to_broad_adjacency[[1L]])) {
    stop("P5 dendrogram-path artifact does not match its frozen reference.")
  }
} else {
  cluster_col <- "q99.9_res0.1_grid30_log1p_freq0.95"
  required <- c(
    "LvsR/LvsR_grid30.png",
    file.path("grid30", "network", paste0("network_", cluster_col, ".png")),
    file.path("grid30", "network", paste0("dendro_network_", cluster_col, ".png")),
    "grid30/grid30_boundary.png",
    "LN_top_pairs_all.tsv", "LN_top_pairs_display_filter.tsv", "LN_Top6.tsv"
  )
  assert_required_figure_outputs(output_root, sample_id, required, 8L)
}

top6 <- utils::read.delim(file.path(output_root, paste0(sample_id, "_Top6.tsv")),
                          stringsAsFactors = FALSE, check.names = FALSE)
assert_reference_top6(script_dir, sample_id, top6)
message(sample_id, " frozen figure bundle verification passed: ", output_root)

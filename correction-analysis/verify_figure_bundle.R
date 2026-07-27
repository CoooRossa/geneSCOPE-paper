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
suppressPackageStartupMessages(library(geneSCOPE))

output_root <- normalizePath(args[[1L]], mustWork = TRUE)
sample_id <- toupper(args[[2L]])
if (!sample_id %in% c("P5", "LN")) stop("Sample must be P5 or LN.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.")
if (!identical(as.character(utils::packageVersion("geneSCOPE")), "1.0.2")) {
  stop("The verifier must run with the frozen geneSCOPE 1.0.2 package.")
}

manifest_path <- file.path(output_root, paste0(sample_id, "_figure_manifest.json"))
if (!file.exists(manifest_path)) stop("Missing figure manifest: ", manifest_path)
manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
scalar <- function(x) as.character(unlist(x, use.names = FALSE)[[1L]])

if (!identical(scalar(manifest$sample_id), sample_id) ||
    !identical(scalar(manifest$package_version), "1.0.2")) {
  stop("Figure manifest sample or package version mismatch.")
}
if (!identical(scalar(manifest$parameters$analysis_mode),
               "render_from_hash_pinned_authoritative_results") ||
    !identical(scalar(manifest$parameters$display_derivations),
               "computeDensity_only")) {
  stop("Figure manifest does not declare the frozen render-only policy.")
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
workflow_text <- paste(readLines(workflow_path, warn = FALSE), collapse = "\n")
forbidden_analysis_calls <- c(
  "computeL", "computeLvsRCurve", "clusterGenes", "getTopLvsR", "computeIDelta"
)
forbidden_pattern <- paste0(
  "\\b(", paste(forbidden_analysis_calls, collapse = "|"), ")\\s*\\("
)
if (grepl(forbidden_pattern, workflow_text, perl = TRUE)) {
  stop("Correction figure workflow contains a forbidden analysis recomputation call.")
}

expected_sources <- if (identical(sample_id, "P5")) {
  list(
    scope = c("MAIN_RESULTS", "P5/P5_scope_shuffle_v102.rds"),
    top_pairs = c("MAIN_RESULTS", "P5/P5_toplvsr_all_shuffleFDR.tsv"),
    analysis_manifest = c("MAIN_RESULTS", "P5/manifest.json")
  )
} else {
  list(
    scope = c("MAIN_RESULTS", "LN/LN_scope_shuffle_v102.rds"),
    top_pairs = c("LN_COMPLETE_RESULTS", "LN_top_pairs_complete_delta_v102.tsv"),
    analysis_manifest = c("MAIN_RESULTS", "LN/manifest.json"),
    delta_manifest = c("LN_COMPLETE_RESULTS", "manifest.json")
  )
}
analysis_sources <- manifest$inputs$analysis_sources
if (!is.list(analysis_sources) ||
    !setequal(names(analysis_sources), names(expected_sources))) {
  stop("Figure manifest is missing the frozen analysis sources.")
}
verified_sources <- lapply(names(expected_sources), function(source_name) {
  source <- analysis_sources[[source_name]]
  expected <- expected_sources[[source_name]]
  if (!identical(scalar(source$root_token), expected[[1L]]) ||
      !identical(scalar(source$relative_path), expected[[2L]])) {
    stop("Unexpected ", source_name, " analysis source in figure manifest.")
  }
  verified <- assert_reference_artifact(
    script_dir, expected[[1L]], expected[[2L]], scalar(source$path)
  )
  if (!identical(verified$sha256, scalar(source$sha256))) {
    stop("Recorded ", source_name, " SHA256 does not match its reference.")
  }
  verified
})
names(verified_sources) <- names(expected_sources)

authoritative_cluster_col <- if (identical(sample_id, "P5")) {
  "shuffle_q95_res0.1_grid30"
} else {
  "shuffle_q99.9_res0.1_grid30"
}
input_root <- scalar(manifest$inputs$xenium_outs)
scope_obj <- readRDS(verified_sources$scope$path)
assert_authoritative_scope(
  scope_obj, sample_id, "grid30", "LR_curve_30_shuffle",
  authoritative_cluster_col
)
assert_scope_xenium_identity(
  scope_obj, input_root, scalar(manifest$inputs$roi_file), sample_id
)
source_membership <- assert_reference_membership(
  script_dir, sample_id, rownames(scope_obj@meta.data),
  scope_obj@meta.data[[authoritative_cluster_col]]
)
membership_fields <- c(
  "exact_gene_module_sha256", "partition_sha256", "assigned", "modules"
)
for (field in membership_fields) {
  if (!identical(as.character(source_membership[[field]]),
                 scalar(manifest$membership[[field]]))) {
    stop("Figure manifest membership does not match its source: ", field)
  }
}
raw_membership <- as.character(scope_obj@meta.data[[authoritative_cluster_col]])
assigned_mask <- !is.na(raw_membership) & nzchar(raw_membership) & raw_membership != "-1"
assigned_genes <- rownames(scope_obj@meta.data)[assigned_mask]
display_mapping <- utils::read.delim(mapping_path, stringsAsFactors = FALSE,
                                     check.names = FALSE)
display_membership <- apply_display_mapping(raw_membership, display_mapping)
display_cluster_ids <- sort(unique(na.omit(display_membership)))

expected_pair_rows <- if (identical(sample_id, "P5")) 1578L else 75405L
source_top_pairs <- read_authoritative_top_pairs(
  verified_sources$top_pairs$path, sample_id, expected_pair_rows
)
assert_scope_pair_table(scope_obj, source_top_pairs, sample_id)
source_display_pairs <- filter_display_pairs(source_top_pairs)
assert_reference_top6(script_dir, sample_id, source_display_pairs)
assert_analysis_provenance(
  sample_id, verified_sources$analysis_manifest$path, expected_pair_rows,
  delta_manifest_path = if (identical(sample_id, "LN")) {
    verified_sources$delta_manifest$path
  } else {
    NULL
  }
)
rm(scope_obj)
invisible(gc())

raw_gate <- assert_reference_raw_inputs(
  script_dir, sample_id, input_root, scalar(manifest$inputs$roi_file)
)
recorded_inputs <- unlist(manifest$inputs$sha256, use.names = TRUE)
observed_inputs <- unlist(raw_gate$sha256, use.names = TRUE)
if (!identical(names(recorded_inputs), names(observed_inputs)) ||
    !identical(unname(observed_inputs), unname(recorded_inputs))) {
  stop("Figure manifest raw-input hash mismatch.")
}

files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
files <- files[file.exists(files) & !dir.exists(files)]
files <- files[grepl("\\.(png|pdf|tsv|csv|rds)$", files, ignore.case = TRUE)]
files <- setdiff(files, manifest_path)
files <- sort(files, method = "radix")
relative <- substring(files, nchar(output_root) + 2L)
keep <- !grepl("^\\.geneSCOPE-v1\\.0\\.2-library/", relative)
files <- files[keep]
relative <- relative[keep]
recorded_outputs <- unlist(manifest$outputs_sha256, use.names = TRUE)
if (!setequal(relative, names(recorded_outputs))) {
  stop("Figure manifest output inventory does not match the bundle.")
}
observed_outputs <- stats::setNames(vapply(files, sha256_file, character(1L)), relative)
if (!identical(unname(observed_outputs[names(recorded_outputs)]),
               unname(recorded_outputs))) {
  stop("Figure manifest output hash mismatch.")
}
all_pairs_relative <- paste0(sample_id, "_top_pairs_all.tsv")
if (!all_pairs_relative %in% names(observed_outputs) ||
    !identical(unname(observed_outputs[[all_pairs_relative]]),
               verified_sources$top_pairs$sha256)) {
  stop("Bundled all-pair table is not the authoritative source table.")
}

assert_same_pair_table <- function(observed, expected, label) {
  observed <- as.data.frame(observed, stringsAsFactors = FALSE, check.names = FALSE)
  expected <- as.data.frame(expected, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(observed), names(expected)) || nrow(observed) != nrow(expected)) {
    stop(label, " schema or row count differs from the authoritative derivation.")
  }
  for (column in names(expected)) {
    if (is.numeric(expected[[column]])) {
      if (!isTRUE(all.equal(as.numeric(observed[[column]]), expected[[column]],
                            tolerance = 1e-12, check.attributes = FALSE))) {
        stop(label, " numeric column differs: ", column)
      }
    } else if (!identical(as.character(observed[[column]]),
                          as.character(expected[[column]]))) {
      stop(label, " text column differs: ", column)
    }
  }
  invisible(TRUE)
}

bundled_display_pairs <- utils::read.delim(
  file.path(output_root, paste0(sample_id, "_top_pairs_display_filter.tsv")),
  stringsAsFactors = FALSE, check.names = FALSE
)
bundled_top6 <- utils::read.delim(
  file.path(output_root, paste0(sample_id, "_Top6.tsv")),
  stringsAsFactors = FALSE, check.names = FALSE
)
assert_same_pair_table(bundled_display_pairs, source_display_pairs,
                       paste0(sample_id, " display-filter table"))
assert_same_pair_table(bundled_top6, utils::head(source_display_pairs, 6L),
                       paste0(sample_id, " Top6 table"))

cluster_col <- if (identical(sample_id, "P5")) {
  "q95_res0.1_grid30_log1p_freq0.95"
} else {
  "q99.9_res0.1_grid30_log1p_freq0.95"
}
network_paths <- c(
  file.path("grid30", "network", paste0("network_", cluster_col, ".png")),
  file.path("grid30", "network", paste0("dendro_network_", cluster_col, ".png"))
)
if (identical(sample_id, "P5")) {
  display_for_assigned <- as.character(display_membership[assigned_mask])
  ord <- order(suppressWarnings(as.integer(display_for_assigned)), assigned_genes,
               method = "radix")
  cluster_genes <- assigned_genes[ord]
  cluster_density_dir <- file.path(
    "grid30", "density",
    paste0("clusters_", paste(as.character(display_cluster_ids), collapse = "_"))
  )
  cluster_density_paths <- file.path(
    cluster_density_dir, paste0(cluster_genes, "_density.png")
  )
  top20 <- utils::head(source_display_pairs, 20L)
  top_delta_paths <- file.path(
    "TopDelta",
    paste0(top20$gene1, "_", top20$gene2, "_grid30_L_0.1.1.png")
  )
  png_paths <- c(
    "LvsR/LvsR_grid30.png", network_paths,
    "grid30/density/CEACAM5_ACTA2_grid30.png",
    "grid30/density/CEACAM5_CEACAM6_grid30.png",
    "cells/density/CEACAM5_CEACAM6_centroids_grid_cells.png",
    "cells/density/CEACAM5_ACTA2_centroids_grid_cells.png",
    "grid30/grid30_boundary.png", cluster_density_paths, top_delta_paths,
    file.path("grid30", "heatmap",
              paste0("LeeL_heatmap_grid30_", cluster_col, ".png")),
    file.path("grid30", "idelta", paste0(cluster_col, "_idelta_by_cluster.png"))
  )
  non_png_paths <- c(
    "P5_top_pairs_all.tsv", "P5_top_pairs_display_filter.tsv", "P5_Top6.tsv",
    "P5_dendro_path_audit.tsv"
  )
} else {
  png_paths <- c(
    "LvsR/LvsR_grid30.png", network_paths,
    "grid30/density/ITGB2_PDGFRA_grid30.png",
    "grid30/density/ITGB2_PTPN6_grid30.png",
    "cells/density/ITGB2_PDGFRA_centroids_grid_cells.png",
    "cells/density/ITGB2_PTPN6_centroids_grid_cells.png",
    "grid30/grid30_boundary.png"
  )
  non_png_paths <- c(
    "LN_top_pairs_all.tsv", "LN_top_pairs_display_filter.tsv", "LN_Top6.tsv"
  )
}
expected_inventory <- sort(c(png_paths, non_png_paths), method = "radix")
if (!identical(sort(relative, method = "radix"), expected_inventory)) {
  stop(
    sample_id, " figure bundle inventory differs from the exact frozen set; missing=",
    paste(setdiff(expected_inventory, relative), collapse = ", "), "; unexpected=",
    paste(setdiff(relative, expected_inventory), collapse = ", ")
  )
}

if (!requireNamespace("png", quietly = TRUE)) stop("png is required for image decoding.")
expected_png_dimensions <- stats::setNames(
  rep(list(c(height = 3000L, width = 3000L)), length(png_paths)), png_paths
)
expected_png_dimensions[["LvsR/LvsR_grid30.png"]] <- c(height = 3600L, width = 3600L)
for (path in network_paths) {
  expected_png_dimensions[[path]] <- c(height = 9000L, width = 9000L)
}
if (identical(sample_id, "P5")) {
  expected_png_dimensions[[file.path(
    "grid30", "heatmap", paste0("LeeL_heatmap_grid30_", cluster_col, ".png")
  )]] <- c(height = 12000L, width = 12000L)
  expected_png_dimensions[[file.path(
    "grid30", "idelta", paste0(cluster_col, "_idelta_by_cluster.png")
  )]] <- c(height = 4200L, width = 9000L)
}
for (path in names(expected_png_dimensions)) {
  decoded <- png::readPNG(file.path(output_root, path), native = TRUE)
  observed_dim <- dim(decoded)
  expected_dim <- unname(expected_png_dimensions[[path]])
  if (length(observed_dim) < 2L ||
      !identical(as.integer(observed_dim[1:2]), as.integer(expected_dim))) {
    stop(path, " did not decode at its frozen dimensions: ",
         paste(observed_dim[1:2], collapse = "x"))
  }
  rm(decoded)
}
invisible(gc())

if (identical(sample_id, "P5")) {
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

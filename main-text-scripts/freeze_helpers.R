genescope_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) return(normalizePath(getwd()))
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
}

required_directory_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set ", name, " to the Xenium outs directory.")
  value <- normalizePath(value, mustWork = TRUE)
  if (!dir.exists(value)) stop(name, " is not a directory: ", value)
  value
}

integer_env <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value) || value < 1L) stop(name, " must be a positive integer.")
  value
}

reset_freeze_rng <- function(seed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed))
}

configure_freeze_runtime <- function() {
  max_bytes <- 500000 * 1024^2
  options(
    future.globals.maxSize = max_bytes,
    genescope.future.globals.maxSize = max_bytes
  )
  invisible(max_bytes)
}

sha256_file <- function(path) {
  if (!file.exists(path)) return("")
  sha256sum <- Sys.which("sha256sum")
  shasum <- Sys.which("shasum")
  if (nzchar(sha256sum)) {
    line <- system2(sha256sum, shQuote(path), stdout = TRUE, stderr = TRUE)
  } else if (nzchar(shasum)) {
    line <- system2(shasum, c("-a", "256", shQuote(path)),
                    stdout = TRUE, stderr = TRUE)
  } else {
    stop("sha256sum or shasum is required for freeze manifests.")
  }
  status <- attr(line, "status", exact = TRUE)
  digest <- if (length(line)) strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]] else ""
  if ((!is.null(status) && status != 0L) || !grepl("^[0-9a-fA-F]{64}$", digest)) {
    stop("Could not hash: ", path)
  }
  tolower(digest)
}

sha256_text <- function(text) {
  path <- tempfile("genescope-freeze-text-")
  on.exit(unlink(path), add = TRUE)
  writeLines(enc2utf8(as.character(text)), path, useBytes = TRUE)
  sha256_file(path)
}

canonical_lee_s2_gate <- function(tolerance = 1e-12) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required for the independent Lee's L gate.")
  }
  Xz <- scale(matrix(c(
    1, 4, 2, 5,
    2, 1, 4, 3,
    5, 2, 1, 4
  ), nrow = 4), center = TRUE, scale = FALSE)
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 1),
    j = c(2, 1, 3, 2, 4, 3, 1, 4),
    x = c(1, 1, 2, 2, 1, 1, 3, 3),
    dims = c(4, 4)
  )
  observed <- getFromNamespace("lee_L", "geneSCOPE")(Xz, W, 1L)
  Wz <- as.matrix(W %*% Xz)
  S2 <- sum(Matrix::rowSums(W)^2)
  norms <- sqrt(colSums(Xz^2))
  expected <- nrow(Xz) / S2 * crossprod(Wz) / outer(norms, norms)
  max_abs_error <- max(abs(observed - expected), na.rm = TRUE)
  if (!is.finite(max_abs_error) || max_abs_error > tolerance) {
    stop("Independent Lee's L/S2 gate failed: max_abs_error=", max_abs_error)
  }
  max_abs_error
}

require_freeze_source_metadata <- function() {
  source_commit <- Sys.getenv("GENESCOPE_SOURCE_COMMIT", unset = "")
  vendor_tree_sha256 <- Sys.getenv("GENESCOPE_VENDOR_TREE_SHA256", unset = "")
  if (!grepl("^[0-9a-f]{40}$", source_commit)) {
    stop("GENESCOPE_SOURCE_COMMIT must identify the frozen 40-character commit.")
  }
  if (!grepl("^[0-9a-f]{64}$", vendor_tree_sha256)) {
    stop("GENESCOPE_VENDOR_TREE_SHA256 must identify the frozen source tree.")
  }
  list(source_commit = source_commit, vendor_tree_sha256 = vendor_tree_sha256)
}

require_paper_commit_metadata <- function() {
  paper_commit <- Sys.getenv("GENESCOPE_PAPER_COMMIT", unset = "")
  if (!grepl("^[0-9a-f]{40}$", paper_commit)) {
    stop("GENESCOPE_PAPER_COMMIT must identify the frozen 40-character paper commit.")
  }
  paper_commit
}

assert_fresh_output_dir <- function(path,
                                    allowed_entries = ".geneSCOPE-v1.0.2-library") {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path <- normalizePath(path, mustWork = TRUE)
  entries <- list.files(path, all.files = TRUE, no.. = TRUE)
  unexpected <- setdiff(entries, allowed_entries)
  if (length(unexpected)) {
    stop(
      "Figure workflow requires a fresh output directory; found: ",
      paste(unexpected, collapse = ", ")
    )
  }
  allowed_present <- intersect(entries, allowed_entries)
  if (length(allowed_present) &&
      any(!dir.exists(file.path(path, allowed_present)))) {
    stop("Reserved freeze-output entries must be directories.")
  }
  path
}

read_display_mapping <- function(script_dir, sample_id) {
  filename <- paste0(sample_id, "_display_mapping.tsv")
  candidates <- c(
    file.path(script_dir, "display-mappings", filename),
    file.path(script_dir, "..", "correction-analysis", "display-mappings", filename)
  )
  candidates <- candidates[file.exists(candidates)]
  path <- if (length(candidates)) candidates[[1L]] else ""
  if (!file.exists(path)) stop("Missing frozen display mapping: ", path)
  mapping <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("display_module", "current_module", "color")
  missing <- setdiff(required, names(mapping))
  if (length(missing)) stop("Display mapping is missing: ", paste(missing, collapse = ", "))
  if (anyDuplicated(mapping$current_module) || anyDuplicated(mapping$display_module)) {
    stop("Display mapping must be one-to-one: ", path)
  }
  mapping
}

display_palette <- function(mapping) {
  display_id <- suppressWarnings(as.integer(mapping$display_module))
  if (any(!is.finite(display_id))) stop("Display module IDs must be integers.")
  ord <- order(display_id, method = "radix")
  stats::setNames(as.character(mapping$color[ord]), as.character(display_id[ord]))
}

apply_display_mapping <- function(membership, mapping) {
  raw <- as.character(membership)
  lookup <- stats::setNames(as.character(mapping$display_module),
                            as.character(mapping$current_module))
  assigned <- !is.na(raw) & nzchar(raw) & raw != "-1"
  unknown <- sort(unique(raw[assigned & !raw %in% names(lookup)]))
  if (length(unknown)) {
    stop("Unmapped current module(s): ", paste(unknown, collapse = ", "))
  }
  display <- rep(NA_character_, length(raw))
  display[assigned] <- unname(lookup[raw[assigned]])
  factor(display, levels = as.character(sort(as.integer(mapping$display_module))))
}

membership_digests <- function(genes, membership) {
  genes <- as.character(genes)
  membership <- as.character(membership)
  if (length(genes) != length(membership) || anyNA(genes) || any(!nzchar(genes))) {
    stop("Membership genes and module labels must be aligned and non-empty.")
  }
  if (anyDuplicated(genes)) stop("Membership genes must be unique.")
  membership[is.na(membership) | !nzchar(membership)] <- "-1"

  gene_order <- order(genes, method = "radix")
  exact_lines <- paste(genes[gene_order], membership[gene_order], sep = "\t")
  assigned <- membership != "-1"
  groups <- split(genes[assigned], membership[assigned])
  group_lines <- if (length(groups)) {
    sort(vapply(
      groups,
      function(x) paste(sort(x, method = "radix"), collapse = ","),
      character(1L)
    ), method = "radix")
  } else {
    character()
  }
  unassigned <- paste(sort(genes[!assigned], method = "radix"), collapse = ",")
  partition_lines <- c(paste0("M\t", group_lines), paste0("U\t", unassigned))

  list(
    exact_gene_module_sha256 = sha256_text(exact_lines),
    partition_sha256 = sha256_text(partition_lines),
    assigned = sum(assigned),
    modules = length(groups)
  )
}

assert_reference_membership <- function(script_dir, sample_id, genes, membership) {
  candidates <- c(
    file.path(script_dir, "reference_membership_digests.tsv"),
    file.path(script_dir, "..", "correction-analysis", "reference_membership_digests.tsv")
  )
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path) || !file.exists(path)) stop("Missing membership digest reference.")
  reference <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  reference <- reference[reference$sample_id == sample_id, , drop = FALSE]
  if (nrow(reference) != 1L) stop("Missing unique membership reference for ", sample_id)
  observed <- membership_digests(genes, membership)
  checks <- c(
    identical(observed$exact_gene_module_sha256, reference$exact_gene_module_sha256[[1L]]),
    identical(observed$partition_sha256, reference$partition_sha256[[1L]]),
    identical(as.integer(observed$assigned), as.integer(reference$assigned[[1L]])),
    identical(as.integer(observed$modules), as.integer(reference$modules[[1L]]))
  )
  if (!all(checks)) {
    stop(
      sample_id, " membership freeze gate failed: exact=", observed$exact_gene_module_sha256,
      ", partition=", observed$partition_sha256,
      ", assigned=", observed$assigned, ", modules=", observed$modules
    )
  }
  observed
}

assert_reference_top6 <- function(script_dir, sample_id, x) {
  path <- file.path(script_dir, "reference_top6.tsv")
  if (!file.exists(path)) stop("Missing Top6 reference: ", path)
  reference <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  reference <- reference[reference$sample_id == sample_id, , drop = FALSE]
  reference <- reference[order(reference$rank), , drop = FALSE]
  if (nrow(reference) != 6L || nrow(x) < 6L) {
    stop(sample_id, " Top6 freeze gate requires six observed and six reference pairs.")
  }
  pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "--")
  observed_key <- pair_key(as.character(x$gene1[seq_len(6L)]),
                           as.character(x$gene2[seq_len(6L)]))
  reference_key <- pair_key(reference$gene1, reference$gene2)
  if (!identical(observed_key, reference_key)) {
    stop(sample_id, " Top6 freeze gate failed: ", paste(observed_key, collapse = ", "))
  }
  invisible(observed_key)
}

write_freeze_output_manifest <- function(output_root, sample_id, freeze_source,
                                         formula_gate, membership_gate,
                                         input_dir, roi_file,
                                         display_mapping_path, workflow_path,
                                         parameters) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to write the figure freeze manifest.")
  }
  output_root <- normalizePath(output_root, mustWork = TRUE)
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  roi_file <- normalizePath(roi_file, mustWork = TRUE)
  display_mapping_path <- normalizePath(display_mapping_path, mustWork = TRUE)
  workflow_path <- normalizePath(workflow_path, mustWork = TRUE)
  paper_commit <- require_paper_commit_metadata()
  raw_inputs <- c(
    roi = roi_file,
    cell_feature_matrix = file.path(input_dir, "cell_feature_matrix.h5"),
    cells_parquet = file.path(input_dir, "cells.parquet"),
    transcripts_parquet = file.path(input_dir, "transcripts.parquet")
  )
  raw_inputs <- raw_inputs[file.exists(raw_inputs)]
  if (!"roi" %in% names(raw_inputs) || length(raw_inputs) < 2L) {
    stop("Could not identify the ROI and at least one frozen Xenium input file.")
  }
  input_hashes <- stats::setNames(
    as.list(vapply(raw_inputs, sha256_file, character(1L))), names(raw_inputs)
  )
  manifest_path <- file.path(output_root, paste0(sample_id, "_figure_manifest.json"))
  files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  files <- files[grepl("\\.(png|pdf|tsv|csv|rds)$", files, ignore.case = TRUE)]
  files <- setdiff(files, manifest_path)
  files <- sort(files, method = "radix")
  relative <- substring(files, nchar(output_root) + 2L)
  hashes <- if (length(files)) {
    stats::setNames(as.list(vapply(files, sha256_file, character(1L))), relative)
  } else {
    list()
  }
  manifest <- list(
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    sample_id = sample_id,
    package_version = as.character(utils::packageVersion("geneSCOPE")),
    package_source_commit = freeze_source$source_commit,
    package_vendor_tree_sha256 = freeze_source$vendor_tree_sha256,
    paper_workflow_commit = paper_commit,
    workflow = list(path = workflow_path, sha256 = sha256_file(workflow_path)),
    gate_max_abs_L_diff = formula_gate,
    membership = membership_gate,
    parameters = parameters,
    inputs = list(
      xenium_outs = input_dir,
      roi_file = roi_file,
      sha256 = input_hashes,
      display_mapping = list(
        path = display_mapping_path,
        sha256 = sha256_file(display_mapping_path)
      )
    ),
    outputs_sha256 = hashes
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE,
                       pretty = TRUE, digits = 17)
  invisible(manifest_path)
}

filter_display_pairs <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  required <- c("gene1", "gene2", "L", "r", "pct1", "pct2")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("getTopLvsR output is missing: ", paste(missing, collapse = ", "))

  q_col <- intersect(c("q_Delta", "delta_fdr", "fdr"), names(x))
  if (!length(q_col)) stop("getTopLvsR output has no Delta-adjusted p-value column.")
  x$q_Delta <- as.numeric(x[[q_col[[1L]]]])
  x$Delta <- as.numeric(x$L) - as.numeric(x$r)

  keep <- is.finite(x$q_Delta) & x$q_Delta < 0.05 &
    is.finite(x$L) & x$L > 0 &
    is.finite(x$r) & x$r < 0.05 &
    is.finite(x$pct1) & x$pct1 > 20 &
    is.finite(x$pct2) & x$pct2 > 20
  x <- x[keep, , drop = FALSE]
  x[order(-x$Delta, x$gene1, x$gene2, method = "radix"), , drop = FALSE]
}

assert_complete_delta_universe <- function(x) {
  provenance <- attr(x, "permutation_provenance")
  ok <- is.list(provenance) && identical(provenance$p_adj_mode, "BH_universe") &&
    identical(provenance$use_blocks, FALSE) &&
    is.finite(provenance$total_universe) && is.finite(provenance$selected_pairs) &&
    provenance$total_universe == provenance$selected_pairs
  if (!ok) {
    stop(
      "Delta inference must cover the complete eligible universe; increase top_n. ",
      "eligible=", provenance$total_universe, ", selected=", provenance$selected_pairs
    )
  }
  invisible(provenance)
}

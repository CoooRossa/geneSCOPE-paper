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

required_file_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set ", name, " to the frozen input file.")
  value <- normalizePath(value, mustWork = TRUE)
  if (dir.exists(value)) stop(name, " must be a file: ", value)
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

bundle_file_inventory <- function(output_root, manifest_path = NULL) {
  output_root <- normalizePath(output_root, mustWork = TRUE)
  entries <- list.files(
    output_root, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE, include.dirs = TRUE
  )
  relative <- substring(entries, nchar(output_root) + 2L)
  keep <- !grepl("^\\.geneSCOPE-v1\\.0\\.2-library(/|$)", relative)
  entries <- entries[keep]
  relative <- relative[keep]
  symbolic_links <- nzchar(Sys.readlink(entries))
  if (any(symbolic_links)) {
    stop("Figure bundle contains symbolic-link outputs: ",
         paste(relative[symbolic_links], collapse = ", "))
  }
  keep <- file.exists(entries) & !dir.exists(entries)
  entries <- entries[keep]
  relative <- relative[keep]
  if (!is.null(manifest_path)) {
    manifest_relative <- substring(
      normalizePath(manifest_path, mustWork = FALSE), nchar(output_root) + 2L
    )
    keep <- relative != manifest_relative
    entries <- entries[keep]
    relative <- relative[keep]
  }
  ord <- order(relative, method = "radix")
  list(files = entries[ord], relative = relative[ord])
}

reference_table_dir <- function(script_dir) {
  candidates <- unique(c(
    script_dir,
    file.path(script_dir, "correction-analysis"),
    file.path(script_dir, "..", "correction-analysis")
  ))
  hit <- candidates[file.exists(file.path(candidates, "reference_artifact_sha256.tsv"))]
  if (!length(hit)) stop("Could not locate the correction-analysis reference tables.")
  normalizePath(hit[[1L]], mustWork = TRUE)
}

assert_reference_artifact <- function(script_dir, root_token, relative_path, path) {
  reference_dir <- reference_table_dir(script_dir)
  tables <- c(
    file.path(reference_dir, "reference_artifact_sha256.tsv"),
    file.path(reference_dir, "reference_external_bundle_hashes.tsv")
  )
  tables <- tables[file.exists(tables)]
  reference <- do.call(rbind, lapply(tables, function(table_path) {
    utils::read.delim(table_path, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  key <- reference$root_token == root_token & reference$relative_path == relative_path
  if (sum(key) != 1L) {
    stop("Expected one artifact reference for ", root_token, "/", relative_path, ".")
  }
  path <- normalizePath(path, mustWork = TRUE)
  if (dir.exists(path)) stop("Frozen artifact must be a file: ", path)
  observed <- sha256_file(path)
  expected <- tolower(reference$sha256[key][[1L]])
  if (!identical(observed, expected)) {
    stop(
      "Frozen artifact SHA256 mismatch for ", root_token, "/", relative_path,
      ": observed=", observed, "; expected=", expected
    )
  }
  list(
    path = path,
    root_token = root_token,
    relative_path = relative_path,
    sha256 = observed
  )
}

assert_reference_raw_inputs <- function(script_dir, sample_id, input_dir, roi_file) {
  reference_path <- file.path(reference_table_dir(script_dir),
                              "reference_raw_input_sha256.tsv")
  if (!file.exists(reference_path)) stop("Missing frozen raw-input SHA256 reference.")
  reference <- utils::read.delim(reference_path, stringsAsFactors = FALSE,
                                 check.names = FALSE)
  reference <- reference[reference$sample_id == sample_id, , drop = FALSE]
  expected_names <- c("roi", "cell_feature_matrix", "cells_parquet",
                      "transcripts_parquet")
  if (nrow(reference) != length(expected_names) ||
      !setequal(reference$input_name, expected_names) ||
      anyDuplicated(reference$input_name)) {
    stop("Raw-input reference is incomplete for ", sample_id, ".")
  }
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  roi_file <- normalizePath(roi_file, mustWork = TRUE)
  paths <- c(
    roi = roi_file,
    cell_feature_matrix = file.path(input_dir, "cell_feature_matrix.h5"),
    cells_parquet = file.path(input_dir, "cells.parquet"),
    transcripts_parquet = file.path(input_dir, "transcripts.parquet")
  )
  if (!all(file.exists(paths)) || any(dir.exists(paths))) {
    stop(sample_id, " requires ROI, H5, cells.parquet, and transcripts.parquet.")
  }
  observed <- vapply(paths, sha256_file, character(1L))
  expected <- stats::setNames(tolower(reference$sha256), reference$input_name)
  if (!identical(unname(observed[expected_names]), unname(expected[expected_names]))) {
    mismatch <- expected_names[observed[expected_names] != expected[expected_names]]
    stop(sample_id, " raw-input SHA256 mismatch: ", paste(mismatch, collapse = ", "))
  }
  list(paths = as.list(paths), sha256 = as.list(observed))
}

assert_authoritative_scope <- function(scope_obj, sample_id, grid_name,
                                       curve_name, cluster_col) {
  if (!inherits(scope_obj, "scope_object")) {
    stop(sample_id, " frozen analysis source is not a scope_object.")
  }
  if (!grid_name %in% names(scope_obj@grid) ||
      is.null(scope_obj@stats[[grid_name]]$LeeStats_Xz)) {
    stop(sample_id, " frozen analysis source is missing ", grid_name, "/LeeStats_Xz.")
  }
  lee <- scope_obj@stats[[grid_name]]$LeeStats_Xz
  required_layers <- c("L", "P", "FDR", curve_name, cluster_col)
  missing_layers <- required_layers[vapply(required_layers, function(layer) {
    is.null(lee[[layer]])
  }, logical(1L))]
  if (length(missing_layers)) {
    stop(sample_id, " frozen analysis source is missing: ",
         paste(missing_layers, collapse = ", "))
  }
  meta_genes <- rownames(scope_obj@meta.data)
  if (!length(meta_genes) || anyNA(meta_genes) || any(!nzchar(meta_genes)) ||
      anyDuplicated(meta_genes)) {
    stop(sample_id, " frozen analysis source has invalid gene metadata.")
  }
  if (!cluster_col %in% colnames(scope_obj@meta.data)) {
    stop(sample_id, " frozen analysis source is missing membership: ", cluster_col)
  }
  if (!inherits(lee[[cluster_col]], "igraph")) {
    stop(sample_id, " frozen analysis source is missing the authoritative graph: ",
         cluster_col)
  }
  meta <- lee$meta
  formula_ok <- is.list(meta) && length(meta$formula_id) == 1L &&
    meta$formula_id %in% c("Lee2009_S2_v1", "Lee_S2_v1")
  permutation_ok <- is.list(meta) && is.character(meta$permutation_scheme) &&
    length(meta$permutation_scheme) == 1L &&
    grepl("all-grid joint shuffle|global_joint_shuffle", meta$permutation_scheme)
  if (!formula_ok || !permutation_ok) {
    stop(sample_id, " frozen analysis source failed formula/permutation provenance.")
  }
  genes <- rownames(lee$L)
  if (!length(genes) || anyNA(genes) || anyDuplicated(genes) ||
      !all(genes %in% meta_genes)) {
    stop(sample_id, " frozen Lee matrices have invalid gene names.")
  }
  matrices <- lapply(c("L", "P", "FDR"), function(layer) lee[[layer]])
  matrix_ok <- vapply(matrices, function(x) {
    identical(dim(x), c(length(genes), length(genes))) &&
      identical(rownames(x), genes) && identical(colnames(x), genes)
  }, logical(1L))
  if (!all(matrix_ok)) {
    stop(sample_id, " frozen L/P/FDR matrices are not aligned to gene metadata.")
  }
  curve <- lee[[curve_name]]
  if (!is.data.frame(curve) || nrow(curve) != 1000L ||
      !all(c("Pear", "fit", "lo95", "hi95") %in% names(curve))) {
    stop(sample_id, " frozen L-vs-r curve is malformed: ", curve_name)
  }
  required_cells <- c("counts", "logCPM")
  if (!all(required_cells %in% names(scope_obj@cells))) {
    stop(sample_id, " frozen analysis source is missing cell layers: ",
         paste(setdiff(required_cells, names(scope_obj@cells)), collapse = ", "))
  }
  if (is.null(scope_obj@grid[[grid_name]]$counts) ||
      is.null(scope_obj@grid[[grid_name]]$Xz) ||
      is.null(scope_obj@grid[[grid_name]]$W)) {
    stop(sample_id, " frozen analysis source is missing grid counts/Xz/W.")
  }
  if (is.null(scope_obj@stats[[grid_name]]$iDeltaStats) ||
      !paste0(grid_name, "_iDelta") %in% colnames(scope_obj@meta.data)) {
    stop(sample_id, " frozen analysis source is missing I-delta results.")
  }
  raw_membership <- as.character(scope_obj@meta.data[[cluster_col]])
  assigned <- !is.na(raw_membership) & nzchar(raw_membership) & raw_membership != "-1"
  graph_genes <- igraph::V(lee[[cluster_col]])$name
  if (!setequal(graph_genes, meta_genes[assigned])) {
    stop(sample_id, " authoritative graph vertices do not match assigned genes.")
  }
  invisible(list(
    formula_id = meta$formula_id,
    permutation_scheme = meta$permutation_scheme,
    curve_name = curve_name,
    cluster_column = cluster_col,
    genes = length(genes),
    graph_vertices = igraph::vcount(lee[[cluster_col]]),
    graph_edges = igraph::ecount(lee[[cluster_col]])
  ))
}

assert_scope_xenium_identity <- function(scope_obj, input_dir, roi_file, sample_id,
                                         tolerance = 1e-6) {
  if (!requireNamespace("arrow", quietly = TRUE) ||
      !requireNamespace("rhdf5", quietly = TRUE) ||
      !requireNamespace("sp", quietly = TRUE)) {
    stop("arrow, rhdf5, and sp are required for the scope-to-Xenium identity gate.")
  }
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  roi_file <- normalizePath(roi_file, mustWork = TRUE)
  cells_path <- file.path(input_dir, "cells.parquet")
  matrix_path <- file.path(input_dir, "cell_feature_matrix.h5")
  if (!file.exists(cells_path) || !file.exists(matrix_path)) {
    stop(sample_id, " Xenium identity gate requires cells.parquet and cell_feature_matrix.h5.")
  }
  raw_cells <- arrow::read_parquet(
    cells_path,
    col_select = c("cell_id", "x_centroid", "y_centroid"),
    as_data_frame = TRUE
  )
  required_raw <- c("cell_id", "x_centroid", "y_centroid")
  if (!identical(names(raw_cells), required_raw) || anyDuplicated(raw_cells$cell_id)) {
    stop(sample_id, " raw cells.parquet has invalid cell identity columns.")
  }
  centroids <- scope_obj@coord$centroids
  if (!is.data.frame(centroids) ||
      !all(c("cell", "x", "y") %in% names(centroids)) ||
      anyDuplicated(centroids$cell)) {
    stop(sample_id, " frozen scope has invalid cell centroids.")
  }
  matched <- match(centroids$cell, raw_cells$cell_id)
  if (anyNA(matched)) {
    stop(sample_id, " frozen scope contains cells absent from the supplied Xenium data.")
  }
  x_error <- max(abs(as.numeric(centroids$x) -
                     as.numeric(raw_cells$x_centroid[matched])), na.rm = TRUE)
  y_direct_error <- max(abs(as.numeric(centroids$y) -
                            as.numeric(raw_cells$y_centroid[matched])), na.rm = TRUE)
  y_sum <- as.numeric(centroids$y) + as.numeric(raw_cells$y_centroid[matched])
  y_flip_constant <- stats::median(y_sum, na.rm = TRUE)
  y_flip_error <- max(abs(y_sum - y_flip_constant), na.rm = TRUE)
  coordinate_transform <- if (is.finite(y_direct_error) && y_direct_error <= tolerance) {
    "identity"
  } else if (is.finite(y_flip_error) && y_flip_error <= tolerance) {
    "x_identity_y_reflection"
  } else {
    stop(
      sample_id, " frozen scope coordinates do not match supplied Xenium data: ",
      "x_error=", x_error, ", y_direct_error=", y_direct_error,
      ", y_flip_error=", y_flip_error
    )
  }
  if (!is.finite(x_error) || x_error > tolerance) {
    stop(sample_id, " frozen scope x coordinates do not match supplied Xenium data.")
  }
  cell_matrix <- scope_obj@cells$logCPM
  scope_cells <- colnames(cell_matrix)
  scope_genes <- rownames(cell_matrix)
  if (!identical(as.character(scope_cells), as.character(centroids$cell))) {
    stop(sample_id, " frozen scope cell matrix and centroids are not identically ordered.")
  }
  raw_barcodes <- as.character(rhdf5::h5read(matrix_path, "/matrix/barcodes"))
  raw_features <- as.character(rhdf5::h5read(matrix_path, "/matrix/features/name"))
  if (anyDuplicated(raw_barcodes) || anyDuplicated(raw_features) ||
      !setequal(raw_barcodes, as.character(raw_cells$cell_id)) ||
      !all(scope_cells %in% raw_barcodes) || !all(scope_genes %in% raw_features)) {
    stop(sample_id, " frozen scope cells/genes do not match the supplied Xenium matrix.")
  }
  roi <- utils::read.csv(roi_file, comment.char = "#", check.names = FALSE)
  if (!identical(names(roi), c("X", "Y")) || nrow(roi) < 4L ||
      any(!is.finite(as.matrix(roi)))) {
    stop(sample_id, " ROI coordinate file is malformed.")
  }
  inside <- sp::point.in.polygon(
    raw_cells$x_centroid, raw_cells$y_centroid, roi$X, roi$Y
  ) > 0L
  roi_cells <- as.character(raw_cells$cell_id[inside])
  if (!setequal(roi_cells, as.character(centroids$cell))) {
    stop(
      sample_id, " frozen scope cell set does not equal the supplied ROI selection: ",
      "ROI=", length(roi_cells), ", scope=", nrow(centroids)
    )
  }
  invisible(list(
    raw_cells = length(raw_barcodes),
    raw_features = length(raw_features),
    scope_cells = length(scope_cells),
    scope_genes = length(scope_genes),
    roi_selected_cells = length(roi_cells),
    coordinate_transform = coordinate_transform,
    max_abs_x_error = x_error,
    max_abs_y_error = if (identical(coordinate_transform, "identity")) {
      y_direct_error
    } else {
      y_flip_error
    },
    y_reflection_constant = if (identical(coordinate_transform,
                                           "x_identity_y_reflection")) {
      y_flip_constant
    } else {
      NA_real_
    }
  ))
}

assert_scope_pair_table <- function(scope_obj, pairs, sample_id,
                                    numeric_tolerance = 1e-12,
                                    pct_tolerance = 5.01e-4) {
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  required <- c("gene1", "gene2", "L", "r", "Delta", "pct1", "pct2")
  missing <- setdiff(required, names(pairs))
  if (length(missing)) {
    stop(sample_id, " pair/scope gate is missing: ", paste(missing, collapse = ", "))
  }
  lee <- scope_obj@stats$grid30$LeeStats_Xz$L
  pearson <- scope_obj@cells$.pearson_cor
  if (is.null(lee) || is.null(pearson)) {
    stop(sample_id, " pair/scope gate requires Lee L and cell Pearson matrices.")
  }
  i <- match(pairs$gene1, rownames(lee))
  j <- match(pairs$gene2, colnames(lee))
  ir <- match(pairs$gene1, rownames(pearson))
  jr <- match(pairs$gene2, colnames(pearson))
  if (anyNA(c(i, j, ir, jr))) {
    stop(sample_id, " pair table contains genes absent from its frozen scope.")
  }
  observed_L <- as.numeric(lee[cbind(i, j)])
  observed_r <- as.numeric(pearson[cbind(ir, jr)])
  L_error <- max(abs(as.numeric(pairs$L) - observed_L), na.rm = TRUE)
  r_error <- max(abs(as.numeric(pairs$r) - observed_r), na.rm = TRUE)
  delta_error <- max(abs(as.numeric(pairs$Delta) -
                         (as.numeric(pairs$L) - as.numeric(pairs$r))), na.rm = TRUE)
  if (any(!is.finite(c(L_error, r_error, delta_error))) ||
      any(c(L_error, r_error, delta_error) > numeric_tolerance)) {
    stop(
      sample_id, " pair table does not match its frozen scope: L=", L_error,
      ", r=", r_error, ", Delta=", delta_error
    )
  }
  pair_genes <- unique(c(as.character(pairs$gene1), as.character(pairs$gene2)))
  grid_layer <- scope_obj@grid$grid30
  counts <- grid_layer$counts
  if (!is.data.frame(counts) ||
      !all(c("gene", "grid_id") %in% names(counts)) ||
      is.null(grid_layer$grid_info) || nrow(grid_layer$grid_info) < 1L) {
    stop(sample_id, " pair/scope gate requires long grid counts and grid_info.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required for the pair prevalence gate.")
  }
  counts <- data.table::as.data.table(counts)
  gene_coverage <- counts[
    gene %in% pair_genes,
    .(detected_grids = data.table::uniqueN(grid_id)),
    by = gene
  ]
  pct_lookup <- stats::setNames(rep(0, length(pair_genes)), pair_genes)
  pct_lookup[gene_coverage$gene] <- round(
    100 * gene_coverage$detected_grids / nrow(grid_layer$grid_info), 3
  )
  pct1_error <- max(abs(as.numeric(pairs$pct1) - pct_lookup[pairs$gene1]), na.rm = TRUE)
  pct2_error <- max(abs(as.numeric(pairs$pct2) - pct_lookup[pairs$gene2]), na.rm = TRUE)
  if (any(!is.finite(c(pct1_error, pct2_error))) ||
      any(c(pct1_error, pct2_error) > pct_tolerance)) {
    stop(sample_id, " pair-table prevalence does not match its frozen scope: pct1=",
         pct1_error, ", pct2=", pct2_error)
  }
  invisible(list(
    pairs = nrow(pairs),
    pair_genes = length(pair_genes),
    max_abs_L_error = L_error,
    max_abs_r_error = r_error,
    max_abs_Delta_error = delta_error,
    max_abs_pct1_error = pct1_error,
    max_abs_pct2_error = pct2_error
  ))
}

generator_call_name <- function(call) {
  if (!is.call(call) || !length(call)) return(NA_character_)
  head <- call[[1L]]
  if (is.symbol(head)) return(as.character(head))
  if (is.call(head) && length(head) == 3L &&
      as.character(head[[1L]]) %in% c("::", ":::")) {
    return(as.character(head[[3L]]))
  }
  NA_character_
}

collect_generator_calls <- function(node, target) {
  found <- list()
  if (is.call(node)) {
    if (identical(generator_call_name(node), target)) found <- list(node)
    for (index in seq_along(node)[-1L]) {
      if (is.symbol(node[[index]]) && !nzchar(as.character(node[[index]]))) next
      found <- c(found, collect_generator_calls(node[[index]], target))
    }
  } else if (is.expression(node) || is.pairlist(node) ||
             (is.list(node) && !is.object(node))) {
    for (index in seq_along(node)) {
      if (is.symbol(node[[index]]) && !nzchar(as.character(node[[index]]))) next
      found <- c(found, collect_generator_calls(node[[index]], target))
    }
  }
  found
}

generator_argument <- function(call, name) {
  args <- as.list(call)[-1L]
  hit <- which(names(args) == name)
  if (length(hit) != 1L) return(NA_character_)
  gsub("[[:space:]]+", "", paste(deparse(args[[hit]]), collapse = ""))
}

assert_generator_provenance <- function(sample_id, generator_path,
                                        top_pair_rows) {
  generator_path <- normalizePath(generator_path, mustWork = TRUE)
  parsed <- parse(generator_path, keep.source = FALSE)
  top_calls <- collect_generator_calls(parsed, "getTopLvsR")
  if (length(top_calls) != 1L) {
    stop(sample_id, " generator must contain exactly one getTopLvsR call.")
  }
  top_call <- top_calls[[1L]]
  common_ok <- identical(generator_argument(top_call, "use_blocks"), "FALSE") &&
    identical(generator_argument(top_call, "direction"), '"largest"') &&
    identical(generator_argument(top_call, "pval_mode"), '"uniform"')
  if (!isTRUE(common_ok)) {
    stop(sample_id, " generator does not declare the frozen shuffle contract.")
  }

  assignments <- collect_generator_calls(parsed, "<-")
  seed_assignment_ok <- any(vapply(assignments, function(call) {
    length(call) == 3L && identical(as.character(call[[2L]]), "SEED") &&
      identical(gsub("[[:space:]]+", "", paste(deparse(call[[3L]]),
                                               collapse = "")), "1L")
  }, logical(1L)))
  set_seed_calls <- collect_generator_calls(parsed, "set.seed")
  rng_calls <- collect_generator_calls(parsed, "RNGkind")

  if (sample_id %in% c("P1", "P2", "P5")) {
    top_n <- 40000L
    seed_ok <- seed_assignment_ok && any(vapply(set_seed_calls, function(call) {
      length(call) >= 2L && identical(as.character(call[[2L]]), "SEED")
    }, logical(1L)))
    rng_ok <- any(vapply(rng_calls, function(call) {
      length(call) >= 2L &&
        identical(as.character(call[[2L]]), "L'Ecuyer-CMRG")
    }, logical(1L)))
    contract_ok <- generator_argument(top_call, "top_n") %in% c("40000", "40000L") &&
      identical(generator_argument(top_call, "curve_layer"), "curve_name") &&
      is.na(generator_argument(top_call, "p_adj_mode")) &&
      is.na(generator_argument(top_call, "perms"))
  } else if (identical(sample_id, "LN")) {
    top_n <- 100000L
    seed_ok <- any(vapply(set_seed_calls, function(call) {
      length(call) >= 2L &&
        gsub("[[:space:]]+", "", paste(deparse(call[[2L]]), collapse = "")) == "1L"
    }, logical(1L)))
    rng_ok <- any(vapply(rng_calls, function(call) {
      length(call) >= 2L &&
        identical(as.character(call[[2L]]), "L'Ecuyer-CMRG")
    }, logical(1L)))
    contract_ok <- generator_argument(top_call, "top_n") %in%
      c("100000", "100000L") &&
      generator_argument(top_call, "perms") %in% c("1000", "1000L") &&
      identical(generator_argument(top_call, "p_adj_mode"), '"BH_universe"') &&
      identical(generator_argument(top_call, "curve_layer"),
                '"LR_curve_30_shuffle"')
  } else {
    stop("Unsupported correction-render sample: ", sample_id)
  }
  if (!isTRUE(seed_ok) || !isTRUE(rng_ok) || !isTRUE(contract_ok) ||
      top_pair_rows >= top_n) {
    stop(sample_id, " generator failed its seed/RNG/complete-universe contract.")
  }
  invisible(list(
    path = generator_path,
    sha256 = sha256_file(generator_path),
    top_n = top_n,
    total_universe = as.integer(top_pair_rows),
    seed = 1L,
    rng = "L'Ecuyer-CMRG",
    use_blocks = FALSE
  ))
}

assert_analysis_provenance <- function(sample_id, analysis_manifest_path,
                                       top_pair_rows,
                                       delta_manifest_path = NULL,
                                       generator_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required for analysis-provenance gates.")
  }
  analysis <- jsonlite::read_json(analysis_manifest_path, simplifyVector = TRUE)
  expected_cluster <- if (identical(sample_id, "P5")) {
    "shuffle_q95_res0.1_grid30"
  } else if (identical(sample_id, "LN")) {
    "shuffle_q99.9_res0.1_grid30"
  } else {
    stop("Unsupported correction-render sample: ", sample_id)
  }
  common_ok <- identical(as.character(analysis$dataset), sample_id) &&
    identical(as.character(analysis$package_version), "1.0.2") &&
    identical(as.integer(analysis$permutations), 1000L) &&
    identical(as.integer(analysis$seed), 1L) &&
    is.character(analysis$permutation) &&
    grepl("all-grid joint shuffle", analysis$permutation) &&
    is.character(analysis$main_FDR) && grepl("BH", analysis$main_FDR) &&
    is.character(analysis$toplvsr) &&
    grepl("use_blocks=FALSE", analysis$toplvsr, fixed = TRUE) &&
    grepl("Delta FDR", analysis$toplvsr, fixed = TRUE) &&
    identical(as.character(analysis$cluster_column), expected_cluster)
  if (!isTRUE(common_ok)) {
    stop(sample_id, " authoritative analysis manifest failed its frozen contract.")
  }
  generator_gate <- assert_generator_provenance(
    sample_id, generator_path, top_pair_rows
  )
  summary_rows <- as.integer(analysis$summary$candidates_after_curve[[1L]])
  if (!identical(summary_rows, as.integer(if (identical(sample_id, "P5")) {
    top_pair_rows
  } else {
    40000L
  }))) {
    stop(sample_id, " authoritative analysis manifest has an unexpected candidate count.")
  }

  if (identical(sample_id, "P5")) {
    if (!is.null(delta_manifest_path)) {
      stop("P5 must not use the LN complete-Delta manifest.")
    }
    return(invisible(list(
      permutations = 1000L,
      seed = 1L,
      permutation_scheme = "global_joint_shuffle",
      use_blocks = FALSE,
      delta_adjustment = "BH over the complete eligible universe",
      total_universe = as.integer(top_pair_rows),
      selected_pairs = as.integer(top_pair_rows),
      top_pair_fdr_semantics = "Delta permutation FDR",
      generator = generator_gate
    )))
  }

  if (is.null(delta_manifest_path) || !file.exists(delta_manifest_path)) {
    stop("LN requires the complete-Delta provenance manifest.")
  }
  delta <- jsonlite::read_json(delta_manifest_path, simplifyVector = TRUE)
  permutation <- delta$delta_permutation
  delta_ok <- identical(as.character(delta$sample_id), "LN") &&
    identical(as.character(delta$package_version), "1.0.2") &&
    identical(as.character(permutation$schema), "geneSCOPE_delta_permutation_v1") &&
    isTRUE(permutation$performed) &&
    identical(as.integer(permutation$permutations), 1000L) &&
    identical(permutation$use_blocks, FALSE) &&
    identical(as.character(permutation$pval_mode), "uniform") &&
    identical(as.character(permutation$p_adj_mode), "BH_universe") &&
    identical(as.integer(permutation$total_universe), as.integer(top_pair_rows)) &&
    identical(as.integer(permutation$selected_pairs), as.integer(top_pair_rows))
  if (!isTRUE(delta_ok)) {
    stop("LN complete-Delta manifest failed its frozen contract.")
  }
  invisible(list(
    permutations = 1000L,
    seed = 1L,
    permutation_scheme = "global_joint_shuffle",
    use_blocks = FALSE,
    delta_adjustment = "BH_universe",
    total_universe = as.integer(permutation$total_universe),
    selected_pairs = as.integer(permutation$selected_pairs),
    top_pair_fdr_semantics = "Delta permutation FDR",
    generator = generator_gate
  ))
}

read_authoritative_top_pairs <- function(path, sample_id, expected_rows) {
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("gene1", "gene2", "L", "r", "pct1", "pct2", "Delta")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(sample_id, " authoritative pair table is missing: ", paste(missing, collapse = ", "))
  }
  q_col <- intersect(c("q_Delta", "delta_fdr", "fdr"), names(x))
  if (!length(q_col)) stop(sample_id, " authoritative pair table has no adjusted P column.")
  if (nrow(x) != as.integer(expected_rows)) {
    stop(sample_id, " authoritative pair-table row count changed: observed=", nrow(x),
         "; expected=", expected_rows)
  }
  key <- paste(pmin(x$gene1, x$gene2), pmax(x$gene1, x$gene2), sep = "--")
  if (anyNA(x$gene1) || anyNA(x$gene2) || any(x$gene1 == x$gene2) ||
      anyDuplicated(key)) {
    stop(sample_id, " authoritative pair table contains invalid or duplicate pairs.")
  }
  delta_error <- max(abs(as.numeric(x$Delta) -
                         (as.numeric(x$L) - as.numeric(x$r))), na.rm = TRUE)
  if (!is.finite(delta_error) || delta_error > 1e-12) {
    stop(sample_id, " authoritative pair table failed Delta=L-r: ", delta_error)
  }
  q <- as.numeric(x[[q_col[[1L]]]])
  if (any(!is.finite(q)) || any(q < 0 | q > 1)) {
    stop(sample_id, " authoritative pair table has invalid adjusted P values.")
  }
  x
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

audit_p5_dendrogram_path <- function(dnet_obj, genes, raw_membership,
                                     display_mapping, reference_path,
                                     output_path) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required for the P5 dendrogram-path freeze gate.")
  }
  graph <- if (is.list(dnet_obj)) dnet_obj$graph else NULL
  if (is.null(graph) || !inherits(graph, "igraph") ||
      !igraph::is_connected(graph) || !igraph::is_tree(graph)) {
    stop("P5 dendrogram-path freeze gate requires one connected igraph tree.")
  }

  genes <- as.character(genes)
  raw_membership <- as.character(raw_membership)
  if (length(genes) != length(raw_membership) || anyNA(genes) ||
      any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("Invalid gene membership supplied to the P5 dendrogram-path gate.")
  }
  gene_module <- stats::setNames(raw_membership, genes)
  endpoints <- igraph::as_data_frame(graph, what = "edges")
  endpoints$module_from <- unname(gene_module[endpoints$from])
  endpoints$module_to <- unname(gene_module[endpoints$to])
  cross <- endpoints[
    !is.na(endpoints$module_from) & nzchar(endpoints$module_from) &
      endpoints$module_from != "-1" &
      !is.na(endpoints$module_to) & nzchar(endpoints$module_to) &
      endpoints$module_to != "-1" &
      endpoints$module_from != endpoints$module_to,
    c("module_from", "module_to"), drop = FALSE
  ]
  if (!nrow(cross)) stop("P5 dendrogram has no between-module edges.")
  cross$key <- apply(cross, 1L, function(z) paste(sort(z), collapse = "--"))
  cross <- cross[!duplicated(cross$key), , drop = FALSE]
  module_edges <- do.call(rbind, strsplit(cross$key, "--", fixed = TRUE))
  module_graph <- igraph::graph_from_data_frame(
    data.frame(from = module_edges[, 1L], to = module_edges[, 2L]),
    directed = FALSE
  )

  start_module <- unname(gene_module[["C3"]])
  end_module <- unname(gene_module[["GPX2"]])
  if (is.null(start_module) || is.null(end_module) ||
      anyNA(c(start_module, end_module))) {
    stop("C3 or GPX2 is absent from the P5 membership.")
  }
  module_path <- names(igraph::shortest_paths(
    module_graph, from = start_module, to = end_module
  )$vpath[[1L]])
  if (!length(module_path)) stop("No C3-to-GPX2 module path was found.")

  path_result <- getDendroWalkPaths(
    dnet_obj, gene = "C3", gene2 = "GPX2", cutoff = 20L,
    max_paths = 50000L, verbose = FALSE
  )
  endpoint_paths <- Filter(function(path) {
    if (length(path) < 2L) return(FALSE)
    ends <- c(path[[1L]], path[[length(path)]])
    identical(ends, c("C3", "GPX2")) || identical(ends, c("GPX2", "C3"))
  }, path_result$paths)
  if (length(endpoint_paths) != 1L) {
    stop("Expected one C3--GPX2 dendrogram path; found ", length(endpoint_paths), ".")
  }
  gene_path <- endpoint_paths[[1L]]
  if (!identical(gene_path[[1L]], "C3")) gene_path <- rev(gene_path)

  display_lookup <- stats::setNames(
    as.character(display_mapping$display_module),
    as.character(display_mapping$current_module)
  )
  display_path <- unname(display_lookup[module_path])
  if (anyNA(display_path)) stop("P5 dendrogram path contains an unmapped module.")
  raw_path_string <- paste(module_path, collapse = "->")
  display_path_string <- paste(display_path, collapse = "->")
  internal_display <- if (length(display_path) > 2L) {
    display_path[seq.int(2L, length(display_path) - 1L)]
  } else {
    character()
  }
  broad_module <- display_path[[length(display_path)]]
  consecutive <- if (length(display_path) > 1L) {
    paste(display_path[-length(display_path)], display_path[-1L], sep = "--")
  } else {
    character()
  }
  direct_stem_to_broad <- any(consecutive %in% c(
    paste("3", broad_module, sep = "--"),
    paste(broad_module, "3", sep = "--")
  ))
  observed <- data.frame(
    gene_query = "C3--GPX2",
    gene_path = paste(gene_path, collapse = "->"),
    raw_module_path = raw_path_string,
    display_module_path = display_path_string,
    stem_positioned_between = "3" %in% internal_display,
    direct_stem_to_broad_adjacency = direct_stem_to_broad,
    enumerated_paths_cutoff20 = length(path_result$paths),
    stringsAsFactors = FALSE
  )

  reference <- utils::read.delim(reference_path, stringsAsFactors = FALSE,
                                 check.names = FALSE)
  if (nrow(reference) != 1L ||
      !identical(observed$gene_query[[1L]], reference$gene_query[[1L]]) ||
      !identical(observed$raw_module_path[[1L]], reference$raw_module_path[[1L]]) ||
      !identical(observed$display_module_path[[1L]], reference$display_module_path[[1L]]) ||
      !identical(observed$stem_positioned_between[[1L]],
                 as.logical(reference$stem_positioned_between[[1L]])) ||
      !identical(observed$direct_stem_to_broad_adjacency[[1L]],
                 as.logical(reference$direct_stem_to_broad_adjacency[[1L]]))) {
    stop(
      "P5 dendrogram-path freeze gate failed: raw=", raw_path_string,
      "; display=", display_path_string,
      "; gene_path=", observed$gene_path[[1L]]
    )
  }
  utils::write.table(observed, output_path, sep = "\t", row.names = FALSE,
                     quote = FALSE)
  invisible(observed)
}

assert_required_figure_outputs <- function(output_root, sample_id,
                                           required_relative,
                                           expected_png_count) {
  output_root <- normalizePath(output_root, mustWork = TRUE)
  required_relative <- sub("^\\./+", "", unique(as.character(required_relative)))
  if (any(!nzchar(required_relative)) ||
      any(grepl("^/|(^|/)\\.\\.(/|$)", required_relative))) {
    stop(sample_id, " figure bundle has an unsafe required-output path.")
  }
  required <- file.path(output_root, required_relative)
  missing <- required[!file.exists(required) | dir.exists(required)]
  if (length(missing)) {
    stop(sample_id, " figure bundle is missing required outputs: ",
         paste(basename(missing), collapse = ", "))
  }
  sizes <- file.info(required)$size
  if (anyNA(sizes) || any(sizes <= 0)) {
    stop(sample_id, " figure bundle contains an empty required output.")
  }

  files <- list.files(
    output_root, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE
  )
  files <- files[file.exists(files) & !dir.exists(files)]
  relative <- substring(files, nchar(output_root) + 2L)
  keep <- !grepl("^\\.geneSCOPE-v1\\.0\\.2-library/", relative)
  files <- files[keep]
  relative <- relative[keep]
  png_files <- files[grepl("\\.png$", relative, ignore.case = TRUE)]
  if (length(png_files) != as.integer(expected_png_count)) {
    stop(sample_id, " figure bundle PNG count changed: observed=",
         length(png_files), ", expected=", expected_png_count)
  }
  png_ok <- vapply(png_files, function(path) {
    size <- file.info(path)$size
    if (is.na(size) || size < 1000) return(FALSE)
    signature <- readBin(path, what = "raw", n = 8L)
    identical(as.integer(signature), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))
  }, logical(1L))
  if (!all(png_ok)) {
    stop(sample_id, " figure bundle contains an invalid PNG: ",
         paste(basename(png_files[!png_ok]), collapse = ", "))
  }
  invisible(list(
    required_files = required_relative,
    expected_png_count = as.integer(expected_png_count),
    observed_png_count = length(png_files)
  ))
}

write_freeze_output_manifest <- function(output_root, sample_id, freeze_source,
                                         formula_gate, membership_gate,
                                         input_dir, roi_file,
                                         display_mapping_path, workflow_path,
                                         analysis_sources, parameters, output_gate) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to write the figure freeze manifest.")
  }
  output_root <- normalizePath(output_root, mustWork = TRUE)
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  roi_file <- normalizePath(roi_file, mustWork = TRUE)
  display_mapping_path <- normalizePath(display_mapping_path, mustWork = TRUE)
  workflow_path <- normalizePath(workflow_path, mustWork = TRUE)
  paper_commit <- require_paper_commit_metadata()
  if (!is.list(analysis_sources) ||
      !all(c("scope", "top_pairs", "analysis_manifest") %in%
           names(analysis_sources)) || anyDuplicated(names(analysis_sources))) {
    stop("Figure manifest requires named scope, top_pairs, and analysis_manifest sources.")
  }
  analysis_sources <- lapply(analysis_sources, function(source) {
    required <- c("path", "root_token", "relative_path", "sha256")
    if (!is.list(source) || !all(required %in% names(source))) {
      stop("Malformed frozen analysis-source provenance.")
    }
    source$path <- normalizePath(source$path, mustWork = TRUE)
    if (!identical(sha256_file(source$path), tolower(source$sha256))) {
      stop("Frozen analysis source changed while rendering: ", source$path)
    }
    source
  })
  raw_gate <- assert_reference_raw_inputs(
    dirname(workflow_path), sample_id, input_dir, roi_file
  )
  raw_inputs <- unlist(raw_gate$paths, use.names = TRUE)
  input_hashes <- raw_gate$sha256
  manifest_path <- file.path(output_root, paste0(sample_id, "_figure_manifest.json"))
  inventory <- bundle_file_inventory(output_root, manifest_path)
  files <- inventory$files
  relative <- inventory$relative
  hashes <- if (length(files)) {
    stats::setNames(as.list(vapply(files, sha256_file, character(1L))), relative)
  } else {
    list()
  }
  loaded_namespaces <- sort(unique(loadedNamespaces()), method = "radix")
  direct_packages <- c(
    "geneSCOPE", "arrow", "rhdf5", "sp", "data.table", "Matrix",
    "jsonlite", "png", "ggplot2", "ggraph", "igraph", "scales"
  )
  if (identical(sample_id, "P5")) {
    direct_packages <- c(direct_packages, "ComplexHeatmap", "circlize")
  }
  runtime_packages <- sort(
    unique(c(loaded_namespaces, direct_packages)), method = "radix"
  )
  runtime_versions <- vapply(runtime_packages, function(package) {
    suppressWarnings(tryCatch(
      as.character(utils::packageVersion(package)),
      error = function(...) NA_character_
    ))
  }, character(1L))
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
    output_gate = output_gate,
    parameters = parameters,
    runtime = list(
      R_version = R.version.string,
      platform = R.version$platform,
      R_executable = file.path(R.home("bin"), "R"),
      package_inventory_policy =
        "loaded_namespaces_plus_direct_workflow_and_verifier_dependencies",
      loaded_namespaces = loaded_namespaces,
      packages = as.list(runtime_versions)
    ),
    inputs = list(
      xenium_outs = input_dir,
      roi_file = roi_file,
      sha256 = input_hashes,
      display_mapping = list(
        path = display_mapping_path,
        sha256 = sha256_file(display_mapping_path)
      ),
      analysis_sources = analysis_sources
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

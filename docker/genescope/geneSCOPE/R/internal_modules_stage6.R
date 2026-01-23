#' Convert Flag
#' @description
#' Internal helper for `.convert_flag`.
#' @param value Parameter value.
#' @return Return value used internally.
#' @keywords internal
.convert_flag <- function(value) {
  if (identical(value, TRUE)) return(TRUE)
  if (identical(value, FALSE)) return(FALSE)
  if (is.character(value)) {
    if (tolower(value[1]) %in% c("true", "t", "1")) return(TRUE)
    if (tolower(value[1]) %in% c("false", "f", "0")) return(FALSE)
  }
  if (is.numeric(value)) {
    return(!is.na(value) && value != 0)
  }
  FALSE
}

#' Select Candidate
#' @description
#' Internal helper for `.select_candidate`.
#' @param candidates Parameter value.
#' @return Return value used internally.
#' @keywords internal
.select_candidate <- function(candidates) {
  if (!length(candidates)) return(NULL)
  candidates[order(
    !grepl("\\.gz$", candidates, ignore.case = TRUE),
    candidates
  )][1]
}

#' Contract Trace Enabled
#' @description
#' Internal helper for `.contract_trace_enabled`.
#' @return Return value used internally.
#' @keywords internal
.contract_trace_enabled <- function() {
  flag <- Sys.getenv("GENESCOPE_CONTRACT_TRACE", unset = "")
  if (!nzchar(flag)) {
    return(FALSE)
  }
  flag <- tolower(flag)
  flag %in% c("1", "true", "yes", "on")
}

#' Trace Contract Field
#' @description
#' Internal helper for `.trace_contract_field`.
#' @param container Parameter value.
#' @param field Parameter value.
#' @param where Parameter value.
#' @param preview Parameter value.
#' @return Return value used internally.
#' @keywords internal
.trace_contract_field <- function(container,
                                  field,
                                  where = NA_character_,
                                  preview = 3L) {
  if (!.contract_trace_enabled()) {
    return(invisible(NULL))
  }
  loc <- if (!is.na(where) && nzchar(where)) where else "unknown_context"
  value <- NULL
  has_field <- FALSE
  try({
    if (!is.null(container) && !is.null(field) && nzchar(field) && !is.environment(container)) {
      if (is.list(container) || is.environment(container)) {
        has_field <- field %in% names(container)
      }
    }
  }, silent = TRUE)
  if (is.list(container) || is.environment(container)) {
    value <- tryCatch(container[[field]], error = function(e) NULL)
  } else if (!is.null(container) && !is.na(field)) {
    value <- tryCatch(container[[field]], error = function(e) NULL)
  }
  type_info <- if (is.null(value)) "NULL" else paste(class(value), collapse = "/")
  len_info <- if (is.null(value)) 0L else {
    if (is.environment(value)) NA_integer_ else length(value)
  }
  preview_vec <- NULL
  if (is.atomic(value) || is.list(value)) {
    vec <- suppressWarnings(as.character(head(value, preview)))
    preview_vec <- paste(vec, collapse = ", ")
  } else if (is.null(value)) {
    preview_vec <- "<missing>"
  } else {
    preview_vec <- "<non-atomic>"
  }
  msg <- sprintf("[contract_trace] %s$%s exists=%s type=%s length=%s preview=%s",
                 loc,
                 field,
                 ifelse(isTRUE(has_field), "TRUE", "FALSE"),
                 type_info,
                 ifelse(is.na(len_info), "NA", as.character(len_info)),
                 preview_vec)
  message(msg)
  invisible(value)
}

#' Scope Input Validate Inputs
#' @description
#' Internal helper for `.scope_input_validate_inputs`.
#' @param input_dir Filesystem path.
#' @param grid_sizes Parameter value.
#' @param segmentation_strategy Parameter value.
#' @param allow_flatfile_generation Parameter value.
#' @param derive_cells_from_polygons Parameter value.
#' @param data_type Parameter value.
#' @param pixel_size_um Parameter value.
#' @param dataset_id Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.scope_input_validate_inputs <- function(input_dir,
                                         grid_sizes,
                                         segmentation_strategy,
                                         allow_flatfile_generation,
                                         derive_cells_from_polygons,
                                         data_type,
                                         pixel_size_um,
                                         dataset_id,
                                         verbose) {
  if (!dir.exists(input_dir)) {
    stop("input_dir not found: ", input_dir)
  }
  if (!is.numeric(grid_sizes) || length(grid_sizes) == 0L || any(!is.finite(grid_sizes))) {
    stop("grid_sizes must be a non-empty numeric vector.")
  }
  if (!is.numeric(pixel_size_um) || length(pixel_size_um) != 1L || !is.finite(pixel_size_um) || pixel_size_um <= 0) {
    stop("pixel_size_um must be a single positive finite number.")
  }
  res <- list(
    input_dir = input_dir,
    grid_sizes = as.numeric(grid_sizes),
    segmentation_strategy = match.arg(segmentation_strategy, c("cell", "nucleus", "both", "none")),
    allow_flatfile_generation = .convert_flag(allow_flatfile_generation),
    derive_cells_from_polygons = .convert_flag(derive_cells_from_polygons),
    data_type = if (!is.null(data_type)) as.character(data_type)[1] else "cosmx",
    pixel_size_um = as.numeric(pixel_size_um),
    dataset_id = if (!is.null(dataset_id) && nzchar(as.character(dataset_id))) as.character(dataset_id)[1] else NULL,
    verbose = .convert_flag(verbose)
  )
  .trace_contract_field(res, "grid_sizes", "scope_input.validate_inputs")
  res
}

#' Scope Input Spec Build
#' @description
#' Internal helper for `.scope_input_spec_build`.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.scope_input_spec_build <- function(validated_inputs) {
  input_dir <- validated_inputs$input_dir
  ff_dir <- file.path(input_dir, "flatFiles")
  transcripts_path <- file.path(input_dir, "transcripts.parquet")
  cells_path <- file.path(input_dir, "cells.parquet")
  seg_candidates <- c(
    file.path(input_dir, "cell_boundaries.parquet"),
    file.path(input_dir, "segmentation_boundaries.parquet"),
    file.path(input_dir, "nucleus_boundaries.parquet")
  )
  tx_candidates <- unique(c(
    Sys.glob(file.path(ff_dir, "*/*_tx_file.csv.gz")),
    Sys.glob(file.path(ff_dir, "*/*_tx_file.csv"))
  ))
  polygon_candidates <- unique(c(
    Sys.glob(file.path(ff_dir, "*/*-polygons.csv.gz")),
    Sys.glob(file.path(ff_dir, "*/*_polygons.csv.gz"))
  ))
  flatfile_dirs <- character(0)
  if (dir.exists(ff_dir)) {
    flatfile_dirs <- sort(list.dirs(ff_dir, recursive = FALSE, full.names = TRUE))
  }
  existing_segments <- seg_candidates[file.exists(seg_candidates)]
  list(
    module_name = ".scope_input_module",
    input_dir = input_dir,
    grid_length = validated_inputs$grid_sizes,
    segmentation_strategy = validated_inputs$segmentation_strategy,
    allow_flatfile_generation = validated_inputs$allow_flatfile_generation,
    derive_cells_from_polygons = validated_inputs$derive_cells_from_polygons,
    data_type = validated_inputs$data_type,
    pixel_size_um = validated_inputs$pixel_size_um,
    dataset_id = validated_inputs$dataset_id,
    transcripts_present = file.exists(transcripts_path),
    cells_present = file.exists(cells_path),
    flatfile_dirs = flatfile_dirs,
    tx_candidates = sort(tx_candidates),
    polygon_candidates = sort(polygon_candidates),
    segmentation_candidates = sort(existing_segments),
    ff_has_candidates = length(tx_candidates) > 0 || length(polygon_candidates) > 0
  )
}

#' Scope Input Materialize
#' @description
#' Internal helper for `.scope_input_materialize`.
#' @param spec Parameter value.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.scope_input_materialize <- function(spec, validated_inputs) {
  ff_dir <- file.path(validated_inputs$input_dir, "flatFiles")
  transcripts_path <- file.path(validated_inputs$input_dir, "transcripts.parquet")
  cells_path <- file.path(validated_inputs$input_dir, "cells.parquet")
  seg_type <- validated_inputs$segmentation_strategy
  work_dir <- validated_inputs$input_dir
  need_temp <- FALSE
  generated_transcripts <- FALSE
  generated_cells <- FALSE
  polygon_path_used <- NULL
  tx_file_used <- NULL
  transcripts_source <- "existing"
  cells_source <- "existing"

  ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg, " is required for .scope_input_module")
    }
  }
  ensure_package("data.table")

  if (!file.exists(transcripts_path)) {
    if (!validated_inputs$allow_flatfile_generation) {
      stop("transcripts.parquet not found and allow_flatfile_generation is FALSE.")
    }
    if (validated_inputs$verbose) {
      .log_info("createSCOPE", "S04", "transcripts.parquet not found; generating from flatFiles ...", validated_inputs$verbose)
    }
    selected_tx <- .select_candidate(spec$tx_candidates)
    if (is.null(selected_tx)) {
      stop("No transcript candidate found under flatFiles/")
    }
    txt <- fread(selected_tx)
    cn <- tolower(names(txt))
    names(txt) <- cn
    headerless <- all(grepl("^v[0-9]+$", cn))
    if (!headerless && all(c("x_global_px", "y_global_px", "target") %in% cn)) {
      xg <- txt[["x_global_px"]]
      yg <- txt[["y_global_px"]]
      tgt <- txt[["target"]]
      qv <- if ("qv" %in% cn) txt[["qv"]] else Inf
      nd <- if ("nucleus_distance" %in% cn) txt[["nucleus_distance"]] else 0
    } else {
      if (ncol(txt) < 9) stop("headerless tx_file has insufficient columns (<9)")
      if (validated_inputs$verbose) {
        .log_info("createSCOPE", "S04", "tx_file headerless - using positional mapping col6/col7/col9 (and col8 as nucleus_distance)", validated_inputs$verbose)
      }
      xg <- txt[[6]]
      yg <- txt[[7]]
      tgt <- txt[[9]]
      qv <- if (ncol(txt) >= 10 && is.numeric(txt[[10]])) txt[[10]] else Inf
      nd <- if (ncol(txt) >= 8 && is.numeric(txt[[8]])) txt[[8]] else 0
    }
    tx_out <- data.table(
      x_location = as.numeric(xg) * validated_inputs$pixel_size_um,
      y_location = as.numeric(yg) * validated_inputs$pixel_size_um,
      feature_name = as.character(tgt),
      qv = as.numeric(qv),
      nucleus_distance = as.numeric(nd)
    )
    temp_dir <- file.path(tempdir(), paste0("cosmx_scope_", as.integer(Sys.time())))
    if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
    work_dir <- temp_dir
    need_temp <- TRUE
    transcripts_path <- file.path(work_dir, "transcripts.parquet")
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("arrow is required to emit transcripts.parquet")
    }
    if (validated_inputs$verbose) {
      .log_info("createSCOPE", "S04", paste0("Writing ", transcripts_path), validated_inputs$verbose)
    }
    arrow::write_parquet(tx_out, transcripts_path)
    generated_transcripts <- TRUE
    transcripts_source <- "generated"
    tx_file_used <- basename(selected_tx)
    if (is.null(validated_inputs$dataset_id)) {
      derived <- sub("_tx_file\\.csv(\\.gz)?$", "", basename(selected_tx))
      validated_inputs$dataset_id <- if (nzchar(derived)) derived else NULL
    }
  } else if (validated_inputs$verbose) {
    .log_info("createSCOPE", "S04", "Using existing transcripts.parquet", validated_inputs$verbose)
  }

  if (!file.exists(cells_path) &&
      validated_inputs$derive_cells_from_polygons &&
      validated_inputs$allow_flatfile_generation) {
    polygon_candidates <- spec$polygon_candidates
    selected_poly <- .select_candidate(polygon_candidates)
    if (!is.null(selected_poly)) {
      if (validated_inputs$verbose) {
        .log_info(
          "createSCOPE",
          "S04",
          paste0("cells.parquet not found; deriving centroids from polygons: ", basename(selected_poly)),
          validated_inputs$verbose
        )
      }
      pol <- fread(selected_poly)
      cn <- tolower(names(pol))
      names(pol) <- cn
      headerless <- all(grepl("^v[0-9]+$", cn))
      if (headerless) {
        if (ncol(pol) < 7) stop("headerless polygons has insufficient columns (<7)")
        fov <- pol[[1]]
        cid <- pol[[2]]
        xg <- pol[[6]]
        yg <- pol[[7]]
      } else {
        fov_val <- if ("fov" %in% cn) pol[["fov"]] else pol[[1]]
        cid_val <- if ("cellid" %in% cn) pol[["cellid"]] else if ("cell_id" %in% cn) pol[["cell_id"]] else pol[[2]]
        xg <- if ("x_global_px" %in% cn) pol[["x_global_px"]] else pol[[which.max(cn == "x")]]
        yg <- if ("y_global_px" %in% cn) pol[["y_global_px"]] else pol[[which.max(cn == "y")]]
        fov <- fov_val
        cid <- cid_val
      }
      dt <- data.table(
        fov = as.integer(fov),
        cell_id = as.character(cid),
        x = as.numeric(xg) * validated_inputs$pixel_size_um,
        y = as.numeric(yg) * validated_inputs$pixel_size_um
      )
      cent <- dt[, .(x_centroid = mean(x, na.rm = TRUE),
                     y_centroid = mean(y, na.rm = TRUE)),
                  by = .(cell_id, fov)]
      if (!need_temp) {
        temp_dir <- file.path(tempdir(), paste0("cosmx_scope_", as.integer(Sys.time())))
        if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)
        work_dir <- temp_dir
        need_temp <- TRUE
        if (!generated_transcripts) {
          transcripts_path <- file.path(work_dir, "transcripts.parquet")
        }
        cells_path <- file.path(work_dir, "cells.parquet")
      }
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("arrow is required to emit cells.parquet")
      }
      if (validated_inputs$verbose) {
        .log_info("createSCOPE", "S04", paste0("Writing ", cells_path), validated_inputs$verbose)
      }
      arrow::write_parquet(cent, file.path(work_dir, "cells.parquet"))
      cells_path <- file.path(work_dir, "cells.parquet")
      generated_cells <- TRUE
      cells_source <- "polygon"
      polygon_path_used <- basename(selected_poly)
    }
  }

  have_cell_seg <- file.exists(file.path(validated_inputs$input_dir, "cell_boundaries.parquet")) ||
    file.exists(file.path(validated_inputs$input_dir, "segmentation_boundaries.parquet"))
  if (seg_type %in% c("cell", "both") && !have_cell_seg) {
    if (validated_inputs$verbose) {
      .log_info("createSCOPE", "S11", "Cell segmentation parquet not found; downgrade segmentation strategy", validated_inputs$verbose)
    }
    seg_type <- if (file.exists(file.path(validated_inputs$input_dir, "nucleus_boundaries.parquet"))) "nucleus" else "none"
  }
  if (seg_type %in% c("nucleus", "both") && !file.exists(file.path(validated_inputs$input_dir, "nucleus_boundaries.parquet"))) {
    if (validated_inputs$verbose) {
      .log_info("createSCOPE", "S11", "Nucleus segmentation parquet not found; downgrade segmentation strategy", validated_inputs$verbose)
    }
    seg_type <- if (have_cell_seg) "cell" else "none"
  }

  if (need_temp) {
    link_or_copy <- function(src, dst) {
      if (file.exists(src) && !file.exists(dst)) {
        if (!suppressWarnings(file.symlink(src, dst))) {
          file.copy(src, dst, overwrite = FALSE)
        }
      }
    }
    link_or_copy(file.path(validated_inputs$input_dir, "segmentation_boundaries.parquet"),
                 file.path(work_dir, "segmentation_boundaries.parquet"))
    link_or_copy(file.path(validated_inputs$input_dir, "cell_boundaries.parquet"),
                 file.path(work_dir, "cell_boundaries.parquet"))
    link_or_copy(file.path(validated_inputs$input_dir, "nucleus_boundaries.parquet"),
                 file.path(work_dir, "nucleus_boundaries.parquet"))
  }

  if (generated_cells) {
    cells_path <- file.path(work_dir, "cells.parquet")
  } else {
    cells_path <- file.path(validated_inputs$input_dir, "cells.parquet")
  }
  if (!file.exists(cells_path)) {
    cells_source <- "missing"
  }
  payload <- list(
    work_dir = work_dir,
    grid_length = validated_inputs$grid_sizes,
    segmentation_strategy = seg_type,
    dataset_id = validated_inputs$dataset_id,
    data_type = validated_inputs$data_type,
    transcripts_path = transcripts_path,
    cells_path = cells_path,
    flatfile_dir = if (length(spec$flatfile_dirs)) spec$flatfile_dirs[1] else NA_character_,
    generated_transcripts = generated_transcripts,
    transcripts_source = transcripts_source,
    generated_cells = generated_cells,
    cells_source = cells_source,
    tx_file_used = tx_file_used,
    polygon_file_used = polygon_path_used
  )
  .trace_contract_field(payload, "grid_length", "scope_input.materialize")
  payload
}

#' Scope Input Validate Outputs
#' @description
#' Internal helper for `.scope_input_validate_outputs`.
#' @param payload Parameter value.
#' @param spec Parameter value.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.scope_input_validate_outputs <- function(payload, spec, validated_inputs) {
  if (!dir.exists(payload$work_dir)) {
    stop("Final work_dir does not exist: ", payload$work_dir)
  }
  if (!file.exists(payload$transcripts_path)) {
    stop("transcripts.parquet missing after module execution.")
  }
  if (!file.exists(payload$cells_path)) {
    stop("cells.parquet missing after module execution.")
  }
  TRUE
}

#' Run Scope Input Module
#' @description
#' Internal helper for `.run_scope_input_module`.
#' @param fixture Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_scope_input_module <- function(fixture) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for .scope_input_module")
  }
  validated <- .scope_input_validate_inputs(
    input_dir = fixture$input_dir,
    grid_sizes = fixture$grid_sizes,
    segmentation_strategy = if (!is.null(fixture$segmentation_strategy)) fixture$segmentation_strategy else "cell",
    allow_flatfile_generation = fixture$allow_flatfile_generation %||% TRUE,
    derive_cells_from_polygons = fixture$derive_cells_from_polygons %||% TRUE,
    data_type = fixture$data_type %||% "cosmx",
    pixel_size_um = if (!is.null(fixture$pixel_size_um)) fixture$pixel_size_um else 0.120280945,
    dataset_id = fixture$dataset_id,
    verbose = fixture$verbose %||% FALSE
  )
  spec <- .scope_input_spec_build(validated)
  payload <- .scope_input_materialize(spec, validated)
  .trace_contract_field(payload, "grid_length", "run_scope_input_module.materialize")
  .scope_input_validate_outputs(payload, spec, validated)
  fixture_id <- fixture$fixture_id
  if (is.null(fixture_id) || !nzchar(fixture_id)) {
    fixture_id <- basename(spec$input_dir)
  }
  digest_fields <- c(
    "grid_length",
    "segmentation_strategy",
    "dataset_id",
    "data_type",
    "flatfile_dir",
    "generated_transcripts",
    "transcripts_source",
    "generated_cells",
    "cells_source",
    "tx_file_used",
    "polygon_file_used"
  )
  payload_for_digest <- payload[intersect(names(payload), digest_fields)]
  res <- payload
  res$module_name <- ".scope_input_module"
  res$fixture_id <- fixture_id
  res$spec_fields <- paste(sort(names(spec)), collapse = ",")
  res$payload_fields <- paste(sort(names(payload_for_digest)), collapse = ",")
  res$spec_digest <- digest::digest(spec, algo = "sha256")
  res$payload_digest <- digest::digest(payload_for_digest, algo = "sha256")
  .trace_contract_field(res, "grid_length", "run_scope_input_module.output")
  res
}

#' Fdr Runner Validate Inputs
#' @description
#' Internal helper for `.fdr_runner_validate_inputs`.
#' @param Z_mat Parameter value.
#' @param P Parameter value.
#' @param perms Parameter value.
#' @param chunk_size Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param p_values Parameter value.
#' @param p_adj_mode Parameter value.
#' @param total_universe Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fdr_runner_validate_inputs <- function(Z_mat = NULL,
                                        P = NULL,
                                        perms = 0L,
                                        chunk_size = 32L,
                                        verbose = FALSE,
                                        p_values = NULL,
                                        p_adj_mode = "bh",
                                        total_universe = NULL) {
  if (!is.null(Z_mat) && !is.null(p_values)) {
    stop("Provide either Z_mat/P (matrix mode) or p_values (vector mode), not both.")
  }
  if (is.null(Z_mat) && is.null(p_values)) {
    stop("Either Z_mat (matrix mode) or p_values (vector mode) must be supplied.")
  }
  mode <- if (!is.null(Z_mat)) "matrix" else "vector"
  if (mode == "matrix") {
    if (!is.matrix(Z_mat) && !inherits(Z_mat, "big.matrix")) {
      stop("Z_mat must be a matrix or big.matrix.")
    }
    dims <- list(nrow = nrow(Z_mat), ncol = ncol(Z_mat))
  } else {
    if (!is.numeric(p_values)) stop("p_values must be numeric.")
    dims <- list(length = length(p_values))
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  list(
    mode = mode,
    Z_mat = if (mode == "matrix") Z_mat else NULL,
    P = if (mode == "matrix") P else NULL,
    perms = as.integer(max(0L, perms)),
    chunk_size = chunk_size,
    verbose = isTRUE(verbose),
    p_values = if (mode == "vector") as.numeric(p_values) else NULL,
    p_adj_mode = tolower(p_adj_mode),
    total_universe = total_universe,
    dims = dims,
    permutation_exists = !is.null(P)
  )
}

#' Fdr Runner Spec Build
#' @description
#' Internal helper for `.fdr_runner_spec_build`.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fdr_runner_spec_build <- function(validated_inputs) {
  list(
    module_name = ".fdr_runner_module",
    mode = validated_inputs$mode,
    dims = validated_inputs$dims,
    chunk_size = validated_inputs$chunk_size,
    permutation = list(
      perms = validated_inputs$perms,
      has_permutations = validated_inputs$permutation_exists
    ),
    p_adj_mode = validated_inputs$p_adj_mode,
    universe = validated_inputs$total_universe
  )
}

#' Fdr Runner Materialize
#' @description
#' Internal helper for `.fdr_runner_materialize`.
#' @param spec Parameter value.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fdr_runner_materialize <- function(spec, validated_inputs) {
  adjust_universe <- function(base_len, universe) {
    if (!is.null(universe) && is.finite(universe) && universe > 0) {
      universe
    } else {
      base_len
    }
  }

  adjust_mode_label <- function(mode) {
    switch(mode,
      "bh" = "BH",
      "by" = "BY",
      "bh_universe" = "BH (universe)",
      "by_universe" = "BY (universe)",
      "bonferroni" = "Bonferroni",
      toupper(mode)
    )
  }

  apply_adjustment <- function(values, mode, universe) {
    switch(mode,
      "by" = p.adjust(values, "BY"),
      "bh_universe" = p.adjust(values, "BH", n = adj_universe),
      "by_universe" = p.adjust(values, "BY", n = adj_universe),
      "bonferroni" = p.adjust(values, "bonferroni", n = adj_universe),
      p.adjust(values, "BH")
    )
  }

  mode <- spec$mode
  perms <- validated_inputs$perms
  min_p_possible <- if (validated_inputs$permutation_exists) 1 / (perms + 1) else NA_real_
  if (mode == "matrix") {
    Z_mat <- validated_inputs$Z_mat
    P <- validated_inputs$P
    chunk_size <- validated_inputs$chunk_size
    verbose <- validated_inputs$verbose
    is_big <- inherits(Z_mat, "big.matrix")
    n_genes <- ncol(Z_mat)
    idx_chunks <- seq(1, n_genes, by = chunk_size)
    allocate <- function() {
      bigmemory::big.matrix(
        nrow = nrow(Z_mat),
        ncol = n_genes,
        type = "double",
        init = NA_real_
      )
    }
    if (is_big && !requireNamespace("bigmemory", quietly = TRUE)) {
      stop("bigmemory is required for chunked matrix FDR smoothing")
    }
    FDR_disc <- if (is_big) allocate() else matrix(NA_real_, nrow = nrow(Z_mat), ncol = n_genes)
    FDR_beta <- if (is_big) allocate() else matrix(NA_real_, nrow = nrow(Z_mat), ncol = n_genes)
    FDR_mid <- if (is_big) allocate() else matrix(NA_real_, nrow = nrow(Z_mat), ncol = n_genes)
    FDR_uniform <- if (is_big) allocate() else matrix(NA_real_, nrow = nrow(Z_mat), ncol = n_genes)
    for (start in idx_chunks) {
      end <- min(start + chunk_size - 1L, n_genes)
      if (!is.null(P)) {
        p_disc <- P[, start:end, drop = FALSE]
        k_est <- round(p_disc * (perms + 1) - 1)
        k_est[k_est < 0] <- 0
        k_est[k_est > perms] <- perms
        p_beta <- (k_est + 1) / (perms + 2)
        p_mid <- (k_est + 0.5) / (perms + 1)
        p_uniform <- (k_est + matrix(runif(length(k_est)), nrow = nrow(k_est))) / (perms + 1)
      } else {
        Z_chunk <- Z_mat[, start:end, drop = FALSE]
        p_disc <- 2 * pnorm(-abs(Z_chunk))
        p_beta <- p_mid <- p_uniform <- p_disc
      }
      flat <- function(x) matrix(p.adjust(c(x), "BH"), nrow = nrow(x))
      if (is_big) {
        FDR_disc[, start:end] <- flat(p_disc)
        FDR_beta[, start:end] <- flat(p_beta)
        FDR_mid[, start:end] <- flat(p_mid)
        FDR_uniform[, start:end] <- flat(p_uniform)
      } else {
        FDR_disc[, start:end] <- flat(p_disc)
        FDR_beta[, start:end] <- flat(p_beta)
        FDR_mid[, start:end] <- flat(p_mid)
        FDR_uniform[, start:end] <- flat(p_uniform)
      }
    }
    FDR_storey <- NULL
    pi0_hat <- NA_real_
    FDR_main_method <- if (is_big) "BH(beta p) (chunked, no q-value)" else "Storey q-value (beta p)"
    if (!is_big) {
      reference <- if (!is.null(P)) p.adjust(P, "BH") * NA_real_ else NA_real_
      p_vec <- as.numeric(if (!is.null(P)) P else 2 * pnorm(-abs(Z_mat)))
      m_tot <- length(p_vec)
      if (m_tot > 0) {
        lambdas <- seq(0.5, 0.95, by = 0.05)
        pi0_vals <- sapply(lambdas, function(lam) {
          mean(p_vec > lam) / (1 - lam)
        })
        pi0_hat <- min(1, min(pi0_vals, na.rm = TRUE))
        if (!is.finite(pi0_hat) || pi0_hat <= 0) pi0_hat <- 1
        o <- order(p_vec, na.last = NA)
        ro <- integer(m_tot); ro[o] <- seq_along(o)
        q_raw <- pi0_hat * m_tot * p_vec / pmax(ro, 1)
        q_ord <- q_raw[o]
        for (i in seq_len(length(q_ord) - 1)) {
          idx <- length(q_ord) - i
          if (q_ord[idx] > q_ord[idx + 1]) q_ord[idx] <- q_ord[idx + 1]
        }
        q_final <- numeric(m_tot)
        q_final[o] <- pmin(q_ord, 1)
        FDR_storey <- matrix(q_final, nrow = nrow(FDR_beta), dimnames = dimnames(FDR_beta))
      }
    }
    FDR_main <- if (!is.null(FDR_storey)) FDR_storey else FDR_beta
    n_sig_005 <- sum(FDR_main < 0.05, na.rm = TRUE)
    return(list(
      FDR_out_disc = FDR_disc,
      FDR_out_beta = FDR_beta,
      FDR_out_mid = FDR_mid,
      FDR_out_uniform = FDR_uniform,
      FDR_storey = FDR_storey,
      FDR_main = FDR_main,
      FDR_main_method = FDR_main_method,
      pi0_hat = pi0_hat,
      n_sig_005 = n_sig_005,
      min_p_possible = min_p_possible,
      raw_p = if (!is.null(P)) P else 2 * pnorm(-abs(Z_mat))
    ))
  } else {
    p_values <- validated_inputs$p_values
    return(list(
      FDR_main = apply_adjustment(p_values, validated_inputs$p_adj_mode, validated_inputs$total_universe),
      FDR_main_method = adjust_mode_label(validated_inputs$p_adj_mode),
      pi0_hat = NA_real_,
      n_sig_005 = sum(apply_adjustment(p_values, validated_inputs$p_adj_mode, validated_inputs$total_universe) < 0.05, na.rm = TRUE),
      min_p_possible = min_p_possible,
      raw_p = p_values
    ))
  }
}

#' Fdr Runner Validate Outputs
#' @description
#' Internal helper for `.fdr_runner_validate_outputs`.
#' @param payload Parameter value.
#' @param spec Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fdr_runner_validate_outputs <- function(payload, spec) {
  if (!is.list(payload)) stop("Module output must be a list.")
  if (spec$mode == "matrix") {
    dims <- spec$dims
    check <- function(mat, name) {
      if (inherits(mat, "big.matrix")) return()
      if (!is.matrix(mat)) stop(name, " must be a matrix.")
      if (!identical(dim(mat), c(dims$nrow, dims$ncol))) stop(name, " dimensions must match spec.")
    }
    check(payload$FDR_out_disc, "FDR_out_disc")
    check(payload$FDR_out_beta, "FDR_out_beta")
    check(payload$FDR_out_mid, "FDR_out_mid")
    check(payload$FDR_out_uniform, "FDR_out_uniform")
  } else {
    if (!is.numeric(payload$FDR_main)) stop("FDR_main must be numeric.")
    if (length(payload$raw_p) != length(payload$FDR_main)) stop("p length mismatch.")
  }
  TRUE
}

#' Run Fdr Runner Module
#' @description
#' Internal helper for `.run_fdr_runner_module`.
#' @param fixture Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_fdr_runner_module <- function(fixture) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for .fdr_runner_module")
  }
  object_path <- fixture$input_rds
  if (is.null(object_path) || !nzchar(object_path) || !file.exists(object_path)) {
    stop("fixture input_rds must point to an existing RDS file.")
  }
  data <- readRDS(object_path)
  validated <- do.call(.fdr_runner_validate_inputs, data$inputs)
  spec <- .fdr_runner_spec_build(validated)
  payload <- .fdr_runner_materialize(spec, validated)
  .fdr_runner_validate_outputs(payload, spec)
  digest_fields <- names(payload)
  list(
    module_name = ".fdr_runner_module",
    fixture_id = fixture$fixture_id %||% "fdr_runner",
    spec_fields = paste(sort(names(spec)), collapse = ","),
    payload_fields = paste(sort(digest_fields), collapse = ","),
    spec_digest = digest::digest(spec, algo = "sha256"),
    payload_digest = digest::digest(payload, algo = "sha256"),
    pos_path = object_path
  )
}

#' Fdr Payload Finalize Validate Inputs
#' @description
#' Internal helper for `.fdr_payload_finalize_validate_inputs`.
#' @param mode Parameter value.
#' @param fdr_spec Parameter value.
#' @param fdr_payload Parameter value.
#' @param include_matrices Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fdr_payload_finalize_validate_inputs <- function(mode,
                                                  fdr_spec,
                                                  fdr_payload,
                                                  include_matrices = NULL) {
  mode <- match.arg(mode, c("matrix", "vector"))
  if (!is.list(fdr_spec) || fdr_spec$module_name != ".fdr_runner_module") {
    stop("fdr_spec must originate from .fdr_runner_module")
  }
  if (!is.list(fdr_payload)) stop("fdr_payload must be a list")
  base_required <- c("FDR_main", "FDR_main_method", "pi0_hat", "n_sig_005", "min_p_possible", "raw_p")
  missing <- setdiff(base_required, names(fdr_payload))
  if (length(missing)) stop("Missing ModuleF payload fields: ", paste(missing, collapse = ", "))
  include_matrices <- if (mode == "matrix") {
    TRUE
  } else {
    if (isTRUE(include_matrices)) stop("Vector mode cannot request matrix outputs")
    FALSE
  }
  list(
    mode = mode,
    include_matrices = include_matrices,
    fdr_spec = fdr_spec,
    fdr_payload = fdr_payload
  )
}

#' Fdr Payload Finalize Spec Build
#' @description
#' Internal helper for `.fdr_payload_finalize_spec_build`.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fdr_payload_finalize_spec_build <- function(validated_inputs) {
  list(
    module_name = ".fdr_payload_finalize_module",
    mode = validated_inputs$mode,
    include_matrices = validated_inputs$include_matrices,
    dims = validated_inputs$fdr_spec$dims
  )
}

#' Fdr Payload Finalize Materialize
#' @description
#' Internal helper for `.fdr_payload_finalize_materialize`.
#' @param spec Parameter value.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fdr_payload_finalize_materialize <- function(spec, validated_inputs) {
  payload <- validated_inputs$fdr_payload
  result <- list(
    FDR_main = payload$FDR_main,
    FDR_main_method = payload$FDR_main_method,
    pi0_hat = payload$pi0_hat,
    n_sig_005 = payload$n_sig_005,
    min_p_possible = payload$min_p_possible,
    raw_p = payload$raw_p
  )
  if (spec$include_matrices) {
    result$FDR_out_disc <- payload$FDR_out_disc
    result$FDR_out_beta <- payload$FDR_out_beta
    result$FDR_out_mid <- payload$FDR_out_mid
    result$FDR_out_uniform <- payload$FDR_out_uniform
    result$FDR_storey <- payload$FDR_storey
  }
  result
}

#' Fdr Payload Finalize Validate Outputs
#' @description
#' Internal helper for `.fdr_payload_finalize_validate_outputs`.
#' @param payload Parameter value.
#' @param spec Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fdr_payload_finalize_validate_outputs <- function(payload, spec) {
  if (!is.list(payload)) stop("Module output must be a list.")
  required <- c("FDR_main", "FDR_main_method", "pi0_hat", "n_sig_005", "min_p_possible", "raw_p")
  missing <- setdiff(required, names(payload))
  if (length(missing)) stop("Missing required output elements: ", paste(missing, collapse = ", "))

  if (!is.character(payload$FDR_main_method) || length(payload$FDR_main_method) != 1) {
    stop("FDR_main_method must be a single character string.")
  }
  if (!is.numeric(payload$pi0_hat) || length(payload$pi0_hat) != 1) {
    stop("pi0_hat must be a numeric scalar.")
  }
  if (!is.numeric(payload$n_sig_005) || length(payload$n_sig_005) != 1) {
    stop("n_sig_005 must be a numeric scalar.")
  }
  if (!is.numeric(payload$min_p_possible) || length(payload$min_p_possible) != 1) {
    stop("min_p_possible must be a numeric scalar.")
  }

  dims <- spec$dims
  if (spec$mode == "matrix") {
    expected_dims <- c(dims$nrow, dims$ncol)
    matrix_check <- function(mat, name) {
      if (!is.matrix(mat)) stop(name, " must be a matrix.")
      if (!identical(dim(mat), expected_dims)) stop(name, " dimensions must match spec.")
    }
    matrix_check(payload$FDR_main, "FDR_main")
    matrix_check(payload$raw_p, "raw_p")
    matrix_fields <- c("FDR_out_disc", "FDR_out_beta", "FDR_out_mid", "FDR_out_uniform")
    missing_mat <- setdiff(matrix_fields, names(payload))
    if (length(missing_mat)) stop("Missing matrix fields: ", paste(missing_mat, collapse = ", "))
    for (field in matrix_fields) {
      matrix_check(payload[[field]], field)
    }
    if (!is.null(payload$FDR_storey)) {
      matrix_check(payload$FDR_storey, "FDR_storey")
    }
  } else {
    if (!is.numeric(payload$FDR_main)) stop("FDR_main must be numeric.")
    if (!is.numeric(payload$raw_p)) stop("raw_p must be numeric.")
    if (length(payload$FDR_main) != dims$length) stop("FDR_main length must match spec.")
    if (length(payload$raw_p) != dims$length) stop("raw_p length must match spec.")
  }
  TRUE
}

#' Run Fdr Payload Finalize Module
#' @description
#' Internal helper for `.run_fdr_payload_finalize_module`.
#' @param fixture Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_fdr_payload_finalize_module <- function(fixture) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for .fdr_payload_finalize_module")
  }
  object_path <- fixture$input_rds
  if (is.null(object_path) || !nzchar(object_path) || !file.exists(object_path)) {
    stop("fixture input_rds must point to an existing RDS file.")
  }
  data <- readRDS(object_path)
  validated <- do.call(.fdr_payload_finalize_validate_inputs, data$inputs)
  spec <- .fdr_payload_finalize_spec_build(validated)
  payload <- .fdr_payload_finalize_materialize(spec, validated)
  .fdr_payload_finalize_validate_outputs(payload, spec)
  digest_fields <- names(payload)
  list(
    module_name = ".fdr_payload_finalize_module",
    fixture_id = fixture$fixture_id %||% "fdr_payload_finalize",
    spec_fields = paste(sort(names(spec)), collapse = ","),
    payload_fields = paste(sort(digest_fields), collapse = ","),
    spec_digest = digest::digest(spec, algo = "sha256"),
    payload_digest = digest::digest(payload, algo = "sha256"),
    pos_path = object_path
  )
}

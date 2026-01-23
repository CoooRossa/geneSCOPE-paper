#' Visium Scope Validate Inputs
#' @description
#' Internal helper for `.visium_scope_validate_inputs`.
#' @param input_dir Filesystem path.
#' @param grid_multiplier Parameter value.
#' @param include_in_tissue_spots Parameter value.
#' @param flip_y Parameter value.
#' @param data_type Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_scope_validate_inputs <- function(input_dir,
                                          grid_multiplier,
                                          include_in_tissue_spots,
                                          flip_y,
                                          data_type) {
  stopifnot(dir.exists(input_dir))
  if (!is.numeric(grid_multiplier) || length(grid_multiplier) != 1L || !is.finite(grid_multiplier) || grid_multiplier <= 0) {
    stop("grid_multiplier must be a single positive finite numeric value.")
  }
  list(
    input_dir = input_dir,
    grid_multiplier = as.numeric(grid_multiplier),
    include_in_tissue_spots = isTRUE(include_in_tissue_spots),
    flip_y = isTRUE(flip_y),
    data_type = data_type
  )
}

#' Visium Scope Spec Build
#' @description
#' Internal helper for `.visium_scope_spec_build`.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_scope_spec_build <- function(validated_inputs) {
  spatial_dir <- file.path(validated_inputs$input_dir, "spatial")
  if (!dir.exists(spatial_dir)) stop("spatial/ directory not found under ", validated_inputs$input_dir)
  candidate_dirs <- unique(file.path(validated_inputs$input_dir, c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")))
  list(
    input_dir = validated_inputs$input_dir,
    data_type = validated_inputs$data_type,
    spatial_dir = spatial_dir,
    matrix_candidates = candidate_dirs,
    pos_candidates = file.path(spatial_dir, c("tissue_positions_list.csv", "tissue_positions.csv")),
    scalefactor_candidates = file.path(spatial_dir, c("scalefactors_json.json", "scalefactors_json_fullres.json")),
    include_in_tissue_spots = validated_inputs$include_in_tissue_spots,
    flip_y = validated_inputs$flip_y,
    grid_multiplier = validated_inputs$grid_multiplier
  )
}

#' Visium Scope Payload Build
#' @description
#' Internal helper for `.visium_scope_payload_build`.
#' @param spec Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_scope_payload_build <- function(spec) {
  pos_path <- spec$pos_candidates[file.exists(spec$pos_candidates)][1]
  if (is.na(pos_path)) stop("Could not find tissue_positions_list.csv or tissue_positions.csv under ", spec$spatial_dir)
  if (is.na(spec$scalefactor_candidates[file.exists(spec$scalefactor_candidates)][1])) stop("scalefactors_json.json not found under ", spec$spatial_dir)
  matrix_dir <- spec$matrix_candidates[dir.exists(spec$matrix_candidates)][1]
  feature_file <- NA_character_
  if (!is.na(matrix_dir)) {
    gz_path <- file.path(matrix_dir, "features.tsv.gz")
    txt_path <- file.path(matrix_dir, "features.tsv")
    if (file.exists(gz_path)) {
      feature_file <- gz_path
    } else if (file.exists(txt_path)) {
      feature_file <- txt_path
    }
  }
  list(
    pos_path = pos_path,
    scalefactor_path = spec$scalefactor_candidates[file.exists(spec$scalefactor_candidates)][1],
    matrix_dir = matrix_dir,
    feature_file = feature_file
  )
}

#' Normalize Flag Value
#' @description
#' Internal helper for `.normalize_flag_value`.
#' @param value Parameter value.
#' @return Return value used internally.
#' @keywords internal
.normalize_flag_value <- function(value) {
  if (isTRUE(value)) return(TRUE)
  if (isFALSE(value)) return(FALSE)
  if (is.character(value)) {
    if (tolower(value) %in% c("true", "t", "1")) return(TRUE)
    if (tolower(value) %in% c("false", "f", "0")) return(FALSE)
  }
  if (is.numeric(value)) {
    return(!is.na(value) && value != 0)
  }
  NA
}

#' Run Visium Scope Module
#' @description
#' Internal helper for `.run_visium_scope_module`.
#' @param fixture Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_visium_scope_module <- function(fixture) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for visium scope module validation")
  }
  validated <- .visium_scope_validate_inputs(
    input_dir = fixture$input_dir,
    grid_multiplier = if (!is.null(fixture$grid_multiplier)) as.numeric(fixture$grid_multiplier) else 1,
    include_in_tissue_spots = .normalize_flag_value(fixture$include_in_tissue_spots),
    flip_y = .normalize_flag_value(fixture$flip_y),
    data_type = if (!is.null(fixture$data_type) && nzchar(fixture$data_type)) fixture$data_type else "visium"
  )
  spec <- .visium_scope_spec_build(validated)
  payload <- .visium_scope_payload_build(spec)
  fixture_id <- fixture$fixture_id
  if (is.null(fixture_id) || !nzchar(fixture_id)) {
    fixture_id <- basename(spec$input_dir)
  }
  list(
    module_name = ".visium_scope_module",
    fixture_id = fixture_id,
    spec_fields = paste(sort(names(spec)), collapse = ","),
    payload_fields = paste(sort(names(payload)), collapse = ","),
    spec_digest = digest::digest(spec, algo = "sha256"),
    payload_digest = digest::digest(payload, algo = "sha256"),
    pos_path = payload$pos_path,
    scalefactor_path = payload$scalefactor_path,
    matrix_dir = payload$matrix_dir,
    feature_file = payload$feature_file
  )
}

# ModuleD - stats table postprocess helpers (PURE layer, used by .compute_l)
#' Stats Table Validate Inputs
#' @description
#' Internal helper for `.stats_table_validate_inputs`.
#' @param L Parameter value.
#' @param X_full Parameter value.
#' @param X_used Parameter value.
#' @param W Parameter value.
#' @param grid_inf Parameter value.
#' @param cells Parameter value.
#' @param genes Parameter value.
#' @param within Parameter value.
#' @param chunk_size Parameter value.
#' @param L_min Numeric threshold.
#' @param block_id Parameter value.
#' @param perms Parameter value.
#' @param block_size Parameter value.
#' @param current_cores Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.stats_table_validate_inputs <- function(L,
                                         X_full,
                                         X_used,
                                         W,
                                         grid_inf,
                                         cells,
                                         genes,
                                         within,
                                         chunk_size,
                                         L_min,
                                         block_id,
                                         perms,
                                         block_size,
                                         current_cores,
                                         verbose) {
  if (!is.matrix(L) && !inherits(L, "big.matrix")) {
    stop("L must be a matrix or big.matrix.")
  }
  if (nrow(L) != ncol(L)) stop("L must be square.")
  if (!is.matrix(X_full)) stop("X_full must be a matrix.")
  if (!is.matrix(X_used)) stop("X_used must be a matrix.")

  n_cells <- nrow(X_full)
  n_genes <- ncol(L)
  if (ncol(X_used) != n_genes) stop("X_used column count must match L dimension.")
  if (ncol(X_full) != n_genes) stop("X_full column count must match L dimension.")

  if (!is(W, "Matrix") && !is.matrix(W)) {
    stop("W must be a matrix or Matrix object.")
  }
  if (nrow(W) != n_cells || ncol(W) != n_cells) stop("W must be square with rows equal to cells.")

  if (length(block_id) != n_cells) stop("block_id length must equal number of cells.")
  if (!is.numeric(chunk_size) || chunk_size <= 0) stop("chunk_size must be a positive number.")
  if (!is.numeric(perms) || perms < 0) stop("perms must be a non-negative number.")
  if (!is.numeric(block_size) || block_size <= 0) stop("block_size must be a positive number.")
  if (!is.numeric(current_cores) || current_cores < 1) stop("current_cores must be >= 1.")

  required_grid <- c("grid_id", "xmin", "xmax", "ymin", "ymax")
  if (!is.data.frame(grid_inf) || !all(required_grid %in% colnames(grid_inf))) {
    stop("grid_inf must be a data.frame with grid_id, xmin/xmax, and ymin/ymax.")
  }
  if (length(cells) != n_cells) stop("cells length must equal number of rows in X_full.")
  if (!all(cells %in% grid_inf$grid_id)) stop("cells must be a subset of grid_inf$grid_id.")

  gene_names <- if (!is.null(genes)) {
    as.character(genes)
  } else if (!is.null(dimnames(L)[[1]])) {
    as.character(dimnames(L)[[1]])
  } else if (!is.null(colnames(X_full))) {
    as.character(colnames(X_full))
  } else {
    seq_len(n_genes)
  }

  sample_names <- if (!is.null(rownames(X_full))) {
    rownames(X_full)
  } else {
    seq_len(n_cells)
  }

  list(
    L = L,
    X_full = X_full,
    X_used = X_used,
    W = W,
    grid_inf = grid_inf,
    cells = as.character(cells),
    genes = if (is.null(genes)) NULL else as.character(genes),
    within = isTRUE(within),
    chunk_size = as.integer(chunk_size),
    L_min = as.numeric(L_min),
    block_id = as.integer(block_id),
    perms = as.integer(perms),
    block_size = as.integer(block_size),
    current_cores = max(1L, as.integer(current_cores)),
    verbose = isTRUE(verbose),
    gene_order = gene_names,
    sample_order = as.character(sample_names),
    n_cells = n_cells,
    n_genes = n_genes
  )
}

#' Stats Table Postprocess Spec Build
#' @description
#' Internal helper for `.stats_table_postprocess_spec_build`.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stats_table_postprocess_spec_build <- function(validated_inputs) {
  grid_inf <- validated_inputs$grid_inf
  list(
    module_name = ".stats_table_postprocess_module",
    dims = list(cells = validated_inputs$n_cells, genes = validated_inputs$n_genes),
    gene_order = validated_inputs$gene_order,
    sample_order = validated_inputs$sample_order,
    within_mode = validated_inputs$within,
    chunk_size = validated_inputs$chunk_size,
    L_min = validated_inputs$L_min,
    permutation = list(perms = validated_inputs$perms, block_size = validated_inputs$block_size),
    block_summary = list(
      unique_blocks = length(unique(validated_inputs$block_id)),
      block_size = validated_inputs$block_size
    ),
    grid_bounds = list(
      xmin = min(grid_inf$xmin, na.rm = TRUE),
      xmax = max(grid_inf$xmax, na.rm = TRUE),
      ymin = min(grid_inf$ymin, na.rm = TRUE),
      ymax = max(grid_inf$ymax, na.rm = TRUE)
    ),
    invariant_columns = c("Z_mat", "P", "FDR_out_disc", "FDR_out_beta", "FDR_out_mid", "FDR_out_uniform"),
    metadata_columns = c("FDR_storey", "FDR_main", "FDR_main_method", "pi0_hat", "n_sig_005", "min_p_possible", "betas", "qc"),
    fdr_adjustments = c("BH", "Storey q")
  )
}

#' Stats Table Postprocess Materialize
#' @description
#' Internal helper for `.stats_table_postprocess_materialize`.
#' @param spec Parameter value.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stats_table_postprocess_materialize <- function(spec, validated_inputs) {
  L <- validated_inputs$L
  X_full <- validated_inputs$X_full
  X_used <- validated_inputs$X_used
  W <- validated_inputs$W
  grid_inf <- validated_inputs$grid_inf
  cells <- validated_inputs$cells
  genes <- validated_inputs$genes
  within <- validated_inputs$within
  chunk_size <- validated_inputs$chunk_size
  L_min <- validated_inputs$L_min
  block_id <- validated_inputs$block_id
  perms <- validated_inputs$perms
  block_size <- validated_inputs$block_size
  current_cores <- validated_inputs$current_cores
  verbose <- validated_inputs$verbose

  n <- nrow(X_full)

  ## --- 3. Analytical Z computation with memory optimization ---
  if (within || is.null(genes)) {
    S0 <- sum(W)
    EZ <- -1 / (n - 1)
    Var <- (n^2 * (n - 2)) / ((n - 1)^2 * (n - 3) * S0)

    if (inherits(L, "big.matrix")) {
      if (verbose) .log_info("computeL", "S04", "Computing Z-scores in chunks", verbose)
      Z_mat <- L
      n_genes <- ncol(L)
      chunk_genes <- seq(1, n_genes, by = chunk_size)
      for (start in chunk_genes) {
        end <- min(start + chunk_size - 1, n_genes)
        L_chunk <- L[, start:end]
        Z_chunk <- (L_chunk - EZ) / sqrt(Var)
        Z_mat[, start:end] <- Z_chunk
      }
    } else {
      Z_mat <- (L - EZ) / sqrt(Var)
    }
  } else {
    if (inherits(L, "big.matrix")) {
      stop("Asymmetric Z-score computation not yet supported for big.matrix")
    }
    Z_mat <- t(apply(L, 1, function(v) {
      sdv <- sd(v)
      if (sdv == 0) rep(0, length(v)) else (v - mean(v)) / sdv
    }))
    dimnames(Z_mat) <- dimnames(L)
  }

  ## --- 4. Monte Carlo p-values with BLAS control and error recovery ---
  P <- if (perms > 0) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
    }

    L_reg <- if (inherits(L, "big.matrix")) {
      if (verbose) .log_info("computeL", "S05", "Converting to regular matrix for permutation tests", verbose)
      as.matrix(L)
    } else {
      L
    }

    if (verbose) .log_info("computeL", "S05", "Running permutation tests", verbose)
    perm_success <- FALSE
    perm_cores <- current_cores
    perm_attempt <- 1

    while (!perm_success && perm_cores >= 1) {
      if (verbose && perm_attempt > 1) {
        .log_info("computeL", "S05", paste0("Retry #", perm_attempt, " with ", perm_cores, " cores"), verbose)
      }

      if (ncol(X_used) != nrow(L_reg) || ncol(X_used) != ncol(L_reg)) {
        stop(sprintf("Permutation input mismatch: ncol(X_used)=%d but L_ref is %dx%d",
          ncol(X_used), nrow(L_reg), ncol(L_reg)))
      }

      perm_result <- tryCatch(
        {
          t_start <- Sys.time()
          p_result <- .lee_l_perm_block(X_used, W, L_reg,
            block_id = block_id,
            perms = perms,
            block_size = block_size,
            n_threads = perm_cores
          )
          t_end <- Sys.time()
          if (verbose) {
            time_msg <- if (perm_attempt == 1) "completed" else "retry successful"
            .log_info("computeL", "S05", paste0("Permutation test ", time_msg, " (", format(t_end - t_start), ")"), verbose)
          }
          list(success = TRUE, object = p_result)
        },
        error = function(e) {
          if (verbose && perm_attempt > 1) {
            .log_info("computeL", "S05", paste0("Test failed: ", conditionMessage(e)), verbose)
          }
          list(success = FALSE, error = e)
        }
      )

      if (perm_result$success) {
        perm_success <- TRUE
        P <- perm_result$object
      } else {
        perm_attempt <- perm_attempt + 1
        perm_cores <- max(floor(perm_cores / 2), 1)
        if (verbose) {
          .log_info("computeL", "S05", paste0("Reducing cores to ", perm_cores, " and retrying"), verbose)
        }
        Sys.sleep(2)
        gc(verbose = FALSE)
      }
    }

    if (!perm_success) {
      if (verbose) {
        .log_info("computeL", "S05", "!!! Warning: Permutation test failed; setting P = NULL !!!", verbose)
      }
      NULL
    } else {
      P
    }
  } else {
    NULL
  }

  ## --- 5. FDR: Three Monte Carlo smoothing strategies as in .get_top_l_vs_r ---
  ## --- 5. FDR smoothing via ModuleF ---
  fdr_validated <- .fdr_runner_validate_inputs(
    Z_mat = Z_mat,
    P = P,
    perms = perms,
    chunk_size = chunk_size,
    verbose = verbose
  )
  fdr_spec <- .fdr_runner_spec_build(fdr_validated)
  fdr_payload <- .fdr_runner_materialize(fdr_spec, fdr_validated)
  .fdr_runner_validate_outputs(fdr_payload, fdr_spec)
  fdr_finalize_validated <- .fdr_payload_finalize_validate_inputs(
    mode = fdr_spec$mode,
    fdr_spec = fdr_spec,
    fdr_payload = fdr_payload,
    include_matrices = TRUE
  )
  fdr_finalize_spec <- .fdr_payload_finalize_spec_build(fdr_finalize_validated)
  fdr_final <- .fdr_payload_finalize_materialize(fdr_finalize_spec, fdr_finalize_validated)
  .fdr_payload_finalize_validate_outputs(fdr_final, fdr_finalize_spec)
  FDR_out_disc <- fdr_final$FDR_out_disc
  FDR_out_beta <- fdr_final$FDR_out_beta
  FDR_out_mid <- fdr_final$FDR_out_mid
  FDR_out_uniform <- fdr_final$FDR_out_uniform
  FDR_storey <- fdr_final$FDR_storey
  FDR_main <- fdr_final$FDR_main
  FDR_main_method <- fdr_final$FDR_main_method
  pi0_hat <- fdr_final$pi0_hat
  n_sig_005 <- fdr_final$n_sig_005
  min_p_possible <- fdr_final$min_p_possible

  ## --- 6. betax / betay ---
  centres <- with(
    grid_inf[match(cells, grid_inf$grid_id), ],
    data.frame(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2
    )
  )
  betas <- t(apply(X_full, 2, function(v) {
    c(beta_x = coef(lm(v ~ centres$x + centres$y))[2:3])
  }))
  if (!is.null(genes)) betas <- betas[genes, , drop = FALSE]
  if (!is.null(spec$gene_order) && length(spec$gene_order) == nrow(betas) &&
      all(spec$gene_order %in% rownames(betas))) {
    betas <- betas[spec$gene_order, , drop = FALSE]
    if (!is.null(rownames(betas))) rownames(betas) <- spec$gene_order
  }

  ## --- 7. QC computation ---
  qc <- NULL
  if (within || is.null(genes)) {
    L_qc <- if (inherits(L, "big.matrix")) {
      .log_info("computeL", "S07", "Converting subset of big.matrix for QC computation...", TRUE)
      n_sample <- min(1000, nrow(L))
      sample_idx <- sample(nrow(L), n_sample)
      as.matrix(L[sample_idx, sample_idx])
    } else {
      L
    }

    A_bin <- abs(L_qc) >= L_min
    A_bin <- A_bin | t(A_bin)
    diag(A_bin) <- FALSE

    A_num <- (A_bin | t(A_bin)) * 1
    g_tmp <- igraph::simplify(
      igraph::graph_from_adjacency_matrix(A_num,
        mode = "undirected",
        diag = FALSE
      ),
      remove.multiple = TRUE,
      remove.loops = TRUE
    )

    deg <- igraph::degree(g_tmp)
    memb <- tryCatch(
      igraph::cluster_louvain(g_tmp)$membership,
      error = function(e) {
        igraph::cluster_leiden(g_tmp)$membership
      }
    )

    qc <- list(
      edge_density = 2 * sum(A_bin[upper.tri(A_bin)]) /
        (ncol(A_bin) * (ncol(A_bin) - 1)),
      components = igraph::components(g_tmp)$no,
      modularity_Q = igraph::modularity(g_tmp, membership = memb),
      mean_degree = mean(deg),
      sd_degree = sd(deg),
      hub_ratio = mean(deg > 2 * median(deg)),
      sig_edge_frac = if (is.null(P)) {
        NA
      } else {
        if (inherits(P, "big.matrix")) {
          p_sample <- as.matrix(P[sample_idx, sample_idx])
          mean(p_sample[upper.tri(p_sample)] < 0.05, na.rm = TRUE)
        } else {
          mean(P[upper.tri(P)] < 0.05, na.rm = TRUE)
        }
      }
    )
  }

  list(
    Z_mat = Z_mat,
    P = P,
    FDR_out_disc = FDR_out_disc,
    FDR_out_beta = FDR_out_beta,
    FDR_out_mid = FDR_out_mid,
    FDR_out_uniform = FDR_out_uniform,
    FDR_storey = FDR_storey,
    FDR_main = FDR_main,
    FDR_main_method = FDR_main_method,
    pi0_hat = pi0_hat,
    n_sig_005 = n_sig_005,
    min_p_possible = min_p_possible,
    betas = betas,
    qc = qc
  )
}

#' Stats Table Validate Outputs
#' @description
#' Internal helper for `.stats_table_validate_outputs`.
#' @param out_tbl Parameter value.
#' @param spec Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stats_table_validate_outputs <- function(out_tbl, spec) {
  if (!is.list(out_tbl)) stop("Output must be a list.")
  required <- c("Z_mat", "FDR_out_disc", "FDR_out_beta", "FDR_out_mid", "FDR_out_uniform", "FDR_main", "FDR_main_method", "pi0_hat", "n_sig_005", "min_p_possible", "betas", "qc")
  if (length(setdiff(required, names(out_tbl)))) stop("Missing required output elements: ", paste(setdiff(required, names(out_tbl)), collapse = ", "))

  check_matrix <- function(mat, name) {
    if (!is.matrix(mat)) stop(name, " must be a matrix.")
    if (!identical(
        dim(mat),
        c(spec$dims$genes, spec$dims$genes)
    )) stop(
        name,
        " must be ",
        spec$dims$genes,
        "x",
        spec$dims$genes
    )
    if (!is.null(spec$gene_order) && !is.null(rownames(mat)) &&
        !identical(rownames(mat), spec$gene_order)) {
      stop(name, " rownames must match gene_order.")
    }
    if (!is.null(spec$gene_order) && !is.null(colnames(mat)) &&
        !identical(colnames(mat), spec$gene_order)) {
      stop(name, " colnames must match gene_order.")
    }
  }

  check_matrix(out_tbl$Z_mat, "Z_mat")
  check_matrix(out_tbl$FDR_out_disc, "FDR_out_disc")
  check_matrix(out_tbl$FDR_out_beta, "FDR_out_beta")
  check_matrix(out_tbl$FDR_out_mid, "FDR_out_mid")
  check_matrix(out_tbl$FDR_out_uniform, "FDR_out_uniform")
  check_matrix(out_tbl$FDR_main, "FDR_main")

  if (spec$permutation$perms > 0) {
    if (is.null(out_tbl$P)) stop("P must be present when perms > 0.")
    check_matrix(out_tbl$P, "P")
  } else if (!is.null(out_tbl$P)) {
    stop("P must be NULL when perms == 0.")
  }

  if (!is.null(out_tbl$FDR_storey)) {
    check_matrix(out_tbl$FDR_storey, "FDR_storey")
  }

  if (!is.matrix(out_tbl$betas)) stop("betas must be a matrix.")
  if (nrow(out_tbl$betas) != spec$dims$genes) stop("betas row count must match gene dimension.")
  if (!is.null(spec$gene_order) && !is.null(rownames(out_tbl$betas)) &&
      !identical(rownames(out_tbl$betas), spec$gene_order)) {
    stop("betas rownames must match gene_order.")
  }

  if (!is.character(out_tbl$FDR_main_method) || length(out_tbl$FDR_main_method) != 1) {
    stop("FDR_main_method must be a single character string.")
  }
  if (!is.numeric(out_tbl$pi0_hat) || length(out_tbl$pi0_hat) != 1) stop("pi0_hat must be a numeric scalar.")
  if (!is.numeric(out_tbl$n_sig_005) || length(out_tbl$n_sig_005) != 1) stop("n_sig_005 must be numeric scalar.")
  if (!is.numeric(out_tbl$min_p_possible) || length(out_tbl$min_p_possible) != 1) stop("min_p_possible must be numeric scalar.")

  if (!is.null(out_tbl$qc) && !is.list(out_tbl$qc)) stop("qc must be NULL or a list.")
  TRUE
}

#' Run Stats Table Postprocess Module
#' @description
#' Internal helper for `.run_stats_table_postprocess_module`.
#' @param fixture Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_stats_table_postprocess_module <- function(fixture) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for stats table postprocess module validation")
  }
  fixture_id <- fixture$fixture_id
  if (is.null(fixture_id) || !nzchar(fixture_id)) {
    fixture_id <- "sample_stats_postprocess"
  }
  input_rds <- fixture$input_rds
  if (is.null(input_rds) || !nzchar(input_rds) || !file.exists(input_rds)) {
    stop("input_rds must point to an existing file.")
  }
  data <- readRDS(input_rds)
  validated <- do.call(.stats_table_validate_inputs, data)
  spec <- .stats_table_postprocess_spec_build(validated)
  payload <- .stats_table_postprocess_materialize(spec, validated)
  .stats_table_validate_outputs(payload, spec)
  list(
    module_name = ".stats_table_postprocess_module",
    fixture_id = fixture_id,
    spec_fields = paste(sort(names(spec)), collapse = ","),
    payload_fields = paste(sort(names(payload)), collapse = ","),
    spec_digest = digest::digest(spec, algo = "sha256"),
    payload_digest = digest::digest(payload, algo = "sha256"),
    pos_path = input_rds
  )
}

#' Stage2 Native Inputs Validate Inputs
#' @description
#' Internal helper for `.stage2_native_inputs_validate_inputs`.
#' @param kept_genes Parameter value.
#' @param similarity_matrix Parameter value.
#' @param config Parameter value.
#' @param stage1_membership_labels Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage2_native_inputs_validate_inputs <- function(kept_genes, similarity_matrix, config, stage1_membership_labels) {
  if (length(kept_genes) == 0L) stop("kept_genes must contain at least one entry.")
  if (!any(inherits(similarity_matrix, "Matrix")) && !is.matrix(similarity_matrix)) stop("similarity_matrix must be a Matrix or matrix.")
  matrix_names <- rownames(similarity_matrix)
  if (is.null(matrix_names)) stop("similarity_matrix must have rownames.")
  if (!all(kept_genes %in% matrix_names)) stop("kept_genes must be a subset of similarity_matrix rownames.")
  list(
    kept_genes = as.character(kept_genes),
    similarity_dim = dim(similarity_matrix),
    stage1_membership_labels = stage1_membership_labels,
    config_dimensions = length(config)
  )
}

#' Stage2 Native Inputs Spec Build
#' @description
#' Internal helper for `.stage2_native_inputs_spec_build`.
#' @param validated_inputs Parameter value.
#' @param config Parameter value.
#' @param stage1_membership_labels Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage2_native_inputs_spec_build <- function(validated_inputs, config, stage1_membership_labels) {
  stage1_summary <- list(
    total_labels = length(validated_inputs$stage1_membership_labels),
    unique_labels = sort(unique(validated_inputs$stage1_membership_labels))
  )
  list(
    kept_genes = validated_inputs$kept_genes,
    config_summary = list(
      min_cutoff = as.numeric(config$min_cutoff),
      pct_min = config$pct_min,
      use_significance = isTRUE(config$use_significance),
      significance_max = as.numeric(config$significance_max),
      CI95_filter = isTRUE(config$CI95_filter),
      CI_rule = config$CI_rule,
      use_mh_weight = isTRUE(config$use_mh_weight),
      use_log1p_weight = isTRUE(config$use_log1p_weight),
      post_smooth = isTRUE(config$post_smooth),
      post_smooth_quant = config$post_smooth_quant,
      post_smooth_power = config$post_smooth_power,
      keep_stage1_backbone = isTRUE(config$keep_stage1_backbone),
      backbone_floor_q = config$backbone_floor_q
    ),
    stage1_summary = stage1_summary
  )
}

#' Stage2 Native Inputs Materialize
#' @description
#' Internal helper for `.stage2_native_inputs_materialize`.
#' @param spec Parameter value.
#' @param similarity_matrix Parameter value.
#' @param FDR Parameter value.
#' @param aux_stats Parameter value.
#' @param pearson_matrix Parameter value.
#' @param mh_object Parameter value.
#' @param stage1_membership_labels Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.stage2_native_inputs_materialize <- function(spec,
                                              similarity_matrix,
                                              FDR,
                                              aux_stats,
                                              pearson_matrix,
                                              mh_object,
                                              stage1_membership_labels,
                                              config,
                                              verbose) {
  kept_genes <- spec$kept_genes
  L_post <- similarity_matrix[kept_genes, kept_genes, drop = FALSE]
  L_post[abs(L_post) < config$min_cutoff] <- 0
  if (config$use_significance && !is.null(FDR)) {
    L_post <- tryCatch(
      .align_and_filter_fdr(L_post, similarity_matrix, FDR, config$significance_max),
      error = function(e) {
        .cluster_message(verbose, "[cluster]   significance subset disabled: ", conditionMessage(e))
        L_post
      }
    )
  }
  L_post[L_post < 0] <- 0
  diag(L_post) <- 0
  L_post <- .filter_matrix_by_quantile(L_post, config$pct_min, "q100")
  if (!inherits(L_post, "sparseMatrix")) {
    L_post <- drop0(Matrix(L_post, sparse = TRUE))
  }

  if (config$CI95_filter) {
    if (is.null(aux_stats)) stop("CI95_filter requires aux_stats input.")
    curve_obj <- aux_stats[[config$curve_layer]]
    if (is.null(curve_obj)) stop("CI95 curve missing in stats layer.")
    curve_mat <- if (inherits(curve_obj, "big.matrix")) {
      bigmemory::as.matrix(curve_obj)
    } else if (is.matrix(curve_obj)) curve_obj else as.matrix(as.data.frame(curve_obj))
    xp <- as.numeric(curve_mat[, "Pear"])
    lo <- as.numeric(curve_mat[, "lo95"])
    hi <- as.numeric(curve_mat[, "hi95"])
    if (anyNA(c(xp, lo, hi))) stop("LR curve contains NA")
    f_lo <- approxfun(xp, lo, rule = 2)
    f_hi <- approxfun(xp, hi, rule = 2)
    if (is.null(pearson_matrix)) stop("CI95_filter requires pearson_matrix input.")
    rMat <- as.matrix(pearson_matrix[kept_genes, kept_genes, drop = FALSE])
    L_post[
        if (config$CI_rule == "remove_within") (L_post <= f_hi(rMat)) else (L_post < f_lo(rMat) | L_post > f_hi(rMat))
    ] <- 0
  }

  if (config$use_mh_weight) {
    if (is.null(mh_object)) stop("MH weighting requested but mh_object is NULL.")
    if (inherits(mh_object, "igraph")) {
      ed_sim <- igraph::as_data_frame(mh_object, what = "edges")
      wcol <- if ("MH" %in% names(ed_sim)) "MH" else if ("CMH" %in% names(ed_sim)) "CMH" else if ("weight" %in% names(ed_sim)) "weight" else stop("MH graph lacks weight attribute.")
      key_sim <- paste(pmin(ed_sim$from, ed_sim$to), pmax(ed_sim$from, ed_sim$to), sep = "|")
      sim_map <- setNames(ed_sim[[wcol]], key_sim)
      idx <- .arrind_from_matrix_predicate(
        L_post,
        op = "gt",
        cutoff = 0,
        triangle = "upper",
        keep_diag = FALSE
      )
      if (nrow(idx)) {
        g1_names <- kept_genes[idx[, 1]]
        g2_names <- kept_genes[idx[, 2]]
        k <- paste(pmin(g1_names, g2_names), pmax(g1_names, g2_names), sep = "|")
        s <- sim_map[k]
        s[is.na(s)] <- median(sim_map, na.rm = TRUE)
        L_post[cbind(idx[, 1], idx[, 2])] <- L_post[cbind(idx[, 1], idx[, 2])] * s
        L_post[cbind(idx[, 2], idx[, 1])] <- L_post[cbind(idx[, 2], idx[, 1])] * s
      }
    } else if (inherits(mh_object, "Matrix") || is.matrix(mh_object)) {
      M <- mh_object
      if (inherits(M, "Matrix") && !inherits(M, "CsparseMatrix")) M <- as(M, "CsparseMatrix")
      if (!(identical(rownames(M), rownames(similarity_matrix)) && identical(colnames(M), colnames(similarity_matrix)))) {
        if (!is.null(rownames(M)) && !is.null(colnames(M)) &&
            all(rownames(similarity_matrix) %in% rownames(M)) &&
            all(colnames(similarity_matrix) %in% colnames(M))) {
          M <- M[rownames(similarity_matrix), colnames(similarity_matrix), drop = FALSE]
        } else {
          stop("MH matrix cannot be aligned to similarity matrix dimensions.")
        }
      }
      ridx <- match(kept_genes, rownames(M))
      M_sub <- M[ridx, ridx, drop = FALSE]
      idx <- .arrind_from_matrix_predicate(
        L_post,
        op = "gt",
        cutoff = 0,
        triangle = "upper",
        keep_diag = FALSE
      )
      if (nrow(idx)) {
        s <- M_sub[cbind(idx[, 1], idx[, 2])]
        vec_all <- if (inherits(M_sub, "sparseMatrix")) M_sub@x else as.numeric(M_sub)
        vec_all <- vec_all[is.finite(vec_all)]
        s[is.na(s)] <- median(vec_all, na.rm = TRUE)
        L_post[cbind(idx[, 1], idx[, 2])] <- L_post[cbind(idx[, 1], idx[, 2])] * s
        L_post[cbind(idx[, 2], idx[, 1])] <- L_post[cbind(idx[, 2], idx[, 1])] * s
      }
    } else {
      stop("Unsupported mh_object class: ", paste(class(mh_object), collapse = ","))
    }
  }

  L_post <- pmax(L_post, t(L_post))
  diag(L_post) <- 0

  if (isTRUE(config$keep_stage1_backbone)) {
    pos_vals <- as.numeric(L_post[L_post > 0])
    if (length(pos_vals) == 0) {
      w_floor <- 1e-6
    } else {
      w_floor <- max(
        1e-8,
        as.numeric(quantile(pos_vals, probs = min(max(config$backbone_floor_q, 0), 0.25), na.rm = TRUE))
      )
    }
    for (cid in sort(na.omit(unique(stage1_membership_labels)))) {
      genes_c <- kept_genes[stage1_membership_labels[kept_genes] == cid]
      if (length(genes_c) < 2) next
      M <- L_post[genes_c, genes_c, drop = FALSE]
      if (!any(M > 0, na.rm = TRUE)) {
        M <- pmax(similarity_matrix, 0)[genes_c, genes_c, drop = FALSE]
      }
      idx_back <- .arrind_from_matrix_predicate(
        M,
        op = "gt",
        cutoff = 0,
        triangle = "upper",
        keep_diag = FALSE
      )
      if (nrow(idx_back) == 0) next
      mst_edge_table <- data.frame(from = genes_c[idx_back[, 1]], to = genes_c[idx_back[, 2]], w = M[idx_back], stringsAsFactors = FALSE)
      gtmp <- igraph::graph_from_data_frame(mst_edge_table, directed = FALSE, vertices = genes_c)
      if (igraph::ecount(gtmp) == 0) next
      mst_c <- igraph::mst(gtmp, weights = 1 / (igraph::E(gtmp)$w + 1e-9))
      if (igraph::ecount(mst_c) == 0) next
      ep <- igraph::as_data_frame(mst_c, what = "edges")
      for (k in seq_len(nrow(ep))) {
        i <- ep$from[k]; j <- ep$to[k]
        L_post[i, j] <- max(L_post[i, j], w_floor)
        L_post[j, i] <- L_post[i, j]
      }
    }
  }

  ed1 <- summary(L_post)
  ed1 <- ed1[ed1$i < ed1$j & ed1$x > 0, , drop = FALSE]
  if (!nrow(ed1)) stop("No edges remain after Stage-2 corrections.")
  edges_corr <- data.frame(
    from = kept_genes[ed1$i],
    to = kept_genes[ed1$j],
    L_corr = ed1$x,
    stringsAsFactors = FALSE
  )
  use_log1p_weight <- isTRUE(config$use_log1p_weight)
  edges_corr$w_raw <- if (use_log1p_weight) log1p(edges_corr$L_corr) else edges_corr$L_corr
  if (isTRUE(config$post_smooth)) {
    q <- quantile(edges_corr$w_raw, probs = config$post_smooth_quant, na.rm = TRUE)
    w_clipped <- pmin(pmax(edges_corr$w_raw, q[1]), q[2])
    rng <- q[2] - q[1]
    if (!is.finite(rng) || rng <= 0) rng <- max(1e-8, max(w_clipped) - min(w_clipped))
    w01 <- (w_clipped - min(w_clipped)) / (rng + 1e-12)
    if (!is.null(config$post_smooth_power) && is.finite(config$post_smooth_power) && config$post_smooth_power != 1) {
      w01 <- w01 ^ config$post_smooth_power
    }
    edges_corr$weight <- w01
  } else {
    edges_corr$weight <- edges_corr$w_raw
  }
  corrected_similarity_graph <- igraph::graph_from_data_frame(
    edges_corr[, c("from", "to", "weight")],
    directed = FALSE,
    vertices = kept_genes
  )
  W <- as.matrix(igraph::as_adjacency_matrix(corrected_similarity_graph, attr = "weight", sparse = TRUE))
  rownames(W) <- colnames(W) <- kept_genes

  list(
    L_post = L_post,
    corrected_similarity_graph = corrected_similarity_graph,
    edges_corr = edges_corr,
    W = W,
    use_log1p_weight = use_log1p_weight
  )
}

#' Stage2 Native Inputs Validate Outputs
#' @description
#' Internal helper for `.stage2_native_inputs_validate_outputs`.
#' @param payload Parameter value.
#' @param spec Parameter value.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage2_native_inputs_validate_outputs <- function(payload, spec, validated_inputs) {
  kept_genes <- spec$kept_genes
  if (!identical(rownames(payload$L_post), kept_genes)) stop("L_post rownames do not match kept_genes.")
  if (!identical(rownames(payload$W), kept_genes)) stop("W rownames do not match kept_genes.")
  graph_genes <- igraph::V(payload$corrected_similarity_graph)$name
  if (!identical(sort(graph_genes), sort(kept_genes))) stop("Graph vertex names drift.")
  if (any(!payload$edges_corr$from %in% kept_genes) || any(!payload$edges_corr$to %in% kept_genes)) {
    stop("Edge table contains genes outside kept_genes.")
  }
  TRUE
}

#' Run Stage2 Native Inputs Module
#' @description
#' Internal helper for `.run_stage2_native_inputs_module`.
#' @param fixture Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_stage2_native_inputs_module <- function(fixture) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for stage2 native inputs validation")
  }
  data <- readRDS(fixture$object_path)
  validated <- .stage2_native_inputs_validate_inputs(
    kept_genes = data$kept_genes,
    similarity_matrix = data$similarity_matrix,
    config = data$config,
    stage1_membership_labels = data$stage1_membership_labels
  )
  spec <- .stage2_native_inputs_spec_build(validated, data$config, data$stage1_membership_labels)
  payload <- .stage2_native_inputs_materialize(
    spec,
    similarity_matrix = data$similarity_matrix,
    FDR = data$FDR,
    aux_stats = data$aux_stats,
    pearson_matrix = data$pearson_matrix,
    mh_object = data$mh_object,
    stage1_membership_labels = data$stage1_membership_labels,
    config = data$config,
    verbose = FALSE
  )
  .stage2_native_inputs_validate_outputs(payload, spec, validated)
  list(
    module_name = ".stage2_native_inputs_module",
    fixture_id = fixture$fixture_id,
    spec_fields = paste("kept_genes", length(spec$kept_genes), sep = "="),
    payload_fields = paste(sort(names(payload)), collapse = ","),
    spec_digest = digest::digest(spec, algo = "sha256"),
    payload_digest = digest::digest(list(L_post = payload$L_post, edges_corr = payload$edges_corr, W = payload$W), algo = "sha256"),
    kept_genes = paste(spec$kept_genes, collapse = ","),
    edge_count = nrow(payload$edges_corr)
  )
}

#' Stage2 Refine Blocks Validate Inputs
#' @description
#' Internal helper for `.stage2_refine_blocks_validate_inputs`.
#' @param kept_genes Parameter value.
#' @param stage1_membership_labels Parameter value.
#' @param corrected_similarity_graph Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage2_refine_blocks_validate_inputs <- function(kept_genes,
                                                 stage1_membership_labels,
                                                 corrected_similarity_graph) {
  if (length(kept_genes) == 0L) stop("kept_genes must contain at least one entry.")
  if (!inherits(corrected_similarity_graph, "igraph")) stop("corrected_similarity_graph must be an igraph object.")
  graph_names <- igraph::V(corrected_similarity_graph)$name
  if (is.null(graph_names) || length(graph_names) == 0L) {
    stop("corrected_similarity_graph must have vertex names.")
  }
  if (!all(kept_genes %in% graph_names)) {
    stop("kept_genes must appear in corrected_similarity_graph vertex names.")
  }
  if (is.null(names(stage1_membership_labels))) {
    names(stage1_membership_labels) <- kept_genes
  }
  if (!all(kept_genes %in% names(stage1_membership_labels))) {
    stop("stage1_membership_labels must cover all kept_genes.")
  }
  stage1_membership_labels <- stage1_membership_labels[kept_genes]
  list(
    kept_genes = as.character(kept_genes),
    stage1_membership_labels = stage1_membership_labels,
    graph_names = graph_names
  )
}

#' Stage2 Refine Blocks Spec Build
#' @description
#' Internal helper for `.stage2_refine_blocks_spec_build`.
#' @param validated_inputs Parameter value.
#' @param stage2_algo_final Parameter value.
#' @param stage2_backend Parameter value.
#' @param runtime_cfg Parameter value.
#' @param config Parameter value.
#' @param use_log1p_weight Logical flag.
#' @return Return value used internally.
#' @keywords internal
.stage2_refine_blocks_spec_build <- function(validated_inputs,
                                            stage2_algo_final,
                                            stage2_backend,
                                            runtime_cfg,
                                            config,
                                            use_log1p_weight) {
  stage1_summary <- list(
    total_labels = length(validated_inputs$stage1_membership_labels),
    unique_labels = sort(unique(validated_inputs$stage1_membership_labels))
  )
  config_summary <- list(
    use_consensus = isTRUE(runtime_cfg$use_consensus),
    prefer_fast = isTRUE(config$prefer_fast),
    enable_subcluster = isTRUE(config$enable_subcluster),
    hotspot_k = config$hotspot_k,
    hotspot_min_module_size = config$hotspot_min_module_size,
    use_log1p_weight = isTRUE(use_log1p_weight)
  )
  list(
    kept_genes = validated_inputs$kept_genes,
    stage1_summary = stage1_summary,
    stage2_algo_final = stage2_algo_final,
    stage2_backend = stage2_backend,
    runtime_cfg = runtime_cfg,
    config_summary = config_summary
  )
}

#' Stage2 Refine Blocks Materialize
#' @description
#' Internal helper for `.stage2_refine_blocks_materialize`.
#' @param spec Parameter value.
#' @param corrected_similarity_graph Parameter value.
#' @param stage1_membership_labels Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param use_log1p_weight Logical flag.
#' @return Return value used internally.
#' @keywords internal
.stage2_refine_blocks_materialize <- function(spec,
                                             corrected_similarity_graph,
                                             stage1_membership_labels,
                                             config,
                                             verbose,
                                             use_log1p_weight) {
  runtime_cfg <- spec$runtime_cfg
  stage2_algo_final <- spec$stage2_algo_final
  stage2_backend <- spec$stage2_backend
  kept_genes <- spec$kept_genes

  use_consensus <- runtime_cfg$use_consensus
  consensus_thr <- runtime_cfg$consensus_thr
  n_restart <- runtime_cfg$n_restart
  n_threads <- runtime_cfg$n_threads
  objective <- runtime_cfg$objective
  res_param_base <- runtime_cfg$res_param_base
  mode_setting <- runtime_cfg$mode_setting
  large_n_threshold <- runtime_cfg$large_n_threshold
  aggr_future_workers <- runtime_cfg$aggr_future_workers
  aggr_batch_size <- runtime_cfg$aggr_batch_size
  nk_leiden_iterations <- runtime_cfg$nk_leiden_iterations
  nk_leiden_randomize <- runtime_cfg$nk_leiden_randomize
  stage1_mode_selected <- runtime_cfg$stage1_mode_selected

  ng <- length(kept_genes)
  memb_final <- rep(NA_integer_, ng)
  names(memb_final) <- kept_genes
  cons <- Matrix(0, ng, ng, sparse = TRUE, dimnames = list(kept_genes, kept_genes))

  next_global_id <- 1L
  for (cid in sort(na.omit(unique(stage1_membership_labels)))) {
    idx <- which(stage1_membership_labels == cid)
    genes_c <- kept_genes[idx]
    if (!length(genes_c)) next
    g_sub <- igraph::induced_subgraph(corrected_similarity_graph, vids = genes_c)
    if (length(genes_c) <= 1 || igraph::ecount(g_sub) == 0) {
      memb_final[genes_c] <- next_global_id
      cons[genes_c, genes_c] <- 0
      next_global_id <- next_global_id + 1L
      next
    }

    if (identical(stage2_algo_final, "hotspot-like")) {
      S_sub <- igraph::as_adjacency_matrix(g_sub, attr = "weight", sparse = TRUE)
      S_sub <- as(S_sub, "dgCMatrix")
      hotspot_k <- config$hotspot_k
      if (is.null(hotspot_k) || !is.finite(hotspot_k)) hotspot_k <- length(genes_c)
      hotspot_res <- .run_hotspot_clustering(S_sub, as.integer(max(1L, min(hotspot_k, nrow(S_sub)))), config$hotspot_min_module_size, use_log1p_weight)
      memb_sub <- hotspot_res$membership
      if (is.null(names(memb_sub))) names(memb_sub) <- genes_c
      cons_sub <- hotspot_res$consensus
      if (is.null(cons_sub)) {
        cons_sub <- Matrix(0, length(memb_sub), length(memb_sub), sparse = TRUE,
          dimnames = list(names(memb_sub), names(memb_sub)))
      }
    } else {
      mode_sub <- if (identical(mode_setting, "auto")) {
        if (length(genes_c) > large_n_threshold) "aggressive" else stage1_mode_selected
      } else stage1_mode_selected
      algo_per_run_sub <- stage2_algo_final
      nk_inputs_sub <- .nk_prepare_graph_input(g_sub)
      g_sub <- nk_inputs_sub$graph
      ew <- nk_inputs_sub$edge_weights
      el <- nk_inputs_sub$edge_table
      vertex_names_sub <- nk_inputs_sub$vertex_names
      res_param_sub <- res_param_base

      backend_sub <- if (identical(stage2_backend, "networkit")) "networkit" else "igraph"
      if (backend_sub == "networkit" && !stage2_algo_final %in% c("leiden", "louvain")) {
        backend_sub <- "igraph"
      }
      if (isTRUE(config$prefer_fast)) backend_sub <- "igraph"
      if (backend_sub == "networkit" && !.nk_available()) {
        stop("Stage-2 NetworKit backend requested but NetworKit is unavailable.", call. = FALSE)
      }

      if (backend_sub == "networkit") {
        seeds <- seq_len(n_restart)
        chunk_size <- aggr_batch_size
        if (is.null(chunk_size) || chunk_size <= 0) {
          chunk_size <- ceiling(n_restart / max(1L, aggr_future_workers))
        }
        split_chunks <- function(x, k) {
          if (k <= 0) return(list(x))
          split(x, ceiling(seq_along(x) / k))
        }
        chunks <- split_chunks(seeds, chunk_size)
        threads_per_worker <- max(1L, floor(n_threads / max(1L, aggr_future_workers)))
        worker_fun <- function(ss) {
          if (identical(stage2_algo_final, "leiden")) {
            .nk_parallel_leiden_mruns(
              vertex_names = vertex_names_sub,
              edge_table = el,
              edge_weights = ew,
              gamma = as.numeric(res_param_sub),
              iterations = nk_leiden_iterations,
              randomize = nk_leiden_randomize,
              threads = threads_per_worker,
              seeds = as.integer(ss)
            )
          } else {
            .nk_plm_mruns(
              vertex_names = vertex_names_sub,
              edge_table = el,
              edge_weights = ew,
              gamma = as.numeric(res_param_sub),
              refine = TRUE,
              threads = threads_per_worker,
              seeds = as.integer(ss)
            )
          }
        }
        memb_blocks <- if (aggr_future_workers > 1L) {
          future.apply::future_lapply(chunks, worker_fun, future.seed = TRUE)
        } else {
          lapply(chunks, worker_fun)
        }
        memb_mat_sub <- do.call(cbind, memb_blocks)
      } else {
        memb_mat_sub <- future.apply::future_sapply(
          seq_len(n_restart),
          function(ii) {
            g_local <- igraph::graph_from_data_frame(
              cbind(el, weight = ew),
              directed = FALSE,
              vertices = vertex_names_sub
            )
            res <- try(.run_graph_algorithm(g_local, algo_per_run_sub, res_param_sub, objective), silent = TRUE)
            if (inherits(res, "try-error")) {
              rep.int(1L, length(genes_c))
            } else {
              mm <- try(igraph::membership(res), silent = TRUE)
              if (inherits(mm, "try-error") || is.null(mm)) rep.int(1L, length(genes_c)) else mm
            }
          },
          future.seed = TRUE
        )
      }
      if (is.null(dim(memb_mat_sub))) {
        memb_mat_sub <- matrix(memb_mat_sub, ncol = 1)
      }
      rownames(memb_mat_sub) <- igraph::V(g_sub)$name

      if (!use_consensus) {
        .cluster_message(verbose, "[cluster]   Stage-2 consensus backend: single-run (no aggregation)")
        memb_sub <- memb_mat_sub[, 1]
        labs <- as.integer(memb_sub)
        n_loc <- length(labs)
        split_idx <- split(seq_len(n_loc), labs)
        if (!length(split_idx)) {
          cons_sub <- Matrix(0, n_loc, n_loc, sparse = TRUE, dimnames = list(names(memb_sub), names(memb_sub)))
        } else {
          ii <- integer(0); jj <- integer(0)
          for (idv in split_idx) {
            s <- length(idv)
            if (s > 1) {
              comb <- combn(idv, 2)
              ii <- c(ii, comb[1, ])
              jj <- c(jj, comb[2, ])
            }
          }
          cons_sub <- if (length(ii)) {
            sparseMatrix(i = ii, j = jj, x = rep(1, length(ii)), dims = c(n_loc, n_loc), symmetric = TRUE)
          } else {
            Matrix(0, n_loc, n_loc, sparse = TRUE)
          }
          if (!is.null(names(memb_sub)) && length(names(memb_sub)) == n_loc) {
            dimnames(cons_sub) <- list(names(memb_sub), names(memb_sub))
          }
        }
      } else {
        cons_sub <- try(consensus_sparse(memb_mat_sub,
          thr = consensus_thr,
          n_threads = n_threads
        ), silent = TRUE)
        if (!inherits(cons_sub, "try-error")) {
          .cluster_message(verbose, "[cluster]   Stage-2 consensus backend: C++ (consensus_sparse)")
        } else {
          .cluster_message(verbose, "[cluster]   Stage-2 consensus backend: R fallback (.compute_consensus_sparse_r)")
          cons_sub <- .compute_consensus_sparse_r(memb_mat_sub, thr = consensus_thr)
        }
        if (length(cons_sub@i)) diag(cons_sub) <- 0
        if (nnzero(cons_sub) == 0) {
          memb_sub <- memb_mat_sub[, 1]
        } else {
          g_cons_sub <- igraph::graph_from_adjacency_matrix(cons_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
          backend_consensus_sub <- backend_sub
          if (!backend_consensus_sub %in% c("networkit", "igraph")) backend_consensus_sub <- "igraph"
          backend_label <- if (identical(backend_consensus_sub, "networkit")) {
            paste0("NetworKit (", stage2_algo_final, ")")
          } else {
            paste0("igraph (", stage2_algo_final, ")")
          }
          .cluster_message(verbose, "[cluster]   Stage-2 consensus graph backend: ", backend_label)
          nk_opts_cons <- list(
            iterations = config$nk_leiden_iterations,
            randomize = config$nk_leiden_randomize,
            refine = TRUE,
            seeds = 1L
          )
          memb_sub <- .run_single_partition(
            g_cons_sub,
            backend = backend_consensus_sub,
            algo_name = stage2_algo_final,
            res_param = res_param_base,
            objective = objective,
            threads = n_threads,
            nk_opts = nk_opts_cons
          )
        }
      }
    }

    sub_ids <- sort(unique(memb_sub))
    if (length(sub_ids) <= 1) {
      memb_final[genes_c] <- next_global_id
      next_global_id <- next_global_id + 1L
    } else {
      for (sid in sub_ids) {
        genes_sid <- names(memb_sub)[memb_sub == sid]
        memb_final[genes_sid] <- next_global_id
        next_global_id <- next_global_id + 1L
      }
    }

    block_names <- intersect(rownames(cons_sub), genes_c)
    if (length(block_names)) {
      cons[block_names, block_names] <- as.matrix(cons_sub[block_names, block_names])
    } else {
      cons[genes_c, genes_c] <- 0
    }
  }

  list(
    memb_final = memb_final,
    cons = cons
  )
}

#' Stage2 Refine Blocks Validate Outputs
#' @description
#' Internal helper for `.stage2_refine_blocks_validate_outputs`.
#' @param payload Parameter value.
#' @param spec Parameter value.
#' @param validated_inputs Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage2_refine_blocks_validate_outputs <- function(payload, spec, validated_inputs) {
  kept_genes <- spec$kept_genes
  if (!identical(names(payload$memb_final), kept_genes)) stop("memb_final names do not match kept_genes.")
  if (!identical(rownames(payload$cons), kept_genes)) stop("cons rownames do not match kept_genes.")
  if (!identical(colnames(payload$cons), kept_genes)) stop("cons colnames do not match kept_genes.")
  TRUE
}

#' Run Stage2 Refine Blocks Module
#' @description
#' Internal helper for `.run_stage2_refine_blocks_module`.
#' @param fixture Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_stage2_refine_blocks_module <- function(fixture) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for stage2 refine blocks validation")
  }
  data <- readRDS(fixture$object_path)
  validated <- .stage2_refine_blocks_validate_inputs(
    kept_genes = data$kept_genes,
    stage1_membership_labels = data$stage1_membership_labels,
    corrected_similarity_graph = data$corrected_similarity_graph
  )
  spec <- .stage2_refine_blocks_spec_build(
    validated,
    stage2_algo_final = data$stage2_algo_final,
    stage2_backend = data$stage2_backend,
    runtime_cfg = data$runtime_cfg,
    config = data$config,
    use_log1p_weight = data$use_log1p_weight
  )
  payload <- .stage2_refine_blocks_materialize(
    spec,
    corrected_similarity_graph = data$corrected_similarity_graph,
    stage1_membership_labels = data$stage1_membership_labels,
    config = data$config,
    verbose = FALSE,
    use_log1p_weight = data$use_log1p_weight
  )
  .stage2_refine_blocks_validate_outputs(payload, spec, validated)
  list(
    module_name = ".stage2_refine_blocks_module",
    fixture_id = fixture$fixture_id,
    spec_fields = paste(sort(names(spec)), collapse = ","),
    payload_fields = paste(sort(names(payload)), collapse = ","),
    spec_digest = digest::digest(spec, algo = "sha256"),
    payload_digest = digest::digest(list(memb_final = payload$memb_final, cons = payload$cons), algo = "sha256"),
    memb_final = paste(payload$memb_final, collapse = ","),
    cons_nnz = nnzero(payload$cons)
  )
}

#' Module Contract Placeholder
#' @description
#' Internal helper for `.module_contract_placeholder`.
#' @return Return value used internally.
#' @keywords internal
.module_contract_placeholder <- function() {
  # Placeholder for future pure spec helpers; keep implementations deterministic and side-effect free.
  NULL
}

#' Module Runner Placeholder
#' @description
#' Internal helper for `.module_runner_placeholder`.
#' @return Return value used internally.
#' @keywords internal
.module_runner_placeholder <- function() {
  # Placeholder for impure runners; keep I/O/RNG/parallel side effects confined here.
  invisible(NULL)
}

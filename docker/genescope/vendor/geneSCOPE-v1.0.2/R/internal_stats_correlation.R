
#' Cmh Lookup Pairs
#' @description
#' Internal helper for `.cmh_lookup_pairs`.
#' @param kept_genes Parameter value.
#' @param ei Parameter value.
#' @param ej Parameter value.
#' @param g_sim Parameter value.
#' @param weight_col Parameter value.
#' @param n_threads Number of threads to use.
#' @return Return value used internally.
#' @keywords internal
.cmh_lookup_pairs <- function(kept_genes, ei, ej, g_sim, weight_col = c("CMH", "weight"), n_threads = NULL) {
    weight_col <- match.arg(weight_col)
    if (is.null(n_threads)) {
        # Robust threads selection across versions
        n_threads <- tryCatch({
            if (exists(".safe_thread_count", mode = "function", inherits = TRUE)) {
                .safe_thread_count()
            } else if (exists(".get_safe_thread_count", mode = "function", inherits = TRUE)) {
                .get_safe_thread_count(default = 8L)
            } else if (exists(".get_safe_thread_count_v2", mode = "function", inherits = TRUE)) {
                f <- get(".get_safe_thread_count_v2", mode = "function")
                an <- tryCatch(names(formals(f)), error = function(e) character(0))
                if ("max_requested" %in% an) f(max_requested = 8L)
                else if ("requested_cores" %in% an) f(requested_cores = 8L)
                else f(8L)
            } else {
                max(1L, min(8L, .detect_cores_safe(logical = TRUE) - 1L))
            }
        }, error = function(e) 1L)
    }
    ed_sim <- igraph::as_data_frame(g_sim, what = "edges")
    wcol <- if (!is.null(ed_sim[[weight_col]])) weight_col else if (!is.null(ed_sim$weight)) "weight" else stop("CMH graph lacks weight column")
    si <- match(ed_sim$from, kept_genes)
    sj <- match(ed_sim$to, kept_genes)
    ok <- which(!is.na(si) & !is.na(sj))
    if (!length(ok)) {
        return(rep_len(1, length(ei)))
    }
    si <- as.integer(si[ok])
    sj <- as.integer(sj[ok])
    sw <- as.numeric(ed_sim[[wcol]][ok])
    fallback <- median(sw, na.rm = TRUE)
    if (!is.finite(fallback)) fallback <- 1
    if (exists("cmh_lookup_rcpp", mode = "function")) {
        return(cmh_lookup_rcpp(as.integer(ei), as.integer(ej), si, sj, sw, fallback, as.integer(n_threads)))
    }
    # R fallback: direct keyed lookup to maintain length(ei)
    sim_names <- paste(
        pmin(kept_genes[si], kept_genes[sj]),
        pmax(kept_genes[si], kept_genes[sj]),
        sep = "|"
    )
    w_map <- sw
    names(w_map) <- sim_names
    edge_keys <- paste(
        pmin(kept_genes[ei], kept_genes[ej]),
        pmax(kept_genes[ei], kept_genes[ej]),
        sep = "|"
    )
    ww <- w_map[edge_keys]
    ww[is.na(ww)] <- fallback
    as.numeric(ww)
}

#' Compute correlation between two sparse matrix blocks
#' @description
#' Internal helper for `.compute_block_correlation`.
#' @param block_i First block
#' @param block_j Second block
#' @param means_i Column means of first block
#' @param means_j Column means of second block
#' @param sds_i Column standard deviations of first block
#' @param sds_j Column standard deviations of second block
#' @return Correlation matrix block
#' @keywords internal
.compute_block_correlation <- function(block_i, block_j, means_i, means_j, sds_i, sds_j) {
  n_obs <- nrow(block_i)
  n_i <- ncol(block_i)
  n_j <- ncol(block_j)

  cor_block <- matrix(0, n_i, n_j)

  for (i in seq_len(n_i)) {
    for (j in seq_len(n_j)) {
      if (sds_i[i] == 0 || sds_j[j] == 0) {
        cor_block[i, j] <- 0
        next
      }

      # Safely extract column vectors and convert to numeric vectors
      x_col <- block_i[, i, drop = FALSE]
      y_col <- block_j[, j, drop = FALSE]

      # Convert to numeric vectors
      x_vec <- as.numeric(x_col)
      y_vec <- as.numeric(y_col)

      # Compute covariance - use standard formula
      if (all(x_vec == 0) && all(y_vec == 0)) {
        # Both columns are zero
        cov_xy <- 0
      } else {
        # Use standard covariance formula
        mean_x <- means_i[i]
        mean_y <- means_j[j]
        cov_xy <- sum((x_vec - mean_x) * (y_vec - mean_y)) / (n_obs - 1)
      }

      # Compute correlation coefficient
      cor_block[i, j] <- cov_xy / (sds_i[i] * sds_j[j])
    }
  }

  return(cor_block)
}

#' Sparse-safe Pearson correlation for dgCMatrix inputs
#' @description Computes Pearson correlation using sparse tcrossprod without densifying the full samples × genes matrix.
#' @param mat dgCMatrix genes × samples.
#' @param verbose Logical.
#' @return Dense Pearson correlation matrix (genes × genes).
#' @keywords internal
.compute_correlation_sparse_pearson <- function(mat, verbose = TRUE) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  n_genes <- nrow(mat)
  n_samples <- ncol(mat)

  if (n_genes > 50000L) {
    stop("[geneSCOPE::.compute_correlation] Sparse Pearson path is intended for moderate gene counts; reduce genes or increase memory_limit_gb.")
  }

  .log_info("computeCorrelation", "S04",
    "Computing Pearson correlation via sparse tcrossprod...", verbose
  )

  mu <- rowMeans(mat)
  xxT <- tcrossprod(mat)
  xxT_dense <- as.matrix(xxT)

  cov_mat <- (xxT_dense - n_samples * (mu %o% mu)) / max(1, n_samples - 1)
  cov_mat[!is.finite(cov_mat)] <- 0

  sd_vec <- sqrt(pmax(diag(cov_mat), 0))
  denom <- sd_vec %o% sd_vec
  cor_mat <- cov_mat
  pos <- denom > 0
  cor_mat[pos] <- cov_mat[pos] / denom[pos]
  cor_mat[!pos] <- 0
  diag(cor_mat) <- 1

  rownames(cor_mat) <- rownames(mat)
  colnames(cor_mat) <- rownames(mat)
  cor_mat
}

#' Extract Pearson correlations for specific edge pairs without densifying.
#' @description
#' Internal helper for `.get_pairwise_cor_for_edges`.
#' Computes only the requested pairwise correlations, keeping large Pearson matrices in sparse form.
#' @param scope_obj scope_object
#' @param grid_name character grid layer name
#' @param level "grid" or "cell" (default used by cluster code is "cell")
#' @param kept_genes character vector of genes (in the same order used by edge_i/j)
#' @param edge_i integer vectors (1-based) indexing into kept_genes
#' @param edge_j Parameter value.
#' @return numeric vector of correlations for each pair (edge_i[k], edge_j[k]).
#' @keywords internal
.get_pairwise_cor_for_edges <- function(scope_obj, grid_name, level = c("cell", "grid"),
                                        kept_genes, edge_i, edge_j) {
    level <- match.arg(level)
    rMat <- .get_pearson_matrix(scope_obj, grid_name = grid_name, level = level)
    if (!is.matrix(rMat)) {
        # Best-effort: try to coerce only if not file-backed; avoid huge subsetting
        rMat <- try(as.matrix(rMat), silent = TRUE)
        if (inherits(rMat, "try-error") || !is.matrix(rMat)) {
            stop("Pearson correlation matrix is not a base matrix; cannot extract pairs efficiently.")
        }
    }
    ridx <- match(kept_genes, rownames(rMat))
    if (anyNA(ridx)) {
        stop("Some kept_genes not present in Pearson matrix.")
    }
    rMat[cbind(ridx[edge_i], ridx[edge_j])]
}

#' Get Lee Matrix
#' @description
#' Internal helper for `.get_lee_matrix`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_layer Layer name.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.get_lee_matrix <- function(scope_obj, grid_name = NULL, lee_layer = NULL, verbose = FALSE) {
    ## ---- 0. Select grid sub-layer ------------------------------------------------
    g_layer <- .select_grid_layer(scope_obj, grid_name, verbose = verbose)
    if (is.null(grid_name)) { # Write back the actual name
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
    }

    ## ---- 1. Auto-detect LeeStats layer name ----------------------------------------
    if (is.null(lee_layer)) {
        cand <- character(0)
        if (!is.null(scope_obj@stats[[grid_name]])) {
            cand <- names(scope_obj@stats[[grid_name]])
        }
        cand <- c(cand, names(g_layer))
        cand <- unique(cand[grepl("^LeeStats_", cand)])

        if (length(cand) == 0L) {
            stop(
                "No layer starting with 'LeeStats_' found for grid '",
                grid_name, "'."
            )
        }
        if (length(cand) == 1L) {
            lee_layer <- cand
        } else if ("LeeStats_Xz" %in% cand) { # Prefer default naming
            lee_layer <- "LeeStats_Xz"
        } else {
            stop(
                "Multiple LeeStats layers detected (",
                paste(cand, collapse = ", "),
                "); please specify `lee_layer` explicitly."
            )
        }
    }

    ## ---- 2. Search @stats -> @grid in order ------------------------------------
    leeStat <- NULL
    if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_layer]])) {
        leeStat <- scope_obj@stats[[grid_name]][[lee_layer]]
    }

    if (is.null(leeStat) && !is.null(g_layer[[lee_layer]])) {
        leeStat <- g_layer[[lee_layer]]
    }

    if (is.null(leeStat) || is.null(leeStat$L)) {
        stop(
            "Layer '", lee_layer, "' in grid '", grid_name,
            "' does not contain a valid Lee's L matrix."
        )
    }

    Lmat <- leeStat$L

    ## ---- 3. If Pearson correlation matrix exists, take intersection ---------------------------------
    # Note: do NOT intersect Lee's L with Pearson gene set here.
    # This keeps the full L matrix as provided by the LeeStats layer
    # and avoids unintended shrinkage of the gene set during Stage-1.

    Lmat
}

#' Get Pearson Matrix
#' @description
#' Internal helper for `.get_pearson_matrix`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param level Parameter value.
#' @return Return value used internally.
#' @keywords internal
.get_pearson_matrix <- function(scope_obj,
                              grid_name = NULL,
                              level = c("grid", "cell")) {
    level <- match.arg(level)

    ## ---------- 1. Set layer name & final target ------------------------------------
    corr_name <- ".pearson_cor"
    f_cell_suf <- "_cell" # single-cell layer suffix

    ## ---------- 2. Get matrix --------------------------------------------------
    if (level == "grid") {
        g_layer <- .select_grid_layer(scope_obj, grid_name)
        if (is.null(grid_name)) {
            grid_name <- names(scope_obj@grid)[
                vapply(scope_obj@grid, identical, logical(1), g_layer)
            ]
        }

        ##   2a. New version: @stats[[grid_name]]
        rmat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
            !is.null(scope_obj@stats[[grid_name]][[corr_name]])) {
            scope_obj@stats[[grid_name]][[corr_name]]
        } else {
            NULL
        }

        ##   2b. Fallback: @grid[[grid_name]]
        if (is.null(rmat) && !is.null(g_layer[[corr_name]])) {
            rmat <- g_layer[[corr_name]]
        }

        if (is.null(rmat)) {
            stop("Pearson matrix not found for grid layer '", grid_name, "'.")
        }

        ## Lee's L alignment
        Lmat <- tryCatch(
            .get_lee_matrix(scope_obj, grid_name = grid_name),
            error = function(e) NULL
        )
        if (!is.null(Lmat)) {
            common <- intersect(rownames(rmat), rownames(Lmat))
            rmat <- rmat[common, common, drop = FALSE]
        }
    } else { # ---------- single-cell ----------

        ##   2a. New version: @stats[["cell"]]
        rmat <- if (!is.null(scope_obj@stats[["cell"]]) &&
            !is.null(scope_obj@stats[["cell"]][[paste0(corr_name, f_cell_suf)]])) {
            scope_obj@stats[["cell"]][[paste0(corr_name, f_cell_suf)]]
        } else {
            NULL
        }

        ##   2b. Fallback: @cells$.pearson_cor
        if (is.null(rmat) && !is.null(scope_obj@cells[[corr_name]])) {
            rmat <- scope_obj@cells[[corr_name]]
        }

        if (is.null(rmat)) {
            stop("Pearson matrix at cell level not found in scope_obj.")
        }

        ## Lee's L alignment (if single-cell layer exists)
        Lmat <- NULL
        if (!is.null(scope_obj@stats[["cell"]])) {
            layer_L <- intersect(
                grep("^LeeStats_", names(scope_obj@stats[["cell"]]),
                    value = TRUE
                ),
                paste0("LeeStats_Xz", f_cell_suf)
            )
            if (length(layer_L)) {
                Lmat <- scope_obj@stats[["cell"]][[layer_L[1]]]$L
            }
        }
        if (is.null(Lmat) && !is.null(scope_obj@cells$LeeStats_Xz)) {
            Lmat <- scope_obj@cells$LeeStats_Xz$L
        }

        if (!is.null(Lmat)) {
            common <- intersect(rownames(rmat), rownames(Lmat))
            rmat <- rmat[common, common, drop = FALSE]
        }
    }

    return(rmat)
}

#' Intelligent gene subset selection
#' @description Select gene subset intelligently based on dataset size and available memory
#' @param mat Expression matrix
#' @param max_genes Maximum number of genes
#' @param method Selection method: "variance", "mean", "dropout", "random"
#' @param verbose Logical.
#' @return Selected gene names vector
#' @keywords internal
.select_top_genes <- function(mat, max_genes = 50000, method = "variance",
                              verbose = TRUE) {
  if (ncol(mat) <= max_genes) {
    return(colnames(mat))
  }

  .log_info("computeCorrelation", "S03",
    paste0("Selecting top ", max_genes, " genes by ", method, "..."), verbose
  )

  if (method == "variance") {
    # Compute gene variance
    gene_vars <- colMeans(mat^2) - colMeans(mat)^2
    top_idx <- order(gene_vars, decreasing = TRUE)[1:max_genes]
  } else if (method == "mean") {
    # Select by mean expression
    gene_means <- colMeans(mat)
    top_idx <- order(gene_means, decreasing = TRUE)[1:max_genes]
  } else if (method == "dropout") {
    # Select by dropout rate (proportion of cells expressing)
    gene_dropout <- colMeans(mat > 0)
    top_idx <- order(gene_dropout, decreasing = TRUE)[1:max_genes]
  } else if (method == "random") {
    # Random selection
    top_idx <- sample(ncol(mat), max_genes)
  } else {
    stop("Unknown selection method: ", method)
  }

  return(colnames(mat)[top_idx])
}

#' Chunked conversion of sparse matrix to dense matrix
#' @description
#' Internal helper for `.sparse_to_dense_chunked`.
#' @param sparse_mat Sparse matrix
#' @param chunk_size Chunk size
#' @return Dense matrix
#' @keywords internal
.sparse_to_dense_chunked <- function(sparse_mat, chunk_size) {
  n_cols <- ncol(sparse_mat)
  dense_mat <- matrix(0, nrow(sparse_mat), n_cols)

  for (start in seq(1, n_cols, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n_cols)
    dense_mat[, start:end] <- as.matrix(sparse_mat[, start:end])
  }

  return(dense_mat)
}

#' Compute pairwise correlation matrices
#' @description
#' Internal helper for `.compute_correlation`.
#' @param scope_obj A `scope_object` containing expression layers.
#' @param level Correlation level (`cell` or `grid`).
#' @param grid_name Optional grid name when operating on grid-level data.
#' @param layer Expression layer name to correlate (e.g., `logCPM`).
#' @param method Correlation method (`pearson`, `spearman`, `kendall`).
#' @param ncores Number of threads to use for computation.
#' @param store_layer Target slot name for storing the correlation matrix.
#' @param compute_fdr Whether to compute FDR alongside correlations.
#' @param verbose Emit progress messages when TRUE.
#' @param blocksize Parameter value.
#' @param chunk_size Parameter value.
#' @param memory_limit_gb Parameter value.
#' @param use_bigmemory Logical flag.
#' @param backing_path Filesystem path.
#' @param force_compute Parameter value.
#' @param fdr_method Parameter value.
#' @return The updated scope object (invisibly).
#' @keywords internal
.compute_correlation <- function(scope_obj,
                               level = c("cell", "grid"),
                               grid_name = NULL,
                               layer = "logCPM",
                               method = c("pearson", "spearman", "kendall"),
                               blocksize = 2000,
                               ncores = 16,
                               chunk_size = 1000,
                               memory_limit_gb = 16,
                               use_bigmemory = TRUE,
                               backing_path = tempdir(),
                               force_compute = FALSE,
                               store_layer = ".pearson_cor",
                               compute_fdr = TRUE,
                               fdr_method = "BH",
                               verbose = TRUE) {
  level <- match.arg(level)
  method <- match.arg(method)

  step_s01 <- .log_step("computeCorrelation", "S01", "resolve runtime config", verbose)
  summary_parts <- c(
    paste0("level=", level),
    if (!is.null(grid_name)) paste0("grid=", grid_name) else NULL,
    paste0("layer=", layer),
    paste0("method=", method),
    paste0("ncores=", ncores),
    paste0("memory_limit_gb=", memory_limit_gb),
    paste0("use_bigmemory=", use_bigmemory),
    paste0("compute_fdr=", compute_fdr)
  )
  step_s01$enter(extra = paste(summary_parts, collapse = " "))

  runtime_cfg <- .compute_correlation_resolve_runtime_config(ncores, verbose)
  ncores_safe <- runtime_cfg$ncores_safe
  thread_source <- runtime_cfg$thread_source
  old_blas_env <- runtime_cfg$old_blas_env
  old_mkl_env <- runtime_cfg$old_mkl_env

  on.exit(.compute_correlation_restore_env(old_blas_env, old_mkl_env))

  inputs <- .compute_correlation_validate_inputs(
    scope_obj = scope_obj,
    level = level,
    grid_name = grid_name,
    layer = layer,
    verbose = verbose
  )
  expr_mat <- inputs$expr_mat
  grid_name <- inputs$grid_name
  sample_label <- inputs$sample_label

  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)

  step_s01$done(extra = paste0("n_genes=", n_genes, " n_samples=", n_samples))

  ## ---- Memory check and chunking strategy -------------------------------
  step_s02 <- .log_step("computeCorrelation", "S02", "memory guard and mode choice", verbose)
  step_s02$enter(extra = paste0("n_genes=", n_genes, " n_samples=", n_samples, " sample=", sample_label))
  .log_info("computeCorrelation", "S02",
    paste0("Data dimensions: ", n_genes, " genes x ", n_samples, " ", sample_label),
    verbose
  )
  # Estimate correlation matrix size (n_genes x n_genes x 8 bytes)
  cor_matrix_gb <- (n_genes^2 * 8) / (1024^3)

  # Memory guard: assume each thread needs at least single-core footprint
  sys_mem_gb <- .get_system_memory_gb()
  est_total_gb <- cor_matrix_gb * ncores_safe
  if (est_total_gb > sys_mem_gb) {
    stop(
      "[geneSCOPE::.compute_correlation] Estimated memory requirement (",
      round(est_total_gb, 1), " GB) exceeds system capacity (",
      round(sys_mem_gb, 1), " GB). Reduce ncores or gene count."
    )
  }

  .log_info("computeCorrelation", "S02",
    paste0("Estimated correlation matrix size: ", round(cor_matrix_gb, 2), " GB"),
    verbose
  )

  cor_matrix <- NULL
  bigmemory_mode <- FALSE
  use_sparse_path <- FALSE
  expr_centered <- NULL
  if (cor_matrix_gb > memory_limit_gb) {
    if (!use_bigmemory) {
      stop(
        "Correlation matrix too large (", round(cor_matrix_gb, 1), " GB) exceeds limit (",
        memory_limit_gb, " GB). Please reduce gene number or increase memory limit or set use_bigmemory=TRUE."
      )
    }
    bigmemory_mode <- TRUE
    .log_backend("computeCorrelation", "S02", "mode", "bigmemory",
      reason = "cor_matrix_gb>memory_limit_gb", verbose = verbose
    )
    .log_info("computeCorrelation", "S02",
      "Dense correlation exceeds memory_limit_gb; switching to bigmemory-backed chunked computation.",
      verbose
    )
  } else {
    .log_backend("computeCorrelation", "S02", "mode", "in_memory",
      reason = "within_limit", verbose = verbose
    )
  }
  step_s02$done(extra = paste0("mode=", ifelse(bigmemory_mode, "bigmemory", "in_memory")))

  ## ---- Data preprocessing --------------------------------------
  step_s03 <- .log_step("computeCorrelation", "S03", "preprocess expression matrix", verbose)
  if (bigmemory_mode) {
    step_s03$enter(extra = "skip (bigmemory mode)")
    .log_info("computeCorrelation", "S03", "Skipping preprocessing in bigmemory mode.", verbose)
    step_s03$done(extra = "skip")
  } else {
    step_s03$enter()
    .log_info("computeCorrelation", "S03", "Preprocessing expression matrix...", verbose)

    # If densifying a huge sparse matrix would exceed memory, use sparse Pearson path.
    dense_input_gb <- (as.double(n_genes) * as.double(n_samples) * 8) / (1024^3)
    densify_risky <- inherits(expr_mat, "dgCMatrix") &&
      identical(method, "pearson") &&
      dense_input_gb * 2 > min(memory_limit_gb, sys_mem_gb * 0.8)

    if (densify_risky) {
      .log_info("computeCorrelation", "S03",
        paste0("Input densification estimated at ~", round(dense_input_gb, 1), " GB; using sparse Pearson path."),
        verbose
      )
      .log_backend("computeCorrelation", "S03", "R",
        paste0("sparse_tcrossprod matrix=dgCMatrix spmm=FALSE center=implicit scale=FALSE densify_gb=", round(dense_input_gb, 1)),
        verbose = verbose
      )
      use_sparse_path <- TRUE
    } else {
      .log_backend("computeCorrelation", "S03", "native",
        "pearson_cor matrix=dense spmm=FALSE center=TRUE scale=FALSE",
        verbose = verbose
      )
      # Transpose matrix to cells x genes (.pearson_cor expects this)
      expr_dense <- as.matrix(t(expr_mat))

      # Center data (.pearson_cor expects centered data)
      .log_info("computeCorrelation", "S03", "Centering data...", verbose)
      expr_centered <- scale(expr_dense, center = TRUE, scale = FALSE)

      # Remove NaN/Inf
      expr_centered[!is.finite(expr_centered)] <- 0

      .log_info("computeCorrelation", "S03", "Data preprocessing completed", verbose)
    }
    step_s03$done(extra = if (use_sparse_path) "path=sparse_tcrossprod" else "path=dense")
  }

  ## ---- Call C++ correlation function --------------------------------
  step_s04 <- .log_step("computeCorrelation", "S04", "compute correlation matrix", verbose)
  step_s04$enter(extra = paste0("blocksize=", blocksize, " threads=", ncores_safe))
  if (bigmemory_mode) {
    cor_matrix <- .compute_correlation_big_matrix(
      mat = expr_mat,
      method = method,
      chunk_size = chunk_size,
      ncores = ncores_safe,
      force_compute = force_compute,
      memory_limit_gb = memory_limit_gb,
      backing_path = backing_path,
      verbose = verbose,
      thread_source = thread_source
    )
  } else if (use_sparse_path) {
    .log_backend("computeCorrelation", "S04", "R",
      "sparse_tcrossprod matrix=dgCMatrix threads=1 thread_source=single spmm=FALSE",
      verbose = verbose
    )
    cor_matrix <- .compute_correlation_sparse_pearson(
      expr_mat,
      verbose = verbose
    )
  } else {
    .log_info("computeCorrelation", "S04",
      paste0("Starting correlation matrix computation (using ", ncores_safe, " threads)..."),
      verbose
    )
    .log_backend("computeCorrelation", "S04", "native",
      paste0(
        "pearson_cor threads=", ncores_safe,
        " thread_source=", thread_source,
        " blocksize=", blocksize
      ),
      verbose = verbose
    )

    start_time <- Sys.time()

    # Use safe thread count to call C++ function
    cor_matrix <- tryCatch(
      {
        .pearson_cor(
          X = expr_centered,
          bs = blocksize,
          n_threads = ncores_safe
        )
      },
      error = function(e) {
        # If error, try single thread
        .log_backend("computeCorrelation", "S04", "native",
          "pearson_cor threads=1 thread_source=fallback",
          reason = "fallback_single_thread", verbose = verbose
        )
        .log_info("computeCorrelation", "S04",
          paste0("!!! Warning: Multithreaded computation failed, trying single thread: ", e$message, " !!!"),
          verbose
        )
        .pearson_cor(
          X = expr_centered,
          bs = blocksize,
          n_threads = 1
        )
      }
    )

    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

    .log_info("computeCorrelation", "S04",
      paste0("Correlation matrix computation completed, time elapsed: ", round(elapsed, 2), " minutes"),
      verbose
    )
  }
  step_s04$done(extra = paste0("mode=", ifelse(bigmemory_mode, "bigmemory", "in_memory")))

  ## ---- Compute FDR (if needed) ----------------------------------
  fdr_matrix <- NULL
  step_s05 <- .log_step("computeCorrelation", "S05", "compute FDR (optional)", verbose)
  step_s05$enter(extra = paste0("requested=", compute_fdr))
  if (!compute_fdr) {
    .log_backend("computeCorrelation", "S05", "fdr", "skipped",
      reason = "disabled", verbose = verbose
    )
    step_s05$done(extra = "skipped")
  } else if (bigmemory_mode) {
    compute_fdr <- FALSE
    .log_backend("computeCorrelation", "S05", "fdr", "skipped",
      reason = "bigmemory", verbose = verbose
    )
    .log_info("computeCorrelation", "S05",
      "Skipping FDR computation in bigmemory mode (file-backed matrix).",
      verbose
    )
    step_s05$done(extra = "skipped")
  } else {
    .log_info("computeCorrelation", "S05", "Computing FDR...", verbose)

    # Avoid large conversion from sparse to dense matrix
    # Directly extract upper triangle values from correlation matrix, avoid as.vector()
    upper_tri <- upper.tri(cor_matrix, diag = FALSE)

    # Use matrix indexing instead of converting entire matrix
    upper_indices <- which(upper_tri, arr.ind = TRUE)
    cor_values <- cor_matrix[upper_indices]

    # Remove NaN and Inf values
    valid_mask <- is.finite(cor_values) & !is.na(cor_values)
    cor_values_clean <- cor_values[valid_mask]
    valid_indices <- upper_indices[valid_mask, , drop = FALSE]

    if (length(cor_values_clean) == 0) {
      .log_info("computeCorrelation", "S05",
        "!!! Warning: All correlation coefficients are NaN or Inf, cannot compute FDR !!!",
        verbose
      )
      fdr_matrix <- matrix(1, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
    } else {
      # Compute two-tailed p-values (only for valid values)
      n_samples <- n_samples
      df <- n_samples - 2

      if (df > 0) {
        # Compute t-statistic (avoid division by zero)
        cor_values_safe <- pmax(pmin(cor_values_clean, 0.9999), -0.9999)
        t_stat <- cor_values_safe * sqrt(df / (1 - cor_values_safe^2))

        # Compute two-tailed p-values
        p_values <- 2 * pt(-abs(t_stat), df = df)

        # Apply FDR correction
        # Reconstruct symmetric matrix (initialize to 1)
        fdr_matrix <- matrix(1, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))

        # Fill only valid positions
        fdr_matrix[valid_indices] <- p.adjust(p_values, method = fdr_method)
        fdr_matrix[valid_indices[, c(2, 1)]] <- p.adjust(p_values, method = fdr_method) # Symmetric fill
        diag(fdr_matrix) <- 0 # Set diagonal to 0

        # Set row and column names
        rownames(fdr_matrix) <- rownames(cor_matrix)
        colnames(fdr_matrix) <- colnames(cor_matrix)

        .log_info("computeCorrelation", "S05",
          paste0("FDR computation completed, using method: ", fdr_method),
          verbose
        )
        .log_info("computeCorrelation", "S05",
          paste0("Valid correlation coefficients: ", length(cor_values_clean), " / ", length(cor_values)),
          verbose
        )
      } else {
        .log_info("computeCorrelation", "S05",
          "!!! Warning: Insufficient sample size, cannot compute FDR !!!",
          verbose
        )
        fdr_matrix <- matrix(1, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
      }
    }
    step_s05$done(extra = paste0("method=", fdr_method))
  }

  ## ---- Set gene names ------------------------------------
  gene_names <- rownames(expr_mat)
  rownames(cor_matrix) <- gene_names
  colnames(cor_matrix) <- gene_names

  step_s06 <- .log_step("computeCorrelation", "S06", "store outputs", verbose)
  step_s06$enter(extra = paste0("store_layer=", store_layer))
  scope_obj <- .compute_correlation_assemble_outputs(
    scope_obj = scope_obj,
    level = level,
    grid_name = grid_name,
    store_layer = store_layer,
    cor_matrix = cor_matrix,
    fdr_matrix = fdr_matrix,
    verbose = verbose
  )
  step_s06$done(extra = paste0("fdr=", !is.null(fdr_matrix)))

  invisible(scope_obj)
}

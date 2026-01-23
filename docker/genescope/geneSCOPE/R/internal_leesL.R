#' Add Lee's L statistics to a scope_object
#' @description
#'   High-level wrapper that computes Lee's L, analytical Z-scores, empirical
#'   p-values via block permutations, FDR, spatial gradients, quality-control
#'   metrics, and stores everything under a new layer in \code{@grid}.
#' @param scope_obj A \code{scope_object} with at least one populated \code{@grid} slot.
#' @param grid_name Character. Name of the grid sub-layer to process. If
#'   \code{NULL} and only one sub-layer exists, it is selected automatically.
#' @param genes Optional character vector of genes to include; if \code{NULL}
#'   all genes are used.
#' @param within Logical. If \code{TRUE} restrict analysis to the selected
#'   gene set on both axes (default). Otherwise compute gene × all.
#' @param ncores Integer. Number of cores for parallel processing (default 1).
#' @param block_side Integer. Number of grid cells per side for block partitioning (default 8).
#' @param perms Integer. Number of permutations for Monte-Carlo p-values (default 1000).
#' @param block_size Integer. Number of permutations processed per batch (default 64).
#' @param L_min Numeric threshold used when building the QC similarity graph.
#' @param norm_layer Character. Name of the normalised expression layer (default "Xz").
#' @param lee_stats_layer_name Character. Name for the output statistics layer.
#' @param legacy_formula Logical. Use legacy denominator for compatibility.
#' @param mem_limit_GB Numeric. RAM threshold that triggers streaming mode (default 8).
#' @param chunk_size Integer. Number of columns processed per chunk (default 128).
#' @param use_bigmemory Logical. Use file-backed matrices for large computations (default TRUE).
#' @param backing_path Character. Directory for temporary files (default tempdir()).
#' @param cache_inputs Logical. Cache preprocessed X/Z/W and block IDs for reuse across calls (default TRUE).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @param ncore Deprecated. Use \code{ncores} instead.
#' @return The modified \code{scope_object}.
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom stats sd pnorm p.adjust coef lm
#' @importFrom igraph graph_from_adjacency_matrix simplify degree cluster_louvain cluster_leiden components modularity
#' @keywords internal
.compute_l <- function(scope_obj,
                     grid_name = NULL,
                     genes = NULL,
                     within = TRUE,
                     ncores = 1,
                     block_side = 8,
                     perms = 1000,
                     block_size = 64,
                     L_min = 0,
                     norm_layer = "Xz",
                     lee_stats_layer_name = NULL,
                     legacy_formula = FALSE,
                     mem_limit_GB = 2,
                     chunk_size = 32L,
                     use_bigmemory = TRUE,
                     backing_path = tempdir(),
                     cache_inputs = TRUE,
                     verbose = TRUE,
                     ncore = NULL) {
    parent <- "computeL"

    ## --- 0. Thread-safe preprocessing: automatic thread management and error recovery ---

    # Handle deprecated parameter
    if (!is.null(ncore)) {
        warning("'ncore' is deprecated, please use 'ncores' instead. Will use its value this time.",
            call. = FALSE, immediate. = TRUE
        )
        ncores <- ncore
    }

    user_set_bigmemory <- !missing(use_bigmemory)
    user_set_chunk <- !missing(chunk_size)

    # Get system information and aggressive thread count (use all visible cores up to request)
    os_type <- .detect_os()
    avail_cores <- max(1L, detectCores(logical = TRUE))
    ncores <- max(1L, min(ncores, avail_cores))
    grid_label <- if (is.null(grid_name)) "auto" else as.character(grid_name)[1]
    step01 <- .log_step(parent, "S01", "resolve inputs and weights", verbose)
    step01$enter(paste0("grid_name=", grid_label, " ncores=", ncores))

    # Use thread configuration for mixed OpenMP/BLAS operations
    thread_config <- .configure_threads_for("mixed", ncores, restore_after = TRUE)
    on.exit({
        restore_fn <- attr(thread_config, "restore_function")
        if (!is.null(restore_fn)) restore_fn()
    })

    # Use OpenMP threads for C++ operations
    ncores_cpp <- thread_config$openmp_threads
    .log_info(parent, "S01", paste0(
        "ncores_use=", ncores, "/", avail_cores,
        " openmp_threads=", ncores_cpp,
        " os=", os_type
    ), verbose)

    ## --- 1. Check if spatial weights exist, auto-compute if missing ---
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else {
        grid_name
    }

    # Check if spatial weight matrix exists
    if (is.null(g_layer$W) || sum(g_layer$W) == 0) {
        .log_info(parent, "S01", "W missing; running computeWeights(style=B)", verbose)
        .log_backend(parent, "S01", "weights", "computeWeights(style=B)", verbose = verbose)
        scope_obj <- .compute_weights(scope_obj,
            grid_name = grid_name,
            style = "B",
            store_mat = TRUE,
            store_listw = FALSE,
            verbose = verbose
        )
        # Re-extract layer after computing weights
        g_layer <- .select_grid_layer(scope_obj, grid_name)
    }

    step01$done(paste0("grid_name=", grid_name))

    ## --- 2. Memory guard and mode selection ---
    step02 <- .log_step(parent, "S02", "memory guard and mode selection", verbose)
    step02$enter(paste0(
        "mem_limit_GB=", mem_limit_GB,
        " use_bigmemory=", use_bigmemory,
        " chunk_size=", chunk_size
    ))
    bigmemory_available <- requireNamespace("bigmemory", quietly = TRUE)
    bigmemory_forced_off <- use_bigmemory && !bigmemory_available
    if (bigmemory_forced_off) {
        if (verbose) {
            .log_info(parent, "S02", "bigmemory not available; forcing use_bigmemory=FALSE", verbose)
        }
        use_bigmemory <- FALSE
    }

    # Cross-platform memory calculation
    all_genes <- if (!is.null(g_layer[[norm_layer]])) {
        ncol(g_layer[[norm_layer]])
    } else {
        stop("Normalized layer '", norm_layer, "' not found")
    }

    n_genes_use <- if (!is.null(genes)) length(genes) else all_genes
    matrix_size_gb <- (n_genes_use^2 * 8) / (1024^3)

    # Memory guard: estimate per-thread usage as single-core size × threads
    sys_mem_gb <- .get_system_memory_gb()
    est_total_gb <- matrix_size_gb * ncores
    .log_info(parent, "S02", paste0(
        "matrix_size_gb=", round(matrix_size_gb, 1),
        " est_total_gb=", round(est_total_gb, 1),
        " sys_mem_gb=", round(sys_mem_gb, 1)
    ), verbose)
    if (est_total_gb > sys_mem_gb) {
        stop(
            "[geneSCOPE::.compute_l] Estimated memory requirement (",
            round(est_total_gb, 1), " GB) exceeds system capacity (",
            round(sys_mem_gb, 1), " GB). Reduce ncores or gene set size."
        )
    }

    mem_reason <- NULL
    # Respect user bigmemory/chunk overrides while still hinting when streaming is advisable
    if (matrix_size_gb > mem_limit_GB) {
        if (use_bigmemory) {
            mem_reason <- "matrix_size_gb>mem_limit_gb"
            if (verbose) {
                .log_info(parent, "S02", paste0(
                    "large matrix detected (", round(matrix_size_gb, 1),
                    " GB); staying in bigmemory/streaming mode.",
                    if (!user_set_chunk) " You may tune chunk_size to trade IO vs RAM." else ""
                ), verbose)
            }
        } else if (user_set_bigmemory && !bigmemory_forced_off) {
            mem_reason <- "user_forced_inmemory"
            if (verbose) .log_info(parent, "S02", paste0(
                "Warning: requested use_bigmemory=FALSE with large matrix (",
                round(matrix_size_gb, 1), " GB); proceeding in-memory as requested."
            ), verbose)
        } else if (!bigmemory_available) {
            mem_reason <- "bigmemory_unavailable"
            if (verbose) .log_info(parent, "S02",
                "!!! Warning: bigmemory not available; using regular matrices !!!",
                verbose
            )
            use_bigmemory <- FALSE
        } else {
            mem_reason <- "auto_enable"
            use_bigmemory <- TRUE
            if (verbose) .log_info(parent, "S02", paste0(
                "large matrix detected (", round(matrix_size_gb, 1),
                " GB); enabling bigmemory/streaming."
            ), verbose)
        }
    } else {
        mem_reason <- "within_limit"
    }

    mem_mode <- if (use_bigmemory) "bigmemory" else "inmemory"
    .log_backend(parent, "S02", "mem_mode", mem_mode, reason = mem_reason, verbose = verbose)

    # Cross-platform temporary directory
    if (use_bigmemory) {
        if (backing_path == tempdir()) {
            # Use platform-specific temporary directory
            backing_path <- switch(os_type,
                windows = Sys.getenv("TEMP", tempdir()),
                linux = Sys.getenv("TMPDIR", "/tmp"),
                macos = Sys.getenv("TMPDIR", "/tmp")
            )
        }

        if (!dir.exists(backing_path)) {
            dir.create(backing_path, recursive = TRUE)
        }
    }

    step02$done(paste0("mem_mode=", mem_mode))

    ## --- 3. Lee's L computation with error recovery ---
    step03 <- .log_step(parent, "S03", "compute Lee's L", verbose)
    step03$enter(paste0("mem_mode=", mem_mode, " ncores=", ncores))
    .log_backend(parent, "S03", "lee_L", paste0("C++ lee_L threads=", ncores_cpp), verbose = verbose)

    # Multi-retry mechanism with reduced thread count
    current_cores <- ncores
    min_cores <- 1
    success <- FALSE
    attempt <- 1

    while (!success && current_cores >= min_cores) {
        if (verbose && attempt > 1) {
            .log_info(parent, "S03", paste0("retry #", attempt, " with ", current_cores, " cores"), verbose)
            .log_backend(parent, "S03", "lee_L",
                paste0("C++ lee_L threads=", current_cores),
                reason = "retry",
                verbose = verbose
            )
        }

        result <- tryCatch(
            {
                # Execute computation
                t_start <- Sys.time()
                res <- .compute_lee_l(scope_obj,
                    grid_name = grid_name,
                    norm_layer = norm_layer,
                    genes = genes,
                    within = within,
                    ncores = current_cores, # Use current core count
                    mem_limit_GB = mem_limit_GB,
                    chunk_size = chunk_size,
                    use_bigmemory = use_bigmemory,
                    backing_path = backing_path,
                    block_side = block_side,
                    cache_inputs = cache_inputs,
                    input_cache = if (cache_inputs && !is.null(scope_obj@stats[[grid_name]])) {
                        scope_obj@stats[[grid_name]][["_lee_input_cache"]]
                    } else {
                        NULL
                    }
                )
                t_end <- Sys.time()

                if (verbose) {
                    time_msg <- if (attempt == 1) "completed" else "retry successful"
                    .log_info(parent, "S03", paste0(
                        "Lee's L ", time_msg, " (", format(t_end - t_start), ")"
                    ), verbose)
                }

                list(success = TRUE, object = res)
            },
            error = function(e) {
                if (verbose && attempt > 1) {
                    .log_info(parent, "S03", paste0("attempt failed: ", conditionMessage(e)), verbose)
                }
                list(success = FALSE, error = e)
            }
        )

        # Update status
        success <- result$success

        if (success) {
            res <- result$object
        } else {
            attempt <- attempt + 1
            # Reduce thread count
            current_cores <- max(floor(current_cores / 2), min_cores)
            if (verbose) {
                .log_info(parent, "S03", paste0("reducing cores to ", current_cores, " and retrying"), verbose)
                .log_backend(parent, "S03", "lee_L",
                    paste0("C++ lee_L threads=", current_cores),
                    reason = "retry",
                    verbose = verbose
                )
            }
            # Give system some time to recover
            Sys.sleep(3)
            # Force garbage collection
            gc(verbose = FALSE)
        }
    }

    if (!success) {
        stop("Unable to compute Lee's L statistics even with minimal thread count")
    }

    step03$done(paste0("threads_final=", current_cores))

    L <- res$Lmat
    X_full <- res$X_full
    X_used <- res$X_used
    W <- res$W
    grid_inf <- res$grid_info
    gname <- res$grid_name
    n <- nrow(X_full)

    ## --- 4. Analytical Z computation with memory optimization ---
    step04 <- .log_step(parent, "S04", "compute Z-scores", verbose)
    step04$enter(paste0("within=", within, " mem_mode=", mem_mode))
    if (within || is.null(genes)) {
        S0 <- sum(W)
        EZ <- -1 / (n - 1)
        Var <- (n^2 * (n - 2)) / ((n - 1)^2 * (n - 3) * S0)

        # Handle large matrices more carefully
        if (inherits(L, "big.matrix")) {
            if (verbose) {
                .log_info(parent, "S04", "computing Z-scores in chunks", verbose)
            }
            Z_mat <- L # Use the same backing store
            # Process in chunks to avoid memory issues
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

    step04$done()

    block_id <- res$block_id
    input_cache <- res$input_cache

    ## --- 5. Monte Carlo p-values with BLAS control and error recovery ---
    step05 <- .log_step(parent, "S05", "permutation p-values", verbose)
    step05$enter(paste0("perms=", perms, " block_size=", block_size))
    P <- if (perms > 0) {
        # Temporarily disable BLAS threads for permutation tests
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
            RhpcBLASctl::blas_set_num_threads(1)
        }

        if (inherits(L, "big.matrix")) {
            if (verbose) .log_info(parent, "S05", "converting to regular matrix for permutation tests", verbose)
            L_reg <- as.matrix(L)
        } else {
            L_reg <- L
        }

        # Permutation test with error recovery
        if (verbose) .log_info(parent, "S05", "running permutation tests", verbose)
        perm_success <- FALSE
        perm_cores <- current_cores
        perm_attempt <- 1
        .log_backend(parent, "S05", "permutation",
            paste0("C++ lee_l_perm_block threads=", perm_cores),
            verbose = verbose
        )

        while (!perm_success && perm_cores >= 1) {
            if (verbose && perm_attempt > 1) {
                .log_info(parent, "S05", paste0("retry #", perm_attempt, " with ", perm_cores, " cores"), verbose)
                .log_backend(parent, "S05", "permutation",
                    paste0("C++ lee_l_perm_block threads=", perm_cores),
                    reason = "retry",
                    verbose = verbose
                )
            }

            # Sanity check to avoid C++ OOB when a gene subset is used
            if (ncol(X_used) != nrow(L_reg) || ncol(X_used) != ncol(L_reg)) {
                stop(sprintf("Permutation input mismatch: ncol(X_used)=%d but L_ref is %dx%d",
                    ncol(X_used), nrow(L_reg), ncol(L_reg)))
            }

            perm_result <- tryCatch(
                {
                    t_start <- Sys.time()
                    # Use X_used to match the dimensions of L_reg when a gene subset is used
                    p_result <- .lee_l_perm_block(X_used, W, L_reg,
                        block_id   = block_id,
                        perms      = perms,
                        block_size = block_size,
                        n_threads  = perm_cores
                    )
                    t_end <- Sys.time()
                    if (verbose) {
                        time_msg <- if (perm_attempt == 1) "completed" else "retry successful"
                        .log_info(parent, "S05", paste0(
                            "permutation test ", time_msg, " (", format(t_end - t_start), ")"
                        ), verbose)
                    }
                    list(success = TRUE, object = p_result)
                },
                error = function(e) {
                    if (verbose && perm_attempt > 1) {
                        .log_info(parent, "S05", paste0("test failed: ", conditionMessage(e)), verbose)
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
                    .log_info(parent, "S05", paste0("reducing cores to ", perm_cores, " and retrying"), verbose)
                    .log_backend(parent, "S05", "permutation",
                        paste0("C++ lee_l_perm_block threads=", perm_cores),
                        reason = "retry",
                        verbose = verbose
                    )
                }
                Sys.sleep(2)
                gc(verbose = FALSE)
            }
        }

        if (!perm_success) {
            if (verbose) .log_info(parent, "S05",
                "!!! Warning: Permutation test failed; setting P = NULL !!!",
                verbose
            )
            .log_backend(parent, "S05", "permutation", "skipped", reason = "failed", verbose = verbose)
            NULL
        } else {
            P
        }
    } else {
        .log_backend(parent, "S05", "permutation", "skipped", reason = "perms<=0", verbose = verbose)
        NULL
    }
    step05$done()

    ## --- 6. FDR: Three Monte Carlo smoothing strategies as in .get_top_l_vs_r ---
    step06 <- .log_step(parent, "S06", "FDR and q-values", verbose)
    step06$enter(paste0("mem_mode=", mem_mode, " perms=", perms))
    # Explanation:
    #   If P exists: P = (k+1)/(perms+1); infer k = round(P*(perms+1)-1)
    #   p_beta    = (k+1)/(perms+2)        (Jeffreys, default more robust/slightly conservative)
    #   p_mid     = (k+0.5)/(perms+1)      (mid-p, less conservative)
    #   p_uniform = (k+U)/(perms+1)        (random jitter, expected equal to discrete grid, breaks staircase)
    #   Main output FDR = BH(p_beta)
    if (inherits(Z_mat, "big.matrix")) {
        if (verbose) .log_info(parent, "S06", "computing FDR corrections", verbose)
        n_genes <- ncol(Z_mat)
        idx_chunks <- seq(1, n_genes, by = chunk_size)

        FDR_disc <- bigmemory::big.matrix(
            nrow = nrow(Z_mat), ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_beta <- bigmemory::big.matrix(
            nrow = nrow(Z_mat), ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_mid <- bigmemory::big.matrix(
            nrow = nrow(Z_mat), ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_uniform <- bigmemory::big.matrix(
            nrow = nrow(Z_mat), ncol = n_genes,
            type = "double", init = NA_real_
        )

        for (start in idx_chunks) {
            end <- min(start + chunk_size - 1L, n_genes)
            if (!is.null(P)) {
                p_disc <- P[, start:end]
                k_est <- round(p_disc * (perms + 1) - 1)
                k_est[k_est < 0] <- 0
                k_est[k_est > perms] <- perms
                p_beta <- (k_est + 1) / (perms + 2)
                p_mid <- (k_est + 0.5) / (perms + 1)
                p_uniform <- (k_est + matrix(runif(length(k_est)), nrow = nrow(k_est))) / (perms + 1)
            } else {
                # No permutation → analytical p
                Z_chunk <- Z_mat[, start:end]
                p_disc <- 2 * pnorm(-abs(Z_chunk))
                p_beta <- p_mid <- p_uniform <- p_disc
            }
            FDR_disc[, start:end] <- matrix(p.adjust(c(p_disc), "BH"), nrow = nrow(p_disc))
            FDR_beta[, start:end] <- matrix(p.adjust(c(p_beta), "BH"), nrow = nrow(p_beta))
            FDR_mid[, start:end] <- matrix(p.adjust(c(p_mid), "BH"), nrow = nrow(p_mid))
            FDR_uniform[, start:end] <- matrix(p.adjust(c(p_uniform), "BH"), nrow = nrow(p_uniform))
        }
        FDR_out_disc <- FDR_disc
        FDR_out_beta <- FDR_beta
        FDR_out_mid <- FDR_mid
        FDR_out_uniform <- FDR_uniform
    } else {
        # In-memory mode: compute in one go (exact BH)
        p_disc <- if (!is.null(P)) P else 2 * pnorm(-abs(Z_mat))
        if (!is.null(P)) {
            k_est <- round(p_disc * (perms + 1) - 1)
            k_est[k_est < 0] <- 0
            k_est[k_est > perms] <- perms
            p_beta <- (k_est + 1) / (perms + 2)
            p_mid <- (k_est + 0.5) / (perms + 1)
            p_uniform <- (k_est + runif(length(k_est))) / (perms + 1)
        } else {
            p_beta <- p_mid <- p_uniform <- p_disc
        }
        # BH on flattened vectors, then reshape
        shp <- dim(p_disc)
        FDR_out_disc <- matrix(p.adjust(p_disc, "BH"), nrow = shp[1], dimnames = dimnames(p_disc))
        FDR_out_beta <- matrix(p.adjust(p_beta, "BH"), nrow = shp[1], dimnames = dimnames(p_disc))
        FDR_out_mid <- matrix(p.adjust(p_mid, "BH"), nrow = shp[1], dimnames = dimnames(p_disc))
        FDR_out_uniform <- matrix(p.adjust(p_uniform, "BH"), nrow = shp[1], dimnames = dimnames(p_disc))
    }

    ## --- 5.x Extra: Storey q-value main FDR (in-memory mode only, big.matrix fallback) ---
    FDR_storey <- NULL
    pi0_hat <- NA_real_
    FDR_main_method <- if (inherits(Z_mat, "big.matrix")) "BH(beta p) (chunked, no q-value)" else "Storey q-value (beta p)"
    if (!inherits(Z_mat, "big.matrix")) {
        # Flatten to vector
        p_vec <- as.numeric(if (!is.null(P)) {
            # p_beta already computed above
            p_beta
        } else {
            p_disc
        })
        m_tot <- length(p_vec)
        if (m_tot > 0) {
            # 1) Estimate pi0
            lambdas <- seq(0.5, 0.95, by = 0.05)
            pi0_vals <- sapply(lambdas, function(lam) {
                mean(p_vec > lam) / (1 - lam)
            })
            pi0_hat <- min(1, min(pi0_vals, na.rm = TRUE))
            if (!is.finite(pi0_hat) || pi0_hat <= 0) pi0_hat <- 1
            # 2) Compute raw q
            o <- order(p_vec, na.last = NA)
            ro <- integer(m_tot)
            ro[o] <- seq_along(o)
            q_raw <- pi0_hat * m_tot * p_vec / pmax(ro, 1)
            # 3) Monotonic decreasing adjustment (tail to head)
            q_ord <- q_raw[o]
            for (i in seq_along(q_ord)[(length(q_ord) - 1):1]) {
                if (q_ord[i] > q_ord[i + 1]) q_ord[i] <- q_ord[i + 1]
            }
            q_final <- numeric(m_tot)
            q_final[o] <- pmin(q_ord, 1)
            # 4) Restore matrix
            FDR_storey <- matrix(q_final, nrow = nrow(FDR_out_beta), dimnames = dimnames(FDR_out_beta))
        }
    }
    # Main FDR choice: use FDR_storey if available, otherwise fallback to FDR_out_beta
    FDR_main <- if (!is.null(FDR_storey)) FDR_storey else FDR_out_beta
    n_sig_005 <- sum(FDR_main < 0.05, na.rm = TRUE)
    min_p_possible <- if (!is.null(P)) 1 / (perms + 1) else NA_real_

    step06$done(paste0("FDR_main_method=", FDR_main_method))

    ## --- 7. gradients and QC metrics ---
    step07 <- .log_step(parent, "S07", "gradients and QC metrics", verbose)
    step07$enter(paste0("within=", within, " L_min=", L_min))
    .log_backend(parent, "S07", "qc_graph", "igraph louvain (fallback=leiden)", verbose = verbose)

    ## --- 7. βx / βy ---
    centres <- with(
        grid_inf[match(res$cells, grid_inf$grid_id), ],
        data.frame(
            x = (xmin + xmax) / 2,
            y = (ymin + ymax) / 2
        )
    )
    betas <- t(apply(X_full, 2, function(v) {
        c(beta_x = coef(lm(v ~ centres$x + centres$y))[2:3])
    }))
    if (!is.null(genes)) betas <- betas[genes, , drop = FALSE]

    ## --- 7. QC computation ---
    qc <- NULL
    if (within || is.null(genes)) {
        # Convert to regular matrix for QC if needed
        L_qc <- if (inherits(L, "big.matrix")) {
            .log_info(parent, "S07", "converting subset of big.matrix for QC computation", verbose)
            # Sample subset for QC to avoid memory issues
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
                    # Sample for significance fraction
                    p_sample <- as.matrix(P[sample_idx, sample_idx])
                    mean(p_sample[upper.tri(p_sample)] < 0.05, na.rm = TRUE)
                } else {
                    mean(P[upper.tri(P)] < 0.05, na.rm = TRUE)
                }
            }
        )
    }

    step07$done()

    ## --- 8. Write to stats (main FDR replacement; retain old matrices; add FDR_storey) ---
    step08 <- .log_step(parent, "S08", "store outputs", verbose)
    layer_name <- if (is.null(lee_stats_layer_name)) paste0("LeeStats_", norm_layer) else lee_stats_layer_name
    step08$enter(paste0("layer_name=", layer_name))
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()

    scope_obj@stats[[gname]][[layer_name]] <- list(
        L = L,
        Z = Z_mat,
        P = P,
        grad = betas,
        L_min = L_min,
        qc = qc,
        FDR = FDR_main, # Main recommendation (q-value or beta-BH fallback)
        FDR_storey = FDR_storey, # Available in in-memory mode only
        FDR_disc = FDR_out_disc, # Original discrete BH
        FDR_beta = FDR_out_beta,
        FDR_mid = FDR_out_mid,
        FDR_uniform = FDR_out_uniform,
        meta = list(
            perms = perms,
            block_size = block_size,
            ncores = ncores,
            mem_mode = if (inherits(L, "big.matrix")) "bigmemory" else "inmemory",
            p_source = if (!is.null(P)) "permutation" else "analytic",
            p_resolution = if (!is.null(P)) paste0("~", format(1 / (perms + 1), scientific = TRUE)) else "analytic",
            pi0_hat = pi0_hat,
            n_tests = length(FDR_main),
            n_sig_FDR_lt_0.05 = n_sig_005,
            min_p_possible = min_p_possible,
            FDR_main_method = FDR_main_method,
            pval_modes = c("disc=(k+1)/(B+1)", "beta=(k+1)/(B+2)", "mid=(k+0.5)/(B+1)", "uniform=(k+U)/(B+1)"),
            note = "FDR: main=Storey q (beta p) if possible; fallback=BH(beta). All matrices retained for diagnostics."
        )
    )

    if (cache_inputs) {
        if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()
        scope_obj@stats[[gname]][["_lee_input_cache"]] <- input_cache
    }

    step08$done(paste0("layer_name=", layer_name))

    invisible(scope_obj)
}

#' Bootstrap Lee's L vs Pearson r relationship
#' @description
#' Internal helper for `.compute_l_vs_r_curve`.
#' Computes a LOESS-based relationship between Lee's L and Pearson correlation
#' with bootstrap confidence intervals.
#' @param scope_obj A `scope_object` containing Lee's L and correlation matrices.
#' @param grid_name Grid layer name to operate on.
#' @param level Correlation level (`cell` or `grid`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param span LOESS span parameter.
#' @param B Number of bootstrap iterations.
#' @param deg Degree for the LOESS fit.
#' @param ncores Number of threads to use.
#' @param length_out Grid size for the fitted curve.
#' @param downsample Optional downsampling rate for points.
#' @param n_strata Number of strata used when drawing bootstrap samples.
#' @param jitter_eps Optional jitter applied to correlation values.
#' @param ci_method Confidence interval method to use.
#' @param ci_adjust Optional analytic adjustment for CI width.
#' @param verbose Emit progress messages when TRUE.
#' @param k_max Numeric threshold.
#' @param min_rel_width Numeric threshold.
#' @param widen_span Parameter value.
#' @param curve_name Parameter value.
#' @return The modified `scope_object` with stored curve data.
#' @keywords internal
.compute_l_vs_r_curve <- function(scope_obj,
                        grid_name,
                        level = c("grid", "cell"),
                        lee_stats_layer = "LeeStats_Xz",
                        span = 0.45,
                        B = 1000,
                        deg = 1,
                        ncores = max(1, detectCores() - 1),
                        length_out = 1000,
                        downsample = 1,
                        n_strata = 50,
                        k_max = Inf,
                        jitter_eps = 0,
                        ci_method = c("percentile", "basic", "bc"),
                        ci_adjust = c("none", "analytic"),
                        min_rel_width = 0,
                        widen_span = 0.1,
                        curve_name = "LR_curve2",
                        verbose = TRUE) {
  ci_method <- match.arg(ci_method)
  ci_adjust <- match.arg(ci_adjust)
  level <- match.arg(level)
  ncores <- max(1L, min(as.integer(ncores), detectCores(logical = TRUE)))
  if (min_rel_width < 0) stop("min_rel_width cannot be negative")
  parent <- "computeLvsRCurve"
  if (B < 20 && verbose) {
    .log_info(parent, "S01", "Warning: B < 20 may be unstable", verbose)
  }

  step01 <- .log_step(parent, "S01", "extract L and correlation", verbose)
  step01$enter(paste0("grid_name=", grid_name, " level=", level))
  .log_info(parent, "S01", paste0(
    "bootstrap_iterations=", B,
    " ci_method=", ci_method,
    " ncores=", ncores
  ), verbose)

  # 1. Extract matrices
  if (verbose) .log_info(parent, "S01", "extracting Lee's L and Pearson correlation matrices", verbose)
  g_layer <- .select_grid_layer(scope_obj, grid_name)
  grid_name <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
  Lmat <- .get_lee_matrix(scope_obj, grid_name, lee_layer = lee_stats_layer)
  rmat <- .get_pearson_matrix(scope_obj, grid_name, level = ifelse(level == "grid", "grid", "cell"))

  common <- intersect(rownames(Lmat), rownames(rmat))
  if (length(common) < 2) stop("Insufficient common genes")

  if (verbose) {
    .log_info(parent, "S01", paste0("common_genes=", length(common)), verbose)
    .log_info(parent, "S01", paste0("matrix_dim=", nrow(Lmat), "x", ncol(Lmat)), verbose)
  }

  step01$done(paste0("common_genes=", length(common)))

  step02 <- .log_step(parent, "S02", "prepare data vectors", verbose)
  step02$enter(paste0("downsample=", downsample, " jitter_eps=", jitter_eps))

  # Memory guard: estimate per-thread footprint (two dense matrices) and stop if over system RAM
  sys_mem_gb <- .get_system_memory_gb()
  per_thread_gb <- (length(common)^2 * 8 * 2) / (1024^3)
  est_total_gb <- per_thread_gb * ncores
  if (verbose) {
    .log_info(parent, "S02", paste0(
      "est_total_gb=", round(est_total_gb, 1),
      " sys_mem_gb=", round(sys_mem_gb, 1)
    ), verbose)
  }
  if (est_total_gb > sys_mem_gb) {
    stop(
      "[geneSCOPE::.compute_l_vs_r_curve] Estimated memory requirement (",
      round(est_total_gb, 1), " GB) exceeds system capacity (",
      round(sys_mem_gb, 1), " GB). Reduce ncores, downsample, or gene set size."
    )
  }

  Lmat <- Lmat[common, common, drop = FALSE]
  rmat <- rmat[common, common, drop = FALSE]
  ut <- upper.tri(Lmat, diag = FALSE)
  Lv <- Lmat[ut]
  rv <- rmat[ut]

  # 2. Clean / Downsample
  if (verbose) .log_info(parent, "S02", "cleaning and preprocessing data points", verbose)
  ok <- is.finite(Lv) & is.finite(rv)
  Lv <- Lv[ok]
  rv <- rv[ok]
  if (!length(Lv)) stop("No valid points")

  if (verbose) .log_info(parent, "S02", paste0("initial_points=", length(Lv)), verbose)

  if (is.numeric(downsample) && downsample < 1) {
    keep <- sample.int(length(Lv), max(1L, floor(downsample * length(Lv))))
    Lv <- Lv[keep]
    rv <- rv[keep]
    if (verbose) .log_info(parent, "S02", paste0(
      "downsampled_points=", length(Lv), " ratio=", downsample
    ), verbose)
  } else if (is.numeric(downsample) && downsample >= 1 && length(Lv) > downsample) {
    keep <- sample.int(length(Lv), downsample)
    Lv <- Lv[keep]
    rv <- rv[keep]
    if (verbose) .log_info(parent, "S02", paste0(
      "downsampled_points=", length(Lv), " target=", downsample
    ), verbose)
  }
  if (jitter_eps > 0) {
    rv <- jitter(rv, factor = jitter_eps)
    if (verbose) .log_info(parent, "S02", paste0("jitter_eps=", jitter_eps), verbose)
  }

  step02$done(paste0("points=", length(Lv)))

  # 3. Stratify (by r quantiles) — fallback if failed
  step03 <- .log_step(parent, "S03", "stratify and build grid", verbose)
  step03$enter(paste0("n_strata=", n_strata, " length_out=", length_out))
  if (verbose) .log_info(parent, "S03", "setting up stratified sampling", verbose)
  uniq_r <- sort(unique(rv))
  if (length(uniq_r) < 3) stop("r values too discrete")
  n_strata_eff <- min(n_strata, length(uniq_r) - 1)

  if (verbose) {
    .log_info(parent, "S03", paste0("unique_r=", length(uniq_r)), verbose)
    .log_info(parent, "S03", paste0("effective_strata=", n_strata_eff), verbose)
  }

  probs <- seq(0, 1, length.out = n_strata_eff + 1)
  brks <- unique(quantile(rv, probs, na.rm = TRUE))
  if (length(brks) < 2) {
    # Fallback: uniform slices
    brks <- seq(min(rv), max(rv), length.out = n_strata_eff + 1)
    brks <- unique(brks)
    if (verbose) .log_info(parent, "S03", "using fallback uniform stratification", verbose)
  }
  if (length(brks) < 2) stop("Cannot establish strata")
  strat <- cut(rv, breaks = brks, include.lowest = TRUE, labels = FALSE)
  ok2 <- !is.na(strat)
  rv <- rv[ok2]
  Lv <- Lv[ok2]
  strat <- strat[ok2]

  # 4. xgrid
  if (verbose) .log_info(parent, "S03", "preparing analysis grid and fitting LOESS model", verbose)
  xr <- range(rv)
  if (diff(xr) <= 0) stop("r has no span")
  xgrid <- seq(xr[1], xr[2], length.out = length_out)

  if (verbose) {
    .log_info(parent, "S03", paste0(
      "analysis_range=[", round(xr[1], 3), ", ", round(xr[2], 3), "]"
    ), verbose)
    .log_info(parent, "S03", paste0("grid_points=", length_out), verbose)
    .log_info(parent, "S03", paste0("bootstrap_cores=", ncores), verbose)
  }

  step03$done(paste0("xgrid_points=", length_out))

  # Directly call unified interface (no longer check old version)
  step04 <- .log_step(parent, "S04", "bootstrap LOESS curve", verbose)
  step04$enter(paste0("B=", B, " span=", span, " deg=", deg))
  .log_backend(parent, "S04", "loess_bootstrap",
    paste0("C++ loess_residual_bootstrap threads=", ncores),
    verbose = verbose
  )
  if (verbose) .log_info(parent, "S04", "running LOESS residual bootstrap analysis", verbose)
  keep_boot <- TRUE
  adjust_mode <- if (ci_method == "percentile" && ci_adjust == "analytic") 1L else 0L
  res <- .loess_residual_bootstrap(
    x = rv, y = Lv, strat = as.integer(strat),
    grid = xgrid,
    B = as.integer(B),
    span = span,
    deg = as.integer(deg),
    n_threads = as.integer(max(1, ncores)),
    k_max = if (is.finite(k_max)) as.integer(k_max) else -1L,
    keep_boot = keep_boot,
    adjust_mode = adjust_mode,
    ci_type = 0L,
    level = 0.95
  )

  fit <- res$fit
  lo <- res$lo
  hi <- res$hi
  if (ci_method == "basic") {
    lo <- res$lo_basic
    hi <- res$hi_basic
  } else if (ci_method == "bc") {
    lo <- res$lo_bc
    hi <- res$hi_bc
  }

  step04$done(paste0("B=", res$B))

  step05 <- .log_step(parent, "S05", "adjust confidence intervals", verbose)
  step05$enter(paste0("ci_method=", ci_method, " ci_adjust=", ci_adjust))
  if (verbose) {
    .log_info(parent, "S05", "bootstrap completed, processing confidence intervals", verbose)
    .log_info(parent, "S05", paste0("ci_method_applied=", ci_method), verbose)
  }

  # 6. floor widen
  if (verbose && min_rel_width > 0) {
    .log_info(parent, "S05", "applying minimum relative width constraints", verbose)
  }
  if (min_rel_width > 0) {
    rng_fit <- diff(range(fit, finite = TRUE))
    if (!is.finite(rng_fit) || rng_fit <= 0) rng_fit <- 1
    width <- hi - lo
    target_w <- min_rel_width * rng_fit
    width <- pmax(width, target_w)
    # Smooth width (LOESS)
    if (is.finite(widen_span) && widen_span > 0 && length(width) > 10) {
      sm <- tryCatch(loess(width ~ xgrid, span = widen_span), error = function(e) NULL)
      if (!is.null(sm)) {
        wp <- tryCatch(predict(sm, xgrid), error = function(e) NULL)
        if (!is.null(wp) && all(is.finite(wp))) {
          width <- pmax(wp, target_w)
        }
      }
    }
    center <- (lo + hi) / 2
    lo <- center - width / 2
    hi <- center + width / 2
  }

  # ---- New: Local residual scale adaptive lower bound (no new parameters) ----------------------------
  # Purpose: Avoid "points more dispersed but CI narrower"; force width >= 2*1.96*local MAD
  if (verbose) .log_info(parent, "S05", "applying local residual scale adaptation", verbose)
  {
    if (length(Lv) > 20) {
      # Fit values interpolated to original r
      fit_at_rv <- tryCatch(approx(xgrid, fit, xout = rv, rule = 2, ties = "ordered")$y,
        error = function(e) rep(mean(fit), length(rv))
      )
      res_raw <- Lv - fit_at_rv
      ord_r <- order(rv)
      rv_sorted <- rv[ord_r]
      res_sorted <- res_raw[ord_r]
      n_pts <- length(rv_sorted)
      # Fixed K (internal constant, no new parameters)
      K <- min(200L, n_pts)
      if (K > 5) {
        # Nearest neighbor index function (bidirectional expansion)
        get_knn_idx <- function(x0) {
          pos <- findInterval(x0, rv_sorted)
          l <- pos
          r <- pos + 1
          out <- integer(0)
          while (length(out) < K && (l >= 1 || r <= n_pts)) {
            dl <- if (l >= 1) abs(rv_sorted[l] - x0) else Inf
            dr <- if (r <= n_pts) abs(rv_sorted[r] - x0) else Inf
            if (dl <= dr) {
              if (l >= 1) {
                out <- c(out, l)
                l <- l - 1
              }
            } else {
              if (r <= n_pts) {
                out <- c(out, r)
                r <- r + 1
              }
            }
          }
          out
        }
        # Calculate local MAD
        loc_mad <- vapply(xgrid, function(xx) {
          idx <- get_knn_idx(xx)
          if (!length(idx)) {
            return(NA_real_)
          }
          rr <- res_sorted[idx]
          med <- median(rr)
          1.4826 * median(abs(rr - med))
        }, numeric(1))
        # Global benchmark (prevent all NA)
        med_mad <- median(loc_mad[is.finite(loc_mad) & loc_mad > 0], na.rm = TRUE)
        if (is.finite(med_mad) && med_mad > 0) {
          width <- hi - lo
          # Required minimum width
          w_need <- 2 * 1.96 * loc_mad
          # If local MAD is NA or 0, no restriction
          bad <- !is.finite(w_need) | w_need <= 0
          if (any(!bad)) {
            width_new <- width
            width_new[!bad] <- pmax(width[!bad], w_need[!bad])
            if (!identical(width_new, width)) {
              center <- (lo + hi) / 2
              lo <- center - width_new / 2
              hi <- center + width_new / 2
              # Record diagnostics (add to meta)
              local_mad_diag <- list(
                K = K,
                mad_global_med = med_mad,
                frac_expanded = mean(width_new > width)
              )
            }
          }
        }
      }
    }
  }

  # ---- New: CI edge smoothing (post-processing only, does not change algorithm) ---------------------------
  if (verbose) .log_info(parent, "S05", "applying confidence interval edge smoothing", verbose)
  {
    if (length(lo) >= 15 && all(is.finite(xgrid))) {
      span_s <- 0.06 + 10 / length(lo) # Adaptive small span, smaller with longer length
      lo_s <- tryCatch(
        predict(loess(lo ~ xgrid,
          span = span_s, degree = 1,
          surface = "direct", family = "gaussian"
        )),
        error = function(e) lo
      )
      hi_s <- tryCatch(
        predict(loess(hi ~ xgrid,
          span = span_s, degree = 1,
          surface = "direct", family = "gaussian"
        )),
        error = function(e) hi
      )
      # Use smoothed width as main, keep original center to avoid overall drift
      center_old <- (lo + hi) / 2
      width_new <- pmax(hi_s - lo_s, 0)
      lo_new <- center_old - width_new / 2
      hi_new <- center_old + width_new / 2
      # Prevent numerical error causing reversal
      swap_idx <- which(hi_new < lo_new)
      if (length(swap_idx)) {
        tmp <- lo_new[swap_idx]
        lo_new[swap_idx] <- hi_new[swap_idx]
        hi_new[swap_idx] <- tmp
      }
      # If NA produced, fallback
      if (all(is.finite(lo_new)) && all(is.finite(hi_new))) {
        edge_smooth_info <- list(
          edge_smooth = TRUE,
          edge_smooth_span = span_s,
          edge_smooth_expanded = mean(width_new > (hi - lo))
        )
        lo <- lo_new
        hi <- hi_new
      } else {
        edge_smooth_info <- list(edge_smooth = FALSE)
      }
    } else {
      edge_smooth_info <- list(edge_smooth = FALSE)
    }
  }

  # ---- CI edge smoothing completed, enter length check ----
  step05$done()

  step06 <- .log_step(parent, "S06", "store curve", verbose)
  step06$enter(paste0("curve_name=", curve_name))
  if (verbose) .log_info(parent, "S06", "finalizing curve data and storing results", verbose)
  ## ---- 6. Write back to @stats -------------------------------------------------------
  n_grid <- length(xgrid)
  if (!all(
    length(fit) == n_grid,
    length(lo) == n_grid,
    length(hi) == n_grid
  )) {
    stop(sprintf(
      "Internal length mismatch: fit=%d lo=%d hi=%d xgrid=%d",
      length(fit), length(lo), length(hi), n_grid
    ))
  }

  if (is.null(scope_obj@stats)) scope_obj@stats <- list()
  if (is.null(scope_obj@stats[[grid_name]])) {
    scope_obj@stats[[grid_name]] <- list()
  }
  if (is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
    scope_obj@stats[[grid_name]][[lee_stats_layer]] <- list()
  }

  ## ---- meta update ----
  meta_obj <- list(
    B = res$B,
    span = span,
    deg = deg,
    ci_method = ci_method,
    ci_adjust = ci_adjust,
    min_rel_width = min_rel_width,
    n_strata = n_strata,
    k_max = k_max,
    edf = res$edf,
    sigma2_raw = res$sigma2_raw,
    sigma2_edf = res$sigma2_edf,
    resid_global_mad = res$resid_mad,
    adjust_mode = adjust_mode,
    note = "Generated by .compute_l_vs_r_curve with residual-MAD floor + edge smoothing",
    local_mad_diag = if (exists("local_mad_diag")) local_mad_diag else NULL,
    edge_smooth = if (exists("edge_smooth_info")) edge_smooth_info else NULL
  )
  meta_col <- rep(list(meta_obj), length(fit))

  scope_obj@stats[[grid_name]][[lee_stats_layer]][[curve_name]] <-
    data.frame(Pear = xgrid, fit = fit, lo95 = lo, hi95 = hi, meta = I(meta_col))

  if (verbose) {
    .log_info(parent, "S06", "analysis completed successfully", verbose)
    .log_info(parent, "S06", paste0("curve_name=", curve_name), verbose)
    .log_info(parent, "S06", paste0("curve_points=", length(xgrid)), verbose)
    .log_info(parent, "S06", paste0("bootstrap_iterations=", res$B), verbose)
  }

  step06$done(paste0("curve_name=", curve_name))

  invisible(scope_obj)
}

#' Rank gene pairs by Lee's L vs Pearson r
#' @description
#' Internal helper for `.get_top_l_vs_r`.
#' Retrieves top gene pairs with large deviations between Lee's L and Pearson
#' correlation, with optional confidence interval filtering and permutation
#' testing.
#' @param scope_obj A `scope_object` containing Lee's L and correlation matrices.
#' @param grid_name Grid layer name.
#' @param pear_level Correlation level (`cell` or `grid`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param curve_layer Optional precomputed curve layer from `.compute_l_vs_r_curve`.
#' @param direction Direction of Delta selection (`both`, `pos`, `neg`).
#' @param top_n Number of top pairs to return.
#' @param ncores Number of threads to use for permutation testing.
#' @param perms Number of permutations for p-value estimation.
#' @param verbose Emit progress messages when TRUE.
#' @param expr_layer Layer name.
#' @param pear_range Parameter value.
#' @param L_range Parameter value.
#' @param do_perm Parameter value.
#' @param block_side Parameter value.
#' @param use_blocks Logical flag.
#' @param clamp_mode Parameter value.
#' @param p_adj_mode Parameter value.
#' @param mem_limit_GB Parameter value.
#' @param pval_mode Parameter value.
#' @param CI_rule Parameter value.
#' @return A data.frame of top pairs with statistics.
#' @keywords internal
.get_top_l_vs_r <- function(scope_obj,
                       grid_name,
                       pear_level = c("cell", "grid"),
                       lee_stats_layer = "LeeStats_Xz",
                       expr_layer = NULL,
                       pear_range = c(-1, 1),
                       L_range = c(-1, 1),
                       top_n = 10,
                       direction = c("largest", "smallest", "both"),
                       do_perm = TRUE,
                       perms = 1000,
                       block_side = 8,
                       use_blocks = TRUE,
                       ncores = 1,
                       clamp_mode = c("none", "ref_only", "both"),
                       p_adj_mode = c("BH", "BY", "BH_universe", "BY_universe", "bonferroni"),
                       mem_limit_GB = 2,
                       pval_mode = c("beta", "mid", "uniform"),
                       curve_layer = NULL,
                       CI_rule = c("remove_within", "remove_outside", "none"),
                       verbose = TRUE) {
  pear_level <- match.arg(pear_level)
  direction <- match.arg(direction)
  p_adj_mode <- match.arg(p_adj_mode)
  clamp_mode <- match.arg(clamp_mode)
  pval_mode <- match.arg(pval_mode)
  CI_rule <- match.arg(CI_rule)

  parent <- "getTopLvsR"
  step01 <- .log_step(parent, "S01", "load matrices and compute Delta", verbose)
  step01$enter(paste0("grid_name=", grid_name, " pear_level=", pear_level))
  if (verbose) {
    .log_info(parent, "S01", "starting top Delta analysis", verbose)
    .log_info(parent, "S01", paste0("direction=", direction, " top_n=", top_n), verbose)
  }

  logi <- detectCores(TRUE)
  safe_cores <- max(1L, min(ncores, logi))
  if (ncores > safe_cores) {
    .log_info(parent, "S01", paste0(
      "adjusting ncores: requested=", ncores,
      " capped at available=", safe_cores
    ), verbose)
    ncores <- safe_cores
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
  Sys.setenv(OMP_NUM_THREADS = ncores)

  L_mat <- .get_lee_matrix(scope_obj, grid_name, lee_layer = lee_stats_layer)
  r_mat <- .get_pearson_matrix(scope_obj, grid_name, level = pear_level)
  common <- intersect(rownames(L_mat), rownames(r_mat))
  if (length(common) < 2) stop("Insufficient common genes")

  if (verbose) {
    .log_info(parent, "S01", paste0("common_genes=", length(common)), verbose)
    .log_info(parent, "S01", paste0("matrix_dim=", nrow(L_mat), "x", ncol(L_mat)), verbose)
  }

  L_mat <- L_mat[common, common]
  r_mat <- r_mat[common, common]
  diag(L_mat) <- NA
  diag(r_mat) <- NA
  ut <- upper.tri(L_mat)
  LeesL_vec <- L_mat[ut]
  Pear_vec <- r_mat[ut]

  if (verbose) .log_info(parent, "S01", "computing Delta values and applying filters", verbose)
  internal_clamp_mode <- clamp_mode
  if (clamp_mode == "both") {
    .log_info(parent, "S01",
      "clamp_mode='both' currently equivalent to ref_only (reference truncation only)",
      verbose
    )
    internal_clamp_mode <- "ref_only"
  }
  Pear_for_delta <- if (internal_clamp_mode == "ref_only") pmax(Pear_vec, 0) else Pear_vec

  df <- data.frame(
    gene1 = rep(common, each = length(common))[ut],
    gene2 = rep(common, length(common))[ut],
    LeesL = LeesL_vec,
    Pear  = Pear_vec,
    Delta = LeesL_vec - Pear_for_delta
  )

  # Expression coverage percentages
  {
    expr_pct_map <- setNames(rep(0, length(common)), common)
    g_layer_try <- tryCatch(.select_grid_layer(scope_obj, grid_name), error = function(e) NULL)
    coverage_done <- FALSE

    if (!is.null(g_layer_try) && !is.null(g_layer_try$counts)) {
      ct <- g_layer_try$counts
      if (is.data.frame(ct) && all(c("gene", "grid_id") %in% colnames(ct))) {
        total_cells <- if (!is.null(g_layer_try$grid_info)) nrow(g_layer_try$grid_info) else length(unique(ct$grid_id))
        if (total_cells <= 0) total_cells <- NA_real_
        if (inherits(ct, "data.table")) {
          gene_cells <- ct[gene %in% common, .(cells = uniqueN(grid_id)), by = gene]
        } else if (requireNamespace("data.table", quietly = TRUE)) {
          dct <- as.data.table(ct)
          gene_cells <- dct[gene %in% common, .(cells = uniqueN(grid_id)), by = gene]
        } else {
          keep <- ct$gene %in% common
          if (any(keep)) {
            spl <- split(ct$grid_id[keep], ct$gene[keep])
            gene_cells <- data.frame(
              gene = names(spl),
              cells = vapply(spl, function(v) length(unique(v)), integer(1)),
              row.names = NULL
            )
          } else {
            gene_cells <- data.frame(gene = character(0), cells = integer(0))
          }
        }
        if (nrow(gene_cells)) {
          expr_pct_map[gene_cells$gene] <- (gene_cells$cells / total_cells) * 100
        }
        coverage_done <- TRUE
      }
    }

    if (!coverage_done) {
      pick_matrix <- function(layer) {
        if (is.null(layer)) {
          return(NULL)
        }
        for (nm in c("counts", "raw_counts", "expr", "X", "data", "logCPM", "Xz")) {
          m <- layer[[nm]]
          if (!is.null(m) && (is.matrix(m) || inherits(m, "dgCMatrix"))) {
            return(m)
          }
        }
        NULL
      }
      expr_mat <- pick_matrix(g_layer_try)
      if (is.null(expr_mat) && pear_level == "cell") {
        cell_env <- tryCatch(scope_obj@cell, error = function(e) NULL)
        if (is.null(cell_env)) cell_env <- tryCatch(scope_obj@cells, error = function(e) NULL)
        expr_mat <- pick_matrix(cell_env)
      }
      if (!is.null(expr_mat)) {
        rn <- rownames(expr_mat)
        cn <- colnames(expr_mat)
        inter_col <- intersect(cn, common)
        inter_row <- intersect(rn, common)
        if (length(inter_col) == 0 && length(inter_row) == 0) {
          # keep zeros
        } else {
          gene_in_col <- (length(inter_col) >= length(inter_row))
          if (!gene_in_col) {
            expr_mat <- if (inherits(expr_mat, "dgCMatrix")) t(expr_mat) else t(expr_mat)
            cn <- colnames(expr_mat)
            inter_col <- intersect(cn, common)
          }
          if (length(inter_col)) {
            nz_counts <- if (inherits(expr_mat, "dgCMatrix")) {
              colSums(expr_mat[, inter_col, drop = FALSE] > 0)
            } else {
              colSums(expr_mat[, inter_col, drop = FALSE] > 0)
            }
            expr_pct_map[inter_col] <- (nz_counts / nrow(expr_mat)) * 100
          }
        }
      }
    }

    df$gene1_expr_pct <- round(expr_pct_map[df$gene1], 3)
    df$gene2_expr_pct <- round(expr_pct_map[df$gene2], 3)
  }

  step01$done(paste0("pairs_total=", nrow(df)))

  # Curve filtering
  step02 <- .log_step(parent, "S02", "apply curve/CI filters", verbose)
  step02$enter(paste0("curve_layer=", if (is.null(curve_layer)) "none" else curve_layer, " CI_rule=", CI_rule))
  if (!is.null(curve_layer) || CI_rule != "none") {
    if (is.null(curve_layer)) stop("curve_layer must be provided when CI_rule != 'none'")
    curve_obj <- scope_obj@stats[[grid_name]][[lee_stats_layer]][[curve_layer]]
    if (is.null(curve_obj)) stop("Curve layer '", curve_layer, "' not found in stats")
    required_cols <- c("Pear", "lo95", "hi95")
    if (!all(required_cols %in% colnames(curve_obj))) {
      stop("Curve layer must contain columns: ", paste(required_cols, collapse = ", "))
    }
    lo_fun <- approxfun(curve_obj$Pear, curve_obj$lo95, rule = 2)
    hi_fun <- approxfun(curve_obj$Pear, curve_obj$hi95, rule = 2)
    df$curve_lo <- lo_fun(df$Pear)
    df$curve_hi <- hi_fun(df$Pear)
    df$outside_ci <- is.na(df$curve_lo) | is.na(df$curve_hi) | df$LeesL < df$curve_lo | df$LeesL > df$curve_hi
    if (CI_rule == "remove_within") {
      df <- df[df$outside_ci, , drop = FALSE]
    } else if (CI_rule == "remove_outside") {
      df <- df[!df$outside_ci, , drop = FALSE]
    }
    if (!nrow(df)) stop("No gene pairs satisfy the requested CI_rule")
  }

  if (verbose && !is.null(df$outside_ci) && CI_rule != "none") {
    if (CI_rule == "remove_within") {
      .log_info(parent, "S02", paste0("pairs_outside_ci95=", sum(df$outside_ci)), verbose)
    } else if (CI_rule == "remove_outside") {
      .log_info(parent, "S02", paste0("pairs_inside_ci95=", sum(!df$outside_ci)), verbose)
    }
  }

  step02$done(paste0("pairs_after_ci=", nrow(df)))

  step03 <- .log_step(parent, "S03", "filter ranges and select pairs", verbose)
  step03$enter(paste0(
    "pear_range=[", pear_range[1], ",", pear_range[2], "]",
    " L_range=[", L_range[1], ",", L_range[2], "]"
  ))
  if (verbose) .log_info(parent, "S03", "applying threshold filters", verbose)
  df <- df[df$Pear >= pear_range[1] & df$Pear <= pear_range[2] &
    df$LeesL >= L_range[1] & df$LeesL <= L_range[2], ]
  if (!nrow(df)) stop("Thresholds remove all pairs")
  total_universe <- nrow(df)

  if (verbose) {
    .log_info(parent, "S03", paste0("pairs_after_filtering=", total_universe), verbose)
    .log_info(parent, "S03", paste0(
      "pear_range=[", pear_range[1], ", ", pear_range[2], "]",
      " L_range=[", L_range[1], ", ", L_range[2], "]"
    ), verbose)
  }

  if (verbose) .log_info(parent, "S03", "selecting top gene pairs by Delta values", verbose)
  sel <- switch(direction,
    largest = slice_max(df, Delta, n = top_n),
    smallest = slice_min(df, Delta, n = top_n),
    both = bind_rows(
      slice_max(df, Delta, n = top_n),
      slice_min(df, Delta, n = top_n)
    )
  )
  sel <- distinct(sel, gene1, gene2, .keep_all = TRUE)
  rownames(sel) <- NULL

  if (verbose) {
    .log_info(parent, "S03", paste0("selected_pairs=", nrow(sel)), verbose)
    if (nrow(sel) > 0) {
      delta_range <- range(sel$Delta)
      .log_info(parent, "S03", paste0(
        "delta_range=[", round(delta_range[1], 4), ", ", round(delta_range[2], 4), "]"
      ), verbose)
    }
  }

  step03$done(paste0("selected_pairs=", nrow(sel)))

  if (!nrow(sel) || !do_perm) {
    step04 <- .log_step(parent, "S04", "permutation testing", verbose)
    step04$enter(paste0("perms=", perms, " do_perm=", do_perm))
    .log_backend(parent, "S04", "permutation", "skipped",
      reason = if (!nrow(sel)) "no_pairs" else "do_perm=FALSE",
      verbose = verbose
    )
    step04$done("skipped")

    step05 <- .log_step(parent, "S05", "return results", verbose)
    step05$enter("no permutation results")
    sel$FDR <- NA_real_
    out <- transmute(
      sel,
      gene1,
      gene2,
      L = LeesL,
      r = Pear,
      pct1 = gene1_expr_pct,
      pct2 = gene2_expr_pct,
      fdr = FDR
    )
    step05$done(paste0("pairs_returned=", nrow(out)))
    return(out)
  }

  step04 <- .log_step(parent, "S04", "permutation testing", verbose)
  step04$enter(paste0("perms=", perms, " pval_mode=", pval_mode))
  if (verbose) {
    .log_info(parent, "S04", "preparing permutation analysis", verbose)
    .log_info(parent, "S04", paste0("permutations=", perms), verbose)
    .log_info(parent, "S04", paste0("pval_mode=", pval_mode), verbose)
    .log_info(parent, "S04", paste0("p_adj_mode=", p_adj_mode), verbose)
  }
  g_layer <- .select_grid_layer(scope_obj, grid_name)
  # Infer expression layer from lee_stats_layer if not specified
  if (is.null(expr_layer) || !nzchar(expr_layer)) {
    inferred <- sub("^LeeStats_", "", lee_stats_layer)
    expr_layer <- if (identical(inferred, lee_stats_layer)) "Xz" else inferred
    if (verbose) {
      .log_info(parent, "S04", paste0(
        "expr_layer=", expr_layer, " (inferred from lee_stats_layer)"
      ), verbose)
    }
  }
  Xcand <- g_layer[[expr_layer]]
  if (is.null(Xcand) && expr_layer != "Xz") Xcand <- g_layer$Xz
  Xz <- Xcand
  W <- g_layer$W
  grid_info <- g_layer$grid_info
  if (is.null(Xz)) stop("Grid layer missing expression layer '", expr_layer, "' and 'Xz'")
  if (is.null(W)) stop("Grid layer missing spatial weights W; run .compute_weights() first")
  # Align rows of Xz to grid_info if possible
  if (!is.null(rownames(Xz)) && !is.null(grid_info$grid_id)) {
    ord <- match(grid_info$grid_id, rownames(Xz))
    if (!any(is.na(ord))) Xz <- Xz[ord, , drop = FALSE]
  }
  if (!is.matrix(Xz)) Xz <- as.matrix(Xz)
  genes_top <- unique(c(sel$gene1, sel$gene2))
  gene_map <- match(genes_top, colnames(Xz))
  if (any(is.na(gene_map))) stop("Selected genes not found in Xz")
  gene_pairs <- cbind(
    match(sel$gene1, genes_top) - 1L,
    match(sel$gene2, genes_top) - 1L
  )
  delta_ref <- sel$Delta

  backend_label <- if (use_blocks) "C++ delta_lr_perm_csr_block" else "C++ delta_lr_perm_csr"
  .log_backend(parent, "S04", "permutation_backend", paste0(
    backend_label,
    " threads=", ncores,
    " mem_limit_GB=", mem_limit_GB
  ), verbose = verbose)

  block_id <- if (use_blocks) {
    bx <- (grid_info$gx - 1L) %/% block_side
    by <- (grid_info$gy - 1L) %/% block_side
    max_by <- max(by)
    bx * (max_by + 1L) + by + 1L
  } else {
    seq_len(nrow(Xz))
  }
  split_rows <- split(seq_along(block_id), block_id)

  W_rows <- lapply(seq_len(nrow(W)), function(i) {
    nz <- which(W[i, ] != 0)
    if (length(nz)) list(indices = nz - 1L, values = as.numeric(W[i, nz])) else list(indices = integer(0), values = numeric(0))
  })
  W_row_lengths <- vapply(W_rows, function(x) length(x$indices), integer(1))
  W_row_ptr <- c(0L, cumsum(W_row_lengths))
  W_indices <- unlist(lapply(W_rows, `[[`, "indices"), use.names = FALSE)
  W_values <- unlist(lapply(W_rows, `[[`, "values"), use.names = FALSE)
  if (!length(W_indices)) stop("Weight matrix has no non-zero entries")
  Xz_sub <- Xz[, gene_map, drop = FALSE]
  n_cells <- nrow(Xz_sub)

  target_batch <- min(100L, perms)
  max_idx_bytes <- mem_limit_GB * 1024^3 * 0.30
  est_bytes <- function(bs) n_cells * bs * 4
  while (target_batch > 1L && est_bytes(target_batch) > max_idx_bytes) {
    target_batch <- max(1L, floor(target_batch / 2))
  }
  if (verbose) {
    .log_info(parent, "S04", paste0(
      "planned_batch_size=", target_batch,
      " (est_idx_mat_mb=", sprintf("%.2f", est_bytes(target_batch) / 1024^2),
      " limit_mb=", sprintf("%.2f", max_idx_bytes / 1024^2), ")"
    ), verbose)
    .log_info(parent, "S04", "starting permutation testing loop", verbose)
  }
  remaining <- perms
  exceed_count <- rep(0L, nrow(sel))
  perm_threads <- ncores
  attempt <- 1
  while (remaining > 0) {
    bsz <- min(target_batch, remaining)
    success <- FALSE
    while (!success && perm_threads >= 1) {
      if (verbose) .log_info(parent, "S04", sprintf(
        "permutation attempt #%d threads=%d batch=%d remain=%d clamp_mode=%s",
        attempt, perm_threads, bsz, remaining, clamp_mode
      ), verbose)
      attempt <- attempt + 1
      idx_mat <- matrix(integer(0), nrow = n_cells, ncol = bsz)
      for (p in seq_len(bsz)) {
        new_order <- sample(length(split_rows))
        idx_mat[, p] <- unlist(split_rows[new_order], use.names = FALSE) - 1L
      }
      res <- tryCatch(
        {
          if (use_blocks) {
            .delta_lr_perm_csr_block(
              Xz_sub, W_indices, W_values, W_row_ptr, idx_mat,
              as.integer(block_id) - 1L, gene_pairs, delta_ref,
              perm_threads
            )
          } else {
            .delta_lr_perm_csr(
              Xz_sub, W_indices, W_values, W_row_ptr, idx_mat,
              gene_pairs, delta_ref,
              perm_threads
            )
          }
        },
        error = function(e) e
      )
      if (inherits(res, "error")) {
        if (verbose) .log_info(parent, "S04", paste0("batch failed: ", conditionMessage(res)), verbose)
        if (perm_threads > 1) {
          perm_threads <- max(1, floor(perm_threads / 2))
          .log_backend(parent, "S04", "permutation_backend", paste0(
            backend_label,
            " threads=", perm_threads,
            " mem_limit_GB=", mem_limit_GB
          ), reason = "retry", verbose = verbose)
          next
        } else {
          stop("Permutation failed at single-thread: ", conditionMessage(res))
        }
      } else {
        exceed_count <- exceed_count + res
        success <- TRUE
      }
    }
    remaining <- remaining - bsz
  }

  if (verbose) .log_info(parent, "S04", "computing p-values and applying multiple testing correction", verbose)
  N <- perms
  k <- exceed_count
  p_values <- switch(pval_mode,
    beta    = (k + 1) / (N + 2),
    mid     = (k + 0.5) / (N + 1),
    uniform = (k + runif(length(k))) / (N + 1)
  )
  p_values[p_values > 1] <- 1
  mc_se <- sqrt(p_values * (1 - p_values) / N)
  p_ci_lo <- qbeta(0.025, k + 1, N - k + 1)
  p_ci_hi <- qbeta(0.975, k + 1, N - k + 1)

  FDR <- switch(p_adj_mode,
    BH          = p.adjust(p_values, "BH"),
    BY          = p.adjust(p_values, "BY"),
    BH_universe = p.adjust(p_values, "BH", n = total_universe),
    BY_universe = p.adjust(p_values, "BY", n = total_universe),
    bonferroni  = p.adjust(p_values, "bonferroni", n = total_universe)
  )

  if (verbose) {
    sig_count <- sum(p_values < 0.05, na.rm = TRUE)
    fdr_sig_count <- sum(FDR < 0.05, na.rm = TRUE)
    .log_info(parent, "S04", paste0(
      "significant_pairs_p_lt_0.05=", sig_count, "/", length(p_values)
    ), verbose)
    .log_info(parent, "S04", paste0(
      "significant_pairs_fdr_lt_0.05=", fdr_sig_count, "/", length(FDR)
    ), verbose)
  }

  sel$raw_p <- p_values
  sel$mc_se <- mc_se
  sel$p_ci_lo <- p_ci_lo
  sel$p_ci_hi <- p_ci_hi
  sel$FDR <- FDR
  sel$stat_type <- switch(clamp_mode,
    none     = "Delta",
    ref_only = "Delta_refClamp",
    both     = "Delta_clamp_Ronly"
  )
  sel$pval_mode <- pval_mode
  if (p_adj_mode == "bonferroni") {
    sel <- sel[sel$FDR < 0.05, , drop = FALSE]
    if (!nrow(sel)) {
      .log_info(parent, "S04", "no Bonferroni-significant pairs (FDR < 0.05)", verbose)
    }
  }

  step04$done(paste0("pairs_with_fdr=", nrow(sel)))

  step05 <- .log_step(parent, "S05", "return results", verbose)
  step05$enter(paste0("pairs_final=", nrow(sel)))
  out <- transmute(
    sel,
    gene1,
    gene2,
    L = LeesL,
    r = Pear,
    pct1 = gene1_expr_pct,
    pct2 = gene2_expr_pct,
    fdr = FDR
  )
  step05$done(paste0("pairs_returned=", nrow(out)))
  out
}

#' Compute Lee's L for Visium (convenience wrapper)
#' @description
#' Performs Visium-specific prescreening and setup (spatial weights, optional Iδ
#' prescreening, and gene filtering) before dispatching to `computeL()`.
#' @param scope_obj A `scope_object` containing Visium data.
#' @param grid_name Grid layer name to operate on.
#' @param norm_layer Normalised expression layer name.
#' @param use_idelta Use Iδ prescreening when TRUE.
#' @param S_target Target number of genes after prescreening.
#' @param min_detect Minimum detection rate when filtering genes.
#' @param winsor_high Winsorization quantile for expression scaling.
#' @param ncores Number of threads to use (defaults to safe thread count).
#' @param verbose Emit progress messages when TRUE.
#' @param ... Additional arguments (currently unused).
#' @return The modified `scope_object`.
#' @seealso `computeL()`, `computeWeights()`, `computeIDelta()`
#' @keywords internal
computeL_visium <- function(
    scope_obj,
    grid_name = NULL,
    norm_layer = "Xz",
    use_idelta = TRUE,
    S_target = 8000L,
    min_detect = 0.10,
    winsor_high = 0.95,
    ncores = NULL,
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...) {

    parent <- "computeL_visium"

    ## ---- 0. Thread count and grid layer selection ----
    if (is.null(ncores)) {
        ncores <- .get_safe_thread_count(default = 8L)
    }
    grid_label <- if (is.null(grid_name)) "auto" else as.character(grid_name)[1]
    step01 <- .log_step(parent, "S01", "select grid layer and validate", verbose)
    step01$enter(paste0("grid_name=", grid_label, " norm_layer=", norm_layer))

    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else grid_name

    if (is.null(g_layer[[norm_layer]])) {
        stop("computeL_visium: normalized layer '", norm_layer, "' not found.")
    }
    X <- g_layer[[norm_layer]]  # n × g
    if (!is.matrix(X)) X <- as.matrix(X)
    n <- nrow(X); G <- ncol(X)
    if (n < 2L || G < 2L) stop("Insufficient data size to compute Lee's L.")
    .log_info(parent, "S01", paste0("cells=", n, " genes=", G, " ncores=", ncores), verbose)
    step01$done(paste0("grid_name=", grid_name, " genes=", G))

    ## ---- 1. Check/compute spatial weights W ----
    step02 <- .log_step(parent, "S02", "ensure spatial weights", verbose)
    step02$enter(paste0("grid_name=", grid_name))
    if (is.null(g_layer$W) || sum(g_layer$W) == 0) {
        .log_info(parent, "S02", "W missing; running computeWeights(style=B)", verbose)
        .log_backend(parent, "S02", "weights", "computeWeights(style=B)", verbose = verbose)
        scope_obj <- .compute_weights(
            scope_obj,
            grid_name = grid_name,
            style = "B",
            store_mat = TRUE,
            verbose = verbose
        )
        g_layer <- .select_grid_layer(scope_obj, grid_name)
    }
    step02$done(paste0("W=", if (is.null(g_layer$W)) "missing" else "ready"))
    W <- g_layer$W  # dgCMatrix (ensures W exists for .compute_l later)

    ## ---- 2. Iδ and gene prescreening ----
    step03 <- .log_step(parent, "S03", "optional I-delta prescreen", verbose)
    step03$enter(paste0("use_idelta=", use_idelta))
    if (use_idelta) {
        meta_col <- paste0(grid_name, "_iDelta")
        if (is.null(scope_obj@meta.data) || !(meta_col %in% colnames(scope_obj@meta.data))) {
            .log_info(parent, "S03", "computing I-delta", verbose)
            .log_backend(parent, "S03", "idelta", "C++ idelta_sparse_cpp", verbose = verbose)
            scope_obj <- .compute_idelta(
                scope_obj,
                grid_name = grid_name,
                level = "grid",
                ncores = min(8L, ncores),
                verbose = FALSE
            )
        } else {
            .log_info(parent, "S03", "I-delta already present; skipping recompute", verbose)
        }
    } else {
        .log_info(parent, "S03", "I-delta prescreen disabled", verbose)
    }
    step03$done()
    genes_all <- colnames(X)
    if (is.null(genes_all)) genes_all <- rownames(scope_obj@meta.data)
    if (is.null(genes_all)) genes_all <- as.character(seq_len(G))

    ## 2.1 Detection rate (based on counts)
    detect_rate <- NULL
    total_count <- NULL
    if (!is.null(g_layer$counts)) {
        dt <- g_layer$counts
        dt <- dt[dt$count > 0, c("gene", "grid_id", "count")]
        # Fast path when data.table is available
        if (requireNamespace("data.table", quietly = TRUE)) {
            dtt <- as.data.table(dt)
            n_spots <- length(g_layer$grid_info$grid_id)
            dr <- dtt[, .N, by = gene]
            detect_rate <- setNames(as.numeric(dr$N) / n_spots, dr$gene)
            total_count <- dtt[, sum(count), by = gene]
            total_count <- setNames(as.numeric(total_count$V1), total_count$gene)
        } else {
            n_spots <- length(unique(dt$grid_id))
            detect_rate <- tapply(dt$grid_id, dt$gene, function(x) length(unique(x)) / n_spots)
            total_count <- tapply(dt$count, dt$gene, sum)
        }
    }

    ## 2.2 Determine candidate gene set
    keep <- rep(TRUE, length(genes_all))
    names(keep) <- genes_all
    if (!is.null(detect_rate)) {
        low <- detect_rate < min_detect
        keep[names(low)] <- !low
    }
    if (!is.null(total_count)) {
        total_count <- total_count[names(keep)]
        total_count[is.na(total_count)] <- 0
    }
    genes_keep <- names(keep)[keep]

    if (!is.null(total_count)) {
        if (length(genes_keep) > S_target) {
            ord <- order(total_count[genes_keep], decreasing = TRUE)
            genes_keep <- genes_keep[ord][seq_len(S_target)]
        }
    }

    ## ---- 3. Winsorize & compute L ----
    if (!is.null(g_layer[[norm_layer]])) {
        X_use <- g_layer[[norm_layer]]
        if (!is.matrix(X_use)) X_use <- as.matrix(X_use)
        if (is.numeric(winsor_high) && winsor_high > 0 && winsor_high < 1) {
            qv <- apply(X_use, 2, stats::quantile, probs = winsor_high, na.rm = TRUE)
            X_use <- sweep(X_use, 2, qv, pmin)
            g_layer[[norm_layer]] <- X_use
            scope_obj@grid[[grid_name]][[norm_layer]] <- X_use
        }
    }

    scope_obj <- .compute_l(
        scope_obj,
        grid_name = grid_name,
        genes = genes_keep,
        within = TRUE,
        ncores = ncores,
        verbose = verbose
    )
    scope_obj
}

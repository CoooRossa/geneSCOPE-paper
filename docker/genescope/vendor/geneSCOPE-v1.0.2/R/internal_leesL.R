#' Adjust a pairwise p-value matrix over unique unordered gene pairs.
#' @keywords internal
.lee_pairwise_p_adjust <- function(p_mat, method = "BH") {
    p_mat <- as.matrix(p_mat)
    out <- matrix(NA_real_, nrow(p_mat), ncol(p_mat), dimnames = dimnames(p_mat))
    if (nrow(p_mat) == ncol(p_mat)) {
        ut <- upper.tri(p_mat, diag = FALSE)
        vals <- p_mat[ut]
        good <- is.finite(vals)
        adj <- rep(NA_real_, length(vals))
        if (any(good)) adj[good] <- p.adjust(vals[good], method = method)
        out[ut] <- adj
        out[lower.tri(out)] <- t(out)[lower.tri(out)]
        diag(out) <- NA_real_
    } else {
        vals <- as.numeric(p_mat)
        good <- is.finite(vals)
        adj <- rep(NA_real_, length(vals))
        if (any(good)) adj[good] <- p.adjust(vals[good], method = method)
        out[] <- adj
    }
    out
}

#' Convert exact two-sided permutation p-values to signed normal scores.
#' @keywords internal
.lee_empirical_signed_z <- function(L, P) {
    L <- as.matrix(L)
    P <- as.matrix(P)
    if (!identical(dim(L), dim(P))) stop("L and P dimensions must match.", call. = FALSE)
    out <- sign(L) * stats::qnorm(pmax(P / 2, .Machine$double.xmin), lower.tail = FALSE)
    out[!is.finite(P)] <- NA_real_
    out[P >= 1] <- 0
    dimnames(out) <- dimnames(L)
    out
}

#' Add Lee's L statistics to a scope_object
#' @description
#'   High-level wrapper that computes Lee's L, empirical p-values and signed
#'   Z-scores via joint row permutations, FDR, spatial gradients, and quality-control
#'   metrics, and stores everything under a new layer in \code{@stats}.
#' @param scope_obj A \code{scope_object} with at least one populated \code{@grid} slot.
#' @param grid_name Character. Name of the grid sub-layer to process. If
#'   \code{NULL} and only one sub-layer exists, it is selected automatically.
#' @param genes Optional character vector of genes to include; if \code{NULL}
#'   all genes are used.
#' @param within A single logical value. If \code{TRUE} restrict analysis to the selected
#'   gene set on both axes (default). Otherwise compute gene × all.
#' @param ncores Positive integer. Number of cores for parallel processing (default 1).
#' @param block_side Positive integer. Number of grid cells per side for block
#'   partitioning (default 8); used only when `use_blocks = TRUE`.
#' @param perms Non-negative integer. Number of permutations for Monte-Carlo p-values (default 1000).
#' @param block_size Positive integer. Number of permutations processed per batch (default 64).
#' @param L_min Numeric threshold used when building the QC similarity graph.
#' @param norm_layer Character. Name of the normalized expression layer (default
#'   "Xz"). Custom layers must be numeric, finite, and column-centred relative
#'   to their RMS magnitude; unit variance is not required.
#' @param lee_stats_layer_name Character. Name for the output statistics layer.
#' @param legacy_formula Must remain `FALSE`; the non-canonical legacy
#'   denominator is disabled.
#' @param mem_limit_GB Numeric. Maximum allowed dense Lee result size in GiB
#'   (default 2); oversized results fail closed.
#' @param chunk_size Retained for API compatibility; the unsafe chunk route is disabled.
#' @param use_bigmemory Must remain `FALSE`; the former chunk route was
#'   RAM-backed rather than file-backed and is disabled.
#' @param backing_path Character. Directory for temporary files (default tempdir()).
#' @param cache_inputs Logical. Cache preprocessed X/Z/W and block IDs for reuse across calls (default TRUE).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @param ncore Deprecated. Use \code{ncores} instead.
#' @param use_blocks Logical. `FALSE` (the default) applies one joint
#'   unrestricted row shuffle over all grid locations. `TRUE` constrains that
#'   joint shuffle within spatial blocks defined by `block_side`.
#' @return The modified \code{scope_object}.
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
                     use_bigmemory = FALSE,
                     backing_path = tempdir(),
                     cache_inputs = TRUE,
                     verbose = TRUE,
                     ncore = NULL,
                     use_blocks = FALSE) {
    parent <- "computeL"

    ## --- 0. Thread-safe preprocessing: automatic thread management and error recovery ---

    # Handle deprecated parameter
    if (!is.null(ncore)) {
        warning("'ncore' is deprecated, please use 'ncores' instead. Will use its value this time.",
            call. = FALSE, immediate. = TRUE
        )
        ncores <- ncore
    }
    within <- .lee_assert_flag(within, "within")
    ncores <- .lee_assert_positive_integer(ncores, "ncores")
    block_side <- .lee_assert_positive_integer(block_side, "block_side")
    perms <- .lee_assert_nonnegative_integer(perms, "perms")
    block_size <- .lee_assert_positive_integer(block_size, "block_size")
    mem_limit_GB <- .lee_assert_positive_finite(mem_limit_GB, "mem_limit_GB")
    chunk_size <- .lee_assert_positive_integer(chunk_size, "chunk_size")
    legacy_formula <- .lee_assert_flag(legacy_formula, "legacy_formula")
    use_bigmemory <- .lee_assert_flag(use_bigmemory, "use_bigmemory")
    cache_inputs <- .lee_assert_flag(cache_inputs, "cache_inputs")
    verbose <- .lee_assert_flag(verbose, "verbose")
    use_blocks <- .lee_assert_flag(use_blocks, "use_blocks")
    norm_layer <- .lee_assert_layer_name(norm_layer, "norm_layer")
    if (!is.null(grid_name)) grid_name <- .lee_assert_layer_name(grid_name, "grid_name")
    if (!is.null(lee_stats_layer_name)) {
        lee_stats_layer_name <- .lee_assert_layer_name(lee_stats_layer_name, "lee_stats_layer_name")
    }
    if (!is.numeric(L_min) || length(L_min) != 1L || is.na(L_min) || !is.finite(L_min)) {
        stop("L_min must be a single finite number.", call. = FALSE)
    }
    if (!is.null(genes) && (!is.character(genes) || !length(genes) || anyNA(genes) ||
        any(!nzchar(genes)) || anyDuplicated(genes))) {
        stop("genes must be NULL or a non-empty vector of unique gene names.", call. = FALSE)
    }
    if (!is.character(backing_path) || length(backing_path) != 1L || is.na(backing_path) || !nzchar(backing_path)) {
        stop("backing_path must be a single non-empty character string.", call. = FALSE)
    }
    if (isTRUE(legacy_formula)) {
        stop(
            "legacy_formula=TRUE is no longer supported because it selects a non-canonical Lee's L denominator. ",
            "Use the canonical default (legacy_formula=FALSE).",
            call. = FALSE
        )
    }
    if (isTRUE(use_bigmemory)) {
        stop(
            "use_bigmemory=TRUE is disabled because the former Lee chunk route was RAM-backed, not file-backed.",
            call. = FALSE
        )
    }
    if (!within && !is.null(genes) && perms > 0L) {
        stop(
            "Permutation inference for within=FALSE with a gene subset is rectangular and is not supported. ",
            "Use perms=0 for observed gene-by-all Lee's L, or use within=TRUE for permutation inference.",
            call. = FALSE
        )
    }

    # Get system information and resolve the explicit request consistently.
    # When core detection is unavailable, a valid explicit request remains
    # authoritative; when detection succeeds, it remains an upper bound.
    os_type <- .detect_os()
    avail_cores <- .detect_cores_safe(logical = TRUE, fallback = ncores)
    ncores <- .clamp_ncores_safe(ncores, logical = TRUE)
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

    # Cross-platform memory calculation
    all_genes <- if (!is.null(g_layer[[norm_layer]])) {
        ncol(g_layer[[norm_layer]])
    } else {
        stop("Normalized layer '", norm_layer, "' not found")
    }

    n_genes_use <- if (!is.null(genes) && isTRUE(within)) length(genes) else all_genes
    matrix_size_gb <- (as.double(n_genes_use)^2 * 8) / (1024^3)
    input_size_gb <- (as.double(nrow(g_layer[[norm_layer]])) * n_genes_use * 8) / (1024^3)

    # The observed Lee kernel shares X/WZ across workers.  The permutation
    # kernel deliberately keeps one g x g counter per OpenMP worker so that
    # updates are race-free and worker allocation failures can be reported
    # without returning partial counts.  Include those buffers in the guard.
    sys_mem_gb <- .get_system_memory_gb()
    observed_est_gb <- matrix_size_gb + 2 * input_size_gb
    permutation_est_gb <- if (perms > 0L) {
        (2 + ncores_cpp) * matrix_size_gb +
            (1 + 2 * ncores_cpp) * input_size_gb
    } else {
        0
    }
    est_total_gb <- max(observed_est_gb, permutation_est_gb)
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

    mem_reason <- "within_limit"
    if (matrix_size_gb > mem_limit_GB) {
        stop(
            "Lee L output exceeds mem_limit_GB, and the former bigmemory route is disabled because it was RAM-backed. ",
            "Reduce the gene set or raise mem_limit_GB only after confirming available RAM.",
            call. = FALSE
        )
    }

    mem_mode <- "inmemory"
    .log_backend(parent, "S02", "mem_mode", mem_mode, reason = mem_reason, verbose = verbose)

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
                    use_blocks = use_blocks,
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
            if (current_cores <= min_cores) {
                stop(
                    "Lee's L computation failed at one thread: ",
                    conditionMessage(result$error),
                    call. = FALSE
                )
            }
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
    if (perms > 0L && inherits(L, "big.matrix")) {
        stop(
            "Permutation/FDR inference is not supported for a streaming big.matrix Lee result because ",
            "global BH must be computed over all unique unordered pairs. Rerun with use_bigmemory=FALSE ",
            "and sufficient memory, reduce the gene set, or use perms=0 for observed L only.",
            call. = FALSE
        )
    }

    ## --- 4. Inference placeholder ---
    ## Moran's-I moments are not the moments of Lee's L.  Z is therefore
    ## derived from the exact two-sided permutation P matrix after Step 5.
    step04 <- .log_step(parent, "S04", "prepare empirical inference", verbose)
    step04$enter(paste0("within=", within, " mem_mode=", mem_mode))
    Z_mat <- if (perms > 0L) {
        matrix(NA_real_, nrow(L), ncol(L), dimnames = dimnames(L))
    } else {
        NULL
    }
    step04$done("analytic Moran moments disabled")

    block_id <- res$block_id
    input_cache <- res$input_cache

    ## --- 5. Monte Carlo p-values with BLAS control and error recovery ---
    step05 <- .log_step(parent, "S05", "permutation p-values", verbose)
    step05$enter(paste0("perms=", perms, " block_size=", block_size))
    permutation_backend <- if (perms < 1L) {
        "skipped"
    } else if (use_blocks) {
        "C++ lee_perm_block (within-block joint shuffle)"
    } else {
        "C++ lee_perm (all-grid joint shuffle)"
    }
    P <- if (perms > 0) {
        # Temporarily disable BLAS threads for permutation tests
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
            try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
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
            paste0(permutation_backend, " threads=", perm_cores),
            verbose = verbose
        )

        while (!perm_success && perm_cores >= 1) {
            if (verbose && perm_attempt > 1) {
                .log_info(parent, "S05", paste0("retry #", perm_attempt, " with ", perm_cores, " cores"), verbose)
                .log_backend(parent, "S05", "permutation",
                    paste0(permutation_backend, " threads=", perm_cores),
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
                    p_result <- if (use_blocks) {
                        .lee_l_perm_block(X_used, W, L_reg,
                            block_id   = block_id,
                            perms      = perms,
                            block_size = block_size,
                            n_threads  = perm_cores
                        )
                    } else {
                        .lee_l_perm_global(X_used, W, L_reg,
                            perms      = perms,
                            block_size = block_size,
                            n_threads  = perm_cores
                        )
                    }
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
                if (perm_cores <= 1L) {
                    stop(
                        "Lee's L permutation test failed at one thread: ",
                        conditionMessage(perm_result$error),
                        call. = FALSE
                    )
                }
                perm_attempt <- perm_attempt + 1
                perm_cores <- max(floor(perm_cores / 2), 1)
                if (verbose) {
                    .log_info(parent, "S05", paste0("reducing cores to ", perm_cores, " and retrying"), verbose)
                    .log_backend(parent, "S05", "permutation",
                        paste0(permutation_backend, " threads=", perm_cores),
                        reason = "retry",
                        verbose = verbose
                    )
                }
                Sys.sleep(2)
                gc(verbose = FALSE)
            }
        }

        if (!perm_success) stop("Lee's L permutation test failed.", call. = FALSE)
        P
    } else {
        .log_backend(parent, "S05", "permutation", "skipped", reason = "perms<=0", verbose = verbose)
        NULL
    }
    if (!is.null(P)) {
        dimnames(P) <- dimnames(L_reg)
        zero_var <- colSums(X_used^2) <= .Machine$double.eps
        if (any(zero_var)) {
            P[zero_var, ] <- NA_real_
            P[, zero_var] <- NA_real_
        }
        Z_mat <- .lee_empirical_signed_z(L_reg, P)
    }
    step05$done(if (!is.null(P)) "empirical P and signed Z ready" else "no inference")

    ## --- 6. FDR: Three Monte Carlo smoothing strategies as in .get_top_l_vs_r ---
    step06 <- .log_step(parent, "S06", "FDR and q-values", verbose)
    step06$enter(paste0("mem_mode=", mem_mode, " perms=", perms))
    if (is.null(P)) {
        FDR_out_disc <- FDR_out_beta <- FDR_out_mid <- FDR_out_uniform <- NULL
        FDR_storey <- FDR_main <- NULL
        pi0_hat <- NA_real_
        n_sig_005 <- NA_integer_
        min_p_possible <- NA_real_
        FDR_main_method <- "unavailable (perms=0)"
    } else {
        # All adjustments use the intended family: unique unordered,
        # off-diagonal pairs.  No per-chunk or duplicated-triangle BH is allowed.
        p_disc <- P
        k_est <- round(p_disc * (perms + 1) - 1)
        k_est[k_est < 0] <- 0
        k_est[k_est > perms] <- perms
        p_beta <- (k_est + 1) / (perms + 2)
        p_mid <- (k_est + 0.5) / (perms + 1)
        p_uniform <- (k_est + matrix(stats::runif(length(k_est)), nrow = nrow(k_est))) / (perms + 1)
        FDR_out_disc <- .lee_pairwise_p_adjust(p_disc, "BH")
        FDR_out_beta <- .lee_pairwise_p_adjust(p_beta, "BH")
        FDR_out_mid <- .lee_pairwise_p_adjust(p_mid, "BH")
        FDR_out_uniform <- .lee_pairwise_p_adjust(p_uniform, "BH")

        # Exploratory Storey q-values use that same unique-pair family.
        FDR_storey <- matrix(NA_real_, nrow(P), ncol(P), dimnames = dimnames(P))
        ut <- upper.tri(P, diag = FALSE)
        p_vec_all <- p_beta[ut]
        good <- is.finite(p_vec_all)
        p_vec <- p_vec_all[good]
        pi0_hat <- NA_real_
        if (length(p_vec)) {
            lambdas <- seq(0.5, 0.95, by = 0.05)
            pi0_vals <- vapply(lambdas, function(lam) mean(p_vec > lam) / (1 - lam), numeric(1))
            pi0_hat <- min(1, min(pi0_vals, na.rm = TRUE))
            if (!is.finite(pi0_hat) || pi0_hat <= 0) pi0_hat <- 1
            o <- order(p_vec)
            q_ord <- pi0_hat * length(p_vec) * p_vec[o] / seq_along(o)
            if (length(q_ord) > 1L) {
                for (i in (length(q_ord) - 1L):1L) q_ord[i] <- min(q_ord[i], q_ord[i + 1L])
            }
            q_vec <- rep(NA_real_, length(p_vec_all))
            q_good <- numeric(length(p_vec))
            q_good[o] <- pmin(q_ord, 1)
            q_vec[good] <- q_good
            FDR_storey[ut] <- q_vec
            FDR_storey[lower.tri(FDR_storey)] <- t(FDR_storey)[lower.tri(FDR_storey)]
        }
        diag(FDR_storey) <- NA_real_
        FDR_main <- FDR_out_disc
        FDR_main_method <- "BH(exact empirical P; unique unordered off-diagonal pairs)"
        n_sig_005 <- sum(FDR_main[upper.tri(FDR_main)] < 0.05, na.rm = TRUE)
        min_p_possible <- 1 / (perms + 1)
    }

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
                mode = "max",
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
        FDR = FDR_main, # Authoritative global BH on exact empirical P
        FDR_storey = FDR_storey, # Available in in-memory mode only
        FDR_disc = FDR_out_disc, # Original discrete BH
        FDR_beta = FDR_out_beta,
        FDR_mid = FDR_out_mid,
        FDR_uniform = FDR_out_uniform,
        meta = list(
            formula_id = "Lee2009_S2_v1",
            norm_layer = norm_layer,
            input_fingerprint = res$input_fingerprint,
            data_fingerprint = res$input_fingerprint$data,
            permutation_fingerprint = res$input_fingerprint$permutation,
            weight_style = res$weight_style,
            S2 = res$S2,
            perms = perms,
            block_side = block_side,
            effective_block_side = if (use_blocks) block_side else NA_integer_,
            block_size = block_size,
            use_blocks = use_blocks,
            permutation_scheme = if (use_blocks) "spatial_block_joint" else "global_joint_shuffle",
            permutation_backend = permutation_backend,
            ncores = ncores,
            mem_mode = if (inherits(L, "big.matrix")) "bigmemory" else "inmemory",
            p_source = if (is.null(P)) {
                "unavailable"
            } else if (use_blocks) {
                "within-block joint row permutation"
            } else {
                "all-grid joint row permutation"
            },
            z_source = if (!is.null(P)) "signed inverse-normal transform of exact two-sided permutation P" else "unavailable",
            p_resolution = if (!is.null(P)) paste0("~", format(1 / (perms + 1), scientific = TRUE)) else "unavailable",
            pi0_hat = pi0_hat,
            n_tests = if (is.null(FDR_main)) 0L else if (nrow(FDR_main) == ncol(FDR_main)) choose(nrow(FDR_main), 2) else length(FDR_main),
            n_sig_FDR_lt_0.05 = n_sig_005,
            min_p_possible = min_p_possible,
            FDR_main_method = FDR_main_method,
            pval_modes = c("exact=(k+1)/(B+1)", "beta=(k+1)/(B+2)", "mid=(k+0.5)/(B+1)", "uniform=(k+U)/(B+1)"),
            note = "FDR main is BH on exact empirical P over unique unordered off-diagonal pairs; beta/mid/uniform/Storey matrices are exploratory diagnostics."
        )
    )

    if (cache_inputs) {
        if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()
        scope_obj@stats[[gname]][["_lee_input_cache"]] <- input_cache
    }

    step08$done(paste0("layer_name=", layer_name))

    invisible(scope_obj)
}

#' Reapply pointwise confidence-interval width floors
#'
#' The L-vs-r curve applies minimum-width constraints before optional edge
#' smoothing.  Smoothing the two interval boundaries can narrow the resulting
#' interval again, so this helper reapplies the accumulated pointwise floor
#' while preserving each interval's post-processing centre.
#' @param lo Numeric lower confidence bound.
#' @param hi Numeric upper confidence bound.
#' @param floor_width Numeric pointwise minimum interval width.
#' @return A list containing corrected bounds and floor diagnostics.
#' @keywords internal
.reapply_ci_width_floor <- function(lo, hi, floor_width) {
  if (length(lo) != length(hi) || length(lo) != length(floor_width)) {
    stop("lo, hi, and floor_width must have identical lengths", call. = FALSE)
  }

  width_before <- hi - lo
  active <- is.finite(floor_width) & floor_width > 0
  finite_width <- is.finite(width_before)
  reapplied <- active & finite_width & width_before < floor_width
  width_after <- width_before
  width_after[reapplied] <- floor_width[reapplied]

  center <- (lo + hi) / 2
  lo_out <- lo
  hi_out <- hi
  changed <- which(reapplied)
  if (length(changed)) {
    lo_out[changed] <- center[changed] - width_after[changed] / 2
    hi_out[changed] <- center[changed] + width_after[changed] / 2
  }

  shortfall <- rep(0, length(width_before))
  comparable <- active & finite_width
  shortfall[comparable] <- pmax(
    floor_width[comparable] - width_before[comparable],
    0
  )

  list(
    lo = lo_out,
    hi = hi_out,
    floor_width = floor_width,
    active_count = sum(active),
    reapplied = any(reapplied),
    reapplied_count = sum(reapplied),
    reapplied_fraction = if (any(active)) sum(reapplied) / sum(active) else 0,
    max_shortfall_before = if (any(comparable)) max(shortfall[comparable]) else 0
  )
}

#' Fingerprint the aligned inputs used to construct an L-vs-r curve
#'
#' The fingerprint records both gene identity/order and the exact aligned Lee
#' and Pearson matrices.  It therefore distinguishes a change in values from a
#' change in gene ordering, while retaining the upstream Lee input fingerprint.
#' @param Lmat Aligned Lee's L matrix.
#' @param rmat Aligned Pearson correlation matrix.
#' @param common Character vector giving the common genes in analysis order.
#' @param pearson_level Pearson source level (`grid` or `cell`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param lee_input_fingerprint Upstream Lee input fingerprint.
#' @return A nested SHA-256 provenance record.
#' @keywords internal
.lvsr_curve_source_fingerprint <- function(Lmat, rmat, common,
                                          pearson_level,
                                          lee_stats_layer,
                                          lee_input_fingerprint) {
  pearson_level <- match.arg(pearson_level, c("grid", "cell"))
  common <- unname(as.character(common))
  expected_dim <- c(length(common), length(common))
  if (!identical(dim(Lmat), expected_dim) || !identical(dim(rmat), expected_dim)) {
    stop("Lmat and rmat must be square matrices aligned to common", call. = FALSE)
  }
  if (!identical(unname(rownames(Lmat)), common) ||
      !identical(unname(colnames(Lmat)), common) ||
      !identical(unname(rownames(rmat)), common) ||
      !identical(unname(colnames(rmat)), common)) {
    stop("Lmat and rmat dimnames must match common in order", call. = FALSE)
  }

  # Matrix constructors may label the two dimname components (for example,
  # names(dimnames(Lmat)) == c("row", "col")) while otherwise carrying the
  # same ordered genes.  Those container labels have no biological or pair-
  # universe meaning.  Canonicalize them before hashing, but retain the strict
  # row/column identity-and-order checks above.
  canonical_dimnames <- list(common, common)
  dimnames(Lmat) <- canonical_dimnames
  dimnames(rmat) <- canonical_dimnames

  hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
  matrix_payload <- function(x) {
    list(class = class(x), dim = dim(x), dimnames = dimnames(x), values = x)
  }
  ut <- upper.tri(Lmat, diag = FALSE)
  common_sorted <- sort(unique(common), method = "radix")

  list(
    schema = "lvsr_curve_source_fingerprint_v1",
    hash_algorithm = "SHA-256 over R serialization",
    common_genes = list(
      count = length(common),
      order_sha256 = hash(common),
      set_sha256 = hash(common_sorted)
    ),
    pearson = list(
      level = pearson_level,
      dimensions = dim(rmat),
      aligned_matrix_sha256 = hash(matrix_payload(rmat))
    ),
    lee = list(
      stats_layer = lee_stats_layer,
      aligned_matrix_sha256 = hash(matrix_payload(Lmat)),
      upstream_input_fingerprint = lee_input_fingerprint,
      upstream_input_fingerprint_sha256 = hash(lee_input_fingerprint)
    ),
    unordered_pair_universe = list(
      count = sum(ut),
      aligned_L_and_Pearson_sha256 = hash(list(
        LeesL = Lmat[ut],
        Pearson = rmat[ut]
      ))
    )
  )
}

#' Snapshot the R RNG without advancing it
#' @return RNG kind and a SHA-256 fingerprint of `.Random.seed`, when present.
#' @keywords internal
.lvsr_r_rng_snapshot <- function() {
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  state <- if (present) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  list(
    kind = unname(RNGkind()),
    state_present = present,
    state_sha256 = if (present) {
      digest::digest(state, algo = "sha256", serialize = TRUE)
    } else {
      NULL
    }
  )
}

#' Describe the deterministic C++ bootstrap RNG contract
#' @param B Number of bootstrap iterations.
#' @param ncores_requested Thread count requested by the caller.
#' @param ncores_passed Thread count passed to the C++ implementation.
#' @return Bootstrap RNG provenance.
#' @keywords internal
.lvsr_bootstrap_rng_provenance <- function(B, ncores_requested, ncores_passed) {
  list(
    schema = "lvsr_cpp_bootstrap_rng_v1",
    engine = "std::mt19937_64",
    seed_derivation = "SplitMix64 finalizer(base_seed + bootstrap_iteration + golden_gamma)",
    base_seed_hex = "0x00000000B5297A4D",
    golden_gamma_hex = "0x9E3779B97F4A7C15",
    splitmix64_multiplier_1_hex = "0xBF58476D1CE4E5B9",
    splitmix64_multiplier_2_hex = "0x94D049BB133111EB",
    bootstrap_iteration_origin = 0L,
    streams = "one deterministic RNG stream per bootstrap iteration",
    B = as.integer(B),
    ncores_requested = suppressWarnings(as.integer(ncores_requested)[1L]),
    ncores_passed_to_cpp = as.integer(ncores_passed),
    thread_invariant = TRUE,
    thread_invariant_scope = paste(
      "bitwise invariant across thread counts and OpenMP schedules",
      "given identical x, y, strata, grid, and bootstrap parameters"
    ),
    uses_R_rng = FALSE
  )
}

#' Bootstrap Lee's L vs Pearson r relationship
#' @description
#' Internal helper for `.compute_l_vs_r_curve`.
#'
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
                        ncores = max(1L, .detect_cores_safe(logical = TRUE) - 1L),
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
  ncores_requested <- ncores
  ncores <- .clamp_ncores_safe(ncores, logical = TRUE)
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
  lee_layer_obj <- if (!is.null(scope_obj@stats[[grid_name]])) {
    scope_obj@stats[[grid_name]][[lee_stats_layer]]
  } else {
    NULL
  }
  lee_meta <- if (is.list(lee_layer_obj)) lee_layer_obj$meta else NULL
  if (!is.list(lee_meta) ||
      !identical(lee_meta$formula_id, "Lee2009_S2_v1") ||
      !is.list(lee_meta$input_fingerprint) ||
      !identical(lee_meta$input_fingerprint$schema, "lee_input_fingerprint_v2")) {
    stop(
      "LeeStats provenance is missing or incompatible. Rerun computeL() with the canonical Lee's L/S2 implementation before computeLvsRCurve().",
      call. = FALSE
    )
  }
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
  source_fingerprint <- .lvsr_curve_source_fingerprint(
    Lmat = Lmat,
    rmat = rmat,
    common = common,
    pearson_level = level,
    lee_stats_layer = lee_stats_layer,
    lee_input_fingerprint = lee_meta$input_fingerprint
  )
  ut <- upper.tri(Lmat, diag = FALSE)
  Lv <- Lmat[ut]
  rv <- rmat[ut]
  pair_points_total <- length(Lv)
  preprocessing_rng_before <- .lvsr_r_rng_snapshot()
  downsample_keep <- NULL

  # 2. Clean / Downsample
  if (verbose) .log_info(parent, "S02", "cleaning and preprocessing data points", verbose)
  ok <- is.finite(Lv) & is.finite(rv)
  Lv <- Lv[ok]
  rv <- rv[ok]
  if (!length(Lv)) stop("No valid points")
  pair_points_finite <- length(Lv)

  if (verbose) .log_info(parent, "S02", paste0("initial_points=", length(Lv)), verbose)

  if (is.numeric(downsample) && downsample < 1) {
    keep <- sample.int(length(Lv), max(1L, floor(downsample * length(Lv))))
    downsample_keep <- keep
    Lv <- Lv[keep]
    rv <- rv[keep]
    if (verbose) .log_info(parent, "S02", paste0(
      "downsampled_points=", length(Lv), " ratio=", downsample
    ), verbose)
  } else if (is.numeric(downsample) && downsample >= 1 && length(Lv) > downsample) {
    keep <- sample.int(length(Lv), downsample)
    downsample_keep <- keep
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
  pair_points_after_downsample <- length(Lv)
  preprocessing_rng_after <- .lvsr_r_rng_snapshot()

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
  pair_points_after_stratification <- length(Lv)

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
  cpp_ci_type <- 0L
  confidence_level <- 0.95
  k_max_cpp <- if (is.finite(k_max)) as.integer(k_max) else -1L
  bootstrap_rng <- .lvsr_bootstrap_rng_provenance(
    B = B,
    ncores_requested = ncores_requested,
    ncores_passed = ncores
  )
  res <- .loess_residual_bootstrap(
    x = rv, y = Lv, strat = as.integer(strat),
    grid = xgrid,
    B = as.integer(B),
    span = span,
    deg = as.integer(deg),
    n_threads = as.integer(max(1, ncores)),
    k_max = k_max_cpp,
    keep_boot = keep_boot,
    adjust_mode = adjust_mode,
    ci_type = cpp_ci_type,
    level = confidence_level
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

  # Accumulate every mathematical lower bound on CI width.  This is reapplied
  # after edge smoothing, which otherwise can undo either constraint.
  ci_width_floor <- rep(0, length(hi))

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
    ci_width_floor <- pmax(ci_width_floor, target_w)
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
            ci_width_floor[!bad] <- pmax(ci_width_floor[!bad], w_need[!bad])
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

  # Edge smoothing is allowed to regularise the boundaries, but it must not
  # invalidate either the local-MAD or minimum-relative-width lower bound.
  post_smooth_floor <- .reapply_ci_width_floor(lo, hi, ci_width_floor)
  lo <- post_smooth_floor$lo
  hi <- post_smooth_floor$hi

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

  hash_sha256 <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
  downsample_mode <- if (is.null(downsample_keep)) {
    "none"
  } else if (is.numeric(downsample) && downsample < 1) {
    "random_fraction_without_replacement"
  } else {
    "random_maximum_count_without_replacement"
  }
  preprocessing_provenance <- list(
    schema = "lvsr_curve_preprocessing_v1",
    unordered_pair_points = list(
      total = pair_points_total,
      finite = pair_points_finite,
      after_downsample_and_jitter = pair_points_after_downsample,
      after_stratification = pair_points_after_stratification
    ),
    downsample = list(
      requested = downsample,
      mode = downsample_mode,
      selection_index_basis = "finite unordered-pair vector before downsampling",
      selected_indices_sha256 = if (is.null(downsample_keep)) {
        NULL
      } else {
        hash_sha256(as.integer(downsample_keep))
      }
    ),
    jitter = list(
      factor = jitter_eps,
      applied_to = if (jitter_eps > 0) "Pearson pair vector" else "none"
    ),
    R_rng = list(
      used = !is.null(downsample_keep) || jitter_eps > 0,
      functions_called = c(
        if (!is.null(downsample_keep)) "sample.int" else character(0),
        if (jitter_eps > 0) "jitter" else character(0)
      ),
      before = preprocessing_rng_before,
      after = preprocessing_rng_after
    ),
    stratification = list(
      variable = "Pearson_r",
      requested_strata = n_strata,
      effective_strata_limit = n_strata_eff,
      realized_strata = length(unique(strat)),
      breakpoints_sha256 = hash_sha256(brks),
      assignments_sha256 = hash_sha256(as.integer(strat))
    ),
    processed_inputs = list(
      pair_vectors_sha256 = hash_sha256(list(LeesL = Lv, Pearson = rv)),
      xgrid_sha256 = hash_sha256(xgrid)
    )
  )
  bootstrap_settings <- list(
    algorithm = "stratified residual bootstrap with local weighted regression",
    B_requested = as.integer(B),
    B_completed = as.integer(res$B),
    span = span,
    degree = as.integer(deg),
    k_max_requested = k_max,
    k_max_passed_to_cpp = k_max_cpp,
    keep_boot = keep_boot,
    adjust_mode = adjust_mode,
    cpp_ci_type = cpp_ci_type,
    confidence_level = confidence_level,
    ncores_requested = suppressWarnings(as.integer(ncores_requested)[1L]),
    ncores_passed_to_cpp = as.integer(ncores),
    length_out = as.integer(length_out),
    rng = bootstrap_rng
  )
  curve_provenance <- list(
    schema = "lvsr_curve_provenance_v1",
    grid_name = grid_name,
    pearson_level = level,
    lee_stats_layer = lee_stats_layer,
    source = source_fingerprint,
    preprocessing = preprocessing_provenance,
    bootstrap = bootstrap_settings,
    confidence_interval = list(
      method = ci_method,
      adjust = ci_adjust,
      confidence_level = confidence_level,
      min_rel_width = min_rel_width,
      widen_span = widen_span,
      local_mad_floor = "2 * 1.96 * local_MAD",
      edge_smoothing_floor_reapplied = TRUE
    )
  )

  ## ---- meta update ----
  meta_obj <- list(
    schema = "lvsr_curve_meta_v2",
    level = level,
    pearson_level = level,
    lee_stats_layer = lee_stats_layer,
    B = res$B,
    span = span,
    deg = deg,
    ncores = ncores,
    ncores_requested = suppressWarnings(as.integer(ncores_requested)[1L]),
    length_out = length_out,
    downsample = downsample,
    jitter_eps = jitter_eps,
    ci_method = ci_method,
    ci_adjust = ci_adjust,
    confidence_level = confidence_level,
    min_rel_width = min_rel_width,
    widen_span = widen_span,
    n_strata = n_strata,
    n_strata_effective = n_strata_eff,
    n_strata_realized = length(unique(strat)),
    k_max = k_max,
    edf = res$edf,
    sigma2_raw = res$sigma2_raw,
    sigma2_edf = res$sigma2_edf,
    resid_global_mad = res$resid_mad,
    adjust_mode = adjust_mode,
    note = paste(
      "Generated by .compute_l_vs_r_curve with residual-MAD floor + edge smoothing;",
      "all CI width floors are reapplied after smoothing"
    ),
    local_mad_diag = if (exists("local_mad_diag")) local_mad_diag else NULL,
    edge_smooth = if (exists("edge_smooth_info")) edge_smooth_info else NULL,
    post_smooth_floor_reapplied = post_smooth_floor$reapplied,
    post_smooth_floor_reapplied_count = post_smooth_floor$reapplied_count,
    post_smooth_floor_diag = post_smooth_floor[c(
      "active_count", "reapplied", "reapplied_count",
      "reapplied_fraction", "max_shortfall_before"
    )],
    common_genes_fingerprint = source_fingerprint$common_genes,
    pearson_fingerprint = source_fingerprint$pearson,
    source_fingerprint = source_fingerprint,
    preprocessing_provenance = preprocessing_provenance,
    bootstrap_settings = bootstrap_settings,
    bootstrap_rng = bootstrap_rng,
    provenance = curve_provenance
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
#' @param expr_layer Optional expression layer name; when supplied it must match
#'   the `norm_layer` recorded in LeeStats metadata.
#' @param pear_range Parameter value.
#' @param L_range Parameter value.
#' @param do_perm Parameter value.
#' @param block_side Positive integer block side for the current permutation
#'   analysis, used only when `use_blocks = TRUE`; it may differ from the value
#'   used by `computeL()`.
#' @param use_blocks A single logical flag. `FALSE` (the default) applies one
#'   joint unrestricted row shuffle over all grid locations. `TRUE` constrains
#'   that joint shuffle within spatial blocks defined by `block_side`.
#' @param clamp_mode Parameter value.
#' @param p_adj_mode Multiple-testing adjustment mode. `"BH"` (the default)
#'   adjusts the selected top pairs. The `_universe` variants scale by the
#'   eligible finite-pair count after confidence-interval and range filters. No
#'   p-values are computed for unselected pairs.
#' @param mem_limit_GB Positive finite memory budget for permutation batches.
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
                       use_blocks = FALSE,
                       ncores = 1,
                       clamp_mode = c("none", "ref_only", "both"),
                       p_adj_mode = c("BH", "BY", "BH_universe", "BY_universe", "bonferroni"),
                       mem_limit_GB = 2,
                       pval_mode = c("exact", "beta", "mid", "uniform"),
                       curve_layer = NULL,
                       CI_rule = c("remove_within", "remove_outside", "none"),
                       verbose = TRUE) {
  pear_level <- match.arg(pear_level)
  direction <- match.arg(direction)
  p_adj_mode <- match.arg(p_adj_mode)
  clamp_mode <- match.arg(clamp_mode)
  pval_mode <- match.arg(pval_mode)
  CI_rule <- match.arg(CI_rule)
  grid_name <- .lee_assert_layer_name(grid_name, "grid_name")
  lee_stats_layer <- .lee_assert_layer_name(lee_stats_layer, "lee_stats_layer")
  do_perm <- .lee_assert_flag(do_perm, "do_perm")
  use_blocks <- .lee_assert_flag(use_blocks, "use_blocks")
  verbose <- .lee_assert_flag(verbose, "verbose")
  ncores <- .lee_assert_positive_integer(ncores, "ncores")
  ncores_requested <- ncores
  block_side <- .lee_assert_positive_integer(block_side, "block_side")
  mem_limit_GB <- .lee_assert_positive_finite(mem_limit_GB, "mem_limit_GB")
  top_n <- .lee_assert_positive_integer(top_n, "top_n")
  perms <- .lee_assert_nonnegative_integer(perms, "perms")
  if (isTRUE(do_perm) && perms < 1L) {
    stop("When do_perm=TRUE, perms must be a positive integer.", call. = FALSE)
  }
  if (!is.null(expr_layer)) expr_layer <- .lee_assert_layer_name(expr_layer, "expr_layer")
  if (!is.null(curve_layer)) curve_layer <- .lee_assert_layer_name(curve_layer, "curve_layer")
  validate_range <- function(x, name) {
    if (!is.numeric(x) || length(x) != 2L || anyNA(x) || any(!is.finite(x)) || x[1L] > x[2L]) {
      stop(name, " must contain two finite numbers in non-decreasing order.", call. = FALSE)
    }
    as.numeric(x)
  }
  pear_range <- validate_range(pear_range, "pear_range")
  L_range <- validate_range(L_range, "L_range")

  parent <- "getTopLvsR"
  step01 <- .log_step(parent, "S01", "load matrices and compute Delta", verbose)
  step01$enter(paste0("grid_name=", grid_name, " pear_level=", pear_level))
  if (verbose) {
    .log_info(parent, "S01", "starting top Delta analysis", verbose)
    .log_info(parent, "S01", paste0("direction=", direction, " top_n=", top_n), verbose)
  }

  ncores <- .get_top_lvs_r_resolve_runtime_config(ncores, verbose)
  ncores_resolved <- ncores

  ## The observed L matrix and every later permutation must be tied to the
  ## exact normalized data and weights that produced it. Never infer a layer
  ## from the LeeStats name and never fall back to Xz.
  provenance_grid <- .select_grid_layer(scope_obj, grid_name)
  lee_stats_obj <- NULL
  if (!is.null(scope_obj@stats[[grid_name]])) {
    lee_stats_obj <- scope_obj@stats[[grid_name]][[lee_stats_layer]]
  }
  if (is.null(lee_stats_obj)) lee_stats_obj <- provenance_grid[[lee_stats_layer]]
  lee_meta <- if (is.list(lee_stats_obj)) lee_stats_obj$meta else NULL
  if (!is.list(lee_meta) || !identical(lee_meta$formula_id, "Lee2009_S2_v1")) {
    stop(
      "LeeStats layer '", lee_stats_layer,
      "' has missing or incompatible Lee formula provenance. Rerun computeL() with geneSCOPE >= 1.0.1.",
      call. = FALSE
    )
  }
  observed_norm_layer <- lee_meta$norm_layer
  if (!is.character(observed_norm_layer) || length(observed_norm_layer) != 1L ||
      is.na(observed_norm_layer) || !nzchar(observed_norm_layer)) {
    stop(
      "LeeStats metadata has no valid norm_layer. Rerun computeL(); getTopLvsR() will not fall back to Xz.",
      call. = FALSE
    )
  }
  if (!is.null(expr_layer) && !identical(expr_layer, observed_norm_layer)) {
    stop(
      "expr_layer ('", expr_layer, "') does not match the norm_layer ('",
      observed_norm_layer, "') that produced observed Lee's L. Rerun computeL() or omit expr_layer.",
      call. = FALSE
    )
  }
  expr_layer <- observed_norm_layer
  X_current_raw <- provenance_grid[[observed_norm_layer]]
  if (is.null(X_current_raw)) {
    stop(
      "Observed Lee norm_layer '", observed_norm_layer,
      "' is missing from the current grid. Rerun computeL(); fallback to Xz is prohibited.",
      call. = FALSE
    )
  }
  X_current <- .validate_lee_norm_layer(
    X_current_raw, observed_norm_layer,
    context = "[geneSCOPE::getTopLvsR]"
  )
  W_current <- provenance_grid$W
  grid_info_current <- provenance_grid$grid_info
  if (is.null(W_current)) {
    stop("Grid layer missing spatial weights W; rerun computeWeights() and computeL().", call. = FALSE)
  }
  target_ids <- as.character(grid_info_current$grid_id)
  if (!length(target_ids) || anyNA(target_ids) || any(!nzchar(target_ids)) || anyDuplicated(target_ids)) {
    stop("grid_info$grid_id must be complete, non-empty, and unique.", call. = FALSE)
  }
  align_ids <- function(source_ids, label) {
    if (is.null(source_ids) || anyNA(source_ids) || anyDuplicated(source_ids)) {
      stop(label, " grid IDs must be complete and unique.", call. = FALSE)
    }
    idx <- match(target_ids, as.character(source_ids))
    if (anyNA(idx) || length(source_ids) != length(target_ids)) {
      stop(label, " grid ID set does not match grid_info$grid_id.", call. = FALSE)
    }
    idx
  }
  X_current <- X_current[align_ids(rownames(X_current), observed_norm_layer), , drop = FALSE]
  current_weight_style <- .lee_weight_style(W_current)
  W_current <- W_current[
    align_ids(rownames(W_current), "W rows"),
    align_ids(colnames(W_current), "W columns"),
    drop = FALSE
  ]
  attr(W_current, "weight_style") <- current_weight_style
  rownames(X_current) <- target_ids
  dimnames(W_current) <- list(target_ids, target_ids)
  current_fingerprint <- .lee_input_fingerprint(
    X_current, W_current, grid_info_current,
    norm_layer = observed_norm_layer,
    block_side = block_side,
    use_blocks = use_blocks
  )
  if (!is.list(lee_meta$input_fingerprint) ||
      !identical(lee_meta$data_fingerprint, lee_meta$input_fingerprint$data) ||
      !.lee_same_observed_data(lee_meta$input_fingerprint, current_fingerprint)) {
    stop(
      "Current X/W/grid data do not match the inputs that produced observed Lee's L. Rerun computeL() before getTopLvsR().",
      call. = FALSE
    )
  }
  if (!identical(lee_meta$weight_style, current_weight_style)) {
    stop("Current weight_style does not match observed Lee metadata; rerun computeL().", call. = FALSE)
  }
  current_S2 <- .lee_s2_value(W_current)
  if (!is.numeric(lee_meta$S2) || length(lee_meta$S2) != 1L ||
      is.na(lee_meta$S2) || !is.finite(lee_meta$S2) || lee_meta$S2 <= 0 ||
      abs(current_S2 - lee_meta$S2) > 1e-12 * max(1, abs(lee_meta$S2))) {
    stop("Current Lee S2 does not match observed Lee metadata; rerun computeL().", call. = FALSE)
  }

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

  finite_pair <- is.finite(df$LeesL) & is.finite(df$Pear) & is.finite(df$Delta)
  nonfinite_pairs <- sum(!finite_pair)
  if (nonfinite_pairs > 0L && verbose) {
    .log_info(parent, "S01", paste0(
      "excluding_nonfinite_pairs=", nonfinite_pairs
    ), verbose)
  }
  df <- df[finite_pair, , drop = FALSE]
  if (!nrow(df)) {
    stop("No finite Lee's L/Pearson/Delta gene pairs are available.", call. = FALSE)
  }

  # Expression coverage percentages
  {
    expr_pct_map <- setNames(rep(0, length(common)), common)
    g_layer_try <- tryCatch(.select_grid_layer(scope_obj, grid_name), error = function(e) NULL)
    coverage_done <- FALSE

    if (!is.null(g_layer_try) && !is.null(g_layer_try$counts)) {
      ct <- g_layer_try$counts
      # data.table queries can repair/add internal self-reference attributes
      # by reference.  getTopLvsR() is analytically read-only, so detach the
      # counts table before coverage summaries to keep the input scope object
      # byte-stable.
      if (inherits(ct, "data.table")) ct <- data.table::copy(ct)
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
  range_keep <- df$Pear >= pear_range[1] & df$Pear <= pear_range[2] &
    df$LeesL >= L_range[1] & df$LeesL <= L_range[2]
  df <- df[range_keep, , drop = FALSE]
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
      fdr = FDR,
      Delta,
      delta_fdr = FDR
    )
    attr(out, "permutation_provenance") <- list(
      schema = "geneSCOPE_delta_permutation_v1",
      performed = FALSE,
      requested_threads = as.integer(ncores_requested),
      resolved_threads = as.integer(ncores_resolved),
      effective_threads = NA_integer_,
      fallback_occurred = FALSE,
      permutations = as.integer(perms),
      pval_mode = pval_mode,
      p_adj_mode = p_adj_mode,
      total_universe = as.integer(total_universe),
      selected_pairs = as.integer(nrow(out)),
      selection_scope = "top_n_by_observed_delta"
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
  Xz <- X_current
  W <- W_current
  grid_info <- grid_info_current
  genes_top <- unique(c(sel$gene1, sel$gene2))
  gene_map <- match(genes_top, colnames(Xz))
  if (any(is.na(gene_map))) stop("Selected genes not found in Xz")
  gene_pairs <- cbind(
    match(sel$gene1, genes_top) - 1L,
    match(sel$gene2, genes_top) - 1L
  )
  delta_ref <- sel$Delta
  # A common row permutation is applied to every gene, so Pearson r is
  # invariant.  Use exactly the observed reference (including any clamp) in
  # every null Delta instead of recomputing grid-level Pearson in C++.
  pearson_ref <- sel$LeesL - delta_ref
  if (any(!is.finite(pearson_ref)) || any(!is.finite(delta_ref))) {
    stop("Selected Delta/Pearson references must be finite.", call. = FALSE)
  }

  backend_label <- if (use_blocks) "C++ delta_l_fixed_r_perm_csr_block" else "C++ delta_l_fixed_r_perm_csr"
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
  retry_count <- 0L
  successful_batches <- 0L
  while (remaining > 0) {
    bsz <- min(target_batch, remaining)
    success <- FALSE
    while (!success && perm_threads >= 1) {
      if (verbose) .log_info(parent, "S04", sprintf(
        "permutation attempt #%d threads=%d batch=%d remain=%d clamp_mode=%s",
        attempt, perm_threads, bsz, remaining, clamp_mode
      ), verbose)
      attempt <- attempt + 1
      idx_mat <- matrix(NA_integer_, nrow = n_cells, ncol = bsz)
      for (p in seq_len(bsz)) {
        if (use_blocks) {
          idx <- seq_len(n_cells)
          for (rows in split_rows) {
            idx[rows] <- rows[sample.int(length(rows), length(rows), replace = FALSE)]
          }
          if (!identical(sort(idx), seq_len(n_cells)) || any(block_id[idx] != block_id)) {
            stop("Internal block permutation invariant failed", call. = FALSE)
          }
        } else {
          idx <- sample.int(n_cells, n_cells, replace = FALSE)
        }
        idx_mat[, p] <- idx - 1L
      }
      res <- tryCatch(
        {
          if (use_blocks) {
            .delta_l_fixed_r_perm_csr_block(
              Xz_sub, W_indices, W_values, W_row_ptr, idx_mat,
              as.integer(block_id) - 1L, gene_pairs, delta_ref, pearson_ref,
              perm_threads
            )
          } else {
            .delta_l_fixed_r_perm_csr(
              Xz_sub, W_indices, W_values, W_row_ptr, idx_mat,
              gene_pairs, delta_ref, pearson_ref,
              perm_threads
            )
          }
        },
        error = function(e) e
      )
      if (inherits(res, "error")) {
        retry_count <- retry_count + 1L
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
        successful_batches <- successful_batches + 1L
        success <- TRUE
      }
    }
    remaining <- remaining - bsz
  }

  if (verbose) .log_info(parent, "S04", "computing p-values and applying multiple testing correction", verbose)
  N <- perms
  k <- as.integer(exceed_count)
  p_values <- as.numeric(switch(pval_mode,
    exact   = (k + 1) / (N + 1),
    beta    = (k + 1) / (N + 2),
    mid     = (k + 0.5) / (N + 1),
    uniform = (k + runif(length(k))) / (N + 1)
  ))
  p_values[p_values > 1] <- 1
  mc_se <- as.numeric(sqrt(p_values * (1 - p_values) / N))
  p_ci_lo <- as.numeric(qbeta(0.025, k + 1, N - k + 1))
  p_ci_hi <- as.numeric(qbeta(0.975, k + 1, N - k + 1))

  FDR <- as.numeric(switch(p_adj_mode,
    BH          = p.adjust(p_values, "BH"),
    BY          = p.adjust(p_values, "BY"),
    BH_universe = p.adjust(p_values, "BH", n = total_universe),
    BY_universe = p.adjust(p_values, "BY", n = total_universe),
    bonferroni  = p.adjust(p_values, "bonferroni", n = total_universe)
  ))

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
  sel$exceed_count <- as.integer(exceed_count)
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
    fdr = FDR,
    Delta,
    delta_fdr = FDR
  )
  diagnostics <- data.frame(
    gene1 = as.character(sel$gene1),
    gene2 = as.character(sel$gene2),
    exceed_count = as.integer(sel$exceed_count),
    raw_p = as.numeric(sel$raw_p),
    mc_se = as.numeric(sel$mc_se),
    p_ci_lo = as.numeric(sel$p_ci_lo),
    p_ci_hi = as.numeric(sel$p_ci_hi),
    fdr = as.numeric(sel$FDR),
    stringsAsFactors = FALSE
  )
  attr(out, "permutation_diagnostics") <- diagnostics
  attr(out, "permutation_provenance") <- list(
    schema = "geneSCOPE_delta_permutation_v1",
    performed = TRUE,
    requested_threads = as.integer(ncores_requested),
    resolved_threads = as.integer(ncores_resolved),
    effective_threads = as.integer(perm_threads),
    fallback_occurred = isTRUE(retry_count > 0L) ||
      !identical(as.integer(perm_threads), as.integer(ncores_resolved)),
    retry_count = as.integer(retry_count),
    successful_batches = as.integer(successful_batches),
    permutations = as.integer(perms),
    block_side = as.integer(block_side),
    use_blocks = isTRUE(use_blocks),
    pval_mode = pval_mode,
    p_adj_mode = p_adj_mode,
    total_universe = as.integer(total_universe),
    selected_pairs = as.integer(nrow(out)),
    selection_scope = "top_n_by_observed_delta",
    backend = backend_label
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

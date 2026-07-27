#' Adjust Future Globals Maxsize
#' @description
#' Internal helper for `.adjust_future_globals_maxsize`.
#' @param min_bytes Numeric threshold.
#' @return Return value used internally.
#' @keywords internal
.adjust_future_globals_maxsize <- function(min_bytes = 2 * 1024^3) {
    override <- .future_maxsize_env_override()
    target <- if (is.finite(override)) max(min_bytes, override) else min_bytes
    old_value <- getOption("future.globals.maxSize")
    needs_update <- is.null(old_value) || (!is.infinite(old_value) && old_value < target)
    if (needs_update) {
        options(future.globals.maxSize = target)
    }
    function() {
        if (needs_update) {
            options(future.globals.maxSize = old_value)
        }
        invisible(NULL)
    }
}

#' Autoconfigure Networkit Runtime
#' @description
#' Internal helper for `.autoconfigure_networkit_runtime`.
#' @return Return value used internally.
#' @keywords internal
.autoconfigure_networkit_runtime <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))

    python_opt <- getOption(".cluster_genes.python", default = getOption("geneSCOPE.python", default = NULL))
    if (is.character(python_opt) && nzchar(python_opt)) {
        try(.use_python_interpreter(python_opt, required = FALSE), silent = TRUE)
    } else {
        conda_opt <- getOption(".cluster_genes.condaenv", default = getOption("geneSCOPE.condaenv", default = NULL))
        if (is.character(conda_opt) && nzchar(conda_opt)) {
            try(.use_conda_environment(conda_opt, required = FALSE), silent = TRUE)
        }
        venv_opt <- getOption(".cluster_genes.virtualenv", default = getOption("geneSCOPE.virtualenv", default = NULL))
        if (is.character(venv_opt) && nzchar(venv_opt)) {
            try(.use_virtual_environment(venv_opt, required = FALSE), silent = TRUE)
        }
    }

    env_python <- Sys.getenv("NETWORKIT_PYTHON", unset = Sys.getenv("NK_PYTHON", unset = Sys.getenv("RETICULATE_PYTHON", unset = "")))
    if (nzchar(env_python)) {
        try(.use_python_interpreter(env_python, required = FALSE), silent = TRUE)
    }
    invisible(TRUE)
}

#' Build a binary co-membership consensus matrix.
#' @description
#' Internal helper for `.build_membership_consensus`.
#' @param membership Integer vector of cluster labels (NA allowed).
#' @param gene_names Gene identifiers.
#' @return Symmetric sparseMatrix with 1 for within-cluster pairs.
#' @keywords internal
.build_membership_consensus <- function(membership, gene_names) {
    n <- length(membership)
    if (n == 0) {
        return(Matrix(0, 0, 0, sparse = TRUE))
    }
    membership <- as.integer(membership)
    groups <- split(seq_len(n), membership)
    ii <- jj <- integer(0)
    for (idx in groups) {
        idx <- idx[!is.na(idx)]
        s <- length(idx)
        if (s > 1) {
            cmb <- combn(idx, 2)
            ii <- c(ii, cmb[1, ])
            jj <- c(jj, cmb[2, ])
        }
    }
    consensus <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n), symmetric = TRUE)
    diag(consensus) <- 0
    consensus <- drop0(consensus)
    if (!is.null(gene_names) && length(gene_names) == n) {
        dimnames(consensus) <- list(gene_names, gene_names)
    }
    consensus
}

#' Build default clustering configuration.
#' @description
#' Internal helper for `.build_pipeline_config`.
#' @param min_cutoff Numeric Lee's L cutoff.
#' @param use_significance Logical; apply FDR alignment.
#' @param significance_max Maximum allowed FDR value.
#' @param pct_min Quantile filter specification.
#' @param drop_isolated Logical; drop zero-degree nodes.
#' @param algo Stage-1 algorithm name.
#' @param stage2_algo Optional Stage-2 algorithm name.
#' @param resolution CPM resolution parameter.
#' @param gamma Modularity resolution parameter.
#' @param objective Objective name (`"CPM"` or `"modularity"`).
#' @param use_log1p_weight Logical; log1p transform weights.
#' @param use_consensus Logical; enable consensus stage.
#' @param consensus_thr Numeric consensus threshold.
#' @param n_restart Integer number of restarts.
#' @param n_threads Integer thread count.
#' @param mode Character mode flag.
#' @param large_n_threshold Integer node threshold for aggressive mode.
#' @param nk_leiden_iterations Integer NetworKit iteration count.
#' @param nk_leiden_randomize Logical; randomize NetworKit initialisation.
#' @param aggr_future_workers Integer future workers.
#' @param aggr_batch_size Integer restart batch size.
#' @param enable_subcluster Logical; attempt Stage-2 subclustering.
#' @param sub_min_size Minimum Stage-2 cluster size.
#' @param sub_min_child_size Minimum Stage-2 child size.
#' @param sub_resolution_factor Multiplicative Stage-2 resolution factor.
#' @param sub_within_cons_max Maximum within-consensus for subclustering.
#' @param sub_conductance_min Minimum conductance for subclustering.
#' @param sub_improve_within_cons_min Minimum within-consensus gain.
#' @param sub_max_groups Maximum Stage-2 splits.
#' @param enable_qc_filter Logical; enable QC pruning.
#' @param qc_gene_intra_cons_min Minimum within-consensus per gene.
#' @param qc_gene_best_out_cons_min Maximum best-out consensus per gene.
#' @param qc_gene_intra_weight_q Quantile threshold for intra weights.
#' @param keep_cross_stable Logical; retain cross-cluster stable edges.
#' @param min_cluster_size Minimum cluster size to keep.
#' @param CI_rule Character CI95 filtering rule; set `NULL` (default) to skip CI filtering.
#' @param curve_layer Name of curve layer used when CI filtering is enabled.
#' @param similarity_slot Similarity slot name.
#' @param significance_slot Significance slot name.
#' @param use_mh_weight Logical; apply Morisita-Horn weights (when `NULL`, inferred from `mh_slot`/`cmh_slot`).
#' @param mh_slot MH slot name; provide `NULL` to disable MH-based corrections.
#' @param post_smooth Logical; smooth post Stage-2 weights.
#' @param post_smooth_quant Numeric vector of smoothing quantiles.
#' @param post_smooth_power Numeric smoothing power.
#' @param keep_stage1_backbone Logical; keep Stage-1 backbone.
#' @param backbone_floor_q Quantile floor for backbone.
#' @param CI95_filter Parameter value.
#' @param hotspot_k Parameter value.
#' @param hotspot_min_module_size Parameter value.
#' @param future_globals_min_bytes Parameter value.
#' @param system_memory_bytes Parameter value.
#' @return A named configuration list consumed by downstream helpers.
#' @keywords internal
.build_pipeline_config <- function(
    min_cutoff = 0,
    use_significance = TRUE,
    significance_max = 0.05,
    pct_min = "q0",
    drop_isolated = TRUE,
    algo = "leiden",
    stage2_algo = NULL,
    resolution = 1,
    gamma = 1,
    objective = "CPM",
    use_log1p_weight = TRUE,
    use_consensus = TRUE,
    consensus_thr = 0.95,
    n_restart = 100,
    n_threads = NULL,
    mode = "auto",
    large_n_threshold = 1000L,
    nk_leiden_iterations = 10L,
    nk_leiden_randomize = TRUE,
    aggr_future_workers = 2L,
    aggr_batch_size = NULL,
    enable_subcluster = TRUE,
    sub_min_size = 10,
    sub_min_child_size = 5,
    sub_resolution_factor = 1.3,
    sub_within_cons_max = 0.5,
    sub_conductance_min = 0.6,
    sub_improve_within_cons_min = 0.07,
    sub_max_groups = 3,
    enable_qc_filter = TRUE,
    qc_gene_intra_cons_min = 0.25,
    qc_gene_best_out_cons_min = 0.6,
    qc_gene_intra_weight_q = 0.05,
    keep_cross_stable = TRUE,
    min_cluster_size = 2L,
    CI95_filter = NULL,
    CI_rule = NULL,
    curve_layer = NULL,
    similarity_slot = "L",
    significance_slot = "FDR",
    use_mh_weight = FALSE,
    mh_slot = NULL,
    post_smooth = TRUE,
    post_smooth_quant = c(0.05, 0.95),
    post_smooth_power = 0.5,
    keep_stage1_backbone = TRUE,
    backbone_floor_q = 0.02,
    hotspot_k = NULL,
    hotspot_min_module_size = 5L,
    future_globals_min_bytes = 2 * 1024^3,
    system_memory_bytes = NA_real_) {

    algo <- match.arg(algo, c("leiden", "louvain", "hotspot-like", "hotspot"))
    if (identical(algo, "hotspot")) algo <- "hotspot-like"
    if (!is.null(stage2_algo)) {
        stage2_algo <- match.arg(stage2_algo, c("leiden", "louvain", "hotspot-like", "hotspot"))
        if (identical(stage2_algo, "hotspot")) stage2_algo <- "hotspot-like"
    }

    valid_ci_rules <- c("remove_within", "remove_outside")
    if (!is.null(CI_rule)) {
        CI_rule <- match.arg(CI_rule, valid_ci_rules)
    }
    if (is.null(CI_rule) && !is.null(CI95_filter) && isTRUE(CI95_filter)) {
        CI_rule <- "remove_within"
    }
    ci_filter_flag <- !is.null(CI_rule)

    list(
        min_cutoff = min_cutoff,
        use_significance = use_significance,
        significance_max = significance_max,
        pct_min = pct_min,
        drop_isolated = drop_isolated,
        algo = algo,
        stage2_algo = stage2_algo,
        resolution = resolution,
        gamma = gamma,
        objective = objective,
        use_log1p_weight = use_log1p_weight,
        use_consensus = use_consensus,
        consensus_thr = consensus_thr,
        n_restart = n_restart,
        n_threads = n_threads,
        mode = mode,
        large_n_threshold = as.integer(large_n_threshold),
        nk_leiden_iterations = nk_leiden_iterations,
        nk_leiden_randomize = nk_leiden_randomize,
        aggr_future_workers = aggr_future_workers,
        aggr_batch_size = aggr_batch_size,
        enable_subcluster = enable_subcluster,
        sub_min_size = sub_min_size,
        sub_min_child_size = sub_min_child_size,
        sub_resolution_factor = sub_resolution_factor,
        sub_within_cons_max = sub_within_cons_max,
        sub_conductance_min = sub_conductance_min,
        sub_improve_within_cons_min = sub_improve_within_cons_min,
        sub_max_groups = sub_max_groups,
        enable_qc_filter = enable_qc_filter,
        qc_gene_intra_cons_min = qc_gene_intra_cons_min,
        qc_gene_best_out_cons_min = qc_gene_best_out_cons_min,
        qc_gene_intra_weight_q = qc_gene_intra_weight_q,
        keep_cross_stable = keep_cross_stable,
        min_cluster_size = as.integer(min_cluster_size),
        CI95_filter = ci_filter_flag,
        CI_rule = CI_rule,
        curve_layer = curve_layer,
        similarity_slot = similarity_slot,
        significance_slot = significance_slot,
        use_mh_weight = use_mh_weight,
        mh_slot = mh_slot,
        post_smooth = post_smooth,
        post_smooth_quant = post_smooth_quant,
        post_smooth_power = post_smooth_power,
        keep_stage1_backbone = keep_stage1_backbone,
        backbone_floor_q = backbone_floor_q,
        hotspot_k = if (is.null(hotspot_k)) NULL else as.integer(hotspot_k),
        hotspot_min_module_size = as.integer(hotspot_min_module_size),
        future_globals_min_bytes = as.numeric(future_globals_min_bytes),
        system_memory_bytes = as.numeric(system_memory_bytes)
    )
}

#' Cluster Message
#' @description
#' Internal helper for `.cluster_message`.
#' @param verbose Logical; whether to emit progress messages.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.cluster_message <- function(verbose, ..., parent = NULL, step = NULL) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    txt <- paste0(...)
    txt <- sub("^\\[cluster\\]\\s*", "", txt)
    if (!is.null(parent) && !is.null(step)) {
        .log_info(parent, step, txt, verbose)
    }
    invisible(NULL)
}

#' Cluster Metrics
#' @description
#' Internal helper for `.cluster_metrics`.
#' @param members Parameter value.
#' @param cons_mat Parameter value.
#' @param W_mat Parameter value.
#' @return Return value used internally.
#' @keywords internal
.cluster_metrics <- function(members, cons_mat, W_mat) {
    cl_ids <- sort(na.omit(unique(members)))
    res <- lapply(cl_ids, function(cid) {
        idx <- which(members == cid)
        n <- length(idx)
        if (n <= 1) return(data.frame(cluster = cid, size = n, within_cons = NA_real_, conductance = NA_real_))
        C <- as.matrix(cons_mat[idx, idx])
        within_cons <- if (n > 1) mean(C[upper.tri(C)], na.rm = TRUE) else NA_real_
        Win <- sum(W_mat[idx, idx]) / 2
        Wcut <- sum(W_mat[idx, -idx, drop = FALSE])
        conductance <- Wcut / (Wcut + Win + 1e-12)
        data.frame(cluster = cid, size = n, within_cons = within_cons, conductance = conductance)
    })
    do.call(rbind, res)
}

#' Cluster Metrics Sparse
#' @description
#' Internal helper for `.cluster_metrics_sparse`.
#' @param members Parameter value.
#' @param cons_mat Parameter value.
#' @param W_mat Parameter value.
#' @param genes Parameter value.
#' @return Return value used internally.
#' @keywords internal
.cluster_metrics_sparse <- function(members, cons_mat, W_mat, genes) {
    cl_ids <- sort(na.omit(unique(members)))
    if (!inherits(cons_mat, "sparseMatrix")) cons_mat <- as(cons_mat, "dgCMatrix")
    if (!inherits(W_mat, "sparseMatrix")) W_mat <- as(W_mat, "dgCMatrix")
    TTc <- as(cons_mat, "TsparseMatrix")
    maskU_c <- (TTc@i < TTc@j)
    ei_c <- TTc@i[maskU_c] + 1L; ej_c <- TTc@j[maskU_c] + 1L; ex_c <- TTc@x[maskU_c]
    TTw <- as(W_mat, "TsparseMatrix")
    maskU_w <- (TTw@i < TTw@j)
    ei_w <- TTw@i[maskU_w] + 1L; ej_w <- TTw@j[maskU_w] + 1L; ex_w <- TTw@x[maskU_w]
    res <- lapply(cl_ids, function(cid) {
        idx <- which(members == cid)
        n <- length(idx)
        if (n <= 1) return(data.frame(cluster = cid, size = n, within_cons = NA_real_, conductance = NA_real_))
        in_cl <- logical(length(members)); in_cl[idx] <- TRUE
        m_in_c <- in_cl[ei_c] & in_cl[ej_c]
        within_cons <- if (any(m_in_c)) mean(ex_c[m_in_c], na.rm = TRUE) else NA_real_
        m_in_w <- in_cl[ei_w] & in_cl[ej_w]
        m_cut_w <- xor(in_cl[ei_w], in_cl[ej_w])
        Win <- if (any(m_in_w)) sum(ex_w[m_in_w]) else 0
        Wcut <- if (any(m_cut_w)) sum(ex_w[m_cut_w]) else 0
        conductance <- Wcut / (Wcut + Win + 1e-12)
        data.frame(cluster = cid, size = n, within_cons = within_cons, conductance = conductance)
    })
    do.call(rbind, res)
}

#' Cluster Timestamp
#' @description
#' Internal helper for `.cluster_timestamp`.
#' @return Return value used internally.
#' @keywords internal
.cluster_timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

#' Compute consensus matrix in pure R (sparse).
#' @description
#' Internal helper for `.compute_consensus_sparse_r`.
#' @param memb Integer matrix of memberships (genes x runs).
#' @param thr Numeric consensus threshold.
#' @return Symmetric sparseMatrix of mean co-occurrence scores.
#' @keywords internal
.compute_consensus_sparse_r <- function(memb, thr = 0) {
    if (is.null(dim(memb))) memb <- matrix(as.integer(memb), nrow = length(memb))
    if (!is.matrix(memb)) memb <- as.matrix(memb)
    storage.mode(memb) <- "integer"
    n <- nrow(memb); r <- ncol(memb)
    if (!is.finite(r) || r < 1L) return(Matrix(0, n, n, sparse = TRUE))
    acc <- Matrix(0, n, n, sparse = FALSE)
    for (jj in seq_len(r)) {
        labs <- as.integer(memb[, jj])
        keep <- !is.na(labs)
        if (!any(keep)) next
        labs2 <- labs[keep]
        K <- suppressWarnings(max(labs2))
        if (!is.finite(K) || K < 1L) next
        ind <- sparseMatrix(i = which(keep), j = labs2, x = 1, dims = c(n, K))
        acc <- acc + as.matrix(ind %*% t(ind))
    }
    acc <- acc / max(1L, r)
    diag(acc) <- 0
    M <- drop0(Matrix(acc, sparse = TRUE))
    rn <- rownames(memb)
    if (!is.null(rn) && length(rn) == nrow(M)) {
        try({ dimnames(M) <- list(rn, rn) }, silent = TRUE)
    }
    if (thr > 0) {
        TT <- as(M, "TsparseMatrix")
        if (length(TT@x)) {
            keep <- TT@x >= thr
            M <- sparseMatrix(i = TT@i[keep] + 1L, j = TT@j[keep] + 1L, x = TT@x[keep], dims = dim(M), symmetric = TRUE)
            try({ dimnames(M) <- list(rn, rn) }, silent = TRUE)
        }
    }
    M
}

#' Configure Networkit Python
#' @description
#' Internal helper for `.configure_networkit_python`.
#' @param conda_env Parameter value.
#' @param conda_bin Parameter value.
#' @param python_path Filesystem path.
#' @param required Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.configure_networkit_python <- function(conda_env = NULL,
                                        conda_bin = NULL,
                                        python_path = NULL,
                                        required = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))

    if (!is.null(python_path) && nzchar(python_path)) {
        .use_python_interpreter(python_path, required = required)
        return(invisible(TRUE))
    }

    if (!is.null(conda_env) && nzchar(conda_env)) {
        .use_conda_environment(conda_env, conda_bin = conda_bin, required = required)
        return(invisible(TRUE))
    }

    # fallback: try CONDA_PREFIX
    conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
    if (nzchar(conda_prefix)) {
        default_python <- file.path(conda_prefix, "bin", "python")
        try(.use_python_interpreter(default_python, required = FALSE), silent = TRUE)
        return(invisible(TRUE))
    }
    invisible(FALSE)
}

#' Detect System Memory Bytes
#' @description
#' Internal helper for `.detect_system_memory_bytes`.
#' @return Return value used internally.
#' @keywords internal
.detect_system_memory_bytes <- function() {
    total <- NA_real_
    if (requireNamespace("ps", quietly = TRUE)) {
        total <- tryCatch(as.numeric(ps::ps_system_memory()["total"]), error = function(e) NA_real_)
    }
    if (!is.na(total) && is.finite(total) && total > 0) return(total)
    sysname <- tolower(Sys.info()[["sysname"]] %||% "")
    if (identical(.Platform$OS.type, "windows")) {
        lim <- suppressWarnings(memory.limit(size = NA))
        if (is.finite(lim) && lim > 0) {
            total <- lim * 1024^2
        }
    } else if (sysname == "darwin") {
        out <- suppressWarnings(system("/usr/sbin/sysctl -n hw.memsize", intern = TRUE))
        total <- suppressWarnings(as.numeric(out))
    } else {
        if (file.exists("/proc/meminfo")) {
            meminfo <- tryCatch(readLines("/proc/meminfo"), error = function(e) character(0))
            line <- meminfo[grepl("^MemTotal", meminfo)][1]
            if (!is.na(line)) {
                value <- as.numeric(gsub("[^0-9]", "", line))
                total <- value * 1024
            }
        }
    }
    if (!is.na(total) && is.finite(total) && total > 0) total else NA_real_
}

#' Estimate Fast Mode Memory
#' @description
#' Internal helper for `.estimate_fast_mode_memory`.
#' @param similarity_matrix Parameter value.
#' @param config Parameter value.
#' @return Return value used internally.
#' @keywords internal
.estimate_fast_mode_memory <- function(similarity_matrix, config) {
    if (is.null(similarity_matrix)) return(NA_real_)
    base <- suppressWarnings(as.numeric(object.size(similarity_matrix)))
    if (!is.finite(base)) base <- 0
    n <- nrow(similarity_matrix)
    restarts <- as.numeric(config$n_restart %||% 1L)
    if (!is.finite(restarts) || restarts < 1) restarts <- 1
    seed_buf <- as.numeric(n) * restarts * 16
    consensus_buf <- if (isTRUE(config$use_consensus)) {
        as.numeric(n) * as.numeric(n) * 8
    } else 0
    overhead <- base * 0.5
    total <- base + seed_buf + consensus_buf + overhead
    if (!is.finite(total)) NA_real_ else total
}

#' Extract primary matrices from a CosMx scope object.
#' @description
#' Internal helper for `.extract_scope_layers`.
#' @param scope_obj CosMx scope object.
#' @param grid_name Character grid identifier.
#' @param stats_layer Statistics layer name.
#' @param similarity_slot Slot holding the similarity matrix.
#' @param significance_slot Slot holding the significance matrix.
#' @param use_significance Logical; whether to return significance matrix.
#' @param verbose Logical; emit progress messages.
#' @return A list containing similarity and significance matrices plus
#'   metadata handles required by downstream stages.
#' @keywords internal
.extract_scope_layers <- function(scope_obj,
                                        grid_name,
                                        stats_layer = "LeeStats_Xz",
                                        similarity_slot = "L",
                                        significance_slot = "FDR",
                                        use_significance = TRUE,
                                        verbose = TRUE) {
    .cluster_message(verbose, "[cluster] Prep: extracting container layers…",
        parent = "clusterGenes", step = "S03"
    )

    g_layer <- .select_grid_layer(scope_obj, grid_name)
    gname <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    .check_grid_content(scope_obj, gname)

    if (is.null(scope_obj@stats[[gname]]) ||
        is.null(scope_obj@stats[[gname]][[stats_layer]])) {
        stop("Statistics layer missing: ", grid_name, "/", stats_layer)
    }

    stats_obj <- scope_obj@stats[[gname]][[stats_layer]]
    if (identical(similarity_slot, "L")) {
        lee_meta <- stats_obj$meta
        provenance_ok <- is.list(lee_meta) &&
            identical(lee_meta$formula_id, "Lee2009_S2_v1") &&
            is.list(lee_meta$input_fingerprint) &&
            identical(lee_meta$input_fingerprint$schema, "lee_input_fingerprint_v2") &&
            is.numeric(lee_meta$S2) && length(lee_meta$S2) == 1L &&
            is.finite(lee_meta$S2) && lee_meta$S2 > 0
        if (!isTRUE(provenance_ok)) {
            stop(
                "Lee similarity provenance is missing or incompatible in ",
                grid_name, "/", stats_layer,
                ". Rerun computeL() with the canonical Lee's L/S2 implementation before clusterGenes().",
                call. = FALSE
            )
        }
    }
    sim_mat <- stats_obj[[similarity_slot]]
    if (is.null(sim_mat)) stop("Similarity slot missing in stats layer: ", similarity_slot)
    sig_mat <- if (use_significance) stats_obj[[significance_slot]] else NULL
    if (use_significance && is.null(sig_mat)) {
        stop("Significance slot missing in stats layer: ", significance_slot, call. = FALSE)
    }

    list(
        grid_name = gname,
        stats_layer = stats_layer,
        similarity = sim_mat,
        significance = sig_mat,
        stats_object = stats_obj,
        scope_stats = scope_obj@stats[[gname]],
        L_raw = sim_mat,          # backward compatibility
        FDR = sig_mat,            # backward compatibility
        lee_stats = stats_obj     # backward compatibility
    )
}

#' Future Maxsize Env Override
#' @description
#' Internal helper for `.future_maxsize_env_override`.
#' @return Return value used internally.
#' @keywords internal
.future_maxsize_env_override <- function() {
    env_bytes <- suppressWarnings(as.numeric(Sys.getenv("FUTURE_GLOBALS_MAXSIZE_BYTES", "")))
    env_gb <- suppressWarnings(as.numeric(Sys.getenv("FUTURE_GLOBALS_MAXSIZE_GB", "")))
    opt_bytes <- getOption("genescope.future.globals.maxSize", default = NA_real_)
    vals <- c(env_bytes, env_gb * 1024^3, opt_bytes)
    vals <- vals[is.finite(vals) & vals > 0]
    if (!length(vals)) return(NA_real_)
    max(vals, na.rm = TRUE)
}

#' Merge small hotspot clusters into nearest neighbours.
#' @description
#' Internal helper for `.merge_hotspot_clusters`.
#' @param membership Integer vector of cluster assignments.
#' @param similarity_dense Dense similarity matrix.
#' @param min_size Minimum allowed cluster size.
#' @return Updated membership vector.
#' @keywords internal
.merge_hotspot_clusters <- function(membership, similarity_dense, min_size) {
    membership <- as.integer(membership)
    names(membership) <- rownames(similarity_dense)
    repeat {
        freq <- table(membership, useNA = "no")
        small_clusters <- as.integer(names(freq)[freq < min_size])
        if (!length(small_clusters)) break
        small_clusters <- small_clusters[order(as.integer(freq[as.character(small_clusters)]))]
        changed <- FALSE
        for (cid in small_clusters) {
            current_idx <- which(membership == cid)
            other_ids <- setdiff(as.integer(names(freq)), cid)
            if (!length(other_ids)) next
            best_target <- NA_integer_
            best_value <- -Inf
            for (target in other_ids) {
                target_idx <- which(membership == target)
                if (!length(target_idx)) next
                avg_sim <- mean(similarity_dense[current_idx, target_idx, drop = FALSE], na.rm = TRUE)
                if (!is.finite(avg_sim)) avg_sim <- 0
                if (avg_sim > best_value) {
                    best_value <- avg_sim
                    best_target <- target
                }
            }
            if (!is.na(best_target)) {
                membership[current_idx] <- best_target
                changed <- TRUE
            }
        }
        if (!changed) break
    }
    unique_ids <- sort(na.omit(unique(membership)))
    if (length(unique_ids)) {
        id_map <- setNames(seq_along(unique_ids), unique_ids)
        idx <- !is.na(membership)
        membership[idx] <- id_map[as.character(membership[idx])]
    }
    membership
}

#' Prune unsupported or too-small clusters after Stage-2 correction.
#' @description
#' Internal helper for `.prune_unstable_memberships`.
#' @param memberships Integer vector of cluster assignments (NA allowed).
#' @param L_matrix Square matrix of canonical Lee's L weights (gene x gene).
#' @param min_cluster_size Minimum cluster size to retain (default 2).
#' @return A list with elements `membership` (pruned vector) and
#'   `stable_clusters` (character vector of surviving cluster ids).
#' @keywords internal
.prune_unstable_memberships <- function(memberships,
                                      L_matrix,
                                      min_cluster_size = 2L) {
    if (!length(memberships)) {
        return(list(membership = memberships, stable_clusters = character(0)))
    }
    if (is.null(dim(L_matrix)) ||
        nrow(L_matrix) != length(memberships) ||
        ncol(L_matrix) != length(memberships)) {
        stop("`L_matrix` must be square and aligned with `memberships`.", call. = FALSE)
    }
    min_cluster_size <- max(1L, as.integer(min_cluster_size))
    logical_L <- L_matrix > 0
    if (!is.matrix(logical_L)) logical_L <- as.matrix(logical_L)
    same_cluster <- outer(memberships, memberships, "==")
    di <- rowSums(logical_L & same_cluster, na.rm = TRUE)
    memberships[di == 0] <- NA
    tbl <- table(memberships, useNA = "no")
    stable_clusters <- character(0)
    if (length(tbl)) {
        if (min_cluster_size > 1L) {
            small <- names(tbl)[tbl < min_cluster_size]
            if (length(small)) {
                memb_chr <- as.character(memberships)
                memberships[memb_chr %in% small] <- NA
            }
        }
        tbl <- table(memberships, useNA = "no")
        stable_clusters <- names(tbl)[tbl >= min_cluster_size]
    }
    list(membership = memberships, stable_clusters = stable_clusters)
}

#' Perform hotspot-like hierarchical clustering on a similarity matrix.
#' @description
#' Internal helper for `.run_hotspot_clustering`.
#' @param similarity_matrix Sparse similarity matrix (already thresholded).
#' @param cluster_count Desired cluster count.
#' @param min_module_size Minimum module size after merge.
#' @param apply_log1p Whether to apply log1p before normalisation.
#' @return List with `membership` and `consensus` entries.
#' @keywords internal
.run_hotspot_clustering <- function(similarity_matrix,
                                    cluster_count,
                                    min_module_size,
                                    apply_log1p = TRUE) {
    if (cluster_count <= 0) cluster_count <- 1L
    if (!inherits(similarity_matrix, "Matrix") && !is.matrix(similarity_matrix)) {
        similarity_matrix <- as.matrix(similarity_matrix)
    }
    matrix_work <- similarity_matrix
    if (apply_log1p && inherits(matrix_work, "sparseMatrix")) {
        matrix_work <- as(matrix_work, "dgCMatrix")
        if (length(matrix_work@x)) matrix_work@x <- log1p(matrix_work@x)
    } else if (apply_log1p) {
        matrix_work <- log1p(matrix_work)
    }
    dense_matrix <- if (inherits(matrix_work, "Matrix")) as.matrix(matrix_work) else matrix_work
    if (!is.double(dense_matrix)) dense_matrix <- matrix(as.numeric(dense_matrix), nrow = nrow(dense_matrix), dimnames = dimnames(dense_matrix))
    max_val <- suppressWarnings(max(dense_matrix, na.rm = TRUE))
    if (!is.finite(max_val) || max_val <= 0) max_val <- 1
    dense_matrix <- dense_matrix / max_val
    diag(dense_matrix) <- 1
    gene_names <- rownames(dense_matrix)
    cluster_count <- max(1L, min(as.integer(cluster_count), nrow(dense_matrix)))
    dist_matrix <- as.dist(1 - dense_matrix)
    hc <- try(if (requireNamespace("fastcluster", quietly = TRUE)) fastcluster::hclust(dist_matrix, method = "ward.D2") else hclust(dist_matrix, method = "ward.D2"), silent = TRUE)
    if (inherits(hc, "try-error")) hc <- hclust(dist_matrix, method = "average")
    membership <- cutree(hc, k = cluster_count)
    if (is.null(names(membership))) names(membership) <- gene_names
    if (!is.null(min_module_size) && is.finite(min_module_size) && min_module_size > 1) {
        membership <- .merge_hotspot_clusters(membership, dense_matrix, as.integer(min_module_size))
    }
    consensus <- .build_membership_consensus(membership, gene_names)
    list(membership = membership, consensus = consensus)
}

#' Run Single Partition
#' @description
#' Internal helper for `.run_single_partition`.
#' @param graph Parameter value.
#' @param backend Parameter value.
#' @param algo_name Parameter value.
#' @param res_param Parameter value.
#' @param objective Parameter value.
#' @param threads Number of threads to use.
#' @param nk_opts Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_single_partition <- function(graph,
                                  backend = c("igraph", "networkit"),
                                  algo_name = c("leiden", "louvain"),
                                  res_param,
                                  objective,
                                  threads = 1L,
                                  nk_opts = list()) {
    backend <- match.arg(backend)
    algo_name <- match.arg(algo_name)
    vertex_names <- igraph::V(graph)$name

    if (backend == "networkit") {
        nk_inputs <- .nk_prepare_graph_input(graph)
        graph <- nk_inputs$graph
        edge_df <- nk_inputs$edge_table
        weights <- nk_inputs$edge_weights
        vertex_names <- nk_inputs$vertex_names
        seeds <- nk_opts$seeds %||% 1L
        if (identical(algo_name, "leiden")) {
            memb_mat <- .networkit_parallel_leiden_runs(
                vertex_names = vertex_names,
                edge_table = edge_df,
                edge_weights = weights,
                gamma = as.numeric(res_param),
                iterations = nk_opts$iterations %||% 10L,
                randomize = nk_opts$randomize %||% TRUE,
                threads = threads,
                seeds = as.integer(seeds)
            )
        } else {
            memb_mat <- .networkit_plm_runs(
                vertex_names = vertex_names,
                edge_table = edge_df,
                edge_weights = weights,
                gamma = as.numeric(res_param),
                refine = nk_opts$refine %||% TRUE,
                threads = threads,
                seeds = as.integer(seeds)
            )
        }
        membership <- as.integer(memb_mat[, 1])
        names(membership) <- vertex_names
        return(membership)
    }

    clustering <- .run_graph_algorithm(graph, algo_name, res_param, objective)
    membership <- igraph::membership(clustering)
    if (is.null(names(membership))) names(membership) <- vertex_names
    membership
}

#' Select Grid Layer
#' @description
#' Internal helper for `.select_grid_layer`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.select_grid_layer <- function(scope_obj, grid_name = NULL, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    if (is.null(scope_obj@grid) || length(scope_obj@grid) == 0) {
        stop("No grid layers found in scope_obj")
    }

    if (is.null(grid_name)) {
        if (length(scope_obj@grid) == 1) {
            # Removed verbose message to avoid redundancy
            return(scope_obj@grid[[1]])
        } else {
            stop("Multiple grid layers found. Please specify grid_name.")
        }
    }

    if (!grid_name %in% names(scope_obj@grid)) {
        stop("Grid layer '", grid_name, "' not found.")
    }

    # Removed verbose message to avoid redundancy with main functions
    return(scope_obj@grid[[grid_name]])
}

#' Stage1 Consensus Workflow V2
#' @description
#' Internal helper for `.stage1_consensus_workflow_v2`.
#' @param similarity_matrix Parameter value.
#' @param threshold Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param networkit_config Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage1_consensus_workflow_v2 <- function(similarity_matrix,
                                      threshold,
                                      config,
                                      verbose = TRUE,
                                      networkit_config = list()) {
    .cluster_message(verbose, "[cluster] Step2: running Stage-1 consensus...",
        parent = "clusterGenes", step = "S05"
    )

    validated <- .stage1_validate_inputs(threshold, config, networkit_config)
    kept_genes <- validated$kept_genes
    similarity_sub <- validated$similarity_sub
    genes_all <- validated$genes_all

    list2env(
        .stage1_prepare_consensus_inputs(similarity_sub, kept_genes, config),
        envir = environment()
    )

    runtime_cfg <- .stage1_configure_runtime_threads(config, kept_genes, verbose)
    if (!is.null(runtime_cfg$future_guard)) {
        on.exit(runtime_cfg$future_guard(), add = TRUE)
    }

    run_res <- .stage1_run_consensus_and_collect(
        similarity_sub = similarity_sub,
        nk_inputs_stage1 = nk_inputs_stage1,
        similarity_graph = similarity_graph,
        runtime_cfg = runtime_cfg,
        config = config,
        verbose = verbose
    )
    stage1_membership_labels <- run_res$stage1_membership_labels
    stage1_consensus_matrix <- run_res$stage1_consensus_matrix
    restart_memberships <- run_res$restart_memberships
    stage1_backend <- run_res$stage1_backend
    stage1_algo_per_run <- run_res$stage1_algo_per_run

    .stage1_postprocess_and_assemble(
        kept_genes = kept_genes,
        genes_all = genes_all,
        threshold = threshold,
        similarity_sub = similarity_sub,
        similarity_graph = similarity_graph,
        mode_selected = runtime_cfg$mode_selected,
        restart_memberships = restart_memberships,
        stage1_consensus_matrix = stage1_consensus_matrix,
        stage1_membership_labels = stage1_membership_labels,
        stage1_backend = stage1_backend,
        stage1_algo = runtime_cfg$algo,
        stage1_algo_per_run = stage1_algo_per_run
    )
}

#' Stage1 Validate Inputs
#' @description
#' Internal helper for `.stage1_validate_inputs`.
#' @param threshold Parameter value.
#' @param config Parameter value.
#' @param networkit_config Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage1_validate_inputs <- function(threshold, config, networkit_config) {
    kept_genes <- threshold$kept_genes
    similarity_sub <- threshold$similarity_sub
    genes_all <- threshold$genes_all

    if (!is.null(networkit_config) && length(networkit_config)) {
        config_conda_env <- networkit_config[["conda_env"]]
        if (is.null(config_conda_env)) config_conda_env <- networkit_config[["nk_condaenv"]]
        config_conda_bin <- networkit_config[["conda_bin"]]
        if (is.null(config_conda_bin)) config_conda_bin <- networkit_config[["nk_conda_bin"]]
        config_python_path <- networkit_config[["python_path"]]
        if (is.null(config_python_path)) config_python_path <- networkit_config[["nk_python"]]
        try(.configure_networkit_python(
            conda_env = config_conda_env,
            conda_bin = config_conda_bin,
            python_path = config_python_path,
            required = FALSE
        ), silent = TRUE)
    }

    list(
        kept_genes = kept_genes,
        similarity_sub = similarity_sub,
        genes_all = genes_all
    )
}

#' Stage1 Prepare Consensus Inputs
#' @description
#' Internal helper for `.stage1_prepare_consensus_inputs`.
#' @param similarity_sub Parameter value.
#' @param kept_genes Parameter value.
#' @param config Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage1_prepare_consensus_inputs <- function(similarity_sub, kept_genes, config) {
    adjacency_triplet <- as.data.frame(summary(similarity_sub))
    if (!all(c("i", "j", "x") %in% colnames(adjacency_triplet))) {
        colnames(adjacency_triplet) <- c("i", "j", "x")[seq_len(ncol(adjacency_triplet))]
    }
    adjacency_triplet <- adjacency_triplet[adjacency_triplet$i < adjacency_triplet$j & adjacency_triplet$x > 0, , drop = FALSE]
    if (!nrow(adjacency_triplet)) stop("All edges filtered out at Stage-1.")

    edge_table <- data.frame(
        from = kept_genes[adjacency_triplet$i],
        to = kept_genes[adjacency_triplet$j],
        sim_edge = adjacency_triplet$x,
        stringsAsFactors = FALSE
    )
    edge_table$weight <- if (config$use_log1p_weight) log1p(edge_table$sim_edge) else edge_table$sim_edge
    similarity_graph <- igraph::graph_from_data_frame(edge_table[, c("from", "to", "weight")],
        directed = FALSE, vertices = kept_genes
    )

    nk_inputs_stage1 <- .nk_prepare_graph_input(similarity_graph)
    similarity_graph <- nk_inputs_stage1$graph

    list(
        similarity_graph = similarity_graph,
        nk_inputs_stage1 = nk_inputs_stage1
    )
}

#' Stage1 Configure Runtime Threads
#' @description
#' Internal helper for `.stage1_configure_runtime_threads`.
#' @param config Parameter value.
#' @param kept_genes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.stage1_configure_runtime_threads <- function(config, kept_genes, verbose) {
    algo <- match.arg(config$algo, c("leiden", "louvain", "hotspot-like"))
    objective <- match.arg(config$objective, c("CPM", "modularity"))
    mode_selected <- switch(config$mode,
        fast = "fast",
        auto = if (isTRUE(config$prefer_fast)) "fast" else if (length(kept_genes) > config$large_n_threshold) "aggressive" else "safe",
        safe = "safe",
        aggressive = "aggressive"
    )
    backend_mode <- if (identical(mode_selected, "fast")) "safe" else mode_selected
    .cluster_message(verbose, sprintf("[cluster]   Stage-1 mode=%s (nodes=%d)", mode_selected, length(kept_genes)),
        parent = "clusterGenes", step = "S05"
    )
    if (identical(backend_mode, "aggressive") && identical(objective, "CPM")) {
        warning("[cluster] CPM objective not supported in aggressive mode; falling back to 'modularity'.")
        objective <- "modularity"
    }

    res_param <- if (identical(objective, "modularity")) as.numeric(config$gamma) else as.numeric(config$resolution)
    algo_per_run <- if (identical(backend_mode, "aggressive") && identical(algo, "leiden")) "louvain" else algo
    n_threads <- .clamp_ncores_safe(config$n_threads %||% 1L, logical = TRUE)
    future_guard <- .adjust_future_globals_maxsize(
        config$future_globals_min_bytes %||% (2 * 1024^3)
    )

    list(
        algo = algo,
        objective = objective,
        mode_selected = mode_selected,
        backend_mode = backend_mode,
        res_param = res_param,
        algo_per_run = algo_per_run,
        n_threads = n_threads,
        future_guard = future_guard
    )
}

#' Stage1 Run Consensus And Collect
#' @description
#' Internal helper for `.stage1_run_consensus_and_collect`.
#' @param similarity_sub Parameter value.
#' @param nk_inputs_stage1 Parameter value.
#' @param similarity_graph Parameter value.
#' @param runtime_cfg Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.stage1_run_consensus_and_collect <- function(similarity_sub,
                                              nk_inputs_stage1,
                                              similarity_graph,
                                              runtime_cfg,
                                              config,
                                              verbose) {
    algo <- runtime_cfg$algo
    objective <- runtime_cfg$objective
    backend_mode <- runtime_cfg$backend_mode
    res_param <- runtime_cfg$res_param
    algo_per_run <- runtime_cfg$algo_per_run
    n_threads <- runtime_cfg$n_threads

    kept_genes <- igraph::V(similarity_graph)$name
    stage1_consensus_matrix <- NULL
    stage1_membership_labels <- NULL
    restart_memberships <- NULL
    stage1_backend <- "igraph"
    stage1_backend <- "igraph"

    if (identical(algo, "hotspot-like")) {
        .log_backend("clusterGenes", "S05", "stage1_backend", "hotspot-like",
            reason = "algo=hotspot-like", verbose = verbose
        )
        hotspot_k <- config$hotspot_k
        if (is.null(hotspot_k) || !is.finite(hotspot_k)) {
            hotspot_k <- max(1L, min(20L, length(kept_genes)))
        }
        list2env(
            .run_hotspot_clustering(similarity_sub, hotspot_k, config$hotspot_min_module_size, config$use_log1p_weight),
            envir = environment()
        )
        stage1_membership_labels <- membership
        stage1_consensus_matrix <- consensus
        restart_memberships <- matrix(as.integer(stage1_membership_labels), ncol = 1,
            dimnames = list(kept_genes, "hotspot"))
        stage1_backend <- "hotspot"
    } else {
        edge_weights <- nk_inputs_stage1$edge_weights
        vertex_names <- nk_inputs_stage1$vertex_names
        edge_pairs <- nk_inputs_stage1$edge_table

        aggressive_requested <- identical(backend_mode, "aggressive")
        if (aggressive_requested && !.networkit_available()) {
            stop("Aggressive mode requires the NetworKit backend; please configure reticulate/networkit before running clusterGenes_new.", call. = FALSE)
        }
        use_networkit_leiden <- aggressive_requested && identical(algo_per_run, "leiden")
        use_networkit_plm <- aggressive_requested && identical(algo_per_run, "louvain")

        if (use_networkit_leiden || use_networkit_plm) {
            backend_label <- if (use_networkit_leiden) "Leiden" else "PLM/louvain"
            .log_backend("clusterGenes", "S05", "stage1_backend", "networkit",
                reason = paste0("algo=", backend_label), verbose = verbose
            )
            stage1_backend <- "networkit"
            seeds <- seq_len(config$n_restart)
            chunk_size <- config$aggr_batch_size
            if (is.null(chunk_size) || chunk_size <= 0) {
                chunk_size <- ceiling(config$n_restart / max(1L, as.integer(config$aggr_future_workers)))
            }
            split_chunks <- function(x, k) {
                if (k <= 0) return(list(x))
                split(x, ceiling(seq_along(x) / k))
            }
            chunks <- split_chunks(seeds, chunk_size)
            threads_per_worker <- max(1L, floor(n_threads / max(1L, as.integer(config$aggr_future_workers))))

            worker_fun <- function(ss) {
                if (use_networkit_leiden) {
                    .networkit_parallel_leiden_runs(
                        vertex_names = vertex_names,
                        edge_table = edge_pairs,
                        edge_weights = edge_weights,
                        gamma = as.numeric(res_param),
                        iterations = as.integer(config$nk_leiden_iterations),
                        randomize = isTRUE(config$nk_leiden_randomize),
                        threads = threads_per_worker,
                        seeds = as.integer(ss)
                    )
                } else {
                    .networkit_plm_runs(
                        vertex_names = vertex_names,
                        edge_table = edge_pairs,
                        edge_weights = edge_weights,
                        gamma = as.numeric(res_param),
                        refine = TRUE,
                        threads = threads_per_worker,
                        seeds = as.integer(ss)
                    )
                }
            }
            membership_blocks <- if (config$aggr_future_workers > 1L) {
                future.apply::future_lapply(chunks, worker_fun, future.seed = TRUE)
            } else {
                lapply(chunks, worker_fun)
            }
            restart_memberships <- do.call(cbind, membership_blocks)
        } else {
            .log_backend("clusterGenes", "S05", "stage1_backend", "igraph",
                reason = paste0("algo=", algo_per_run), verbose = verbose
            )
            # Aggressive mode no longer falls back to igraph; the original igraph fallback is retained below for reference.
            restart_memberships <- future.apply::future_sapply(
                seq_len(config$n_restart),
                function(ii) {
                    g_local <- igraph::graph_from_data_frame(
                        cbind(edge_pairs, weight = edge_weights),
                        directed = FALSE,
                        vertices = vertex_names
                    )
                    res <- try(.run_graph_algorithm(g_local, algo_per_run, res_param, objective), silent = TRUE)
                    if (inherits(res, "try-error")) {
                        stop("Stage-1 graph algorithm failed: ", as.character(res), call. = FALSE)
                    } else {
                        mm <- try(igraph::membership(res), silent = TRUE)
                        if (inherits(mm, "try-error") || is.null(mm)) {
                            stop("Stage-1 membership extraction failed: ", as.character(mm), call. = FALSE)
                        }
                        mm
                    }
                },
                future.seed = TRUE
            )
        }
        if (is.null(dim(restart_memberships))) restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
        if (!is.matrix(restart_memberships) && !inherits(restart_memberships, "Matrix")) {
            restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
        }
        if (is.null(nrow(restart_memberships)) || nrow(restart_memberships) != length(kept_genes)) {
            len <- length(restart_memberships)
            if (len %% length(kept_genes) == 0) {
                restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
            } else {
                stop("Stage-1 membership shape mismatch: ", len, " elements vs ", length(kept_genes))
            }
        }
        rownames(restart_memberships) <- kept_genes

        if (!config$use_consensus) {
            stage1_membership_labels <- restart_memberships[, 1]
            stage1_consensus_matrix <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                dims = c(length(kept_genes), length(kept_genes)),
                dimnames = list(kept_genes, kept_genes))
        } else {
            stage1_consensus_matrix <- try(consensus_sparse(restart_memberships,
                thr = config$consensus_thr,
                n_threads = n_threads
            ), silent = TRUE)
            if (!inherits(stage1_consensus_matrix, "try-error")) {
                .log_backend("clusterGenes", "S05", "stage1_consensus_backend", "cpp",
                    reason = "consensus_sparse", verbose = verbose
                )
            } else {
                .log_backend("clusterGenes", "S05", "stage1_consensus_backend", "r",
                    reason = "consensus_sparse_failed", verbose = verbose
                )
                stage1_consensus_matrix <- .compute_consensus_sparse_r(restart_memberships, thr = config$consensus_thr)
            }
            try({ if (length(stage1_consensus_matrix@i)) diag(stage1_consensus_matrix) <- 0 }, silent = TRUE)
            if (nnzero(stage1_consensus_matrix) == 0) {
                stage1_membership_labels <- restart_memberships[, 1]
            } else {
                consensus_graph_view <- igraph::graph_from_adjacency_matrix(stage1_consensus_matrix,
                    mode = "undirected", weighted = TRUE, diag = FALSE)
                consensus_backend <- if (identical(stage1_backend, "networkit")) "networkit" else "igraph"
                backend_label <- if (identical(consensus_backend, "networkit")) {
                    paste0("NetworKit (", algo, ")")
                } else {
                    paste0("igraph (", algo, ")")
                }
                .log_backend("clusterGenes", "S05", "stage1_consensus_graph_backend", backend_label,
                    reason = paste0("stage1_backend=", stage1_backend), verbose = verbose
                )
                nk_opts <- list(
                    iterations = config$nk_leiden_iterations,
                    randomize = config$nk_leiden_randomize,
                    refine = TRUE,
                    seeds = 1L
                )
                stage1_membership_labels <- .run_single_partition(
                    consensus_graph_view,
                    backend = consensus_backend,
                    algo_name = algo,
                    res_param = res_param,
                    objective = objective,
                    threads = n_threads,
                    nk_opts = nk_opts
                )
            }
        }
    }

    list(
        stage1_membership_labels = stage1_membership_labels,
        stage1_consensus_matrix = stage1_consensus_matrix,
        restart_memberships = restart_memberships,
        stage1_backend = stage1_backend,
        stage1_algo_per_run = algo_per_run
    )
}

#' Stage1 Postprocess And Assemble
#' @description
#' Internal helper for `.stage1_postprocess_and_assemble`.
#' @param kept_genes Parameter value.
#' @param genes_all Parameter value.
#' @param threshold Parameter value.
#' @param similarity_sub Parameter value.
#' @param similarity_graph Parameter value.
#' @param mode_selected Parameter value.
#' @param restart_memberships Parameter value.
#' @param stage1_consensus_matrix Parameter value.
#' @param stage1_membership_labels Parameter value.
#' @param stage1_backend Parameter value.
#' @param stage1_algo Parameter value.
#' @param stage1_algo_per_run Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage1_postprocess_and_assemble <- function(kept_genes,
                                             genes_all,
                                             threshold,
                                             similarity_sub,
                                             similarity_graph,
                                             mode_selected,
                                             restart_memberships,
                                             stage1_consensus_matrix,
                                             stage1_membership_labels,
                                             stage1_backend,
                                             stage1_algo,
                                             stage1_algo_per_run) {
    stage1_edge_triplet <- summary(stage1_consensus_matrix)
    if (!all(c("i", "j", "x") %in% colnames(stage1_edge_triplet))) {
        colnames(stage1_edge_triplet) <- c("i", "j", "x")[seq_len(ncol(stage1_edge_triplet))]
    }
    stage1_edge_triplet <- stage1_edge_triplet[stage1_edge_triplet$i < stage1_edge_triplet$j & stage1_edge_triplet$x > 0, , drop = FALSE]
    stage1_cons_graph <- NULL
    if (nrow(stage1_edge_triplet)) {
        consensus_edge_df <- data.frame(
            from = kept_genes[stage1_edge_triplet$i],
            to = kept_genes[stage1_edge_triplet$j],
            weight = stage1_edge_triplet$x,
            stringsAsFactors = FALSE
        )
        stage1_cons_graph <- igraph::graph_from_data_frame(
            consensus_edge_df,
            directed = FALSE,
            vertices = kept_genes
        )
    } else {
        stage1_cons_graph <- igraph::make_empty_graph(n = length(kept_genes), directed = FALSE)
        stage1_cons_graph <- igraph::set_vertex_attr(stage1_cons_graph, "name", value = kept_genes)
    }

    list(
        kept_genes = kept_genes,
        genes_all = genes_all,
        similarity_full = threshold$similarity,
        similarity_sub = similarity_sub,
        A = threshold$similarity,          # backward compatibility
        A_sub = similarity_sub,            # backward compatibility
        g_stage1 = similarity_graph,
        mode_selected = mode_selected,
        stage1_membership_matrix = restart_memberships,
        stage1_consensus = stage1_consensus_matrix,
        stage1_membership = stage1_membership_labels,
        stage1_consensus_graph = stage1_cons_graph,
        stage1_backend = stage1_backend,
        stage1_algo = stage1_algo,
        stage1_algo_per_run = stage1_algo_per_run %||% stage1_algo
    )
}


#' Stage2 Validate Inputs
#' @description
#' Internal helper for `.stage2_validate_inputs`.
#' @param stage1 Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.stage2_validate_inputs <- function(stage1, config, verbose) {
    kept_genes <- stage1$kept_genes
    genes_all <- stage1$genes_all
    stage1_membership_labels <- stage1$stage1_membership
    stage1_consensus_matrix <- stage1$stage1_consensus
    stage1_similarity_graph <- stage1$g_stage1
    stage1_mode_selected <- stage1$mode_selected
    if (is.null(stage1_mode_selected)) stage1_mode_selected <- "safe"
    stage1_membership_matrix <- stage1$stage1_membership_matrix
    stage1_backend <- stage1$stage1_backend %||% "igraph"
    stage1_algo_effective <- stage1$stage1_algo_per_run %||% stage1$stage1_algo %||% config$algo

    stage2_algo_final <- tryCatch({
        if (is.null(config$stage2_algo)) config$algo else match.arg(config$stage2_algo, c("leiden", "louvain", "hotspot-like"))
    }, error = function(e) config$algo)
    .cluster_message(verbose, "[cluster]   Stage-2 base algorithm: ", stage2_algo_final,
        parent = "clusterGenes", step = "S07"
    )
    stage2_backend <- if (!is.null(config$stage2_algo)) {
        if (identical(stage2_algo_final, "hotspot-like")) "hotspot" else "igraph"
    } else {
        if (identical(stage2_algo_final, "hotspot-like")) {
            "hotspot"
        } else if (identical(stage1_backend, "networkit")) {
            "networkit"
        } else {
            "igraph"
        }
    }
    if (isTRUE(config$prefer_fast) && !identical(stage2_backend, "hotspot")) {
        stage2_backend <- "igraph"
    }
    reason <- if (!is.null(config$stage2_algo)) {
        paste0("stage2_algo=", stage2_algo_final)
    } else {
        paste0("stage1_backend=", stage1_backend)
    }
    .log_backend("clusterGenes", "S07", "stage2_backend", stage2_backend,
        reason = reason, verbose = verbose
    )

    list(
        kept_genes = kept_genes,
        genes_all = genes_all,
        stage1_membership_labels = stage1_membership_labels,
        stage1_consensus_matrix = stage1_consensus_matrix,
        stage1_similarity_graph = stage1_similarity_graph,
        stage1_mode_selected = stage1_mode_selected,
        stage1_membership_matrix = stage1_membership_matrix,
        stage1_backend = stage1_backend,
        stage1_algo_effective = stage1_algo_effective,
        stage2_algo_final = stage2_algo_final,
        stage2_backend = stage2_backend
    )
}

#' Stage2 Prepare Paths And Io
#' @description
#' Internal helper for `.stage2_prepare_paths_and_io`.
#' @param stage1_membership_labels Parameter value.
#' @param kept_genes Parameter value.
#' @param stage1_consensus_matrix Parameter value.
#' @param stage1 Parameter value.
#' @param genes_all Parameter value.
#' @param stage1_similarity_graph Parameter value.
#' @param similarity_matrix Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.stage2_prepare_paths_and_io <- function(stage1_membership_labels,
                                         kept_genes,
                                         stage1_consensus_matrix,
                                         stage1,
                                         genes_all,
                                         stage1_similarity_graph,
                                         similarity_matrix,
                                         config,
                                         verbose) {
    if (!config$use_mh_weight && !config$CI95_filter) {
        .log_info("clusterGenes", "S07", "no Stage-2 corrections requested; returning Stage-1 result", verbose)

        # --- Stage-1 membership cleanup (keeps behaviour identical to old pipeline) ----
        if (length(stage1_membership_labels) && length(kept_genes)) {
            stage1_mem <- stage1_membership_labels
            if (!is.null(names(stage1_mem))) {
                stage1_mem <- stage1_mem[kept_genes]
            } else if (length(stage1_mem) == length(kept_genes)) {
                names(stage1_mem) <- kept_genes
            } else {
                stage1_mem <- setNames(rep(NA_integer_, length(kept_genes)), kept_genes)
            }
            stage1_mem <- as.integer(stage1_mem)
            names(stage1_mem) <- kept_genes
            sim_stage1 <- stage1$similarity_sub %||% stage1$A_sub
            if (!is.null(sim_stage1)) {
                if (!identical(rownames(sim_stage1), kept_genes)) {
                    sim_stage1 <- sim_stage1[kept_genes, kept_genes, drop = FALSE]
                }
                min_cluster_keep <- max(1L, as.integer(config$min_cluster_size %||% 2L))
                pruned <- .prune_unstable_memberships(
                    memberships = stage1_mem,
                    L_matrix = sim_stage1,
                    min_cluster_size = min_cluster_keep
                )
                stage1_mem <- pruned$membership
                idx_valid <- which(!is.na(stage1_mem))
                if (length(idx_valid)) {
                    unique_ids <- sort(unique(stage1_mem[idx_valid]))
                    id_map <- setNames(seq_along(unique_ids), unique_ids)
                    stage1_mem[idx_valid] <- as.integer(id_map[as.character(stage1_mem[idx_valid])])
                }
                names(stage1_mem) <- kept_genes
                stage1_membership_labels <- stage1_mem
            }
        }

        # Rebuild only from the already FDR/q95-filtered Stage-1 matrix.
        A_sub_local <- stage1$similarity_sub %||% stage1$A_sub
        if (is.null(A_sub_local)) {
            stop("Stage-1 did not retain its filtered similarity matrix.", call. = FALSE)
        }
        if (!identical(rownames(A_sub_local), kept_genes)) {
            A_sub_local <- A_sub_local[kept_genes, kept_genes, drop = FALSE]
        }
        assigned_genes <- kept_genes[!is.na(stage1_membership_labels[kept_genes])]
        A_final <- A_sub_local[assigned_genes, assigned_genes, drop = FALSE]
        diag(A_final) <- 0
        edges_keep <- .arrind_from_matrix_predicate(
            A_final, op = "gt", cutoff = 0,
            triangle = "upper", keep_diag = FALSE
        )
        if (nrow(edges_keep)) {
            base_weight <- as.numeric(A_final[edges_keep])
            edge_weight <- if (isTRUE(config$use_log1p_weight)) log1p(base_weight) else base_weight
            valid_edge <- is.finite(edge_weight) & edge_weight > 0
            edges_keep <- edges_keep[valid_edge, , drop = FALSE]
            base_weight <- base_weight[valid_edge]
            edge_weight <- edge_weight[valid_edge]
            edges_df_final <- data.frame(
                from = assigned_genes[edges_keep[, 1]],
                to = assigned_genes[edges_keep[, 2]],
                weight = edge_weight,
                sign = ifelse(base_weight < 0, "neg", "pos"),
                stringsAsFactors = FALSE
            )
        } else {
            edges_df_final <- data.frame(
                from = character(0), to = character(0),
                weight = numeric(0), sign = character(0),
                stringsAsFactors = FALSE
            )
        }
        g_cons <- igraph::graph_from_data_frame(
            edges_df_final, directed = FALSE, vertices = assigned_genes
        )

        return(list(
            early_return = TRUE,
            result = list(
                membership = stage1_membership_labels,
                genes_all = genes_all,
                kept_genes = assigned_genes,
                consensus_graph = g_cons,
                membership_matrix = stage1$stage1_membership_matrix,
                stage1_consensus = stage1_consensus_matrix,
                metrics_before = NULL,
                metrics_after = NULL,
                genes_before = NULL,
                genes_after = NULL,
                W = igraph::as_adjacency_matrix(g_cons, attr = "weight", sparse = TRUE)
            )
        ))
    }

    list(
        early_return = FALSE,
        stage1_membership_labels = stage1_membership_labels
    )
}

#' Stage2 Configure Runtime Threads
#' @description
#' Internal helper for `.stage2_configure_runtime_threads`.
#' @param config Parameter value.
#' @param stage1_mode_selected Parameter value.
#' @param kept_genes Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage2_configure_runtime_threads <- function(config, stage1_mode_selected, kept_genes) {
    use_consensus <- isTRUE(config$use_consensus)
    consensus_thr <- config$consensus_thr
    n_restart <- max(1L, as.integer(config$n_restart))
    n_threads <- if (is.null(config$n_threads)) 1L else max(1L, as.integer(config$n_threads))
    objective <- config$objective
    res_param_base <- if (identical(objective, "modularity")) as.numeric(config$gamma) else as.numeric(config$resolution)
    if (!is.finite(res_param_base)) res_param_base <- 1
    large_n_threshold <- if (is.null(config$large_n_threshold)) length(kept_genes) else as.integer(config$large_n_threshold)
    aggr_future_workers <- if (is.null(config$aggr_future_workers)) 1L else max(1L, as.integer(config$aggr_future_workers))
    aggr_batch_size <- config$aggr_batch_size
    nk_leiden_iterations <- if (is.null(config$nk_leiden_iterations)) 10L else as.integer(config$nk_leiden_iterations)
    nk_leiden_randomize <- if (is.null(config$nk_leiden_randomize)) TRUE else isTRUE(config$nk_leiden_randomize)
    min_cluster_keep <- if (is.null(config$min_cluster_size)) 2L else max(1L, as.integer(config$min_cluster_size))
    sub_min_size <- max(1L, as.integer(config$sub_min_size))
    sub_min_child_size <- max(1L, as.integer(config$sub_min_child_size))
    sub_resolution_factor <- if (is.null(config$sub_resolution_factor)) 1.3 else config$sub_resolution_factor
    sub_within_cons_max <- config$sub_within_cons_max
    sub_conductance_min <- config$sub_conductance_min
    sub_improve_within_cons_min <- config$sub_improve_within_cons_min
    sub_max_groups <- max(1L, as.integer(config$sub_max_groups))
    enable_qc_filter <- isTRUE(config$enable_qc_filter)
    qc_gene_intra_cons_min <- config$qc_gene_intra_cons_min
    qc_gene_best_out_cons_min <- config$qc_gene_best_out_cons_min
    qc_gene_intra_weight_q <- config$qc_gene_intra_weight_q
    keep_cross_stable <- isTRUE(config$keep_cross_stable)
    min_cutoff <- config$min_cutoff

    list(
        use_consensus = use_consensus,
        consensus_thr = consensus_thr,
        n_restart = n_restart,
        n_threads = n_threads,
        objective = objective,
        res_param_base = res_param_base,
        mode_setting = config$mode,
        large_n_threshold = large_n_threshold,
        aggr_future_workers = aggr_future_workers,
        aggr_batch_size = aggr_batch_size,
        nk_leiden_iterations = nk_leiden_iterations,
        nk_leiden_randomize = nk_leiden_randomize,
        min_cluster_keep = min_cluster_keep,
        sub_min_size = sub_min_size,
        sub_min_child_size = sub_min_child_size,
        sub_resolution_factor = sub_resolution_factor,
        sub_within_cons_max = sub_within_cons_max,
        sub_conductance_min = sub_conductance_min,
        sub_improve_within_cons_min = sub_improve_within_cons_min,
        sub_max_groups = sub_max_groups,
        enable_qc_filter = enable_qc_filter,
        qc_gene_intra_cons_min = qc_gene_intra_cons_min,
        qc_gene_best_out_cons_min = qc_gene_best_out_cons_min,
        qc_gene_intra_weight_q = qc_gene_intra_weight_q,
        keep_cross_stable = keep_cross_stable,
        min_cutoff = min_cutoff,
        stage1_mode_selected = stage1_mode_selected
    )
}

#' Stage2 Run Refine Blocks
#' @description
#' Internal helper for `.stage2_run_refine_blocks`.
#' @param kept_genes Parameter value.
#' @param stage1_membership_labels Parameter value.
#' @param corrected_similarity_graph Parameter value.
#' @param stage2_algo_final Parameter value.
#' @param stage2_backend Parameter value.
#' @param runtime_cfg Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param use_log1p_weight Logical flag.
#' @return Return value used internally.
#' @keywords internal
.stage2_run_refine_blocks <- function(kept_genes,
                                      stage1_membership_labels,
                                      corrected_similarity_graph,
                                      stage2_algo_final,
                                      stage2_backend,
                                      runtime_cfg,
                                      config,
                                      verbose,
                                      use_log1p_weight) {
    validated <- .stage2_refine_blocks_validate_inputs(
        kept_genes = kept_genes,
        stage1_membership_labels = stage1_membership_labels,
        corrected_similarity_graph = corrected_similarity_graph
    )
    spec <- .stage2_refine_blocks_spec_build(
        validated,
        stage2_algo_final = stage2_algo_final,
        stage2_backend = stage2_backend,
        runtime_cfg = runtime_cfg,
        config = config,
        use_log1p_weight = use_log1p_weight
    )
    payload <- .stage2_refine_blocks_materialize(
        spec,
        corrected_similarity_graph = corrected_similarity_graph,
        stage1_membership_labels = stage1_membership_labels,
        config = config,
        verbose = verbose,
        use_log1p_weight = use_log1p_weight
    )
    .stage2_refine_blocks_validate_outputs(payload, spec, validated)
    payload
}

#' Stage2 Apply Qc And Subcluster
#' @description
#' Internal helper for `.stage2_apply_qc_and_subcluster`.
#' @param memb_final Parameter value.
#' @param cons Parameter value.
#' @param metrics_before Parameter value.
#' @param genes_before Parameter value.
#' @param corrected_similarity_graph Parameter value.
#' @param runtime_cfg Parameter value.
#' @param config Parameter value.
#' @param stage2_algo_final Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage2_apply_qc_and_subcluster <- function(memb_final,
                                            cons,
                                            metrics_before,
                                            genes_before,
                                            corrected_similarity_graph,
                                            runtime_cfg,
                                            config,
                                            stage2_algo_final) {
    enable_qc_filter <- runtime_cfg$enable_qc_filter
    qc_gene_intra_weight_q <- runtime_cfg$qc_gene_intra_weight_q
    qc_gene_intra_cons_min <- runtime_cfg$qc_gene_intra_cons_min
    qc_gene_best_out_cons_min <- runtime_cfg$qc_gene_best_out_cons_min
    sub_min_size <- runtime_cfg$sub_min_size
    sub_min_child_size <- runtime_cfg$sub_min_child_size
    sub_resolution_factor <- runtime_cfg$sub_resolution_factor
    sub_within_cons_max <- runtime_cfg$sub_within_cons_max
    sub_conductance_min <- runtime_cfg$sub_conductance_min
    sub_improve_within_cons_min <- runtime_cfg$sub_improve_within_cons_min
    sub_max_groups <- runtime_cfg$sub_max_groups
    consensus_thr <- runtime_cfg$consensus_thr
    n_restart <- runtime_cfg$n_restart
    objective <- runtime_cfg$objective
    res_param_base <- runtime_cfg$res_param_base

    kept_genes <- names(memb_final)
    if (enable_qc_filter) {
        bad <- which(!is.na(genes_before$p_in) &
            genes_before$p_in < qc_gene_intra_cons_min &
            genes_before$p_best_out > qc_gene_best_out_cons_min &
            genes_before$w_in <= quantile(genes_before$w_in, probs = qc_gene_intra_weight_q, na.rm = TRUE))
        if (length(bad)) {
            memb_final[match(genes_before$gene[bad], kept_genes)] <- NA
        }
    }

    if (config$enable_subcluster && !identical(stage2_algo_final, "hotspot-like")) {
        next_id <- if (length(na.omit(memb_final))) max(na.omit(memb_final)) + 1L else 1L
        for (row_idx in seq_len(nrow(metrics_before))) {
            cid <- metrics_before$cluster[row_idx]
            size <- metrics_before$size[row_idx]
            wcons <- metrics_before$within_cons[row_idx]
            condc <- metrics_before$conductance[row_idx]
            if (is.na(size) || size < sub_min_size || is.na(wcons) || is.na(condc)) next
            if (!(wcons < sub_within_cons_max && condc > sub_conductance_min)) next

            genes_c <- kept_genes[memb_final == cid]
            if (!length(genes_c)) next
            g_sub <- igraph::induced_subgraph(corrected_similarity_graph, vids = genes_c)
            if (igraph::ecount(g_sub) == 0) next
            ew <- igraph::E(g_sub)$weight
            el <- igraph::as_data_frame(g_sub, what = "edges")[, c("from", "to")]
            memb_mat_sub <- future.apply::future_sapply(
                seq_len(min(n_restart, 50)),
                function(ii) {
                    g_local <- igraph::graph_from_data_frame(cbind(el, weight = ew), directed = FALSE, vertices = igraph::V(g_sub)$name)
                    res <- .run_graph_algorithm(
                        g_local,
                        stage2_algo_final,
                        if (identical(objective, "modularity")) res_param_base * sub_resolution_factor else res_param_base * sub_resolution_factor,
                        objective
                    )
                    igraph::membership(res)
                },
                future.seed = TRUE
            )
            if (is.null(dim(memb_mat_sub))) memb_mat_sub <- matrix(memb_mat_sub, ncol = 1)
            rownames(memb_mat_sub) <- igraph::V(g_sub)$name
            cons_sub <- .compute_consensus_sparse_r(memb_mat_sub, thr = consensus_thr)
            if (nnzero(cons_sub) == 0) next
            g_cons_sub <- igraph::graph_from_adjacency_matrix(cons_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
            comm_sub <- .run_graph_algorithm(
                g_cons_sub,
                stage2_algo_final,
                if (identical(objective, "modularity")) res_param_base * sub_resolution_factor else res_param_base * sub_resolution_factor,
                objective
            )
            memb_sub <- igraph::membership(comm_sub)

            sub_ids <- sort(unique(memb_sub))
            if (length(sub_ids) <= 1 || length(sub_ids) > sub_max_groups) next
            sizes <- vapply(sub_ids, function(sid) sum(memb_sub == sid), integer(1))
            if (any(sizes < sub_min_child_size)) next

            wcons_sub <- vapply(sub_ids, function(sid) {
                idx_local <- which(memb_sub == sid)
                if (length(idx_local) > 1) mean(cons_sub[idx_local, idx_local][upper.tri(cons_sub[idx_local, idx_local])], na.rm = TRUE) else 0
            }, numeric(1))
            wcons_old <- wcons
            wcons_new <- sum(wcons_sub * sizes) / sum(sizes)
            if (!(is.finite(wcons_new) && wcons_new >= (wcons_old + sub_improve_within_cons_min))) next

            for (sid in sub_ids) {
                genes_sid <- names(memb_sub)[memb_sub == sid]
                memb_final[genes_sid] <- next_id
                next_id <- next_id + 1L
            }
            block_names <- intersect(rownames(cons_sub), genes_c)
            if (length(block_names)) {
                cons[block_names, block_names] <- as.matrix(cons_sub[block_names, block_names])
            }
        }
    }

    list(
        memb_final = memb_final,
        cons = cons
    )
}

#' Stage2 Postprocess And Assemble
#' @description
#' Internal helper for `.stage2_postprocess_and_assemble`.
#' @param memb_final Parameter value.
#' @param cons Parameter value.
#' @param metrics_before Parameter value.
#' @param genes_before Parameter value.
#' @param edges_corr Parameter value.
#' @param corrected_similarity_graph Parameter value.
#' @param W Parameter value.
#' @param kept_genes Parameter value.
#' @param genes_all Parameter value.
#' @param stage1_membership_matrix Parameter value.
#' @param stage1_consensus_matrix Parameter value.
#' @param stage1_similarity_graph Parameter value.
#' @param similarity_matrix Parameter value.
#' @param use_log1p_weight Logical flag.
#' @param runtime_cfg Parameter value.
#' @param config Parameter value.
#' @param L_post Parameter value.
#' @param stage2_backend Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage2_postprocess_and_assemble <- function(memb_final,
                                             cons,
                                             metrics_before,
                                             genes_before,
                                             edges_corr,
                                             corrected_similarity_graph,
                                             W,
                                             kept_genes,
                                             genes_all,
                                             stage1_membership_matrix,
                                             stage1_consensus_matrix,
                                             stage1_similarity_graph,
                                             similarity_matrix,
                                             use_log1p_weight,
                                             runtime_cfg,
                                             config,
                                             L_post,
                                             stage2_backend) {
    min_cluster_keep <- runtime_cfg$min_cluster_keep
    keep_cross_stable <- runtime_cfg$keep_cross_stable
    min_cutoff <- runtime_cfg$min_cutoff

    di <- rowSums((L_post > 0)[kept_genes, kept_genes, drop = FALSE] &
        outer(memb_final, memb_final, "=="), na.rm = TRUE)
    memb_final[di == 0] <- NA
    tbl <- table(memb_final, useNA = "no")
    small <- names(tbl)[tbl < min_cluster_keep]
    if (length(small)) {
        memb_final[memb_final %in% small] <- NA
    }
    tbl <- table(memb_final, useNA = "no")
    if (length(tbl)) {
        unique_ids <- sort(as.integer(names(tbl)))
        id_map <- setNames(seq_along(unique_ids), unique_ids)
        if (length(which(!is.na(memb_final)))) {
            memb_final[which(!is.na(memb_final))] <- as.integer(id_map[as.character(memb_final[which(!is.na(memb_final))])])
        }
        tbl <- table(memb_final, useNA = "no")
        stable_cl <- as.integer(names(tbl)[tbl >= min_cluster_keep])
    } else {
        stable_cl <- integer(0)
    }
    names(memb_final) <- kept_genes

    cons_ts <- as(cons, "TsparseMatrix")
    if (length(cons_ts@x)) {
        sm2 <- data.frame(
            i = cons_ts@i + 1L,
            j = cons_ts@j + 1L,
            x = cons_ts@x,
            stringsAsFactors = FALSE
        )
    } else {
        sm2 <- data.frame(i = integer(0), j = integer(0), x = numeric(0), stringsAsFactors = FALSE)
    }
    sm2 <- sm2[sm2$i < sm2$j & sm2$x > 0, , drop = FALSE]

    key_corr <- paste(pmin(edges_corr$from, edges_corr$to), pmax(edges_corr$from, edges_corr$to), sep = "|")
    intra_keys <- paste(
        pmin(kept_genes[sm2$i], kept_genes[sm2$j]),
        pmax(kept_genes[sm2$i], kept_genes[sm2$j]),
        sep = "|"
    )
    keep_mask <- (function(vec_tmp) !is.na(vec_tmp) & vec_tmp > 0)(
        as.numeric(setNames(edges_corr$weight, key_corr)[intra_keys])
    )
    intra_df <- data.frame(
        from = kept_genes[sm2$i[keep_mask]],
        to = kept_genes[sm2$j[keep_mask]],
        weight = as.numeric(setNames(edges_corr$weight, key_corr)[intra_keys])[keep_mask],
        sign = ifelse(similarity_matrix[cbind(
            match(kept_genes[sm2$i[keep_mask]], rownames(similarity_matrix)),
            match(kept_genes[sm2$j[keep_mask]], colnames(similarity_matrix))
        )] < 0, "neg", "pos"),
        stringsAsFactors = FALSE
    )

    cross_df <- NULL
    if (keep_cross_stable && length(stable_cl) > 1) {
        genes_stable <- kept_genes[memb_final %in% stable_cl]
        if (length(genes_stable) > 1) {
	            L_sub_raw <- L_post[genes_stable, genes_stable, drop = FALSE]
	            ij <- .arrind_from_matrix_predicate(
	                L_sub_raw,
	                op = "ge",
	                cutoff = min_cutoff,
	                triangle = "upper",
	                keep_diag = FALSE
	            )
            if (nrow(ij)) {
                g1 <- genes_stable[ij[, 1]]
                g2 <- genes_stable[ij[, 2]]
                cl1 <- memb_final[g1]
                cl2 <- memb_final[g2]
                keep_cross <- cl1 != cl2 & !is.na(cl1) & !is.na(cl2)
                if (any(keep_cross)) {
                    g1 <- g1[keep_cross]
                    g2 <- g2[keep_cross]
                    Lval <- L_sub_raw[cbind(match(g1, genes_stable), match(g2, genes_stable))]
                    mask <- (if (use_log1p_weight) log1p(Lval) else Lval) > 0
                    if (any(mask)) {
                        cross_df <- data.frame(
                            from = g1[mask],
                            to = g2[mask],
                            weight = (if (use_log1p_weight) log1p(Lval) else Lval)[mask],
                            sign = ifelse(Lval[mask] < 0, "neg", "pos")
                        )
                    }
                }
            }
        }
    }

    edges_df_final <- if (!is.null(cross_df)) unique(rbind(intra_df, cross_df)) else intra_df
    if (!nrow(edges_df_final)) {
        cons_graph_final <- igraph::make_empty_graph(n = length(kept_genes), directed = FALSE)
        cons_graph_final <- igraph::set_vertex_attr(cons_graph_final, "name", value = kept_genes)
    } else {
        cons_graph_final <- igraph::graph_from_data_frame(edges_df_final,
            directed = FALSE, vertices = kept_genes
        )
    }

    list(
        membership = memb_final,
        kept_genes = kept_genes,
        genes_all = genes_all,
        consensus_graph = cons_graph_final,
        membership_matrix = stage1_membership_matrix,
        stage1_consensus = stage1_consensus_matrix,
        metrics_before = metrics_before,
        metrics_after = .summarise_cluster_metrics(memb_final, cons, W),
        genes_before = genes_before,
        genes_after = .summarise_gene_membership(memb_final, cons, W, kept_genes),
        corrected_graph = corrected_similarity_graph,
        weights_matrix = W,
        stage2_backend = stage2_backend
    )
}

# Compatibility shims for legacy stage workflow entrypoints.
if (!exists(".stage1_consensus_workflow", inherits = FALSE) && exists(".stage1_consensus_workflow_v2", inherits = FALSE)) {
    .stage1_consensus_workflow <- .stage1_consensus_workflow_v2
}
#' Summarise Cluster Metrics
#' @description
#' Internal helper for `.summarise_cluster_metrics`.
#' @param memberships Parameter value.
#' @param consensus_matrix Parameter value.
#' @param weight_matrix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.summarise_cluster_metrics <- function(memberships,
                                    consensus_matrix,
                                    weight_matrix) {
    if (!length(memberships)) {
        return(data.frame(
            cluster = integer(0),
            size = integer(0),
            within_cons = numeric(0),
            conductance = numeric(0)
        ))
    }
    cl_ids <- sort(na.omit(unique(memberships)))
    if (!length(cl_ids)) {
        return(data.frame(
            cluster = integer(0),
            size = integer(0),
            within_cons = numeric(0),
            conductance = numeric(0)
        ))
    }
    res <- lapply(cl_ids, function(cid) {
        idx <- which(memberships == cid)
        n <- length(idx)
        if (n <= 1) {
            return(data.frame(
                cluster = cid,
                size = n,
                within_cons = NA_real_,
                conductance = NA_real_
            ))
        }
        cons_block <- as.matrix(consensus_matrix[idx, idx])
        within_cons <- mean(cons_block[upper.tri(cons_block)], na.rm = TRUE)
        W_block <- as.matrix(weight_matrix[idx, idx])
        Win <- sum(W_block) / 2
        Wcut <- sum(weight_matrix[idx, -idx, drop = FALSE])
        conductance <- Wcut / (Wcut + Win + 1e-12)
        data.frame(
            cluster = cid,
            size = n,
            within_cons = within_cons,
            conductance = conductance
        )
    })
    do.call(rbind, res)
}

#' Summarise Gene Membership
#' @description
#' Internal helper for `.summarise_gene_membership`.
#' @param memberships Parameter value.
#' @param consensus_matrix Parameter value.
#' @param weight_matrix Parameter value.
#' @param genes Parameter value.
#' @return Return value used internally.
#' @keywords internal
.summarise_gene_membership <- function(memberships,
                                         consensus_matrix,
                                         weight_matrix,
                                         genes = NULL) {
    n <- length(memberships)
    if (is.null(genes)) {
        genes <- rownames(consensus_matrix)
        if (is.null(genes)) genes <- as.character(seq_len(n))
    }
    if (length(genes) != n) {
        stop("Length of `genes` must match length of `memberships`.", call. = FALSE)
    }
    cl_ids <- sort(na.omit(unique(memberships)))
    p_in <- numeric(n); p_best_out <- numeric(n); w_in <- numeric(n)
    for (i in seq_len(n)) {
        cid <- memberships[i]
        if (is.na(cid)) {
            p_in[i] <- NA_real_
            p_best_out[i] <- NA_real_
            w_in[i] <- NA_real_
            next
        }
        idx <- which(memberships == cid)
        idx_noi <- setdiff(idx, i)
        if (length(idx_noi)) {
            p_in[i] <- mean(consensus_matrix[i, idx_noi], na.rm = TRUE)
            w_in[i] <- sum(weight_matrix[i, idx_noi])
        } else {
            p_in[i] <- 0
            w_in[i] <- 0
        }
        other <- setdiff(cl_ids, cid)
        if (length(other)) {
            means <- vapply(other, function(cj) {
                jj <- which(memberships == cj)
                if (length(jj)) mean(consensus_matrix[i, jj], na.rm = TRUE) else 0
            }, numeric(1))
            p_best_out[i] <- max(means, na.rm = TRUE)
        } else {
            p_best_out[i] <- 0
        }
    }
    data.frame(
        gene = genes,
        cluster = memberships,
        p_in = p_in,
        p_best_out = p_best_out,
        w_in = w_in,
        stringsAsFactors = FALSE
    )
}

#' Use Conda Environment
#' @description
#' Internal helper for `.use_conda_environment`.
#' @param env_name Parameter value.
#' @param conda_bin Parameter value.
#' @param required Parameter value.
#' @return Return value used internally.
#' @keywords internal
.use_conda_environment <- function(env_name, conda_bin = NULL, required = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))
    if (is.null(conda_bin)) {
        reticulate::use_condaenv(env_name, required = required)
    } else {
        reticulate::use_condaenv(env_name, conda = conda_bin, required = required)
    }
    assign("configured", TRUE, envir = .networkit_runtime_state)
    invisible(TRUE)
}

#' Use Python Interpreter
#' @description
#' Internal helper for `.use_python_interpreter`.
#' @param python_path Filesystem path.
#' @param required Parameter value.
#' @return Return value used internally.
#' @keywords internal
.use_python_interpreter <- function(python_path, required = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))
    reticulate::use_python(python_path, required = required)
    assign("configured", TRUE, envir = .networkit_runtime_state)
    invisible(TRUE)
}

#' Use Virtual Environment
#' @description
#' Internal helper for `.use_virtual_environment`.
#' @param env_name Parameter value.
#' @param required Parameter value.
#' @return Return value used internally.
#' @keywords internal
.use_virtual_environment <- function(env_name, required = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))
    reticulate::use_virtualenv(env_name, required = required)
    assign("configured", TRUE, envir = .networkit_runtime_state)
    invisible(TRUE)
}

`%||%` <- function(x, y) {
    if (!is.null(x)) x else y
}

#' Stage-1 consensus workflow (baseline)
#' @description
#' Executes multiple community-detection runs on the similarity graph, aggregates
#' them via consensus (or single run when consensus disabled), and returns the
#' stage-1 memberships plus consensus graph artifacts.
#' @param similarity_matrix Full similarity matrix (dgCMatrix).
#' @param threshold List from thresholding stage containing `kept_genes`,
#'   `similarity_sub`, and `genes_all`.
#' @param config List of clustering parameters (algo/objective/mode, consensus,
#'   restart counts, NetworKit options, etc.).
#' @param verbose Logical; emit progress logs.
#' @param networkit_config Optional reticulate configuration for NetworKit backend.
#' @return List with stage-1 membership labels/matrix, consensus matrix/graph,
#'   and bookkeeping fields used by later stages.
#' @keywords internal
# Baseline overrides for Stage-1/Stage-2 consensus workflows to preserve algorithm equivalence
#' Stage1 Consensus Workflow V2
#' @description
#' Internal helper for `.stage1_consensus_workflow_v2`.
#' @param similarity_matrix Parameter value.
#' @param threshold Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param networkit_config Parameter value.
#' @return Return value used internally.
#' @keywords internal
.stage1_consensus_workflow_v2 <- function(similarity_matrix,
                                      threshold,
                                      config,
                                      verbose = TRUE,
                                      networkit_config = list()) {
    .cluster_message(verbose, "[cluster] Step2: running Stage-1 consensus…",
        parent = "clusterGenes", step = "S05"
    )

    kept_genes <- threshold$kept_genes
    similarity_sub <- threshold$similarity_sub
    genes_all <- threshold$genes_all

    if (!is.null(networkit_config) && length(networkit_config)) {
        config_conda_env <- networkit_config[["conda_env"]]
        if (is.null(config_conda_env)) config_conda_env <- networkit_config[["nk_condaenv"]]
        config_conda_bin <- networkit_config[["conda_bin"]]
        if (is.null(config_conda_bin)) config_conda_bin <- networkit_config[["nk_conda_bin"]]
        config_python_path <- networkit_config[["python_path"]]
        if (is.null(config_python_path)) config_python_path <- networkit_config[["nk_python"]]
        try(.configure_networkit_python(
            conda_env = config_conda_env,
            conda_bin = config_conda_bin,
            python_path = config_python_path,
            required = FALSE
        ), silent = TRUE)
    }

    adjacency_triplet <- as.data.frame(summary(similarity_sub))
    if (!all(c("i", "j", "x") %in% colnames(adjacency_triplet))) {
        colnames(adjacency_triplet) <- c("i", "j", "x")[seq_len(ncol(adjacency_triplet))]
    }
    adjacency_triplet <- adjacency_triplet[adjacency_triplet$i < adjacency_triplet$j & adjacency_triplet$x > 0, , drop = FALSE]
    if (!nrow(adjacency_triplet)) stop("All edges filtered out at Stage-1.")

    edge_table <- data.frame(
        from = kept_genes[adjacency_triplet$i],
        to = kept_genes[adjacency_triplet$j],
        sim_edge = adjacency_triplet$x,
        stringsAsFactors = FALSE
    )
    edge_table$weight <- if (config$use_log1p_weight) log1p(edge_table$sim_edge) else edge_table$sim_edge
    similarity_graph <- igraph::graph_from_data_frame(edge_table[, c("from", "to", "weight")],
        directed = FALSE, vertices = kept_genes
    )

    nk_inputs_stage1 <- .nk_prepare_graph_input(similarity_graph)
    similarity_graph <- nk_inputs_stage1$graph

    algo <- match.arg(config$algo, c("leiden", "louvain", "hotspot-like"))
    objective <- match.arg(config$objective, c("CPM", "modularity"))
    prefer_fast <- isTRUE(config$prefer_fast)
    mode_selected <- switch(config$mode,
        fast = "fast",
        auto = if (prefer_fast) "fast" else if (length(kept_genes) > config$large_n_threshold) "aggressive" else "safe",
        safe = "safe",
        aggressive = "aggressive"
    )
    backend_mode <- if (identical(mode_selected, "fast")) "safe" else mode_selected
    .cluster_message(verbose, sprintf("[cluster]   Stage-1 mode=%s (nodes=%d)", mode_selected, length(kept_genes)),
        parent = "clusterGenes", step = "S05"
    )
    if (identical(backend_mode, "aggressive") && identical(objective, "CPM")) {
        warning("[cluster] CPM objective not supported in aggressive mode; falling back to 'modularity'.")
        objective <- "modularity"
    }

    stage1_consensus_matrix <- NULL
    stage1_membership_labels <- NULL
    restart_memberships <- NULL
    stage1_backend <- "igraph"
    stage1_backend <- "igraph"

    if (identical(algo, "hotspot-like")) {
        .log_backend("clusterGenes", "S05", "stage1_backend", "hotspot-like",
            reason = "algo=hotspot-like", verbose = verbose
        )
        hotspot_k <- config$hotspot_k
        if (is.null(hotspot_k) || !is.finite(hotspot_k)) {
            hotspot_k <- max(1L, min(20L, length(kept_genes)))
        }
        hotspot_result <- .run_hotspot_clustering(similarity_sub, hotspot_k, config$hotspot_min_module_size, config$use_log1p_weight)
        stage1_membership_labels <- hotspot_result$membership
        stage1_consensus_matrix <- hotspot_result$consensus
        restart_memberships <- matrix(as.integer(stage1_membership_labels), ncol = 1,
            dimnames = list(kept_genes, "hotspot"))
        stage1_backend <- "hotspot"
    } else {
        n_threads <- .clamp_ncores_safe(config$n_threads %||% 1L, logical = TRUE)

        edge_weights <- nk_inputs_stage1$edge_weights
        vertex_names <- nk_inputs_stage1$vertex_names
        edge_pairs <- nk_inputs_stage1$edge_table

        res_param <- if (identical(objective, "modularity")) as.numeric(config$gamma) else as.numeric(config$resolution)
        algo_per_run <- if (identical(backend_mode, "aggressive") && identical(algo, "leiden")) "louvain" else algo

        future_guard <- .adjust_future_globals_maxsize(config$future_globals_min_bytes %||% (2 * 1024^3))
        on.exit(future_guard(), add = TRUE)

        aggressive_requested <- identical(backend_mode, "aggressive")
        if (aggressive_requested && !.networkit_available()) {
            stop("Aggressive mode requires the NetworKit backend; please configure reticulate/networkit before running clusterGenes_new.", call. = FALSE)
        }
        use_networkit_leiden <- aggressive_requested && identical(algo_per_run, "leiden")
        use_networkit_plm <- aggressive_requested && identical(algo_per_run, "louvain")

        if (use_networkit_leiden || use_networkit_plm) {
            backend_label <- if (use_networkit_leiden) "Leiden" else "PLM/louvain"
            .log_backend("clusterGenes", "S05", "stage1_backend", "networkit",
                reason = paste0("algo=", backend_label), verbose = verbose
            )
            stage1_backend <- "networkit"
            seeds <- seq_len(config$n_restart)
            chunk_size <- config$aggr_batch_size
            if (is.null(chunk_size) || chunk_size <= 0) {
                chunk_size <- ceiling(config$n_restart / max(1L, as.integer(config$aggr_future_workers)))
            }
            split_chunks <- function(x, k) {
                if (k <= 0) return(list(x))
                split(x, ceiling(seq_along(x) / k))
            }
            chunks <- split_chunks(seeds, chunk_size)
            threads_per_worker <- max(1L, floor(n_threads / max(1L, as.integer(config$aggr_future_workers))))

            worker_fun <- function(ss) {
                if (use_networkit_leiden) {
                    .networkit_parallel_leiden_runs(
                        vertex_names = vertex_names,
                        edge_table = edge_pairs,
                        edge_weights = edge_weights,
                        gamma = as.numeric(res_param),
                        iterations = as.integer(config$nk_leiden_iterations),
                        randomize = isTRUE(config$nk_leiden_randomize),
                        threads = threads_per_worker,
                        seeds = as.integer(ss)
                    )
                } else {
                    .networkit_plm_runs(
                        vertex_names = vertex_names,
                        edge_table = edge_pairs,
                        edge_weights = edge_weights,
                        gamma = as.numeric(res_param),
                        refine = TRUE,
                        threads = threads_per_worker,
                        seeds = as.integer(ss)
                    )
                }
            }
            membership_blocks <- if (config$aggr_future_workers > 1L) {
                future.apply::future_lapply(chunks, worker_fun, future.seed = TRUE)
            } else {
                lapply(chunks, worker_fun)
            }
            restart_memberships <- do.call(cbind, membership_blocks)
        } else {
            .log_backend("clusterGenes", "S05", "stage1_backend", "igraph",
                reason = paste0("algo=", algo_per_run), verbose = verbose
            )
            restart_memberships <- future.apply::future_sapply(
                seq_len(config$n_restart),
                function(ii) {
                    g_local <- igraph::graph_from_data_frame(
                        cbind(edge_pairs, weight = edge_weights),
                        directed = FALSE,
                        vertices = vertex_names
                    )
                    res <- try(.run_graph_algorithm(g_local, algo_per_run, res_param, objective), silent = TRUE)
                    if (inherits(res, "try-error")) {
                        stop("Stage-1 graph algorithm failed: ", as.character(res), call. = FALSE)
                    } else {
                        mm <- try(igraph::membership(res), silent = TRUE)
                        if (inherits(mm, "try-error") || is.null(mm)) {
                            stop("Stage-1 membership extraction failed: ", as.character(mm), call. = FALSE)
                        }
                        mm
                    }
                },
                future.seed = TRUE
            )
        }
        if (is.null(dim(restart_memberships))) restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
        if (!is.matrix(restart_memberships) && !inherits(restart_memberships, "Matrix")) {
            restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
        }
        if (is.null(nrow(restart_memberships)) || nrow(restart_memberships) != length(kept_genes)) {
            len <- length(restart_memberships)
            if (len %% length(kept_genes) == 0) {
                restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
            } else {
                stop("Stage-1 membership shape mismatch: ", len, " elements vs ", length(kept_genes))
            }
        }
        rownames(restart_memberships) <- kept_genes

        if (!config$use_consensus) {
            stage1_membership_labels <- restart_memberships[, 1]
            stage1_consensus_matrix <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                dims = c(length(kept_genes), length(kept_genes)),
                dimnames = list(kept_genes, kept_genes))
        } else {
            stage1_consensus_matrix <- try(consensus_sparse(restart_memberships,
                thr = config$consensus_thr,
                n_threads = n_threads
            ), silent = TRUE)
            if (!inherits(stage1_consensus_matrix, "try-error")) {
                .log_backend("clusterGenes", "S05", "stage1_consensus_backend", "cpp",
                    reason = "consensus_sparse", verbose = verbose
                )
            } else {
                .log_backend("clusterGenes", "S05", "stage1_consensus_backend", "r",
                    reason = "consensus_sparse_failed", verbose = verbose
                )
                stage1_consensus_matrix <- .compute_consensus_sparse_r(restart_memberships, thr = config$consensus_thr)
            }
            try({ if (length(stage1_consensus_matrix@i)) diag(stage1_consensus_matrix) <- 0 }, silent = TRUE)
            if (nnzero(stage1_consensus_matrix) == 0) {
                stage1_membership_labels <- restart_memberships[, 1]
            } else {
                consensus_graph_view <- igraph::graph_from_adjacency_matrix(stage1_consensus_matrix,
                    mode = "undirected", weighted = TRUE, diag = FALSE)
                consensus_backend <- if (identical(stage1_backend, "networkit")) "networkit" else "igraph"
                backend_label <- if (identical(consensus_backend, "networkit")) {
                    paste0("NetworKit (", algo, ")")
                } else {
                    paste0("igraph (", algo, ")")
                }
                .log_backend("clusterGenes", "S05", "stage1_consensus_graph_backend", backend_label,
                    reason = paste0("stage1_backend=", stage1_backend), verbose = verbose
                )
                nk_opts <- list(
                    iterations = config$nk_leiden_iterations,
                    randomize = config$nk_leiden_randomize,
                    refine = TRUE,
                    seeds = 1L
                )
                stage1_membership_labels <- .run_single_partition(
                    consensus_graph_view,
                    backend = consensus_backend,
                    algo_name = algo,
                    res_param = res_param,
                    objective = objective,
                    threads = n_threads,
                    nk_opts = nk_opts
                )
            }
        }
    }

    stage1_edge_triplet <- summary(stage1_consensus_matrix)
    if (!all(c("i", "j", "x") %in% colnames(stage1_edge_triplet))) {
        colnames(stage1_edge_triplet) <- c("i", "j", "x")[seq_len(ncol(stage1_edge_triplet))]
    }
    stage1_edge_triplet <- stage1_edge_triplet[stage1_edge_triplet$i < stage1_edge_triplet$j & stage1_edge_triplet$x > 0, , drop = FALSE]
    stage1_cons_graph <- NULL
    if (nrow(stage1_edge_triplet)) {
        consensus_edge_df <- data.frame(
            from = kept_genes[stage1_edge_triplet$i],
            to = kept_genes[stage1_edge_triplet$j],
            weight = stage1_edge_triplet$x,
            stringsAsFactors = FALSE
        )
        stage1_cons_graph <- igraph::graph_from_data_frame(
            consensus_edge_df,
            directed = FALSE,
            vertices = kept_genes
        )
    } else {
        stage1_cons_graph <- igraph::make_empty_graph(n = length(kept_genes), directed = FALSE)
        stage1_cons_graph <- igraph::set_vertex_attr(stage1_cons_graph, "name", value = kept_genes)
    }

    list(
        kept_genes = kept_genes,
        genes_all = genes_all,
        similarity_full = threshold$similarity,
        similarity_sub = similarity_sub,
        A = threshold$similarity,
        A_sub = similarity_sub,
        g_stage1 = similarity_graph,
        mode_selected = mode_selected,
        stage1_membership_matrix = restart_memberships,
        stage1_consensus = stage1_consensus_matrix,
        stage1_membership = stage1_membership_labels,
        stage1_consensus_graph = stage1_cons_graph,
        stage1_backend = stage1_backend,
        stage1_algo = algo,
        stage1_algo_per_run = if (exists("algo_per_run")) algo_per_run else algo
    )
}

#' Build non-destructive consensus graph views
#' @description
#' Adds provenance to the legacy Stage-1 filtered L backbone without changing
#' its vertices, edge order, or `weight` values.  It also creates two explicit
#' derived views: the pure co-clustering-frequency graph and the intersection
#' of that graph with the legacy L backbone.
#' @param backbone_graph Legacy graph returned by Stage 2.
#' @param stage1_consensus_matrix Sparse Stage-1 co-clustering matrix.  This is
#'   already thresholded by `consensus_thr`.
#' @param stage1_membership_matrix Gene-by-restart membership matrix, used to
#'   recover the exact co-clustering frequency for every backbone edge,
#'   including edges below `consensus_thr`.
#' @param raw_similarity Untransformed similarity matrix used to annotate
#'   backbone edges with `raw_L`.
#' @param consensus_thr Co-clustering-frequency threshold.
#' @param use_log1p_weight Whether legacy backbone weights were log1p
#'   transformed.
#' @return A list with `backbone_graph`, `consensus_frequency_graph`, and
#'   `paper_consensus_graph`.
#' @keywords internal
.build_consensus_graph_views <- function(backbone_graph,
                                         stage1_consensus_matrix,
                                         stage1_membership_matrix,
                                         raw_similarity,
                                         consensus_thr,
                                         use_log1p_weight) {
    stopifnot(igraph::is_igraph(backbone_graph))

    threshold <- as.numeric(consensus_thr)[1]
    if (!is.finite(threshold) || threshold < 0 || threshold > 1) {
        stop("consensus_thr must be a finite value in [0, 1].", call. = FALSE)
    }
    vertices <- igraph::V(backbone_graph)$name
    if (is.null(vertices)) vertices <- as.character(seq_len(igraph::vcount(backbone_graph)))
    transform_name <- if (isTRUE(use_log1p_weight)) "log1p" else "identity"

    edge_table <- igraph::as_data_frame(backbone_graph, what = "edges")
    if (nrow(edge_table)) {
        raw_i <- match(edge_table$from, rownames(raw_similarity))
        raw_j <- match(edge_table$to, colnames(raw_similarity))
        if (anyNA(raw_i) || anyNA(raw_j)) {
            stop("Backbone graph edges cannot be aligned to raw_similarity.", call. = FALSE)
        }
        raw_L <- as.numeric(raw_similarity[cbind(raw_i, raw_j)])

        frequency <- rep(NA_real_, nrow(edge_table))
        membership <- stage1_membership_matrix
        if (!is.null(membership)) {
            if (is.null(dim(membership))) membership <- matrix(membership, ncol = 1L)
            membership <- as.matrix(membership)
            membership_names <- rownames(membership)
            if (is.null(membership_names) && nrow(membership) == nrow(raw_similarity)) {
                membership_names <- rownames(raw_similarity)
            }
            if (!is.null(membership_names)) {
                mi <- match(edge_table$from, membership_names)
                mj <- match(edge_table$to, membership_names)
                aligned <- !is.na(mi) & !is.na(mj)
                if (any(aligned)) {
                    frequency[aligned] <- vapply(which(aligned), function(k) {
                        lhs <- membership[mi[k], ]
                        rhs <- membership[mj[k], ]
                        mean(!is.na(lhs) & !is.na(rhs) & lhs == rhs)
                    }, numeric(1))
                }
            }
        }

        # The thresholded sparse matrix is a safe fallback for old objects
        # that did not retain restart memberships.  Missing entries remain NA
        # because their exact below-threshold frequency is not recoverable.
        unresolved <- which(!is.finite(frequency))
        if (length(unresolved) && !is.null(stage1_consensus_matrix)) {
            cons_names <- rownames(stage1_consensus_matrix)
            if (!is.null(cons_names)) {
                ci <- match(edge_table$from[unresolved], cons_names)
                cj <- match(edge_table$to[unresolved], cons_names)
                ok <- !is.na(ci) & !is.na(cj)
                if (any(ok)) {
                    fallback_frequency <- as.numeric(
                        stage1_consensus_matrix[cbind(ci[ok], cj[ok])]
                    )
                    known <- is.finite(fallback_frequency) & fallback_frequency > 0
                    if (any(known)) {
                        frequency[unresolved[ok][known]] <- fallback_frequency[known]
                    }
                }
            }
        }

        backbone_graph <- igraph::set_edge_attr(backbone_graph, "raw_L", value = raw_L)
        backbone_graph <- igraph::set_edge_attr(
            backbone_graph, "weight_transform", value = rep(transform_name, nrow(edge_table))
        )
        backbone_graph <- igraph::set_edge_attr(
            backbone_graph, "consensus_frequency", value = frequency
        )
        backbone_graph <- igraph::set_edge_attr(
            backbone_graph, "is_consensus_edge",
            value = is.finite(frequency) & frequency >= threshold
        )
    } else {
        backbone_graph <- igraph::set_edge_attr(backbone_graph, "raw_L", value = numeric(0))
        backbone_graph <- igraph::set_edge_attr(backbone_graph, "weight_transform", value = character(0))
        backbone_graph <- igraph::set_edge_attr(backbone_graph, "consensus_frequency", value = numeric(0))
        backbone_graph <- igraph::set_edge_attr(backbone_graph, "is_consensus_edge", value = logical(0))
    }
    backbone_graph <- igraph::set_graph_attr(backbone_graph, "graph_role", "thresholded_L_backbone")
    backbone_graph <- igraph::set_graph_attr(backbone_graph, "consensus_threshold", threshold)
    backbone_graph <- igraph::set_graph_attr(backbone_graph, "weight_transform", transform_name)
    backbone_graph <- igraph::set_graph_attr(backbone_graph, "weight_semantics",
        if (isTRUE(use_log1p_weight)) "log1p(raw_L)" else "raw_L")
    backbone_graph <- igraph::set_graph_attr(backbone_graph, "graph_schema_version", "1")

    consensus_sub <- NULL
    if (!is.null(stage1_consensus_matrix)) {
        cons_names <- rownames(stage1_consensus_matrix)
        if (!is.null(cons_names) && all(vertices %in% cons_names)) {
            consensus_sub <- stage1_consensus_matrix[vertices, vertices, drop = FALSE]
        }
    }
    if (is.null(consensus_sub)) {
        consensus_sub <- Matrix::Matrix(0, length(vertices), length(vertices), sparse = TRUE,
            dimnames = list(vertices, vertices))
    }
    consensus_sub <- Matrix::drop0(Matrix::Matrix(consensus_sub, sparse = TRUE))
    if (nrow(consensus_sub)) diag(consensus_sub) <- 0
    cons_triplet <- as.data.frame(summary(consensus_sub))
    if (nrow(cons_triplet)) {
        cons_triplet <- cons_triplet[
            cons_triplet$i < cons_triplet$j &
                is.finite(cons_triplet$x) & cons_triplet$x >= threshold,
            , drop = FALSE
        ]
    }
    if (nrow(cons_triplet)) {
        consensus_edges <- data.frame(
            from = vertices[cons_triplet$i],
            to = vertices[cons_triplet$j],
            weight = cons_triplet$x,
            consensus_frequency = cons_triplet$x,
            is_consensus_edge = TRUE,
            stringsAsFactors = FALSE
        )
        consensus_graph <- igraph::graph_from_data_frame(
            consensus_edges, directed = FALSE, vertices = vertices
        )
    } else {
        consensus_graph <- igraph::make_empty_graph(n = length(vertices), directed = FALSE)
        consensus_graph <- igraph::set_vertex_attr(consensus_graph, "name", value = vertices)
    }
    consensus_graph <- igraph::set_graph_attr(consensus_graph, "graph_role", "consensus_frequency_graph")
    consensus_graph <- igraph::set_graph_attr(consensus_graph, "consensus_threshold", threshold)
    consensus_graph <- igraph::set_graph_attr(consensus_graph, "weight_semantics", "co_clustering_frequency")
    consensus_graph <- igraph::set_graph_attr(consensus_graph, "graph_schema_version", "1")

    keep_ids <- if (igraph::ecount(backbone_graph)) {
        which(igraph::E(backbone_graph)$is_consensus_edge %in% TRUE)
    } else {
        integer(0)
    }
    paper_graph <- igraph::subgraph.edges(backbone_graph, eids = keep_ids, delete.vertices = FALSE)
    paper_graph <- igraph::set_graph_attr(paper_graph, "graph_role", "paper_consensus_network")
    paper_graph <- igraph::set_graph_attr(paper_graph, "consensus_threshold", threshold)
    paper_graph <- igraph::set_graph_attr(
        paper_graph, "edge_definition", "FDR/quantile L backbone intersect co-clustering frequency threshold"
    )
    paper_graph <- igraph::set_graph_attr(paper_graph, "graph_schema_version", "1")

    list(
        backbone_graph = backbone_graph,
        consensus_frequency_graph = consensus_graph,
        paper_consensus_graph = paper_graph
    )
}

#' Stage-2 refinement workflow (baseline)
#' @description
#' Internal helper for `.stage2_refine_workflow_v2`.
#' Applies MH weighting/CI95 filtering and optional sub-clustering to refine
#' stage-1 results, producing corrected memberships and non-destructive graph
#' views. The legacy graph remains the filtered L backbone; consensus-frequency
#' views are returned separately.
#' @param stage1 Output list from `.stage1_consensus_workflow_v2()`.
#' @param similarity_matrix Full similarity matrix.
#' @param config List of refinement parameters (algo selections, thresholds,
#'   consensus options, minimum cluster sizes, etc.).
#' @param FDR Optional FDR matrix aligned to similarity_matrix when CI95 filtering.
#' @param mh_object Optional MH weighting object (igraph or matrix).
#' @param aux_stats Optional auxiliary stats (e.g., curve layers for CI95).
#' @param pearson_matrix Optional Pearson correlation matrix for CI95 filtering.
#' @param verbose Logical; emit progress logs.
#' @return List containing refined membership, the backward-compatible L
#'   backbone, pure consensus-frequency graph, paper-consensus intersection,
#'   and associated metrics/intermediates.
#' @keywords internal
.stage2_refine_workflow_v2 <- function(stage1,
                                   similarity_matrix,
                                   config,
                                   FDR = NULL,
                                   mh_object = NULL,
                                   aux_stats = NULL,
                                   pearson_matrix = NULL,
                                   verbose = TRUE) {
    .cluster_message(verbose, "[cluster] Step3: Stage-2 correction & clustering…",
        parent = "clusterGenes", step = "S07"
    )

    kept_genes <- stage1$kept_genes
    genes_all <- stage1$genes_all
    stage1_membership_labels <- stage1$stage1_membership
    stage1_consensus_matrix <- stage1$stage1_consensus
    stage1_similarity_graph <- stage1$g_stage1
    stage1_mode_selected <- stage1$mode_selected
    if (is.null(stage1_mode_selected)) stage1_mode_selected <- "safe"
    stage1_membership_matrix <- stage1$stage1_membership_matrix
    stage1_backend <- stage1$stage1_backend %||% "igraph"
    stage1_algo_effective <- stage1$stage1_algo_per_run %||% stage1$stage1_algo %||% config$algo

    stage2_algo_final <- tryCatch({
        if (is.null(config$stage2_algo)) config$algo else match.arg(config$stage2_algo, c("leiden", "louvain", "hotspot-like"))
    }, error = function(e) config$algo)
    .cluster_message(verbose, "[cluster]   Stage-2 base algorithm: ", stage2_algo_final,
        parent = "clusterGenes", step = "S07"
    )
    stage2_backend <- if (!is.null(config$stage2_algo)) {
        if (identical(stage2_algo_final, "hotspot-like")) "hotspot" else "igraph"
    } else {
        if (identical(stage2_algo_final, "hotspot-like")) {
            "hotspot"
        } else if (identical(stage1_backend, "networkit")) {
            "networkit"
        } else {
            "igraph"
        }
    }
    if (isTRUE(config$prefer_fast) && !identical(stage2_backend, "hotspot")) {
        stage2_backend <- "igraph"
    }
    reason <- if (!is.null(config$stage2_algo)) {
        paste0("stage2_algo=", stage2_algo_final)
    } else {
        paste0("stage1_backend=", stage1_backend)
    }
    .log_backend("clusterGenes", "S07", "stage2_backend", stage2_backend,
        reason = reason, verbose = verbose
    )

    if (!config$use_mh_weight && !config$CI95_filter) {
        .log_info("clusterGenes", "S07", "no Stage-2 corrections requested; returning Stage-1 result", verbose)

        if (length(stage1_membership_labels) && length(kept_genes)) {
            stage1_mem <- stage1_membership_labels
            if (!is.null(names(stage1_mem))) {
                stage1_mem <- stage1_mem[kept_genes]
            } else if (length(stage1_mem) == length(kept_genes)) {
                names(stage1_mem) <- kept_genes
            } else {
                stage1_mem <- setNames(rep(NA_integer_, length(kept_genes)), kept_genes)
            }
            stage1_mem <- as.integer(stage1_mem)
            names(stage1_mem) <- kept_genes
            sim_stage1 <- stage1$similarity_sub %||% stage1$A_sub
            if (!is.null(sim_stage1)) {
                if (!identical(rownames(sim_stage1), kept_genes)) {
                    sim_stage1 <- sim_stage1[kept_genes, kept_genes, drop = FALSE]
                }
                min_cluster_keep <- max(1L, as.integer(config$min_cluster_size %||% 2L))
                pruned <- .prune_unstable_memberships(
                    memberships = stage1_mem,
                    L_matrix = sim_stage1,
                    min_cluster_size = min_cluster_keep
                )
                stage1_mem <- pruned$membership
                idx_valid <- which(!is.na(stage1_mem))
                if (length(idx_valid)) {
                    unique_ids <- sort(unique(stage1_mem[idx_valid]))
                    id_map <- setNames(seq_along(unique_ids), unique_ids)
                    stage1_mem[idx_valid] <- as.integer(id_map[as.character(stage1_mem[idx_valid])])
                }
                names(stage1_mem) <- kept_genes
                stage1_membership_labels <- stage1_mem
            }
        }

        # The stored graph must be a view of the already thresholded Stage-1
        # similarity, not a second pass over raw L.  In particular, do not add
        # cross-module edges that bypass FDR or q95 filtering.
        A_sub_local <- stage1$similarity_sub %||% stage1$A_sub
        if (is.null(A_sub_local)) {
            stop("Stage-1 did not retain its filtered similarity matrix.", call. = FALSE)
        }
        if (!identical(rownames(A_sub_local), kept_genes)) {
            A_sub_local <- A_sub_local[kept_genes, kept_genes, drop = FALSE]
        }
        assigned_genes <- kept_genes[!is.na(stage1_membership_labels[kept_genes])]
        A_final <- A_sub_local[assigned_genes, assigned_genes, drop = FALSE]
        diag(A_final) <- 0
        edges_keep <- .arrind_from_matrix_predicate(
            A_final,
            op = "gt",
            cutoff = 0,
            triangle = "upper",
            keep_diag = FALSE
        )
        if (nrow(edges_keep)) {
            base_weight <- as.numeric(A_final[edges_keep])
            edge_weight <- if (isTRUE(config$use_log1p_weight)) log1p(base_weight) else base_weight
            valid_edge <- is.finite(edge_weight) & edge_weight > 0
            edges_keep <- edges_keep[valid_edge, , drop = FALSE]
            base_weight <- base_weight[valid_edge]
            edge_weight <- edge_weight[valid_edge]
            edges_df_final <- data.frame(
                from = assigned_genes[edges_keep[, 1]],
                to = assigned_genes[edges_keep[, 2]],
                weight = edge_weight,
                sign = ifelse(base_weight < 0, "neg", "pos"),
                stringsAsFactors = FALSE
            )
        } else {
            edges_df_final <- data.frame(
                from = character(0), to = character(0),
                weight = numeric(0), sign = character(0),
                stringsAsFactors = FALSE
            )
        }
        g_cons <- igraph::graph_from_data_frame(edges_df_final,
            directed = FALSE, vertices = assigned_genes
        )

        graph_views <- .build_consensus_graph_views(
            backbone_graph = g_cons,
            stage1_consensus_matrix = stage1_consensus_matrix,
            stage1_membership_matrix = stage1_membership_matrix,
            raw_similarity = similarity_matrix,
            consensus_thr = config$consensus_thr %||% 0.95,
            use_log1p_weight = isTRUE(config$use_log1p_weight)
        )

        return(list(
            membership = stage1_membership_labels,
            genes_all = genes_all,
            kept_genes = assigned_genes,
            consensus_graph = graph_views$backbone_graph,
            consensus_frequency_graph = graph_views$consensus_frequency_graph,
            paper_consensus_graph = graph_views$paper_consensus_graph,
            membership_matrix = stage1$stage1_membership_matrix,
            stage1_consensus = stage1_consensus_matrix,
            metrics_before = NULL,
            metrics_after = NULL,
            genes_before = NULL,
            genes_after = NULL,
            W = igraph::as_adjacency_matrix(g_cons, attr = "weight", sparse = TRUE)
        ))
    }

    L_post <- similarity_matrix[kept_genes, kept_genes, drop = FALSE]
    L_post[abs(L_post) < config$min_cutoff] <- 0
    if (config$use_significance) {
        if (is.null(FDR)) stop("Significance filtering requested but FDR is missing.", call. = FALSE)
        L_post <- .align_and_filter_fdr(L_post, similarity_matrix, FDR, config$significance_max)
    }
    L_post[L_post < 0] <- 0
    L_post <- pmax(L_post, t(L_post))
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
        mask <- if (config$CI_rule == "remove_within") (L_post <= f_hi(rMat)) else (L_post < f_lo(rMat) | L_post > f_hi(rMat))
        L_post[mask] <- 0
    }

    if (config$use_mh_weight) {
        if (is.null(mh_object)) stop("MH weighting requested but mh_object is NULL.")
        if (inherits(mh_object, "igraph")) {
            ed_sim <- igraph::as_data_frame(mh_object, what = "edges")
            wcol <- if ("MH" %in% names(ed_sim)) "MH" else if ("CMH" %in% names(ed_sim)) "CMH" else if ("weight" %in% names(ed_sim)) "weight" else stop("MH graph lacks weight attribute.")
            key_sim <- paste(pmin(ed_sim$from, ed_sim$to), pmax(ed_sim$from, ed_sim$to), sep = "|")
            sim_map <- setNames(ed_sim[[wcol]], key_sim)
            idx <- if (inherits(L_post, "sparseMatrix")) {
                ed <- summary(L_post)
                ed <- ed[ed$i < ed$j & ed$x > 0, , drop = FALSE]
                if (!nrow(ed)) matrix(integer(0), ncol = 2) else as.matrix(ed[, c("i", "j"), drop = FALSE])
            } else {
                which(upper.tri(L_post) & L_post > 0, arr.ind = TRUE)
            }
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
            idx <- if (inherits(L_post, "sparseMatrix")) {
                ed <- summary(L_post)
                ed <- ed[ed$i < ed$j & ed$x > 0, , drop = FALSE]
                if (!nrow(ed)) matrix(integer(0), ncol = 2) else as.matrix(ed[, c("i", "j"), drop = FALSE])
            } else {
                which(upper.tri(L_post) & L_post > 0, arr.ind = TRUE)
            }
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
        L_post <- .stage2_restore_filtered_backbone(
            L_post = L_post,
            stage1_similarity = stage1$similarity_sub %||% stage1$A_sub,
            kept_genes = kept_genes,
            stage1_membership_labels = stage1_membership_labels,
            backbone_floor_q = config$backbone_floor_q
        )
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

    use_consensus <- isTRUE(config$use_consensus)
    consensus_thr <- config$consensus_thr
    n_restart <- max(1L, as.integer(config$n_restart))
    n_threads <- if (is.null(config$n_threads)) 1L else max(1L, as.integer(config$n_threads))
    objective <- config$objective
    res_param_base <- if (identical(objective, "modularity")) as.numeric(config$gamma) else as.numeric(config$resolution)
    if (!is.finite(res_param_base)) res_param_base <- 1
    mode_setting <- config$mode
    large_n_threshold <- if (is.null(config$large_n_threshold)) length(kept_genes) else as.integer(config$large_n_threshold)
    aggr_future_workers <- if (is.null(config$aggr_future_workers)) 1L else max(1L, as.integer(config$aggr_future_workers))
    aggr_batch_size <- config$aggr_batch_size
    nk_leiden_iterations <- if (is.null(config$nk_leiden_iterations)) 10L else as.integer(config$nk_leiden_iterations)
    nk_leiden_randomize <- if (is.null(config$nk_leiden_randomize)) TRUE else isTRUE(config$nk_leiden_randomize)
    min_cluster_keep <- if (is.null(config$min_cluster_size)) 2L else max(1L, as.integer(config$min_cluster_size))
    sub_min_size <- max(1L, as.integer(config$sub_min_size))
    sub_min_child_size <- max(1L, as.integer(config$sub_min_child_size))
    sub_resolution_factor <- if (is.null(config$sub_resolution_factor)) 1.3 else config$sub_resolution_factor
    sub_within_cons_max <- config$sub_within_cons_max
    sub_conductance_min <- config$sub_conductance_min
    sub_improve_within_cons_min <- config$sub_improve_within_cons_min
    sub_max_groups <- max(1L, as.integer(config$sub_max_groups))
    enable_qc_filter <- isTRUE(config$enable_qc_filter)
    qc_gene_intra_cons_min <- config$qc_gene_intra_cons_min
    qc_gene_best_out_cons_min <- config$qc_gene_best_out_cons_min
    qc_gene_intra_weight_q <- config$qc_gene_intra_weight_q
    keep_cross_stable <- isTRUE(config$keep_cross_stable)
    min_cutoff <- config$min_cutoff

    ng <- length(kept_genes)
    memb_final <- rep(NA_integer_, ng)
    names(memb_final) <- kept_genes
    cons <- Matrix(0, ng, ng, sparse = TRUE, dimnames = list(kept_genes, kept_genes))

    cl_ids_stage1 <- sort(na.omit(unique(stage1_membership_labels)))
    next_global_id <- 1L
    for (cid in cl_ids_stage1) {
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
                            stop("Stage-2 subcluster graph algorithm failed: ", as.character(res), call. = FALSE)
                        } else {
                            mm <- try(igraph::membership(res), silent = TRUE)
                            if (inherits(mm, "try-error") || is.null(mm)) {
                                stop("Stage-2 subcluster membership extraction failed: ", as.character(mm), call. = FALSE)
                            }
                            mm
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
                .log_backend("clusterGenes", "S07", "stage2_consensus_backend", "single-run",
                    reason = "use_consensus=FALSE", verbose = verbose
                )
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
                    .log_backend("clusterGenes", "S07", "stage2_consensus_backend", "cpp",
                        reason = "consensus_sparse", verbose = verbose
                    )
                } else {
                    .log_backend("clusterGenes", "S07", "stage2_consensus_backend", "r",
                        reason = "consensus_sparse_failed", verbose = verbose
                    )
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
                    .log_backend("clusterGenes", "S07", "stage2_consensus_graph_backend", backend_label,
                        reason = paste0("stage2_backend=", backend_sub), verbose = verbose
                    )
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

    metrics_before <- .summarise_cluster_metrics(memb_final, cons, W)
    genes_before <- .summarise_gene_membership(memb_final, cons, W, kept_genes)

    if (enable_qc_filter) {
        q_w <- quantile(genes_before$w_in, probs = qc_gene_intra_weight_q, na.rm = TRUE)
        bad <- which(!is.na(genes_before$p_in) &
            genes_before$p_in < qc_gene_intra_cons_min &
            genes_before$p_best_out > qc_gene_best_out_cons_min &
            genes_before$w_in <= q_w)
        if (length(bad)) {
            memb_final[match(genes_before$gene[bad], kept_genes)] <- NA
        }
    }

    if (config$enable_subcluster && !identical(stage2_algo_final, "hotspot-like")) {
        next_id <- if (length(na.omit(memb_final))) max(na.omit(memb_final)) + 1L else 1L
        for (row_idx in seq_len(nrow(metrics_before))) {
            cid <- metrics_before$cluster[row_idx]
            size <- metrics_before$size[row_idx]
            wcons <- metrics_before$within_cons[row_idx]
            condc <- metrics_before$conductance[row_idx]
            if (is.na(size) || size < sub_min_size || is.na(wcons) || is.na(condc)) next
            if (!(wcons < sub_within_cons_max && condc > sub_conductance_min)) next

            genes_c <- kept_genes[memb_final == cid]
            if (!length(genes_c)) next
            g_sub <- igraph::induced_subgraph(corrected_similarity_graph, vids = genes_c)
            if (igraph::ecount(g_sub) == 0) next
            ew <- igraph::E(g_sub)$weight
            el <- igraph::as_data_frame(g_sub, what = "edges")[, c("from", "to")]
            memb_mat_sub <- future.apply::future_sapply(
                seq_len(min(n_restart, 50)),
                function(ii) {
                    g_local <- igraph::graph_from_data_frame(cbind(el, weight = ew), directed = FALSE, vertices = igraph::V(g_sub)$name)
                    res <- .run_graph_algorithm(
                        g_local,
                        stage2_algo_final,
                        if (identical(objective, "modularity")) res_param_base * sub_resolution_factor else res_param_base * sub_resolution_factor,
                        objective
                    )
                    igraph::membership(res)
                },
                future.seed = TRUE
            )
            if (is.null(dim(memb_mat_sub))) memb_mat_sub <- matrix(memb_mat_sub, ncol = 1)
            rownames(memb_mat_sub) <- igraph::V(g_sub)$name
            cons_sub <- .compute_consensus_sparse_r(memb_mat_sub, thr = consensus_thr)
            if (nnzero(cons_sub) == 0) next
            g_cons_sub <- igraph::graph_from_adjacency_matrix(cons_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
            comm_sub <- .run_graph_algorithm(
                g_cons_sub,
                stage2_algo_final,
                if (identical(objective, "modularity")) res_param_base * sub_resolution_factor else res_param_base * sub_resolution_factor,
                objective
            )
            memb_sub <- igraph::membership(comm_sub)

            sub_ids <- sort(unique(memb_sub))
            if (length(sub_ids) <= 1 || length(sub_ids) > sub_max_groups) next
            sizes <- vapply(sub_ids, function(sid) sum(memb_sub == sid), integer(1))
            if (any(sizes < sub_min_child_size)) next

            wcons_sub <- vapply(sub_ids, function(sid) {
                idx_local <- which(memb_sub == sid)
                if (length(idx_local) > 1) mean(cons_sub[idx_local, idx_local][upper.tri(cons_sub[idx_local, idx_local])], na.rm = TRUE) else 0
            }, numeric(1))
            wcons_old <- wcons
            wcons_new <- sum(wcons_sub * sizes) / sum(sizes)
            if (!(is.finite(wcons_new) && wcons_new >= (wcons_old + sub_improve_within_cons_min))) next

            for (sid in sub_ids) {
                genes_sid <- names(memb_sub)[memb_sub == sid]
                memb_final[genes_sid] <- next_id
                next_id <- next_id + 1L
            }
            block_names <- intersect(rownames(cons_sub), genes_c)
            if (length(block_names)) {
                cons[block_names, block_names] <- as.matrix(cons_sub[block_names, block_names])
            }
        }
    }

    di <- rowSums((L_post > 0)[kept_genes, kept_genes, drop = FALSE] &
        outer(memb_final, memb_final, "=="), na.rm = TRUE)
    memb_final[di == 0] <- NA
    tbl <- table(memb_final, useNA = "no")
    small <- names(tbl)[tbl < min_cluster_keep]
    if (length(small)) {
        memb_final[memb_final %in% small] <- NA
    }
    tbl <- table(memb_final, useNA = "no")
    if (length(tbl)) {
        unique_ids <- sort(as.integer(names(tbl)))
        id_map <- setNames(seq_along(unique_ids), unique_ids)
        idx_keep <- which(!is.na(memb_final))
        if (length(idx_keep)) {
            memb_final[idx_keep] <- as.integer(id_map[as.character(memb_final[idx_keep])])
        }
        tbl <- table(memb_final, useNA = "no")
        stable_cl <- as.integer(names(tbl)[tbl >= min_cluster_keep])
    } else {
        stable_cl <- integer(0)
    }
    names(memb_final) <- kept_genes

    metrics_after <- .summarise_cluster_metrics(memb_final, cons, W)
    genes_after <- .summarise_gene_membership(memb_final, cons, W, kept_genes)

    cons_ts <- as(cons, "TsparseMatrix")
    if (length(cons_ts@x)) {
        sm2 <- data.frame(
            i = cons_ts@i + 1L,
            j = cons_ts@j + 1L,
            x = cons_ts@x,
            stringsAsFactors = FALSE
        )
    } else {
        sm2 <- data.frame(i = integer(0), j = integer(0), x = numeric(0), stringsAsFactors = FALSE)
    }
    sm2 <- sm2[sm2$i < sm2$j & sm2$x > 0, , drop = FALSE]

    key_corr <- paste(pmin(edges_corr$from, edges_corr$to), pmax(edges_corr$from, edges_corr$to), sep = "|")
    w_map <- setNames(edges_corr$weight, key_corr)
    intra_keys <- paste(
        pmin(kept_genes[sm2$i], kept_genes[sm2$j]),
        pmax(kept_genes[sm2$i], kept_genes[sm2$j]),
        sep = "|"
    )
    w_vec <- as.numeric(w_map[intra_keys])
    keep_mask <- !is.na(w_vec) & w_vec > 0
    intra_df <- data.frame(
        from = kept_genes[sm2$i[keep_mask]],
        to = kept_genes[sm2$j[keep_mask]],
        weight = w_vec[keep_mask],
        sign = ifelse(similarity_matrix[cbind(
            match(kept_genes[sm2$i[keep_mask]], rownames(similarity_matrix)),
            match(kept_genes[sm2$j[keep_mask]], colnames(similarity_matrix))
        )] < 0, "neg", "pos"),
        stringsAsFactors = FALSE
    )

    cross_df <- NULL
    if (keep_cross_stable && length(stable_cl) > 1) {
        idx_stable <- memb_final %in% stable_cl
        genes_stable <- kept_genes[idx_stable]
        if (length(genes_stable) > 1) {
	            L_sub_raw <- L_post[genes_stable, genes_stable, drop = FALSE]
	            ij <- .arrind_from_matrix_predicate(
	                L_sub_raw,
	                op = "ge",
	                cutoff = min_cutoff,
	                triangle = "upper",
	                keep_diag = FALSE
	            )
            if (nrow(ij)) {
                g1 <- genes_stable[ij[, 1]]
                g2 <- genes_stable[ij[, 2]]
                cl1 <- memb_final[g1]
                cl2 <- memb_final[g2]
                keep_cross <- cl1 != cl2 & !is.na(cl1) & !is.na(cl2)
                if (any(keep_cross)) {
                    g1 <- g1[keep_cross]
                    g2 <- g2[keep_cross]
                    Lval <- L_sub_raw[cbind(match(g1, genes_stable), match(g2, genes_stable))]
                    wL <- if (use_log1p_weight) log1p(Lval) else Lval
                    mask <- wL > 0
                    if (any(mask)) {
                        cross_df <- data.frame(
                            from = g1[mask],
                            to = g2[mask],
                            weight = wL[mask],
                            sign = ifelse(Lval[mask] < 0, "neg", "pos")
                        )
                    }
                }
            }
        }
    }

    edges_df_final <- if (!is.null(cross_df)) unique(rbind(intra_df, cross_df)) else intra_df
    if (!nrow(edges_df_final)) {
        cons_graph_final <- igraph::make_empty_graph(n = length(kept_genes), directed = FALSE)
        cons_graph_final <- igraph::set_vertex_attr(cons_graph_final, "name", value = kept_genes)
    } else {
        cons_graph_final <- igraph::graph_from_data_frame(edges_df_final,
            directed = FALSE, vertices = kept_genes
        )
    }

    graph_views <- .build_consensus_graph_views(
        backbone_graph = cons_graph_final,
        stage1_consensus_matrix = stage1_consensus_matrix,
        stage1_membership_matrix = stage1_membership_matrix,
        raw_similarity = similarity_matrix,
        consensus_thr = consensus_thr %||% 0.95,
        use_log1p_weight = use_log1p_weight
    )

    list(
        membership = memb_final,
        kept_genes = kept_genes,
        genes_all = genes_all,
        consensus_graph = graph_views$backbone_graph,
        consensus_frequency_graph = graph_views$consensus_frequency_graph,
        paper_consensus_graph = graph_views$paper_consensus_graph,
        membership_matrix = stage1_membership_matrix,
        stage1_consensus = stage1_consensus_matrix,
        metrics_before = metrics_before,
        metrics_after = metrics_after,
        genes_before = genes_before,
        genes_after = genes_after,
        corrected_graph = corrected_similarity_graph,
        weights_matrix = W,
        stage2_backend = stage2_backend
    )
}

# Compatibility shim kept as a dynamic forwarder so it cannot retain a stale
# closure if the active Stage-2 implementation is changed during package load.
.stage2_refine_workflow <- function(...) {
    .stage2_refine_workflow_v2(...)
}

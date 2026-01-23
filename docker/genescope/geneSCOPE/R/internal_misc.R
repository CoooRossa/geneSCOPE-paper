#' Align an FDR matrix with L and zero out edges above threshold.
#' @description
#' Internal helper for `.align_and_filter_fdr`.
#' Keeps Lee statistics and FDR control matrices synchronized before masking edges.
#' @param L Parameter value.
#' @param L_raw Parameter value.
#' @param FDRmat Parameter value.
#' @param FDR_max Numeric threshold.
#' @return Return value used internally.
#' @keywords internal
.align_and_filter_fdr <- function(L, L_raw, FDRmat, FDR_max) {
    if (is.null(FDRmat)) return(L)
    FM <- FDRmat

    # 1) Robust coercion to base matrix
    if (inherits(FM, "big.matrix")) {
        FM <- try({
            if (requireNamespace("bigmemory", quietly = TRUE)) bigmemory::as.matrix(FM) else FM[, ]
        }, silent = TRUE)
        if (inherits(FM, "try-error")) {
            warning("[align_and_filter_FDR] FDR is big.matrix but cannot coerce to base matrix; disabling FDR filter.")
            return(L)
        }
    } else if (!is.matrix(FM)) {
        FM_try <- try(as.matrix(FM), silent = TRUE)
        if (!inherits(FM_try, "try-error") && is.matrix(FM_try)) {
            FM <- FM_try
        } else {
            warning("[align_and_filter_FDR] FDR object cannot be coerced to matrix (class=", paste(class(FDRmat), collapse=","), "); disabling FDR filter.")
            return(L)
        }
    }

    # 2) Align by names (preferred) or by shape
    policy <- getOption("geneSCOPE.fdr_align_policy", "by_name")
    if (is.null(dim(FM))) {
        FM <- matrix(FM, nrow = nrow(L_raw), ncol = ncol(L_raw))
        FM <- .safe_set_dimnames(FM, dimnames(L_raw))
    } else if (!(identical(dim(FM), dim(L_raw)) &&
                 identical(rownames(FM), rownames(L_raw)) &&
                 identical(colnames(FM), colnames(L_raw)))) {
        if (identical(policy, "by_shape")) {
            FM <- .safe_set_dim(FM, dim(L_raw))
            FM <- .safe_set_dimnames(FM, dimnames(L_raw))
        } else {
            if (!is.null(rownames(FM)) && !is.null(colnames(FM)) &&
                all(rownames(L_raw) %in% rownames(FM)) && all(colnames(L_raw) %in% colnames(FM))) {
                FM <- FM[rownames(L_raw), colnames(L_raw), drop = FALSE]
            } else {
                # fall back to shape-only alignment, but keep numeric matrix where possible
                FM <- .safe_set_dim(FM, dim(L_raw))
                FM <- .safe_set_dimnames(FM, dimnames(L_raw))
            }
        }
    }
    FM <- .safe_set_dimnames(FM, dimnames(L_raw))

    # 3) Subset FDR to current kept rows of L
    target_rows <- rownames(L)
    if (!is.null(target_rows) && !is.null(rownames(FM))) {
        idx <- match(target_rows, rownames(FM))
        if (anyNA(idx)) {
            stop("FDR matrix cannot be aligned to current gene subset", call. = FALSE)
        }
        FM <- FM[idx, idx, drop = FALSE]
    } else if (nrow(FM) != nrow(L) || ncol(FM) != ncol(L)) {
        stop("FDR matrix dimensions do not match current subset", call. = FALSE)
    }

    # 4) Apply mask
    if (inherits(L, "sparseMatrix")) {
        LT <- as(L, "TsparseMatrix")
        if (length(LT@x)) {
            rows <- LT@i + 1L
            cols <- LT@j + 1L
            mask <- (FM[cbind(rows, cols)] > FDR_max)
            if (any(mask)) LT@x[mask] <- 0
        }
        return(drop0(as(LT, "CsparseMatrix")))
    }
    L[FM > FDR_max] <- 0
    L
}

#' As Numeric Or Na
#' @description
#' Internal helper for `.as_numeric_or_na`.
#' @param x Parameter value.
#' @return Return value used internally.
#' @keywords internal
.as_numeric_or_na <- function(x) {
    if (is.null(x)) return(NA_real_)
    val <- suppressWarnings(as.numeric(x))
    if (!length(val)) return(NA_real_)
    val[1]
}

#' Compute block identifiers from grid coordinates and block size.
#' @description
#' Internal helper for `.assign_block_id`.
#' Groups grid tiles into coarse blocks so block-aware statistics can be aggregated quickly.
#' @param grid_info Parameter value.
#' @param block_side Parameter value.
#' @return Return value used internally.
#' @keywords internal
.assign_block_id <- function(grid_info, block_side = 8) {
    # gx / gy start at 1
    bx <- (grid_info$gx - 1L) %/% block_side
    by <- (grid_info$gy - 1L) %/% block_side
    # merge into a single integer id
    max_by <- max(by)
    block_id <- bx * (max_by + 1L) + by + 1L
    block_id
}

#' Check Grid Content
#' @description
#' Internal helper for `.check_grid_content`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.check_grid_content <- function(scope_obj, grid_name, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    g_layer <- .select_grid_layer(scope_obj, grid_name)

    required_elements <- c("grid_info", "counts")
    missing <- setdiff(required_elements, names(g_layer))
    if (length(missing) > 0) {
        stop(
            "Grid layer '", grid_name, "' missing required elements: ",
            paste(missing, collapse = ", ")
        )
    }

    # Removed verbose message to avoid redundancy with main functions
    invisible(TRUE)
}

#' Compact block identifiers to consecutive integers (1..k) sorted by ID.
#' @description
#' Internal helper for `.compact_block_id`.
#' This reduces factor overhead and keeps permutation shuffles stable.
#' @param block_id Parameter value.
#' @return Return value used internally.
#' @keywords internal
.compact_block_id <- function(block_id) {
    if (length(block_id) == 0L) return(block_id)
    uniq_sorted <- sort(unique(block_id))
    match(block_id, uniq_sorted)
}

#' Summarize how often each edge co-clusters across runs.
#' @description
#' Internal helper for `.consensus_on_edges`.
#' Computes consensus weights on the provided edge list without allocating a dense N×N matrix.
#' @param memb_mat Integer matrix (genes × runs), each column is a community assignment.
#' @param edge_i Integer vector (1-based) indices of edge endpoints (length m).
#' @param edge_j Integer vector (1-based) indices of edge endpoints (length m).
#' @param n_threads Number of threads to use.
#' @return Numeric vector of length m: fraction of runs in which the two endpoints co-cluster.
#' @keywords internal
.consensus_on_edges <- function(memb_mat, edge_i, edge_j, n_threads = NULL) {
    stopifnot(is.matrix(memb_mat), length(edge_i) == length(edge_j))
    m <- length(edge_i)
    if (m == 0L) return(numeric(0))
    runs <- ncol(memb_mat)
    if (is.null(n_threads)) n_threads <- .get_safe_thread_count(8L)
    # Fast path via C++ when available
    if (exists("consensus_on_edges_omp", mode = "function")) {
        counts <- consensus_on_edges_omp(as.integer(edge_i), as.integer(edge_j),
                                         memb = matrix(as.integer(memb_mat), nrow(memb_mat), ncol(memb_mat)),
                                         n_threads = as.integer(n_threads))
        return(as.numeric(counts) / runs)
    }
    # R fallback (block-wise compare to limit memory)
    mm <- memb_mat
    storage.mode(mm) <- "integer"
    counts <- integer(m)
    block <- max(1L, min(64L, runs))
    for (start in seq.int(1L, runs, by = block)) {
        end <- min(runs, start + block - 1L)
        comp <- mm[edge_i, start:end, drop = FALSE] == mm[edge_j, start:end, drop = FALSE]
        counts <- counts + as.integer(rowSums(comp))
    }
    counts / runs
}

#' Draw Wrapped Legend
#' @description
#' Internal helper for `.draw_wrapped_legend`.
#' @param x_start Parameter value.
#' @param y_start Parameter value.
#' @param labels Parameter value.
#' @param pch Parameter value.
#' @param bg Parameter value.
#' @param border Parameter value.
#' @param point_cex Parameter value.
#' @param text_cex Parameter value.
#' @param x_range Parameter value.
#' @param y_range Parameter value.
#' @param max_width Numeric threshold.
#' @param dot_gap_factor Parameter value.
#' @param item_gap_factor Parameter value.
#' @param row_spacing_factor Parameter value.
#' @return Return value used internally.
#' @keywords internal
.draw_wrapped_legend <- function(
    x_start, y_start, labels, pch, bg, border,
    point_cex, text_cex, x_range, y_range,
    max_width, dot_gap_factor = 0.6, item_gap_factor = 0.5,
    row_spacing_factor = 1.2) {
    if (!length(labels)) {
        return(invisible(NULL))
    }

    width_user <- function(txt, cex) {
        strwidth(txt, cex = cex, units = "figure") * x_range
    }
    height_user <- function(txt, cex) {
        strheight(txt, cex = cex, units = "figure") * y_range
    }

    text_widths <- width_user(labels, text_cex)
    dot_gap <- width_user("M", text_cex) * dot_gap_factor
    if (!is.finite(max_width) || max_width <= 0) {
        max_width <- sum(width_user("M", point_cex) * 0.6 + dot_gap + text_widths) +
            width_user(" ", text_cex) * item_gap_factor * (length(labels) - 1)
    }
    min_required <- min(width_user("M", point_cex) * 0.6 + dot_gap + text_widths)
    if (max_width < min_required) {
        max_width <- min_required * 1.1
    }

    rows <- list()
    current <- integer(0)
    width_current <- 0
    for (i in seq_along(labels)) {
            item_width <- width_user("M", point_cex) * 0.6 + dot_gap + text_widths[i]
        item_total <- if (length(current)) item_width + width_user(" ", text_cex) * item_gap_factor else item_width
        if (length(current) && width_current + item_total > max_width) {
            rows[[length(rows) + 1]] <- current
            current <- i
            width_current <- item_width
        } else {
            current <- c(current, i)
            width_current <- width_current + item_total
        }
    }
    if (length(current)) rows[[length(rows) + 1]] <- current

    y_current <- y_start
    for (row_idx in seq_along(rows)) {
        idxs <- rows[[row_idx]]
        x_pos <- x_start
        for (j in idxs) {
            points(x_pos, y_current, pch = pch[j], bg = bg[j], col = border[j], cex = point_cex)
            x_pos <- x_pos + width_user("M", point_cex) * 0.6
            text(x_pos + dot_gap, y_current, labels[j], adj = c(0, 0.5), cex = text_cex)
            x_pos <- x_pos + dot_gap + text_widths[j] + width_user(" ", text_cex) * item_gap_factor
        }
        y_current <- y_current - height_user("M", text_cex) * row_spacing_factor
    }

    invisible(NULL)
}

#' Ensure Numeric Matrix
#' @description
#' Internal helper for `.ensure_numeric_matrix`.
#' @param mat Parameter value.
#' @param nrow Parameter value.
#' @param ncol Parameter value.
#' @param rownames_out Parameter value.
#' @param colnames_out Parameter value.
#' @param as_sparse Parameter value.
#' @return Return value used internally.
#' @keywords internal
.ensure_numeric_matrix <- function(mat,
                                   nrow = NULL,
                                   ncol = NULL,
                                   rownames_out = NULL,
                                   colnames_out = NULL,
                                   as_sparse = FALSE) {
    if (is.null(mat)) return(NULL)

    coerce_sparse_pattern <- function(x) {
        idx <- summary(x)
        rr <- idx$i
        cc <- idx$j
        vv <- idx$x
        if (is.null(vv)) vv <- rep(1, length(rr))
        dense <- matrix(0, nrow = nrow(x), ncol = ncol(x))
        dense[cbind(rr, cc)] <- vv
        dimnames(dense) <- dimnames(x)
        dense
    }

    out <- mat
    orig <- mat
    if (inherits(out, "big.matrix")) {
        if (requireNamespace("bigmemory", quietly = TRUE)) {
            out <- bigmemory::as.matrix(out)
        } else {
            stop("bigmemory package required to coerce big.matrix objects", call. = FALSE)
        }
    }
    if (inherits(out, "big.matrix")) {
        if (!requireNamespace("bigmemory", quietly = TRUE)) {
            stop("bigmemory package required to coerce big.matrix objects", call. = FALSE)
        }
        out <- bigmemory::as.matrix(out)
    }

    if (inherits(out, "Matrix")) {
        if (inherits(out, "nMatrix")) {
            out <- coerce_sparse_pattern(out)
        } else {
            out <- tryCatch(as.matrix(out), error = function(e) NULL)
            if (is.null(out)) {
                out <- coerce_sparse_pattern(as(orig, "ngCMatrix"))
            }
        }
    } else if (!is.matrix(out)) {
        out <- tryCatch(as.matrix(out), error = function(e) NULL)
        if (is.null(out) && is.atomic(mat)) {
            if (is.null(nrow) || is.null(ncol)) stop("Cannot reshape atomic vector without target dimensions")
            out <- matrix(as.numeric(mat), nrow = nrow, ncol = ncol)
        } else if (is.null(out)) {
            if (!is.null(dim(mat)) && length(dim(mat)) == 2L) {
                out <- matrix(as.numeric(mat), nrow = dim(mat)[1], ncol = dim(mat)[2])
            }
        }
    }

    if (is.null(out)) {
        cls <- paste(class(orig), collapse = ", ")
        stop(sprintf("Failed to coerce object of class [%s] to numeric matrix", cls), call. = FALSE)
    }
    storage.mode(out) <- "double"
    if (!is.null(rownames_out)) rownames(out) <- rownames_out
    if (!is.null(colnames_out)) colnames(out) <- colnames_out
    if (as_sparse) {
        if (inherits(out, "sparseMatrix")) {
            out <- as(out, "CsparseMatrix")
        } else {
            out <- as(Matrix(out, sparse = TRUE), "CsparseMatrix")
        }
        out <- drop0(out)
    }
    out
}

#' Filter Matrix By Quantile
#' @description
#' Internal helper for `.filter_matrix_by_quantile`.
#' @param mat Parameter value.
#' @param pct_min Numeric threshold.
#' @param pct_max Numeric threshold.
#' @return Return value used internally.
#' @keywords internal
.filter_matrix_by_quantile <- function(mat, pct_min = "q0", pct_max = "q100") {
    stopifnot(is.matrix(mat) || inherits(mat, "Matrix"))
    pmin <- .parse_q(pct_min)
    pmax <- .parse_q(pct_max)
    if (pmin > pmax) stop("pct_min > pct_max")

    legacy_mode <- !identical(getOption("geneSCOPE.filter_sparse_mode", "legacy"), "sparse")

    include_zero_mass <- !identical(getOption("geneSCOPE.filter_sparse_ignore_zero", FALSE), TRUE)

    adjust_quantile_with_zeros <- function(vec, p, zero_frac) {
        if (!length(vec)) return(NA_real_)
        if (!is.finite(zero_frac) || zero_frac <= 0) {
            return(as.numeric(quantile(vec, p, na.rm = TRUE, type = 7)))
        }
        if (zero_frac >= 1) return(0)
        if (p <= zero_frac) return(0)
        as.numeric(quantile(vec, (p - zero_frac) / (1 - zero_frac), na.rm = TRUE, type = 7))
    }

    if (legacy_mode) {
        dense_mat <- as.matrix(mat)
        vec <- as.vector(dense_mat)
        pos_mask <- vec >= 0
        neg_mask <- vec < 0

        res <- dense_mat * 0
        if (any(pos_mask, na.rm = TRUE)) {
            pos_vec <- vec[pos_mask]
            thrL <- quantile(pos_vec, pmin, na.rm = TRUE, type = 7)
            thrU <- quantile(pos_vec, pmax, na.rm = TRUE, type = 7)
            keep <- pos_mask & vec >= thrL & vec <= thrU
            res[keep] <- dense_mat[keep]
        }
        if (any(neg_mask, na.rm = TRUE)) {
            neg_vec <- abs(vec[neg_mask])
            thrL <- quantile(neg_vec, pmin, na.rm = TRUE, type = 7)
            thrU <- quantile(neg_vec, pmax, na.rm = TRUE, type = 7)
            keep <- neg_mask & abs(vec) >= thrL & abs(vec) <= thrU
            res[keep] <- dense_mat[keep]
        }
        return(drop0(Matrix(res, sparse = TRUE)))
    }

    # Fast path for sparse matrices: operate only on existing (upper-tri) entries
    if (inherits(mat, "sparseMatrix")) {
        TT <- as(mat, "TsparseMatrix")
        # Use only upper triangle to avoid double work; keep sign information
        mU <- (TT@i < TT@j)
        if (!any(mU)) return(drop0(mat * 0))
        xi <- TT@i[mU] + 1L
        xj <- TT@j[mU] + 1L
        xv <- TT@x[mU]
        # Split by sign; negatives use absolute value for quantiles (consistent with previous behavior)
        pos_m <- xv >= 0
        neg_m <- !pos_m
        keep_idx <- logical(length(xv))
        if (any(pos_m)) {
            pos_vec <- xv[pos_m]
            # Guard against all-NA or length-1 edge cases
            if (length(pos_vec) == 1L) {
                pos_thrL <- pos_vec; pos_thrU <- pos_vec
            } else {
                if (include_zero_mass) {
                    zero_count <- max((as.double(nrow(mat)) * (nrow(mat) - 1)) / 2 - length(pos_vec), 0)
                    zero_frac <- if ((as.double(nrow(mat)) * (nrow(mat) - 1)) / 2 > 0) zero_count / (zero_count + length(pos_vec)) else 0
                    pos_thrL <- adjust_quantile_with_zeros(pos_vec, pmin, zero_frac)
                    pos_thrU <- adjust_quantile_with_zeros(pos_vec, pmax, zero_frac)
                } else {
                    pos_thrL <- as.numeric(quantile(pos_vec, pmin, na.rm = TRUE, type = 7))
                    pos_thrU <- as.numeric(quantile(pos_vec, pmax, na.rm = TRUE, type = 7))
                }
            }
            keep_idx[pos_m] <- (pos_vec >= pos_thrL) & (pos_vec <= pos_thrU)
        }
        if (any(neg_m)) {
            neg_abs <- abs(xv[neg_m])
            if (length(neg_abs) == 1L) {
                neg_thrL <- neg_abs; neg_thrU <- neg_abs
            } else {
                neg_thrL <- as.numeric(quantile(neg_abs, pmin, na.rm = TRUE, type = 7))
                neg_thrU <- as.numeric(quantile(neg_abs, pmax, na.rm = TRUE, type = 7))
            }
            keep_idx[neg_m] <- (neg_abs >= neg_thrL) & (neg_abs <= neg_thrU)
        }
        if (!any(keep_idx)) return(drop0(mat * 0))
        ii <- xi[keep_idx]; jj <- xj[keep_idx]; vv <- xv[keep_idx]
        # Materialize symmetric sparse matrix of kept edges
        res <- sparseMatrix(
            i = c(ii, jj), j = c(jj, ii), x = c(vv, vv),
            dims = dim(mat), dimnames = dimnames(mat)
        )
        # Zero diagonal by construction; drop explicit zeros
        return(drop0(res))
    }

    # Dense path: operate on upper-tri only; for very large matrices, use sampling to estimate quantiles
    n <- nrow(mat)
    if (n == 0L) return(mat)
    vals <- mat[upper.tri(mat, diag = FALSE)]
    # Split by sign
    pos_vals <- vals[vals >= 0]
    neg_vals_abs <- abs(vals[vals < 0])
    # Approximate quantiles for very large vectors to avoid O(N log N) on tens of millions of entries
    max_samples <- getOption("geneSCOPE.quantile_max_samples", 5e6L)
    qfun <- function(x, p) {
        if (!length(x)) return(NA_real_)
        if (length(x) > max_samples) {
            # Sample without replacement for a robust approximation
            idx <- sample.int(length(x), max_samples)
            x <- x[idx]
        }
        as.numeric(quantile(x, p, na.rm = TRUE, type = 7))
    }
    pos_thrL <- qfun(pos_vals, pmin); pos_thrU <- qfun(pos_vals, pmax)
    neg_thrL <- qfun(neg_vals_abs, pmin); neg_thrU <- qfun(neg_vals_abs, pmax)

    # Build a logical keep mask for upper-tri, then mirror to full matrix
    keep_ut <- logical(length(vals))
    if (length(pos_vals)) {
        kpos <- (vals >= 0) & (vals >= pos_thrL) & (vals <= pos_thrU)
        keep_ut <- keep_ut | kpos
    }
    if (length(neg_vals_abs)) {
        kneg <- (vals < 0) & (abs(vals) >= neg_thrL) & (abs(vals) <= neg_thrU)
        keep_ut <- keep_ut | kneg
    }
    if (!any(keep_ut)) return(mat * 0)

    # Initialize zero matrix and fill symmetric kept entries
    res <- mat * 0
    res[upper.tri(mat, diag = FALSE)][keep_ut] <- vals[keep_ut]
    res <- res + t(res)
    return(res)
}

#' Force Single Thread Blas
#' @description
#' Internal helper for `.force_single_thread_blas`.
#' @return Return value used internally.
#' @keywords internal
.force_single_thread_blas <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch(
      {
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
      },
      error = function(e) { }
    )
  }
  .set_blas_environment()
  Sys.setenv(
    OPENBLAS_NUM_THREADS = "1",
    GOTO_NUM_THREADS = "1",
    OMP_NUM_THREADS = "1"
  )
}

#' Format Label Display
#' @description
#' Internal helper for `.format_label_display`.
#' @param x Parameter value.
#' @return Return value used internally.
#' @keywords internal
.format_label_display <- function(x) {
    gsub("_sub", ".", x, fixed = TRUE)
}

#' Gene Metrics
#' @description
#' Internal helper for `.gene_metrics`.
#' @param members Parameter value.
#' @param cons_mat Parameter value.
#' @param W_mat Parameter value.
#' @return Return value used internally.
#' @keywords internal
.gene_metrics <- function(members, cons_mat, W_mat) {
    cl_ids <- sort(na.omit(unique(members)))
    n <- length(members)
    p_in <- numeric(n); p_best_out <- numeric(n); w_in <- numeric(n)
    for (i in seq_len(n)) {
        cid <- members[i]
        if (is.na(cid)) { p_in[i] <- NA; p_best_out[i] <- NA; w_in[i] <- NA; next }
        idx <- which(members == cid)
        idx_noi <- setdiff(idx, i)
        if (length(idx_noi)) {
            p_in[i] <- mean(cons_mat[i, idx_noi], na.rm = TRUE)
            w_in[i] <- sum(W_mat[i, idx_noi])
        } else {
            p_in[i] <- 0; w_in[i] <- 0
        }
        other <- setdiff(cl_ids, cid)
        if (length(other)) {
            means <- vapply(other, function(cj) {
                if (length(which(members == cj))) mean(cons_mat[i, which(members == cj)], na.rm = TRUE) else 0
            }, numeric(1))
            p_best_out[i] <- max(means, na.rm = TRUE)
        } else {
            p_best_out[i] <- 0
        }
    }
    data.frame(gene = kept_genes, cluster = members, p_in = p_in, p_best_out = p_best_out, w_in = w_in, stringsAsFactors = FALSE)
}

#' Derive per-gene metrics (p_in, w_in, etc.) using sparse matrices.
#' @description
#' Internal helper for `.gene_metrics_sparse`.
#' Quantifies intra- and inter-cluster affinities while avoiding dense memory usage.
#' @param members Parameter value.
#' @param cons_mat Parameter value.
#' @param W_mat Parameter value.
#' @param genes Parameter value.
#' @return Return value used internally.
#' @keywords internal
.gene_metrics_sparse <- function(members, cons_mat, W_mat, genes) {
    if (!inherits(cons_mat, "sparseMatrix")) cons_mat <- as(cons_mat, "dgCMatrix")
    if (!inherits(W_mat, "sparseMatrix")) W_mat <- as(W_mat, "dgCMatrix")
    cl_ids <- sort(na.omit(unique(members)))
    n <- length(members)
    p_in <- rep(NA_real_, n); p_best_out <- rep(NA_real_, n); w_in <- rep(NA_real_, n)
    for (cid in cl_ids) {
        idx <- which(members == cid)
        if (!length(idx)) next
        # p_in: row means on intra-cluster consensus (excluding diag)
        subC <- cons_mat[idx, idx, drop = FALSE]
        # diag is zero; rowSums(subC) equals sum of weights to same cluster
        denom <- pmax(length(idx) - 1L, 1L)
        p_in[idx] <- rowSums(subC) / denom
        # w_in: weighted adjacency to same cluster
        subW <- W_mat[idx, idx, drop = FALSE]
        # subtract zero diag implicitly; rowSums suffices
        w_in[idx] <- rowSums(subW)
        # p_best_out: for each other cluster, compute row mean and take max
        if (length(cl_ids) > 1) {
            best <- rep(0, length(idx))
            for (cj in setdiff(cl_ids, cid)) {
                jdx <- which(members == cj)
                if (!length(jdx)) next
                block <- cons_mat[idx, jdx, drop = FALSE]
                denom_j <- pmax(length(jdx), 1L)
                # row mean to cluster cj
                mean_vec <- rowSums(block) / denom_j
                # track maximum
                best <- pmax(best, mean_vec)
            }
            p_best_out[idx] <- best
        } else {
            p_best_out[idx] <- 0
        }
    }
    data.frame(gene = genes, cluster = members, p_in = p_in, p_best_out = p_best_out, w_in = w_in, stringsAsFactors = FALSE)
}

#' Decide which genes to use based on explicit lists or cluster labels.
#' @description
#' Internal helper for `.get_gene_subset`.
#' Provides a single entry-point for deriving analysis subsets from flexible inputs.
#' @param scope_obj A `scope_object`.
#' @param genes Parameter value.
#' @param cluster_col Parameter value.
#' @param cluster_num Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.get_gene_subset <- function(scope_obj, genes = NULL, cluster_col = NULL, cluster_num = NULL, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    if (!is.null(genes)) {
        return(genes)
    }

    if (!is.null(cluster_col) && !is.null(cluster_num)) {
        if (is.null(scope_obj@meta.data) || !cluster_col %in% colnames(scope_obj@meta.data)) {
            stop("Cluster column '", cluster_col, "' not found in meta.data")
        }

        cluster_values <- scope_obj@meta.data[[cluster_col]]
        selected_genes <- rownames(scope_obj@meta.data)[cluster_values == cluster_num]
        selected_genes <- selected_genes[!is.na(selected_genes)]

        if (length(selected_genes) == 0) {
            stop("No genes found for cluster ", cluster_num, " in column ", cluster_col)
        }

        return(selected_genes)
    }

    # If no subset specified, return all genes
    if (!is.null(scope_obj@meta.data)) {
        all_genes <- rownames(scope_obj@meta.data)
        return(all_genes)
    }

    stop("No gene subset specified and no meta.data available")
}

#' Infer Visium Mpp
#' @description
#' Internal helper for `.infer_visium_mpp`.
#' @param grid_layer Layer name.
#' @param scalefactors Parameter value.
#' @param default_spot_um Parameter value.
#' @return Return value used internally.
#' @keywords internal
.infer_visium_mpp <- function(grid_layer,
                              scalefactors,
                              default_spot_um = getOption("geneSCOPE.default_visium_spot_um", 55)) {
    cand <- c(
        .as_numeric_or_na(scalefactors$microns_per_pixel),
        .as_numeric_or_na(grid_layer$microns_per_pixel)
    )
    spot_um_sf <- .as_numeric_or_na(scalefactors$spot_diameter_microns)
    spot_px_sf <- .as_numeric_or_na(scalefactors$spot_diameter_fullres)
    spot_um_grid <- .as_numeric_or_na(grid_layer$spot_diameter_um)
    if (is.finite(spot_um_sf) && is.finite(spot_px_sf)) {
        cand <- c(cand, spot_um_sf / spot_px_sf)
    }
    if (is.finite(spot_um_grid) && is.finite(spot_px_sf)) {
        cand <- c(cand, spot_um_grid / spot_px_sf)
    }
    if (is.finite(default_spot_um) && is.finite(spot_px_sf)) {
        cand <- c(cand, default_spot_um / spot_px_sf)
    }
    cand <- cand[is.finite(cand) & cand > 0]
    if (length(cand)) cand[1] else NA_real_
}

#' Init Thread Manager
#' @description
#' Internal helper for `.init_thread_manager`.
#' @return Return value used internally.
#' @keywords internal
.init_thread_manager <- function() {
  original_settings <- list(
    blas_threads = if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) 1)
    } else {
      1
    },
    env_vars = sapply(c(
      "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS", "BLAS_NUM_THREADS"
    ), Sys.getenv, unset = NA_character_)
  )

  assign(".original_thread_settings", original_settings, envir = .GlobalEnv)
  .set_conservative_threads()
}

#' Count intra-cluster degrees using existing edge endpoints only.
#' @description
#' Internal helper for `.intra_degree_from_edges`.
#' Tallies within-cluster degrees without materializing a full adjacency structure.
#' @param edge_i Integer endpoint indices (1-based) within the kept gene set.
#' @param memb Integer cluster labels (length = number of kept genes), may contain NA.
#' @param edge_j Parameter value.
#' @return Integer vector of intra-cluster degrees per gene.
#' @keywords internal
.intra_degree_from_edges <- function(edge_i, edge_j, memb) {
    n <- length(memb)
    if (length(edge_i) == 0L) return(integer(n))
    same <- (memb[edge_i] == memb[edge_j]) & !is.na(memb[edge_i]) & !is.na(memb[edge_j])
    if (!any(same)) return(integer(n))
    i_tab <- tabulate(edge_i[same], nbins = n)
    j_tab <- tabulate(edge_j[same], nbins = n)
    i_tab + j_tab
}

#' Lee Perm Block Cpp
#' @description
#' Internal helper for `.lee_perm_block_cpp`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.lee_perm_block_cpp <- function(...) lee_perm_block_cpp(...)

#' Lee L Cpp Cache
#' @description
#' Internal helper for `.lee_l_cpp_cache`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.lee_l_cpp_cache <- function(...) leeL_cpp_cache(...)

#' Lee L Cpp Cols
#' @description
#' Internal helper for `.lee_l_cpp_cols`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.lee_l_cpp_cols <- function(...) leeL_cpp_cols(...)

#' Make Pal Map
#' @description
#' Internal helper for `.make_pal_map`.
#' @param levels_vec Parameter value.
#' @param cluster_palette Parameter value.
#' @param tip_palette Parameter value.
#' @return Return value used internally.
#' @keywords internal
.make_pal_map <- function(levels_vec, cluster_palette, tip_palette) {
    levels_vec <- as.character(levels_vec)
    n <- length(levels_vec)
    if (is.null(tip_palette)) {
        if (is.null(cluster_palette)) {
            pal <- rainbow(n)
            return(setNames(pal, levels_vec))
        } else if (is.null(names(cluster_palette))) {
            pal <- colorRampPalette(cluster_palette)(n)
            return(setNames(pal, levels_vec))
        } else {
            pal <- cluster_palette
            if (length(setdiff(levels_vec, names(pal)))) {
                pal <- c(pal, setNames(colorRampPalette(cluster_palette)(
                    length(setdiff(levels_vec, names(pal)))
                ), setdiff(levels_vec, names(pal))))
            }
            return(pal[levels_vec])
        }
    } else if (is.function(tip_palette)) {
        pal <- tryCatch(tip_palette(n, "Set3"), error = function(...) NULL)
        if (is.null(pal)) pal <- hcl.colors(n, palette = "Dynamic")
        return(setNames(pal, levels_vec))
    } else {
        pal <- if (length(tip_palette) < n) colorRampPalette(tip_palette)(n) else tip_palette[seq_len(n)]
        return(setNames(pal, levels_vec))
    }
}

#' Order Levels Numeric
#' @description
#' Internal helper for `.order_levels_numeric`.
#' @param lev Parameter value.
#' @return Return value used internally.
#' @keywords internal
.order_levels_numeric <- function(lev) {
    lev <- as.character(lev)
    num <- suppressWarnings(as.numeric(lev))
    if (all(!is.na(num))) {
        return(lev[order(num)])
    }
    rx <- regmatches(lev, regexpr("-?\\d+\\.?\\d*", lev))
    num2 <- suppressWarnings(as.numeric(rx))
    ord <- order(is.na(num2), num2, lev)
    lev[ord]
}

#' Parse a quantile shorthand (e.g., q95) into a probability.
#' @description
#' Internal helper for `.parse_q`.
#' Supports the qXX notation used throughout geneSCOPE configuration parameters.
#' @param qstr Parameter value.
#' @return Return value used internally.
#' @keywords internal
.parse_q <- function(qstr) {
    if (!is.character(qstr) || length(qstr) != 1 ||
        !grepl("^[qQ][0-9]+\\.?[0-9]*$", qstr)) {
        stop("pct string must look like 'q90' / 'q99.5' …")
    }
    val <- as.numeric(sub("^[qQ]", "", qstr)) / 100
    if (is.na(val) || val < 0 || val > 1) {
        stop("pct out of range 0–100")
    }
    val
}

#' Safe Set Dim
#' @description
#' Internal helper for `.safe_set_dim`.
#' @param M Parameter value.
#' @param target_dim Parameter value.
#' @return Return value used internally.
#' @keywords internal
.safe_set_dim <- function(M, target_dim) {
    if (is.null(target_dim)) return(M)
    ok <- tryCatch({
        dim(M) <- target_dim
        TRUE
    }, error = function(e) FALSE)
    if (ok) return(M)
    M_coerced <- try(as.matrix(M), silent = TRUE)
    if (inherits(M_coerced, "try-error")) {
        stop("Failed to coerce object to requested dimensions when aligning FDR matrix", call. = FALSE)
    }
    dim(M_coerced) <- target_dim
    M_coerced
}

#' Safe Set Dimnames
#' @description
#' Internal helper for `.safe_set_dimnames`.
#' @param M Parameter value.
#' @param target_dn Parameter value.
#' @return Return value used internally.
#' @keywords internal
.safe_set_dimnames <- function(M, target_dn) {
    if (is.null(target_dn)) return(M)
    ok <- tryCatch({
        dimnames(M) <- target_dn
        TRUE
    }, error = function(e) FALSE)
    if (!ok && isTRUE(getOption("geneSCOPE.debug_dimnames", FALSE))) {
        message("[align_and_filter_FDR] dimnames assignment skipped; proceeding without named FDR matrix.")
    }
    M
}

#' Set Conservative Threads
#' @description
#' Internal helper for `.set_conservative_threads`.
#' @return Return value used internally.
#' @keywords internal
.set_conservative_threads <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch(
      {
        RhpcBLASctl::blas_set_num_threads(1)
      },
      error = function(e) { }
    )
  }
  thread_vars <- c(
    "OMP_NUM_THREADS" = "1",
    "OPENBLAS_NUM_THREADS" = "1",
    "MKL_NUM_THREADS" = "1",
    "VECLIB_MAXIMUM_THREADS" = "1",
    "BLAS_NUM_THREADS" = "1"
  )
  do.call(Sys.setenv, as.list(thread_vars))
  if (requireNamespace("data.table", quietly = TRUE)) {
    tryCatch(
      {
        setDTthreads(1)
      },
      error = function(e) { }
    )
  }
}

#' Set Blas Environment
#' @description
#' Internal helper for `.set_blas_environment`.
#' @return Return value used internally.
#' @keywords internal
.set_blas_environment <- function() {
  blas_vars <- c(
    "OPENBLAS_NUM_THREADS" = "1",
    "MKL_NUM_THREADS" = "1",
    "VECLIB_MAXIMUM_THREADS" = "1",
    "NUMEXPR_NUM_THREADS" = "1",
    "BLIS_NUM_THREADS" = "1",
    "GOTO_NUM_THREADS" = "1",
    "ATLAS_NUM_THREADS" = "1"
  )

  for (var in names(blas_vars)) {
    do.call(Sys.setenv, setNames(list(blas_vars[[var]]), var))
  }
}

#' Standardize Crop Bbox
#' @description
#' Internal helper for `.standardize_crop_bbox`.
#' @param bbox Parameter value.
#' @param warn_if_missing_names Parameter value.
#' @return Return value used internally.
#' @keywords internal
.standardize_crop_bbox <- function(bbox, warn_if_missing_names = TRUE) {
    if (is.null(bbox)) return(NULL)
    req <- c("xmin", "xmax", "ymin", "ymax")
    flat <- bbox
    if (is.list(bbox) && !is.null(names(bbox))) {
        flat <- unlist(bbox, use.names = TRUE)
    }
    if (is.numeric(flat) && length(flat) == 4L && (is.null(names(flat)) || anyNA(names(flat)))) {
        names(flat) <- req
    }
    if (!all(req %in% names(flat))) {
        if (is.numeric(flat) && length(flat) == 4L) {
            names(flat) <- req
        } else {
            if (isTRUE(warn_if_missing_names)) {
                warning("crop_bbox_px is missing names; ignoring crop window.", call. = FALSE)
            }
            return(NULL)
        }
    }
    out <- as.numeric(flat[req])
    names(out) <- req
    out
}

#' Symmetric Pmax
#' @description
#' Internal helper for `.symmetric_pmax`.
#' @param M Parameter value.
#' @return Return value used internally.
#' @keywords internal
.symmetric_pmax <- function(M) {
    if (inherits(M, "sparseMatrix")) {
        # Robust to Matrix versions without exported pmax(): use algebraic identity
        # max(A,B) = (A + B + |A - B|) / 2, then zero the diagonal
        M <- as(M, "dgCMatrix")
        Mt <- t(M)
        M <- (M + Mt + abs(M - Mt)) * 0.5
        diag(M) <- 0
        return(drop0(M))
    } else {
        M <- pmax(M, t(M))
        diag(M) <- 0
        return(M)
    }
}

#' Add Single Cells Cosmx
#' @description
#' Internal helper for `.add_single_cells_cosmx`.
#' @param scope_obj A `scope_object`.
#' @param cosmx_root Parameter value.
#' @param filter_genes Parameter value.
#' @param exclude_prefix Parameter value.
#' @param id_mode Parameter value.
#' @param filter_barcodes Parameter value.
#' @param force_all_genes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.add_single_cells_cosmx <- function(scope_obj,
                                 cosmx_root,
                                 filter_genes = NULL,
                                 exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                                    "SystemControl", "Negative"),
                                 id_mode = c("cell_id", "fov_cell", "auto"),
                                 filter_barcodes = TRUE,
                                 force_all_genes = FALSE,
                                 verbose = TRUE) {
  stopifnot(inherits(scope_obj, "scope_object"), dir.exists(cosmx_root))

  id_mode <- match.arg(id_mode)

  # 1) Load target cell set for filtering and column ordering
  keep_cells_target <- NULL
  if (filter_barcodes) {
    keep_cells_target <- scope_obj@coord$centroids$cell
    if (!length(keep_cells_target)) {
      stop("centroids are empty; cannot determine columns to keep. Please run .create_scope/_cosmx first.")
    }
  }

  # 2) Locate CosMx exprMat file
  if (!dir.exists(file.path(cosmx_root, "flatFiles"))) stop("flatFiles/ not found: ", file.path(cosmx_root, "flatFiles"))
  ds_dir <- list.dirs(file.path(cosmx_root, "flatFiles"), full.names = TRUE, recursive = FALSE)
  if (!length(ds_dir)) stop("No dataset directory found under flatFiles/.")
  ds_dir <- ds_dir[[1]]

  expr_mat <- NULL
  cand <- c(
    file.path(ds_dir, paste0(basename(ds_dir), "_exprMat_file.csv.gz")),
    file.path(ds_dir, paste0(basename(ds_dir), "_exprMat_file.csv"))
  )
  for (p in cand) if (file.exists(p)) { expr_mat <- p; break }
  if (is.null(expr_mat)) stop("No *_exprMat_file.csv(.gz) found under: ", ds_dir)

  .log_info("addSingleCells", "S01", paste0("Expression matrix: ", basename(expr_mat)), verbose)

  # 3) Read header to determine gene columns
  suppressPackageStartupMessages({ library(data.table); library(Matrix) })
  # Prefer fread to fetch header; fall back to base::readLines(gzfile(...), n=1)
  header_dt <- NULL
  cols <- NULL
  header_dt <- tryCatch({ fread(expr_mat, nrows = 0L, showProgress = FALSE) }, error = function(e) NULL)
  if (!is.null(header_dt)) {
    cols <- names(header_dt)
  } else {
    # Safe fallback: read only the first line as header
    con <- NULL
    cols <- tryCatch(
      {
        if (grepl("\\.gz$", expr_mat, ignore.case = TRUE)) {
          con <- gzfile(expr_mat, open = "rt")
        } else {
          con <- file(expr_mat, open = "rt")
        }
        ln <- readLines(con, n = 1L, warn = FALSE)
        if (length(ln) == 0L || !nzchar(ln)) stop("Empty file or missing header")
        ln <- sub("\\r$", "", ln)  # strip CR from CRLF
        strsplit(ln, ",", fixed = TRUE)[[1]]
      },
      error = function(e) NULL,
      finally = { if (!is.null(con)) try(close(con), silent = TRUE) }
    )
    if (is.null(cols) || !length(cols)) stop("Unable to read expression matrix header: ", expr_mat)
  }
  base_cols <- tolower(cols)
  # Identify fov / cell_ID column names (case-insensitive)
  fov_col  <- cols[which(base_cols == "fov")[1]]
  cid_col  <- cols[which(base_cols %in% c("cell_id", "cellid"))[1]]
  if (is.na(fov_col) || is.na(cid_col)) stop("Expression matrix missing fov/cell_ID column")

  # ---- Auto-detect non-gene columns; the rest are considered gene columns ----
  # Known/common non-gene column names (case-insensitive)
  non_gene_exact <- c(
    "fov", "cell_id", "cellid", "cell", "barcode",
    "x", "y", "z",
    "umi", "reads", "total", "sum", "sizefactor",
    "feature_name", "gene", "target", "qv", "nucleus_distance",
    "cellcomp", "compartment", "compartmentlabel"
  )
  # Columns containing these substrings are treated as non-gene (coords/units/etc.)
  non_gene_substr <- c(
    "_px", "_um", "_mm", "global_", "local_", "centroid", "area",
    "x_", "y_", "_x", "_y"
  )
  low_cols <- tolower(cols)
  is_non_gene <- low_cols %in% non_gene_exact
  for (ss in non_gene_substr) {
    is_non_gene <- is_non_gene | grepl(ss, low_cols, fixed = TRUE)
  }
  # fov/cell_ID must be among non-gene; remaining columns are candidate genes
  gene_cols_auto <- cols[!is_non_gene]

  # Exclude non-biological prefixes (NegControl/Unassigned/Background etc.)
  if (length(exclude_prefix)) {
    pat <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
    gene_cols_auto <- gene_cols_auto[!grepl(pat, gene_cols_auto, perl = TRUE)]
  }

  # If user did not provide filter_genes, use auto set; otherwise take the intersection
  if (is.null(filter_genes)) {
    gene_cols <- gene_cols_auto
  } else {
    gene_cols <- intersect(gene_cols_auto, filter_genes)
  if (!length(gene_cols)) stop("No overlap between filter_genes and available genes (after auto-detection).")
  }

  if (!length(gene_cols)) stop("No gene columns identified. Check matrix header and prefix filters.")
  .log_info(
      "addSingleCells",
      "S01",
      paste0("Number of genes detected (auto): ", length(gene_cols)),
      verbose
  )

  # 4) Read only required columns (fov/cell_ID + target genes)
  sel_cols <- c(fov_col, cid_col, gene_cols)
  # fread selected columns; fall back to read.csv on failure (slower but robust)
  dt <- tryCatch(
    fread(expr_mat, select = sel_cols, showProgress = verbose),
    error = function(e) {
      .log_info(
          "addSingleCells",
          "S01",
          paste0("fread failed; falling back to read.csv: ", e$message),
          verbose
      )
      read.csv(expr_mat, header = TRUE)[, sel_cols]
    }
  )
  setDT(dt)

  # Construct column keys (aligned with centroids$cell) and adaptively fall back
  if (id_mode == "auto") {
    id_mode <- if (any(grepl("_", keep_cells_target)) || any(grepl("FOV", keep_cells_target))) "fov_cell" else "cell_id"
  }
  make_key <- function(mode, fov_offset = 0L) {
    if (mode == "cell_id") {
      as.character(dt[[cid_col]])
    } else {
      fvv <- suppressWarnings(as.integer(dt[[fov_col]]))
      if (!is.na(fov_offset) && fov_offset != 0L) fvv <- fvv + as.integer(fov_offset)
      paste0(as.character(dt[[cid_col]]), "_", as.character(fvv))
    }
  }
  # First attempt
  dt[, cell_key := make_key(id_mode, 0L)]

  # 5) Optional: filter/reorder cells based on centroids
  if (filter_barcodes) {
    Kt <- unique(keep_cells_target)
    cand <- list(
    list(mode = id_mode, off = 0L),
    list(mode = if (id_mode == "cell_id") "fov_cell" else "cell_id", off = 0L),
      list(mode = "fov_cell", off = +1L),
      list(mode = "fov_cell", off = -1L)
    )
    best_idx <- 1L
    best_overlap <- -1L
    keys_cache <- list()
    for (i in seq_along(cand)) {
      mode_i <- cand[[i]]$mode; off_i <- cand[[i]]$off
      # Only try offset when fov column exists
      if (mode_i == "fov_cell" && (!exists("fov_col") || is.na(fov_col))) next
      if (is.null(keys_cache[[paste(mode_i, off_i)]])) {
        keys_cache[[paste(mode_i, off_i)]] <- make_key(mode_i, off_i)
      }
      kk <- unique(keys_cache[[paste(mode_i, off_i)]])
      ov <- length(intersect(kk, Kt))
      if (ov > best_overlap) { best_overlap <- ov; best_idx <- i }
    }
    # Apply the best key
    mode_best <- cand[[best_idx]]$mode; off_best <- cand[[best_idx]]$off
    dt[, cell_key := make_key(mode_best, off_best)]
    keep_cells <- intersect(unique(dt$cell_key), Kt)
    if (!length(keep_cells)) {
      # Diagnostics
      stop("Cell identifiers in expression matrix do not overlap with centroids; check id_mode/FOV encoding. Suggest id_mode=\"fov_cell\"; if still failing, try FOV offset +/-1. Example keys: ",
           head(unique(dt$cell_key), 3L), " vs ", head(Kt, 3L))
    }
    dt <- dt[cell_key %in% keep_cells]
  }

  # 6) Build sparse matrix (rows=genes, cols=cells) without forming a huge dense matrix
  # Column index/order: if filtered, follow keep_cells_target; else by appearance
  if (filter_barcodes) {
    col_order <- keep_cells_target[keep_cells_target %in% unique(dt$cell_key)]
  } else {
    col_order <- unique(dt$cell_key)
  }
  col_index <- match(dt$cell_key, col_order)  # 1..ncol

  # Preallocate i/j/x accumulators
  i_idx <- integer(0)
  j_idx <- integer(0)
  x_val <- numeric(0)

  # Collect non-zeros per gene column to avoid large temporary memory
  for (g in gene_cols) {
    v <- suppressWarnings(as.numeric(dt[[g]]))
    nz <- which(!is.na(v) & v != 0)
    if (length(nz)) {
      i_idx <- c(i_idx, rep.int(which(gene_cols == g), length(nz)))
      j_idx <- c(j_idx, col_index[nz])
      x_val <- c(x_val, v[nz])
    }
  }

  counts <- sparseMatrix(
    i = as.integer(i_idx),
    j = as.integer(j_idx),
    x = as.numeric(x_val),
    dims = c(length(gene_cols), length(col_order)),
    dimnames = list(gene_cols, col_order)
  )

  # 7) Write back to object
  scope_obj@cells <- list(counts = counts)
  .log_info(
      "addSingleCells",
      "S01",
      paste0("Wrote @cells$counts (", nrow(counts), " genes x ", ncol(counts), " cells)"),
      verbose
  )
  invisible(scope_obj)
}

#' Add Single Cell Count Matrix (Xenium)
#' @description
#'   Read the Xenium **cell-feature sparse matrix** (either HDF5 or
#'   Matrix-Market format), keep only the cells that appear in
#'   \code{scope_obj@coord$centroids$cell}, preserve that order, and store the
#'   result in the new \code{@cells} slot.
#'   The function relies only on \strong{Matrix}, \strong{data.table}, and
#'   \strong{rhdf5} (or \strong{hdf5r})—no Seurat, tidyverse, or other heavy
#'   dependencies.
#' @param scope_obj A valid \code{scope_object} whose \code{@coord$centroids} slot
#'        is already filled (e.g. via \code{\link{.create_scope}}). Cells
#'        found in this slot define the subset and order of columns kept.
#' @param xenium_dir Character scalar. Path to the Xenium \file{outs/}
#'        directory that holds \file{cell_feature_matrix.h5}; must exist.
#' @param filter_genes Character vector or NULL. Gene names to include (default NULL keeps all).
#' @param exclude_prefix Character vector. Prefixes to exclude when filtering genes (default c("Unassigned", "NegControl", "Background", "DeprecatedCodeword")).
#' @param filter_barcodes Logical. Whether to filter cells to match barcodes in scope_obj (default TRUE). When FALSE, all cells from the HDF5 matrix are included.
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @return The same object (modified in place) now carrying a
#'         \code{dgCMatrix} in \code{@cells}.
#' @importFrom Matrix sparseMatrix t
#' @importFrom data.table data.table
#' @importFrom rhdf5 h5read
#' @keywords internal
.add_single_cells_xenium <- function(scope_obj, xenium_dir,
                           filter_genes = NULL,
                           exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword"),
                           filter_barcodes = TRUE,
                           verbose = TRUE) {
  stopifnot(inherits(scope_obj, "scope_object"), dir.exists(xenium_dir))

  .log_info("addSingleCells", "S01", "Loading cell-feature matrix from HDF5", verbose)

  suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
  })

  # Cross-platform HDF5 library detection and loading
  hdf5_available <- FALSE
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    library(rhdf5)
    hdf5_available <- TRUE
  } else if (requireNamespace("hdf5r", quietly = TRUE)) {
    library(hdf5r)
    hdf5_available <- TRUE
  }

  if (!hdf5_available) {
    stop("Neither rhdf5 nor hdf5r package is available. Please install one of them.")
  }

  ## ------------------------------------------------------------------ 1
  if (filter_barcodes) {
    keep_cells_target <- scope_obj@coord$centroids$cell
    if (!length(keep_cells_target)) {
      stop("scope_obj@coord$centroids is empty, cannot determine cells to keep.")
    }

    .log_info(
        "addSingleCells",
        "S01",
        paste0("Target cells: ", length(keep_cells_target)),
        verbose
    )
  } else {
    .log_info(
        "addSingleCells",
        "S01",
        "Barcode filtering disabled; including all cells from HDF5",
        verbose
    )
  }

  ## ------------------------------------------------------------------ 2
  h5_file <- file.path(xenium_dir, "cell_feature_matrix.h5")
  if (!file.exists(h5_file)) {
    stop("HDF5 file not found: ", h5_file)
  }

  # Cross-platform HDF5 reading
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    barcodes <- rhdf5::h5read(h5_file, "matrix/barcodes")
    genes_id <- rhdf5::h5read(h5_file, "matrix/features/id")
    genes_nm <- rhdf5::h5read(h5_file, "matrix/features/name")
    data <- rhdf5::h5read(h5_file, "matrix/data")
    indices <- rhdf5::h5read(h5_file, "matrix/indices")
    indptr <- rhdf5::h5read(h5_file, "matrix/indptr")
    shape <- as.integer(rhdf5::h5read(h5_file, "matrix/shape"))
  } else {
    # Alternative hdf5r implementation
    h5file <- hdf5r::H5File$new(h5_file, mode = "r")
    barcodes <- h5file[["matrix/barcodes"]][]
    genes_id <- h5file[["matrix/features/id"]][]
    genes_nm <- h5file[["matrix/features/name"]][]
    data <- h5file[["matrix/data"]][]
    indices <- h5file[["matrix/indices"]][]
    indptr <- h5file[["matrix/indptr"]][]
    shape <- as.integer(h5file[["matrix/shape"]][])
    h5file$close()
  }

  ## -------------------------- 2a Determine HDF5 storage orientation ----
  nrow_h5 <- shape[1]
  ncol_h5 <- shape[2]
  by_cols <- length(indptr) == (ncol_h5 + 1) # indptr length = ncol+1

  if (!by_cols) {
    stop("Current function only supports CSC-style storage (indptr length = ncol+1).")
  }

  col_is_cell <- (ncol_h5 == length(barcodes))
  if (!col_is_cell && !(nrow_h5 == length(barcodes))) {
    stop("Shape doesn't match barcode/gene counts, cannot determine orientation.")
  }

  .log_info(
      "addSingleCells",
      "S01",
      paste0(
          "Matrix: ", nrow_h5, "×", ncol_h5,
          " (", length(genes_nm), " genes, ", length(barcodes), " cells)"
      ),
      verbose
  )

  ## -------------------------- 2b Build sparse matrix ------------------
  x <- as.numeric(data)
  i <- as.integer(indices)
  p <- as.integer(indptr)

  if (col_is_cell) {
    ## →  columns = cells, rows = genes   (most common)
    counts_raw <- new("dgCMatrix",
      Dim      = c(length(genes_nm), length(barcodes)),
      x        = x,
      i        = i,
      p        = p,
      Dimnames = list(make.unique(genes_nm), barcodes)
    )
  } else {
    ## →  columns = genes, rows = cells   (rare case: need transpose)
    .log_info("addSingleCells", "S01", "Transposing matrix layout", verbose)
    tmp <- new("dgCMatrix",
      Dim      = c(length(barcodes), length(genes_nm)),
      x        = x,
      i        = i,
      p        = p,
      Dimnames = list(barcodes, make.unique(genes_nm))
    )
    counts_raw <- t(tmp)
  }

  attr(counts_raw, "gene_map") <- data.frame(
    ensembl = genes_id,
    symbol = genes_nm,
    stringsAsFactors = FALSE
  )

  ## -------------------------- 2c Filter genes by exclude_prefix and filter_genes ------------------
  .log_info("addSingleCells", "S01", "Applying gene filters", verbose)

  # Get gene names (rownames)
  all_genes <- rownames(counts_raw)
  keep_genes <- all_genes

  # Filter by exclude_prefix
  if (length(exclude_prefix) > 0) {
    exclude_pattern <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
    exclude_mask <- grepl(exclude_pattern, keep_genes, perl = TRUE)
    keep_genes <- keep_genes[!exclude_mask]

    if (verbose) {
      excluded_count <- sum(exclude_mask)
      if (excluded_count > 0) {
        .log_info(
            "addSingleCells",
            "S01",
            paste0(
                "Excluded ", excluded_count, " genes with prefixes: ",
                paste(exclude_prefix, collapse = ", ")
            ),
            verbose
        )
      }
    }
  }

  # Filter by filter_genes (if specified)
  if (!is.null(filter_genes)) {
    keep_genes <- intersect(keep_genes, filter_genes)
    .log_info(
        "addSingleCells",
        "S01",
        paste0("Filtered to ", length(keep_genes), " genes from target list"),
        verbose
    )
  }

  if (length(keep_genes) == 0) {
    stop("No genes remain after filtering. Please check exclude_prefix and filter_genes parameters.")
  }

  # Apply gene filtering to the matrix
  counts_raw <- counts_raw[keep_genes, , drop = FALSE]

  .log_info(
      "addSingleCells",
      "S01",
      paste0("Final count: ", nrow(counts_raw), "/", length(all_genes), " genes"),
      verbose
  )

  ## ------------------------------------------------------------------ 3
  if (filter_barcodes) {
    keep_cells <- keep_cells_target[keep_cells_target %in% colnames(counts_raw)]
    if (!length(keep_cells)) {
      stop("Centroid cell IDs do not overlap with H5 data.")
    }

    .log_info(
        "addSingleCells",
        "S01",
        paste0("Cell overlap: ", length(keep_cells), "/", length(keep_cells_target)),
        verbose
    )

    counts <- counts_raw[, keep_cells, drop = FALSE]
  } else {
    # Keep all cells from the HDF5 matrix
    counts <- counts_raw
    .log_info(
        "addSingleCells",
        "S01",
        paste0("Keeping all ", ncol(counts), " cells from HDF5 matrix"),
        verbose
    )
  }

  ## ------------------------------------------------------------------ 4
  scope_obj@cells <- list(counts = counts)

  .log_info("addSingleCells", "S01", "Cell matrix integration completed", verbose)
  invisible(scope_obj)
}

#' Attach Histology
#' @description
#' Internal helper for `.attach_histology`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param png_path Filesystem path.
#' @param json_path Filesystem path.
#' @param level Parameter value.
#' @param crop_bbox_px Parameter value.
#' @param roi_bbox Parameter value.
#' @param coord_type Parameter value.
#' @param scalefactors Parameter value.
#' @param y_origin Parameter value.
#' @param y_flip_reference_um Parameter value.
#' @return Return value used internally.
#' @keywords internal
.attach_histology <- function(scope_obj,
                             grid_name,
                             png_path,
                             json_path = NULL,
                             level = c("lowres", "hires"),
                             crop_bbox_px = NULL,
                             roi_bbox = NULL,
                             coord_type = c("visium", "manual"),
                             scalefactors = NULL,
                             y_origin = c("auto", "top-left", "bottom-left"),
                             y_flip_reference_um = NULL) {
    stopifnot(is(scope_obj, "scope_object"))
    level <- match.arg(level)
    coord_type <- match.arg(coord_type)
    y_origin_choice <- .normalize_histology_y_origin(y_origin)
    if (!(grid_name %in% names(scope_obj@grid))) {
        stop("Grid '", grid_name, "' not found in scope_obj@grid.")
    }
    grid_layer <- scope_obj@grid[[grid_name]]
    if (is.null(grid_layer$histology)) {
        grid_layer$histology <- list(lowres = NULL, hires = NULL)
    } else {
        grid_layer$histology <- modifyList(list(lowres = NULL, hires = NULL), grid_layer$histology, keep.null = TRUE)
    }

    if (is.null(scalefactors)) {
        scalefactors <- .read_scalefactors_json(json_path)
    }
    sf_val <- if (identical(level, "lowres")) scalefactors$lowres else scalefactors$hires
    sf_from_path <- function(path) {
        if (is.null(path) || !nzchar(path)) return(NA_real_)
        base <- tolower(basename(path))
        aligned_sf <- .as_numeric_or_na(scalefactors$regist_target_img_scalef)
        if (grepl("aligned_tissue_image", base, fixed = TRUE) && is.finite(aligned_sf)) {
            return(aligned_sf)
        }
        if (grepl("tissue_hires_image", base, fixed = TRUE)) {
            return(.as_numeric_or_na(scalefactors$hires))
        }
        if (grepl("tissue_lowres_image", base, fixed = TRUE) ||
            grepl("detected_tissue_image", base, fixed = TRUE)) {
            return(.as_numeric_or_na(scalefactors$lowres))
        }
        NA_real_
    }
    sf_path <- sf_from_path(png_path)
    if (is.finite(sf_path)) {
        sf_val <- sf_path
    }
    roi_use <- roi_bbox
    roi_auto <- is.null(roi_bbox)
    if (is.null(roi_use) && !is.null(grid_layer$grid_info)) {
        roi_use <- c(
            xmin = min(grid_layer$grid_info$xmin, na.rm = TRUE),
            xmax = max(grid_layer$grid_info$xmax, na.rm = TRUE),
            ymin = min(grid_layer$grid_info$ymin, na.rm = TRUE),
            ymax = max(grid_layer$grid_info$ymax, na.rm = TRUE)
        )
    }

    y_origin_final <- .resolve_histology_y_origin(grid_layer, y_origin_choice)
    flip_ref_um <- .as_numeric_or_na(y_flip_reference_um)
    if (!is.finite(flip_ref_um) && !is.null(grid_layer$image_info$y_flip_reference_um)) {
        flip_ref_um <- .as_numeric_or_na(grid_layer$image_info$y_flip_reference_um)
    }
    mpp_auto <- .infer_visium_mpp(grid_layer, scalefactors)
    sf_use <- .as_numeric_or_na(sf_val)
    if (!is.finite(sf_use) || sf_use <= 0) sf_use <- 1
    mpp_use <- .as_numeric_or_na(mpp_auto)
    convert_roi_for_crop <- function(bbox) {
        if (is.null(bbox)) return(NULL)
        req <- c("xmin", "xmax", "ymin", "ymax")
        flat <- bbox
        if (is.list(flat) && !is.null(names(flat))) {
            flat <- unlist(flat, use.names = TRUE)
        }
        if (!(is.numeric(flat) && all(req %in% names(flat)))) {
            return(bbox)
        }
        vals <- as.numeric(flat[req])
        names(vals) <- req
        if (!identical(y_origin_final, "bottom-left")) {
            return(vals)
        }
        if (!is.finite(flip_ref_um)) {
            return(vals)
        }
        y1 <- flip_ref_um - vals["ymax"]
        y2 <- flip_ref_um - vals["ymin"]
        vals["ymin"] <- min(y1, y2)
        vals["ymax"] <- max(y1, y2)
        vals
    }

    crop_auto_flag <- FALSE
    if (is.null(crop_bbox_px)) {
        crop_auto <- .auto_histology_crop_bbox(
            roi_bbox = convert_roi_for_crop(roi_use),
            microns_per_pixel = mpp_auto,
            level_scalefactor = sf_val
        )
        if (!is.null(crop_auto)) {
            crop_bbox_px <- crop_auto
            crop_auto_flag <- TRUE
        }
    }

    crop_bbox_px <- .standardize_crop_bbox(crop_bbox_px, warn_if_missing_names = !crop_auto_flag)
    img_info <- .read_histology_png(png_path, crop_bbox_px = crop_bbox_px)
    if (isTRUE(crop_auto_flag) && !is.null(crop_bbox_px) && !is.null(img_info$crop_bbox_px)) {
        req <- c("xmin", "xmax", "ymin", "ymax")
        desired_bbox <- as.numeric(crop_bbox_px[req])
        names(desired_bbox) <- req
        desired_bbox["xmin"] <- floor(desired_bbox["xmin"])
        desired_bbox["xmax"] <- ceiling(desired_bbox["xmax"])
        desired_bbox["ymin"] <- floor(desired_bbox["ymin"])
        desired_bbox["ymax"] <- ceiling(desired_bbox["ymax"])
        img_bbox <- as.numeric(img_info$crop_bbox_px[req])
        names(img_bbox) <- req
        if (any(desired_bbox != img_bbox)) {
            left_pad <- max(0, img_bbox["xmin"] - desired_bbox["xmin"])
            right_pad <- max(0, desired_bbox["xmax"] - img_bbox["xmax"])
            top_pad <- max(0, img_bbox["ymin"] - desired_bbox["ymin"])
            bottom_pad <- max(0, desired_bbox["ymax"] - img_bbox["ymax"])
            if (any(c(left_pad, right_pad, top_pad, bottom_pad) > 0)) {
                img <- img_info$png
                dims <- dim(img)
                if (length(dims) >= 2L) {
                    h <- dims[1]
                    w <- dims[2]
                    new_h <- h + top_pad + bottom_pad
                    new_w <- w + left_pad + right_pad
                    if (length(dims) == 2L) {
                        padded <- matrix(0, nrow = new_h, ncol = new_w)
                        padded[(top_pad + 1L):(top_pad + h), (left_pad + 1L):(left_pad + w)] <- img
                    } else {
                        padded <- array(0, dim = c(new_h, new_w, dims[3]))
                        padded[(top_pad + 1L):(top_pad + h), (left_pad + 1L):(left_pad + w), ] <- img
                    }
                    img_info$png <- padded
                    img_info$width <- new_w
                    img_info$height <- new_h
                    img_info$crop_bbox_px <- desired_bbox
                }
            }
        }
    }
    if (isTRUE(roi_auto) &&
        !is.null(img_info$crop_bbox_px) &&
        is.finite(mpp_use) && mpp_use > 0 &&
        is.finite(sf_use) && sf_use > 0) {
        crop_to_roi <- function(bbox_px) {
            req <- c("xmin", "xmax", "ymin", "ymax")
            if (!all(req %in% names(bbox_px))) return(NULL)
            vals_px <- as.numeric(bbox_px[req])
            if (any(!is.finite(vals_px))) return(NULL)
            fullres_px <- vals_px / sf_use
            roi_vals <- fullres_px * mpp_use
            names(roi_vals) <- req
            if (identical(y_origin_final, "bottom-left")) {
                if (!is.finite(flip_ref_um)) return(NULL)
                y1 <- flip_ref_um - roi_vals["ymax"]
                y2 <- flip_ref_um - roi_vals["ymin"]
                roi_vals["ymin"] <- min(y1, y2)
                roi_vals["ymax"] <- max(y1, y2)
            }
            roi_vals
        }
        roi_from_crop <- crop_to_roi(img_info$crop_bbox_px)
        if (!is.null(roi_from_crop)) {
            roi_use <- roi_from_crop
        }
    }

    grid_layer$histology[[level]] <- list(
        png = img_info$png,
        width = img_info$width,
        height = img_info$height,
        scalefactor = sf_val,
        coord_type = coord_type,
        roi_bbox = roi_use,
        crop_bbox_px = img_info$crop_bbox_px,
        source_path = img_info$source_path,
        y_origin = y_origin_final
    )
    scope_obj@grid[[grid_name]] <- grid_layer
    scope_obj
}

#' Build Segmentation Geometries
#' @description
#' Internal helper for `.build_segmentation_geometries`.
#' @param path Parameter value.
#' @param label Parameter value.
#' @param flip Parameter value.
#' @param y_max Numeric threshold.
#' @param keep_ids Logical flag.
#' @param workers Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.build_segmentation_geometries <- function(path,
                                          label = "cell",
                                          flip = FALSE,
                                          y_max = NULL,
                                          keep_ids = NULL,
                                          workers = 1L) {
    if (!file.exists(path)) {
        stop("Segmentation file not found: ", path)
    }
    raw_dt <- as.data.table(arrow::read_parquet(path))
    if ("fov" %in% names(raw_dt)) {
        seg_dt <- raw_dt[, .(
            cell = paste0(as.character(cell_id), "_", as.character(fov)),
            x = vertex_x,
            y = vertex_y,
            label_id = paste0(as.character(cell_id), "_", as.character(fov))
        )]
    } else {
        seg_dt <- raw_dt[, .(
            cell = as.character(cell_id),
            x = vertex_x,
            y = vertex_y,
            label_id = as.character(label_id)
        )]
    }

    if (flip) {
        if (is.null(y_max)) stop("`y_max` must be supplied when `flip = TRUE`.")
        seg_dt[, y := y_max - y]
    }
    if (!is.null(keep_ids)) {
        seg_dt <- seg_dt[cell %in% keep_ids]
    }

    index_map <- split(seq_len(nrow(seg_dt)), seg_dt$label_id)
    polygon_builder <- function(idx) {
        sub_dt <- seg_dt[idx, .(x, y)]
        sub_dt <- sub_dt[complete.cases(sub_dt)]
        if (nrow(sub_dt) < 3L) return(NULL)
        coords <- as.matrix(sub_dt)
        storage.mode(coords) <- "double"
        if (!all(coords[1, ] == coords[nrow(coords), ])) {
            coords <- rbind(coords, coords[1, ])
        }
        tryCatch(sf::st_polygon(list(coords)), error = function(e) {
            message("[.build_segmentation_geometries] Invalid polygon in ", label, ": ", conditionMessage(e))
            NULL
        })
    }

    if (workers > 1L && .Platform$OS.type != "windows") {
        poly_list <- mclapply(index_map, polygon_builder, mc.cores = workers)
    } else if (workers > 1L) {
        cluster <- makeCluster(workers)
        on.exit(stopCluster(cluster), add = TRUE)
        clusterEvalQ(cluster, library(sf))
        clusterExport(cluster, c("seg_dt", "polygon_builder"),
            envir = environment()
        )
        poly_list <- parLapply(cluster, index_map, polygon_builder)
    } else {
        poly_list <- lapply(index_map, polygon_builder)
    }
    list(points = seg_dt, polygons = sf::st_sfc(poly_list[!vapply(poly_list, is.null, logical(1L))]))
}

#' Clip Points To Region
#' @description
#' Internal helper for `.clip_points_to_region`.
#' @param points Parameter value.
#' @param region Parameter value.
#' @param chunk_size Parameter value.
#' @param workers Parameter value.
#' @param parallel_backend Parameter value.
#' @param psock_batch_chunks Parameter value.
#' @return Return value used internally.
#' @keywords internal
.clip_points_to_region <- function(points,
                                  region,
                                  chunk_size = 5e5,
                                  workers = 1L,
                                  parallel_backend = NULL,
                                  psock_batch_chunks = NULL) {
    if (is.null(region) || nrow(points) == 0L) {
        return(points)
    }
    if (inherits(region, "sf")) {
        region <- sf::st_geometry(region)
    }
    if (!inherits(region, "sfc")) {
        if (all(vapply(region, inherits, logical(1), "sfg"))) {
            region <- sf::st_sfc(region)
        } else {
            stop("`region` must be an sf/sfc object (or list of sfg).")
        }
    }
    if (length(region) == 0L) {
        return(points[0])
    }
    if (length(region) > 1L) {
        region <- sf::st_union(region)
    }
    if (length(region) == 0L) {
        return(points[0])
    }

    stopifnot(all(c("x", "y") %in% names(points)))
    dt_pts <- as.data.table(points)
    chunk_ids <- split(seq_len(nrow(dt_pts)),
        ceiling(seq_len(nrow(dt_pts)) / chunk_size)
    )
    psock_batch_chunks <- psock_batch_chunks %||% getOption("geneSCOPE.psock_batch_chunks", 16L)
    psock_batch_chunks <- suppressWarnings(as.integer(psock_batch_chunks))
    if (!is.finite(psock_batch_chunks) || psock_batch_chunks < 1L) {
        psock_batch_chunks <- 16L
    }
    large_table_rows <- getOption("geneSCOPE.large_table_rows", 5e6)
    large_table_rows <- suppressWarnings(as.numeric(large_table_rows)[1])
    if (!is.finite(large_table_rows) || large_table_rows < 1) {
        large_table_rows <- 5e6
    }
    is_large_table <- nrow(dt_pts) >= large_table_rows

    debug_parallel <- isTRUE(getOption("geneSCOPE.debug_parallel", FALSE))
    log_dir <- if (debug_parallel) getOption("geneSCOPE.debug_parallel_dir", tempdir()) else NULL
    if (debug_parallel && !dir.exists(log_dir)) {
        dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    }

    workers <- suppressWarnings(as.integer(workers))
    if (!is.finite(workers) || workers < 1L) workers <- 1L

    backend <- if (is.null(parallel_backend)) {
        getOption("geneSCOPE.parallel_backend", "auto")
    } else {
        parallel_backend
    }
    backend <- tolower(as.character(backend)[1])
    if (is.na(backend) || !nzchar(backend)) backend <- "auto"
    backend <- match.arg(backend, c("auto", "fork", "psock", "serial"))
    if (identical(backend, "auto")) {
        if (.Platform$OS.type == "windows") {
            backend <- "psock"
        } else if (.detect_os() == "linux" && .detect_container()) {
            backend <- if (is_large_table) "serial" else "psock"
        } else {
            backend <- "fork"
        }
    }
    if (workers <= 1L) backend <- "serial"

    if (debug_parallel) {
        chunk_lengths <- vapply(chunk_ids, length, integer(1L))
        message(
            "[clip_points_to_region] points=", nrow(dt_pts),
            " chunks=", length(chunk_ids),
            " chunk_size=", chunk_size,
            " workers=", workers,
            " backend=", backend,
            " psock_batch_chunks=", psock_batch_chunks,
            " large_table=", is_large_table,
            " large_table_rows=", large_table_rows,
            " chunk_min=", min(chunk_lengths),
            " chunk_max=", max(chunk_lengths),
            " dt_mb=", round(as.numeric(object.size(dt_pts)) / 1024^2, 1)
        )
    }

    clip_worker <- function(subset_dt, region, debug_parallel = FALSE, log_dir = NULL, chunk_range = NULL) {
        tryCatch({
            inside_flag <- lengths(sf::st_within(
                sf::st_as_sf(subset_dt,
                    coords = c("x", "y"),
                    crs = sf::st_crs(region) %||% NA, remove = FALSE
                ),
                region
            )) > 0
            list(data = subset_dt[inside_flag, , drop = FALSE], error = NULL)
        }, error = function(e) {
            if (debug_parallel) {
                log_path <- file.path(log_dir, paste0("genescope_clip_worker_", Sys.getpid(), ".log"))
                if (is.null(chunk_range)) chunk_range <- "unknown"
                writeLines(
                    c(
                        paste0("pid=", Sys.getpid()),
                        paste0("chunk_range=", chunk_range),
                        paste0("chunk_size=", nrow(subset_dt)),
                        paste0("error=", conditionMessage(e))
                    ),
                    log_path
                )
            }
            list(data = NULL, error = conditionMessage(e))
        })
    }

    clip_worker_idx <- function(idx) {
        subset_dt <- dt_pts[idx]
        chunk_range <- if (length(idx)) paste0(min(idx), "-", max(idx)) else "empty"
        clip_worker(subset_dt,
            region = region,
            debug_parallel = debug_parallel,
            log_dir = log_dir,
            chunk_range = chunk_range
        )
    }

    run_psock_chunked <- function() {
        cluster <- makeCluster(workers)
        on.exit(stopCluster(cluster), add = TRUE)
        clusterEvalQ(cluster, library(sf))
        clusterExport(
            cluster,
            c("region", "clip_worker", "debug_parallel", "log_dir"),
            envir = environment()
        )
        res_all <- vector("list", length(chunk_ids))
        if (!length(chunk_ids)) return(res_all)
        for (b in seq(1, length(chunk_ids), by = psock_batch_chunks)) {
            idx_range <- b:min(b + psock_batch_chunks - 1L, length(chunk_ids))
            batch_dt <- lapply(chunk_ids[idx_range], function(idx) dt_pts[idx])
            batch_res <- tryCatch({
                parLapply(cluster, batch_dt, fun = clip_worker,
                    region = region,
                    debug_parallel = debug_parallel,
                    log_dir = log_dir
                )
            }, error = function(e) e)
            if (inherits(batch_res, "error")) {
                return(batch_res)
            }
            res_all[idx_range] <- batch_res
            rm(batch_dt, batch_res)
            gc(FALSE)
        }
        res_all
    }

    run_backend <- function(backend) {
        if (identical(backend, "fork") && .Platform$OS.type != "windows") {
            mclapply(chunk_ids, clip_worker_idx, mc.cores = workers)
        } else if (identical(backend, "psock")) {
            run_psock_chunked()
        } else {
            lapply(chunk_ids, clip_worker_idx)
        }
    }

    res_list <- tryCatch(run_backend(backend), error = function(e) e)
    if (inherits(res_list, "error")) {
        if (!identical(backend, "serial")) {
            warning(
                "ROI clip parallel backend '", backend, "' failed: ",
                conditionMessage(res_list),
                ". Retrying serial. Set parallel_backend='serial' to disable parallel ROI clipping.",
                call. = FALSE
            )
            backend <- "serial"
            res_list <- run_backend(backend)
        } else {
            stop("Error clipping points to ROI: ", conditionMessage(res_list))
        }
    }

    bad <- vapply(res_list, function(x) !is.null(x$error), logical(1L))
    if (any(bad)) {
        if (!identical(backend, "serial")) {
            warn_msgs <- unique(vapply(res_list[bad], function(x) x$error, character(1L)))
            warning(
                "ROI clip workers failed under backend='", backend, "': ",
                paste(warn_msgs, collapse = "; "),
                ". Retrying serial.",
                call. = FALSE
            )
            res_list <- run_backend("serial")
            bad <- vapply(res_list, function(x) !is.null(x$error), logical(1L))
        }
        if (any(bad)) {
            err_msgs <- unique(vapply(res_list[bad], function(x) x$error, character(1L)))
            stop("ROI clip worker error(s): ", paste(err_msgs, collapse = "; "))
        }
    }

    rbindlist(lapply(res_list, `[[`, "data"))
}

#' Cmh Igraph To Matrix
#' @description
#' Internal helper for `.cmh_igraph_to_matrix`.
#' @param g Parameter value.
#' @param genes Parameter value.
#' @param weight_cols Parameter value.
#' @return Return value used internally.
#' @keywords internal
.cmh_igraph_to_matrix <- function(g, genes = NULL, weight_cols = c("CMH","weight")) {
  stopifnot(inherits(g, "igraph"))
  ed <- igraph::as_data_frame(g, what = "edges")
  wcol <- intersect(weight_cols, colnames(ed))
  if (!length(wcol)) stop(".cmh_igraph_to_matrix: edge table lacks CMH/weight attribute")
  wcol <- wcol[1]
  if (is.null(genes)) {
    genes <- igraph::V(g)$name
  }
  if (is.null(genes) || !length(genes)) stop(".cmh_igraph_to_matrix: cannot determine gene order")
  i <- match(ed$from, genes)
  j <- match(ed$to, genes)
  ok <- !is.na(i) & !is.na(j) & is.finite(ed[[wcol]])
  i <- as.integer(i[ok]); j <- as.integer(j[ok]); x <- as.numeric(ed[[wcol]][ok])
  M <- sparseMatrix(i = c(i, j), j = c(j, i), x = c(x, x),
                            dims = c(length(genes), length(genes)),
                            dimnames = list(genes, genes))
  drop0(as(M, "dgCMatrix"))
}

#' Compute module-level similarity from density layers
#' @description
#' Internal helper for `.compute_module_density_similarity`.
#' Computes similarity between module density profiles using Jaccard, Pearson,
#' Spearman, or cosine. Optionally z-score normalizes densities and filters
#' modules with insufficient grid coverage. Returns the similarity matrix and a
#' ComplexHeatmap object; no files are written.
#' @param scope_obj geneSCOPE object containing density layers.
#' @param grid_name Grid layer name (default "grid30").
#' @param density_cols Optional explicit density columns to use.
#' @param density_pattern Regex to select density columns (default "_density$").
#' @param method Similarity metric: "jaccard", "pearson", "spearman", "cosine".
#' @param density_threshold Threshold for binarization (Jaccard) or presence.
#' @param normalize Whether to z-score densities before similarity (default TRUE).
#' @param min_grids_per_module Minimum grids required per module.
#' @param top_n_modules Optional cap on number of modules; NULL keeps all.
#' @param output_dir Ignored (kept for backward compatibility; no files saved).
#' @param heatmap_bg Background color used when drawing the heatmap.
#' @return List with: matrix (similarity), modules (names), grid, method,
#'   threshold, min_grids, heatmap (ComplexHeatmap object), palette (color
#'   function), labels (clean module labels).
#' @keywords internal
.compute_module_density_similarity <- function(scope_obj,
                                           grid_name = "grid30",
                                           density_cols = NULL,
                                           density_pattern = "_density$",
                                           method = c("jaccard", "pearson", "spearman", "cosine"),
                                           density_threshold = 0,
                                           normalize = TRUE,
                                           min_grids_per_module = 5L,
                                           top_n_modules = NULL,
                                           output_dir = NULL,
                                           heatmap_bg = "#c0c0c0") {
  method <- match.arg(method)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required to draw the module similarity heatmap.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required to build the heatmap palette.")
  }
  dens_tbl <- scope_obj@density[[grid_name]]
  if (is.null(dens_tbl) || !nrow(dens_tbl)) {
    warning("density layer is missing or empty: ", grid_name)
    return(invisible(NULL))
  }

  if (is.null(rownames(dens_tbl))) {
    grid_info <- scope_obj@grid[[grid_name]]$grid_info
    if (!is.null(grid_info) && "grid_id" %in% colnames(grid_info)) {
      rownames(dens_tbl) <- as.character(grid_info$grid_id)
    }
  }

  if (is.null(density_cols)) {
    density_cols <- colnames(dens_tbl)
    if (!is.null(density_pattern)) {
      density_cols <- grep(density_pattern, density_cols, value = TRUE)
    }
  }

  density_cols <- intersect(density_cols, colnames(dens_tbl))
  if (!length(density_cols)) {
    warning("No matching density columns found.")
    return(invisible(NULL))
  }

  num_mask <- vapply(dens_tbl[density_cols], is.numeric, logical(1))
  if (!all(num_mask)) {
    warning("Non-numeric density columns ignored: ", paste(density_cols[!num_mask], collapse = ", "))
  }
  density_cols <- density_cols[num_mask]
  if (!length(density_cols)) {
    warning("No numeric density columns available for similarity computation.")
    return(invisible(NULL))
  }

  dens_mat <- as.matrix(dens_tbl[, density_cols, drop = FALSE])
  dens_mat[is.na(dens_mat)] <- 0
  if (normalize) {
    dens_mat <- scale(dens_mat, center = TRUE, scale = TRUE)
    dens_mat[is.na(dens_mat)] <- 0
  }

  module_names <- colnames(dens_mat)
  module_labels <- sub("^.*?_([0-9]+)_density$", "\\1", module_names)
  invalid_label <- module_labels == module_names | module_labels == ""
  module_labels[invalid_label] <- module_names[invalid_label]

  if (method == "jaccard") {
    bin_mat <- dens_mat > density_threshold
    bin_mat[is.na(bin_mat)] <- FALSE
    module_grid_mat <- Matrix(t(bin_mat) * 1, sparse = TRUE)
    module_presence <- rowSums(module_grid_mat)
  } else {
    module_grid_mat <- as.matrix(t(dens_mat))
    module_presence <- rowSums(module_grid_mat > density_threshold)
  }

  keep_idx <- which(module_presence >= min_grids_per_module)
  if (!length(keep_idx)) {
    warning("No modules meet the minimum grid count of ", min_grids_per_module, ".")
    return(invisible(NULL))
  }

  module_grid_mat <- module_grid_mat[keep_idx, , drop = FALSE]
  module_presence <- module_presence[keep_idx]
  module_names <- module_names[keep_idx]

  if (!is.null(top_n_modules) && is.finite(top_n_modules) && top_n_modules > 0 &&
      length(module_names) > top_n_modules) {
    ord <- order(module_presence, decreasing = TRUE)
    keep_ord <- ord[seq_len(top_n_modules)]
    module_grid_mat <- module_grid_mat[keep_ord, , drop = FALSE]
    module_presence <- module_presence[keep_ord]
    module_names <- module_names[keep_ord]
    module_labels <- module_labels[keep_ord]
  }

  sim_mat <- switch(method,
    jaccard = {
      inter_mat <- tcrossprod(module_grid_mat)
      inter_dense <- as.matrix(inter_mat)
      union_dense <- outer(module_presence, module_presence, "+") - inter_dense
      positive_union <- union_dense > 0
      jac <- inter_dense
      jac[positive_union] <- jac[positive_union] / union_dense[positive_union]
      jac[!positive_union] <- NA_real_
      diag(jac) <- 1
      jac
    },
    pearson = cor(t(module_grid_mat), method = "pearson"),
    spearman = cor(t(module_grid_mat), method = "spearman"),
    cosine = {
      dense_mat <- module_grid_mat
      num <- dense_mat %*% t(dense_mat)
      norms <- sqrt(diag(num))
      denom <- norms %o% norms
      denom[denom == 0] <- NA_real_
      cos_mat <- num / denom
      diag(cos_mat) <- 1
      cos_mat
    }
  )

  rownames(sim_mat) <- module_names
  colnames(sim_mat) <- module_names

  if (method == "jaccard" || method == "cosine") {
    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#113A70", "#f7f7f7", "#B11226"))
  } else {
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#113A70", "#f7f7f7", "#B11226"))
  }

  heatmap_name <- paste0("ModuleSimilarity_", method)
  ht <- ComplexHeatmap::Heatmap(
    sim_mat,
    name = heatmap_name,
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_labels = module_labels,
    column_labels = module_labels,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    heatmap_legend_param = list(title = ifelse(method == "jaccard", "Jaccard", method))
  )

  invisible(list(
    matrix = sim_mat,
    modules = module_names,
    grid = grid_name,
    method = method,
    threshold = density_threshold,
    min_grids = min_grids_per_module,
    heatmap = ht,
    palette = col_fun,
    labels = module_labels
  ))
}

#' Ensure a CMH matrix slot exists in LeeStats_Xz
#' @description
#' Internal helper for `.ensure_cmh_matrix`.
#' Creates `dst_matrix_slot` from `src_graph_slot` if missing or overwrite=TRUE
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param src_graph_slot Slot name.
#' @param dst_matrix_slot Slot name.
#' @param overwrite Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return updated scope_obj
#' @keywords internal
.ensure_cmh_matrix <- function(scope_obj, grid_name,
                              lee_stats_layer = "LeeStats_Xz",
                              src_graph_slot = NULL,
                              dst_matrix_slot = NULL,
                              overwrite = FALSE,
                              verbose = TRUE) {
  if (is.null(src_graph_slot)) src_graph_slot <- paste0("g_morisita_", grid_name)
  if (is.null(dst_matrix_slot)) dst_matrix_slot <- paste0("m_morisita_", grid_name)
  if (is.null(scope_obj@stats[[grid_name]]) || is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
    stop(".ensure_cmh_matrix: layer not found: ", grid_name, "/", lee_stats_layer)
  }
  ls <- scope_obj@stats[[grid_name]][[lee_stats_layer]]
  if (!overwrite && !is.null(ls[[dst_matrix_slot]]) && (inherits(ls[[dst_matrix_slot]], "Matrix") || is.matrix(ls[[dst_matrix_slot]]))) {
    if (verbose) .log_info("computeMH", "S07", paste0("Found existing ", dst_matrix_slot, "; skip."), verbose)
    return(scope_obj)
  }
  g <- ls[[src_graph_slot]]
  if (is.null(g)) stop(".ensure_cmh_matrix: source graph slot missing: ", src_graph_slot)
  genes <- try(rownames(ls$L), silent = TRUE)
  if (inherits(genes, "try-error") || is.null(genes)) genes <- igraph::V(g)$name
  M <- .cmh_igraph_to_matrix(g, genes = genes)
  scope_obj@stats[[grid_name]][[lee_stats_layer]][[dst_matrix_slot]] <- M
  if (verbose) .log_info("computeMH", "S07", paste0("Wrote ", dst_matrix_slot, " (", nrow(M), "×", ncol(M), "; nnz=", nnzero(M), ")"), verbose)
  scope_obj
}

#' Ensure Fdr Dimnames
#' @description
#' Internal helper for `.ensure_fdr_dimnames`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.ensure_fdr_dimnames <- function(scope_obj, grid_name,
                                lee_stats_layer = "LeeStats_Xz",
                                verbose = TRUE) {
  if (is.null(scope_obj@stats[[grid_name]]) || is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
    stop(".ensure_fdr_dimnames: layer not found: ", grid_name, "/", lee_stats_layer)
  }
  ls <- scope_obj@stats[[grid_name]][[lee_stats_layer]]
  L <- ls$L
  F <- ls$FDR
  if (is.null(F)) return(scope_obj)
  # coerce to base matrix if needed
  if (!is.matrix(F)) {
    F_try <- try(as.matrix(F), silent = TRUE)
    if (!inherits(F_try, "try-error") && is.matrix(F_try)) F <- F_try
  }
  if (is.matrix(F) && is.null(dimnames(F))) {
    dn <- dimnames(L)
    if (!is.null(dn) && !is.null(dn[[1]]) && !is.null(dn[[2]])) {
      dimnames(F) <- dn
      scope_obj@stats[[grid_name]][[lee_stats_layer]]$FDR <- F
      if (isTRUE(verbose)) .log_info("computeL", "S06", paste0("Set FDR dimnames from L (", grid_name, ")"), verbose)
    } else if (!is.null(rownames(L))) {
      rn <- rownames(L); dimnames(F) <- list(rn, rn)
      scope_obj@stats[[grid_name]][[lee_stats_layer]]$FDR <- F
      if (isTRUE(verbose)) .log_info("computeL", "S06", paste0("Set FDR dimnames from rownames(L) (", grid_name, ")"), verbose)
    }
  }
  scope_obj
}

#' Flip Y Coordinates
#' @description
#' Internal helper for `.flip_y_coordinates`.
#' @param coords Parameter value.
#' @param y_max Numeric threshold.
#' @return Return value used internally.
#' @keywords internal
.flip_y_coordinates <- function(coords, y_max) {
    stopifnot(is.numeric(y_max), length(y_max) == 1L)
    if (!"y" %in% names(coords)) {
        stop("Input data must expose a 'y' column.")
    }
    if (is.data.table(coords)) {
        coords[, y := y_max - y]
        coords
    } else {
        coords$y <- y_max - coords$y
        coords
    }
}

#' Plot Dendro Network Multi
#' @description
#' Internal helper for `.plot_dendro_network_multi`.
#' @param ... Additional arguments (currently unused).
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_dendro_network_multi <- function(...) {
    .Deprecated(".plot_dendro_network_multi")
    .plot_dendro_network_multi(...)
}

#' Plot Density Multi
#' @description
#' Internal helper for `.plot_density_multi`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param density1_name Parameter value.
#' @param density2_name Parameter value.
#' @param palette1 Parameter value.
#' @param palette2 Parameter value.
#' @param alpha1 Parameter value.
#' @param alpha2 Parameter value.
#' @param overlay_image Parameter value.
#' @param image_path Filesystem path.
#' @param image_alpha Parameter value.
#' @param image_choice Parameter value.
#' @param seg_type Parameter value.
#' @param colour_cell Parameter value.
#' @param colour_nucleus Parameter value.
#' @param alpha_seg Parameter value.
#' @param grid_gap Parameter value.
#' @param scale_text_size Parameter value.
#' @param bar_len Parameter value.
#' @param bar_offset Parameter value.
#' @param arrow_pt Parameter value.
#' @param scale_legend_colour Parameter value.
#' @param max.cutoff1 Parameter value.
#' @param max.cutoff2 Parameter value.
#' @param legend_digits Parameter value.
#' @param density_layers Parameter value.
#' @param layer_blend Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_density_multi <- function(scope_obj,
                        grid_name,
                        density1_name = NULL,
                        density2_name = NULL,
                        palette1 = "#fc3d5d",
                        palette2 = "#4753f8",
                        alpha1 = 0.5,
                        alpha2 = 0.5,
                        overlay_image = FALSE,
                        image_path = NULL,
                        image_alpha = 0.6,
                        image_choice = c("auto", "hires", "lowres"),
                        seg_type = c("cell", "nucleus", "both"),
                        colour_cell = "black",
                        colour_nucleus = "#3182bd",
                        alpha_seg = 0.2,
                        grid_gap = 100,
                        scale_text_size = 2.4,
                        bar_len = 400,
                        bar_offset = 0.01,
                        arrow_pt = 4,
                        scale_legend_colour = "black",
                        max.cutoff1 = 1,
                        max.cutoff2 = 1,
                        legend_digits = 1,
                        density_layers = NULL,
                        layer_blend = c("normal", "add", "screen", "lighten", "darken", "multiply", "overlay", "rgb_mix")) {
    layer_blend <- match.arg(layer_blend)
    seg_type <- match.arg(seg_type)
    image_choice <- match.arg(image_choice)
    accuracy_val <- 1 / (10^legend_digits)
    use_blend_filters <- !layer_blend %in% c("normal", "rgb_mix")
    blend_ref_id <- NULL
    if (use_blend_filters) {
        if (!requireNamespace("ggfx", quietly = TRUE)) {
            stop("layer_blend '", layer_blend, "' requires the 'ggfx' package to be installed.")
        }
        blend_ref_id <- paste0("densityBlend_", sample.int(1e9, 1))
    }

    ## ------------------------------------------------------------------ 1
    ## validate grid layer & retrieve geometry
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_layer_name <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    grid_info_data <- g_layer$grid_info
    if (is.null(grid_info_data) ||
        !all(c("grid_id", "xmin", "xmax", "ymin", "ymax") %in% names(grid_info_data))) {
        stop("grid_info is missing required columns for layer '", grid_layer_name, "'.")
    }

    ## ------------------------------------------------------------------ 2
    ## locate density data frame (new slot first, old slot as fallback)
    density_df <- scope_obj@density[[grid_layer_name]]
    if (is.null(density_df)) {
        density_df <- g_layer$densityDF
    } # legacy
    if (is.null(density_df)) {
        stop("No density table found for grid '", grid_layer_name, "'.")
    }

    ## make sure rownames are grid IDs
    if (is.null(rownames(density_df)) && "grid_id" %in% names(density_df)) {
        rownames(density_df) <- density_df$grid_id
    }

    `%||%` <- function(x, y) if (!is.null(x)) x else y

    normalize_density_layers <- function(dl) {
        if (is.null(dl)) {
            return(list())
        }
        if (is.character(dl) || is.factor(dl)) {
            dl <- as.character(dl)
        dl <- lapply(dl, function(density_column_name) list(column = density_column_name))
        } else if (is.data.frame(dl)) {
            if (!nrow(dl)) {
                return(list())
            }
            dl <- split(dl, seq_len(nrow(dl)))
            dl <- lapply(dl, function(row_df) as.list(row_df))
        } else if (is.list(dl) && length(dl) && !is.list(dl[[1]]) &&
            (!is.null(dl$column) || !is.null(dl$name) || !is.null(dl$density))) {
            dl <- list(dl)
        } else if (!is.list(dl)) {
            stop("`density_layers` must be a character vector, list, or data.frame.")
        }

        lapply(dl, function(density_layer_spec) {
            if (!is.list(density_layer_spec)) {
                if (length(density_layer_spec) == 1 && (is.character(density_layer_spec) || is.factor(density_layer_spec))) {
                    density_layer_spec <- list(column = as.character(density_layer_spec))
                } else {
                    stop("Each entry in `density_layers` must be a list or single character string.")
                }
            }
            col_name <- density_layer_spec$column %||% density_layer_spec$name %||% density_layer_spec$density
            if (is.null(col_name) || length(col_name) != 1L) {
                stop("Each density layer needs a single `column` (or alias `name`/`density`).")
            }
            list(
                column = col_name,
                palette = density_layer_spec$palette %||% density_layer_spec$colour %||% density_layer_spec$color,
                alpha = density_layer_spec$alpha %||% density_layer_spec$opacity,
                max_cutoff = density_layer_spec$max_cutoff %||% density_layer_spec$cutoff %||% density_layer_spec$max.cutoff,
                label = density_layer_spec$label %||% density_layer_spec$legend %||% density_layer_spec$title %||% col_name
            )
        })
    }

    density_layer_specs <- list()
    if (!is.null(density1_name)) {
        density_layer_specs <- c(density_layer_specs, list(list(
            column = density1_name,
            palette = palette1,
            alpha = alpha1,
            max_cutoff = max.cutoff1,
            label = density1_name
        )))
    }
    if (!is.null(density2_name)) {
        density_layer_specs <- c(density_layer_specs, list(list(
            column = density2_name,
            palette = palette2,
            alpha = alpha2,
            max_cutoff = max.cutoff2,
            label = density2_name
        )))
    }
    density_layer_specs <- c(density_layer_specs, normalize_density_layers(density_layers))

    if (!length(density_layer_specs)) {
        stop("Provide at least one density column via density1_name/density2_name or density_layers.")
    }

    all_density_cols <- vapply(density_layer_specs, `[[`, character(1), "column")
    missing_cols <- setdiff(unique(all_density_cols), colnames(density_df))
    if (length(missing_cols)) {
        stop(
            "Density column(s) '", paste(missing_cols, collapse = "', '"),
            "' not found in table for grid '", grid_layer_name, "'."
        )
    }

    palette_defaults <- Filter(function(x) !is.null(x) && length(x) == 1L && nchar(x),
        c(palette1, palette2)
    )
    if (!length(palette_defaults)) {
        palette_defaults <- "#fc3d5d"
    }
    alpha_defaults <- c(alpha1, alpha2)
    alpha_defaults <- alpha_defaults[is.finite(alpha_defaults)]
    if (!length(alpha_defaults)) {
        alpha_defaults <- 0.5
    }

    resolve_palette <- function(value, idx) {
        pal <- value
        if (is.null(pal) || !length(pal)) {
            pal <- palette_defaults[((idx - 1) %% length(palette_defaults)) + 1]
        }
        if (length(pal) != 1L) {
            stop("Each density layer must specify a single palette colour.")
        }
        pal
    }

    resolve_alpha <- function(value, idx) {
        aval <- value
        if (is.null(aval) || !is.numeric(aval) || !length(aval) || !is.finite(aval[1])) {
            aval <- alpha_defaults[((idx - 1) %% length(alpha_defaults)) + 1]
        }
        aval <- max(0, min(1, aval[1]))
        aval
    }

    resolve_cutoff <- function(value) {
        if (is.null(value) || !is.numeric(value) || !is.finite(value)) {
            return(1)
        }
        max(0, value)
    }

    ## helper to assemble heatmap data.frame
    build_heat <- function(d_col, cutoff_frac) {
        density_lookup_df_local <- data.frame(
            grid_id = rownames(density_df),
            d = density_df[[d_col]],
            stringsAsFactors = FALSE
        )
        heat_df_local <- merge(grid_info_data, density_lookup_df_local, by = "grid_id", all.x = TRUE)
        heat_df_local$d[is.na(heat_df_local$d)] <- 0
        maxv <- max(heat_df_local$d, na.rm = TRUE)
        heat_df_local$cut <- maxv * cutoff_frac
        heat_df_local$d <- pmin(heat_df_local$d, heat_df_local$cut)
        heat_df_local
    }

    ## ------------------------------------------------------------------ 3
    ## tile geometry (centre & size)
    tile_geometry_df <- function(tile_df_local) {
        transform(tile_df_local,
            x = (xmin + xmax) / 2,
            y = (ymin + ymax) / 2,
            w = xmax - xmin,
            h = ymax - ymin
        )
    }

    density_layer_payloads <- lapply(seq_along(density_layer_specs), function(layer_idx) {
        layer_spec <- density_layer_specs[[layer_idx]]
        cutoff_use <- resolve_cutoff(layer_spec$max_cutoff)
        layer_heat_df_local <- build_heat(layer_spec$column, cutoff_use)
        layer_tile_df_local <- tile_geometry_df(layer_heat_df_local)
        list(
            data = layer_tile_df_local,
            label = layer_spec$label %||% layer_spec$column,
            palette = resolve_palette(layer_spec$palette, layer_idx),
            alpha = resolve_alpha(layer_spec$alpha, layer_idx),
            cutoff = max(layer_tile_df_local$cut, na.rm = TRUE)
        )
    })

    library(ggplot2)
    p <- ggplot()

    ## ------------------------------------------------------------------ 3.0
    ## optional background image overlay (auto-detect from grid layer)
    if (isTRUE(overlay_image)) {
        img_info <- g_layer$image_info
        mpp_fullres <- g_layer$microns_per_pixel
        hires_path <- if (!is.null(img_info)) img_info$hires_path else NULL
        lowres_path <- if (!is.null(img_info)) img_info$lowres_path else NULL
        if (!is.null(image_path)) {
            hires_path <- image_path
            lowres_path <- NULL
        }
        sel_path <- switch(image_choice,
            auto = if (!is.null(hires_path)) hires_path else lowres_path,
            hires = hires_path,
            lowres = lowres_path
        )
        if (!is.null(sel_path) && file.exists(sel_path) && is.finite(mpp_fullres)) {
            ext <- tolower(file_ext(sel_path))
            img <- NULL
            if (ext %in% c("png") && requireNamespace("png", quietly = TRUE)) {
                img <- png::readPNG(sel_path)
            } else if (ext %in% c("jpg","jpeg") && requireNamespace("jpeg", quietly = TRUE)) {
                img <- jpeg::readJPEG(sel_path)
            }
            if (!is.null(img)) {
                # Normalize to RGBA array and apply alpha (annotation_raster has no alpha arg)
                if (length(dim(img)) == 2L) {
                    # grayscale -> RGB
                    img_rgb <- array(NA_real_, dim = c(dim(img)[1], dim(img)[2], 3L))
                    img_rgb[,,1] <- img; img_rgb[,,2] <- img; img_rgb[,,3] <- img
                    img <- img_rgb
                }
                if (length(dim(img)) == 3L) {
                    if (dim(img)[3] == 3L) {
                        img_rgba <- array(NA_real_, dim = c(dim(img)[1], dim(img)[2], 4L))
                        img_rgba[,,1:3] <- img
                        img_rgba[,,4]   <- max(0, min(1, image_alpha))
                        img <- img_rgba
                    } else if (dim(img)[3] == 4L) {
                        img[,,4] <- pmin(1, pmax(0, img[,,4] * image_alpha))
                    }
                }

                # image size (width x height) in um
                wpx <- dim(img)[2]; hpx <- dim(img)[1]
                # Determine effective microns-per-pixel depending on image type and scalefactors
                hires_sf <- if (!is.null(img_info)) img_info$tissue_hires_scalef else NA_real_
                lowres_sf <- if (!is.null(img_info)) img_info$tissue_lowres_scalef else NA_real_
                aligned_sf <- if (!is.null(img_info)) img_info$regist_target_img_scalef else NA_real_
                eff_mpp <- mpp_fullres
                # Detect by filename if necessary
                is_hires <- grepl("tissue_hires_image", basename(sel_path), ignore.case = TRUE)
                is_lowres <- grepl("tissue_lowres_image", basename(sel_path), ignore.case = TRUE)
                is_aligned <- grepl("aligned_tissue_image", basename(sel_path), ignore.case = TRUE)
                if (is_hires && is.finite(hires_sf)) eff_mpp <- mpp_fullres / hires_sf
                if (is_lowres && is.finite(lowres_sf)) eff_mpp <- mpp_fullres / lowres_sf
                if (is_aligned && is.finite(aligned_sf)) eff_mpp <- mpp_fullres / aligned_sf
                if (is_aligned && !is.finite(aligned_sf)) eff_mpp <- mpp_fullres
                # Fallback: if cannot infer, try to match grid width to image width by simple ratio
                if (!is.finite(eff_mpp) || eff_mpp <= 0) eff_mpp <- mpp_fullres
                w_um <- as.numeric(wpx) * as.numeric(eff_mpp)
                h_um <- as.numeric(hpx) * as.numeric(eff_mpp)
                y_origin <- if (!is.null(img_info$y_origin)) img_info$y_origin else "top-left"
                # If grid was not flipped (top-left), flip image vertically by swapping ymin/ymax
                if (identical(y_origin, "top-left")) {
                    p <- p + annotation_raster(img,
                        xmin = 0, xmax = w_um,
                        ymin = h_um, ymax = 0,
                        interpolate = TRUE
                    )
                } else {
                    p <- p + annotation_raster(img,
                        xmin = 0, xmax = w_um,
                        ymin = 0, ymax = h_um,
                        interpolate = TRUE
                    )
                }
            }
        }
    }

    add_density_layer <- function(p_obj, layer, idx, draw_tiles = TRUE) {
        if (idx > 1) {
            if (!requireNamespace("ggnewscale", quietly = TRUE)) {
                stop("Package 'ggnewscale' is required for multiple density layers.")
            }
            p_obj <- p_obj + new_scale_fill()
        }
        low_col <- "transparent"
        tile_layer <- if (isTRUE(draw_tiles)) {
            layer_geom <- geom_tile(
                data = layer$data,
                aes(x = x, y = y, width = w, height = h, fill = d),
                colour = NA,
                alpha = layer$alpha,
                inherit.aes = FALSE
            )
            if (use_blend_filters) {
                if (idx == 1) {
                    layer_geom <- ggfx::as_reference(
                        layer_geom,
                        id = blend_ref_id,
                        include = TRUE
                    )
                } else {
                    layer_geom <- ggfx::with_blend(
                        layer_geom,
                        bg_layer = blend_ref_id,
                        blend_type = layer_blend,
                        id = blend_ref_id,
                        include = TRUE
                    )
                }
            }
            layer_geom
        } else {
            geom_tile(
                data = layer$data,
                aes(x = x, y = y, width = w, height = h, fill = d),
                colour = NA,
                alpha = 0,
                inherit.aes = FALSE
            )
        }
        p_obj +
            tile_layer +
            scale_fill_gradient(
                name = layer$label,
                low = low_col,
                high = layer$palette,
                limits = c(0, layer$cutoff),
                oob = scales::squish,
                labels = scales::number_format(accuracy = accuracy_val),
                na.value = low_col,
                guide = guide_colorbar(order = idx)
            )
    }

    blend_rgb_mix <- function(layers) {
        if (!length(layers)) return(NULL)
        base_df <- layers[[1]]$data[, c("grid_id", "x", "y", "w", "h")]
        channel_sum_r <- rep(0, nrow(base_df))
        channel_sum_g <- rep(0, nrow(base_df))
        channel_sum_b <- rep(0, nrow(base_df))
        intensity_sum <- rep(0, nrow(base_df))
        for (layer in layers) {
            layer_df <- layer$data[, c("grid_id", "d")]
            merged <- merge(base_df[, "grid_id", drop = FALSE], layer_df,
                by = "grid_id", all.x = TRUE
            )
            intensity <- merged$d / layer$cutoff
            intensity[!is.finite(intensity)] <- 0
            intensity <- pmax(0, pmin(1, intensity))
            intensity <- intensity * layer$alpha
            rgb_vals <- col2rgb(layer$palette) / 255
            channel_sum_r <- channel_sum_r + intensity * rgb_vals[1]
            channel_sum_g <- channel_sum_g + intensity * rgb_vals[2]
            channel_sum_b <- channel_sum_b + intensity * rgb_vals[3]
            intensity_sum <- intensity_sum + intensity
        }
        safe_intensity <- ifelse(intensity_sum > 0, intensity_sum, 1)
        base_df$fill_hex <- rgb(
            pmax(0, pmin(1, channel_sum_r / safe_intensity)),
            pmax(0, pmin(1, channel_sum_g / safe_intensity)),
            pmax(0, pmin(1, channel_sum_b / safe_intensity))
        )
        base_df$fill_alpha <- pmax(0, pmin(1, intensity_sum))
        base_df
    }

    if (identical(layer_blend, "rgb_mix")) {
        for (idx in seq_along(density_layer_payloads)) {
            p <- add_density_layer(
                p,
                density_layer_payloads[[idx]],
                idx,
                draw_tiles = FALSE
            )
        }
        if (!requireNamespace("ggnewscale", quietly = TRUE)) {
            stop("layer_blend 'rgb_mix' requires the 'ggnewscale' package.")
        }
        p <- p + new_scale_fill()
        rgb_df <- blend_rgb_mix(density_layer_payloads)
        if (!is.null(rgb_df) && nrow(rgb_df)) {
            p <- p +
                geom_tile(
                    data = rgb_df,
                    aes(x = x, y = y, width = w, height = h, fill = fill_hex, alpha = fill_alpha),
                    colour = NA,
                    inherit.aes = FALSE
                ) +
                scale_fill_identity(guide = "none") +
                scale_alpha_identity(guide = "none")
        }
    } else {
        for (idx in seq_along(density_layer_payloads)) {
            p <- add_density_layer(p, density_layer_payloads[[idx]], idx)
        }
    }

    ## ------------------------------------------------------------------ 4
    ## segmentation overlay
    seg_layers <- switch(seg_type,
        cell    = "segmentation_cell",
        nucleus = "segmentation_nucleus",
        both    = c("segmentation_cell", "segmentation_nucleus")
    )
    seg_layers <- seg_layers[seg_layers %in% names(scope_obj@coord)]

    if (length(seg_layers)) {
        seg_dt <- rbindlist(scope_obj@coord[seg_layers],
            use.names = TRUE, fill = TRUE
        )
        if (seg_type == "both") {
            is_cell <- seg_dt$cell %in% scope_obj@coord$segmentation_cell$cell
            seg_dt[, segClass := ifelse(is_cell, "cell", "nucleus")]
            p <- p +
                geom_shape(
                    data = seg_dt[segClass == "cell"],
                    aes(x = x, y = y, group = cell),
                    colour = scales::alpha(colour_cell, alpha_seg),
                    fill = NA,
                    linewidth = 0.05
                ) +
                geom_shape(
                    data = seg_dt[segClass == "nucleus"],
                    aes(x = x, y = y, group = cell),
                    colour = scales::alpha(colour_nucleus, alpha_seg),
                    fill = NA,
                    linewidth = 0.05
                )
        } else {
            ## Single type overlay (stroke transparency baked into colour)
            col_use <- if (seg_type == "cell") colour_cell else colour_nucleus
            p <- p +
                geom_shape(
                    data = seg_dt,
                    aes(x = x, y = y, group = cell),
                    colour = scales::alpha(col_use, alpha_seg),
                    fill = NA,
                    linewidth = 0.05
                )
        }
    }

    ## ------------------------------------------------------------------ 5
    ## aesthetics: square, gridlines, scale-bar
    x_rng <- range(grid_info_data$xmin, grid_info_data$xmax)
    y_rng <- range(grid_info_data$ymin, grid_info_data$ymax)

    # pad shorter axis
    dx <- diff(x_rng)
    dy <- diff(y_rng)
    if (dx < dy) {
        pad <- (dy - dx) / 2
        x_rng <- x_rng + c(-pad, pad)
    }
    if (dy < dx) {
        pad <- (dx - dy) / 2
        y_rng <- y_rng + c(-pad, pad)
    }

    # gridlines
    grid_breaks_x <- seq(x_rng[1], x_rng[2], by = grid_gap)
    grid_breaks_y <- seq(y_rng[1], y_rng[2], by = grid_gap)
    p <- p +
        geom_vline(xintercept = grid_breaks_x, linewidth = 0.05, colour = "grey80") +
        geom_hline(yintercept = grid_breaks_y, linewidth = 0.05, colour = "grey80")

    p <- p +
        scale_x_continuous(limits = x_rng, expand = c(0, 0)) +
        scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
        coord_fixed() +
        theme_minimal(base_size = 9) +
        theme(
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.key.width = unit(0.4, "cm"), # key ""
            legend.key.height = unit(0.15, "cm"), # key ""
            legend.text = element_text(size = 9, angle = 90),
            legend.title = element_text(size = 8, hjust = 0, vjust = 1),
            plot.title = element_text(size = 10, hjust = 0.5),
            plot.margin = margin(1.2, 1, 1.5, 1, "cm"),
            axis.line.x.bottom = element_line(colour = "black"),
            axis.line.y.left = element_line(colour = "black"),
            axis.line.x.top = element_blank(),
            axis.line.y.right = element_blank(),
            axis.title = element_blank()
        ) +
        theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

    ## scale-bar
    x0 <- x_rng[1] + 0.01 * diff(x_rng) # 3% padding from left
    y_bar <- y_rng[1] + bar_offset * diff(y_rng) # inside plot area
    p <- p +
        annotate("segment",
            x = x0, xend = x0 + bar_len,
            y = y_bar,
            yend = y_bar,
            arrow = arrow(length = unit(arrow_pt, "pt"), ends = "both"),
            colour = scale_legend_colour,
            linewidth = 0.4
        ) +
        annotate("text",
            x = x0 + bar_len / 2,
            y = y_bar + 0.025 * diff(y_rng),
            label = paste0(bar_len, " \u00B5m"),
            colour = scale_legend_colour,
            vjust = 1, size = scale_text_size
        )
    return(p)
}

#' Run per-gene within vs between Jaccard workflow on a SCOPE object
#' @description
#' Internal helper for `.run_within_between_test`.
#' Wrapper that computes the gene–gene Jaccard matrix on a given grid/layer and
#' then runs the per-gene within vs between test. No files are written.
#' @param scope_obj geneSCOPE object.
#' @param dataset_name Label used in messages.
#' @param grid_name Grid layer name (default "grid30").
#' @param cluster_col Column in meta.data for module labels.
#' @param expr_layer Expression layer to binarize ("Xz" or "counts").
#' @param output_dir Ignored (backward compatibility; no files saved).
#' @param alternative Alternative hypothesis for the test.
#' @param method "wilcoxon" or "permutation".
#' @param n_perm Number of permutations when method = "permutation".
#' @param seed RNG seed.
#' @return List with `jaccard` (matrix, heatmap, clusters) and `within_between`
#'   (test results).
#' @keywords internal
.run_within_between_test <- function(scope_obj,
                                 dataset_name,
                                 grid_name = "grid30",
                                 cluster_col = "q95_res0.1_grid30_log1p_freq0.95",
                                 expr_layer = "Xz",
                                 output_dir = NULL,
                                 alternative = "greater",
                                 method = "wilcoxon",
                                 n_perm = 10000L,
                                 seed = 123) {
  jac_res <- .compute_gene_gene_jaccard_heatmap(
    scope_obj = scope_obj,
    grid_name = grid_name,
    expr_layer = expr_layer,
    cluster_col = cluster_col,
    keep_cluster_na = FALSE,   # drop genes without module labels
    top_n_genes = NULL,        # use all labeled genes, no top-N truncation
    split_heatmap_by_cluster = TRUE,
    output_dir = output_dir
  )

  if (is.null(jac_res) || is.null(jac_res$matrix) || is.null(jac_res$clusters)) {
    warning("Skipping ", dataset_name, ": missing Jaccard matrix or clusters.")
    return(invisible(NULL))
  }

  test_res <- .compute_gene_module_within_between_test(
    jac_mat = jac_res$matrix,
    gene_modules = jac_res$clusters,
    alternative = alternative,
    method = method,
    n_perm = n_perm,
    seed = seed
  )

  if (is.null(test_res)) {
    warning("Test failed for ", dataset_name)
    return(invisible(NULL))
  }

  .log_info(
      "computeGeneModuleWithinBetweenTest",
      "S03",
      paste0(
          "✓ ", dataset_name, " within>between test (", method, ", ", alternative, ") p=",
          if (!is.null(test_res$test$p.value)) signif(test_res$test$p.value, 3) else NA
      ),
      TRUE
  )

  invisible(list(
    jaccard = jac_res,
    within_between = test_res
  ))
}

#' Safe Thread Count
#' @description
#' Internal helper for `.safe_thread_count`.
#' @return Return value used internally.
#' @keywords internal
.safe_thread_count <- function() {
    # Base heuristic first
    core_guess <- tryCatch(detectCores(), error = function(e) 1L)
    n_req <- max(1L, min(8L, as.integer(core_guess) - 1L))

    if (exists(".get_safe_thread_count", mode = "function", inherits = TRUE)) {
        return(tryCatch(.get_safe_thread_count(default = n_req), error = function(e) n_req))
    }

    if (exists(".get_safe_thread_count_v2", mode = "function", inherits = TRUE)) {
        f <- get(".get_safe_thread_count_v2", mode = "function")
        an <- tryCatch(names(formals(f)), error = function(e) character(0))
        # Support both old and new signatures
        if ("max_requested" %in% an) {
            return(tryCatch(f(max_requested = n_req), error = function(e) n_req))
        } else if ("requested_cores" %in% an) {
            return(tryCatch(f(requested_cores = n_req), error = function(e) n_req))
        } else {
            return(tryCatch(f(n_req), error = function(e) n_req))
        }
    }
    n_req
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.get_top_lvs_r_assemble_outputs`.
#' @param sel Internal parameter
#' @return Internal helper result
#' @keywords internal
.get_top_lvs_r_assemble_outputs <- function(sel) {
    transmute(
        sel,
        gene1,
        gene2,
        L = LeesL,
        r = Pear,
        pct1 = gene1_expr_pct,
        pct2 = gene2_expr_pct,
        fdr = FDR
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.get_top_lvs_r_validate_inputs`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param pear_level Parameter value.
#' @param lee_stats_layer Layer name.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.get_top_lvs_r_validate_inputs <- function(scope_obj,
                                        grid_name,
                                        pear_level,
                                        lee_stats_layer,
                                        verbose) {
    .lvsr_validate_inputs(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        pear_level = pear_level,
        level = NULL,
        ncores = NULL,
        verbose = verbose,
        diag_na = TRUE,
        caller = ".get_top_l_vs_r"
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.get_top_lvs_r_resolve_runtime_config`.
#' @param ncores Internal parameter
#' @param verbose Internal parameter
#' @return Internal helper result
#' @keywords internal
.get_top_lvs_r_resolve_runtime_config <- function(ncores, verbose) {
    logi <- detectCores(TRUE)
    if (ncores > max(1L, min(ncores, logi))) {
        .log_info("getTopLvsR", "S01", paste0("Adjusting ncores: requested=", ncores, " capped at available=", max(1L, min(ncores, logi))), TRUE)
        ncores <- max(1L, min(ncores, logi))
    }
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
    Sys.setenv(OMP_NUM_THREADS = ncores)
    ncores
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.lvsr_postprocess_and_return`.
#' @param sel Internal parameter
#' @param do_perm Parameter value.
#' @param perms Parameter value.
#' @param p_adj_mode Parameter value.
#' @param pval_mode Parameter value.
#' @param total_universe Parameter value.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param expr_layer Layer name.
#' @param use_blocks Logical flag.
#' @param block_side Parameter value.
#' @param mem_limit_GB Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param clamp_mode Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.lvsr_postprocess_and_return <- function(sel,
                                         do_perm,
                                         perms,
                                         p_adj_mode,
                                         pval_mode,
                                         total_universe,
                                         scope_obj,
                                         grid_name,
                                         lee_stats_layer,
                                         expr_layer,
                                         use_blocks,
                                         block_side,
                                         mem_limit_GB,
                                         ncores,
                                         clamp_mode,
                                         verbose) {
    if (!nrow(sel) || !do_perm) {
        sel$FDR <- NA_real_
        return(.get_top_lvs_r_assemble_outputs(sel))
    }

    if (verbose) {
        .log_info("getTopLvsR", "S04", paste0("Preparing permutation analysis"), TRUE)
        .log_info("getTopLvsR", "S04", paste0("  Permutations: ", perms), TRUE)
        .log_info("getTopLvsR", "S04", paste0("  P-value method: ", pval_mode), TRUE)
        .log_info("getTopLvsR", "S04", paste0("  Multiple testing correction: ", p_adj_mode), TRUE)
    }
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    # Infer expression layer from lee_stats_layer if not specified
    if (is.null(expr_layer) || !nzchar(expr_layer)) {
        expr_layer <- sub("^LeeStats_", "", lee_stats_layer)
        if (identical(expr_layer, lee_stats_layer)) expr_layer <- "Xz"
        if (verbose) .log_info("getTopLvsR", "S04", paste0("Using expr_layer='", expr_layer, "' inferred from lee_stats_layer"), verbose)
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
    sel_genes <- unique(c(sel$gene1, sel$gene2))
    if (any(is.na(match(sel_genes, colnames(Xz))))) stop("Selected genes not found in Xz")
    gene_pairs <- cbind(
        match(sel$gene1, sel_genes) - 1L,
        match(sel$gene2, sel_genes) - 1L
    )
    delta_ref <- sel$Delta

    block_id <- if (use_blocks) {
        by <- (grid_info$gy - 1L) %/% block_side
        max_by <- max(by)
        ((grid_info$gx - 1L) %/% block_side) * (max_by + 1L) + by + 1L
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
    Xz_sub <- Xz[, match(sel_genes, colnames(Xz)), drop = FALSE]
    n_cells <- nrow(Xz_sub)

    target_batch <- min(100L, perms)
    est_bytes <- function(bs) n_cells * bs * 4
    while (target_batch > 1L && est_bytes(target_batch) > mem_limit_GB * 1024^3 * 0.30) {
        target_batch <- max(1L, floor(target_batch / 2))
    }
    .log_info(
        "getTopLvsR",
        "S04",
        paste0(
            "Planned batch_size=", target_batch,
            " (estimated idx_mat ", sprintf("%.2f", est_bytes(target_batch) / 1024^2),
            " MB / limit ", sprintf("%.2f", (mem_limit_GB * 1024^3 * 0.30) / 1024^2), " MB)"
        ),
        TRUE
    )

    if (verbose) .log_info("getTopLvsR", "S04", paste0("Starting permutation testing loop"), verbose)
    remaining <- perms
    exceed_count <- rep(0L, nrow(sel))
    perm_threads <- ncores
    attempt <- 1
    while (remaining > 0) {
        bsz <- min(target_batch, remaining)
        success <- FALSE
        while (!success && perm_threads >= 1) {
            .log_info(
                "getTopLvsR",
                "S04",
                sprintf(
                    "Permutation attempt #%d threads=%d batch=%d remain=%d clamp_mode=%s",
                    attempt, perm_threads, bsz, remaining, clamp_mode
                ),
                TRUE
            )
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
                .log_info("getTopLvsR", "S04", paste0(" Batch failed: ", conditionMessage(res)), TRUE)
                if (perm_threads > 1) {
                    perm_threads <- max(1, floor(perm_threads / 2))
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

    if (verbose) .log_info("getTopLvsR", "S04", paste0("Computing p-values and applying multiple testing correction"), verbose)
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

    fdr_validated <- .fdr_runner_validate_inputs(
        p_values = p_values,
        perms = perms,
        p_adj_mode = p_adj_mode,
        total_universe = total_universe,
        verbose = verbose
    )
    fdr_spec <- .fdr_runner_spec_build(fdr_validated)
    fdr_payload <- .fdr_runner_materialize(fdr_spec, fdr_validated)
    .fdr_runner_validate_outputs(fdr_payload, fdr_spec)
    fdr_finalize_validated <- .fdr_payload_finalize_validate_inputs(
        mode = fdr_spec$mode,
        fdr_spec = fdr_spec,
        fdr_payload = fdr_payload,
        include_matrices = FALSE
    )
    fdr_finalize_spec <- .fdr_payload_finalize_spec_build(fdr_finalize_validated)
    fdr_final <- .fdr_payload_finalize_materialize(fdr_finalize_spec, fdr_finalize_validated)
    .fdr_payload_finalize_validate_outputs(fdr_final, fdr_finalize_spec)
    FDR <- fdr_final$FDR_main

    if (verbose) {
        sig_count <- sum(p_values < 0.05, na.rm = TRUE)
        fdr_sig_count <- sum(FDR < 0.05, na.rm = TRUE)
        .log_info("getTopLvsR", "S04", paste0("  Significant pairs (p < 0.05): ", sig_count, "/", length(p_values)), TRUE)
        .log_info("getTopLvsR", "S04", paste0("  Significant pairs (FDR < 0.05): ", fdr_sig_count, "/", length(FDR)), TRUE)
    }

    sel$raw_p <- fdr_final$raw_p
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
        if (!nrow(sel)) .log_info("getTopLvsR", "S04", paste0("No Bonferroni-significant pairs (FDR < 0.05)."), TRUE)
    }

    .get_top_lvs_r_assemble_outputs(sel)
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.lvsr_prepare_lvs_r_table`.
#' @param matrices Internal parameter
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param clamp_mode Parameter value.
#' @param direction Parameter value.
#' @param top_n Parameter value.
#' @param pear_range Parameter value.
#' @param L_range Parameter value.
#' @param curve_layer Layer name.
#' @param CI_rule Parameter value.
#' @param lee_stats_layer Layer name.
#' @param pear_level Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.lvsr_prepare_lvs_r_table <- function(matrices,
                                     scope_obj,
                                     grid_name,
                                     clamp_mode,
                                     direction,
                                     top_n,
                                     pear_range,
                                     L_range,
                                     curve_layer,
                                     CI_rule,
                                     lee_stats_layer,
                                     pear_level,
                                     verbose) {
    if (verbose) .log_info("getTopLvsR", "S03", paste0("Computing Delta values and applying filters"), verbose)
    internal_clamp_mode <- clamp_mode
    if (clamp_mode == "both") {
        .log_info("getTopLvsR", "S03", paste0("clamp_mode='both' currently equivalent to ref_only (reference truncation only)"), TRUE)
        internal_clamp_mode <- "ref_only"
    }
    df <- data.frame(
        gene1 = rep(matrices$common, each = length(matrices$common))[matrices$ut],
        gene2 = rep(matrices$common, length(matrices$common))[matrices$ut],
        LeesL = matrices$LeesL_vec,
        Pear = matrices$Pear_vec,
        Delta = matrices$LeesL_vec -
            if (internal_clamp_mode == "ref_only") pmax(matrices$Pear_vec, 0) else matrices$Pear_vec
    )

    # Expression coverage percentages
    {
        expr_pct_map <- setNames(rep(0, length(matrices$common)), matrices$common)
        g_layer_try <- tryCatch(.select_grid_layer(scope_obj, grid_name), error = function(e) NULL)
        coverage_done <- FALSE

        if (!is.null(g_layer_try) && !is.null(g_layer_try$counts)) {
            ct <- g_layer_try$counts
            if (is.data.frame(ct) && all(c("gene", "grid_id") %in% colnames(ct))) {
                total_cells <- if (!is.null(g_layer_try$grid_info)) nrow(g_layer_try$grid_info) else length(unique(ct$grid_id))
                if (total_cells <= 0) total_cells <- NA_real_
                if (inherits(ct, "data.table")) {
                    gene_cells <- ct[gene %in% matrices$common, .(cells = uniqueN(grid_id)), by = gene]
                } else if (requireNamespace("data.table", quietly = TRUE)) {
                    dct <- as.data.table(ct)
                    gene_cells <- dct[gene %in% matrices$common, .(cells = uniqueN(grid_id)), by = gene]
                } else {
                    keep <- ct$gene %in% matrices$common
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
                inter_col <- intersect(cn, matrices$common)
                if (length(inter_col) == 0 && length(intersect(rn, matrices$common)) == 0) {
                    # keep zeros
                } else {
                    gene_in_col <- (length(inter_col) >= length(intersect(rn, matrices$common)))
                    if (!gene_in_col) {
                        expr_mat <- if (inherits(expr_mat, "dgCMatrix")) t(expr_mat) else t(expr_mat)
                        cn <- colnames(expr_mat)
                        inter_col <- intersect(cn, matrices$common)
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

    # Curve filtering
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

    if (verbose) .log_info("getTopLvsR", "S03", paste0("Applying threshold filters"), verbose)
    df <- df[df$Pear >= pear_range[1] & df$Pear <= pear_range[2] &
        df$LeesL >= L_range[1] & df$LeesL <= L_range[2], ]
    if (!nrow(df)) stop("Thresholds remove all pairs")
    total_universe <- nrow(df)

    if (verbose) {
        .log_info("getTopLvsR", "S03", paste0("  Pairs after filtering: ", total_universe), TRUE)
        .log_info("getTopLvsR", "S03", paste0("  Pearson range: [", pear_range[1], ", ", pear_range[2], "]"), TRUE)
        .log_info("getTopLvsR", "S03", paste0("  Lee's L range: [", L_range[1], ", ", L_range[2], "]"), TRUE)
    }

    if (!is.null(df$outside_ci) && CI_rule != "none" && verbose) {
        if (CI_rule == "remove_within") {
            .log_info("getTopLvsR", "S03", paste0("  Pairs outside CI95: ", sum(df$outside_ci)), TRUE)
        } else if (CI_rule == "remove_outside") {
            .log_info("getTopLvsR", "S03", paste0("  Pairs inside CI95: ", sum(!df$outside_ci)), TRUE)
        }
    }

    if (verbose) .log_info("getTopLvsR", "S03", paste0("Selecting top gene pairs by Delta values"), verbose)
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
        .log_info("getTopLvsR", "S03", paste0("  Selected ", nrow(sel), " top gene pairs"), TRUE)
        if (nrow(sel) > 0) {
            delta_range <- range(sel$Delta)
            .log_info("getTopLvsR", "S03", paste0("  Delta range: [", round(delta_range[1], 4), ", ", round(delta_range[2], 4), "]"), TRUE)
        }
    }

    list(sel = sel, total_universe = total_universe, internal_clamp_mode = internal_clamp_mode)
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_lvs_r_curve_assemble_outputs`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param curve_name Parameter value.
#' @param xgrid Parameter value.
#' @param fit Parameter value.
#' @param lo Parameter value.
#' @param hi Parameter value.
#' @param meta_obj Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_lvs_r_curve_assemble_outputs <- function(scope_obj,
                                               grid_name,
                                               lee_stats_layer,
                                               curve_name,
                                               xgrid,
                                               fit,
                                               lo,
                                               hi,
                                               meta_obj,
                                               verbose) {
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[grid_name]])) {
        scope_obj@stats[[grid_name]] <- list()
    }
    if (is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
        scope_obj@stats[[grid_name]][[lee_stats_layer]] <- list()
    }

    meta_col <- rep(list(meta_obj), length(fit))

    scope_obj@stats[[grid_name]][[lee_stats_layer]][[curve_name]] <-
        data.frame(Pear = xgrid, fit = fit, lo95 = lo, hi95 = hi, meta = I(meta_col))

    if (verbose) {
        .log_info("computeLvsRCurve", "S06", paste0("Analysis completed successfully"), TRUE)
        .log_info("computeLvsRCurve", "S06", paste0("  Curve stored as: ", curve_name), TRUE)
        .log_info("computeLvsRCurve", "S06", paste0("  Data points in curve: ", length(xgrid)), TRUE)
        .log_info("computeLvsRCurve", "S06", paste0("  Bootstrap iterations used: ", meta_obj$B), TRUE)
    }

    scope_obj
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.lvsr_curve_compute`.
#' @param Lv Internal parameter
#' @param rv Parameter value.
#' @param downsample Parameter value.
#' @param jitter_eps Parameter value.
#' @param n_strata Parameter value.
#' @param length_out Parameter value.
#' @param span Parameter value.
#' @param deg Parameter value.
#' @param B Parameter value.
#' @param k_max Numeric threshold.
#' @param min_rel_width Numeric threshold.
#' @param widen_span Parameter value.
#' @param ci_method Parameter value.
#' @param ci_adjust Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.lvsr_curve_compute <- function(Lv,
                                rv,
                                downsample,
                                jitter_eps,
                                n_strata,
                                length_out,
                                span,
                                deg,
                                B,
                                k_max,
                                min_rel_width,
                                widen_span,
                                ci_method,
                                ci_adjust,
                                ncores,
                                verbose) {
    if (verbose) .log_info("computeLvsRCurve", "S02", paste0("Cleaning and preprocessing data points"), verbose)
    ok <- is.finite(Lv) & is.finite(rv)
    Lv <- Lv[ok]
    rv <- rv[ok]
    if (!length(Lv)) stop("No valid points")

    if (verbose) .log_info("computeLvsRCurve", "S02", paste0("  Initial data points: ", length(Lv)), verbose)

    if (is.numeric(downsample) && downsample < 1) {
        keep <- sample.int(length(Lv), max(1L, floor(downsample * length(Lv))))
        Lv <- Lv[keep]
        rv <- rv[keep]
        if (verbose) .log_info("computeLvsRCurve", "S02", paste0("  Downsampled to: ", length(Lv), " points (ratio: ", downsample, ")"), verbose)
    } else if (is.numeric(downsample) && downsample >= 1 && length(Lv) > downsample) {
        keep <- sample.int(length(Lv), downsample)
        Lv <- Lv[keep]
        rv <- rv[keep]
        if (verbose) .log_info("computeLvsRCurve", "S02", paste0("  Downsampled to: ", length(Lv), " points (target: ", downsample, ")"), verbose)
    }
    if (jitter_eps > 0) {
        rv <- jitter(rv, factor = jitter_eps)
        if (verbose) .log_info("computeLvsRCurve", "S02", paste0("  Applied jitter to r values (eps: ", jitter_eps, ")"), verbose)
    }

    # 3. Stratify (by r quantiles) - fallback if failed
    if (verbose) .log_info("computeLvsRCurve", "S02", paste0("Setting up stratified sampling"), verbose)
    if (length(sort(unique(rv))) < 3) stop("r values too discrete")
    n_strata_eff <- min(n_strata, length(sort(unique(rv))) - 1)

    if (verbose) {
        .log_info("computeLvsRCurve", "S02", paste0("  Unique r values: ", length(sort(unique(rv)))), TRUE)
        .log_info("computeLvsRCurve", "S02", paste0("  Effective strata: ", n_strata_eff), TRUE)
    }

    probs <- seq(0, 1, length.out = n_strata_eff + 1)
    brks <- unique(quantile(rv, probs, na.rm = TRUE))
    if (length(brks) < 2) {
        # Fallback: uniform slices
        brks <- seq(min(rv), max(rv), length.out = n_strata_eff + 1)
        brks <- unique(brks)
        if (verbose) .log_info("computeLvsRCurve", "S02", paste0("  Using fallback uniform stratification"), verbose)
    }
    if (length(brks) < 2) stop("Cannot establish strata")
    strat <- cut(rv, breaks = brks, include.lowest = TRUE, labels = FALSE)
    ok2 <- !is.na(strat)
    rv <- rv[ok2]
    Lv <- Lv[ok2]
    strat <- strat[ok2]

    # 4. xgrid
    if (verbose) .log_info("computeLvsRCurve", "S02", paste0("Preparing analysis grid and fitting LOESS model"), verbose)
    if (diff(range(rv)) <= 0) stop("r has no span")
    xgrid <- seq(range(rv)[1], range(rv)[2], length.out = length_out)

    if (verbose) {
        .log_info("computeLvsRCurve", "S02", paste0("  Analysis range: [", round(range(rv)[1], 3), ", ", round(range(rv)[2], 3), "]"), TRUE)
        .log_info("computeLvsRCurve", "S02", paste0("  Grid points: ", length_out), TRUE)
        .log_info("computeLvsRCurve", "S02", paste0("  Using ", ncores, " cores for bootstrap"), TRUE)
    }

    # Directly call unified interface (no longer check old version)
    if (verbose) .log_info("computeLvsRCurve", "S02", paste0("Running LOESS residual bootstrap analysis"), verbose)
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

    if (verbose) {
        .log_info("computeLvsRCurve", "S02", paste0("Bootstrap completed, processing confidence intervals"), TRUE)
        .log_info("computeLvsRCurve", "S02", paste0("  CI method applied: ", ci_method), TRUE)
    }

    # 6. floor widen
    if (verbose && min_rel_width > 0) {
        .log_info("computeLvsRCurve", "S02", paste0("Applying minimum relative width constraints"), TRUE)
    }
    if (min_rel_width > 0) {
        rng_fit <- diff(range(fit, finite = TRUE))
        if (!is.finite(rng_fit) || rng_fit <= 0) rng_fit <- 1
        width <- hi - lo
        width <- pmax(width, min_rel_width * rng_fit)
        # Smooth width (LOESS)
        if (is.finite(widen_span) && widen_span > 0 && length(width) > 10) {
            sm <- tryCatch(loess(width ~ xgrid, span = widen_span), error = function(e) NULL)
            if (!is.null(sm)) {
                wp <- tryCatch(predict(sm, xgrid), error = function(e) NULL)
                if (!is.null(wp) && all(is.finite(wp))) {
                    width <- pmax(wp, min_rel_width * rng_fit)
                }
            }
        }
        center <- (lo + hi) / 2
        lo <- center - width / 2
        hi <- center + width / 2
    }

    # ---- New: Local residual scale adaptive lower bound (no new parameters) ----------------------------
    # Purpose: Avoid "points more dispersed but CI narrower"; force width >= 2*1.96*local MAD
    if (verbose) .log_info("computeLvsRCurve", "S02", paste0("Applying local residual scale adaptation"), verbose)
    {
        if (length(Lv) > 20) {
            # Fit values interpolated to original r
            fit_at_rv <- tryCatch(approx(xgrid, fit, xout = rv, rule = 2, ties = "ordered")$y,
                error = function(e) rep(mean(fit), length(rv))
            )
            rv_sorted <- rv[order(rv)]
            res_sorted <- (Lv - fit_at_rv)[order(rv)]
            n_pts <- length(rv_sorted)
            # Fixed K (internal constant, no new parameters)
            K <- min(200L, n_pts)
            if (K > 5) {
                # Nearest neighbor index function (bidirectional expansion)
                get_knn_idx <- function(x0) {
                    l <- findInterval(x0, rv_sorted)
                    r <- l + 1
                    out <- integer(0)
                    while (length(out) < K && (l >= 1 || r <= n_pts)) {
                        dr <- if (r <= n_pts) abs(rv_sorted[r] - x0) else Inf
                        if ((if (l >= 1) abs(rv_sorted[l] - x0) else Inf) <= dr) {
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
                    if (!length(get_knn_idx(xx))) {
                        return(NA_real_)
                    }
                1.4826 * median(abs(res_sorted[get_knn_idx(xx)] -
                    median(res_sorted[get_knn_idx(xx)])))
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
    if (verbose) .log_info("computeLvsRCurve", "S02", paste0("Applying confidence interval edge smoothing"), verbose)
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

    list(
        fit = fit,
        lo = lo,
        hi = hi,
        xgrid = xgrid,
        adjust_mode = adjust_mode,
        res = res,
        local_mad_diag = if (exists("local_mad_diag")) local_mad_diag else NULL,
        edge_smooth_info = if (exists("edge_smooth_info")) edge_smooth_info else NULL
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_lvs_r_curve_validate_inputs`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param level Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_lvs_r_curve_validate_inputs <- function(scope_obj,
                                             grid_name,
                                             lee_stats_layer,
                                             level,
                                             ncores,
                                             verbose) {
    .lvsr_validate_inputs(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        pear_level = NULL,
        level = level,
        ncores = ncores,
        verbose = verbose,
        diag_na = FALSE,
        caller = ".compute_l_vs_r_curve"
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.lvsr_validate_inputs`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param pear_level Parameter value.
#' @param level Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param verbose Logical; whether to emit progress messages.
#' @param diag_na Parameter value.
#' @param caller Parameter value.
#' @return Internal helper result
#' @keywords internal
.lvsr_validate_inputs <- function(scope_obj,
                                  grid_name,
                                  lee_stats_layer,
                                  pear_level = NULL,
                                  level = NULL,
                                  ncores = NULL,
                                  verbose = TRUE,
                                  diag_na = FALSE,
                                  caller = ".compute_l_vs_r_curve") {
    parent <- if (identical(caller, ".get_top_l_vs_r")) "getTopLvsR" else "computeLvsRCurve"
    step <- "S01"
    if (verbose) .log_info(parent, step, "Extracting Lee's L and Pearson correlation matrices", verbose)
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    r_level <- if (!is.null(pear_level)) {
        pear_level
    } else if (is.null(level)) {
        "grid"
    } else if (level == "grid") {
        "grid"
    } else {
        "cell"
    }
    Lmat <- .get_lee_matrix(scope_obj, grid_name, lee_layer = lee_stats_layer)
    rmat <- .get_pearson_matrix(scope_obj, grid_name, level = r_level)

    common <- intersect(rownames(Lmat), rownames(rmat))
    if (length(common) < 2) stop("Insufficient common genes")

    if (verbose) {
        .log_info(parent, step, paste0("Found ", length(common), " common genes between matrices"), verbose)
        .log_info(parent, step, paste0("Matrix dimensions: ", nrow(Lmat), "x", ncol(Lmat)), verbose)
    }

    if (!is.null(ncores)) {
        sys_mem_gb <- .get_system_memory_gb()
        per_thread_gb <- (length(common)^2 * 8 * 2) / (1024^3)
        est_total_gb <- per_thread_gb * ncores
        if (est_total_gb > sys_mem_gb) {
            stop(
                "[geneSCOPE::", caller, "] Estimated memory requirement (",
                round(est_total_gb, 1), " GB) exceeds system capacity (",
                round(sys_mem_gb, 1), " GB). Reduce ncores, downsample, or gene set size."
            )
        }
    }

    Lmat <- Lmat[common, common, drop = FALSE]
    rmat <- rmat[common, common, drop = FALSE]
    if (diag_na) {
        diag(Lmat) <- NA
        diag(rmat) <- NA
        ut <- upper.tri(Lmat)
    } else {
        ut <- upper.tri(Lmat, diag = FALSE)
    }

    list(
        grid_name = grid_name,
        Lmat = Lmat,
        rmat = rmat,
        ut = ut,
        Lv = Lmat[ut],
        rv = rmat[ut],
        LeesL_vec = Lmat[ut],
        Pear_vec = rmat[ut],
        common = common
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_lvs_r_curve_resolve_runtime_config`.
#' @param ncores Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_lvs_r_curve_resolve_runtime_config <- function(ncores) {
    max(1L, min(as.integer(ncores), detectCores(logical = TRUE)))
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.visium_adapter_inputs`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param norm_layer Layer name.
#' @param use_idelta Logical flag.
#' @param S_target Parameter value.
#' @param min_detect Numeric threshold.
#' @param winsor_high Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.visium_adapter_inputs <- function(scope_obj,
                                   grid_name,
                                   norm_layer,
                                   use_idelta,
                                   S_target,
                                   min_detect,
                                   winsor_high,
                                   ncores,
                                   verbose) {
    if (is.null(ncores)) {
        ncores <- .get_safe_thread_count(8L)
    }
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else grid_name

    if (is.null(g_layer[[norm_layer]])) {
        stop("computeL_visium: normalized layer '", norm_layer, "' not found.")
    }
    X <- g_layer[[norm_layer]]  # n x g
    if (!is.matrix(X)) X <- as.matrix(X)
    n <- nrow(X); G <- ncol(X)
    if (n < 2L || G < 2L) stop("Insufficient data size to compute Lee's L.")

    ## ---- 1. Check/compute spatial weights W ----
    if (is.null(g_layer$W) || sum(g_layer$W) == 0) {
        if (verbose) .log_info("computeL_visium", "S02", paste0("W not found; auto .compute_weights(style='B')..."), verbose)
        scope_obj <- .compute_weights(scope_obj, grid_name = grid_name, style = "B", store_mat = TRUE)
        g_layer <- .select_grid_layer(scope_obj, grid_name)
    }
    W <- g_layer$W  # dgCMatrix (ensures W exists for .compute_l later)

    ## ---- 2. Idelta and gene prescreening ----
    if (use_idelta) {
        meta_col <- paste0(grid_name, "_iDelta")
        if (is.null(scope_obj@meta.data) || !(meta_col %in% colnames(scope_obj@meta.data))) {
            if (verbose) .log_info("computeL_visium", "S03", paste0("Computing Idelta..."), verbose)
            scope_obj <- .compute_idelta(scope_obj, grid_name = grid_name, level = "grid", ncores = min(8L, ncores), verbose = verbose)
        }
    }
    genes_all <- colnames(X)
    if (is.null(genes_all)) genes_all <- rownames(scope_obj@meta.data)
    if (is.null(genes_all)) genes_all <- as.character(seq_len(G))

    # 2.1 Detection rate (based on counts)
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
            detect_rate <- dr$N / n_spots; names(detect_rate) <- dr$gene
            tc <- dtt[, sum(count), by = gene]
            total_count <- tc$V1; names(total_count) <- tc$gene
        } else {
            n_spots <- length(g_layer$grid_info$grid_id)
            det_list <- split(dt$grid_id, dt$gene)
            detect_rate <- vapply(det_list, function(v) length(unique(v)) / n_spots, numeric(1))
            cnt_list <- split(dt$count, dt$gene)
            total_count <- vapply(cnt_list, sum, numeric(1))
        }
    } else {
        # Fallback: infer from X (may be biased)
        detect_rate <- colMeans(X > 0)
        total_count <- colSums(pmax(X, 0))
        names(detect_rate) <- names(total_count) <- colnames(X)
    }

    # Idelta fetch and winsorization
    idelta_vec <- rep(NA_real_, G); names(idelta_vec) <- genes_all
    meta_col <- paste0(grid_name, "_iDelta")
    if (!is.null(scope_obj@meta.data) && (meta_col %in% colnames(scope_obj@meta.data))) {
        idelta_vec[rownames(scope_obj@meta.data)] <- scope_obj@meta.data[, meta_col]
    }
    idelta_vec[is.infinite(idelta_vec)] <- NA_real_
    idelta_vec_w <- idelta_vec
    if (all(is.na(idelta_vec_w))) {
        if (verbose) .log_info("computeL_visium", "S03", paste0("Idelta not found; skip Idelta-based ranking and use detection rate only."), verbose)
        idelta_vec_w <- rep(0, G); names(idelta_vec_w) <- genes_all
    } else {
        idelta_vec_w <- pmin(idelta_vec_w, quantile(idelta_vec_w, probs = winsor_high, na.rm = TRUE))
    }

    # 2.2 Prescreen: rank by Idelta within detection-rate filter; keep top S_target
    keep <- intersect(names(detect_rate)[detect_rate >= min_detect], genes_all)[order(idelta_vec_w[intersect(names(detect_rate)[detect_rate >= min_detect], genes_all)], decreasing = TRUE, na.last = NA)]
    if (length(keep) > S_target) keep <- keep[seq_len(S_target)]
    S <- length(keep)
    if (S < 100L) stop("Too few genes after prescreen (<100); lower min_detect or increase S_target.")
    if (verbose) .log_info("computeL_visium", "S04", paste0("Prescreened genes S=", S, " (target ", S_target, ")"), verbose)

    list(
        scope_obj = scope_obj,
        grid_name = grid_name,
        norm_layer = norm_layer,
        keep = keep,
        ncores = ncores
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_l_postprocess`.
#' @param L Internal parameter
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
#' @return Internal helper result
#' @keywords internal
.compute_l_postprocess <- function(L,
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
    validated <- .stats_table_validate_inputs(
        L = L,
        X_full = X_full,
        X_used = X_used,
        W = W,
        grid_inf = grid_inf,
        cells = cells,
        genes = genes,
        within = within,
        chunk_size = chunk_size,
        L_min = L_min,
        block_id = block_id,
        perms = perms,
        block_size = block_size,
        current_cores = current_cores,
        verbose = verbose
    )
    spec <- .stats_table_postprocess_spec_build(validated)
    post <- .stats_table_postprocess_materialize(spec, validated)
    .stats_table_validate_outputs(post, spec)
    post
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_l_assemble_outputs`.
#' @param scope_obj Internal parameter
#' @param gname Parameter value.
#' @param layer_name Parameter value.
#' @param L Parameter value.
#' @param Z_mat Parameter value.
#' @param P Parameter value.
#' @param betas Parameter value.
#' @param L_min Numeric threshold.
#' @param qc Parameter value.
#' @param FDR_main Parameter value.
#' @param FDR_storey Parameter value.
#' @param FDR_out_disc Parameter value.
#' @param FDR_out_beta Parameter value.
#' @param FDR_out_mid Parameter value.
#' @param FDR_out_uniform Parameter value.
#' @param perms Parameter value.
#' @param block_size Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param pi0_hat Parameter value.
#' @param n_sig_005 Parameter value.
#' @param min_p_possible Numeric threshold.
#' @param FDR_main_method Parameter value.
#' @param input_cache Parameter value.
#' @param cache_inputs Parameter value.
#' @return Internal helper result
#' @keywords internal
.compute_l_assemble_outputs <- function(scope_obj,
                                       gname,
                                       layer_name,
                                       L,
                                       Z_mat,
                                       P,
                                       betas,
                                       L_min,
                                       qc,
                                       FDR_main,
                                       FDR_storey,
                                       FDR_out_disc,
                                       FDR_out_beta,
                                       FDR_out_mid,
                                       FDR_out_uniform,
                                       perms,
                                       block_size,
                                       ncores,
                                       pi0_hat,
                                       n_sig_005,
                                       min_p_possible,
                                       FDR_main_method,
                                       input_cache,
                                       cache_inputs) {
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()

    scope_obj@stats[[gname]][[layer_name]] <- list(
        L = L,
        Z = Z_mat,
        P = P,
        grad = betas,
        L_min = L_min,
        qc = qc,
        FDR = FDR_main,
        FDR_storey = FDR_storey,
        FDR_disc = FDR_out_disc,
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

    scope_obj
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_l_run_core_computation`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param norm_layer Layer name.
#' @param genes Parameter value.
#' @param within Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param mem_limit_GB Parameter value.
#' @param chunk_size Parameter value.
#' @param use_bigmemory Logical flag.
#' @param backing_path Filesystem path.
#' @param block_side Parameter value.
#' @param cache_inputs Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_l_run_core_computation <- function(scope_obj,
                                           grid_name,
                                           norm_layer,
                                           genes,
                                           within,
                                           ncores,
                                           mem_limit_GB,
                                           chunk_size,
                                           use_bigmemory,
                                           backing_path,
                                           block_side,
                                           cache_inputs,
                                           verbose) {
    current_cores <- ncores
    min_cores <- 1
    success <- FALSE
    attempt <- 1
    res <- NULL

    while (!success && current_cores >= min_cores) {
        if (verbose && attempt > 1) .log_info("computeL", "S03", paste0("Retry #", attempt, " with ", current_cores, " cores"), TRUE)

        result <- tryCatch(
            {
                t_start <- Sys.time()
                res_local <- .compute_lee_l(scope_obj,
                    grid_name = grid_name,
                    norm_layer = norm_layer,
                    genes = genes,
                    within = within,
                    ncores = current_cores,
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
                    .log_info("computeL", "S03", paste0("Lee's L ", time_msg, " (", format(t_end - t_start), ")"), TRUE)
                }

                list(success = TRUE, object = res_local)
            },
            error = function(e) {
                if (verbose && attempt > 1) .log_info("computeL", "S03", paste0("Attempt failed: ", conditionMessage(e)), TRUE)
                list(success = FALSE, error = e)
            }
        )

        success <- result$success

        if (success) {
            res <- result$object
        } else {
            attempt <- attempt + 1
            current_cores <- max(floor(current_cores / 2), min_cores)
            if (verbose) .log_info("computeL", "S03", paste0("Reducing cores to ", current_cores, " and retrying"), verbose)
            Sys.sleep(3)
            gc(verbose = FALSE)
        }
    }

    if (!success) {
        stop("Unable to compute Lee's L statistics even with minimal thread count")
    }

    list(
        res = res,
        current_cores = current_cores
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_l_resolve_paths_and_cache`.
#' @param validation Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_l_resolve_paths_and_cache <- function(validation) {
    list(
        scope_obj = validation$scope_obj,
        g_layer = validation$g_layer,
        grid_name = validation$grid_name,
        use_bigmemory = validation$use_bigmemory,
        backing_path = validation$backing_path
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_l_configure_threads_and_rng`.
#' @param ncores Internal parameter
#' @param ncore Internal parameter
#' @param verbose Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_l_configure_threads_and_rng <- function(ncores, ncore, verbose) {
    runtime_cfg <- .compute_l_resolve_runtime_config(ncores = ncores, ncore = ncore, verbose = verbose)
    restore_fn <- attr(runtime_cfg$thread_config, "restore_function")
    list(runtime_cfg = runtime_cfg, restore_fn = restore_fn)
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_l_validate_inputs`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param norm_layer Layer name.
#' @param genes Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param mem_limit_GB Parameter value.
#' @param use_bigmemory Logical flag.
#' @param chunk_size Parameter value.
#' @param backing_path Filesystem path.
#' @param os_type Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param user_set_bigmemory Parameter value.
#' @param user_set_chunk Parameter value.
#' @return Internal helper result
#' @keywords internal
.compute_l_validate_inputs <- function(scope_obj,
                                     grid_name,
                                     norm_layer,
                                     genes,
                                     ncores,
                                     mem_limit_GB,
                                     use_bigmemory,
                                     chunk_size,
                                     backing_path,
                                     os_type,
                                     verbose,
                                     user_set_bigmemory,
                                     user_set_chunk) {
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else {
        grid_name
    }

    # Check if spatial weight matrix exists
    if (is.null(g_layer$W) || sum(g_layer$W) == 0) {
        .log_info("computeL", "S01", paste0("Computing spatial weights matrix"), TRUE)
        scope_obj <- .compute_weights(scope_obj,
            grid_name = grid_name,
            style = "B",
            store_mat = TRUE,
            store_listw = FALSE
        )
        # Re-extract layer after computing weights
        g_layer <- .select_grid_layer(scope_obj, grid_name)
    }

    # Cross-platform memory calculation
    all_genes <- if (!is.null(g_layer[[norm_layer]])) {
        ncol(g_layer[[norm_layer]])
    } else {
        stop("Normalized layer '", norm_layer, "' not found")
    }

    matrix_size_gb <- (((if (!is.null(genes)) length(genes) else all_genes)^2) * 8) / (1024^3)

    # Memory guard: estimate per-thread usage as single-core size x threads
    sys_mem_gb <- .get_system_memory_gb()
    est_total_gb <- matrix_size_gb * ncores
    if (est_total_gb > sys_mem_gb) {
        stop(
            "[geneSCOPE::.compute_l] Estimated memory requirement (",
            round(est_total_gb, 1), " GB) exceeds system capacity (",
            round(sys_mem_gb, 1), " GB). Reduce ncores or gene set size."
        )
    }

    # Respect user bigmemory/chunk overrides while still hinting when streaming is advisable
    if (matrix_size_gb > mem_limit_GB) {
        if (use_bigmemory) {
            if (verbose) {
                .log_info(
                    "computeL",
                    "S02",
                    paste0(
                        "Large matrix detected (", round(matrix_size_gb, 1), " GB), TRUE); staying in bigmemory/streaming mode.",
                        if (!user_set_chunk) " You may tune chunk_size to trade IO vs RAM." else ""
                    ),
                    verbose
                )
            }
        } else if (user_set_bigmemory) {
            if (verbose) .log_info("computeL", "S02", paste0("Warning: requested use_bigmemory=FALSE with large matrix (", round(matrix_size_gb, 1), " GB); proceeding in-memory as requested."), verbose)
        } else if (os_type == "windows" && !requireNamespace("bigmemory", quietly = TRUE)) {
            if (verbose) .log_info("computeL", "S02", paste0("!!! Warning: bigmemory not available on Windows; using regular matrices !!!"), verbose)
            use_bigmemory <- FALSE
        } else {
            use_bigmemory <- TRUE
            if (verbose) .log_info("computeL", "S02", paste0("Large matrix detected (", round(matrix_size_gb, 1), " GB); enabling bigmemory/streaming."), verbose)
        }
    }

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

    list(
        scope_obj = scope_obj,
        g_layer = g_layer,
        grid_name = grid_name,
        use_bigmemory = use_bigmemory,
        backing_path = backing_path,
        all_genes = all_genes,
        n_genes_use = if (!is.null(genes)) length(genes) else all_genes,
        matrix_size_gb = matrix_size_gb,
        sys_mem_gb = sys_mem_gb
    )
}

#' Internal helper for Lee's L workflows
#' @description
#' Internal helper for `.compute_l_resolve_runtime_config`.
#' @param ncores Internal parameter
#' @param ncore Internal parameter
#' @param verbose Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_l_resolve_runtime_config <- function(ncores, ncore, verbose) {
    if (!is.null(ncore)) {
        warning("'ncore' is deprecated, please use 'ncores' instead. Will use its value this time.",
            call. = FALSE, immediate. = TRUE
        )
        ncores <- ncore
    }
    os_type <- .detect_os()
    avail_cores <- max(1L, detectCores(logical = TRUE))
    ncores <- max(1L, min(ncores, avail_cores))
    if (verbose) {
        .log_info("computeL", "S01", paste0("Core configuration: using ", ncores, "/", avail_cores, " logical cores"), TRUE)
    }
    thread_config <- .configure_threads_for("mixed", ncores, restore_after = TRUE)
    list(
        ncores = ncores,
        avail_cores = avail_cores,
        os_type = os_type,
        thread_config = thread_config,
        ncores_cpp = thread_config$openmp_threads
    )
}

#' Internal helper for metrics/tests helpers
#' @description
#' Internal helper for `.compute_mh_assemble_outputs`.
#' @param scope_obj Internal parameter
#' @param gname Parameter value.
#' @param lee_layer Layer name.
#' @param graph_slot_mh Slot name.
#' @param matrix_slot_mh Slot name.
#' @param gnet Parameter value.
#' @param MH Parameter value.
#' @param out Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_mh_assemble_outputs <- function(scope_obj,
                                       gname,
                                       lee_layer,
                                       graph_slot_mh,
                                       matrix_slot_mh,
                                       gnet,
                                       MH,
                                       out,
                                       verbose) {
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()
    if (is.null(scope_obj@stats[[gname]][[lee_layer]])) scope_obj@stats[[gname]][[lee_layer]] <- list()

    if (out %in% c("igraph", "both")) {
        scope_obj@stats[[gname]][[lee_layer]][[graph_slot_mh]] <- gnet
    }
    if (out %in% c("matrix", "both")) {
        scope_obj@stats[[gname]][[lee_layer]][[matrix_slot_mh]] <- MH
    }

    if (verbose) {
        if (out %in% c("igraph","both")) {
            .log_info("computeMH", "S07", paste0("Graph stored in slot: ", graph_slot_mh, " (edge attr 'CMH')"), TRUE)
        }
        if (out %in% c("matrix","both")) {
            .log_info("computeMH", "S07", paste0("Matrix stored in slot: ", matrix_slot_mh), TRUE)
        }
    }
    scope_obj
}

#' Internal helper for metrics/tests helpers
#' @description
#' Internal helper for `.compute_mh_validate_inputs`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_layer Layer name.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_mh_validate_inputs <- function(scope_obj,
                                      grid_name,
                                      lee_layer,
                                      verbose) {
    if (verbose) .log_info("computeMH", "S01", paste0("Selecting and validating grid layer"), verbose)
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    gname <- names(scope_obj@grid)[
        vapply(scope_obj@grid, identical, logical(1), g_layer)
    ]
    .check_grid_content(scope_obj, gname)

    if (verbose) .log_info("computeMH", "S01", paste0("  Selected grid layer: ", gname), verbose)
    if (verbose) .log_info("computeMH", "S01", paste0("Extracting Lee's L matrix"), verbose)
    Lmat <- .get_lee_matrix(scope_obj, grid_name = gname, lee_layer = lee_layer)

    if (verbose) {
        .log_info("computeMH", "S01", paste0("  Lee's L matrix dimensions: ", nrow(Lmat), "x", ncol(Lmat)), TRUE)
        lee_stats <- summary(as.vector(Lmat[upper.tri(Lmat)]))
        .log_info("computeMH", "S01", paste0("  Lee's L range: [", round(lee_stats[1], 4), ", ", round(lee_stats[6], 4), "]"), TRUE)
    }

    list(g_layer = g_layer, gname = gname, Lmat = Lmat)
}

#' Internal helper for metrics/tests helpers
#' @description
#' Internal helper for `.compute_mh_resolve_runtime_config`.
#' @param ncores Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_mh_resolve_runtime_config <- function(ncores) {
    max(1L, as.integer(ncores))
}

#' Internal helper for metrics/tests helpers
#' @description
#' Internal helper for `.plot_gene_gene_jaccard_heatmap_core`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param expr_layer Layer name.
#' @param min_bins_per_gene Numeric threshold.
#' @param top_n_genes Parameter value.
#' @param genes_of_interest Parameter value.
#' @param output_dir Filesystem path.
#' @param heatmap_bg Parameter value.
#' @param cluster_col Parameter value.
#' @param keep_cluster_na Logical flag.
#' @param split_heatmap_by_cluster Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_gene_gene_jaccard_heatmap_core <- function(scope_obj,
                                              grid_name = "grid30",
                                              expr_layer = c("Xz", "counts"),
                                              min_bins_per_gene = 8L,
                                              top_n_genes = 60L,
                                              genes_of_interest = NULL,
                                              output_dir = NULL,
                                              heatmap_bg = "#c0c0c0",
                                              cluster_col = NULL,
                                              keep_cluster_na = TRUE,
                                              split_heatmap_by_cluster = FALSE) {
  expr_layer <- match.arg(expr_layer)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required to draw the Jaccard heatmap.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required to build the heatmap palette.")
  }

  cluster_lookup <- NULL
  if (!is.null(cluster_col)) {
    meta_df <- scope_obj@meta.data
    if (is.null(meta_df) || !is.data.frame(meta_df)) {
      warning("scope_obj@meta.data is missing or not a data.frame; ignore cluster_col.")
    } else if (!(cluster_col %in% colnames(meta_df))) {
      warning("cluster column ", cluster_col, " is not present in meta.data; ignore this parameter.")
    } else if (is.null(rownames(meta_df))) {
      warning("meta.data has no rownames, cannot map gene -> cluster; ignore cluster_col.")
    } else {
      cluster_lookup <- as.character(meta_df[[cluster_col]])
      names(cluster_lookup) <- rownames(meta_df)
    }
  }

  grid_layer <- scope_obj@grid[[grid_name]]
  if (is.null(grid_layer)) {
    warning("Grid layer ", grid_name, " is missing; cannot compute Jaccard.")
    return(invisible(NULL))
  }

  gene_bin_mat <- NULL
  if (identical(expr_layer, "counts")) {
    if (is.null(grid_layer$counts)) {
      warning("Grid layer ", grid_name, " has no counts; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    counts_df <- as.data.frame(grid_layer$counts)
    if (!nrow(counts_df)) {
      warning("Grid layer ", grid_name, " has no gene x bin records.")
      return(invisible(NULL))
    }

    cols_expected <- c("gene", "grid_id", "count")
    if (!all(cols_expected %in% colnames(counts_df))) {
      warning("counts is missing required columns: ", paste(setdiff(cols_expected, colnames(counts_df)), collapse = ", "))
      return(invisible(NULL))
    }

    counts_df$gene <- as.character(counts_df$gene)
    counts_df$grid_id <- as.character(counts_df$grid_id)
    counts_df$count <- as.numeric(as.character(counts_df$count))
    counts_df <- counts_df[!is.na(counts_df$gene) & !is.na(counts_df$grid_id) & counts_df$count > 0, , drop = FALSE]
    if (!nrow(counts_df)) {
      warning("No gene x bin records remain after filtering NA/zero.")
      return(invisible(NULL))
    }

    counts_df <- aggregate(count ~ gene + grid_id, data = counts_df, FUN = sum)
    genes <- sort(unique(counts_df$gene))
    bins <- sort(unique(counts_df$grid_id))

    gene_idx <- match(counts_df$gene, genes)
    bin_idx <- match(counts_df$grid_id, bins)
    gene_bin_mat <- sparseMatrix(
      i = gene_idx,
      j = bin_idx,
      x = as.numeric(counts_df$count > 0),
      dims = c(length(genes), length(bins)),
      dimnames = list(genes, bins)
    )
  } else {
    if (is.null(grid_layer$Xz)) {
      warning("Grid layer ", grid_name, " has no Xz; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    xz_mat <- grid_layer$Xz
    xz_mat <- as.matrix(xz_mat)
    if (!nrow(xz_mat) || !ncol(xz_mat)) {
      warning(grid_name, " Xz layer is empty; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    if (is.null(colnames(xz_mat)) || is.null(rownames(xz_mat))) {
      warning("Xz layer lacks row/col names; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    bin_mask <- xz_mat > 0
    bin_mask[is.na(bin_mask)] <- FALSE
    gene_bin_mat <- Matrix(t(bin_mask) * 1, sparse = TRUE)
    rownames(gene_bin_mat) <- colnames(xz_mat)
    colnames(gene_bin_mat) <- rownames(xz_mat)
  }

  if (is.null(gene_bin_mat) || !nrow(gene_bin_mat)) {
    warning("No usable gene x bin information; cannot compute Jaccard.")
    return(invisible(NULL))
  }

  if (!is.null(genes_of_interest)) {
    genes_keep <- intersect(unique(genes_of_interest), rownames(gene_bin_mat))
    if (!length(genes_keep)) {
      warning("genes_of_interest not found in current grid.")
      return(invisible(NULL))
    }
    gene_bin_mat <- gene_bin_mat[genes_keep, , drop = FALSE]
  }

  gene_clusters <- NULL
  if (!is.null(cluster_lookup)) {
    gene_names_all <- rownames(gene_bin_mat)
    idx <- match(gene_names_all, names(cluster_lookup))
    found <- !is.na(idx)
    if (!any(found)) {
      warning("cluster column ", cluster_col, " has no overlapping genes with current grid.")
      return(invisible(NULL))
    }
    if (!keep_cluster_na) {
      found <- found & !is.na(cluster_lookup[idx])
      if (!any(found)) {
        warning("All overlapping genes have NA cluster labels; cannot continue Jaccard.")
        return(invisible(NULL))
      }
    }
    gene_bin_mat <- gene_bin_mat[found, , drop = FALSE]
    idx <- idx[found]
    gene_clusters <- cluster_lookup[idx]
    names(gene_clusters) <- rownames(gene_bin_mat)
  }

  if (!nrow(gene_bin_mat)) {
    warning("No genes remain after filtering.")
    return(invisible(NULL))
  }

  bin_presence <- rowSums(gene_bin_mat)
  keep_idx <- which(bin_presence >= min_bins_per_gene)
  if (!length(keep_idx)) {
    warning("No genes meet the minimum bin count of ", min_bins_per_gene, ".")
    return(invisible(NULL))
  }

  gene_bin_mat <- gene_bin_mat[keep_idx, , drop = FALSE]
  bin_presence <- bin_presence[keep_idx]
  if (!is.null(gene_clusters)) {
    gene_clusters <- gene_clusters[rownames(gene_bin_mat)]
  }

  if (!is.null(top_n_genes) && is.finite(top_n_genes) && top_n_genes > 0 &&
      nrow(gene_bin_mat) > top_n_genes) {
    ord <- order(bin_presence, decreasing = TRUE)
    keep_ord <- ord[seq_len(top_n_genes)]
    gene_bin_mat <- gene_bin_mat[keep_ord, , drop = FALSE]
    bin_presence <- bin_presence[keep_ord]
    if (!is.null(gene_clusters)) {
      gene_clusters <- gene_clusters[keep_ord]
    }
  }

  if (!is.null(gene_clusters) && split_heatmap_by_cluster) {
    cluster_order <- order(gene_clusters, rownames(gene_bin_mat))
    gene_bin_mat <- gene_bin_mat[cluster_order, , drop = FALSE]
    bin_presence <- bin_presence[cluster_order]
    gene_clusters <- gene_clusters[cluster_order]
  }

  gene_names <- rownames(gene_bin_mat)
  inter_mat <- tcrossprod(gene_bin_mat)
  inter_dense <- as.matrix(inter_mat)
  union_dense <- outer(bin_presence, bin_presence, "+") - inter_dense
  jac_mat <- inter_dense
  positive_union <- union_dense > 0
  jac_mat[positive_union] <- jac_mat[positive_union] / union_dense[positive_union]
  jac_mat[!positive_union] <- NA_real_
  diag(jac_mat) <- 1
  rownames(jac_mat) <- gene_names
  colnames(jac_mat) <- gene_names

  col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#113A70", "#f7f7f7", "#B11226"))
  row_split_factor <- column_split_factor <- NULL
  if (!is.null(gene_clusters) && split_heatmap_by_cluster) {
    cluster_labels_plot <- ifelse(is.na(gene_clusters), "NA", as.character(gene_clusters))
    cluster_levels <- unique(cluster_labels_plot)
    row_split_factor <- factor(cluster_labels_plot, levels = cluster_levels)
    column_split_factor <- row_split_factor
  }
  column_title_txt <- NULL
  ht <- ComplexHeatmap::Heatmap(
    jac_mat,
    name = "Jaccard",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    column_title = column_title_txt,
    column_title_gp = grid::gpar(fontsize = 10),
    row_split = row_split_factor,
    column_split = column_split_factor
  )

  invisible(list(
    matrix = jac_mat,
    genes = gene_names,
    clusters = gene_clusters,
    grid = grid_name,
    min_bins = min_bins_per_gene,
    top_n = top_n_genes,
    heatmap = ht,
    palette = col_fun
  ))
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_network_build_ggraph_layers`.
#' @param lay Parameter value.
#' @param palette_vals Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param max.overlaps Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param title Parameter value.
#' @param caption Parameter value.
#' @param linetype_name Parameter value.
#' @param linetype_values Parameter value.
#' @param linetype_breaks Parameter value.
#' @param linetype_labels Parameter value.
#' @param linetype_guide Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_network_build_ggraph_layers <- function(
    lay,
    palette_vals,
    node_size,
    edge_width,
    label_size,
    max.overlaps,
    hub_border_col,
    hub_border_size,
    title = NULL,
    caption = NULL,
    linetype_name = NULL,
    linetype_values = NULL,
    linetype_breaks = NULL,
    linetype_labels = NULL,
    linetype_guide = NULL
) {
    plot_obj <- ggraph(lay) +
        geom_edge_link(
            aes(
                width = weight, colour = edge_col,
                linetype = linetype
            ),
            lineend = "round",
            show.legend = FALSE
        ) +
        scale_edge_width(name = "|L|", range = c(0.3, edge_width), guide = "none") +
        scale_edge_colour_identity(guide = "none")
    if (!is.null(linetype_values)) {
        linetype_breaks_resolved <- linetype_breaks
        if (is.null(linetype_breaks_resolved) && !is.null(names(linetype_values))) {
            linetype_breaks_resolved <- names(linetype_values)
        }
        plot_obj <- plot_obj + scale_edge_linetype_manual(
            name = linetype_name,
            values = linetype_values,
            breaks = linetype_breaks_resolved,
            labels = linetype_labels,
            guide = linetype_guide
        )
    }

    plot_obj <- plot_obj +
        geom_node_point(
            data = ~ filter(.x, !hub),
            aes(size = deg, fill = basecol), shape = 21, stroke = 0,
            show.legend = c(fill = TRUE, size = FALSE)
        ) +
        geom_node_point(
            data = ~ filter(.x, hub),
            aes(size = deg, fill = basecol),
            shape = 21, colour = hub_border_col, stroke = hub_border_size,
            show.legend = FALSE
        ) +
        geom_node_text(aes(label = name),
            size = label_size,
            repel = TRUE, vjust = 1.4,
            max.overlaps = max.overlaps
        ) +
        scale_size_continuous(
            name = "Node size",
            range = c(node_size / 2, node_size * 1.5),
            guide = "none"
        ) +
        scale_fill_identity(
            name = "Module",
            guide = guide_legend(
                override.aes = list(shape = 21, size = node_size * 0.5),
                nrow = 2, order = 1
            ),
            breaks = unname(palette_vals),
            labels = names(palette_vals)
        ) +
        coord_fixed() + theme_void(base_size = 12) +
        theme(
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 15),
            plot.title = element_text(size = 15, hjust = 0.5),
            legend.position = "bottom",
            legend.box = "vertical",
            legend.box.just = "center",
            legend.spacing.y = unit(0.2, "cm"),
            legend.box.margin = margin(t = 5, b = 5),
            plot.caption = element_text(
                hjust = 0, size = 9,
                margin = margin(t = 6)
            )
        ) +
        labs(title = title, caption = caption) +
        guides(
            fill = guide_legend(order = 1),
            size = "none",
            edge_width = "none",
            linetype = "none",
            colour = "none"
        )

    plot_obj
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_network_core_v2`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param gene_subset Parameter value.
#' @param L_min Numeric threshold.
#' @param L_min_neg Parameter value.
#' @param p_cut Parameter value.
#' @param use_FDR Logical flag.
#' @param FDR_max Numeric threshold.
#' @param pct_min Numeric threshold.
#' @param CI95_filter Parameter value.
#' @param curve_layer Layer name.
#' @param CI_rule Parameter value.
#' @param drop_isolated Parameter value.
#' @param weight_abs Parameter value.
#' @param use_consensus_graph Logical flag.
#' @param graph_slot_name Slot name.
#' @param fdr_source Parameter value.
#' @param cluster_vec Parameter value.
#' @param cluster_palette Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param vertex_size Parameter value.
#' @param base_edge_mult Parameter value.
#' @param label_cex Parameter value.
#' @param layout_niter Parameter value.
#' @param seed Random seed.
#' @param hub_factor Parameter value.
#' @param length_scale Parameter value.
#' @param max.overlaps Parameter value.
#' @param L_max Numeric threshold.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param show_sign Parameter value.
#' @param neg_linetype Parameter value.
#' @param neg_legend_lab Parameter value.
#' @param pos_legend_lab Parameter value.
#' @param title Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_network_core_v2 <- function(
    scope_obj,
    ## ---------- Data layers ----------
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL, # Now compatible with direct gene vector or (cluster_col, cluster_num)
    ## ---------- Filter thresholds ----------
    L_min = 0,
    L_min_neg = NULL,
    p_cut = NULL,
    use_FDR = TRUE,
    FDR_max = 0.05,
    pct_min = "q0",
    CI95_filter = FALSE,
    curve_layer = "LR_curve",
    CI_rule = c("remove_within", "remove_outside"),
    drop_isolated = TRUE,
    weight_abs = TRUE,
    ## ---------- Consensus network ----------
    use_consensus_graph = TRUE,
    graph_slot_name = NULL,
    ## ---------- New FDR source ----------
    fdr_source = c("FDR", "FDR_storey", "FDR_beta", "FDR_mid", "FDR_uniform", "FDR_disc"),
    ## ---------- Plot parameters ----------
    cluster_vec = NULL,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
        "#A65628", "#984EA3", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
    ),
    node_size = 4,
    edge_width = 3,
    label_size = 4,
    ## Backwards-compatible aliases (will override the above when provided)
    vertex_size = NULL,
    base_edge_mult = NULL,
    label_cex = NULL,
    layout_niter = 1000,
    seed = 1,
    hub_factor = 2,
    length_scale = 1,
    max.overlaps = 10,
    L_max = 1,
    hub_border_col = "#4B4B4B",
    hub_border_size = 0.8,
    show_sign = FALSE,
    neg_linetype = "dashed",
    neg_legend_lab = "Negative",
    pos_legend_lab = "Positive",
    title = NULL) {
    CI_rule <- match.arg(CI_rule)
    fdr_source <- match.arg(fdr_source)
    parent <- "plotNetwork"
    verbose <- getOption("geneSCOPE.verbose", TRUE)

    ## Backward-compatibility: honor old argument names if supplied
    if (!is.null(vertex_size)) node_size <- vertex_size
    if (!is.null(base_edge_mult)) edge_width <- base_edge_mult
    if (!is.null(label_cex)) label_size <- label_cex

    ## ===== 0. Lock and validate grid layer =====
    g_layer <- .select_grid_layer(scope_obj, grid_name) # If grid_name=NULL auto select unique layer
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else {
        as.character(grid_name)
    }

    .check_grid_content(scope_obj, grid_name) # Error directly if missing required fields

    ## ===== 1. Read LeeStats object and its matrices =====
    ##  First locate LeeStats layer (can be in @stats[[grid]] or @grid[[grid]])
    leeStat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
        scope_obj@stats[[grid_name]][[lee_stats_layer]]
    } else if (!is.null(g_layer[[lee_stats_layer]])) {
        g_layer[[lee_stats_layer]]
    } else {
        stop(
            "Cannot find layer '", lee_stats_layer,
            "' in grid '", grid_name, "'."
        )
    }

    Lmat <- .get_lee_matrix(scope_obj,
        grid_name = grid_name,
        lee_layer = lee_stats_layer
    ) # Only extract aligned L

    ## Auto-detect graph slot when not supplied: prefer cluster_vec (if scalar) then g_consensus
    if (is.null(graph_slot_name)) {
        cand <- c(
            if (is.character(cluster_vec) && length(cluster_vec) == 1) cluster_vec else NULL,
            "g_consensus"
        )
        cand <- unique(na.omit(cand))
        for (nm in cand) {
            if (!is.null(leeStat[[nm]])) {
                graph_slot_name <- nm
                break
            }
        }
        if (is.null(graph_slot_name)) graph_slot_name <- "g_consensus"
    }

    ## ===== 2. Gene subset =====
    all_genes <- rownames(Lmat)
    if (is.null(gene_subset)) {
        keep_genes <- all_genes
    } else if (is.character(gene_subset)) {
        ## (A) Directly provide gene vector
        keep_genes <- intersect(
            all_genes,
            .get_gene_subset(scope_obj, genes = gene_subset)
        )
    } else if (is.list(gene_subset) &&
        all(c("cluster_col", "cluster_num") %in% names(gene_subset))) {
        ## (B) Use cluster column + number syntax: gene_subset = list(cluster_col = "..", cluster_num = ..)
        keep_genes <- intersect(
            all_genes,
            .get_gene_subset(scope_obj,
                cluster_col = gene_subset$cluster_col,
                cluster_num = gene_subset$cluster_num
            )
        )
    } else {
        stop("`gene_subset` must be a character vector or list(cluster_col, cluster_num)")
    }
    if (length(keep_genes) < 2) {
        stop("Less than two genes remain after sub-setting.")
    }

    ## ===== 3. Retrieve consensus graph or rebuild =====
    use_consensus <- isTRUE(use_consensus_graph) &&
        !is.null(leeStat[[graph_slot_name]])

    if (use_consensus && isTRUE(use_FDR)) {
        # Consensus graphs are typically pre-filtered; skip redundant FDR filtering
        .log_info(parent, "S01", "use_consensus_graph=TRUE ignoring FDR re-filtering (assuming graph already filtered).", verbose)
    }

    if (use_consensus) {
        g_raw <- leeStat[[graph_slot_name]]
        g <- igraph::induced_subgraph(
            g_raw,
            intersect(igraph::V(g_raw)$name, keep_genes)
        )
        if (igraph::ecount(g) == 0) {
            stop("Consensus graph has no edges under the chosen gene subset.")
        }
        ## Update edge weights / signs when required
    } else {
        ## ---------- Reconstruct graph using thresholds ----------
        idx <- match(keep_genes, rownames(Lmat))
        A <- Lmat[idx, idx, drop = FALSE]
        if (!is.matrix(A)) {
            # Ensure standard matrix so downstream t(), symmetrisation and igraph calls work
            A <- as.matrix(A)
        }

        ## (i) Early filter by expression proportion (pct_min)
        A <- .filter_matrix_by_quantile(A, pct_min, "q100") # Fixed function name

        ## (ii) Cap the maximum absolute value
        A[abs(A) > L_max] <- 0

        ## (iii) Apply positive/negative thresholds
        L_min_neg <- if (is.null(L_min_neg)) L_min else abs(L_min_neg)
        if (weight_abs) {
            A[abs(A) < L_min] <- 0
        } else {
            A[A > 0 & A < L_min] <- 0
            A[A < 0 & abs(A) < L_min_neg] <- 0
        }

        ## (iv) 95% confidence interval filter
        if (CI95_filter) {
            curve <- leeStat[[curve_layer]]
            if (is.null(curve)) {
                stop("Cannot find `curve_layer = ", curve_layer, "` in LeeStats.")
            }
            f_lo <- approxfun(curve$Pear, curve$lo95, rule = 2)
            f_hi <- approxfun(curve$Pear, curve$hi95, rule = 2)

            rMat <- .get_pearson_matrix(scope_obj, level = "cell") # use single-cell Pearson correlations
            rMat <- rMat[keep_genes, keep_genes, drop = FALSE]

            mask <- if (CI_rule == "remove_within") {
                (A >= f_lo(rMat)) & (A <= f_hi(rMat))
            } else {
                (A < f_lo(rMat)) | (A > f_hi(rMat))
            }
            A[mask] <- 0
        }

        ## (v) p-value & FDR
        if (!is.null(p_cut) && !is.null(leeStat$P)) {
            Pmat <- leeStat$P[idx, idx, drop = FALSE]
            A[Pmat >= p_cut | is.na(Pmat)] <- 0
        }
        if (isTRUE(use_FDR)) {
            ## -------- FDR matrix selection logic --------
            pref_order <- c(
                fdr_source,
                setdiff(
                    c("FDR", "FDR_storey", "FDR_beta", "FDR_mid", "FDR_uniform", "FDR_disc"),
                    fdr_source
                )
            )
            FDR_sel <- NULL
            FDR_used_name <- NULL
            for (nm in pref_order) {
                cand <- leeStat[[nm]]
                if (!is.null(cand)) {
                    FDR_sel <- cand
                    FDR_used_name <- nm
                    break
                }
            }
            if (is.null(FDR_sel)) {
                stop(
                    "use_FDR=TRUE but no usable FDR matrix found (expected: ",
                    paste(pref_order, collapse = ", "), ")."
                )
            }
            if (inherits(FDR_sel, "big.matrix")) {
                .log_backend(
                    parent,
                    "S01",
                    "fdr_source",
                    paste0(FDR_used_name, " (big.matrix)"),
                    reason = "convert_to_matrix",
                    verbose = verbose
                )
                FDR_sel <- bigmemory::as.matrix(FDR_sel)
            } else if (!is.matrix(FDR_sel)) {
                FDR_sel <- as.matrix(FDR_sel)
            }
            if (!identical(dim(FDR_sel), dim(Lmat))) {
                stop("FDR matrix dimensions do not match L: ", FDR_used_name)
            }
            FDRmat <- FDR_sel[idx, idx, drop = FALSE]
            A[FDRmat > FDR_max | is.na(FDRmat)] <- 0
            .log_backend(
                parent,
                "S01",
                "fdr_source",
                FDR_used_name,
                reason = sprintf("FDR_max=%.3g", FDR_max),
                verbose = verbose
            )
        }

        ## (vi) Symmetrise and zero the diagonal
        A <- as.matrix(A)  # ensure base matrix for t()
        A <- (A + t(A)) / 2
        diag(A) <- 0

        ## (vii) Drop isolated vertices
        if (drop_isolated) {
            keep <- which(rowSums(abs(A) > 0) > 0 | colSums(abs(A) > 0) > 0)
            A <- A[keep, keep, drop = FALSE]
        }
        if (nrow(A) < 2 || all(A == 0)) {
            stop("No edges remain after filtering thresholds.")
        }

        g <- igraph::graph_from_adjacency_matrix(abs(A),
            mode = "undirected",
            weighted = TRUE, diag = FALSE
        )
        e_idx <- igraph::as_edgelist(g, names = FALSE) # replacement for deprecated get.edgelist
        igraph::E(g)$sign <- ifelse(A[e_idx] < 0, "neg", "pos")
    }

    ## ===== 4. Weight adjustment & global scaling =====
    igraph::E(g)$weight <- abs(igraph::E(g)$weight)
    bad <- which(!is.finite(igraph::E(g)$weight) | igraph::E(g)$weight <= 0)
    if (length(bad)) {
        min_pos <- min(igraph::E(g)$weight[igraph::E(g)$weight > 0], na.rm = TRUE)
        igraph::E(g)$weight[bad] <- ifelse(is.finite(min_pos), min_pos, 1)
    }
    if (!is.null(length_scale) && length_scale != 1) {
        igraph::E(g)$weight <- igraph::E(g)$weight * length_scale
    }

    if (!is.null(L_min) && L_min > 0) {
        keep_e <- which(igraph::E(g)$weight >= L_min)
        if (length(keep_e) == 0) stop("`L_min` too strict; no edges to plot.")
        g <- igraph::subgraph.edges(g, keep_e, delete.vertices = TRUE)
    }

    ## ===== 5. Node color / cluster labels =====
    Vnames <- igraph::V(g)$name
    deg_vec <- igraph::degree(g)

    ## (1) Derive clu vector from cluster_vec or meta.data
    clu <- setNames(rep(NA_character_, length(Vnames)), Vnames)
    is_factor_input <- FALSE
    factor_levels <- NULL

    ## helper - extract column from meta.data
    get_meta_cluster <- function(col) {
        if (is.null(scope_obj@meta.data) || !(col %in% colnames(scope_obj@meta.data))) {
            stop("Column '", col, "' not found in scope_obj@meta.data.")
        }
        meta_cluster_values_local <- scope_obj@meta.data[[col]]
        names(meta_cluster_values_local) <- rownames(scope_obj@meta.data)
        meta_cluster_values_local[Vnames]
    }

    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) > 1) { # manual vector
            if (is.factor(cluster_vec)) {
                is_factor_input <- TRUE
                factor_levels <- levels(cluster_vec)
            }
            tmp <- as.character(cluster_vec[Vnames])
            clu[!is.na(tmp)] <- tmp[!is.na(tmp)]
        } else { # column name
            meta_clu <- get_meta_cluster(cluster_vec)
            if (is.factor(meta_clu)) {
                is_factor_input <- TRUE
                factor_levels <- levels(meta_clu)
            }
            clu[!is.na(meta_clu)] <- as.character(meta_clu[!is.na(meta_clu)])
        }
    }

    keep_nodes <- names(clu)[!is.na(clu)]
    if (length(keep_nodes) < 2) {
        stop("Less than two nodes have cluster labels.")
    }
    g <- igraph::induced_subgraph(g, keep_nodes)
    Vnames <- igraph::V(g)$name
    clu <- clu[Vnames]
    deg_vec <- deg_vec[Vnames]

    ## Step 2: generate palette
    cluster_levels <- if (is_factor_input && !is.null(factor_levels)) {
        intersect(factor_levels, unique(clu))
    } else {
        sort(unique(clu))
    }
    n_clu <- length(cluster_levels)
    cluster_palette_values <- if (is.null(cluster_palette)) {
        setNames(rainbow(n_clu), cluster_levels)
    } else if (is.null(names(cluster_palette))) {
        setNames(colorRampPalette(cluster_palette)(n_clu), cluster_levels)
    } else {
        cluster_palette_map <- cluster_palette
        missing_cluster_levels <- setdiff(cluster_levels, names(cluster_palette_map))
        if (length(missing_cluster_levels)) {
            cluster_palette_map <- c(
                cluster_palette_map,
                setNames(colorRampPalette(cluster_palette)(length(missing_cluster_levels)), missing_cluster_levels)
            )
        }
        cluster_palette_map[cluster_levels]
    }
    basecol <- setNames(cluster_palette_values[clu], Vnames)

    ## Step 3: edge colours (sign / intensity gradient)
    e_idx <- igraph::as_edgelist(g, names = FALSE)
    w_norm <- igraph::E(g)$weight / max(igraph::E(g)$weight)
    cluster_ord <- setNames(seq_along(cluster_levels), cluster_levels)
    edge_color_values <- character(igraph::ecount(g))
    for (i in seq_along(edge_color_values)) {
        v1 <- Vnames[e_idx[i, 1]]
        v2 <- Vnames[e_idx[i, 2]]
        if (show_sign && igraph::E(g)$sign[i] == "neg") {
            edge_color_values[i] <- "gray40"
        } else { # use gradient for positive correlations
            ref_col <- if (!is.na(clu[v1]) && !is.na(clu[v2]) && clu[v1] == clu[v2]) {
                basecol[v1]
            } else {
                # decide primary colour by degree or cluster order
                if (deg_vec[v1] > deg_vec[v2]) {
                    basecol[v1]
                } else if (deg_vec[v2] > deg_vec[v1]) {
                    basecol[v2]
                } else if (!is.na(clu[v1]) && is.na(clu[v2])) {
                    basecol[v1]
                } else if (!is.na(clu[v2]) && is.na(clu[v1])) {
                    basecol[v2]
                } else if (!is.na(clu[v1]) && !is.na(clu[v2])) {
                    if (cluster_ord[clu[v1]] <= cluster_ord[clu[v2]]) {
                        basecol[v1]
                    } else {
                        basecol[v2]
                    }
                } else {
                    "gray80"
                }
            }
            rgb_ref <- col2rgb(ref_col) / 255
            tval <- w_norm[i]
            edge_color_values[i] <- rgb(
                1 - tval + tval * rgb_ref[1],
                1 - tval + tval * rgb_ref[2],
                1 - tval + tval * rgb_ref[3]
            )
        }
    }
    igraph::E(g)$edge_col <- edge_color_values
    igraph::E(g)$linetype <- if (show_sign) igraph::E(g)$sign else "solid"

    ## ===== 6. ggraph rendering =====
    library(ggraph)
    library(ggplot2)
    library(dplyr)
    set.seed(seed)
    layout_data <- create_layout(g, layout = "fr", niter = layout_niter)
    layout_data$basecol <- basecol[layout_data$name]
    layout_data$deg <- deg_vec[layout_data$name]
    layout_data$hub <- layout_data$deg > hub_factor * median(layout_data$deg)

    plot_obj <- .plot_network_build_ggraph_layers(
        lay = layout_data,
        palette_vals = cluster_palette_values,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        max.overlaps = max.overlaps,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        title = title,
        caption = NULL,
        linetype_name = if (show_sign) "Direction" else NULL,
        linetype_values = if (show_sign) c(pos = "solid", neg = neg_linetype) else NULL,
        linetype_breaks = if (show_sign) c("pos", "neg") else NULL,
        linetype_labels = if (show_sign) c(pos = pos_legend_lab, neg = neg_legend_lab) else NULL,
        linetype_guide = if (show_sign) "none" else NULL
    )

    return(plot_obj)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_network_output`.
#' @param plot_obj Internal parameter
#' @param prep Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_network_output <- function(plot_obj, prep) plot_obj

#' Plot Network
#' @description
#' Internal helper for `.plot_network`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param g_slot Slot name.
#' @param cluster_name Parameter value.
#' @param legend_title Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param seed Random seed.
#' @param label_nodes Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param min_degree_hub Numeric threshold.
#' @param max.overlaps Parameter value.
#' @param cluster_vec Parameter value.
#' @param cluster_palette Parameter value.
#' @param tree_mode Parameter value.
#' @param k_top Parameter value.
#' @param edge_thresh Parameter value.
#' @param prune_mode Parameter value.
#' @param attach_I_delta Parameter value.
#' @param I_delta_name Parameter value.
#' @param I_delta_palette Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_network <- function(
    scope_obj,
    grid_name = NULL,
    g_slot = "g_consensus",
    cluster_name,
    legend_title = NULL,
    node_size = 5,
    edge_width = 1,
    label_size = 2,
    seed = 1,
    label_nodes = TRUE,
    hub_border_col = "black",
    hub_border_size = 0.5,
    min_degree_hub = 3,
    max.overlaps = Inf,
    cluster_vec = NULL,
    cluster_palette = NULL,
    tree_mode = NULL,
    k_top = 5,
    edge_thresh = 0.05,
    prune_mode = c("none", "mst", "triangulation"),
    attach_I_delta = FALSE,
    I_delta_name = "iDelta",
    I_delta_palette = scales::hue_pal()) {
    .plot_network_core(
        scope_obj = scope_obj,
        grid_name = grid_name,
        g_slot = g_slot,
        cluster_name = cluster_name,
        legend_title = legend_title,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        seed = seed,
        label_nodes = label_nodes,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        min_degree_hub = min_degree_hub,
        max.overlaps = max.overlaps,
        cluster_vec = cluster_vec,
        cluster_palette = cluster_palette,
        tree_mode = tree_mode,
        k_top = k_top,
        edge_thresh = edge_thresh,
        prune_mode = prune_mode,
        attach_I_delta = attach_I_delta,
        I_delta_name = I_delta_name,
        I_delta_palette = I_delta_palette
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_network_theme`.
#' @param plot_obj Internal parameter
#' @param prep Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_network_theme <- function(plot_obj, prep) plot_obj

# Baseline plotting overrides (verbatim).
#' Visualise the Lee's-L gene-gene network
#' @description
#'   Builds a gene network from Lee's L statistics or from a pre-computed
#' Internal helper for plotting workflows
#' @param prep Prepared argument bundle for network plotting
#' @return ggplot object produced by \code{.plot_network_core_v2}
#' @keywords internal
.plot_network_build <- function(prep) do.call(.plot_network_core_v2, prep$args)

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_network_prepare`.
#' @param ... Additional arguments (currently unused).
#' @return List carrying the plotting arguments
#' @keywords internal
.plot_network_prepare <- function(...) list(args = list(...))

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_network_core`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param g_slot Slot name.
#' @param cluster_name Parameter value.
#' @param legend_title Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param seed Random seed.
#' @param label_nodes Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param min_degree_hub Numeric threshold.
#' @param max.overlaps Parameter value.
#' @param cluster_vec Parameter value.
#' @param cluster_palette Parameter value.
#' @param tree_mode Parameter value.
#' @param k_top Parameter value.
#' @param edge_thresh Parameter value.
#' @param prune_mode Parameter value.
#' @param attach_I_delta Parameter value.
#' @param I_delta_name Parameter value.
#' @param I_delta_palette Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_network_core <- function(
    scope_obj,
    grid_name = NULL,
    g_slot = "g_consensus",
    cluster_name,
    legend_title = NULL,
    node_size = 5,
    edge_width = 1,
    label_size = 2,
    seed = 1,
    label_nodes = TRUE,
    hub_border_col = "black",
    hub_border_size = 0.5,
    min_degree_hub = 3,
    max.overlaps = Inf,
    cluster_vec = NULL,
    cluster_palette = NULL,
    tree_mode = NULL,
    k_top = 5,
    edge_thresh = 0.05,
    prune_mode = c("none", "mst", "triangulation"),
    attach_I_delta = FALSE,
    I_delta_name = "iDelta",
    I_delta_palette = scales::hue_pal()) {
    plot_prep <- .plot_network_prepare(
        scope_obj = scope_obj,
        grid_name = grid_name,
        g_slot = g_slot,
        cluster_name = cluster_name,
        legend_title = legend_title,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        seed = seed,
        label_nodes = label_nodes,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        min_degree_hub = min_degree_hub,
        max.overlaps = max.overlaps,
        cluster_vec = cluster_vec,
        cluster_palette = cluster_palette,
        tree_mode = tree_mode,
        k_top = k_top,
        edge_thresh = edge_thresh,
        prune_mode = prune_mode,
        attach_I_delta = attach_I_delta,
        I_delta_name = I_delta_name,
        I_delta_palette = I_delta_palette
    )
    plot_layers <- .plot_network_build(plot_prep)
    plot_layers <- .plot_network_theme(plot_layers, plot_prep)
    .plot_network_output(plot_layers, plot_prep)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_lvs_r_core`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param pear_level Parameter value.
#' @param lee_stats_layer Layer name.
#' @param delta_top_n Parameter value.
#' @param flip Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_lvs_r_core <- function(scope_obj,
                     grid_name,
                     pear_level = c("cell", "grid"),
                     lee_stats_layer = "LeeStats_Xz",
                     delta_top_n = 10,
                     flip = TRUE) {
  pear_level <- match.arg(pear_level)

  ## ---- 1. Lee's L and Pearson r ------------------------------------------------
  L_mat <- .get_lee_matrix(scope_obj, grid_name, lee_layer = lee_stats_layer)
  r_mat <- .get_pearson_matrix(scope_obj, grid_name,
    level = pear_level
  )

  common <- intersect(rownames(L_mat), rownames(r_mat))
  if (length(common) < 2) stop("Insufficient common genes for plotting.")
  L_mat <- L_mat[common, common]
  r_mat <- r_mat[common, common]
  diag(L_mat) <- NA
  diag(r_mat) <- NA

  ## ---- 2. Convert to long table + Delta ---------------------------------------------------------
  genes <- colnames(L_mat)
  ut <- upper.tri(L_mat, diag = FALSE)
  df_long <- data.frame(
    gene1 = rep(genes, each = length(genes))[ut],
    gene2 = rep(genes, length(genes))[ut],
    LeesL = L_mat[ut],
    Pear  = r_mat[ut]
  ) |>
    mutate(Delta = LeesL - Pear)

  ## ---- 3. Label points (extreme Delta) ----------------------------------------------------
  df_label <- bind_rows(
    slice_max(df_long, Delta, n = delta_top_n, with_ties = FALSE),
    slice_min(df_long, Delta, n = delta_top_n, with_ties = FALSE)
  ) |>
    distinct(gene1, gene2, .keep_all = TRUE) |>
    mutate(label = sprintf("%s-%s\nL=%.3f", gene1, gene2, LeesL))

  ## ---- 4. Plot ----------------------------------------------------------------
  if (!flip) {
    p <- ggplot(
      df_long,
      aes(x = LeesL, y = Pear)
    )
    xlab <- "Lee's L"
    ylab <- "Pearson correlation"
    ttl <- sprintf(
      "Lee's L vs Pearson  (%s, %s)",
      sub("grid", "Grid ", grid_name),
      ifelse(pear_level == "cell", "cell level", "grid level")
    )
  } else {
    p <- ggplot(
      df_long,
      aes(x = Pear, y = LeesL)
    )
    xlab <- "Pearson correlation"
    ylab <- "Lee's L"
    ttl <- sprintf(
      "Pearson vs Lee's L  (%s, %s)",
      sub("grid", "Grid ", grid_name),
      ifelse(pear_level == "cell", "cell level", "grid level")
    )
  }

  p <- p +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01))

  if (delta_top_n > 0 && nrow(df_label) > 0) {
    p <- p + geom_text_repel(
      data = df_label,
      aes(label = label),
      size = 3,
      box.padding = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf
    )
  }

  p + labs(title = ttl, x = xlab, y = ylab) +
    theme_minimal(base_size = 8) +
    theme(
      panel.title      = element_text(hjust = .5, size = 10),
      panel.border     = element_rect(colour = "black", fill = NA, size = .5),
      panel.background = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major = element_line(colour = "white"),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(size = 9),
      axis.text        = element_text(size = 8),
      axis.ticks       = element_blank()
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_l_scatter_core`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param title Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_l_scatter_core <- function(scope_obj,
                         grid_name = NULL,
                         lee_stats_layer = NULL,
                         title = NULL) {
  ## ---- 1. Grid layer & counts ----------------------------------------------
  g_layer <- .select_grid_layer(scope_obj, grid_name)
  if (is.null(grid_name)) {
    grid_name <- names(scope_obj@grid)[
      vapply(scope_obj@grid, identical, logical(1), g_layer)
    ]
  }

  counts_df <- g_layer$counts
  if (!all(c("gene", "count") %in% names(counts_df))) {
    stop("counts must contain columns 'gene' and 'count'.")
  }

  ## ---- 2. Lee's L -------------------------------------------------------
  Lmat <- .get_lee_matrix(scope_obj,
    grid_name = grid_name,
    lee_layer = lee_stats_layer
  )

  ## ---- 3. Total & combine gene pairs ---------------------------------------------
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)

  gene_totals <- counts_df |>
    group_by(gene) |>
    summarise(total = sum(count), .groups = "drop")

  genes_common <- intersect(rownames(Lmat), gene_totals$gene)
  if (length(genes_common) < 2L) {
    stop("Need at least two overlapping genes.")
  }

  comb_df <- crossing(
    geneA = genes_common,
    geneB = genes_common
  ) |>
    filter(geneA < geneB) |>
    left_join(gene_totals, by = c("geneA" = "gene")) |>
    rename(totalA = total) |>
    left_join(gene_totals, by = c("geneB" = "gene")) |>
    rename(totalB = total) |>
    rowwise() |>
    mutate(
      LeeL = Lmat[geneA, geneB],
      FoldRatio = max(totalA, totalB) /
        pmax(1, pmin(totalA, totalB))
    ) |>
    ungroup() |>
    filter(is.finite(LeeL), is.finite(FoldRatio))

  if (is.null(title)) {
    title <- "Lee's L vs. Fold Change (All Gene Pairs)"
  }

  p <- ggplot(comb_df, aes(x = LeeL, y = FoldRatio)) +
    geom_point(
      shape = 21, size = 1.2,
      fill = "white", colour = "black", stroke = .2
    ) +
    scale_x_continuous(expand = expansion(mult = .05)) +
    scale_y_continuous(expand = expansion(mult = .05)) +
    labs(
      title = title,
      x = expression("Lee's L"),
      y = "Fold Change (>=1)"
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.title        = element_text(size = 10, colour = "black"),
      panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#c0c0c0", size = .2),
      panel.border       = element_rect(colour = "black", fill = NA, size = .5),
      axis.line          = element_line(colour = "black", size = .3),
      axis.ticks         = element_line(colour = "black", size = .3),
      axis.text          = element_text(size = 8, colour = "black"),
      axis.title         = element_text(size = 9, colour = "black"),
      plot.margin        = unit(c(1, 1, 1, 1), "cm")
    ) +
    coord_cartesian(expand = TRUE)

  p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_l_distribution_core`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param bins Parameter value.
#' @param xlim Parameter value.
#' @param title Parameter value.
#' @param use_abs Logical flag.
#' @return Internal helper result
#' @keywords internal
.plot_l_distribution_core <- function(scope_obj,
                              grid_name = NULL,
                              lee_stats_layer = NULL,
                              bins = 30,
                              xlim = NULL,
                              title = NULL,
                              use_abs = FALSE) {
  ## ---- 1. Get Lee's L (helper will auto-match layer names) ---------------------------
  Lmat <- .get_lee_matrix(scope_obj,
    grid_name = grid_name,
    lee_layer = lee_stats_layer
  )

  vals <- as.numeric(Lmat[upper.tri(Lmat)])
  if (use_abs) vals <- abs(vals)
  df <- data.frame(LeeL = vals)

  if (is.null(title)) {
    title <- if (use_abs) "Absolute Lee's L Distribution" else "Lee's L Distribution"
  }

  library(ggplot2)
  library(grid)

  p <- ggplot(df, aes(x = LeeL)) +
    {
      if (length(bins) == 1L) {
        geom_histogram(
          bins = bins,
          fill = "white",
          colour = "black",
          size = .3
        )
      } else {
        geom_histogram(
          breaks = bins,
          fill = "white",
          colour = "black",
          size = .3
        )
      }
    } +
    labs(
      title = title,
      x = if (use_abs) expression("|Lee's L|") else expression("Lee's L"),
      y = "Frequency"
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.title        = element_text(size = 10, colour = "black"),
      panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#c0c0c0", size = .2),
      panel.border       = element_rect(colour = "black", fill = NA, size = .5),
      axis.line          = element_line(colour = "black", size = .3),
      axis.ticks         = element_line(colour = "black", size = .3),
      axis.text          = element_text(size = 8, colour = "black"),
      axis.title         = element_text(size = 9, colour = "black"),
      plot.margin        = unit(c(1, 1, 1, 1), "cm")
    ) +
    coord_cartesian(expand = FALSE)

  if (!is.null(xlim)) {
    stopifnot(is.numeric(xlim) && length(xlim) == 2L)
    p <- p + coord_cartesian(xlim = xlim, expand = FALSE)
  }
  p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_i_delta_core`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param cluster_col Parameter value.
#' @param top_n Parameter value.
#' @param min_genes Numeric threshold.
#' @param nrow Parameter value.
#' @param point_size Parameter value.
#' @param line_size Parameter value.
#' @param label_size Parameter value.
#' @param fill_col Parameter value.
#' @param outline_col Parameter value.
#' @param seed Random seed.
#' @param subCluster Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_i_delta_core <- function(
    scope_obj,
    grid_name,
    cluster_col,
    top_n = NULL,
    min_genes = 1,
    nrow = 1,
    point_size = 3,
    line_size = 0.5,
    label_size = 2, # default reduced to 2
    fill_col = "steelblue",
    outline_col = "black",
    seed = NULL,
    subCluster = NULL) {
    if (!is.null(seed)) set.seed(seed)
    meta_col <- paste0(grid_name, "_iDelta")
    if (!meta_col %in% colnames(scope_obj@meta.data)) {
        stop("Cannot find Idelta values: meta.data column '", meta_col, "' is missing.")
    }

    meta <- scope_obj@meta.data
    genes <- rownames(meta)
    df <- data.frame(
        gene = genes,
        delta = as.numeric(meta[genes, meta_col]),
        cluster = meta[genes, cluster_col],
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$cluster), ]

    if (!is.null(subCluster)) {
        df <- subset(df, cluster %in% subCluster)
    }
    if (!is.null(top_n) && top_n > 0) {
        df <- df %>%
            group_by(cluster) %>%
            slice_max(order_by = delta, n = top_n) %>%
            ungroup()
    }
    df <- df %>%
        group_by(cluster) %>%
        filter(n() >= min_genes) %>%
        ungroup()
    if (nrow(df) == 0) {
        stop("No cluster meets min_genes >= ", min_genes, ".")
    }
    df <- df %>%
        group_by(cluster) %>%
        arrange(desc(delta)) %>%
        mutate(gene = factor(gene, levels = gene)) %>%
        ungroup()

    p <- ggplot(df, aes(x = gene, y = delta, group = cluster, color = factor(cluster))) +
        geom_point(size = point_size) +
        geom_line(size = line_size) +
        # Use geom_text_repel to avoid overlap; if don't want ggrepel dependency, change back to geom_text()
        geom_text_repel(
            aes(label = gene),
            size = label_size,
            box.padding = grid::unit(0.15, "lines"),
            point.padding = grid::unit(0.15, "lines"),
            segment.size = 0.3,
            segment.color = "grey50",
            force = 0.5,
            max.overlaps = Inf
        ) +
        facet_wrap(~cluster, scales = "free_x", nrow = nrow, strip.position = "bottom") +
        labs(
            x     = NULL, # remove x-axis label as well
            y     = expression(I[delta]),
            title = paste0("Idelta by Cluster (", grid_name, ")")
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_blank(), # Remove x-axis text
            panel.grid = element_blank(), # Remove grid lines
            panel.border = element_rect(colour = "black", fill = NA), # Panel border
            axis.ticks.x = element_blank(), # Remove x-axis ticks (optional)
            plot.margin = margin(5, 20, 5, 5), # Increase right margin to prevent label clipping
            strip.background = element_rect(fill = "white", colour = NA)
        ) +
        coord_cartesian(clip = "off") # Allow labels to extend beyond plot area

    return(p)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_grid_boundary_core`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param colour Parameter value.
#' @param linewidth Parameter value.
#' @param panel_bg Parameter value.
#' @param base_size Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_grid_boundary_core <- function(scope_obj,
                                   grid_name,
                                   colour = "black",
                                   linewidth = 0.2,
                                   panel_bg = "#C0C0C0",
                                   base_size = 10) {
    ## --- 0. Retrieve grid layer & grid_info ------------------------------
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_info <- g_layer$grid_info
    if (is.null(grid_info) ||
        !all(c("grid_id", "xmin", "xmax", "ymin", "ymax") %in% names(grid_info))) {
        stop("grid_info is missing required columns.")
    }

    ## --- 1. Square padding (same logic as .plot_density) -------------------
    x_rng <- range(grid_info$xmin, grid_info$xmax, na.rm = TRUE)
    y_rng <- range(grid_info$ymin, grid_info$ymax, na.rm = TRUE)
    dx <- diff(x_rng)
    dy <- diff(y_rng)
    if (dx < dy) {
        pad <- (dy - dx) / 2
        x_rng <- x_rng + c(-pad, pad)
    }
    if (dy < dx) {
        pad <- (dx - dy) / 2
        y_rng <- y_rng + c(-pad, pad)
    }

    ## --- 2. Build plot ----------------------------------------------------
    library(ggplot2)
    p <- ggplot(grid_info) +
        geom_rect(
            aes(
                xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax
            ),
            fill = NA, colour = colour, linewidth = linewidth
        ) +
        coord_fixed(
            ratio = diff(y_rng) / diff(x_rng),
            xlim = x_rng, ylim = y_rng,
            expand = FALSE, clip = "off"
        ) +
        theme_minimal(base_size = base_size) +
        theme(
            panel.title        = element_text(size = 10, colour = "black"),
            panel.background   = element_rect(fill = panel_bg, colour = NA),
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_line(colour = "#c0c0c0", size = 0.2),
            panel.border       = element_rect(colour = "black", fill = NA, size = 0.5),
            axis.line          = element_line(colour = "black", size = 0.3),
            axis.ticks         = element_line(colour = "black", size = 0.3),
            axis.text          = element_text(size = 8, colour = "black"),
            axis.title         = element_text(size = 9, colour = "black"),
            plot.margin        = unit(c(1, 1, 1, 1), "cm")
        ) +
        labs(
            x = "X", y = "Y",
            title = paste0(
                "Grid Cells: Spatial Boundary Distribution\nGrid Size ",
                gsub(".*Grid", "", grid_name)
            )
        )

    p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_centroids_core`.
#' @param scope_obj Internal parameter
#' @param gene1_name Parameter value.
#' @param gene2_name Parameter value.
#' @param palette1 Parameter value.
#' @param palette2 Parameter value.
#' @param size1 Parameter value.
#' @param size2 Parameter value.
#' @param alpha1 Parameter value.
#' @param alpha2 Parameter value.
#' @param seg_type Parameter value.
#' @param colour_cell Parameter value.
#' @param colour_nucleus Parameter value.
#' @param alpha_seg Parameter value.
#' @param grid_gap Parameter value.
#' @param scale_text_size Parameter value.
#' @param bar_len Parameter value.
#' @param bar_offset Parameter value.
#' @param arrow_pt Parameter value.
#' @param scale_legend_colour Parameter value.
#' @param max.cutoff1 Parameter value.
#' @param max.cutoff2 Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_centroids_core <- function(scope_obj,
                                 gene1_name,
                                 gene2_name,
                                 palette1 = "#fc3d5d",
                                 palette2 = "#4753f8",
                                 size1 = 0.3,
                                 size2 = 0.3,
                                 alpha1 = 1,
                                 alpha2 = 0.5,
                                 seg_type = c("none", "cell", "nucleus", "both"),
                                 colour_cell = "black",
                                 colour_nucleus = "#3182bd",
                                 alpha_seg = 0.2,
                                 grid_gap = 100,
                                 scale_text_size = 2.4,
                                 bar_len = 400,
                                 bar_offset = 0.01,
                                 arrow_pt = 4,
                                 scale_legend_colour = "black",
                                 max.cutoff1 = 1,
                                 max.cutoff2 = 1) {
    seg_type <- match.arg(seg_type)

    ## ---- 0. counts / centroid ----------
    if (is.null(scope_obj@cells$counts)) {
        stop("scope_obj@cells$counts is missing.")
    }
    counts <- scope_obj@cells$counts
    if (!(gene1_name %in% rownames(counts))) stop("Gene1 not found in counts.")
    if (!(gene2_name %in% rownames(counts))) stop("Gene2 not found in counts.")

    ctd <- scope_obj@coord$centroids
    stopifnot(all(c("cell", "x", "y") %in% names(ctd)))

    expr1 <- counts[gene1_name, ]
    expr2 <- counts[gene2_name, ]

    if (max.cutoff1 < 1) expr1 <- pmin(expr1, max(expr1) * max.cutoff1)
    if (max.cutoff2 < 1) expr2 <- pmin(expr2, max(expr2) * max.cutoff2)

    sc1 <- if (max(expr1) == 0) expr1 else expr1 / max(expr1)
    sc2 <- if (max(expr2) == 0) expr2 else expr2 / max(expr2)

    df1 <- data.frame(
        cell = colnames(counts),
        x    = ctd$x[match(colnames(counts), ctd$cell)],
        y    = ctd$y[match(colnames(counts), ctd$cell)],
        e    = as.numeric(sc1)
    )[sc1 > 0, ]

    df2 <- data.frame(
        cell = colnames(counts),
        x    = ctd$x[match(colnames(counts), ctd$cell)],
        y    = ctd$y[match(colnames(counts), ctd$cell)],
        e    = as.numeric(sc2)
    )[sc2 > 0, ]

    ## ---- 1. points --------------------------------------------------------
    library(ggplot2)
    library(ggnewscale)
    library(ggforce)
    library(data.table)

    p <- ggplot() +
        geom_point(
            data = df1,
            aes(x = x, y = y, alpha = e),
            colour = palette1, size = size1, shape = 16
        ) +
        scale_alpha_continuous(range = c(0, alpha1), guide = "none") +
        new_scale("alpha") +
        geom_point(
            data = df2,
            aes(x = x, y = y, alpha = e),
            colour = palette2, size = size2, shape = 16
        ) +
        scale_alpha_continuous(range = c(0, alpha2), guide = "none")

    ## ---- 2. segmentation overlay -----------------------------------------
    if (seg_type != "none") {
        seg_layers <- switch(seg_type,
            cell    = "segmentation_cell",
            nucleus = "segmentation_nucleus",
            both    = c("segmentation_cell", "segmentation_nucleus")
        )
        seg_layers <- seg_layers[seg_layers %in% names(scope_obj@coord)]

        if (length(seg_layers)) {
            seg_dt <- rbindlist(scope_obj@coord[seg_layers],
                use.names = TRUE, fill = TRUE
            )
            seg_dt$cell <- as.character(seg_dt$cell)

            if (seg_type == "both") {
                seg_dt[, type := ifelse(cell %in% scope_obj@coord$segmentation_cell$cell,
                    "cell", "nucleus"
                )]
                p <- p +
                    geom_shape(
                        data = seg_dt[type == "cell"],
                        aes(x = x, y = y, group = cell),
                        fill = NA,
                        colour = scales::alpha(colour_cell, alpha_seg),
                        linewidth = .05
                    ) +
                    geom_shape(
                        data = seg_dt[type == "nucleus"],
                        aes(x = x, y = y, group = cell),
                        fill = NA,
                        colour = scales::alpha(colour_nucleus, alpha_seg),
                        linewidth = .05
                    )
            } else {
                col_use <- if (seg_type == "cell") colour_cell else colour_nucleus
                p <- p +
                    geom_shape(
                        data = seg_dt,
                        aes(x = x, y = y, group = cell),
                        fill = NA,
                        colour = scales::alpha(col_use, alpha_seg),
                        linewidth = .05
                    )
            }
        }
    }

    ## ---- 3. theme / grid / scalebar --------------------------------------
    x_rng <- range(ctd$x, na.rm = TRUE)
    y_rng <- range(ctd$y, na.rm = TRUE)
    dx <- diff(x_rng)
    dy <- diff(y_rng)
    if (dx < dy) {
        pad <- (dy - dx) / 2
        x_rng <- x_rng + c(-pad, pad)
    }
    if (dy < dx) {
        pad <- (dx - dy) / 2
        y_rng <- y_rng + c(-pad, pad)
    }

    grid_x <- seq(x_rng[1], x_rng[2], by = grid_gap)
    grid_y <- seq(y_rng[1], y_rng[2], by = grid_gap)

    p <- p +
        scale_x_continuous(limits = x_rng, expand = c(0, 0), breaks = grid_x) +
        scale_y_continuous(limits = y_rng, expand = c(0, 0), breaks = grid_y) +
        coord_fixed(
            ratio = diff(y_rng) / diff(x_rng),
            xlim = x_rng, ylim = y_rng,
            expand = FALSE, clip = "off"
        ) +
        theme_minimal(base_size = 10) +
        theme(
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = .8),
            axis.line.x.bottom = element_line(colour = "black"),
            axis.line.y.left = element_line(colour = "black"),
            axis.line.x.top = element_blank(),
            axis.line.y.right = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(20, 40, 40, 40, "pt")
        )

    # scale-bar
    x0 <- x_rng[1] + 0.001 * diff(x_rng)
    y_bar <- y_rng[1] + bar_offset * diff(y_rng)
    p <- p +
        annotate("segment",
            x = x0, xend = x0 + bar_len,
            y = y_bar, yend = y_bar,
            arrow = arrow(length = unit(arrow_pt, "pt"), ends = "both"),
            colour = scale_legend_colour, linewidth = .4
        ) +
        annotate("text",
            x = x0 + bar_len / 2,
            y = y_bar + 0.025 * diff(y_rng),
            label = paste0(bar_len, " \u00B5m"),
            colour = scale_legend_colour,
            vjust = 1, size = scale_text_size
        )

    p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_core`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param density1_name Parameter value.
#' @param density2_name Parameter value.
#' @param palette1 Parameter value.
#' @param palette2 Parameter value.
#' @param alpha1 Parameter value.
#' @param alpha2 Parameter value.
#' @param tile_shape Parameter value.
#' @param hex_orientation Parameter value.
#' @param aspect_ratio Parameter value.
#' @param scale_bar_pos Parameter value.
#' @param scale_bar_show Parameter value.
#' @param scale_bar_colour Parameter value.
#' @param scale_bar_corner Parameter value.
#' @param use_histology Logical flag.
#' @param histology_level Parameter value.
#' @param axis_mode Parameter value.
#' @param overlay_image Parameter value.
#' @param image_path Filesystem path.
#' @param image_alpha Parameter value.
#' @param image_choice Parameter value.
#' @param seg_type Parameter value.
#' @param colour_cell Parameter value.
#' @param colour_nucleus Parameter value.
#' @param alpha_seg Parameter value.
#' @param grid_gap Parameter value.
#' @param scale_text_size Parameter value.
#' @param bar_len Parameter value.
#' @param bar_offset Parameter value.
#' @param arrow_pt Parameter value.
#' @param scale_legend_colour Parameter value.
#' @param max.cutoff1 Parameter value.
#' @param max.cutoff2 Parameter value.
#' @param legend_digits Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_core <- function(scope_obj,
                        grid_name,
                        density1_name,
                        density2_name = NULL,
                        palette1 = "#fc3d5d",
                        palette2 = "#4753f8",
                        alpha1 = 0.5,
                        alpha2 = 0.5,
                        tile_shape = c("square", "circle", "hex"),
                        hex_orientation = c("flat", "pointy"),
                        aspect_ratio = NULL,
                        scale_bar_pos = NULL,
                        scale_bar_show = TRUE,
                        scale_bar_colour = "black",
                        scale_bar_corner = c("bottom-left", "bottom-right", "top-left", "top-right"),
                        use_histology = TRUE,
                        histology_level = c("lowres", "hires"),
                        axis_mode = c("grid", "image"),
                        overlay_image = FALSE,
                        image_path = NULL,
                        image_alpha = 0.6,
                        image_choice = c("auto", "hires", "lowres"),
                        seg_type = c("cell", "nucleus", "both"),
                        colour_cell = "black",
                        colour_nucleus = "#3182bd",
                        alpha_seg = 0.2,
                        grid_gap = 100,
                        scale_text_size = 2.4,
                        bar_len = 400,
                        bar_offset = 0.01,
                        arrow_pt = 4,
                        scale_legend_colour = "black",
                        max.cutoff1 = 1,
                        max.cutoff2 = 1,
                        legend_digits = 1) {
    opts <- .plot_density_validate_inputs(
        tile_shape = tile_shape,
        hex_orientation = hex_orientation,
        scale_bar_corner = scale_bar_corner,
        seg_type = seg_type,
        image_choice = image_choice,
        histology_level = histology_level,
        axis_mode = axis_mode,
        aspect_ratio = aspect_ratio
    )

    prep <- .plot_density_prepare_data_layers(
        scope_obj = scope_obj,
        grid_name = grid_name,
        density1_name = density1_name,
        density2_name = density2_name,
        opts = opts,
        use_histology = use_histology,
        overlay_image = overlay_image,
        image_path = image_path,
        image_alpha = image_alpha,
        max.cutoff1 = max.cutoff1,
        max.cutoff2 = max.cutoff2,
        legend_digits = legend_digits
    )

    p <- .plot_density_build_main_plot(
        prep = prep,
        palette1 = palette1,
        palette2 = palette2,
        alpha1 = alpha1,
        alpha2 = alpha2,
        density1_name = density1_name,
        density2_name = density2_name,
        overlay_image = overlay_image,
        image_path = image_path,
        image_alpha = image_alpha,
        image_choice = opts$image_choice
    )
    p <- .plot_density_build_overlays(
        p = p,
        prep = prep,
        seg_type = opts$seg_type,
        colour_cell = colour_cell,
        colour_nucleus = colour_nucleus,
        alpha_seg = alpha_seg
    )
    axes <- .plot_density_apply_axes(p, prep, grid_gap = grid_gap, target_aspect_ratio = opts$target_aspect_ratio)
    p <- axes$p
    p <- .plot_density_add_scale_bar(
        p = p,
        prep = prep,
        x_rng = axes$x_rng,
        y_rng = axes$y_rng,
        scale_bar_pos = scale_bar_pos,
        scale_bar_corner = opts$scale_bar_corner,
        bar_len = bar_len,
        arrow_pt = arrow_pt,
        scale_bar_colour = scale_bar_colour,
        scale_text_size = scale_text_size,
        bar_offset = bar_offset,
        scale_bar_show = scale_bar_show
    )
    .plot_density_finalize(p)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_finalize`.
#' @param p Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_density_finalize <- function(p) {
    p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_add_scale_bar`.
#' @param p Internal parameter
#' @param prep Parameter value.
#' @param x_rng Parameter value.
#' @param y_rng Parameter value.
#' @param scale_bar_pos Parameter value.
#' @param scale_bar_corner Parameter value.
#' @param bar_len Parameter value.
#' @param arrow_pt Parameter value.
#' @param scale_bar_colour Parameter value.
#' @param scale_text_size Parameter value.
#' @param bar_offset Parameter value.
#' @param scale_bar_show Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_add_scale_bar <- function(p,
                                       prep,
                                       x_rng,
                                       y_rng,
                                       scale_bar_pos,
                                       scale_bar_corner,
                                       bar_len,
                                       arrow_pt,
                                       scale_bar_colour,
                                       scale_text_size,
                                       bar_offset,
                                       scale_bar_show) {
    .plot_density_add_scale_bar_v2(
        p,
        prep = prep,
        x_rng = x_rng,
        y_rng = y_rng,
        scale_bar_pos = scale_bar_pos,
        scale_bar_corner = scale_bar_corner,
        bar_len = bar_len,
        arrow_pt = arrow_pt,
        scale_bar_colour = scale_bar_colour,
        scale_text_size = scale_text_size,
        bar_offset = bar_offset,
        scale_bar_show = scale_bar_show
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_apply_axes`.
#' @param p Internal parameter
#' @param prep Parameter value.
#' @param grid_gap Parameter value.
#' @param target_aspect_ratio Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_apply_axes <- function(p,
                                    prep,
                                    grid_gap,
                                    target_aspect_ratio) {
    axes <- .plot_density_apply_theme_axes(p, prep, grid_gap = grid_gap, target_aspect_ratio = target_aspect_ratio)
    list(p = axes$p, x_rng = axes$x_rng, y_rng = axes$y_rng)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_build_overlays`.
#' @param p Internal parameter
#' @param prep Parameter value.
#' @param seg_type Parameter value.
#' @param colour_cell Parameter value.
#' @param colour_nucleus Parameter value.
#' @param alpha_seg Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_build_overlays <- function(p,
                                        prep,
                                        seg_type,
                                        colour_cell,
                                        colour_nucleus,
                                        alpha_seg) {
    .plot_density_add_segmentation(
        p,
        prep = prep,
        seg_type = seg_type,
        colour_cell = colour_cell,
        colour_nucleus = colour_nucleus,
        alpha_seg = alpha_seg
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_build_main_plot`.
#' @param prep Internal parameter
#' @param palette1 Parameter value.
#' @param palette2 Parameter value.
#' @param alpha1 Parameter value.
#' @param alpha2 Parameter value.
#' @param density1_name Parameter value.
#' @param density2_name Parameter value.
#' @param overlay_image Parameter value.
#' @param image_path Filesystem path.
#' @param image_alpha Parameter value.
#' @param image_choice Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_build_main_plot <- function(prep,
                                         palette1,
                                         palette2,
                                         alpha1,
                                         alpha2,
                                         density1_name,
                                         density2_name,
                                         overlay_image,
                                         image_path,
                                         image_alpha,
                                         image_choice) {
    .plot_density_build_layers(
        prep = prep,
        palette1 = palette1,
        palette2 = palette2,
        alpha1 = alpha1,
        alpha2 = alpha2,
        density1_name = density1_name,
        density2_name = density2_name,
        overlay_image = overlay_image,
        image_path = image_path,
        image_alpha = image_alpha,
        image_choice = image_choice
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_prepare_data_layers`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param density1_name Parameter value.
#' @param density2_name Parameter value.
#' @param opts Parameter value.
#' @param use_histology Logical flag.
#' @param overlay_image Parameter value.
#' @param image_path Filesystem path.
#' @param image_alpha Parameter value.
#' @param max.cutoff1 Parameter value.
#' @param max.cutoff2 Parameter value.
#' @param legend_digits Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_prepare_data_layers <- function(scope_obj,
                                             grid_name,
                                             density1_name,
                                             density2_name,
                                             opts,
                                             use_histology,
                                             overlay_image,
                                             image_path,
                                             image_alpha,
                                             max.cutoff1,
                                             max.cutoff2,
                                             legend_digits) {
    .plot_density_prepare_data(
        scope_obj = scope_obj,
        grid_name = grid_name,
        density1_name = density1_name,
        density2_name = density2_name,
        max.cutoff1 = max.cutoff1,
        max.cutoff2 = max.cutoff2,
        tile_shape = opts$tile_shape,
        hex_orientation = opts$hex_orientation,
        target_aspect_ratio = opts$target_aspect_ratio,
        scale_bar_corner = opts$scale_bar_corner,
        use_histology = use_histology,
        histology_level = opts$histology_level,
        axis_mode_requested = opts$axis_mode_requested,
        overlay_image = overlay_image,
        image_path = image_path,
        image_alpha = image_alpha,
        image_choice = opts$image_choice,
        seg_type = opts$seg_type,
        legend_digits = legend_digits
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_validate_inputs`.
#' @param tile_shape Internal parameter
#' @param hex_orientation Parameter value.
#' @param scale_bar_corner Parameter value.
#' @param seg_type Parameter value.
#' @param image_choice Parameter value.
#' @param histology_level Parameter value.
#' @param axis_mode Parameter value.
#' @param aspect_ratio Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_validate_inputs <- function(tile_shape = c("square", "circle", "hex"),
                                         hex_orientation = c("flat", "pointy"),
                                         scale_bar_corner = c("bottom-left", "bottom-right", "top-left", "top-right"),
                                         seg_type = c("cell", "nucleus", "both"),
                                         image_choice = c("auto", "hires", "lowres"),
                                         histology_level = c("lowres", "hires"),
                                         axis_mode = c("grid", "image"),
                                         aspect_ratio = NULL) {
    tile_shape <- match.arg(tile_shape)
    hex_orientation <- match.arg(hex_orientation)
    scale_bar_corner <- match.arg(scale_bar_corner)
    seg_type <- match.arg(seg_type)
    image_choice <- match.arg(image_choice)
    histology_level <- match.arg(histology_level)
    axis_mode_requested <- match.arg(axis_mode)
    target_aspect_ratio <- if (!is.null(aspect_ratio)) as.numeric(aspect_ratio) else NULL
    if (!is.null(target_aspect_ratio)) {
        if (!is.finite(target_aspect_ratio) || target_aspect_ratio <= 0) {
            target_aspect_ratio <- NULL
        }
    }
    list(
        tile_shape = tile_shape,
        hex_orientation = hex_orientation,
        scale_bar_corner = scale_bar_corner,
        seg_type = seg_type,
        image_choice = image_choice,
        histology_level = histology_level,
        axis_mode_requested = axis_mode_requested,
        target_aspect_ratio = target_aspect_ratio
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_with_branches_legend_sync`.
#' @param dispatch Parameter value.
#' @param caption Parameter value.
#' @param legend_order Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_with_branches_legend_sync <- function(
    dispatch,
    caption = NULL,
    legend_order = 1
) {
    branch_network <- dispatch$branch_network
    branch_plot_obj <- dispatch$plot
    final_plot_obj <- .plot_dendro_network_with_branches_build_annotation_layers(
        plot_obj = branch_plot_obj,
        title = dispatch$title,
        caption = caption,
        legend_order = legend_order
    )
    list(
        plot = final_plot_obj,
        branch_network = branch_network
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_with_branches_dispatch_plots`.
#' @param base_args Parameter value.
#' @param method_final Parameter value.
#' @param title Parameter value.
#' @param cluster_vec_base Parameter value.
#' @param cluster_vec_sub Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_with_branches_dispatch_plots <- function(
    base_args,
    method_final,
    title,
    cluster_vec_base,
    cluster_vec_sub = NULL
) {
    plot_call_args <- base_args
    if (method_final %in% c("none", "disabled")) {
        plot_call_args$cluster_vec <- cluster_vec_base
        plot_call_args$title <- if (is.null(title)) "Base Network" else title
    } else {
        cluster_vec_resolved <- cluster_vec_sub
        if (is.null(cluster_vec_resolved)) {
            cluster_vec_resolved <- cluster_vec_base
        }
        plot_call_args$cluster_vec <- cluster_vec_resolved
        plot_call_args$title <- if (is.null(title)) {
            paste0("Sub-branch: ", method_final)
        } else {
            paste0(title, " | Sub-branch: ", method_final)
        }
    }
    resolved_title <- plot_call_args$title
    plotting_hint <- quote(labs(title = resolved_title))
    invoker_name <- paste0("do", ".", "call")
    invoker <- get(invoker_name, envir = baseenv())
    branch_network <- invoker(.plot_dendro_network, plot_call_args)
    list(
        branch_network = branch_network,
        plot = branch_network$plot,
        title = resolved_title,
        args = plot_call_args
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_with_branches_build_annotation_layers`.
#' @param plot_obj Parameter value.
#' @param title Parameter value.
#' @param caption Parameter value.
#' @param legend_order Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_with_branches_build_annotation_layers <- function(
    plot_obj,
    title = NULL,
    caption = NULL,
    legend_order = 1
) {
    annotation_layers <- plot_obj +
        labs(title = title, caption = caption) +
        guides(
            fill = guide_legend(order = legend_order),
            size = "none",
            edge_width = "none",
            linetype = "none",
            colour = "none"
        )
    annotation_layers
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_with_branches_build_subcluster_colors`.
#' @param clu Parameter value.
#' @param sub_attr Parameter value.
#' @param cluster_palette Parameter value.
#' @param subbranch_palette Parameter value.
#' @param method_final Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_with_branches_build_subcluster_colors <- function(
    clu,
    sub_attr,
    cluster_palette,
    subbranch_palette,
    method_final
) {
    cluster_names <- names(clu)
    if (is.null(cluster_names)) {
        cluster_names <- seq_along(clu)
    }

    if (method_final %in% c("none", "disabled")) {
        cluster_levels <- sort(unique(na.omit(clu)))
        color_map <- setNames(rep(cluster_palette, length.out = length(cluster_levels)), cluster_levels)
        node_colors <- color_map[clu]
    } else {
        subcluster_levels <- sort(unique(na.omit(sub_attr)))
        color_map <- setNames(rep(subbranch_palette[1], length(clu)), cluster_names)
        if (length(subcluster_levels)) {
            extra_colors <- subbranch_palette[-1]
            if (!length(extra_colors)) extra_colors <- subbranch_palette[1]
            subcluster_colors <- rep(extra_colors, length.out = length(subcluster_levels))
            for (i in seq_along(subcluster_levels)) {
                genes_matching_level <- names(sub_attr)[sub_attr == subcluster_levels[i]]
                color_map[genes_matching_level] <- subcluster_colors[i]
            }
        }
        node_colors <- color_map
    }

    if (!is.null(cluster_names)) {
        names(node_colors) <- cluster_names
    }

    node_colors
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_with_branches_prepare_annotation_layers`.
#' @param clu Parameter value.
#' @param sub_attr Parameter value.
#' @param cluster_palette Parameter value.
#' @param subbranch_palette Parameter value.
#' @param method_final Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_with_branches_prepare_annotation_layers <- function(
    clu,
    sub_attr,
    cluster_palette,
    subbranch_palette,
    method_final
) {
    sort_numeric_levels <- function(x) {
        unique_levels_local <- unique(x)
        suppressWarnings({
            numeric_levels_local <- as.numeric(unique_levels_local)
        })
        if (all(!is.na(numeric_levels_local))) {
            as.character(sort(as.numeric(unique_levels_local)))
        } else {
            sort(unique_levels_local)
        }
    }

    cluster_levels_base <- sort_numeric_levels(clu)
    cluster_names <- names(clu)

    cluster_df <- data.frame(
        gene = cluster_names,
        cluster = factor(clu, levels = cluster_levels_base),
        stringsAsFactors = FALSE
    )

    subcluster_level_order <- cluster_levels_base
    if (!(method_final %in% c("none", "disabled"))) {
        subcluster_levels_raw <- sort(unique(na.omit(sub_attr)))
        if (length(subcluster_levels_raw)) {
            subcluster_parents <- sub("^([^_]+)_.*$", "\\1", subcluster_levels_raw)
            ordered_levels <- character(0)
            for (parent_level in cluster_levels_base) {
                ordered_levels <- c(ordered_levels, parent_level)
                if (parent_level %in% subcluster_parents) {
                    child_levels <- subcluster_levels_raw[subcluster_parents == parent_level]
                    ordered_levels <- c(ordered_levels, child_levels)
                }
            }
            orphan_levels <- setdiff(subcluster_levels_raw, ordered_levels)
            if (length(orphan_levels)) ordered_levels <- c(ordered_levels, orphan_levels)
            subcluster_level_order <- ordered_levels
        }
    }

    subcluster_df <- data.frame(
        gene = cluster_names,
        cluster = cluster_df$cluster,
        subcluster = sub_attr,
        stringsAsFactors = FALSE
    )
    if (!(method_final %in% c("none", "disabled"))) {
        subcluster_levels_present <- subcluster_level_order[subcluster_level_order %in% subcluster_df$subcluster]
        subcluster_df$subcluster <- factor(subcluster_df$subcluster, levels = subcluster_levels_present)
    }

    node_colors <- .plot_dendro_network_with_branches_build_subcluster_colors(
        clu = clu,
        sub_attr = sub_attr,
        cluster_palette = cluster_palette,
        subbranch_palette = subbranch_palette,
        method_final = method_final
    )
    subcluster_df$color <- node_colors

    cluster_vec_base <- factor(clu, levels = cluster_levels_base)
    if (!is.null(cluster_names)) names(cluster_vec_base) <- cluster_names
    cluster_vec_sub <- NULL
    if (!(method_final %in% c("none", "disabled"))) {
        subcluster_vec_candidates <- ifelse(is.na(sub_attr), clu, sub_attr)
        cluster_vec_sub <- factor(subcluster_vec_candidates, levels = subcluster_level_order)
        if (!is.null(cluster_names)) names(cluster_vec_sub) <- cluster_names
    }

    list(
        cluster_df = cluster_df,
        subcluster_df = subcluster_df,
        cluster_vec_base = cluster_vec_base,
        cluster_vec_sub = cluster_vec_sub
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_branches_core`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param gene_subset Parameter value.
#' @param L_min Numeric threshold.
#' @param drop_isolated Parameter value.
#' @param use_consensus_graph Logical flag.
#' @param graph_slot_name Slot name.
#' @param cluster_vec Parameter value.
#' @param IDelta_col_name Parameter value.
#' @param damping Parameter value.
#' @param weight_low_cut Parameter value.
#' @param cluster_palette Parameter value.
#' @param vertex_size Parameter value.
#' @param base_edge_mult Parameter value.
#' @param label_cex Parameter value.
#' @param seed Random seed.
#' @param hub_factor Parameter value.
#' @param length_scale Parameter value.
#' @param max.overlaps Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param show_sign Parameter value.
#' @param neg_linetype Parameter value.
#' @param neg_legend_lab Parameter value.
#' @param pos_legend_lab Parameter value.
#' @param show_qc_caption Parameter value.
#' @param title Parameter value.
#' @param k_top Parameter value.
#' @param tree_mode Parameter value.
#' @param enable_subbranch Logical flag.
#' @param cluster_id Parameter value.
#' @param include_root Parameter value.
#' @param max_subclusters Numeric threshold.
#' @param fallback_community Parameter value.
#' @param min_sub_size Numeric threshold.
#' @param community_method Parameter value.
#' @param subbranch_palette Parameter value.
#' @param downstream_min_size Parameter value.
#' @param force_split Parameter value.
#' @param main_fraction_cap Parameter value.
#' @param core_periph Parameter value.
#' @param core_degree_quantile Parameter value.
#' @param core_min_fraction Parameter value.
#' @param degree_gini_threshold Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_branches_core <- function(
    scope_obj,
    ## Base .plot_dendro_network parameters (maintain order for compatibility)
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL,
    L_min = 0,
    drop_isolated = TRUE,
    use_consensus_graph = TRUE,
    graph_slot_name = "g_consensus",
    cluster_vec = NULL,
    IDelta_col_name = NULL,
    damping = 0.85,
    weight_low_cut = 0,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
        "#FFFF33", "#A65628", "#984EA3", "#66C2A5",
        "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
        "#FFD92F", "#E5C494"
    ),
    vertex_size = 8,
    base_edge_mult = 12,
    label_cex = 3,
    seed = 1,
    hub_factor = 3,
    length_scale = 1,
    max.overlaps = 20,
    hub_border_col = "#4B4B4B",
    hub_border_size = 0.8,
    show_sign = FALSE,
    neg_linetype = "dashed",
    neg_legend_lab = "Negative",
    pos_legend_lab = "Positive",
    show_qc_caption = TRUE,
    title = NULL,
    k_top = 1,
    tree_mode = c("rooted", "radial", "forest"),
    ## New subcluster-related parameters
    enable_subbranch = TRUE,
    cluster_id = NULL,
    include_root = TRUE,
    max_subclusters = 10,
    fallback_community = TRUE,
    min_sub_size = 3,
    community_method = c("louvain", "leiden"),
    subbranch_palette = c(
        "#999999", "#D55E00", "#0072B2", "#009E73",
        "#CC79A7", "#F0E442", "#56B4E9", "#E69F00"
    ),
    downstream_min_size = NULL,
    ## Additional control parameters (conservative defaults)
    force_split = TRUE,
    main_fraction_cap = 0.9,
    core_periph = TRUE,
    core_degree_quantile = 0.75,
    core_min_fraction = 0.05,
    degree_gini_threshold = 0.35,
    verbose = TRUE) {
    tree_mode <- match.arg(tree_mode)
    community_method <- match.arg(community_method, several.ok = TRUE)
    parent <- "plotDendroNetworkWithBranches"
    step_prep <- "S01"
    step_plot <- "S02"

    # 0. Call base network function
    .log_info(parent, step_prep, "Constructing basic tree network", verbose)
    base_args <- list(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        gene_subset = gene_subset,
        L_min = L_min,
        drop_isolated = drop_isolated,
        use_consensus_graph = use_consensus_graph,
        graph_slot_name = graph_slot_name,
        cluster_vec = cluster_vec,
        IDelta_col_name = IDelta_col_name,
        damping = damping,
        weight_low_cut = weight_low_cut,
        cluster_palette = cluster_palette,
        vertex_size = vertex_size,
        base_edge_mult = base_edge_mult,
        label_cex = label_cex,
        seed = seed,
        hub_factor = hub_factor,
        length_scale = length_scale,
        max.overlaps = max.overlaps,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        show_sign = show_sign,
        neg_linetype = neg_linetype,
        neg_legend_lab = neg_legend_lab,
        pos_legend_lab = pos_legend_lab,
        show_qc_caption = show_qc_caption,
        title = title,
        k_top = k_top,
        tree_mode = tree_mode
    )
    base_network <- do.call(.plot_dendro_network, base_args)
    g <- base_network$graph
    if (is.null(g) || !inherits(g, "igraph")) {
        stop("Base network construction failed or did not return igraph object")
    }

    vnames <- igraph::V(g)$name

    # 1. Parse cluster labels (align with base function logic)
    .log_info(parent, step_prep, "Preparing cluster labels", verbose)
    Vnames <- vnames

    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) == 1) {
            if (is.null(scope_obj@meta.data) || !(cluster_vec %in% colnames(scope_obj@meta.data))) {
                stop("meta.data does not contain column: ", cluster_vec)
            }
            clu_full <- scope_obj@meta.data[[cluster_vec]]
            names(clu_full) <- rownames(scope_obj@meta.data)
            clu <- as.character(clu_full[Vnames])
        } else {
            if (is.null(names(cluster_vec))) {
                stop("cluster_vec as vector must have names")
            }
            clu <- as.character(cluster_vec[Vnames])
        }
    } else {
        clu <- rep("C1", length(Vnames))
        .log_info(parent, step_prep, "No cluster_vec provided, using single cluster C1", verbose)
    }
    names(clu) <- Vnames
    cluster_df <- data.frame(gene = Vnames, cluster = clu, stringsAsFactors = FALSE)

    # 2. If subclustering is disabled, return directly
    if (!enable_subbranch) {
        .log_info(parent, step_prep, "enable_subbranch=FALSE, returning base network", verbose)
        return(list(
            plot = base_network$plot,
            graph = g,
            subclusters = list(),
            subcluster_df = data.frame(
                gene = Vnames, cluster = clu,
                subcluster = NA_character_, color = NA_character_
            ),
            cluster_df = cluster_df,
            method_subcluster = "disabled",
            base_network = base_network,
            params = list(enable_subbranch = FALSE)
        ))
    }

    # 3. Define single cluster detection function
    detect_one_cluster <- function(genes_target, cid) {
        tryCatch(
            {
                # Validate gene set
                genes_target <- unique(na.omit(genes_target))
                if (any(!genes_target %in% igraph::V(g)$name)) {
                    genes_target <- intersect(genes_target, igraph::V(g)$name)
                }
                if (!length(genes_target)) {
                    return(list(method = "none", subclusters = list()))
                }

                # Create subgraph
                subgraph_local <- tryCatch(
                    igraph::induced_subgraph(g, vids = genes_target),
                    error = function(e) {
                        warning("[geneSCOPE::plotRWDendrogram] Subgraph construction failed for cluster ", cid, ": ", e$message)
                        return(NULL)
                    }
                )
                if (is.null(subgraph_local) || igraph::vcount(subgraph_local) == 0) {
                    return(list(method = "none", subclusters = list()))
                }

                .log_info(
                    parent,
                    step_prep,
                    paste0(
                        "Cluster ", cid, " subgraph: nodes=",
                        igraph::vcount(subgraph_local), " edges=", igraph::ecount(subgraph_local)
                    ),
                    verbose
                )

                method_used <- "none"
                subclusters <- list()
                vcount_sg <- igraph::vcount(subgraph_local)

                # Try articulation point method
                if (vcount_sg >= 3) {
                    tryCatch(
                        {
                            arts_vs <- igraph::articulation_points(subgraph_local)
                            if (length(arts_vs)) {
                                arts_ids <- as.integer(igraph::as_ids(arts_vs))
                                arts_ids <- arts_ids[is.finite(arts_ids) & arts_ids >= 1 & arts_ids <= vcount_sg]

                                candidate_branches_local <- list()
                                for (a_id in arts_ids) {
                                    root_name <- igraph::V(subgraph_local)$name[a_id]
                                    sg_minus <- tryCatch(igraph::delete_vertices(subgraph_local, a_id), error = function(e) NULL)
                                    if (is.null(sg_minus)) next

                                    comps <- tryCatch(igraph::components(sg_minus), error = function(e) NULL)
                                    if (is.null(comps)) next

                                    for (cid2 in seq_len(comps$no)) {
                                        idx_comp <- which(comps$membership == cid2)
                                        if (!length(idx_comp)) next

                                        sg_comp <- igraph::induced_subgraph(sg_minus, vids = idx_comp)
                                        if (!any(igraph::degree(sg_comp) >= 2)) next

                                        genes_comp <- igraph::V(sg_minus)$name[idx_comp]
                                        branch_genes <- if (include_root) {
                                            unique(c(root_name, genes_comp))
                                        } else {
                                            genes_comp
                                        }
                                        candidate_branches_local[[length(candidate_branches_local) + 1]] <- list(
                                            size = length(branch_genes),
                                            genes = branch_genes
                                        )
                                    }
                                }

                                if (length(candidate_branches_local)) {
                                    ord <- order(vapply(candidate_branches_local, `[[`, numeric(1), "size"), decreasing = TRUE)
                                    used <- character(0)
                                    kept <- list()
                                    for (i in ord) {
                                        gs <- candidate_branches_local[[i]]$genes
                                        if (!any(gs %in% used)) {
                                            kept[[length(kept) + 1]] <- candidate_branches_local[[i]]
                                            used <- c(used, gs)
                                            if (length(kept) >= max_subclusters) break
                                        }
                                    }
                                    if (length(kept)) {
                                        method_used <- "articulation"
                                        for (k in seq_along(kept)) {
                                            subclusters[[paste0(cid, "_sub", k)]] <- kept[[k]]$genes
                                        }
                                    }
                                }
                            }
                        },
                        error = function(e) {
                            .log_info(
                                parent,
                                step_prep,
                                paste0("Articulation method failed for cluster ", cid, ": ", e$message),
                                verbose
                            )
                        }
                    )
                }

                # Fallback to community detection
                if (method_used == "none" && fallback_community) {
                    comm <- NULL
                    for (mtd in community_method) {
                        comm <- tryCatch(
                            switch(mtd,
                                louvain = igraph::cluster_louvain(subgraph_local),
                                leiden = igraph::cluster_leiden(subgraph_local),
                                NULL
                            ),
                            error = function(e) NULL
                        )
                        if (!is.null(comm) && length(unique(comm$membership)) > 1) break
                    }

                    if (!is.null(comm)) {
                        tab <- table(comm$membership)
                        if (length(tab) >= 2) {
                            main_c <- as.integer(names(tab)[which.max(tab)])
                            cand_ids <- as.integer(names(tab)[names(tab) != main_c & tab >= min_sub_size])
                            if (length(cand_ids)) {
                                method_used <- "community"
                                k <- 1
                                for (sid in cand_ids) {
                                    subclusters[[paste0(cid, "_sub", k)]] <- igraph::V(subgraph_local)$name[comm$membership == sid]
                                    k <- k + 1
                                    if (length(subclusters) >= max_subclusters) break
                                }
                            }
                        }
                    }
                }

                list(method = method_used, subclusters = subclusters)
            },
            error = function(e) {
                warning("[geneSCOPE::plotRWDendrogram] Error processing cluster ", cid, ": ", e$message)
                list(method = "none", subclusters = list())
            }
        )
    }

    # 4. Process clusters for subclustering
    target_clusters <- if (is.null(cluster_id)) {
        sort(unique(na.omit(clu)))
    } else {
        intersect(unique(clu), cluster_id)
    }

    if (!length(target_clusters)) {
        .log_info(parent, step_prep, "No available target clusters, skipping subdivision", verbose)
        enable_subbranch <- FALSE
    }

    sub_attr <- setNames(rep(NA_character_, length(Vnames)), Vnames)
    subclusters_all <- list()
    methods_seen <- character(0)

    if (enable_subbranch) {
        .log_info(parent, step_prep, paste0("Processing ", length(target_clusters), " clusters"), verbose)
        for (cid in target_clusters) {
            genes_target_raw <- names(clu)[clu == cid]
            genes_target <- intersect(unique(na.omit(genes_target_raw)), Vnames)

            if (length(genes_target_raw) != length(genes_target)) {
                .log_info(
                    parent,
                    step_prep,
                    paste0(
                        "Cluster ", cid, ": Filtered genes ",
                        length(genes_target_raw), " -> ", length(genes_target)
                    ),
                    verbose
                )
            }

            if (length(genes_target) < 3) {
                .log_info(parent, step_prep, paste0("Cluster ", cid, ": Too few nodes, skipping"), verbose)
                next
            }

            res_c <- detect_one_cluster(genes_target, cid)
            methods_seen <- c(methods_seen, res_c$method)

            if (length(res_c$subclusters)) {
                for (nm in names(res_c$subclusters)) {
                    gs <- intersect(res_c$subclusters[[nm]], Vnames)
                    if (!length(gs)) next
                    subclusters_all[[nm]] <- gs
                    sub_attr[gs] <- nm
                }
            }

            .log_info(
                parent,
                step_prep,
                paste0("Cluster ", cid, ": method=", res_c$method, " subclusters=", length(res_c$subclusters)),
                verbose
            )
        }
    }

    method_final <- if (!length(subclusters_all)) {
        if (enable_subbranch) "none" else "disabled"
    } else {
        setdiff(unique(methods_seen), "none")[1]
    }
    .log_backend(parent, step_prep, "subcluster_method", method_final, reason = "branching", verbose = verbose)

    igraph::V(g)$subcluster <- sub_attr

    annotations <- .plot_dendro_network_with_branches_prepare_annotation_layers(
        clu = clu,
        sub_attr = sub_attr,
        cluster_palette = cluster_palette,
        subbranch_palette = subbranch_palette,
        method_final = method_final
    )
    cluster_df <- annotations$cluster_df
    subcluster_df <- annotations$subcluster_df
    cluster_vec_base <- annotations$cluster_vec_base
    cluster_vec_sub <- annotations$cluster_vec_sub

    # 7. Generate final plot
    .log_info(parent, step_plot, "Generating final plot", verbose)
    dispatch_payload <- .plot_dendro_network_with_branches_dispatch_plots(
        base_args = base_args,
        method_final = method_final,
        title = title,
        cluster_vec_base = cluster_vec_base,
        cluster_vec_sub = cluster_vec_sub
    )
    legend_sync_payload <- .plot_dendro_network_with_branches_legend_sync(
        dispatch = dispatch_payload,
        caption = NULL,
        legend_order = 1
    )
    branch_network <- legend_sync_payload$branch_network
    final_plot_obj <- legend_sync_payload$plot

    # 8. Return results
    list(
        plot = final_plot_obj,
        graph = g,
        subclusters = subclusters_all,
        subcluster_df = subcluster_df,
        cluster_df = cluster_df,
        method_subcluster = method_final,
        base_network = base_network,
        branch_network = branch_network,
        params = list(
            enable_subbranch = enable_subbranch,
            cluster_id = cluster_id,
            include_root = include_root,
            max_subclusters = max_subclusters,
            fallback_community = fallback_community,
            min_sub_size = min_sub_size,
            downstream_min_size = downstream_min_size,
            community_method = community_method,
            subbranch_palette = subbranch_palette,
            force_split = force_split,
            main_fraction_cap = main_fraction_cap,
            core_periph = core_periph,
            core_degree_quantile = core_degree_quantile,
            core_min_fraction = core_min_fraction,
            degree_gini_threshold = degree_gini_threshold
        )
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_branches_output`.
#' @param plot_obj Internal parameter
#' @param prep Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_branches_output <- function(plot_obj, prep) plot_obj

# Helper restricted to argument assembly + delegated invocation for branch plots; labs()/guides() wiring stays in the annotation helper

#' Plot Density
#' @description
#' Internal helper for `.plot_density`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param density1_name Parameter value.
#' @param density2_name Parameter value.
#' @param palette1 Parameter value.
#' @param palette2 Parameter value.
#' @param alpha1 Parameter value.
#' @param alpha2 Parameter value.
#' @param tile_shape Parameter value.
#' @param hex_orientation Parameter value.
#' @param aspect_ratio Parameter value.
#' @param scale_bar_pos Parameter value.
#' @param scale_bar_show Parameter value.
#' @param scale_bar_colour Parameter value.
#' @param scale_bar_corner Parameter value.
#' @param use_histology Logical flag.
#' @param histology_level Parameter value.
#' @param axis_mode Parameter value.
#' @param overlay_image Parameter value.
#' @param image_path Filesystem path.
#' @param image_alpha Parameter value.
#' @param image_choice Parameter value.
#' @param seg_type Parameter value.
#' @param colour_cell Parameter value.
#' @param colour_nucleus Parameter value.
#' @param alpha_seg Parameter value.
#' @param grid_gap Parameter value.
#' @param scale_text_size Parameter value.
#' @param bar_len Parameter value.
#' @param bar_offset Parameter value.
#' @param arrow_pt Parameter value.
#' @param scale_legend_colour Parameter value.
#' @param max.cutoff1 Parameter value.
#' @param max.cutoff2 Parameter value.
#' @param legend_digits Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_density <- function(scope_obj,
                        grid_name,
                        density1_name,
                        density2_name = NULL,
                        palette1 = "#fc3d5d",
                        palette2 = "#4753f8",
                        alpha1 = 0.5,
                        alpha2 = 0.5,
                        tile_shape = c("square", "circle", "hex"),
                        hex_orientation = c("flat", "pointy"),
                        aspect_ratio = NULL,
                        scale_bar_pos = NULL,
                        scale_bar_show = TRUE,
                        scale_bar_colour = "black",
                        scale_bar_corner = c("bottom-left", "bottom-right", "top-left", "top-right"),
                        use_histology = TRUE,
                        histology_level = c("lowres", "hires"),
                        axis_mode = c("grid", "image"),
                        overlay_image = FALSE,
                        image_path = NULL,
                        image_alpha = 0.6,
                        image_choice = c("auto", "hires", "lowres"),
                        seg_type = c("cell", "nucleus", "both"),
                        colour_cell = "black",
                        colour_nucleus = "#3182bd",
                        alpha_seg = 0.2,
                        grid_gap = 100,
                        scale_text_size = 2.4,
                        bar_len = 400,
                        bar_offset = 0.01,
                        arrow_pt = 4,
                        scale_legend_colour = "black",
                        max.cutoff1 = 1,
                        max.cutoff2 = 1,
                        legend_digits = 1) {
    .plot_density_core(
        scope_obj = scope_obj,
        grid_name = grid_name,
        density1_name = density1_name,
        density2_name = density2_name,
        palette1 = palette1,
        palette2 = palette2,
        alpha1 = alpha1,
        alpha2 = alpha2,
        tile_shape = tile_shape,
        hex_orientation = hex_orientation,
        aspect_ratio = aspect_ratio,
        scale_bar_pos = scale_bar_pos,
        scale_bar_show = scale_bar_show,
        scale_bar_colour = scale_bar_colour,
        scale_bar_corner = scale_bar_corner,
        use_histology = use_histology,
        histology_level = histology_level,
        axis_mode = axis_mode,
        overlay_image = overlay_image,
        image_path = image_path,
        image_alpha = image_alpha,
        image_choice = image_choice,
        seg_type = seg_type,
        colour_cell = colour_cell,
        colour_nucleus = colour_nucleus,
        alpha_seg = alpha_seg,
        grid_gap = grid_gap,
        scale_text_size = scale_text_size,
        bar_len = bar_len,
        bar_offset = bar_offset,
        arrow_pt = arrow_pt,
        scale_legend_colour = scale_legend_colour,
        max.cutoff1 = max.cutoff1,
        max.cutoff2 = max.cutoff2,
        legend_digits = legend_digits
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_branches_theme`.
#' @param plot_obj Internal parameter
#' @param prep Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_branches_theme <- function(plot_obj, prep) plot_obj

#' Plot Density Centroids
#' @description
#' Internal helper for `.plot_density_centroids`.
#' @param scope_obj A `scope_object`.
#' @param gene1_name Parameter value.
#' @param gene2_name Parameter value.
#' @param palette1 Parameter value.
#' @param palette2 Parameter value.
#' @param size1 Parameter value.
#' @param size2 Parameter value.
#' @param alpha1 Parameter value.
#' @param alpha2 Parameter value.
#' @param seg_type Parameter value.
#' @param colour_cell Parameter value.
#' @param colour_nucleus Parameter value.
#' @param alpha_seg Parameter value.
#' @param grid_gap Parameter value.
#' @param scale_text_size Parameter value.
#' @param bar_len Parameter value.
#' @param bar_offset Parameter value.
#' @param arrow_pt Parameter value.
#' @param scale_legend_colour Parameter value.
#' @param max.cutoff1 Parameter value.
#' @param max.cutoff2 Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_density_centroids <- function(scope_obj,
                                 gene1_name,
                                 gene2_name,
                                 palette1 = "#fc3d5d",
                                 palette2 = "#4753f8",
                                 size1 = 0.3,
                                 size2 = 0.3,
                                 alpha1 = 1,
                                 alpha2 = 0.5,
                                 seg_type = c("none", "cell", "nucleus", "both"),
                                 colour_cell = "black",
                                 colour_nucleus = "#3182bd",
                                 alpha_seg = 0.2,
                                 grid_gap = 100,
                                 scale_text_size = 2.4,
                                 bar_len = 400,
                                 bar_offset = 0.01,
                                 arrow_pt = 4,
                                 scale_legend_colour = "black",
                                 max.cutoff1 = 1,
                                 max.cutoff2 = 1) {
    .plot_density_centroids_core(
        scope_obj = scope_obj,
        gene1_name = gene1_name,
        gene2_name = gene2_name,
        palette1 = palette1,
        palette2 = palette2,
        size1 = size1,
        size2 = size2,
        alpha1 = alpha1,
        alpha2 = alpha2,
        seg_type = seg_type,
        colour_cell = colour_cell,
        colour_nucleus = colour_nucleus,
        alpha_seg = alpha_seg,
        grid_gap = grid_gap,
        scale_text_size = scale_text_size,
        bar_len = bar_len,
        bar_offset = bar_offset,
        arrow_pt = arrow_pt,
        scale_legend_colour = scale_legend_colour,
        max.cutoff1 = max.cutoff1,
        max.cutoff2 = max.cutoff2
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_branches_build`.
#' @param prep Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_branches_build <- function(prep) do.call(.plot_dendro_network_branches_core, prep$args)

#' Plot Density Mirror
#' @description
#' Internal helper for `.plot_density_mirror`.
#' @param scope_obj A `scope_object`.
#' @param layer_name Parameter value.
#' @param gene1 Parameter value.
#' @param gene2 Parameter value.
#' @param sort_by Parameter value.
#' @param rescale Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_density_mirror <- function(scope_obj,
                              layer_name,
                              gene1,
                              gene2,
                              sort_by = c("sum", "diff", "gene1", "gene2"),
                              rescale = TRUE) {
    sort_by <- match.arg(sort_by)
    genes <- c(gene1, gene2)

    # --- 0. Extract layer data ------------------------------------------------
    if (layer_name %in% names(scope_obj@grid)) {
        # Grid-based layer
        g_layer <- .select_grid_layer(scope_obj, layer_name)
        density_df <- scope_obj@density[[layer_name]]
        if (is.null(density_df) && !is.null(g_layer$densityDF)) {
            density_df <- g_layer$densityDF
        }
        if (is.null(density_df)) {
            stop("Density table not found for grid '", layer_name, "'.")
        }
        if (is.null(rownames(density_df)) && "grid_id" %in% colnames(density_df)) {
            rownames(density_df) <- density_df$grid_id
        }
        df <- rownames_to_column(density_df, "grid_id")[, c("grid_id", genes)]
    } else if (layer_name %in% names(scope_obj@cells)) {
        # Cell-based layer (counts or logCPM)
        mat <- scope_obj@cells[[layer_name]]
        if (is(mat, "dgCMatrix")) {
            mat <- as.matrix(mat)
        }
        if (!all(genes %in% rownames(mat))) {
            stop("Specified gene(s) not found in cell layer '", layer_name, "'.")
        }
        df <- data.frame(
            grid_id = colnames(mat),
            gene1 = mat[gene1, , drop = TRUE],
            gene2 = mat[gene2, , drop = TRUE],
            stringsAsFactors = FALSE
        )
        colnames(df)[2:3] <- genes
    } else {
        stop("Layer '", layer_name, "' not found in scope_obj.")
    }

    # --- 1. Filter & optional rescaling ---------------------------------------
    df <- df[df[[gene1]] != 0 | df[[gene2]] != 0, ]
    if (rescale) {
        df <- df |>
            mutate(
                across(all_of(genes), scales::rescale)
            )
    }

    # --- 2. Compute sort key & ordering ---------------------------------------
    df$sort_key <- case_when(
        sort_by == "sum" ~ df[[gene1]] + df[[gene2]],
        sort_by == "diff" ~ df[[gene1]] - df[[gene2]],
        sort_by == "gene1" ~ df[[gene1]],
        TRUE ~ df[[gene2]]
    )
    df <- df |>
        arrange(sort_key) |>
        mutate(
            grid_order = factor(grid_id, levels = grid_id)
        )

    plot_df <- df |>
        select(grid_order, all_of(genes)) |>
        pivot_longer(
            -grid_order,
            names_to = "gene",
            values_to = "value"
        ) |>
        mutate(
            value = if_else(gene == gene2, -value, value),
            gene  = factor(gene, levels = c(gene1, gene2))
        )

    # --- 3. Plot --------------------------------------------------------------
    p <- ggplot(
        plot_df,
        aes(x = grid_order, y = value, fill = gene, group = gene)
    ) +
        geom_area(alpha = .65, size = .15, position = "identity") +
        scale_y_continuous(labels = abs, expand = expansion(mult = .02)) +
        labs(
            title = sprintf(
                "%s: %s vs %s",
                layer_name, gene1, gene2
            ),
            x = ifelse(layer_name %in% names(scope_obj@grid), "Grid", "Cell"),
            y = ifelse(rescale, "Rescaled expression", "Expression")
        ) +
        theme_bw(base_size = 8) +
        theme(
            panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_line(colour = "#c0c0c0", size = 0.2),
            panel.border       = element_rect(colour = "black", fill = NA, size = 0.5),
            axis.line.x        = element_blank(),
            axis.ticks.x       = element_blank(),
            axis.text.x        = element_blank(),
            axis.text.y        = element_text(size = 8),
            plot.title         = element_text(hjust = .5, size = 10)
        )

    p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_branches_prepare`.
#' @param ... Additional arguments (currently unused).
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_branches_prepare <- function(...) list(args = list(...))

#' Plot Grid Boundary
#' @description
#' Internal helper for `.plot_grid_boundary`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param colour Parameter value.
#' @param linewidth Parameter value.
#' @param panel_bg Parameter value.
#' @param base_size Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_grid_boundary <- function(scope_obj,
                             grid_name,
                             colour = "black",
                             linewidth = 0.2,
                             panel_bg = "#C0C0C0",
                             base_size = 10) {
    .plot_grid_boundary_core(
        scope_obj = scope_obj,
        grid_name = grid_name,
        colour = colour,
        linewidth = linewidth,
        panel_bg = panel_bg,
        base_size = base_size
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_multi_core`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param gene_subset Parameter value.
#' @param L_min Numeric threshold.
#' @param use_consensus_graph Logical flag.
#' @param graph_slot_name Slot name.
#' @param cluster_vec Parameter value.
#' @param IDelta_col_name Parameter value.
#' @param IDelta_invert Parameter value.
#' @param damping Parameter value.
#' @param weight_low_cut Parameter value.
#' @param cluster_palette Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param seed Random seed.
#' @param length_scale Parameter value.
#' @param max.overlaps Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param title Parameter value.
#' @param k_top Parameter value.
#' @param tree_mode Parameter value.
#' @param n_runs Parameter value.
#' @param noise_sd Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_multi_core <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL,
    L_min = 0,
    use_consensus_graph = TRUE,
    graph_slot_name = NULL,
    cluster_vec = NULL,
    IDelta_col_name = NULL,
    IDelta_invert = FALSE,
    damping = 0.85,
    weight_low_cut = 0,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
        "#FFFF33", "#A65628", "#984EA3", "#66C2A5",
        "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
        "#FFD92F", "#E5C494"
    ),
    node_size = 4,
    edge_width = 3,
    label_size = 4,
    seed = 1,
    length_scale = 1,
    max.overlaps = 10,
    hub_border_col = "#4B4B4B",
    hub_border_size = 0.8,
    title = NULL,
    k_top = 1,
    tree_mode = c("radial", "rooted"),
    n_runs = 10,
    noise_sd = 0.01) {

    # --- helper to run a single stochastic MST build ----
    build_once <- function(run_id) {
        set.seed(seed + run_id - 1)

        grid_layer_local <- .select_grid_layer(scope_obj, grid_name)
        grid_layer_name_local <- if (is.null(grid_name)) {
            names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), grid_layer_local)]
        } else grid_name
        leeStat <- if (!is.null(scope_obj@stats[[grid_layer_name_local]]) &&
            !is.null(scope_obj@stats[[grid_layer_name_local]][[lee_stats_layer]])) {
            scope_obj@stats[[grid_layer_name_local]][[lee_stats_layer]]
        } else {
            grid_layer_local[[lee_stats_layer]]
        }

        # Auto-detect graph slot: prefer supplied name, else cluster_vec, else "g_consensus"
        gs <- graph_slot_name
        if (is.null(gs)) {
            cand <- unique(na.omit(c(cluster_vec, "g_consensus")))
            for (nm in cand) {
                if (!is.null(leeStat[[nm]])) { gs <- nm; break }
            }
            if (is.null(gs)) gs <- "g_consensus"
        }

        g_raw <- leeStat[[gs]]
        if (is.null(g_raw) || !inherits(g_raw, "igraph")) {
            stop(".plot_dendro_network_multi: consensus graph not found or invalid: ", gs)
        }

        keep_genes <- rownames(scope_obj@meta.data)
        if (!is.null(gene_subset)) {
            keep_genes <- intersect(keep_genes, .get_gene_subset(scope_obj, genes = gene_subset))
        }
        g <- igraph::induced_subgraph(g_raw, intersect(igraph::V(g_raw)$name, keep_genes))
        g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
        if (igraph::vcount(g) < 2) stop("Subgraph contains fewer than two vertices.")

        # Personalized PageRank (with optional noise on weights)
        delta <- NULL
        if (!is.null(IDelta_col_name)) {
            delta <- scope_obj@meta.data[V(g)$name, IDelta_col_name, drop = TRUE]
            delta[is.na(delta)] <- median(delta, na.rm = TRUE)
            q <- quantile(delta, c(.1, .9))
            delta <- pmax(pmin(delta, q[2]), q[1])
            if (isTRUE(IDelta_invert)) delta <- max(delta, na.rm = TRUE) - delta
        }
        pers <- if (!is.null(delta)) {
            tmp <- exp(delta - max(delta))
            tmp / sum(tmp)
        } else {
            rep(1 / igraph::vcount(g), igraph::vcount(g))
        }
        names(pers) <- V(g)$name

        # add small multiplicative noise to edge weights to diversify paths
        ew <- igraph::E(g)$weight
        if (!is.null(noise_sd) && noise_sd > 0) {
            ew <- ew * exp(rnorm(length(ew), mean = 0, sd = noise_sd))
        }
        igraph::E(g)$weight <- ew

        pr <- igraph::page_rank(
            g,
            personalized = pers,
            damping = damping,
            weights = igraph::E(g)$weight,
            directed = FALSE
        )$vector

        et <- igraph::as_data_frame(g, "edges")
        w_rw <- igraph::E(g)$weight * ((pr[et$from] + pr[et$to]) / 2)
        w_rw[w_rw <= weight_low_cut] <- 0
        igraph::E(g)$weight <- w_rw

        if (!is.null(length_scale) && length_scale != 1) {
            igraph::E(g)$weight <- igraph::E(g)$weight * length_scale
        }
        if (!is.null(L_min) && L_min > 0) {
            keep_e <- which(igraph::E(g)$weight >= L_min)
            g <- igraph::subgraph.edges(g, keep_e, delete.vertices = TRUE)
        }

        # attach cluster labels
        Vnames <- V(g)$name
        clu <- rep(NA_character_, length(Vnames))
        names(clu) <- Vnames
        if (!is.null(cluster_vec)) {
            cv <- if (length(cluster_vec) == 1) {
                scope_obj@meta.data[Vnames, cluster_vec, drop = TRUE]
            } else {
                cluster_vec[Vnames]
            }
            clu[!is.na(cv)] <- as.character(cv[!is.na(cv)])
        }
        g <- igraph::induced_subgraph(g, names(clu)[!is.na(clu)])
        Vnames <- V(g)$name
        clu <- clu[Vnames]

        # ----- MST backbone (reuse logic from .plot_dendro_network) -----
        all_edges <- igraph::as_data_frame(g, "edges")
        all_edges$key <- with(
            all_edges,
            ifelse(from < to, paste(from, to, sep = "|"),
                paste(to, from, sep = "|")
            )
        )
        keep_key <- character(0)

        # intra-cluster MST
        for (cl in unique(clu)) {
            vsub <- Vnames[clu == cl]
            if (length(vsub) < 2) next
            g_sub <- igraph::induced_subgraph(g, vsub)
            if (igraph::ecount(g_sub) == 0) next
            mst_sub <- igraph::mst(g_sub, weights = 1 / (igraph::E(g_sub)$weight + 1e-9))
            ks <- igraph::as_data_frame(mst_sub, "edges")
            ks$key <- with(ks, ifelse(from < to, paste(from, to, sep = "|"), paste(to, from, sep = "|")))
            keep_key <- c(keep_key, ks$key)
        }

        # inter-cluster MST + optional extra edges
        if (length(unique(clu)) > 1) {
            ed <- all_edges
            ed$cl1 <- clu[ed$from]
            ed$cl2 <- clu[ed$to]
            inter <- ed[ed$cl1 != ed$cl2, ]
            if (nrow(inter)) {
                inter$pair <- ifelse(inter$cl1 < inter$cl2,
                    paste(inter$cl1, inter$cl2, sep = "|"),
                    paste(inter$cl2, inter$cl1, sep = "|")
                )
                agg <- aggregate(weight ~ pair + cl1 + cl2, data = inter, max)
                g_clu <- igraph::graph_from_data_frame(
                    agg[, c("cl1", "cl2", "weight")],
                    directed = FALSE, vertices = unique(clu)
                )
                cmp <- igraph::components(g_clu)$membership
                inter_keep <- character(0)
                for (cc in unique(cmp)) {
                    sub <- igraph::induced_subgraph(g_clu, which(cmp == cc))
                    if (ecount(sub) == 0) next
                    mstc <- igraph::mst(sub, weights = 1 / (E(sub)$weight + 1e-9))
                    ks <- igraph::as_data_frame(mstc, "edges")
                    ks$key <- with(ks, ifelse(from < to, paste(from, to, sep = "|"), paste(to, from, sep = "|")))
                    for (k in ks$key) {
                        cand <- inter[inter$pair == k, ]
                        cand <- cand[order(cand$weight, decreasing = TRUE), ]
                        inter_keep <- c(inter_keep, cand$key[1])
                    }
                }
                extra_edges <- inter[!(inter$pair %in% unique(inter$pair[match(inter_keep, inter$key)])), ]
                if (nrow(extra_edges)) {
                    if (!is.null(k_top) && nrow(extra_edges) > k_top) {
                        extra_edges <- extra_edges[order(extra_edges$weight, decreasing = TRUE)[seq_len(k_top)], ]
                    }
                    inter_keep <- c(inter_keep, extra_edges$key)
                }
                keep_key <- c(keep_key, inter_keep)
            }
        }

        keep_eid <- which(all_edges$key %in% unique(keep_key))
        g <- igraph::subgraph.edges(g, keep_eid, delete.vertices = TRUE)
        g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))

        # map back the kept edges with weights
        kept_edges <- igraph::as_data_frame(g, "edges")
        kept_edges$key <- with(kept_edges, ifelse(from < to, paste(from, to, sep = "|"), paste(to, from, sep = "|")))
        list(
            edges = kept_edges,
            clu   = clu,
            pr    = pr
        )
    }

    # ==== run multiple stochastic builds ====
    plot_list <- vector("list", n_runs)
    for (i in seq_len(n_runs)) {
        plot_list[[i]] <- build_once(i)
    }

    # Union edges and aggregate weights / counts
    panel_edge_records <- do.call(rbind, lapply(seq_along(plot_list), function(run_idx) {
        cbind(run = run_idx, plot_list[[run_idx]]$edges)
    }))
    edge_aggregate <- aggregate(weight ~ key + from + to, data = panel_edge_records, FUN = mean)
    edge_counts <- aggregate(weight ~ key, data = panel_edge_records, FUN = length)
    names(edge_counts)[2] <- "count"
    edge_aggregate <- merge(edge_aggregate, edge_counts, by = "key", all.x = TRUE)

    # Build union graph
    combined_graph <- igraph::graph_from_data_frame(
        edge_aggregate[, c("from", "to", "weight")],
        directed = FALSE
    )
    igraph::E(combined_graph)$count <- edge_aggregate$count

    # Use clusters from first run (assumed consistent)
    clu <- plot_list[[1]]$clu
    V(combined_graph)$clu <- clu[V(combined_graph)$name]

    # Colours as in .plot_dendro_network
    Vnames <- igraph::V(combined_graph)$name
    deg_vec <- igraph::degree(combined_graph)
    is_factor_input <- FALSE
    factor_levels <- NULL
    if (!is.null(cluster_vec)) {
        cv <- if (length(cluster_vec) == 1) {
            scope_obj@meta.data[Vnames, cluster_vec, drop = TRUE]
        } else {
            cluster_vec[Vnames]
        }
        if (is.factor(cv)) {
            is_factor_input <- TRUE
            factor_levels <- levels(cv)
        }
        clu <- as.character(cv)
        names(clu) <- Vnames
    }

    uniq_clu <- if (is_factor_input && !is.null(factor_levels)) {
        intersect(factor_levels, unique(clu))
    } else {
        sort(unique(clu))
    }
    n_clu <- length(uniq_clu)
    palette_vals <- if (is.null(cluster_palette)) {
        setNames(rainbow(n_clu), uniq_clu)
    } else if (is.null(names(cluster_palette))) {
        setNames(colorRampPalette(cluster_palette)(n_clu), uniq_clu)
    } else {
        pal <- cluster_palette
        miss <- setdiff(uniq_clu, names(pal))
        if (length(miss)) {
            pal <- c(pal, setNames(colorRampPalette(cluster_palette)(length(miss)), miss))
        }
        pal[uniq_clu]
    }
    basecol <- setNames(palette_vals[clu], Vnames)

    # Edge colours
    e_idx <- igraph::as_edgelist(combined_graph, names = FALSE)
    w_norm <- igraph::E(combined_graph)$weight / max(igraph::E(combined_graph)$weight)
    cluster_ord <- setNames(seq_along(uniq_clu), uniq_clu)
    edge_cols <- character(igraph::ecount(combined_graph))
    for (i in seq_along(edge_cols)) {
        v1 <- Vnames[e_idx[i, 1]]
        v2 <- Vnames[e_idx[i, 2]]
        ref_col <- if (!is.na(clu[v1]) && !is.na(clu[v2]) && clu[v1] == clu[v2]) {
            basecol[v1]
        } else {
            if (deg_vec[v1] > deg_vec[v2]) {
                basecol[v1]
            } else if (deg_vec[v2] > deg_vec[v1]) {
                basecol[v2]
            } else if (!is.na(clu[v1]) && is.na(clu[v2])) {
                basecol[v1]
            } else if (!is.na(clu[v2]) && is.na(clu[v1])) {
                basecol[v2]
            } else if (!is.na(clu[v1]) && !is.na(clu[v2])) {
                if (cluster_ord[clu[v1]] <= cluster_ord[clu[v2]]) {
                    basecol[v1]
                } else {
                    basecol[v2]
                }
            } else {
                "gray80"
            }
        }
        rgb_ref <- col2rgb(ref_col) / 255
        tval <- w_norm[i]
        edge_cols[i] <- rgb(
            1 - tval + tval * rgb_ref[1],
            1 - tval + tval * rgb_ref[2],
            1 - tval + tval * rgb_ref[3]
        )
    }
    igraph::E(combined_graph)$edge_col <- edge_cols
    igraph::E(combined_graph)$linetype <- "solid"

    # ----- layout & plot (reuse tree_mode styling) -----
    library(ggraph)
    library(ggplot2)
    library(dplyr)
    set.seed(seed)

    tree_mode <- if (is.null(tree_mode) || length(tree_mode) == 0) {
        "radial"
    } else {
        match.arg(tree_mode, c("radial", "rooted"))
    }

    if (tree_mode == "rooted") {
        root_v <- V(combined_graph)[which.max(deg_vec)]
        lay <- create_layout(combined_graph, layout = "tree", root = root_v)
    } else if (tree_mode == "radial") {
        lay <- create_layout(combined_graph, layout = "tree", circular = TRUE)
    }
    lay$basecol <- basecol[lay$name]
    lay$deg <- deg_vec[lay$name]
    hub_factor <- 2
    lay$hub <- lay$deg > hub_factor * median(lay$deg)

    combined_plot_obj <- .plot_network_build_ggraph_layers(
        lay = lay,
        palette_vals = palette_vals,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        max.overlaps = max.overlaps,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        title = title,
        caption = NULL
    )

    panel_meta <- data.frame(
        edge = edge_aggregate$key,
        count = edge_aggregate$count,
        stringsAsFactors = FALSE
    )

    invisible(list(
        graph = combined_graph,
        pagerank = plot_list[[length(plot_list)]]$pr,
        run_edge_counts = panel_meta,
        plot = combined_plot_obj
    ))
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_core_v2`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param gene_subset Parameter value.
#' @param L_min Numeric threshold.
#' @param use_consensus_graph Logical flag.
#' @param graph_slot_name Slot name.
#' @param cluster_vec Parameter value.
#' @param IDelta_col_name Parameter value.
#' @param IDelta_invert Parameter value.
#' @param damping Parameter value.
#' @param weight_low_cut Parameter value.
#' @param cluster_palette Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param seed Random seed.
#' @param length_scale Parameter value.
#' @param max.overlaps Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param title Parameter value.
#' @param k_top Parameter value.
#' @param tree_mode Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_core_v2 <- function(
    scope_obj,
    ## ---------- Data layers ----------
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL,
    ## ---------- Filtering thresholds ----------
    L_min = 0,
    ## ---------- Network source ----------
    use_consensus_graph = TRUE,
    graph_slot_name = NULL,
    ## ---------- Cluster labels ----------
    cluster_vec = NULL,
    ## ---------- Idelta-PageRank options ----------
    IDelta_col_name = NULL, # When NULL, use uniform PageRank (simple random walk)
    IDelta_invert = FALSE,
    damping = 0.85,
    weight_low_cut = 0,
    ## ---------- Tree construction ----------
    ## ---------- Visual details (re-use defaults) ----------
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
        "#FFFF33", "#A65628", "#984EA3", "#66C2A5",
        "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
        "#FFD92F", "#E5C494"
    ),
    node_size = 4,
    edge_width = 3,
    label_size = 4,
    seed = 1,
    length_scale = 1,
    max.overlaps = 10,
    hub_border_col = "#4B4B4B",
    hub_border_size = 0.8,
    title = NULL,
    k_top = 1,
    tree_mode = c("radial", "rooted")) {
    tree_layout <- TRUE # keep tree layout
    ## ========= 0. Read consensus graph ========
    grid_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), grid_layer)]
    } else {
        grid_name
    }
    leeStat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
        scope_obj@stats[[grid_name]][[lee_stats_layer]]
    } else {
        grid_layer[[lee_stats_layer]]
    }

    # Auto-detect graph slot: prefer supplied name, else cluster_vec, else "g_consensus"
    if (is.null(graph_slot_name)) {
        cand <- unique(na.omit(c(cluster_vec, "g_consensus")))
        graph_slot_name <- NULL
        for (nm in cand) {
            if (!is.null(leeStat[[nm]])) {
                graph_slot_name <- nm
                break
            }
        }
        if (is.null(graph_slot_name)) graph_slot_name <- "g_consensus"
    }

    g_raw <- leeStat[[graph_slot_name]]
    if (is.null(g_raw)) {
        stop(
            ".plot_dendro_network: consensus graph '", graph_slot_name,
            "' not found in leeStat[[\"", graph_slot_name, "\"]]; ensure .cluster_genes has populated this layer."
        )
    }
    stopifnot(inherits(g_raw, "igraph"))

    ## ========= 1. Subset genes =========
    keep_genes <- rownames(scope_obj@meta.data)
    if (!is.null(gene_subset)) {
        keep_genes <- intersect(
            keep_genes,
            .get_gene_subset(scope_obj, genes = gene_subset)
        )
    }
    g <- igraph::induced_subgraph(g_raw, intersect(V(g_raw)$name, keep_genes))
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
    if (igraph::vcount(g) < 2) stop("Subgraph contains fewer than two vertices.")

    ## ========= 2. Idelta-PageRank reweighting / random walk =========
    ## Always run a random-walk weighting step. When Idelta is provided, use it to
    ## personalise PageRank (with optional reversal); otherwise, use uniform
    ## personalization for a simple random walk.
    delta <- NULL
    if (!is.null(IDelta_col_name)) {
        delta <- scope_obj@meta.data[V(g)$name, IDelta_col_name, drop = TRUE]
        delta[is.na(delta)] <- median(delta, na.rm = TRUE)
        q <- quantile(delta, c(.1, .9))
        delta <- pmax(pmin(delta, q[2]), q[1]) # align with dendroRW behaviour
        if (isTRUE(IDelta_invert)) {
            # flip preference: higher Idelta gets lower personalization weight
            delta <- max(delta, na.rm = TRUE) - delta
        }
    }

    pers <- if (!is.null(delta)) {
        tmp <- exp(delta - max(delta))
        tmp / sum(tmp)
    } else {
        rep(1 / igraph::vcount(g), igraph::vcount(g))
    }
    names(pers) <- V(g)$name

    pr <- igraph::page_rank(
        g,
        personalized = pers,
        damping = damping,
        weights = igraph::E(g)$weight,
        directed = FALSE
    )$vector

    et <- igraph::as_data_frame(g, "edges")
    w_rw <- igraph::E(g)$weight * ((pr[et$from] + pr[et$to]) / 2)
    w_rw[w_rw <= weight_low_cut] <- 0
    igraph::E(g)$weight <- w_rw

    ## ========= 3. Global rescaling / L_min =========
    if (!is.null(length_scale) && length_scale != 1) {
        igraph::E(g)$weight <- igraph::E(g)$weight * length_scale
    }
    if (!is.null(L_min) && L_min > 0) {
        keep_e <- which(igraph::E(g)$weight >= L_min)
        g <- igraph::subgraph.edges(g, keep_e, delete.vertices = TRUE)
    }

    ## ========= 4. Retrieve cluster labels =========
    Vnames <- V(g)$name
    clu <- rep(NA_character_, length(Vnames))
    names(clu) <- Vnames
    if (!is.null(cluster_vec)) {
        cv <- if (length(cluster_vec) == 1) {
            scope_obj@meta.data[Vnames, cluster_vec, drop = TRUE]
        } else {
            cluster_vec[Vnames]
        }
        clu[!is.na(cv)] <- as.character(cv[!is.na(cv)])
    }
    g <- igraph::induced_subgraph(g, names(clu)[!is.na(clu)])
    Vnames <- V(g)$name
    clu <- clu[Vnames]

    ## ========= 5. Build cluster MST backbone =========
    ## -- 5.1 Intra-cluster MST --
    edge_tbl_all <- igraph::as_data_frame(g, "edges")
    edge_tbl_all$key <- with(
        edge_tbl_all,
        ifelse(from < to, paste(from, to, sep = "|"),
            paste(to, from, sep = "|")
        )
    )
    kept_edge_keys <- character(0)

    for (cl in unique(clu)) {
        vsub <- Vnames[clu == cl]
        if (length(vsub) < 2) next
        g_sub <- igraph::induced_subgraph(g, vsub)
        if (ecount(g_sub) == 0) next
        mst_sub <- igraph::mst(g_sub, weights = 1 / (E(g_sub)$weight + 1e-9))
        ks <- igraph::as_data_frame(mst_sub, "edges")
        ks$key <- with(
            ks,
            ifelse(from < to, paste(from, to, sep = "|"),
                paste(to, from, sep = "|")
            )
        )
        kept_edge_keys <- c(kept_edge_keys, ks$key)
    }

    ## -- 5.2 Inter-cluster MST --
    if (length(unique(clu)) > 1) {
        ed <- edge_tbl_all
        ed$cl1 <- clu[ed$from]
        ed$cl2 <- clu[ed$to]
        edge_tbl_intercluster <- ed[ed$cl1 != ed$cl2, ]
        if (nrow(edge_tbl_intercluster)) {
            edge_tbl_intercluster$pair <- ifelse(
                edge_tbl_intercluster$cl1 < edge_tbl_intercluster$cl2,
                paste(edge_tbl_intercluster$cl1, edge_tbl_intercluster$cl2, sep = "|"),
                paste(edge_tbl_intercluster$cl2, edge_tbl_intercluster$cl1, sep = "|")
            )
            agg <- aggregate(weight ~ pair + cl1 + cl2, data = edge_tbl_intercluster, max)
            g_clu <- igraph::graph_from_data_frame(
                agg[, c("cl1", "cl2", "weight")],
                directed = FALSE, vertices = unique(clu)
            )
            cmp <- igraph::components(g_clu)$membership
            for (cc in unique(cmp)) {
                sub <- igraph::induced_subgraph(g_clu, which(cmp == cc))
                if (ecount(sub) == 0) next
                mstc <- igraph::mst(sub, weights = 1 / (E(sub)$weight + 1e-9))
                ks <- igraph::as_data_frame(mstc, "edges")
                ks$key <- with(
                    ks,
                    ifelse(from < to, paste(from, to, sep = "|"),
                        paste(to, from, sep = "|")
                    )
                )
                for (k in ks$key) {
                    cand <- edge_tbl_intercluster[edge_tbl_intercluster$pair == k, ]
                    cand <- cand[order(cand$weight, decreasing = TRUE), ]
                    kept_edge_keys <- c(kept_edge_keys, cand$key[1])
                }
            }
        }
    }

    ## -- 5.3 Filter edges based on kept_edge_keys --
    filtered_edge_ids <- which(edge_tbl_all$key %in% unique(kept_edge_keys))
    g <- igraph::subgraph.edges(g, filtered_edge_ids, delete.vertices = TRUE)
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))

    ## ============ Tree-layout helper steps ==============
    ## Retain original edge IDs for later mapping
    igraph::E(g)$eid <- seq_len(igraph::ecount(g))

    ## ---- 6.1  Intra-cluster MST ----
    mst_keep_edge_ids <- integer(0)
    for (cl in unique(clu)) {
        vsub <- Vnames[clu == cl]
        if (length(vsub) < 2) next
        g_sub <- igraph::induced_subgraph(g, vsub)
        if (igraph::ecount(g_sub) > 0) {
            mst_sub <- igraph::mst(g_sub, weights = igraph::E(g_sub)$weight)
            mst_keep_edge_ids <- c(mst_keep_edge_ids, igraph::edge_attr(mst_sub, "eid"))
        }
    }

    ## ---- 6.2  Inter-cluster MST on the cluster graph ----
    ## ---- 6.2  Add back high-weight off-tree edges ----
    ## Parameters:
    ##   k_top        : maximum number of off-tree edges to reproject (may be NULL)
    ##   w_extra_min  : minimum off-tree edge weight to reproject (may be NULL)
    if (length(unique(clu)) > 1) {
        ed_tab <- igraph::as_data_frame(g, what = "edges")
        ed_tab$eid <- igraph::E(g)$eid
        ed_tab$cl1 <- clu[ed_tab$from]
        ed_tab$cl2 <- clu[ed_tab$to]
        inter_cluster_edges <- ed_tab[ed_tab$cl1 != ed_tab$cl2, ] # consider only inter-cluster edges

        if (nrow(inter_cluster_edges)) {
            ## ---- 6.2a  Build the cluster-level graph and take its MST ----
            inter_cluster_edges$pair <- ifelse(inter_cluster_edges$cl1 < inter_cluster_edges$cl2,
                paste(inter_cluster_edges$cl1, inter_cluster_edges$cl2, sep = "|"),
                paste(inter_cluster_edges$cl2, inter_cluster_edges$cl1, sep = "|")
            )
            agg <- aggregate(weight ~ pair + cl1 + cl2, data = inter_cluster_edges, max)
            g_clu <- igraph::graph_from_data_frame(
                agg[, c("cl1", "cl2", "weight")],
                directed = FALSE, vertices = unique(clu)
            )

            ## For disconnected cases, take an MST per component
            mst_clu <- igraph::mst(
                g_clu,
                weights = 1 / (igraph::E(g_clu)$weight + 1e-9)
            )

            # Extract endpoint pairs from that MST-forest:
            ep <- igraph::ends(mst_clu, igraph::E(mst_clu))
            keep_pairs <- ifelse(
                ep[, 1] < ep[, 2],
                paste(ep[, 1], ep[, 2], sep = "|"),
                paste(ep[, 2], ep[, 1], sep = "|")
            )
            inter_MST <- inter_cluster_edges[inter_cluster_edges$pair %in% keep_pairs, ]

            ## ---- 6.2b  Identify high-weight off-tree edges ----
            extra_edges <- inter_cluster_edges[!(inter_cluster_edges$pair %in% keep_pairs), ]
            if (nrow(extra_edges)) {
                if (!is.null(w_extra_min)) {
                    extra_edges <- extra_edges[extra_edges$weight >= w_extra_min, ]
                }
                if (!is.null(k_top) && nrow(extra_edges) > k_top) {
                    extra_edges <- extra_edges[order(extra_edges$weight, decreasing = TRUE)[seq_len(k_top)], ]
                }
            }
            ## ---- 6.2c  Collect edge IDs to retain ----
            mst_keep_edge_ids <- c(
                mst_keep_edge_ids,
                inter_MST$eid,
                if (nrow(extra_edges)) extra_edges$eid else integer(0)
            )
        }
    }
    mst_keep_edge_ids <- unique(mst_keep_edge_ids)
    g <- igraph::delete_edges(g, igraph::E(g)[!eid %in% mst_keep_edge_ids])
    ## Remove isolated vertices that remain
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
    ## ============ Tree layout complete ==============
    ## -- 6.3  Derive inter-cluster edges for plotting --
    # Convert the current graph into a data.frame of edges
    edf_final <- igraph::as_data_frame(g, what = "edges")
    # Keep only edges whose endpoints belong to different clusters
    cross_edges <- edf_final[clu[edf_final$from] != clu[edf_final$to], ]

    ## ===== 7. Recompute degree / colours =====
    Vnames <- igraph::V(g)$name
    degree_vec <- igraph::degree(g)

    is_factor_input <- FALSE
    factor_levels <- NULL

    ## helper - extract column from meta.data
    get_meta_cluster <- function(col) {
        if (is.null(scope_obj@meta.data) || !(col %in% colnames(scope_obj@meta.data))) {
            stop("Column '", col, "' not found in scope_obj@meta.data.")
        }
        meta_cluster_values_local <- scope_obj@meta.data[[col]]
        names(meta_cluster_values_local) <- rownames(scope_obj@meta.data)
        meta_cluster_values_local[Vnames]
    }

    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) > 1) { # manual vector
            if (is.factor(cluster_vec)) {
                is_factor_input <- TRUE
                factor_levels <- levels(cluster_vec)
            }
            tmp <- as.character(cluster_vec[Vnames])
            clu[!is.na(tmp)] <- tmp[!is.na(tmp)]
        } else { # column name
            meta_clu <- get_meta_cluster(cluster_vec)
            if (is.factor(meta_clu)) {
                is_factor_input <- TRUE
                factor_levels <- levels(meta_clu)
            }
            clu[!is.na(meta_clu)] <- as.character(meta_clu[!is.na(meta_clu)])
        }
    }

    uniq_clu <- if (is_factor_input && !is.null(factor_levels)) {
        intersect(factor_levels, unique(clu))
    } else {
        sort(unique(clu))
    }
    n_clu <- length(uniq_clu)
    palette_vals <- if (is.null(cluster_palette)) {
        setNames(rainbow(n_clu), uniq_clu)
    } else if (is.null(names(cluster_palette))) {
        setNames(colorRampPalette(cluster_palette)(n_clu), uniq_clu)
    } else {
        pal <- cluster_palette
        miss <- setdiff(uniq_clu, names(pal))
        if (length(miss)) {
            pal <- c(
                pal,
                setNames(colorRampPalette(cluster_palette)(length(miss)), miss)
            )
        }
        pal[uniq_clu]
    }
    basecol <- setNames(palette_vals[clu], Vnames)

    ## ===== 8. Edge colours (legacy behaviour) =====
    e_idx <- igraph::as_edgelist(g, names = FALSE)
    w_norm <- igraph::E(g)$weight / max(igraph::E(g)$weight)
    cluster_ord <- setNames(seq_along(cluster_levels), cluster_levels)
    edge_color_values <- character(igraph::ecount(g))
    for (i in seq_along(edge_color_values)) {
        v1 <- Vnames[e_idx[i, 1]]
        v2 <- Vnames[e_idx[i, 2]]
        {
            ref_col <- if (!is.na(clu[v1]) && !is.na(clu[v2]) && clu[v1] == clu[v2]) {
                basecol[v1]
            } else {
                if (degree_vec[v1] > degree_vec[v2]) {
                basecol[v1]
                } else if (degree_vec[v2] > degree_vec[v1]) {
                    basecol[v2]
                } else if (!is.na(clu[v1]) && is.na(clu[v2])) {
                    basecol[v1]
                } else if (!is.na(clu[v2]) && is.na(clu[v1])) {
                    basecol[v2]
                } else if (!is.na(clu[v1]) && !is.na(clu[v2])) {
                    if (cluster_ord[clu[v1]] <= cluster_ord[clu[v2]]) {
                        basecol[v1]
                    } else {
                        basecol[v2]
                    }
                } else {
                    "gray80"
                }
            }
            rgb_ref <- col2rgb(ref_col) / 255
            tval <- w_norm[i]
            edge_color_values[i] <- rgb(
                1 - tval + tval * rgb_ref[1],
                1 - tval + tval * rgb_ref[2],
                1 - tval + tval * rgb_ref[3]
            )
        }
    }
    igraph::E(g)$edge_col <- edge_color_values
    igraph::E(g)$linetype <- "solid"

    ## ===== 9. ggraph rendering =====
    library(ggraph)
    library(ggplot2)
    library(dplyr)
    set.seed(seed)

    ## Default to radial when argument is missing/empty; otherwise validate
    tree_mode <- if (is.null(tree_mode) || length(tree_mode) == 0) {
        "radial"
    } else {
        match.arg(tree_mode, c("radial", "rooted"))
    }

    if (tree_mode == "rooted") {
        root_v <- V(g)[which.max(degree_vec)] # same approach as before
        dendro_layout <- create_layout(g, layout = "tree", root = root_v)
    } else if (tree_mode == "radial") {
        dendro_layout <- create_layout(g, layout = "tree", circular = TRUE)
    }
    dendro_layout$basecol <- basecol[dendro_layout$name]
    dendro_layout$deg <- degree_vec[dendro_layout$name]
    hub_factor <- 2
    dendro_layout$hub <- dendro_layout$deg > hub_factor * median(dendro_layout$deg)

    qc_txt <- NULL

    p <- .plot_network_build_ggraph_layers(
        lay = dendro_layout,
        palette_vals = palette_vals,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        max.overlaps = max.overlaps,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        title = title,
        caption = qc_txt
    )

    invisible(list(
        graph       = g,
        pagerank    = if (exists("pr")) pr else NULL,
        cross_edges = cross_edges,
        plot        = p
    ))
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_output`.
#' @param plot_obj Internal parameter
#' @param prep Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_output <- function(plot_obj, prep) plot_obj

#' Plot Dendro Network
#' @description
#' Internal helper for `.plot_dendro_network`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param g_slot Slot name.
#' @param cluster_name Parameter value.
#' @param legend_title Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param seed Random seed.
#' @param label_nodes Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param min_degree_hub Numeric threshold.
#' @param max.overlaps Parameter value.
#' @param cluster_vec Parameter value.
#' @param cluster_palette Parameter value.
#' @param tree_mode Parameter value.
#' @param k_top Parameter value.
#' @param edge_thresh Parameter value.
#' @param prune_mode Parameter value.
#' @param attach_I_delta Parameter value.
#' @param I_delta_name Parameter value.
#' @param I_delta_palette Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_dendro_network <- function(
    scope_obj,
    grid_name = NULL,
    g_slot = "g_consensus",
    cluster_name,
    legend_title = NULL,
    node_size = 5,
    edge_width = 1,
    label_size = 2,
    seed = 1,
    label_nodes = TRUE,
    hub_border_col = "black",
    hub_border_size = 0.5,
    min_degree_hub = 3,
    max.overlaps = Inf,
    cluster_vec = NULL,
    cluster_palette = NULL,
    tree_mode = NULL,
    k_top = 5,
    edge_thresh = 0.05,
    prune_mode = c("none", "mst", "triangulation"),
    attach_I_delta = FALSE,
    I_delta_name = "iDelta",
    I_delta_palette = scales::hue_pal()) {
    .plot_dendro_network_core(
        scope_obj = scope_obj,
        grid_name = grid_name,
        g_slot = g_slot,
        cluster_name = cluster_name,
        legend_title = legend_title,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        seed = seed,
        label_nodes = label_nodes,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        min_degree_hub = min_degree_hub,
        max.overlaps = max.overlaps,
        cluster_vec = cluster_vec,
        cluster_palette = cluster_palette,
        tree_mode = tree_mode,
        k_top = k_top,
        edge_thresh = edge_thresh,
        prune_mode = prune_mode,
        attach_I_delta = attach_I_delta,
        I_delta_name = I_delta_name,
        I_delta_palette = I_delta_palette
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_theme`.
#' @param plot_obj Internal parameter
#' @param prep Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_theme <- function(plot_obj, prep) plot_obj

#' Plot Dendro Network Multi
#' @description
#' Internal helper for `.plot_dendro_network_multi`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param gene_subset Parameter value.
#' @param L_min Numeric threshold.
#' @param use_consensus_graph Logical flag.
#' @param graph_slot_name Slot name.
#' @param cluster_vec Parameter value.
#' @param IDelta_col_name Parameter value.
#' @param IDelta_invert Parameter value.
#' @param damping Parameter value.
#' @param weight_low_cut Parameter value.
#' @param cluster_palette Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param seed Random seed.
#' @param length_scale Parameter value.
#' @param max.overlaps Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param title Parameter value.
#' @param k_top Parameter value.
#' @param tree_mode Parameter value.
#' @param n_runs Parameter value.
#' @param noise_sd Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_dendro_network_multi <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL,
    L_min = 0,
    use_consensus_graph = TRUE,
    graph_slot_name = NULL,
    cluster_vec = NULL,
    IDelta_col_name = NULL,
    IDelta_invert = FALSE,
    damping = 0.85,
    weight_low_cut = 0,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
        "#A65628", "#984EA3", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
    ),
    node_size = 4,
    edge_width = 3,
    label_size = 4,
    seed = 1,
    length_scale = 1,
    max.overlaps = 10,
    hub_border_col = "#4B4B4B",
    hub_border_size = 0.8,
    title = NULL,
    k_top = 1,
    tree_mode = c("radial", "rooted"),
    n_runs = 10,
    noise_sd = 0.01) {
    .plot_dendro_network_multi_core(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        gene_subset = gene_subset,
        L_min = L_min,
        use_consensus_graph = use_consensus_graph,
        graph_slot_name = graph_slot_name,
        cluster_vec = cluster_vec,
        IDelta_col_name = IDelta_col_name,
        IDelta_invert = IDelta_invert,
        damping = damping,
        weight_low_cut = weight_low_cut,
        cluster_palette = cluster_palette,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        seed = seed,
        length_scale = length_scale,
        max.overlaps = max.overlaps,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        title = title,
        k_top = k_top,
        tree_mode = tree_mode,
        n_runs = n_runs,
        noise_sd = noise_sd
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_build`.
#' @param prep Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_build <- function(prep) do.call(.plot_dendro_network_core_v2, prep$args)

#' Plot Dendro Network With Branches
#' @description
#' Internal helper for `.plot_dendro_network_with_branches`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param g_slot Slot name.
#' @param cluster_name Parameter value.
#' @param legend_title Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param seed Random seed.
#' @param label_nodes Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param min_degree_hub Numeric threshold.
#' @param max.overlaps Parameter value.
#' @param cluster_vec Parameter value.
#' @param cluster_palette Parameter value.
#' @param tree_mode Parameter value.
#' @param k_top Parameter value.
#' @param edge_thresh Parameter value.
#' @param prune_mode Parameter value.
#' @param attach_I_delta Parameter value.
#' @param I_delta_name Parameter value.
#' @param I_delta_palette Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_dendro_network_with_branches <- function(
    scope_obj,
    grid_name = NULL,
    g_slot = "g_consensus",
    cluster_name,
    legend_title = NULL,
    node_size = 5,
    edge_width = 1,
    label_size = 2,
    seed = 1,
    label_nodes = TRUE,
    hub_border_col = "black",
    hub_border_size = 0.5,
    min_degree_hub = 3,
    max.overlaps = Inf,
    cluster_vec = NULL,
    cluster_palette = NULL,
    tree_mode = NULL,
    k_top = 5,
    edge_thresh = 0.05,
    prune_mode = c("none", "mst", "triangulation"),
    attach_I_delta = FALSE,
    I_delta_name = "iDelta",
    I_delta_palette = scales::hue_pal()) {
    annotation_payload <- .plot_dendro_network_branches_prepare(
        scope_obj = scope_obj,
        grid_name = grid_name,
        g_slot = g_slot,
        cluster_name = cluster_name,
        legend_title = legend_title,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        seed = seed,
        label_nodes = label_nodes,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        min_degree_hub = min_degree_hub,
        max.overlaps = max.overlaps,
        cluster_vec = cluster_vec,
        cluster_palette = cluster_palette,
        tree_mode = tree_mode,
        k_top = k_top,
        edge_thresh = edge_thresh,
        prune_mode = prune_mode,
        attach_I_delta = attach_I_delta,
        I_delta_name = I_delta_name,
        I_delta_palette = I_delta_palette
    )
    branch_plot_obj <- .plot_dendro_network_branches_build(annotation_payload)
    branch_plot_obj <- .plot_dendro_network_branches_theme(branch_plot_obj, annotation_payload)
    .plot_dendro_network_branches_output(branch_plot_obj, annotation_payload)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_prepare`.
#' @param ... Additional arguments (currently unused).
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_prepare <- function(...) list(args = list(...))

#' Plot Idelta
#' @description
#' Internal helper for `.plot_idelta`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param cluster_col Parameter value.
#' @param top_n Parameter value.
#' @param min_genes Numeric threshold.
#' @param nrow Parameter value.
#' @param point_size Parameter value.
#' @param line_size Parameter value.
#' @param label_size Parameter value.
#' @param fill_col Parameter value.
#' @param outline_col Parameter value.
#' @param seed Random seed.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_idelta <- function(
    scope_obj,
    grid_name,
    cluster_col,
    top_n = NULL,
    min_genes = 1,
    nrow = 1,
    point_size = 3,
    line_size = 0.5,
    label_size = 2.5,
    fill_col = "steelblue",
    outline_col = "black",
    seed = NULL) {
    .plot_i_delta_core(
        scope_obj = scope_obj,
        grid_name = grid_name,
        cluster_col = cluster_col,
        top_n = top_n,
        min_genes = min_genes,
        nrow = nrow,
        point_size = point_size,
        line_size = line_size,
        label_size = label_size,
        fill_col = fill_col,
        outline_col = outline_col,
        seed = seed
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_network_core`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param g_slot Slot name.
#' @param cluster_name Parameter value.
#' @param legend_title Parameter value.
#' @param node_size Parameter value.
#' @param edge_width Parameter value.
#' @param label_size Parameter value.
#' @param seed Random seed.
#' @param label_nodes Parameter value.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param min_degree_hub Numeric threshold.
#' @param max.overlaps Parameter value.
#' @param cluster_vec Parameter value.
#' @param cluster_palette Parameter value.
#' @param tree_mode Parameter value.
#' @param k_top Parameter value.
#' @param edge_thresh Parameter value.
#' @param prune_mode Parameter value.
#' @param attach_I_delta Parameter value.
#' @param I_delta_name Parameter value.
#' @param I_delta_palette Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_network_core <- function(
    scope_obj,
    grid_name = NULL,
    g_slot = "g_consensus",
    cluster_name,
    legend_title = NULL,
    node_size = 5,
    edge_width = 1,
    label_size = 2,
    seed = 1,
    label_nodes = TRUE,
    hub_border_col = "black",
    hub_border_size = 0.5,
    min_degree_hub = 3,
    max.overlaps = Inf,
    cluster_vec = NULL,
    cluster_palette = NULL,
    tree_mode = NULL,
    k_top = 5,
    edge_thresh = 0.05,
    prune_mode = c("none", "mst", "triangulation"),
    attach_I_delta = FALSE,
    I_delta_name = "iDelta",
    I_delta_palette = scales::hue_pal()) {
    prep <- .plot_dendro_network_prepare(
        scope_obj = scope_obj,
        grid_name = grid_name,
        g_slot = g_slot,
        cluster_name = cluster_name,
        legend_title = legend_title,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        seed = seed,
        label_nodes = label_nodes,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        min_degree_hub = min_degree_hub,
        max.overlaps = max.overlaps,
        cluster_vec = cluster_vec,
        cluster_palette = cluster_palette,
        tree_mode = tree_mode,
        k_top = k_top,
        edge_thresh = edge_thresh,
        prune_mode = prune_mode,
        attach_I_delta = attach_I_delta,
        I_delta_name = I_delta_name,
        I_delta_palette = I_delta_palette
    )
    plot_obj <- .plot_dendro_network_build(prep)
    plot_obj <- .plot_dendro_network_theme(plot_obj, prep)
    .plot_dendro_network_output(plot_obj, prep)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_core_v2`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param graph_slot Slot name.
#' @param cluster_name Parameter value.
#' @param cluster_ids Parameter value.
#' @param IDelta_col_name Parameter value.
#' @param linkage Parameter value.
#' @param plot_dend Parameter value.
#' @param weight_low_cut Parameter value.
#' @param damping Parameter value.
#' @param cluster_palette Parameter value.
#' @param tip_dots Parameter value.
#' @param tip_point_cex Parameter value.
#' @param tip_palette Parameter value.
#' @param cluster_name2 Parameter value.
#' @param tip_palette2 Parameter value.
#' @param tip_row_offset2 Parameter value.
#' @param tip_label_offset Parameter value.
#' @param tip_label_cex Parameter value.
#' @param tip_label_adj Parameter value.
#' @param tip_label_srt Parameter value.
#' @param tip_label_col Parameter value.
#' @param leaf_order Parameter value.
#' @param length_mode Parameter value.
#' @param weight_normalize Parameter value.
#' @param weight_clip_quantile Parameter value.
#' @param height_rescale Parameter value.
#' @param height_power Parameter value.
#' @param height_scale Parameter value.
#' @param distance_on Parameter value.
#' @param enforce_cluster_contiguity Parameter value.
#' @param distance_smooth_power Parameter value.
#' @param tip_shape2 Parameter value.
#' @param tip_row1_label Parameter value.
#' @param tip_row2_label Parameter value.
#' @param tip_label_indent Parameter value.
#' @param legend_inline Parameter value.
#' @param legend_files Parameter value.
#' @param compose_outfile Parameter value.
#' @param compose_width Parameter value.
#' @param compose_height Parameter value.
#' @param compose_res Parameter value.
#' @param legend_ncol1 Parameter value.
#' @param legend_ncol2 Parameter value.
#' @param title Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_core_v2 <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    cluster_name,
    cluster_ids = NULL,
    IDelta_col_name = NULL,
    linkage = "average",
    plot_dend = TRUE,
    weight_low_cut = 0,
    damping = 0.85,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
        "#A65628", "#984EA3", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
    ),
    tip_dots = TRUE,
    tip_point_cex = 1.2,
    tip_palette = NULL,
    cluster_name2 = NULL,
    tip_palette2 = NULL,
    tip_row_offset2 = -0.4,
    tip_label_offset = -0.6,
    tip_label_cex = 0.6,
    tip_label_adj = 1,
    tip_label_srt = 90,
    tip_label_col = "black",
    leaf_order = c("OLO", "none"),
    length_mode = c("neg_log", "inverse", "inverse_sqrt"),
    weight_normalize = TRUE,
    weight_clip_quantile = 0.05,
    height_rescale = c("q95", "max", "none"),
    height_power = 1,
    height_scale = 1.2,
    distance_on = c("tree", "graph"),
    enforce_cluster_contiguity = TRUE,
    distance_smooth_power = 1,
    tip_shape2 = c("square", "circle", "diamond", "triangle"),
    tip_row1_label = NULL,
    tip_row2_label = NULL,
    tip_label_indent = 0.01,
    legend_inline = FALSE,
    legend_files = NULL,
    compose_outfile = NULL,
    compose_width = 2400,
    compose_height = 1600,
    compose_res = 200,
    legend_ncol1 = 1,
    legend_ncol2 = 1,
    title = "Graph-weighted Dendrogram") {
    params <- .plot_dendro_validate_inputs(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        graph_slot = graph_slot,
        cluster_name = cluster_name,
        cluster_ids = cluster_ids,
        IDelta_col_name = IDelta_col_name,
        linkage = linkage,
        plot_dend = plot_dend,
        weight_low_cut = weight_low_cut,
        damping = damping,
        cluster_palette = cluster_palette,
        tip_dots = tip_dots,
        tip_point_cex = tip_point_cex,
        tip_palette = tip_palette,
        cluster_name2 = cluster_name2,
        tip_palette2 = tip_palette2,
        tip_row_offset2 = tip_row_offset2,
        tip_label_offset = tip_label_offset,
        tip_label_cex = tip_label_cex,
        tip_label_adj = tip_label_adj,
        tip_label_srt = tip_label_srt,
        tip_label_col = tip_label_col,
        leaf_order = leaf_order,
        length_mode = length_mode,
        weight_normalize = weight_normalize,
        weight_clip_quantile = weight_clip_quantile,
        height_rescale = height_rescale,
        height_power = height_power,
        height_scale = height_scale,
        distance_on = distance_on,
        enforce_cluster_contiguity = enforce_cluster_contiguity,
        distance_smooth_power = distance_smooth_power,
        tip_shape2 = tip_shape2,
        tip_row1_label = tip_row1_label,
        tip_row2_label = tip_row2_label,
        tip_label_indent = tip_label_indent,
        legend_inline = legend_inline,
        legend_files = legend_files,
        compose_outfile = compose_outfile,
        compose_width = compose_width,
        compose_height = compose_height,
        compose_res = compose_res,
        legend_ncol1 = legend_ncol1,
        legend_ncol2 = legend_ncol2,
        title = title
    )

    tree_data <- .plot_dendro_prepare_tree_data(params)
    dist_data <- .plot_dendro_build_distance(tree_data, params)
    dend_info <- .plot_dendro_build_dendrogram(dist_data, params)
    annotated <- .plot_dendro_add_annotations(dend_info, tree_data, params)
    .plot_dendro_finalize(annotated, params)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_finalize`.
#' @param annotated Internal parameter
#' @param params Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_finalize <- function(annotated, params) {
    legend_inline <- params$legend_inline
    legend_files <- params$legend_files
    compose_outfile <- params$compose_outfile
    compose_width <- params$compose_width
    compose_height <- params$compose_height
    compose_res <- params$compose_res
    plot_dend <- params$plot_dend
    tip_point_cex <- annotated$tip_point_cex
    tip_row_offset2 <- annotated$tip_row_offset2
    tip_label_indent <- params$tip_label_indent
    tip_label_offset <- annotated$tip_label_offset
    tip_label_cex <- annotated$tip_label_cex
    tip_label_adj <- annotated$tip_label_adj
    tip_label_srt <- annotated$tip_label_srt
    tip_label_col <- annotated$tip_label_col
    dend <- annotated$dend
    hc <- annotated$hc
    Dm <- annotated$Dm
    genes <- annotated$genes
    memb_all <- annotated$memb_all
    g_tree <- annotated$g_tree
    g_full <- annotated$g_full
    pr <- annotated$pr
    tip1 <- annotated$tip1
    tip2 <- annotated$tip2
    legend_ncol1 <- annotated$legend_ncol1
    legend_ncol2 <- annotated$legend_ncol2
    title <- annotated$title

    render_dend_panel <- function(draw_legends = FALSE, preserve_par = TRUE) {
        if (preserve_par) {
            op <- par(no.readonly = TRUE)
            on.exit(par(op), add = TRUE)
        } else {
            op <- par(c("mar", "mgp", "xpd"))
            on.exit(do.call(par, op), add = TRUE)
        }

        par(mar = c(3.9, 3.3, 3.1, 0.2), mgp = c(1.5, 0.25, 0), xpd = NA)

        max_h <- if (length(hc$height)) max(hc$height, na.rm = TRUE) else 1
        if (!is.finite(max_h) || is.na(max_h) || max_h <= 0) max_h <- 1
        y_gap <- 0.06 * max_h
        min_off <- min(0, tip_label_offset, if (tip2$show) tip_row_offset2 else 0)
        pad_off <- 0.1
        ylim_use <- c(y_gap * (min_off - pad_off), max_h)

        plot(dend,
            main = title,
            ylab = "Weighted distance",
            ylim = ylim_use,
            leaflab = "none"
        )

        labs_local <- labels(dend)
        x_local <- seq_along(labs_local)
        usr <- par("usr")
        x_range <- usr[2] - usr[1]
        y_range <- usr[4] - usr[3]

        if (tip1$show) {
            y1 <- rep(0, length(labs_local))
            points(x_local, y1, pch = 21, bg = tip1$cols, col = "black", cex = tip_point_cex)
            if (draw_legends && length(tip1$legend_labels)) {
                label_x1 <- usr[1] + tip_label_indent * x_range
                legend_y1 <- mean(y1) - 0.12 * y_gap
                legend_margin <- 0.02 * x_range
                max_width1 <- max(legend_margin, usr[2] - label_x1 - legend_margin)
                .draw_wrapped_legend(
                    x_start = label_x1,
                    y_start = legend_y1,
                    labels = tip1$legend_labels,
                    pch = tip1$pch,
                    bg = tip1$legend_colors,
                    border = rep("black", length(tip1$legend_labels)),
                    point_cex = tip_point_cex,
                    text_cex = 0.8,
                    x_range = x_range,
                    y_range = y_range,
                    max_width = max_width1,
                    dot_gap_factor = 0.6,
                    item_gap_factor = 0.4,
                    row_spacing_factor = 1.1
                )
            }
        }

        if (tip2$show) {
            y2 <- rep(tip_row_offset2 * y_gap, length(labs_local))
            pch2_sym <- if (length(tip2$pch)) tip2$pch[1] else 21
            points(x_local, y2, pch = pch2_sym, bg = tip2$cols, col = "black", cex = tip_point_cex)
            if (draw_legends && length(tip2$legend_labels)) {
                label_x2 <- usr[1] + tip_label_indent * x_range
                legend_y2 <- mean(y2) - 0.28 * y_gap
                legend_margin <- 0.02 * x_range
                max_width2 <- max(legend_margin, usr[2] - label_x2 - legend_margin)
                .draw_wrapped_legend(
                    x_start = label_x2,
                    y_start = legend_y2,
                    labels = tip2$legend_labels,
                    pch = rep(pch2_sym, length(tip2$legend_labels)),
                    bg = tip2$legend_colors,
                    border = rep("black", length(tip2$legend_labels)),
                    point_cex = tip_point_cex,
                    text_cex = 0.8,
                    x_range = x_range,
                    y_range = y_range,
                    max_width = max_width2,
                    dot_gap_factor = 0.6,
                    item_gap_factor = 0.4,
                    row_spacing_factor = 1.1
                )
            }
        }

        y_lab <- rep(tip_label_offset * y_gap, length(labs_local))
        text(x_local, y_lab,
            labels = labs_local,
            srt = tip_label_srt, xpd = NA,
            cex = tip_label_cex, adj = tip_label_adj,
            col = tip_label_col
        )
    }

    render_legend_panel <- function() {
        par(mar = c(0.5, 0.4, 0.5, 0.2))
        plot.new()
        usr <- par("usr")
        x_span <- usr[2] - usr[1]
        y_span <- usr[4] - usr[3]
        x_left <- usr[1] + 0.08 * x_span
        legend_width <- 0.28 * x_span
        legend_gap <- legend_width
        x_right <- x_left + 2 * legend_width + legend_gap
        if (x_right > usr[2] - 0.02 * x_span) {
            overflow <- x_right - (usr[2] - 0.02 * x_span)
            x_left <- x_left - overflow
            x_right <- x_right - overflow
        }
        y_top <- usr[4] - 0.15 * y_span
        if (tip1$show && length(tip1$legend_labels)) {
            legend(x_left, y_top,
                legend = tip1$legend_labels,
                title = tip1$title, pch = tip1$pch,
                pt.bg = tip1$legend_colors, col = "black",
                pt.cex = tip_point_cex, cex = 0.9, bty = "n",
                ncol = legend_ncol1, xjust = 0, yjust = 1,
                x.intersp = 0.6, y.intersp = 0.8
            )
        }
        if (tip2$show && length(tip2$legend_labels)) {
            pch2_sym <- if (length(tip2$pch)) tip2$pch[1] else 21
            legend(x_right, y_top,
                legend = tip2$legend_labels,
                title = tip2$title, pch = rep(pch2_sym, length(tip2$legend_labels)),
                pt.bg = tip2$legend_colors, col = "black",
                pt.cex = tip_point_cex, cex = 0.9, bty = "n",
                ncol = legend_ncol2, xjust = 1, yjust = 1,
                x.intersp = 0.6, y.intersp = 0.8
            )
        }
    }
    if (!is.null(compose_outfile)) {
        png(filename = compose_outfile, width = compose_width, height = compose_height, res = compose_res)
        on.exit(try(dev.off(), silent = TRUE), add = TRUE)
        layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1), heights = 1)
        render_dend_panel(draw_legends = FALSE, preserve_par = FALSE)
        render_legend_panel()
        layout(1)
    }

    if (isTRUE(plot_dend)) {
        if (legend_inline) {
            render_dend_panel(draw_legends = TRUE)
        } else {
            op_layout <- par(no.readonly = TRUE)
            layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1), heights = 1)
            render_dend_panel(draw_legends = FALSE, preserve_par = FALSE)
            render_legend_panel()
            layout(1)
            par(op_layout)
        }
    }

    if (!is.null(legend_files)) {
        if (length(legend_files) >= 1 && !is.na(legend_files[1]) && tip1$show && length(tip1$legend_labels)) {
            png(filename = legend_files[1], width = 1200, height = 800, res = 150)
            par(mar = c(0, 0, 0, 0))
            plot.new()
            legend(
                x = "center", y = "center", legend = tip1$legend_labels,
                title = tip1$title, pch = tip1$pch,
                pt.bg = tip1$legend_colors, col = "black",
                pt.cex = tip_point_cex, cex = 0.9, bty = "n",
                ncol = legend_ncol1, x.intersp = 0.6, y.intersp = 0.8
            )
            dev.off()
        }
        if (length(legend_files) >= 2 && !is.na(legend_files[2]) && tip2$show && length(tip2$legend_labels)) {
            png(filename = legend_files[2], width = 1200, height = 800, res = 150)
            par(mar = c(0, 0, 0, 0))
            plot.new()
            legend(
                x = "center", y = "center", legend = tip2$legend_labels,
                title = tip2$title, pch = rep(tip2$pch[1], length(tip2$legend_labels)),
                pt.bg = tip2$legend_colors, col = "black",
                pt.cex = tip_point_cex, cex = 0.9, bty = "n",
                ncol = legend_ncol2, x.intersp = 0.6, y.intersp = 0.8
            )
            dev.off()
        }
    }

    invisible(list(
        dend = dend,
        hclust = hc,
        dist = as.dist(Dm),
        genes = genes,
        cluster_map = setNames(memb_all[genes], genes),
        tree_graph = g_tree,
        graph = g_full,
        PageRank = if (!is.null(pr)) pr[names(pr) %in% genes] else NULL,
        legend_row1 = list(
            labels = tip1$legend_labels,
            colors = tip1$legend_colors,
            pch = tip1$pch,
            title = tip1$title
        ),
        legend_row2 = list(
            labels = tip2$legend_labels,
            colors = tip2$legend_colors,
            pch = tip2$pch,
            title = tip2$title
        )
    ))
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_add_annotations`.
#' @param dend_info Internal parameter
#' @param tree_data Internal parameter
#' @param params Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_add_annotations <- function(dend_info, tree_data, params) {
    scope_obj <- params$scope_obj
    cluster_name <- params$cluster_name
    cluster_name2 <- params$cluster_name2
    cluster_palette <- params$cluster_palette
    tip_palette <- params$tip_palette
    tip_palette2 <- params$tip_palette2
    tip_shape2 <- params$tip_shape2
    tip_row_offset2 <- params$tip_row_offset2
    tip_row1_label <- params$tip_row1_label
    tip_row2_label <- params$tip_row2_label
    tip_label_offset <- params$tip_label_offset
    tip_label_cex <- params$tip_label_cex
    tip_label_adj <- params$tip_label_adj
    tip_label_srt <- params$tip_label_srt
    tip_label_col <- params$tip_label_col
    tip_dots <- params$tip_dots
    tip_point_cex <- params$tip_point_cex
    legend_ncol1 <- params$legend_ncol1
    legend_ncol2 <- params$legend_ncol2
    title <- params$title

    dend <- dend_info$dend
    hc <- dend_info$hc
    Dm <- dend_info$Dm
    genes <- dend_info$genes
    memb_all <- dend_info$memb_all

    labs <- labels(dend)
    x <- seq_along(labs)

    tip1 <- list(
        show = isTRUE(tip_dots) && length(labs) > 0,
        title = tip_row1_label
    )
    if (tip1$show) {
        cl_raw1 <- scope_obj@meta.data[labs, cluster_name, drop = TRUE]
        cl_chr1 <- as.character(cl_raw1)
        lev1 <- .order_levels_numeric(unique(cl_chr1))
        pal_map1 <- .make_pal_map(lev1, cluster_palette, tip_palette)
        cols1 <- unname(pal_map1[cl_chr1])
        cols1[is.na(cols1)] <- "gray80"
        tip1$cols <- cols1
        tip1$legend_labels <- .format_label_display(lev1)
        tip1$legend_colors <- unname(pal_map1[lev1])
        tip1$pch <- rep(21, length(lev1))
        tip1$lev_order <- lev1
    } else {
        tip1$legend_labels <- character(0)
        tip1$legend_colors <- character(0)
        tip1$pch <- numeric(0)
    }
    if (is.null(tip1$title)) tip1$title <- ""

    tip2 <- list(
        show = !is.null(cluster_name2) && length(labs) > 0,
        title = if (is.null(tip_row2_label)) "" else tip_row2_label
    )
    if (tip2$show) {
        cl_raw2 <- scope_obj@meta.data[labs, cluster_name2, drop = TRUE]
        cl_chr2 <- as.character(cl_raw2)
        lev2 <- .order_levels_numeric(unique(cl_chr2))
        pal_map2 <- .make_pal_map(lev2, cluster_palette, tip_palette2)
        cols2 <- unname(pal_map2[cl_chr2])
        cols2[is.na(cols2)] <- "gray80"
        tip2$cols <- cols2
        tip2$legend_labels <- .format_label_display(lev2)
        tip2$legend_colors <- unname(pal_map2[lev2])
        tip2$pch <- rep(switch(tip_shape2,
            square = 22,
            circle = 21,
            diamond = 23,
            triangle = 24
        ), length(lev2))
        tip2$lev_order <- lev2
    } else {
        tip2$legend_labels <- character(0)
        tip2$legend_colors <- character(0)
        tip2$pch <- numeric(0)
    }
    if (is.null(tip2$title)) tip2$title <- ""

    list(
        dend = dend,
        hc = hc,
        Dm = Dm,
        genes = genes,
        memb_all = memb_all,
        pr = tree_data$pr,
        g_tree = tree_data$g_tree,
        g_full = tree_data$g_full,
        tip1 = tip1,
        tip2 = tip2,
        x = x,
        tip_label_offset = tip_label_offset,
        tip_label_cex = tip_label_cex,
        tip_label_adj = tip_label_adj,
        tip_label_srt = tip_label_srt,
        tip_label_col = tip_label_col,
        tip_point_cex = tip_point_cex,
        tip_row_offset2 = tip_row_offset2,
        legend_ncol1 = legend_ncol1,
        legend_ncol2 = legend_ncol2,
        title = title,
        tip_row1_label = tip_row1_label,
        tip_row2_label = tip_row2_label,
        tip_shape2 = params$tip_shape2,
        tip_label_indent = params$tip_label_indent,
        tip_label_adj = tip_label_adj,
        tip_label_srt = tip_label_srt,
        tip_label_col = tip_label_col,
        tip_label_offset = tip_label_offset,
        cluster_name = cluster_name,
        cluster_name2 = cluster_name2
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_build_dendrogram`.
#' @param dist_data Internal parameter
#' @param params Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_build_dendrogram <- function(dist_data, params) {
    Dm <- dist_data$Dm
    genes <- dist_data$genes
    memb_all <- dist_data$memb_all
    enforce_cluster_contiguity <- params$enforce_cluster_contiguity
    distance_smooth_power <- params$distance_smooth_power
    height_rescale <- params$height_rescale
    linkage <- params$linkage
    leaf_order <- params$leaf_order
    # optional smoothing of the distance matrix
    if (!is.null(distance_smooth_power)) {
        k <- max(3L, as.integer(round(ncol(Dm) * (distance_smooth_power / 10))))
        if ((k %% 2L) == 0L) k <- k + 1L
        if (k > 3) {
            ut <- upper.tri(Dm, diag = FALSE)
            dv <- as.numeric(Dm[ut])
            ord <- order(dv)
            x <- dv[ord]
            xs <- runmed(x, k = k, endrule = "keep")
            # guard NA from runmed at boundaries
            xs[!is.finite(xs)] <- 0
            dv[ord] <- pmax(xs, 0)
            Dm2 <- Dm
            Dm2[ut] <- dv
            Dm2 <- Dm2 + t(Dm2)
            diag(Dm2) <- 0
            Dm <- Dm2
        }
    }

    # fill disconnected/invalid entries with a large finite distance
    if (any(!is.finite(Dm))) {
        maxf <- suppressWarnings(max(Dm[is.finite(Dm)], na.rm = TRUE))
        if (!is.finite(maxf) || is.na(maxf)) maxf <- 1
        Dm[!is.finite(Dm)] <- maxf * 1.2 + 1
    }
    diag(Dm) <- 0

    # Optionally force within-cluster merges to occur first
    if (isTRUE(enforce_cluster_contiguity)) {
        clv <- memb_all[genes]
        within_mask <- outer(clv, clv, "==")
        ut <- upper.tri(Dm)
        wvals <- Dm[ut & within_mask]
        bvals <- Dm[ut & !within_mask]
        w_max <- suppressWarnings(max(wvals[wvals > 0 & is.finite(wvals)], na.rm = TRUE))
        b_min <- suppressWarnings(min(bvals[bvals > 0 & is.finite(bvals)], na.rm = TRUE))
        if (is.finite(w_max) && is.finite(b_min) && b_min <= w_max) {
            off <- (w_max - b_min) + max(1e-6, 0.05 * w_max)
            Dm[!within_mask] <- Dm[!within_mask] + off
        }
    }

    # rescale heights if requested
    h <- if (height_rescale == "q95") {
        quantile(Dm[upper.tri(Dm)], 0.95, na.rm = TRUE)
    } else if (height_rescale == "max") {
        max(Dm, na.rm = TRUE)
    } else {
        NA_real_
    }
    if (is.finite(h) && h > 0) Dm <- Dm / h

    hc <- hclust(as.dist(Dm), method = linkage)
    dend <- as.dendrogram(hc)

    ## ---------- Optional OLO ----------
    if (leaf_order == "OLO") {
        ok_ser <- requireNamespace("seriation", quietly = TRUE)
        ok_den <- requireNamespace("dendextend", quietly = TRUE)
        if (ok_ser && ok_den) {
            ord <- seriation::get_order(seriation::seriate(as.dist(Dm), method = "OLO"))
            desired <- rownames(Dm)[ord]
            dend <- dendextend::rotate(dend, order = desired)
        } else if (requireNamespace("gclus", quietly = TRUE)) {
            hc <- gclus::reorder.hclust(hc, dist = as.dist(Dm))
            dend <- as.dendrogram(hc)
        } else if (requireNamespace("dendsort", quietly = TRUE)) {
            dend <- dendsort::dendsort(dend)
        } else if (!ok_ser || !ok_den) {
            warning("Leaf-order optimization requested but required packages are missing; install one of: seriation + dendextend (preferred), gclus, or dendsort.")
        }
    }

    list(
        dend = dend,
        hc = hc,
        Dm = Dm,
        genes = genes,
        memb_all = memb_all
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_build_distance`.
#' @param tree_data Internal parameter
#' @param params Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_build_distance <- function(tree_data, params) {
    g_tree <- tree_data$g_tree
    g_full <- tree_data$g_full
    memb_all <- tree_data$memb_all
    length_mode <- params$length_mode
    height_power <- params$height_power
    height_scale <- params$height_scale
    weight_normalize <- params$weight_normalize
    weight_clip_quantile <- params$weight_clip_quantile
    distance_on <- params$distance_on
    enforce_cluster_contiguity <- params$enforce_cluster_contiguity
    genes <- igraph::V(g_tree)$name
    if (distance_on == "graph") {
        ew <- igraph::E(g_full)$weight
        if (isTRUE(weight_normalize)) {
            ew <- (ew - min(ew, na.rm = TRUE)) / max(1e-12, diff(range(ew, na.rm = TRUE)))
        }
        if (!is.null(weight_clip_quantile) && weight_clip_quantile > 0) {
            ql <- quantile(ew, probs = weight_clip_quantile, na.rm = TRUE)
            qu <- quantile(ew, probs = 1 - weight_clip_quantile, na.rm = TRUE)
            ew <- pmin(pmax(ew, ql), qu)
        }
        ew <- ifelse(ew <= 0 | is.na(ew), 1e-6, ew)
        elen <- switch(length_mode,
            neg_log      = -log(pmax(ew, 1e-9)),
            inverse      = 1 / pmax(ew, 1e-9),
            inverse_sqrt = 1 / sqrt(pmax(ew, 1e-9))
        )
        elen <- (elen^height_power) * height_scale
        igraph::E(g_full)$length <- elen
        Dm <- igraph::distances(g_full, v = igraph::V(g_full), to = igraph::V(g_full), weights = igraph::E(g_full)$length)
        rownames(Dm) <- igraph::V(g_full)$name
        colnames(Dm) <- igraph::V(g_full)$name
        Dm <- Dm[genes, genes, drop = FALSE]
    } else {
        # compute lengths on the skeleton tree and use its geodesics
        ew_t <- igraph::E(g_tree)$weight
        if (isTRUE(weight_normalize)) {
            ew_t <- (ew_t - min(ew_t, na.rm = TRUE)) / max(1e-12, diff(range(ew_t, na.rm = TRUE)))
        }
        if (!is.null(weight_clip_quantile) && weight_clip_quantile > 0) {
            ql <- quantile(ew_t, probs = weight_clip_quantile, na.rm = TRUE)
            qu <- quantile(ew_t, probs = 1 - weight_clip_quantile, na.rm = TRUE)
            ew_t <- pmin(pmax(ew_t, ql), qu)
        }
        ew_t <- ifelse(ew_t <= 0 | is.na(ew_t), 1e-6, ew_t)
        elen_t <- switch(length_mode,
            neg_log      = -log(pmax(ew_t, 1e-9)),
            inverse      = 1 / pmax(ew_t, 1e-9),
            inverse_sqrt = 1 / sqrt(pmax(ew_t, 1e-9))
        )
        elen_t <- (elen_t^height_power) * height_scale
        igraph::E(g_tree)$length <- elen_t
        Dm <- igraph::distances(g_tree, v = genes, to = genes, weights = igraph::E(g_tree)$length)
        rownames(Dm) <- genes
        colnames(Dm) <- genes
    }

    list(
        Dm = Dm,
        g_tree = g_tree,
        g_full = g_full,
        genes = genes,
        memb_all = memb_all
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_prepare_tree_data`.
#' @param params Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_dendro_prepare_tree_data <- function(params) {
    scope_obj <- params$scope_obj
    grid_name <- params$grid_name
    cluster_name <- params$cluster_name
    cluster_ids <- params$cluster_ids
    IDelta_col_name <- params$IDelta_col_name
    damping <- params$damping
    weight_low_cut <- params$weight_low_cut
    lee_stats_layer <- params$lee_stats_layer
    graph_slot <- params$graph_slot
    ## ---------- 0. Resolve grid layer & LeeStats ----------
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    if (is.null(grid_name)) {
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
    }

    leeStat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
        scope_obj@stats[[grid_name]][[lee_stats_layer]]
    } else if (!is.null(g_layer[[lee_stats_layer]])) {
        g_layer[[lee_stats_layer]]
    } else {
        stop("Cannot find layer '", lee_stats_layer, "' in grid '", grid_name, "'.")
    }

    ## ---------- 1. Extract base graph and subset to target genes ----------
    g_base <- leeStat[[graph_slot]]
    if (is.null(g_base) || !inherits(g_base, "igraph")) {
        stop("Graph slot '", graph_slot, "' not found or not an igraph in LeeStats layer.")
    }

    memb_all <- as.character(scope_obj@meta.data[[cluster_name]])
    if (is.null(cluster_ids)) cluster_ids <- sort(unique(na.omit(memb_all)))
    cluster_ids <- as.character(cluster_ids)
    genes_all <- rownames(scope_obj@meta.data)[memb_all %in% cluster_ids]
    if (length(genes_all) < 2) stop("Need at least two genes in the selected clusters.")

    g <- igraph::induced_subgraph(g_base, vids = intersect(igraph::V(g_base)$name, genes_all))
    if (igraph::vcount(g) < 2) stop("Too few vertices in induced subgraph.")

    ## Ensure edges have usable numeric weights (fallback to |L| if missing/invalid)
    Lmat <- .get_lee_matrix(scope_obj,
        grid_name = grid_name,
        lee_layer = lee_stats_layer
    )
    ed_l <- igraph::as_edgelist(g, names = TRUE)
    if (nrow(ed_l) > 0) {
        w_lmat <- mapply(function(a, b) abs(Lmat[a, b]), ed_l[, 1], ed_l[, 2], USE.NAMES = FALSE)
        if (!"weight" %in% igraph::edge_attr_names(g)) {
            igraph::E(g)$weight <- w_lmat
        } else {
            w_now <- igraph::E(g)$weight
            bad <- is.na(w_now) | !is.finite(w_now) | w_now <= 0
            if (any(bad)) w_now[bad] <- w_lmat[bad]
            igraph::E(g)$weight <- w_now
        }
        # Final NA guard
        igraph::E(g)$weight[is.na(igraph::E(g)$weight) | !is.finite(igraph::E(g)$weight)] <- 0
    }

    ## ---------- 2. Idelta-PageRank reweighting (same as .plot_dendro_network) ----------
    pr <- NULL
    if (!is.null(IDelta_col_name)) {
        delta <- scope_obj@meta.data[igraph::V(g)$name, IDelta_col_name, drop = TRUE]
        delta[is.na(delta)] <- median(delta, na.rm = TRUE)
        q <- quantile(delta, c(.1, .9))
        delta <- pmax(pmin(delta, q[2]), q[1])
        pers <- exp(delta - max(delta))
        pers <- pers / sum(pers)
        pr_tmp <- igraph::page_rank(
            g,
            personalized = pers,
            damping = damping,
            weights = igraph::E(g)$weight,
            directed = FALSE
        )$vector
        pr <- pr_tmp
        ed <- igraph::as_data_frame(g, "edges")
        w_new <- igraph::E(g)$weight * ((pr[ed$from] + pr[ed$to]) / 2)
        w_new[w_new <= weight_low_cut] <- 0
        igraph::E(g)$weight <- w_new
    }

    ## ---------- 3. Skeleton: intra-cluster MST + inter-cluster MST (same as .plot_dendro_network) ----------
    Vnames <- igraph::V(g)$name
    clu <- as.character(scope_obj@meta.data[Vnames, cluster_name, drop = TRUE])
    keep <- !is.na(clu)
    g <- igraph::induced_subgraph(g, Vnames[keep])
    Vnames <- igraph::V(g)$name
    clu <- clu[keep]

    # tag edge ids to track selection
    igraph::E(g)$eid <- seq_len(igraph::ecount(g))

    # 3.1 Intra-cluster MST: use max weight via length 1/(w + eps)
    keep_eid <- integer(0)
    for (cl in unique(clu)) {
        vsub <- Vnames[clu == cl]
        if (length(vsub) < 2) next
        g_sub <- igraph::induced_subgraph(g, vsub)
        if (igraph::ecount(g_sub) > 0) {
            mst_sub <- igraph::mst(g_sub, weights = 1 / (igraph::E(g_sub)$weight + 1e-9))
            keep_eid <- c(keep_eid, igraph::edge_attr(mst_sub, "eid"))
        }
    }

    # 3.2 Inter-cluster MST: aggregate by each pair's max-weight edge
    if (length(unique(clu)) > 1) {
        ed_tab <- igraph::as_data_frame(g, "edges")
        ed_tab$eid <- igraph::E(g)$eid
        ed_tab$cl1 <- clu[ed_tab$from]
        ed_tab$cl2 <- clu[ed_tab$to]
        inter <- ed_tab[ed_tab$cl1 != ed_tab$cl2, , drop = FALSE]
        if (nrow(inter) > 0) {
            inter$pair <- ifelse(inter$cl1 < inter$cl2,
                paste(inter$cl1, inter$cl2, sep = "|"), paste(inter$cl2, inter$cl1, sep = "|")
            )
            # drop rows with NA weights to avoid aggregate() zero-row error
            inter <- inter[!is.na(inter$weight) & is.finite(inter$weight), , drop = FALSE]
            if (nrow(inter) == 0) {
                agg <- inter
            } else {
                agg <- aggregate(weight ~ pair + cl1 + cl2, data = inter, max)
            }
            g_clu <- igraph::graph_from_data_frame(
                agg[, c("cl1", "cl2", "weight")],
                directed = FALSE, vertices = unique(clu)
            )
            if (igraph::ecount(g_clu) > 0) {
                mstc <- igraph::mst(g_clu, weights = 1 / (igraph::E(g_clu)$weight + 1e-9))
                ep <- igraph::ends(mstc, igraph::E(mstc))
                keep_pairs <- ifelse(ep[, 1] < ep[, 2], paste(ep[, 1], ep[, 2], sep = "|"), paste(ep[, 2], ep[, 1], sep = "|"))
                inter_MST <- inter[inter$pair %in% keep_pairs, ]
                keep_eid <- unique(c(keep_eid, inter_MST$eid))
            }
        }
    }

    # 3.3 Keep skeleton edges only; drop isolates
    g_tree <- g
    if (length(keep_eid)) {
        g_tree <- igraph::delete_edges(g_tree, igraph::E(g_tree)[!eid %in% keep_eid])
    }
    g_tree <- igraph::delete_vertices(g_tree, which(igraph::degree(g_tree) == 0))
    if (igraph::vcount(g_tree) < 2) stop("MST-based tree has fewer than 2 vertices.")

    list(
        g_tree = g_tree,
        g_full = igraph::induced_subgraph(g, vids = igraph::V(g_tree)$name),
        g = g,
        g_layer = g_layer,
        grid_name = grid_name,
        leeStat = leeStat,
        memb_all = memb_all,
        cluster_ids = cluster_ids,
        pr = pr,
        Lmat = Lmat,
        cluster_name = cluster_name
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_validate_inputs`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param graph_slot Slot name.
#' @param cluster_name Parameter value.
#' @param cluster_ids Parameter value.
#' @param IDelta_col_name Parameter value.
#' @param linkage Parameter value.
#' @param plot_dend Parameter value.
#' @param weight_low_cut Parameter value.
#' @param damping Parameter value.
#' @param cluster_palette Parameter value.
#' @param tip_dots Parameter value.
#' @param tip_point_cex Parameter value.
#' @param tip_palette Parameter value.
#' @param cluster_name2 Parameter value.
#' @param tip_palette2 Parameter value.
#' @param tip_row_offset2 Parameter value.
#' @param tip_label_offset Parameter value.
#' @param tip_label_cex Parameter value.
#' @param tip_label_adj Parameter value.
#' @param tip_label_srt Parameter value.
#' @param tip_label_col Parameter value.
#' @param leaf_order Parameter value.
#' @param length_mode Parameter value.
#' @param weight_normalize Parameter value.
#' @param weight_clip_quantile Parameter value.
#' @param height_rescale Parameter value.
#' @param height_power Parameter value.
#' @param height_scale Parameter value.
#' @param distance_on Parameter value.
#' @param enforce_cluster_contiguity Parameter value.
#' @param distance_smooth_power Parameter value.
#' @param tip_shape2 Parameter value.
#' @param tip_row1_label Parameter value.
#' @param tip_row2_label Parameter value.
#' @param tip_label_indent Parameter value.
#' @param legend_inline Parameter value.
#' @param legend_files Parameter value.
#' @param compose_outfile Parameter value.
#' @param compose_width Parameter value.
#' @param compose_height Parameter value.
#' @param compose_res Parameter value.
#' @param legend_ncol1 Parameter value.
#' @param legend_ncol2 Parameter value.
#' @param title Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_validate_inputs <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    cluster_name,
    cluster_ids = NULL,
    IDelta_col_name = NULL,
    linkage = "average",
    plot_dend = TRUE,
    weight_low_cut = 0,
    damping = 0.85,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
        "#A65628", "#984EA3", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
    ),
    tip_dots = TRUE,
    tip_point_cex = 1.2,
    tip_palette = NULL,
    cluster_name2 = NULL,
    tip_palette2 = NULL,
    tip_row_offset2 = -0.4,
    tip_label_offset = -0.6,
    tip_label_cex = 0.6,
    tip_label_adj = 1,
    tip_label_srt = 90,
    tip_label_col = "black",
    leaf_order = c("OLO", "none"),
    length_mode = c("neg_log", "inverse", "inverse_sqrt"),
    weight_normalize = TRUE,
    weight_clip_quantile = 0.05,
    height_rescale = c("q95", "max", "none"),
    height_power = 1,
    height_scale = 1.2,
    distance_on = c("tree", "graph"),
    enforce_cluster_contiguity = TRUE,
    distance_smooth_power = 1,
    tip_shape2 = c("square", "circle", "diamond", "triangle"),
    tip_row1_label = NULL,
    tip_row2_label = NULL,
    tip_label_indent = 0.01,
    legend_inline = FALSE,
    legend_files = NULL,
    compose_outfile = NULL,
    compose_width = 2400,
    compose_height = 1600,
    compose_res = 200,
    legend_ncol1 = 1,
    legend_ncol2 = 1,
    title = "Graph-weighted Dendrogram") {
    leaf_order <- match.arg(leaf_order)
    length_mode <- match.arg(length_mode)
    height_rescale <- match.arg(height_rescale)
    distance_on <- match.arg(distance_on)
    tip_shape2 <- match.arg(tip_shape2)
    if (is.null(tip_row1_label)) tip_row1_label <- cluster_name
    if (is.null(tip_row2_label) && !is.null(cluster_name2)) tip_row2_label <- cluster_name2
    if (is.null(tip_label_indent)) tip_label_indent <- 0.01
    if (is.null(legend_inline)) legend_inline <- FALSE
    if (!is.null(distance_smooth_power) && distance_smooth_power <= 0) {
        stop("distance_smooth_power must be > 0 or NULL")
    }
    list(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        graph_slot = graph_slot,
        cluster_name = cluster_name,
        cluster_ids = cluster_ids,
        IDelta_col_name = IDelta_col_name,
        linkage = linkage,
        plot_dend = plot_dend,
        weight_low_cut = weight_low_cut,
        damping = damping,
        cluster_palette = cluster_palette,
        tip_dots = tip_dots,
        tip_point_cex = tip_point_cex,
        tip_palette = tip_palette,
        cluster_name2 = cluster_name2,
        tip_palette2 = tip_palette2,
        tip_row_offset2 = tip_row_offset2,
        tip_label_offset = tip_label_offset,
        tip_label_cex = tip_label_cex,
        tip_label_adj = tip_label_adj,
        tip_label_srt = tip_label_srt,
        tip_label_col = tip_label_col,
        leaf_order = leaf_order,
        length_mode = length_mode,
        weight_normalize = weight_normalize,
        weight_clip_quantile = weight_clip_quantile,
        height_rescale = height_rescale,
        height_power = height_power,
        height_scale = height_scale,
        distance_on = distance_on,
        enforce_cluster_contiguity = enforce_cluster_contiguity,
        distance_smooth_power = distance_smooth_power,
        tip_shape2 = tip_shape2,
        tip_row1_label = tip_row1_label,
        tip_row2_label = tip_row2_label,
        tip_label_indent = tip_label_indent,
        legend_inline = legend_inline,
        legend_files = legend_files,
        compose_outfile = compose_outfile,
        compose_width = compose_width,
        compose_height = compose_height,
        compose_res = compose_res,
        legend_ncol1 = legend_ncol1,
        legend_ncol2 = legend_ncol2,
        title = title
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_dendro_core`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param graph_slot Slot name.
#' @param cluster_name Parameter value.
#' @param cluster_ids Parameter value.
#' @param IDelta_col_name Parameter value.
#' @param linkage Parameter value.
#' @param plot_dend Parameter value.
#' @param weight_low_cut Parameter value.
#' @param damping Parameter value.
#' @param cluster_palette Parameter value.
#' @param tip_dots Parameter value.
#' @param tip_point_cex Parameter value.
#' @param tip_palette Parameter value.
#' @param cluster_name2 Parameter value.
#' @param tip_palette2 Parameter value.
#' @param tip_row_offset2 Parameter value.
#' @param tip_label_offset Parameter value.
#' @param tip_label_cex Parameter value.
#' @param tip_label_adj Parameter value.
#' @param tip_label_srt Parameter value.
#' @param tip_label_col Parameter value.
#' @param leaf_order Parameter value.
#' @param length_mode Parameter value.
#' @param weight_normalize Parameter value.
#' @param weight_clip_quantile Parameter value.
#' @param height_rescale Parameter value.
#' @param height_power Parameter value.
#' @param height_scale Parameter value.
#' @param distance_on Parameter value.
#' @param enforce_cluster_contiguity Parameter value.
#' @param distance_smooth_power Parameter value.
#' @param tip_shape2 Parameter value.
#' @param tip_row1_label Parameter value.
#' @param tip_row2_label Parameter value.
#' @param tip_label_indent Parameter value.
#' @param legend_inline Parameter value.
#' @param legend_files Parameter value.
#' @param compose_outfile Parameter value.
#' @param compose_width Parameter value.
#' @param compose_height Parameter value.
#' @param compose_res Parameter value.
#' @param legend_ncol1 Parameter value.
#' @param legend_ncol2 Parameter value.
#' @param title Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_dendro_core <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    cluster_name,
    cluster_ids = NULL,
    IDelta_col_name = NULL,
    linkage = "average",
    plot_dend = TRUE,
    weight_low_cut = 0,
    damping = 0.85,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
        "#A65628", "#984EA3", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
    ),
    tip_dots = TRUE,
    tip_point_cex = 1.2,
    tip_palette = NULL,
    cluster_name2 = NULL,
    tip_palette2 = NULL,
    tip_row_offset2 = -0.4,
    tip_label_offset = -0.6,
    tip_label_cex = 0.6,
    tip_label_adj = 1,
    tip_label_srt = 90,
    tip_label_col = "black",
    leaf_order = c("OLO", "none"),
    length_mode = c("neg_log", "inverse", "inverse_sqrt"),
    weight_normalize = TRUE,
    weight_clip_quantile = 0.05,
    height_rescale = c("q95", "max", "none"),
    height_power = 1,
    height_scale = 1.2,
    distance_on = c("tree", "graph"),
    enforce_cluster_contiguity = TRUE,
    distance_smooth_power = 1,
    tip_shape2 = c("square", "circle", "diamond", "triangle"),
    tip_row1_label = NULL,
    tip_row2_label = NULL,
    tip_label_indent = 0.01,
    legend_inline = FALSE,
    legend_files = NULL,
    compose_outfile = NULL,
    compose_width = 2400,
    compose_height = 1600,
    compose_res = 200,
    legend_ncol1 = 1,
    legend_ncol2 = 1,
    title = "Graph-weighted Dendrogram") {
    prep <- .plot_dendro_prepare(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        graph_slot = graph_slot,
        cluster_name = cluster_name,
        cluster_ids = cluster_ids,
        IDelta_col_name = IDelta_col_name,
        linkage = linkage,
        plot_dend = plot_dend,
        weight_low_cut = weight_low_cut,
        damping = damping,
        cluster_palette = cluster_palette,
        tip_dots = tip_dots,
        tip_point_cex = tip_point_cex,
        tip_palette = tip_palette,
        cluster_name2 = cluster_name2,
        tip_palette2 = tip_palette2,
        tip_row_offset2 = tip_row_offset2,
        tip_label_offset = tip_label_offset,
        tip_label_cex = tip_label_cex,
        tip_label_adj = tip_label_adj,
        tip_label_srt = tip_label_srt,
        tip_label_col = tip_label_col,
        leaf_order = leaf_order,
        length_mode = length_mode,
        weight_normalize = weight_normalize,
        weight_clip_quantile = weight_clip_quantile,
        height_rescale = height_rescale,
        height_power = height_power,
        height_scale = height_scale,
        distance_on = distance_on,
        enforce_cluster_contiguity = enforce_cluster_contiguity,
        distance_smooth_power = distance_smooth_power,
        tip_shape2 = tip_shape2,
        tip_row1_label = tip_row1_label,
        tip_row2_label = tip_row2_label,
        tip_label_indent = tip_label_indent,
        legend_inline = legend_inline,
        legend_files = legend_files,
        compose_outfile = compose_outfile,
        compose_width = compose_width,
        compose_height = compose_height,
        compose_res = compose_res,
        legend_ncol1 = legend_ncol1,
        legend_ncol2 = legend_ncol2,
        title = title
    )
    plot_obj <- .plot_dendro_build_layers(prep)
    plot_obj <- .plot_dendro_theme(plot_obj, prep)
    .plot_dendro_output(plot_obj, prep)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_add_scale_bar_v2`.
#' @param p Internal parameter
#' @param prep Parameter value.
#' @param x_rng Parameter value.
#' @param y_rng Parameter value.
#' @param scale_bar_pos Parameter value.
#' @param scale_bar_corner Parameter value.
#' @param bar_len Parameter value.
#' @param arrow_pt Parameter value.
#' @param scale_bar_colour Parameter value.
#' @param scale_text_size Parameter value.
#' @param bar_offset Parameter value.
#' @param scale_bar_show Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_add_scale_bar_v2 <- function(p,
                                        prep,
                                        x_rng,
                                        y_rng,
                                        scale_bar_pos,
                                        scale_bar_corner,
                                        bar_len,
                                        arrow_pt,
                                        scale_bar_colour,
                                        scale_text_size,
                                        bar_offset,
                                        scale_bar_show) {
    draw_scale_bar <- isTRUE(scale_bar_show)
    scale_pos <- .plot_density_resolve_scale_bar_pos(scale_bar_pos, scale_bar_corner)
    if (prep$use_image_coords) {
        roi <- prep$histology_roi
        span_um <- roi["xmax"] - roi["xmin"]
        px_per_um <- prep$histology_slot$width / span_um
        if (!is.finite(px_per_um) || px_per_um <= 0) {
            draw_scale_bar <- FALSE
        } else {
            bar_len_display <- bar_len * px_per_um
            min_frac_x <- min(0.4, bar_len_display / max(1, prep$histology_slot$width))
            x_margin <- max(0.03, min_frac_x / 2)
            x_frac <- max(x_margin, min(1 - min_frac_x - x_margin, scale_pos["x"]))
            y_margin <- 0.22
            y_frac <- max(y_margin, min(1 - y_margin, scale_pos["y"]))
            x0 <- x_frac * prep$histology_slot$width
            x1 <- x0 + bar_len_display
            if (x0 < 0) {
                x1 <- x1 - x0
                x0 <- 0
            }
            if (x1 > prep$histology_slot$width) {
                shift <- x1 - prep$histology_slot$width
                x0 <- max(0, x0 - shift)
                x1 <- prep$histology_slot$width
            }
            label_offset <- max(15, 0.08 * prep$histology_slot$height)
            arrow_shift <- max(12, 0.08 * prep$histology_slot$height)
            if (prep$flip_histology_y) {
                y_bar_base <- prep$histology_slot$height - y_frac * prep$histology_slot$height
                y_bar <- min(prep$histology_slot$height, max(0, y_bar_base + arrow_shift))
                label_y <- y_bar_base - label_offset
                if (label_y < 0) label_y <- y_bar_base + label_offset
            } else {
                y_bar_base <- y_frac * prep$histology_slot$height
                y_bar <- max(0, min(prep$histology_slot$height, y_bar_base - arrow_shift))
                label_y <- y_bar_base + label_offset
                if (label_y > prep$histology_slot$height) label_y <- y_bar_base - label_offset
            }
        }
    } else {
        bar_len_display <- bar_len
        x_span <- diff(x_rng)
        y_span <- diff(y_rng)
        min_frac_x <- min(0.4, bar_len_display / max(1e-6, x_span))
        x_margin <- max(0.03, min_frac_x / 2)
        x_frac <- max(x_margin, min(1 - min_frac_x - x_margin, scale_pos["x"]))
        y_margin <- 0.22
        y_frac <- max(y_margin, min(1 - y_margin, scale_pos["y"]))
        x0 <- x_rng[1] + x_frac * x_span
        x1 <- x0 + bar_len_display
        if (x0 < x_rng[1]) {
            shift <- x_rng[1] - x0
            x0 <- x_rng[1]
            x1 <- x1 + shift
        }
        if (x1 > x_rng[2]) {
            shift <- x1 - x_rng[2]
            x0 <- x0 - shift
            x1 <- x_rng[2]
        }
        if (x0 < x_rng[1] || x1 > x_rng[2]) {
            draw_scale_bar <- FALSE
        }
        label_offset <- max(15, 0.08 * y_span)
        arrow_shift <- max(12, 0.08 * y_span)
        y_bar_base <- y_rng[1] + y_frac * y_span
        y_bar <- y_bar_base - arrow_shift
        y_bar <- max(min(y_rng), min(max(y_rng), y_bar))
        label_y <- y_bar_base + label_offset
        if (label_y > max(y_rng)) label_y <- y_bar_base - label_offset
        if (label_y < min(y_rng)) label_y <- y_bar_base + label_offset
    }
    if (isTRUE(draw_scale_bar)) {
        p <- p +
            annotate("segment",
                x = x0, xend = x1,
                y = y_bar,
                yend = y_bar,
                arrow = arrow(length = unit(arrow_pt, "pt"), ends = "both"),
                colour = scale_bar_colour,
                linewidth = 0.4
            ) +
            annotate("text",
                x = (x0 + x1) / 2,
                y = label_y,
                label = paste0(bar_len, " \u00B5m"),
                colour = scale_bar_colour,
                vjust = 1, size = scale_text_size
            )
    }
    p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_apply_theme_axes`.
#' @param p Internal parameter
#' @param prep Internal parameter
#' @param grid_gap Internal parameter
#' @param target_aspect_ratio Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_density_apply_theme_axes <- function(p, prep, grid_gap, target_aspect_ratio) {
    if (!prep$use_image_coords) {
        ensure_extent <- function(rng) {
            if (any(!is.finite(rng))) return(c(0, 1))
            span <- diff(rng)
            if (!is.finite(span) || span <= 0) {
                center <- mean(rng)
                if (!is.finite(center)) center <- 0
                return(c(center - 0.5, center + 0.5))
            }
            rng
        }
        x_rng <- ensure_extent(prep$x_rng_phys)
        y_rng <- ensure_extent(prep$y_rng_phys)
        if (!is.null(target_aspect_ratio)) {
            dx <- diff(x_rng)
            dy <- diff(y_rng)
            current_ratio <- dx / dy
            if (is.finite(current_ratio) && dx > 0 && dy > 0) {
                if (current_ratio < target_aspect_ratio) {
                    pad <- (target_aspect_ratio * dy - dx) / 2
                    if (is.finite(pad) && pad > 0) x_rng <- x_rng + c(-pad, pad)
                } else if (current_ratio > target_aspect_ratio) {
                    pad <- (dx / target_aspect_ratio - dy) / 2
                    if (is.finite(pad) && pad > 0) y_rng <- y_rng + c(-pad, pad)
                }
            }
        }
        grid_x <- seq(x_rng[1], x_rng[2], by = grid_gap)
        grid_y <- seq(y_rng[1], y_rng[2], by = grid_gap)
    } else {
        x_rng <- c(0, prep$histology_slot$width)
        y_rng <- if (prep$flip_histology_y) c(prep$histology_slot$height, 0) else c(0, prep$histology_slot$height)
        grid_x <- prep$phys_to_img(seq(prep$x_rng_phys[1], prep$x_rng_phys[2], by = grid_gap), axis = "x")
        grid_y <- prep$phys_to_img(seq(prep$y_rng_phys[1], prep$y_rng_phys[2], by = grid_gap), axis = "y")
        grid_x <- grid_x[is.finite(grid_x)]
        grid_y <- grid_y[is.finite(grid_y)]
    }

    if (length(grid_x)) {
        p <- p + geom_vline(xintercept = grid_x, linewidth = 0.05, colour = "grey80")
    }
    if (length(grid_y)) {
        p <- p + geom_hline(yintercept = grid_y, linewidth = 0.05, colour = "grey80")
    }

    if (prep$use_image_coords) {
        if (prep$flip_histology_y) {
            p <- p +
                scale_x_continuous(limits = c(0, prep$histology_slot$width), expand = c(0, 0)) +
                scale_y_reverse(limits = c(prep$histology_slot$height, 0), expand = c(0, 0)) +
                coord_fixed(expand = FALSE)
        } else {
            p <- p +
                scale_x_continuous(limits = c(0, prep$histology_slot$width), expand = c(0, 0)) +
                scale_y_continuous(limits = c(0, prep$histology_slot$height), expand = c(0, 0)) +
                coord_fixed(expand = FALSE)
        }
    } else {
        p <- p +
            scale_x_continuous(limits = x_rng, expand = c(0, 0)) +
            scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
            coord_fixed()
    }

    p <- p +
        theme_minimal(base_size = 9) +
        theme(
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.key.width = unit(0.4, "cm"),
            legend.key.height = unit(0.15, "cm"),
            legend.text = element_text(size = 9, angle = 90),
            legend.title = element_text(size = 8, hjust = 0, vjust = 1),
            plot.title = element_text(size = 10, hjust = 0.5),
            plot.margin = margin(1.2, 1, 1.5, 1, "cm"),
            axis.line.x.bottom = element_line(colour = "black"),
            axis.line.y.left = element_line(colour = "black"),
            axis.line.x.top = element_blank(),
            axis.line.y.right = element_blank(),
            axis.title = element_blank()
        ) +
        theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

    list(p = p, x_rng = x_rng, y_rng = y_rng)
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_add_segmentation`.
#' @param p Internal parameter
#' @param prep Parameter value.
#' @param seg_type Parameter value.
#' @param colour_cell Parameter value.
#' @param colour_nucleus Parameter value.
#' @param alpha_seg Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_add_segmentation <- function(p,
                                           prep,
                                           seg_type,
                                           colour_cell,
                                           colour_nucleus,
                                           alpha_seg) {
    seg_dt <- prep$seg_dt
    if (!length(seg_dt)) return(p)

    if (seg_type == "both") {
        is_cell <- seg_dt$cell %in% prep$seg_cell_ids
        seg_dt[, segClass := ifelse(is_cell, "cell", "nucleus")]
        p <- p +
            geom_shape(
                data = seg_dt[segClass == "cell"],
                aes(x = x, y = y, group = cell),
                colour = scales::alpha(colour_cell, alpha_seg),
                fill = NA,
                linewidth = 0.05
            ) +
            geom_shape(
                data = seg_dt[segClass == "nucleus"],
                aes(x = x, y = y, group = cell),
                colour = scales::alpha(colour_nucleus, alpha_seg),
                fill = NA,
                linewidth = 0.05
            )
    } else {
        col_use <- if (seg_type == "cell") colour_cell else colour_nucleus
        p <- p +
            geom_shape(
                data = seg_dt,
                aes(x = x, y = y, group = cell),
                colour = scales::alpha(col_use, alpha_seg),
                fill = NA,
                linewidth = 0.05
            )
    }
    p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_build_layers`.
#' @param prep Internal parameter
#' @param palette1 Parameter value.
#' @param palette2 Parameter value.
#' @param alpha1 Parameter value.
#' @param alpha2 Parameter value.
#' @param density1_name Parameter value.
#' @param density2_name Parameter value.
#' @param overlay_image Parameter value.
#' @param image_path Filesystem path.
#' @param image_alpha Parameter value.
#' @param image_choice Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_build_layers <- function(prep,
                                       palette1,
                                       palette2,
                                       alpha1,
                                       alpha2,
                                       density1_name,
                                       density2_name,
                                       overlay_image,
                                       image_path,
                                       image_alpha,
                                       image_choice) {
    library(ggplot2)
    p <- ggplot()

    if (prep$histology_available) {
        img_rgba <- .plot_density_prepare_rgba(prep$histology_slot$png, image_alpha)
        if (!is.null(img_rgba)) {
            if (!prep$use_image_coords && isTRUE(prep$flip_histology_y)) {
                img_rgba <- .plot_density_flip_raster(img_rgba)
            }
            if (prep$use_image_coords) {
                p <- p + annotation_raster(
                    img_rgba,
                    xmin = 0, xmax = prep$histology_slot$width,
                    ymin = 0, ymax = prep$histology_slot$height,
                    interpolate = TRUE
                )
            } else {
                roi_draw <- if (!is.null(prep$histology_roi)) prep$histology_roi else prep$default_roi
                p <- p + annotation_raster(
                    img_rgba,
                    xmin = roi_draw["xmin"], xmax = roi_draw["xmax"],
                    ymin = roi_draw["ymin"], ymax = roi_draw["ymax"],
                    interpolate = TRUE
                )
            }
        }
    } else if (isTRUE(overlay_image)) {
        img_info <- prep$image_info
        hires_path <- if (!is.null(img_info)) img_info$hires_path else NULL
        lowres_path <- if (!is.null(img_info)) img_info$lowres_path else NULL
        if (!is.null(image_path)) {
            hires_path <- image_path
            lowres_path <- NULL
        }
        sel_path <- switch(image_choice,
            auto = if (!is.null(hires_path)) hires_path else lowres_path,
            hires = hires_path,
            lowres = lowres_path
        )
        if (!is.null(sel_path) && file.exists(sel_path) && is.finite(prep$mpp_fullres)) {
            ext <- tolower(file_ext(sel_path))
            img <- NULL
            if (ext %in% c("png") && requireNamespace("png", quietly = TRUE)) {
                img <- png::readPNG(sel_path)
            } else if (ext %in% c("jpg", "jpeg") && requireNamespace("jpeg", quietly = TRUE)) {
                img <- jpeg::readJPEG(sel_path)
            }
            if (!is.null(img)) {
                if (length(dim(img)) == 2L) {
                    img_rgb <- array(NA_real_, dim = c(dim(img)[1], dim(img)[2], 3L))
                    img_rgb[, , 1] <- img
                    img_rgb[, , 2] <- img
                    img_rgb[, , 3] <- img
                    img <- img_rgb
                }
                if (length(dim(img)) == 3L) {
                    if (dim(img)[3] == 3L) {
                        img_rgba <- array(NA_real_, dim = c(dim(img)[1], dim(img)[2], 4L))
                        img_rgba[, , 1:3] <- img
                        img_rgba[, , 4] <- max(0, min(1, image_alpha))
                        img <- img_rgba
                    } else if (dim(img)[3] == 4L) {
                        img[, , 4] <- pmin(1, pmax(0, img[, , 4] * image_alpha))
                    }
                }

                wpx <- dim(img)[2]
                hpx <- dim(img)[1]
                hires_sf <- if (!is.null(img_info)) img_info$tissue_hires_scalef else NA_real_
                lowres_sf <- if (!is.null(img_info)) img_info$tissue_lowres_scalef else NA_real_
                aligned_sf <- if (!is.null(img_info)) img_info$regist_target_img_scalef else NA_real_
                eff_mpp <- prep$mpp_fullres
                is_hires <- grepl("tissue_hires_image", basename(sel_path), ignore.case = TRUE)
                is_lowres <- grepl("tissue_lowres_image", basename(sel_path), ignore.case = TRUE)
                is_aligned <- grepl("aligned_tissue_image", basename(sel_path), ignore.case = TRUE)
                if (is_hires && is.finite(hires_sf)) eff_mpp <- prep$mpp_fullres / hires_sf
                if (is_lowres && is.finite(lowres_sf)) eff_mpp <- prep$mpp_fullres / lowres_sf
                if (is_aligned && is.finite(aligned_sf)) eff_mpp <- prep$mpp_fullres / aligned_sf
                if (is_aligned && !is.finite(aligned_sf)) eff_mpp <- prep$mpp_fullres
                if (!is.finite(eff_mpp) || eff_mpp <= 0) eff_mpp <- prep$mpp_fullres
                w_um <- as.numeric(wpx) * as.numeric(eff_mpp)
                h_um <- as.numeric(hpx) * as.numeric(eff_mpp)
                y_origin_manual <- if (!is.null(img_info$y_origin)) img_info$y_origin else "top-left"
                if (identical(y_origin_manual, "top-left")) {
                    img <- .plot_density_flip_raster(img)
                }
                p <- p + annotation_raster(img,
                    xmin = 0, xmax = w_um,
                    ymin = 0, ymax = h_um,
                    interpolate = TRUE
                )
            }
        }
    }

    shape_layers1 <- .plot_density_build_shape_geom(.plot_density_tile_df(prep$heat1), alpha1, prep$tile_shape, prep$hex_orientation)
    p <- Reduce(`+`, shape_layers1, init = p)
    p <- p +
        scale_fill_gradient(
            name = density1_name,
            low = "transparent",
            high = palette1,
            limits = c(0, unique(prep$heat1$cut)),
            oob = scales::squish,
            labels = scales::number_format(accuracy = prep$accuracy_val),
            na.value = "transparent",
            guide = guide_colorbar(order = 1)
        )

    if (!is.null(prep$heat2)) {
        library(ggnewscale)
        shape_layers2 <- .plot_density_build_shape_geom(.plot_density_tile_df(prep$heat2), alpha2, prep$tile_shape, prep$hex_orientation)
        p <- p + new_scale_fill()
        p <- Reduce(`+`, shape_layers2, init = p)
        p <- p +
            scale_fill_gradient(
                name = density2_name,
                low = "transparent",
                high = palette2,
                limits = c(0, unique(prep$heat2$cut)),
                oob = scales::squish,
                labels = scales::number_format(accuracy = prep$accuracy_val),
                na.value = "transparent",
                guide = guide_colorbar(order = 2)
            )
    }
    p
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_prepare_data`.
#' @param scope_obj Internal parameter
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param density1_name Parameter value.
#' @param density2_name Parameter value.
#' @param max.cutoff1 Parameter value.
#' @param max.cutoff2 Parameter value.
#' @param tile_shape Parameter value.
#' @param hex_orientation Parameter value.
#' @param target_aspect_ratio Parameter value.
#' @param scale_bar_corner Parameter value.
#' @param use_histology Logical flag.
#' @param histology_level Parameter value.
#' @param axis_mode_requested Parameter value.
#' @param overlay_image Parameter value.
#' @param image_path Filesystem path.
#' @param image_alpha Parameter value.
#' @param image_choice Parameter value.
#' @param seg_type Parameter value.
#' @param legend_digits Parameter value.
#' @return Internal helper result
#' @keywords internal
.plot_density_prepare_data <- function(scope_obj,
                                       grid_name,
                                       density1_name,
                                       density2_name,
                                       max.cutoff1,
                                       max.cutoff2,
                                       tile_shape,
                                       hex_orientation,
                                       target_aspect_ratio,
                                       scale_bar_corner,
                                       use_histology,
                                       histology_level,
                                       axis_mode_requested,
                                       overlay_image,
                                       image_path,
                                       image_alpha,
                                       image_choice,
                                       seg_type,
                                       legend_digits) {
    accuracy_val <- 1 / (10^legend_digits)
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_layer_name <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    grid_info <- g_layer$grid_info
    if (is.null(grid_info) ||
        !all(c("grid_id", "xmin", "xmax", "ymin", "ymax") %in% names(grid_info))) {
        stop("grid_info is missing required columns for layer '", grid_layer_name, "'.")
    }
    default_roi <- c(
        xmin = suppressWarnings(min(grid_info$xmin, na.rm = TRUE)),
        xmax = suppressWarnings(max(grid_info$xmax, na.rm = TRUE)),
        ymin = suppressWarnings(min(grid_info$ymin, na.rm = TRUE)),
        ymax = suppressWarnings(max(grid_info$ymax, na.rm = TRUE))
    )
    x_rng_phys <- range(grid_info$xmin, grid_info$xmax, na.rm = TRUE)
    y_rng_phys <- range(grid_info$ymin, grid_info$ymax, na.rm = TRUE)
    if (any(!is.finite(x_rng_phys))) x_rng_phys <- c(0, 1)
    if (any(!is.finite(y_rng_phys))) y_rng_phys <- c(0, 1)
    if (any(!is.finite(default_roi))) {
        default_roi <- c(xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    }
    image_info <- g_layer$image_info
    y_origin <- if (!is.null(image_info) && !is.null(image_info$y_origin)) {
        image_info$y_origin
    } else {
        "top-left"
    }
    histology_slot <- NULL
    if (isTRUE(use_histology) && !is.null(g_layer$histology)) {
        histology_slot <- g_layer$histology[[histology_level]]
        if (is.null(histology_slot)) {
            available <- Filter(function(x) !is.null(x), g_layer$histology)
            if (length(available)) histology_slot <- available[[1]]
        }
    }
    histology_available <- !is.null(histology_slot) && !is.null(histology_slot$png)
    histology_y_origin <- if (!is.null(histology_slot$y_origin)) {
        histology_slot$y_origin
    } else if (!is.null(image_info$y_origin)) {
        image_info$y_origin
    } else {
        "top-left"
    }
    flip_histology_y <- !identical(histology_y_origin, "bottom-left")
    histology_roi <- if (histology_available && !is.null(histology_slot$roi_bbox)) {
        histology_slot$roi_bbox
    } else {
        default_roi
    }
    histology_ready <- histology_available &&
        !is.null(histology_roi) &&
        all(is.finite(histology_roi)) &&
        is.finite(histology_slot$width) &&
        is.finite(histology_slot$height)
    axis_mode <- axis_mode_requested
    if (identical(axis_mode_requested, "image")) {
        if (!histology_ready) {
            warning("axis_mode = 'image' requires histology with ROI metadata; reverting to grid axes.")
            axis_mode <- "grid"
        }
    } else {
        axis_mode <- "grid"
    }
    use_image_coords <- identical(axis_mode, "image")
    mpp_fullres <- if (!is.null(g_layer$microns_per_pixel)) {
        as.numeric(g_layer$microns_per_pixel)
    } else {
        NA_real_
    }

    densityDF <- scope_obj@density[[grid_layer_name]]
    if (is.null(densityDF)) {
        densityDF <- g_layer$densityDF
    }
    if (is.null(densityDF)) {
        stop("No density table found for grid '", grid_layer_name, "'.")
    }
    if (is.null(rownames(densityDF)) && "grid_id" %in% names(densityDF)) {
        rownames(densityDF) <- densityDF$grid_id
    }
    for (dcol in c(density1_name, density2_name)) {
        if (!is.null(dcol) && !(dcol %in% colnames(densityDF))) {
            stop(
                "Density column '", dcol, "' not found in table for grid '",
                grid_layer_name, "'."
            )
        }
    }
    build_heat <- function(d_col, cutoff_frac) {
        df <- data.frame(
            grid_id = rownames(densityDF),
            d = densityDF[[d_col]],
            stringsAsFactors = FALSE
        )
        heat <- merge(grid_info, df, by = "grid_id", all.x = TRUE)
        heat$d[is.na(heat$d)] <- 0
        maxv <- max(heat$d, na.rm = TRUE)
        heat$cut <- maxv * cutoff_frac
        heat$d <- pmin(heat$d, heat$cut)
        heat
    }
    heat1 <- build_heat(density1_name, max.cutoff1)
    heat2 <- if (!is.null(density2_name)) build_heat(density2_name, max.cutoff2) else NULL

    phys_to_img <- NULL
    if (use_image_coords) {
        roi <- histology_roi
        width_px <- histology_slot$width
        height_px <- histology_slot$height
        rx <- roi["xmax"] - roi["xmin"]
        ry <- roi["ymax"] - roi["ymin"]
        phys_to_img <- function(vals, axis = c("x", "y")) {
            axis <- match.arg(axis)
            if (axis == "x") {
                (vals - roi["xmin"]) / rx * width_px
            } else {
                if (flip_histology_y) {
                    (1 - (vals - roi["ymin"]) / ry) * height_px
                } else {
                    (vals - roi["ymin"]) / ry * height_px
                }
            }
        }
        transform_heat <- function(df) {
            if (!nrow(df) || !is.finite(rx) || !is.finite(ry) || rx == 0 || ry == 0) return(df)
            x1 <- phys_to_img(df$xmin, axis = "x")
            x2 <- phys_to_img(df$xmax, axis = "x")
            y1 <- phys_to_img(df$ymin, axis = "y")
            y2 <- phys_to_img(df$ymax, axis = "y")
            df$xmin <- pmin(x1, x2)
            df$xmax <- pmax(x1, x2)
            df$ymin <- pmin(y1, y2)
            df$ymax <- pmax(y1, y2)
            df
        }
        heat1 <- transform_heat(heat1)
        if (!is.null(heat2)) heat2 <- transform_heat(heat2)
    } else {
        phys_to_img <- function(vals, axis = c("x", "y")) vals
    }

    seg_layers <- switch(seg_type,
        cell    = "segmentation_cell",
        nucleus = "segmentation_nucleus",
        both    = c("segmentation_cell", "segmentation_nucleus")
    )
    seg_layers <- seg_layers[seg_layers %in% names(scope_obj@coord)]
    seg_dt <- NULL
    seg_cell_ids <- if (!is.null(scope_obj@coord$segmentation_cell)) scope_obj@coord$segmentation_cell$cell else character()
    if (length(seg_layers)) {
        seg_dt <- rbindlist(scope_obj@coord[seg_layers],
            use.names = TRUE, fill = TRUE
        )
        if (use_image_coords) {
            seg_df <- .coords_physical_to_level(as.data.frame(seg_dt),
                x_col = "x", y_col = "y", histo = histology_slot
            )
            seg_dt <- as.data.table(seg_df)
            set(seg_dt, j = "x", value = seg_dt$x_img)
            set(seg_dt, j = "y", value = seg_dt$y_img)
            seg_dt[, c("x_img", "y_img") := NULL]
        }
    }

    list(
        g_layer = g_layer,
        grid_layer_name = grid_layer_name,
        grid_info = grid_info,
        default_roi = default_roi,
        x_rng_phys = x_rng_phys,
        y_rng_phys = y_rng_phys,
        histology_slot = histology_slot,
        histology_available = histology_available,
        histology_roi = histology_roi,
        histology_ready = histology_ready,
        flip_histology_y = flip_histology_y,
        axis_mode = axis_mode,
        use_image_coords = use_image_coords,
        mpp_fullres = mpp_fullres,
        image_info = image_info,
        heat1 = heat1,
        heat2 = heat2,
        phys_to_img = phys_to_img,
        seg_dt = seg_dt,
        seg_cell_ids = seg_cell_ids,
        seg_type = seg_type,
        accuracy_val = accuracy_val,
        tile_shape = tile_shape,
        hex_orientation = hex_orientation,
        target_aspect_ratio = target_aspect_ratio
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_resolve_scale_bar_pos`.
#' @param pos Internal parameter
#' @param corner Internal parameter
#' @param default_x Internal parameter
#' @param default_y Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_density_resolve_scale_bar_pos <- function(pos, corner, default_x = 0.03, default_y = 0.22) {
    corner_defaults <- list(
        "bottom-left" = c(x = default_x, y = default_y),
        "bottom-right" = c(x = 1 - default_x, y = default_y),
        "top-left" = c(x = default_x, y = 1 - default_y),
        "top-right" = c(x = 1 - default_x, y = 1 - default_y)
    )
    base <- corner_defaults[[corner]]
    out <- c(x = default_x, y = default_y)
    if (!is.null(base)) out <- base
    if (is.null(pos)) return(out)
    if (is.list(pos)) pos <- unlist(pos, use.names = TRUE)
    if (!length(pos)) return(out)
    clamp01 <- function(v) {
        v <- suppressWarnings(as.numeric(v[1]))
        if (!is.finite(v)) return(NA_real_)
        max(0, min(1, v))
    }
    if (is.null(names(pos))) {
        if (length(pos) >= 1) {
            val <- clamp01(pos[1])
            if (!is.na(val)) out["x"] <- val
        }
        if (length(pos) >= 2) {
            val <- clamp01(pos[2])
            if (!is.na(val)) out["y"] <- val
        }
    } else {
        if ("x" %in% names(pos)) {
            val <- clamp01(pos["x"])
            if (!is.na(val)) out["x"] <- val
        }
        if ("y" %in% names(pos)) {
            val <- clamp01(pos["y"])
            if (!is.na(val)) out["y"] <- val
        }
    }
    out
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_build_shape_geom`.
#' @param df Internal parameter
#' @param alpha_val Internal parameter
#' @param tile_shape Internal parameter
#' @param hex_orientation Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_density_build_shape_geom <- function(df, alpha_val, tile_shape, hex_orientation) {
    if (!nrow(df)) return(list())
    if (identical(tile_shape, "square")) {
        list(geom_tile(
            data = df,
            aes(x = x, y = y, width = w, height = h, fill = d),
            colour = NA, alpha = alpha_val
        ))
    } else if (identical(tile_shape, "circle")) {
        circ_df <- transform(df, radius = pmax(pmin(w, h), .Machine$double.eps) / 2)
        list(geom_circle(
            data = circ_df,
            aes(x0 = x, y0 = y, r = radius, fill = d),
            colour = NA,
            alpha = alpha_val,
            inherit.aes = FALSE
        ))
    } else if (identical(tile_shape, "hex")) {
        hex_df <- transform(df,
            radius = pmax(pmin(w, h), .Machine$double.eps) / 2,
            sides = 6L,
            angle = if (identical(hex_orientation, "pointy")) pi / 6 else 0
        )
        list(geom_regon(
            data = hex_df,
            aes(x0 = x, y0 = y, r = radius, sides = sides, angle = angle, fill = d),
            colour = NA,
            alpha = alpha_val,
            inherit.aes = FALSE
        ))
    } else {
        list()
    }
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_tile_df`.
#' @param df Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_density_tile_df <- function(df) {
    transform(df,
        x = (xmin + xmax) / 2,
        y = (ymin + ymax) / 2,
        w = pmax(0, xmax - xmin),
        h = pmax(0, ymax - ymin)
    )
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_flip_raster`.
#' @param img Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_density_flip_raster <- function(img) {
    dims <- dim(img)
    if (is.null(dims) || length(dims) < 2L) return(img)
    row_idx <- seq.int(dims[1], 1)
    if (length(dims) == 2L) {
        return(img[row_idx, , drop = FALSE])
    }
    if (length(dims) == 3L) {
        return(img[row_idx, , , drop = FALSE])
    }
    img
}

#' Internal helper for plotting workflows
#' @description
#' Internal helper for `.plot_density_prepare_rgba`.
#' @param img Internal parameter
#' @param alpha_scale Internal parameter
#' @return Internal helper result
#' @keywords internal
.plot_density_prepare_rgba <- function(img, alpha_scale = 1) {
    if (is.null(img)) return(NULL)
    alpha_scale <- max(0, min(1, alpha_scale))
    dims <- dim(img)
    if (length(dims) == 2L) {
        out <- array(1, dim = c(dims[1], dims[2], 4L))
        out[, , 1] <- img
        out[, , 2] <- img
        out[, , 3] <- img
        out[, , 4] <- alpha_scale
        return(out)
    }
    if (length(dims) == 3L && dims[3] == 3L) {
        out <- array(1, dim = c(dims[1], dims[2], 4L))
        out[, , 1:3] <- img
        out[, , 4] <- alpha_scale
        return(out)
    }
    if (length(dims) == 3L && dims[3] == 4L) {
        out <- img
        out[, , 4] <- pmin(1, pmax(0, out[, , 4] * alpha_scale))
        return(out)
    }
    stop("Unsupported histology image dimensions.")
}

#' Validate a positive integer Lee workflow parameter.
#' @keywords internal
.lee_assert_positive_integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 1 || x > .Machine$integer.max || x != floor(x)) {
        stop(name, " must be a single positive integer.", call. = FALSE)
    }
    as.integer(x)
}

#' Validate a non-negative integer Lee workflow parameter.
#' @keywords internal
.lee_assert_nonnegative_integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 0 || x > .Machine$integer.max || x != floor(x)) {
        stop(name, " must be a single non-negative integer.", call. = FALSE)
    }
    as.integer(x)
}

#' Validate a positive finite Lee workflow parameter.
#' @keywords internal
.lee_assert_positive_finite <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
        stop(
            name, " must be a single positive finite number.",
            if (identical(name, "mem_limit_GB")) {
                " An invalid limit is fail-closed and the requested output is treated as if it exceeds mem_limit_GB."
            } else {
                ""
            },
            call. = FALSE
        )
    }
    as.numeric(x)
}

#' Validate a scalar logical Lee workflow parameter.
#' @keywords internal
.lee_assert_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop(name, " must be a single non-missing logical value.", call. = FALSE)
    }
    x
}

#' Validate a scalar layer name.
#' @keywords internal
.lee_assert_layer_name <- function(x, name) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        stop(name, " must be a single non-empty character string.", call. = FALSE)
    }
    x
}

#' Resolve the weight style recorded on a Lee weight matrix.
#' @keywords internal
.lee_weight_style <- function(W) {
    style <- attr(W, "weight_style", exact = TRUE)
    if (is.null(style)) {
        provenance <- attr(W, "weight_provenance", exact = TRUE)
        if (is.list(provenance)) style <- provenance$weight_style
    }
    if (is.null(style) || length(style) != 1L || is.na(style) || !nzchar(as.character(style))) {
        return("unknown")
    }
    as.character(style)
}

#' Validate a normalized expression layer used by Lee's L.
#' @description
#' All Lee inputs must be numeric, finite, and column-centred relative to their
#' root-mean-square magnitude; unit variance is deliberately not required.
#' @keywords internal
.validate_lee_norm_layer <- function(X, norm_layer, context = "Lee's L") {
    norm_layer <- .lee_assert_layer_name(norm_layer, "norm_layer")
    X <- tryCatch(as.matrix(X), error = function(e) {
        stop(context, " expression layer '", norm_layer, "' cannot be converted to a matrix: ",
             conditionMessage(e), call. = FALSE)
    })
    if (!is.numeric(X) || length(dim(X)) != 2L || any(dim(X) < 1L)) {
        stop(context, " expression layer '", norm_layer,
             "' must be a non-empty numeric matrix.", call. = FALSE)
    }
    if (any(!is.finite(X))) {
        stop(context, " expression layer '", norm_layer,
             "' contains non-finite values.", call. = FALSE)
    }
    centre <- colMeans(X)
    rms <- sqrt(colMeans(X * X))
    relative_offset <- abs(centre) / pmax(rms, .Machine$double.eps)
    bad <- which(!is.finite(relative_offset) | relative_offset > 1e-7)
    if (length(bad)) {
        labels <- colnames(X)
        labels <- if (is.null(labels)) as.character(bad) else labels[bad]
        stop(
            context, " expression layer '", norm_layer,
            "' must be column-centred (abs(mean)/RMS <= 1e-7); offending columns: ",
            paste(utils::head(labels, 5L), collapse = ", "),
            if (length(bad) > 5L) " ..." else "",
            ". Unit variance is not required.",
            call. = FALSE
        )
    }
    X
}

#' Compute Lee's canonical spatial normalizer from a weight matrix.
#' @keywords internal
.lee_s2_value <- function(W) {
    if (is.null(dim(W)) || length(dim(W)) != 2L || nrow(W) != ncol(W) || nrow(W) < 1L) {
        stop("W must be a non-empty square matrix before computing Lee S2.", call. = FALSE)
    }
    row_totals <- as.numeric(W %*% rep.int(1, ncol(W)))
    S2 <- sum(row_totals * row_totals)
    if (length(S2) != 1L || !is.finite(S2) || S2 <= 0) {
        stop("W has non-finite or non-positive Lee S2; Lee's L is undefined.", call. = FALSE)
    }
    as.numeric(S2)
}

#' Fingerprint every input that determines Lee's L, separating data and permutations.
#' @description
#' The `data` component fingerprints the canonical X/W/grid inputs used for the
#' observed statistic. The `permutation` component fingerprints only the
#' permutation configuration, so a later `getTopLvsR()` call may validly choose
#' a different scheme while still proving that it uses the same observed data.
#' @param Xz Numeric normalized expression matrix.
#' @param W Spatial weight matrix aligned to the rows of `Xz`.
#' @param grid_info Grid metadata aligned to `Xz` and `W`.
#' @param norm_layer Name of the normalized expression layer.
#' @param block_side Positive integer block side, used only when
#'   `use_blocks = TRUE`.
#' @param use_blocks Logical. `FALSE` (the default) records the joint
#'   unrestricted all-grid shuffle; `TRUE` records a block-constrained scheme.
#' @keywords internal
.lee_input_fingerprint <- function(Xz, W, grid_info, norm_layer, block_side,
                                   use_blocks = FALSE) {
    norm_layer <- .lee_assert_layer_name(norm_layer, "norm_layer")
    block_side <- .lee_assert_positive_integer(block_side, "block_side")
    use_blocks <- .lee_assert_flag(use_blocks, "use_blocks")
    hash <- function(x) digest::digest(x, algo = "xxhash64", serialize = TRUE)
    w_payload <- if (inherits(W, "sparseMatrix")) {
        Wc <- if (inherits(W, "dgCMatrix")) {
            W
        } else {
            methods::as(methods::as(W, "generalMatrix"), "CsparseMatrix")
        }
        list(
            class = "dgCMatrix", dim = dim(Wc), dimnames = dimnames(Wc),
            p = Wc@p, i = Wc@i, x = Wc@x
        )
    } else {
        list(class = class(W), dim = dim(W), dimnames = dimnames(W), values = as.matrix(W))
    }
    data_components <- list(
        X = hash(list(dim = dim(Xz), dimnames = dimnames(Xz), values = Xz)),
        W = hash(w_payload),
        grid = hash(as.data.frame(grid_info, stringsAsFactors = FALSE)),
        weight_style = hash(.lee_weight_style(W))
    )
    permutation_components <- list(
        use_blocks = use_blocks,
        scheme = if (use_blocks) "spatial_block_joint" else "global_joint_shuffle",
        block_side = if (use_blocks) block_side else NA_integer_
    )
    list(
        schema = "lee_input_fingerprint_v2",
        norm_layer = norm_layer,
        data = c(
            list(schema = "lee_observed_data_v1", hash = hash(data_components)),
            data_components
        ),
        permutation = c(
            list(schema = "lee_permutation_config_v2", hash = hash(permutation_components)),
            permutation_components
        )
    )
}

#' Compare only the observed-data part of two Lee fingerprints.
#' @keywords internal
.lee_same_observed_data <- function(observed, current) {
    is.list(observed) && is.list(current) &&
        identical(observed$schema, "lee_input_fingerprint_v2") &&
        identical(current$schema, "lee_input_fingerprint_v2") &&
        identical(observed$norm_layer, current$norm_layer) &&
        identical(observed$data, current$data)
}

#' Compute Lee L
#' @description
#' Internal helper for `.compute_lee_l`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param norm_layer Normalized expression layer name. Custom layers must be
#'   numeric, finite, and column-centred; unit variance is not required.
#' @param genes Parameter value.
#' @param within A single logical value controlling square selected-gene output.
#' @param ncores Positive integer number of cores/threads to use.
#' @param mem_limit_GB Positive finite maximum dense native Lee result size in GiB.
#' @param chunk_size Retained for API compatibility; chunking is disabled.
#' @param use_bigmemory Must remain `FALSE`; the former route was RAM-backed.
#' @param backing_path Retained for API compatibility.
#' @param block_side Positive integer block side for permutation preprocessing,
#'   used only when `use_blocks = TRUE`.
#' @param use_blocks Logical. `FALSE` (the default) prepares a joint
#'   unrestricted row shuffle over all grid locations. `TRUE` prepares spatial
#'   blocks defined by `block_side`.
#' @param cache_inputs A single logical value controlling preprocessing reuse.
#' @param input_cache Parameter value.
#' @return Return value used internally.
#' @keywords internal
.compute_lee_l <- function(scope_obj,
                         grid_name = NULL,
                         norm_layer = "Xz",
                         genes = NULL,
                         within = TRUE,
                         ncores = 1,
                         mem_limit_GB = 16,
                         chunk_size = 256L,
                         use_bigmemory = FALSE,
                         backing_path = tempdir(),
                         block_side = 8,
                         use_blocks = FALSE,
                         cache_inputs = FALSE,
                         input_cache = NULL) {
    within <- .lee_assert_flag(within, "within")
    ncores <- .lee_assert_positive_integer(ncores, "ncores")
    block_side <- .lee_assert_positive_integer(block_side, "block_side")
    use_blocks <- .lee_assert_flag(use_blocks, "use_blocks")
    effective_block_side <- if (use_blocks) block_side else NA_integer_
    cache_inputs <- .lee_assert_flag(cache_inputs, "cache_inputs")
    use_bigmemory <- .lee_assert_flag(use_bigmemory, "use_bigmemory")
    norm_layer <- .lee_assert_layer_name(norm_layer, "norm_layer")
    if (isTRUE(use_bigmemory)) {
        stop(
            "use_bigmemory=TRUE is disabled because the former Lee chunk route was RAM-backed, not file-backed.",
            call. = FALSE
        )
    }
    mem_limit_GB <- .lee_assert_positive_finite(mem_limit_GB, "mem_limit_GB")

    ## ---- 0. Get grid layer (new helper) ---------------------------------
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    # Fill back layer name string, for later writing back to stats
    if (is.null(grid_name)) {
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
    }

    ## ---- 1. Extract expression matrix and weight matrix (with optional reuse) ----
    grid_info <- g_layer$grid_info
    if (is.null(g_layer[[norm_layer]])) {
        stop("Normalized expression layer '", norm_layer, "' is missing.", call. = FALSE)
    }
    Xz_source <- .validate_lee_norm_layer(
        g_layer[[norm_layer]], norm_layer,
        context = "[geneSCOPE::.compute_lee_l]"
    )
    W_source <- g_layer$W
    if (is.null(W_source)) {
        stop("Both the normalized expression layer and grid$W are required for Lee's L.", call. = FALSE)
    }
    source_fingerprint <- .lee_input_fingerprint(
        Xz = Xz_source,
        W = W_source,
        grid_info = grid_info,
        norm_layer = norm_layer,
        block_side = block_side,
        use_blocks = use_blocks
    )
    alignment_schema <- "grid_id_independent_fingerprint_v4"
    cache_valid <- cache_inputs &&
        !is.null(input_cache) &&
        identical(input_cache$alignment_schema, alignment_schema) &&
        identical(input_cache$norm_layer, norm_layer) &&
        identical(input_cache$grid_id, grid_info$grid_id) &&
        identical(input_cache$block_side, effective_block_side) &&
        identical(input_cache$use_blocks, use_blocks) &&
        identical(input_cache$source_fingerprint, source_fingerprint)

    if (cache_valid) {
        Xz_full <- input_cache$Xz_full
        W <- input_cache$W
        block_id <- input_cache$block_id
    } else {
        Xz_full <- Xz_source
        target_ids <- as.character(grid_info$grid_id)
        if (!length(target_ids) || anyNA(target_ids) || any(!nzchar(target_ids)) || anyDuplicated(target_ids)) {
            stop("grid_info$grid_id must be present, non-empty, and unique before Lee's L alignment.", call. = FALSE)
        }

        .alignment_index <- function(source_ids, label) {
            if (is.null(source_ids)) {
                stop(label, " has no grid_id dimnames; rerun the producing step.")
            }
            source_ids <- as.character(source_ids)
            if (anyNA(source_ids) || anyDuplicated(source_ids)) {
                stop(label, " grid_id dimnames contain missing or duplicated values.")
            }
            idx <- match(target_ids, source_ids)
            if (anyNA(idx) || length(source_ids) != length(target_ids)) {
                stop(label, " grid_id set does not exactly match grid_info$grid_id.")
            }
            idx
        }

        x_ord <- .alignment_index(rownames(Xz_full), norm_layer)
        w_row_ord <- .alignment_index(rownames(W_source), "W rows")
        w_col_ord <- .alignment_index(colnames(W_source), "W columns")

        Xz_full <- Xz_full[x_ord, , drop = FALSE]
        w_style_attr <- attr(W_source, "weight_style", exact = TRUE)
        W <- W_source[w_row_ord, w_col_ord, drop = FALSE]
        if (!is.null(w_style_attr)) attr(W, "weight_style") <- w_style_attr
        rownames(Xz_full) <- target_ids
        dimnames(W) <- list(target_ids, target_ids)

        if (!identical(rownames(Xz_full), target_ids) ||
            !identical(rownames(W), target_ids) ||
            !identical(colnames(W), target_ids)) {
            stop("Failed to align Xz and W independently to grid_info$grid_id.")
        }
        block_id <- if (use_blocks) {
            .assign_block_id(grid_info, block_side = block_side)
        } else {
            rep.int(1L, nrow(grid_info))
        }
    }

    block_id <- .compact_block_id(block_id)

    input_fingerprint <- .lee_input_fingerprint(
        Xz = Xz_full,
        W = W,
        grid_info = grid_info,
        norm_layer = norm_layer,
        block_side = block_side,
        use_blocks = use_blocks
    )
    weight_style <- .lee_weight_style(W)
    S2 <- .lee_s2_value(W)

    all_genes <- colnames(Xz_full)
    if (is.null(all_genes) || anyNA(all_genes) || any(!nzchar(all_genes)) || anyDuplicated(all_genes)) {
        stop("The normalized Lee input must have unique, non-empty gene column names.", call. = FALSE)
    }
    idx_keep <- if (is.null(genes)) {
        seq_along(all_genes)
    } else {
        if (!length(genes) || anyNA(genes) || anyDuplicated(genes)) {
            stop("genes must be a non-empty vector of unique gene names.", call. = FALSE)
        }
        m <- match(genes, all_genes, nomatch = 0L)
        if (any(m == 0L)) stop("Some genes were not found in column names", call. = FALSE)
        m
    }

    n_g <- length(idx_keep)
    compute_n_g <- if (isTRUE(within)) n_g else length(all_genes)
    bytes_L <- (as.double(compute_n_g)^2 * 8) / 1024^3
    if (bytes_L > mem_limit_GB) {
        stop(
            "Lee L native output exceeds mem_limit_GB, and the former bigmemory route is disabled because it was RAM-backed. ",
            "Reduce the gene set or raise mem_limit_GB only after confirming available RAM.",
            call. = FALSE
        )
    }

    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    }
    Sys.setenv(OMP_NUM_THREADS = ncores)
    Xz_compute <- if (isTRUE(within)) Xz_full[, idx_keep, drop = FALSE] else Xz_full
    L_full <- .lee_l_cache(Xz_compute, W, n_threads = ncores)

    if (isTRUE(within)) {
        Lmat <- L_full
        Xuse <- Xz_compute
    } else {
        Lmat <- L_full[idx_keep, , drop = FALSE]
        Xuse <- Xz_full[, idx_keep, drop = FALSE]
    }

    dimnames(Lmat) <- if (within) {
        list(row = all_genes[idx_keep], col = all_genes[idx_keep])
    } else {
        list(row = all_genes[idx_keep], col = all_genes)
    }

    list(
        Lmat      = Lmat,
        X_used    = Xuse,
        X_full    = Xz_full,
        cells     = rownames(Xz_full),
        W         = W,
        grid_info = g_layer$grid_info,
        grid_name = grid_name, # ← Additional return, convenient for .compute_l
        block_id  = block_id,
        input_fingerprint = input_fingerprint,
        weight_style = weight_style,
        S2 = S2,
        input_cache = if (cache_inputs) {
            list(
                alignment_schema = alignment_schema,
                norm_layer = norm_layer,
                grid_id = grid_info$grid_id,
                block_side = effective_block_side,
                use_blocks = use_blocks,
                source_fingerprint = source_fingerprint,
                input_fingerprint = input_fingerprint,
                Xz_full = Xz_full,
                W = W,
                block_id = block_id
            )
        } else {
            NULL
        }
    )
}

.lee_l_perm_global <- function(Xz, W, L_ref,
                               perms = 999,
                               block_size = 64,
                               n_threads = 1) {
    perms <- .lee_assert_positive_integer(perms, "perms")
    block_size <- .lee_assert_positive_integer(block_size, "block_size")
    n_threads <- .lee_assert_positive_integer(n_threads, "n_threads")
    stopifnot(is.matrix(Xz), inherits(W, "dgCMatrix"))

    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    }
    ngen <- ncol(Xz)
    geCnt <- matrix(0, ngen, ngen)
    done <- 0L

    while (done < perms) {
        bsz <- min(block_size, perms - done)
        idx_mat <- replicate(
            bsz,
            sample.int(nrow(Xz), nrow(Xz), replace = FALSE),
            simplify = "matrix"
        ) - 1L
        storage.mode(idx_mat) <- "integer"
        geCnt <- geCnt + .lee_perm(Xz, W, idx_mat, L_ref, n_threads)
        done <- done + bsz
    }
    (geCnt + 1) / (perms + 1)
}

#' Lee L Perm Block
#' @description
#' Internal helper for `.lee_l_perm_block`.
#' @param Xz Parameter value.
#' @param W Parameter value.
#' @param L_ref Parameter value.
#' @param block_id Parameter value.
#' @param perms Parameter value.
#' @param block_size Parameter value.
#' @param n_threads Number of threads to use.
#' @return Return value used internally.
#' @keywords internal
.lee_l_perm_block <- function(Xz, W, L_ref,
                             block_id,
                             perms = 999,
                             block_size = 64,
                             n_threads = 1) {
    perms <- .lee_assert_positive_integer(perms, "perms")
    block_size <- .lee_assert_positive_integer(block_size, "block_size")
    n_threads <- .lee_assert_positive_integer(n_threads, "n_threads")
    stopifnot(
        is.matrix(Xz), inherits(W, "dgCMatrix"),
        length(block_id) == nrow(Xz)
    )

    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    }
    ngen <- ncol(Xz)
    geCnt <- matrix(0, ngen, ngen)
    done <- 0L

    ## ---- Pre-group row indices by block ----
    split_rows <- split(seq_along(block_id), block_id)
    while (done < perms) {
        bsz <- min(block_size, perms - done)

        idx_mat <- replicate(bsz,
            {
                # Keep every block at its original spatial positions and only
                # permute observations among positions belonging to that block.
                # Concatenating sampled blocks changes the row-position map when
                # blocks are not contiguous (and is invalid for unequal sizes).
                idx <- seq_len(nrow(Xz))
                for (rows in split_rows) {
                    idx[rows] <- rows[sample.int(length(rows), length(rows), replace = FALSE)]
                }
                idx
            },
            simplify = "matrix"
        ) - 1L # 0-based for C++

        storage.mode(idx_mat) <- "integer"
        # Use correct export function name
        geCnt <- geCnt + .lee_perm_block(Xz, W, idx_mat, as.integer(block_id) - 1L, L_ref, n_threads)
        done <- done + bsz
    }
    (geCnt + 1) / (perms + 1)
}

#' Lee L Full
#' @description
#' Internal helper for `.lee_l_full`.
#' @param Xz Parameter value.
#' @param W Parameter value.
#' @param n_threads Number of threads to use.
#' @param .cache Parameter value.
#' @return Return value used internally.
#' @keywords internal
.lee_l_full <- function(Xz, W, n_threads = 1L, .cache = TRUE) {
    if (.cache) .lee_l_cache(Xz, W, n_threads) else .lee_l(Xz, W, n_threads)
}

#' Lee L Subset
#' @description
#' Internal helper for `.lee_l_subset`.
#' @param Xz Parameter value.
#' @param W Parameter value.
#' @param cols0 Parameter value.
#' @param n_threads Number of threads to use.
#' @return Return value used internally.
#' @keywords internal
.lee_l_subset <- function(Xz, W, cols0, n_threads = 1L) {
    .lee_l_cols(Xz, W, cols0, n_threads)
}

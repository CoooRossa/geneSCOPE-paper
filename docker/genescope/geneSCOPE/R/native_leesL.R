#' Compute Lee L
#' @description
#' Internal helper for `.compute_lee_l`.
#' @param scope_obj A `scope_object`.
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
                         use_bigmemory = TRUE,
                         backing_path = tempdir(),
                         block_side = 8,
                         cache_inputs = FALSE,
                         input_cache = NULL) {
    if (use_bigmemory && !requireNamespace("bigmemory", quietly = TRUE)) {
        use_bigmemory <- FALSE
    }

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
    cache_valid <- cache_inputs &&
        !is.null(input_cache) &&
        identical(input_cache$norm_layer, norm_layer) &&
        identical(input_cache$grid_id, grid_info$grid_id) &&
        identical(input_cache$block_side, block_side)

    if (cache_valid) {
        Xz_full <- input_cache$Xz_full
        W <- input_cache$W
        block_id <- input_cache$block_id
    } else {
        Xz_full <- as.matrix(g_layer[[norm_layer]])
        ord <- match(grid_info$grid_id, rownames(Xz_full))
        Xz_full <- Xz_full[ord, , drop = FALSE]
        W <- g_layer$W[ord, ord]
        block_id <- .assign_block_id(grid_info, block_side = block_side)
    }

    block_id <- .compact_block_id(block_id)

    all_genes <- colnames(Xz_full)
    idx_keep <- if (is.null(genes)) {
        seq_along(all_genes)
    } else {
        m <- match(genes, all_genes, nomatch = 0L)
        if (any(m == 0L)) stop("Some genes were not found in column names")
        m
    }

    n_g <- length(idx_keep)
    bytes_L <- n_g^2 * 8 / 1024^3 # in GB
    need_stream <- use_bigmemory && bytes_L > mem_limit_GB

    ## ======================================================
    ## =============  A. One‑shot computation (fits in RAM) ===============
    ## ======================================================
    if (!need_stream) {
        RhpcBLASctl::blas_set_num_threads(1)
        Sys.setenv(OMP_NUM_THREADS = ncores)
        # Use correct export function name
        L_full <- .lee_l_cache(Xz_full, W, n_threads = ncores)

        if (within) {
            Lmat <- L_full[idx_keep, idx_keep, drop = FALSE]
            Xuse <- Xz_full[, idx_keep, drop = FALSE]
        } else {
            Lmat <- L_full[idx_keep, , drop = FALSE]
            Xuse <- Xz_full[, idx_keep, drop = FALSE]
        }

        ## ======================================================
        ## =============  B. Chunked computation with file mapping ================
        ## ======================================================
    } else {
        ## Use in-memory shared big.matrix to avoid unsupported 'shared' arg on filebacked in some versions
        ## Note: this stores the chunked result in RAM. Ensure sufficient memory is available.
        L_bm <- bigmemory::big.matrix(
            nrow = n_g, ncol = n_g, type = "double",
            init = NA_real_,
            dimnames = NULL,
            shared = TRUE
        )

        RhpcBLASctl::blas_set_num_threads(1)
        Sys.setenv(OMP_NUM_THREADS = ncores)

        ## ---- Process column blocks and write ----
        for (start in seq(1L, n_g, by = chunk_size)) {
            idx_chunk <- idx_keep[start:min(n_g, start + chunk_size - 1L)]
            # Use correct export function name
            L_block <- .lee_l_cols(
                Xz_full, W,
                cols0 = as.integer(idx_chunk - 1L),
                n_threads = ncores
            )

            ## Write column block + symmetric columns
            L_bm[, idx_chunk] <- L_block
            L_bm[idx_chunk, ] <- t(L_block)

            rm(L_block)
            gc(verbose = FALSE)
        }

        ## ---- Add gene names (only two vectors) ----
        dimnames(L_bm) <- list(
            all_genes[idx_keep],
            all_genes[idx_keep]
        )

        Lmat <- L_bm
        Xuse <- Xz_full[, idx_keep, drop = FALSE]
    }

    dimnames(Lmat) <- list(
        row = colnames(Xuse),
        col = colnames(Xuse)
    )

    list(
        Lmat      = Lmat,
        X_used    = Xuse,
        X_full    = Xz_full,
        cells     = rownames(Xz_full),
        W         = W,
        grid_info = g_layer$grid_info,
        grid_name = grid_name, # ← Additional return, convenient for .compute_l
        block_id  = block_id,
        input_cache = if (cache_inputs) {
            list(
                norm_layer = norm_layer,
                grid_id = grid_info$grid_id,
                block_side = block_side,
                Xz_full = Xz_full,
                W = W,
                block_id = block_id
            )
        } else {
            NULL
        }
    )
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
    stopifnot(
        is.matrix(Xz), inherits(W, "dgCMatrix"),
        length(block_id) == nrow(Xz)
    )

    RhpcBLASctl::blas_set_num_threads(1)
    ngen <- ncol(Xz)
    geCnt <- matrix(0, ngen, ngen)
    done <- 0L

    ## ---- Pre-group row indices by block ----
    split_rows <- split(seq_along(block_id), block_id)
    blk_keys <- names(split_rows)
    n_blk <- length(split_rows)

    while (done < perms) {
        bsz <- min(block_size, perms - done)

        idx_mat <- replicate(bsz,
            {
                new_order <- sample(n_blk) # shuffle block order
                unlist(split_rows[new_order], use.names = FALSE)
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

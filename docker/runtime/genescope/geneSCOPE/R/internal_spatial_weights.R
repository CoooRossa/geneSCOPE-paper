
#' Spw Trace Enabled
#' @description
#' Internal helper for `.spw_trace_enabled`.
#' @return Return value used internally.
#' @keywords internal
.spw_trace_enabled <- local({
  .cache <- NULL
  function() {
    if (is.null(.cache)) {
      flag <- Sys.getenv("GENESCOPE_SPW_TRACE", "")
      .cache <<- tolower(flag) %in% c("1", "true", "yes", "on")
    }
    .cache
  }
})

#' Spw Trace
#' @description
#' Internal helper for `.spw_trace`.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.spw_trace <- function(..., parent = "computeWeights", step = "S03", verbose = TRUE) {
  if (!.spw_trace_enabled() || !.log_enabled(verbose)) return(invisible(NULL))
  msg <- paste0(..., collapse = "")
  .log_info(parent, step, msg, verbose)
  invisible(NULL)
}

# Internal aliases for native helpers (RcppExports uses non-dot names).
.grid_nb_omp <- function(nrow, ncol, queen = TRUE) grid_nb_omp(nrow, ncol, queen)
.grid_nb <- function(nrow, ncol, queen = TRUE) grid_nb(nrow, ncol, queen)
.grid_nb_hex_omp <- function(nrow, ncol, oddr = TRUE) grid_nb_hex_omp(nrow, ncol, oddr)
.grid_nb_hex <- function(nrow, ncol, oddr = TRUE) grid_nb_hex(nrow, ncol, oddr)
.grid_nb_hexq_omp <- function(nrow, ncol, oddq = TRUE) grid_nb_hexq_omp(nrow, ncol, oddq)
.grid_nb_hexq <- function(nrow, ncol, oddq = TRUE) grid_nb_hexq(nrow, ncol, oddq)
.grid_weights_kernel_rect_omp <- function(nrow, ncol, gx, gy, queen = TRUE, radius = 2L, kernel = "gaussian", sigma = 1.0) {
  grid_weights_kernel_rect_omp(nrow, ncol, gx, gy, queen, radius, kernel, sigma)
}
.grid_weights_kernel_hexr_omp <- function(nrow, ncol, gx, gy, oddr = TRUE, radius = 2L, kernel = "gaussian", sigma = 1.0) {
  grid_weights_kernel_hexr_omp(nrow, ncol, gx, gy, oddr, radius, kernel, sigma)
}
.grid_weights_kernel_hexq_omp <- function(nrow, ncol, gx, gy, oddq = TRUE, radius = 2L, kernel = "gaussian", sigma = 1.0) {
  grid_weights_kernel_hexq_omp(nrow, ncol, gx, gy, oddq, radius, kernel, sigma)
}
.listw_b_omp <- function(nb) listw_B_omp(nb)

#' Builds spatial neighbourhood structures, converts them to listw/sparse
#' @description
#' Internal helper for `.compute_weights`.
#' matrices, and stores results under the grid layer (including fuzzy weights).
#' @param scope_obj A `scope_object` with grid layers.
#' @param grid_name Grid layer name (defaults to active layer).
#' @param style Style for `nb2listw` (default `B`).
#' @param store_mat Store the adjacency matrix in the grid layer when TRUE.
#' @param zero_policy Whether zero-neighbour cells are allowed.
#' @param store_listw Store the `listw` object when TRUE.
#' @param verbose Emit progress messages when TRUE.
#' @param topology Topology choice (`auto`, `queen`, `rook`, `hex`, `fuzzy_queen`, `fuzzy_hex`).
#' @param nb_order Optional higher-order neighbour expansion.
#' @param min_neighbors Minimum neighbour count when expanding (non-kernel styles).
#' @param max_order Maximum neighbour order (non-kernel styles).
#' @param kernel_radius Kernel radius for `style="kernel_*"` (overrides `max_order`).
#' @param kernel_sigma Kernel sigma for `style="kernel_*"` (overrides `min_neighbors`).
#' @param ncores Number of threads to use.
#' @return The modified `scope_object`.
#' @keywords internal
.compute_weights <- function(scope_obj,
                           grid_name = NULL,
                           style = "B",
                           store_mat = TRUE,
                           zero_policy = TRUE,
                           store_listw = TRUE,
                           verbose = TRUE,
                           topology = c("auto", "queen", "rook", "hex", "fuzzy_queen", "fuzzy_hex"),
                           nb_order = 1L,
                           min_neighbors = NULL,
                           max_order = NULL,
                           kernel_radius = NULL,
                           kernel_sigma = NULL,
                           ncores = detectCores(logical = TRUE)) {
  topology <- match.arg(topology)
  use_fuzzy <- grepl("^fuzzy", topology)
  topo_arg <- if (use_fuzzy) {
    if (identical(topology, "fuzzy_hex")) "hex" else "queen"
  } else {
    topology
  }
  parent <- "computeWeights"
  step01 <- .log_step(parent, "S01", "select grid layer and validate", verbose)
  grid_label <- if (is.null(grid_name)) "auto" else as.character(grid_name)[1]
  step01$enter(paste0("grid_name=", grid_label))

  guard <- .spw_thread_guard_v2()
  on.exit(guard$restore(), add = TRUE)

  # 1) Select grid layer and basic metadata
  layer_info <- .spw_select_layer(scope_obj, grid_name)
  grid_name <- layer_info$grid_name
  layer <- layer_info$layer
  meta <- .spw_validate_grid_metadata(layer, grid_name)
  ncores_use <- max(1L, min(ncores, detectCores(logical = TRUE)))
  .log_info(parent, "S01", paste0(
    "grid_name=", grid_name,
    " cells=", length(meta$grid_id),
    " dims=", meta$ybins_eff, "x", meta$xbins_eff,
    " ncores_use=", ncores_use
  ), verbose)
  step01$done(paste0("grid_name=", grid_name, " cells=", length(meta$grid_id)))

  # 2) Determine lattice topology (auto or forced)
  step02 <- .log_step(parent, "S02", "resolve topology and kernel params", verbose)
  step02$enter(paste0("topology=", topology, " style=", style))
  topo_info <- .spw_choose_topology(scope_obj, layer, topo_arg, verbose,
    parent = parent, step = "S02"
  )

  # Kernel-smoothed spatial weights (main-2 parity): enable via style = "kernel_*".
  kernel_style <- is.character(style) && length(style) == 1L && nzchar(style) && grepl("^kernel", style)
  kernel_token <- NULL
  kernel_radius_use <- NULL
  kernel_sigma_use <- NULL
  store_mat_effective <- store_mat
  if (kernel_style) {
    if (use_fuzzy) {
      stop("Kernel weights are incompatible with fuzzy_* topology.")
    }

    kernel_token <- tolower(sub("^kernel[_-]?", "", style))
    if (!nzchar(kernel_token)) kernel_token <- "gaussian"
    if (!kernel_token %in% c("gaussian", "flat")) {
      stop("Unknown kernel style: ", style, " (supported: kernel_gaussian, kernel_flat)")
    }

    # Defaults mirror main-2 (radius=2, sigma=1.0) without changing exported formals.
    kernel_radius_use <- if (is.null(kernel_radius)) {
      if (!is.null(max_order)) {
        suppressWarnings(as.integer(max_order))
      } else {
        ord <- suppressWarnings(as.integer(nb_order))
        if (length(ord) && is.finite(ord) && ord >= 2L) ord else 2L
      }
    } else {
      suppressWarnings(as.integer(kernel_radius))
    }
    kernel_sigma_use <- if (is.null(kernel_sigma)) {
      if (is.null(min_neighbors)) 1.0 else suppressWarnings(as.numeric(min_neighbors))
    } else {
      suppressWarnings(as.numeric(kernel_sigma))
    }

    if (!is.finite(kernel_radius_use) || kernel_radius_use < 1L) {
      stop("kernel_radius must be >= 1 (set `max_order` when using style='kernel_*').")
    }
    if (!is.finite(kernel_sigma_use) || kernel_sigma_use <= 0) {
      stop("kernel_sigma must be > 0 (set `min_neighbors` when using style='kernel_*').")
    }

    if (!store_mat && verbose) {
      .log_info(parent, "S02", paste0("style='", style, "' forces store_mat=TRUE; storing W."), verbose)
    }
    store_mat_effective <- TRUE
    .log_info(parent, "S02", paste0(
      "kernel=", kernel_token,
      " radius=", kernel_radius_use,
      " sigma=", kernel_sigma_use
    ), verbose)
  }
  step02$done(paste0(
    "topology=", topo_info$topology,
    " kernel=", if (kernel_style) kernel_token else "none"
  ))

  step03 <- .log_step(parent, "S03", "build spatial weights", verbose)
  step03$enter(paste0("path=", if (kernel_style) "kernel" else "adjacency"))
  .log_backend(parent, "S03", "weights_path",
    if (kernel_style) "kernel" else "adjacency",
    reason = paste0("style=", style),
    verbose = verbose
  )

  openmp_info <- tryCatch(.native_openmp_info(), error = function(e) NULL)
  compiled_omp <- FALSE
  omp_threads <- NA_integer_
  if (is.list(openmp_info) && !is.null(openmp_info$compiled_with_openmp)) {
    compiled_omp <- isTRUE(openmp_info$compiled_with_openmp)
    omp_threads <- openmp_info$omp_max_threads
  }
  if (isTRUE(compiled_omp)) {
    .log_backend(parent, "S03", "openmp",
      paste0("compiled_with_openmp=TRUE omp_max_threads=", omp_threads),
      verbose = verbose
    )
  } else {
    .log_backend(parent, "S03", "openmp",
      "compiled_with_openmp=FALSE",
      reason = "fallback",
      verbose = verbose
    )
  }

  W <- NULL
  listw_obj <- NULL
  relabel_nb <- NULL
  nb_backend <- NULL
  zero_idx <- integer(0)

  if (kernel_style) {
    topology_resolved <- topo_info$topology
    .log_backend(parent, "S03", "kernel_builder", "C++", verbose = verbose)
    .log_info(parent, "S03", paste0(
      "topology=", topology_resolved,
      " kernel=", kernel_token,
      " radius=", kernel_radius_use,
      " sigma=", kernel_sigma_use
    ), verbose)
    W <- if (topology_resolved %in% c("queen", "rook")) {
      .grid_weights_kernel_rect_omp(
        nrow = meta$ybins_eff,
        ncol = meta$xbins_eff,
        gx = meta$gx,
        gy = meta$gy,
        queen = identical(topology_resolved, "queen"),
        radius = kernel_radius_use,
        kernel = kernel_token,
        sigma = kernel_sigma_use
      )
    } else if (topology_resolved %in% c("hex-oddr", "hex-evenr")) {
      .grid_weights_kernel_hexr_omp(
        nrow = meta$ybins_eff,
        ncol = meta$xbins_eff,
        gx = meta$gx,
        gy = meta$gy,
        oddr = identical(topology_resolved, "hex-oddr"),
        radius = kernel_radius_use,
        kernel = kernel_token,
        sigma = kernel_sigma_use
      )
    } else if (topology_resolved %in% c("hex-oddq", "hex-evenq")) {
      .grid_weights_kernel_hexq_omp(
        nrow = meta$ybins_eff,
        ncol = meta$xbins_eff,
        gx = meta$gx,
        gy = meta$gy,
        oddq = identical(topology_resolved, "hex-oddq"),
        radius = kernel_radius_use,
        kernel = kernel_token,
        sigma = kernel_sigma_use
      )
    } else {
      stop("Unknown topology for kernel weights: ", topology_resolved)
    }
    diag(W) <- 0
    dimnames(W) <- list(meta$grid_id, meta$grid_id)
  } else {
    # 3) Build the full neighbourhood structure and relabel to active cells
    nb_bundle <- .spw_build_neighbourhood(meta, topo_info,
      ncores = ncores_use, verbose = verbose, parent = parent, step = "S03"
    )
    relabel_nb <- nb_bundle$nb
    nb_backend <- nb_bundle$backend
    nb_label <- if (identical(nb_backend, "R")) {
      "R .grid_nb_r"
    } else if (identical(nb_backend, "omp")) {
      "C++ .grid_nb_omp"
    } else if (identical(nb_backend, "serial") &&
      topo_info$topology %in% c("hex-oddr", "hex-evenr", "hex-oddq", "hex-evenq")) {
      "grid_nb_hex (serial)"
    } else if (identical(nb_backend, "serial")) {
      "C++ .grid_nb"
    } else {
      nb_backend
    }
    .log_backend(parent, "S03", "nb_builder", nb_label, verbose = verbose)

  # 3.1) Optionally expand to higher-order neighbourhoods for broader context
  nb_expand <- .spw_expand_neighbourhood(
    relabel_nb,
    nb_order = nb_order,
    min_neighbors = min_neighbors,
    max_order = max_order,
    verbose = verbose,
    parent = parent,
    step = "S03"
  )
  relabel_nb <- nb_expand$nb
  attr(relabel_nb, "avg_neighbors") <- nb_expand$avg_neighbors

  # 3.2) Handle zero-neighbour regions after any expansion
  zero_idx <- integer(0)
  if (zero_policy) {
    zero_idx <- which(vapply(relabel_nb, length, integer(1)) == 0L)
    if (length(zero_idx) > 0L) {
      for (j in zero_idx) {
        relabel_nb[[j]] <- j
      }
    }
  }

  # 4) Convert to requested outputs
    if (store_mat_effective) {
      W <- .listw_b_omp(relabel_nb)
      diag(W) <- 0
      W@x[] <- 1
      dimnames(W) <- list(meta$grid_id, meta$grid_id)
      .log_backend(parent, "S03", "matrix_builder", "C++ .listw_b_omp", verbose = verbose)
    }
  }

  w_nnz_step3 <- if (!is.null(W) && inherits(W, "Matrix")) {
    length(W@x)
  } else if (!is.null(W) && is.matrix(W)) {
    sum(W != 0)
  } else {
    NA_integer_
  }
  step03$done(paste0(
    "path=", if (kernel_style) "kernel" else "adjacency",
    if (isTRUE(is.finite(w_nnz_step3))) paste0(" nnz=", w_nnz_step3) else ""
  ))

  step04 <- .log_step(parent, "S04", "convert/store listw and fuzzy outputs", verbose)
  step04$enter(paste0(
    "store_mat=", store_mat_effective,
    " store_listw=", store_listw,
    if (use_fuzzy) paste0(" fuzzy=", topo_arg) else ""
  ))

  if (kernel_style) {
    if (store_mat_effective) {
      scope_obj@grid[[grid_name]]$W <- W
    }
    if (store_listw) {
      scope_obj@grid[[grid_name]]$listw <- spdep::mat2listw(W,
        style = "W",
        zero.policy = zero_policy
      )
      .log_backend(parent, "S04", "listw_builder", "spdep::mat2listw", verbose = verbose)
    }
  } else {
    if (store_listw || store_mat_effective) {
      listw_obj <- spdep::nb2listw(relabel_nb,
        style = style,
        zero.policy = zero_policy
      )
    }
    if (store_mat_effective) {
      scope_obj@grid[[grid_name]]$W <- W
    }
    if (store_listw) {
      scope_obj@grid[[grid_name]]$listw <- listw_obj
      .log_backend(parent, "S04", "listw_builder", "spdep::nb2listw", verbose = verbose)
    }
  }

  # Optionally restore empty neighbour vectors for downstream consumers
  if (zero_policy && length(zero_idx)) {
    for (j in zero_idx) relabel_nb[[j]] <- integer(0)
  }

  if (use_fuzzy) {
    if (!store_mat) {
      stop("Fuzzy ", topo_arg, " topology requires store_mat = TRUE to retain the adjacency matrix.")
    }
    W_base <- scope_obj@grid[[grid_name]]$W
    scope_obj@grid[[grid_name]]$fuzzy <- list(
      base_topology = topo_arg,
      W_binary = W_base,
      W_row_normalised = .row_normalize_sparse(W_base)
    )
    .log_info(parent, "S04", paste0(
      "stored fuzzy weights topology=", topo_arg,
      " row_normalised=TRUE"
    ), verbose)
  }
  step04$done("stored outputs")

  step05 <- .log_step(parent, "S05", "finalize and return", verbose)
  step05$enter("summary")
  w_nnz <- if (!is.null(W) && inherits(W, "Matrix")) {
    length(W@x)
  } else if (!is.null(W) && is.matrix(W)) {
    sum(W != 0)
  } else {
    NA_integer_
  }
  w_dim <- if (!is.null(W)) paste0(nrow(W), "x", ncol(W)) else "NA"
  diag_zero <- if (!is.null(W) && !is.matrix(W)) {
    all(diag(W) == 0)
  } else if (!is.null(W) && is.matrix(W)) {
    all(diag(W) == 0)
  } else {
    NA
  }
  w_binary <- if (!is.null(W) && inherits(W, "Matrix")) {
    if (!length(W@x)) TRUE else all(W@x %in% c(0, 1))
  } else if (!is.null(W) && is.matrix(W)) {
    if (!length(W)) TRUE else all(W[W != 0] %in% c(0, 1))
  } else {
    NA
  }
  row_norm_flag <- if (use_fuzzy) "fuzzy" else "FALSE"
  summary_msg <- paste0(
    "grid_name=", grid_name,
    " W=", if (!is.null(W)) "TRUE" else "FALSE",
    " dims=", w_dim,
    " nnz=", w_nnz,
    " diag_zero=", diag_zero,
    " binary=", w_binary,
    " row_normalized=", row_norm_flag,
    " listw=", store_listw,
    if (use_fuzzy) " fuzzy=TRUE" else ""
  )
  step05$done(summary_msg)
  invisible(scope_obj)
}

#' Validate grid dimensions for neighbourhood builders
#' @description
#' Internal helper for `.spw_validate_grid_shape`.
#' Ensures grid dimensions are scalar, integerish, and within safe bounds for
#' native neighbourhood builders.
#' @param xbins_eff ybins_eff Effective grid dimensions.
#' @param context Character string used for error messages.
#' @param ybins_eff Parameter value.
#' @return List with validated integer `xbins_eff`, `ybins_eff`, and `n_cells`.
#' @keywords internal
.spw_validate_grid_shape <- function(xbins_eff, ybins_eff, context = "[geneSCOPE::.compute_weights]") {
  validate_dim <- function(val, label) {
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
      stop(context, " Grid metadata field '", label, "' must be a finite numeric scalar.")
    }
    if (!.spw_is_integerish(val)) {
      stop(context, " Grid metadata field '", label, "' must be an integer.")
    }
    iv <- as.integer(round(val))
    if (iv <= 0L) {
      stop(context, " Grid metadata field '", label, "' must be a positive integer (got ", val, ").")
    }
    iv
  }

  x_int <- validate_dim(xbins_eff, "xbins_eff")
  y_int <- validate_dim(ybins_eff, "ybins_eff")
  n <- as.double(x_int) * as.double(y_int)
  if (!is.finite(n)) {
    stop(context, " Grid resolution produced a non-finite cell count.")
  }
  if (n > .Machine$integer.max) {
    stop(context, " Grid resolution ", x_int, " x ", y_int, " exceeds integer limits for neighbourhood builder.")
  }
  if (n > 5e8) {
    stop(context, " Grid resolution ", x_int, " x ", y_int, " exceeds supported upper bound (5e8 cells).")
  }

  list(xbins_eff = x_int, ybins_eff = y_int, n_cells = n)
}

#' Pure R rectangular neighbourhood builder (queen/rook)
#' @description
#' Internal helper for `.grid_nb_r`.
#' Provides a safe fallback when native .grid_nb is unavailable or unstable.
#' @param nrow Grid dimensions.
#' @param ncol Grid dimensions.
#' @param queen Logical, TRUE for queen (8-neighbour), FALSE for rook (4-neighbour).
#' @return List of integer neighbour vectors with `nb` attributes set.
#' @keywords internal
.grid_nb_r <- function(nrow, ncol, queen = TRUE) {
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  n <- nrow * ncol
  mat <- matrix(seq_len(n), nrow = nrow, ncol = ncol, byrow = TRUE)

  offsets <- list(
    c(0L, -1L), c(0L, 1L), c(-1L, 0L), c(1L, 0L)
  )
  if (isTRUE(queen)) {
    offsets <- c(offsets, list(
      c(-1L, -1L), c(-1L, 1L), c(1L, -1L), c(1L, 1L)
    ))
  }

  edges_src <- integer(0)
  edges_tgt <- integer(0)
  for (off in offsets) {
    dr <- off[[1]]
    dc <- off[[2]]
    r_src <- seq.int(max(1L, 1L - dr), min(nrow, nrow - dr))
    c_src <- seq.int(max(1L, 1L - dc), min(ncol, ncol - dc))
    if (!length(r_src) || !length(c_src)) next
    src <- mat[r_src, c_src, drop = FALSE]
    tgt <- mat[r_src + dr, c_src + dc, drop = FALSE]
    edges_src <- c(edges_src, as.integer(src))
    edges_tgt <- c(edges_tgt, as.integer(tgt))
  }

  split_tgt <- split(edges_tgt, edges_src)
  nb <- vector("list", n)
  if (length(split_tgt)) {
    for (nm in names(split_tgt)) {
      idx <- as.integer(nm)
      nb[[idx]] <- as.integer(split_tgt[[nm]])
    }
  }
  empty_idx <- which(vapply(nb, is.null, logical(1L)))
  if (length(empty_idx)) {
    nb[empty_idx] <- replicate(length(empty_idx), integer(0), simplify = FALSE)
  }

  attr(nb, "class") <- "nb"
  attr(nb, "region.id") <- seq_len(n)
  attr(nb, "queen") <- isTRUE(queen)
  attr(nb, "topology") <- if (isTRUE(queen)) "queen" else "rook"
  nb
}

#' Safe grid neighbourhood dispatcher with platform-aware fallback
#' @description
#' Internal helper for `.grid_nb_safe`.
#' Validates dimensions, honours env override GENESCOPE_FORCE_R_GRID_NB, and
#' uses the native C++ grid builders by default.
#' @param nrow Grid dimensions (validated externally).
#' @param ncol Grid dimensions (validated externally).
#' @param topology Topology label.
#' @param queen Logical, TRUE for queen adjacency.
#' @param use_parallel Logical, whether OpenMP backend was requested.
#' @param verbose Logical, emit selection logs.
#' @return List with `nb` and `backend` fields.
#' @keywords internal
.grid_nb_safe <- function(nrow,
                         ncol,
                         topology = "queen",
                         queen = TRUE,
                         use_parallel = FALSE,
                         verbose = FALSE,
                         parent = "computeWeights",
                         step = "S03") {
  force_r <- tolower(Sys.getenv("GENESCOPE_FORCE_R_GRID_NB", "")) %in% c("1", "true", "yes", "on")
  sysname <- tolower(Sys.info()[["sysname"]])
  fallback_supported <- topology %in% c("queen", "rook")
  fallback_reason <- NULL

  if (force_r) {
    fallback_reason <- "GENESCOPE_FORCE_R_GRID_NB=1"
  }

  if (!is.null(fallback_reason) && fallback_supported) {
    .log_info(parent, step, paste0(
      "Using R fallback grid builder (reason=", fallback_reason, ")."
    ), verbose)
    return(list(nb = .grid_nb_r(nrow, ncol, queen = queen), backend = "R"))
  }
  if (!is.null(fallback_reason) && !fallback_supported && verbose) {
    .log_info(parent, step, paste0(
      "R fallback requested (", fallback_reason,
      ") but topology=", topology,
      " is not supported; using native builder."
    ), verbose)
  }

  backend_label <- if (use_parallel) "omp" else "serial"
  if (verbose || force_r) {
    .log_info(parent, step, paste0(
      "Using native .grid_nb (backend=", backend_label,
      ", reason=", if (!is.null(fallback_reason)) "fallback_not_supported" else "default",
      ")."
    ), verbose)
  }
  nb <- tryCatch(
    {
      if (use_parallel) {
        .grid_nb_omp(nrow, ncol, queen = queen)
      } else {
        .grid_nb(nrow, ncol, queen = queen)
      }
    },
    error = function(e) {
      if (fallback_supported) {
        .log_info(parent, step, paste0(
          "Native .grid_nb failed (sysname=",
          sysname,
          "); falling back to R (error=",
          conditionMessage(e),
          ")."
        ), verbose)
        return(.grid_nb_r(nrow, ncol, queen = queen))
      }
      stop(e)
    }
  )
  list(nb = nb, backend = backend_label)
}

#' Spw Is Integerish
#' @description
#' Internal helper for `.spw_is_integerish`.
#' @param value Parameter value.
#' @param tol Parameter value.
#' @return Return value used internally.
#' @keywords internal
.spw_is_integerish <- function(value, tol = 1e-6) {
  is.numeric(value) &&
    length(value) == 1L &&
    is.finite(value) &&
    abs(value - round(value)) <= tol
}

#' Diffuse Binary Pattern
#' @description
#' Internal helper for `.diffuse_binary_pattern`.
#' @param x Parameter value.
#' @param W Parameter value.
#' @param alpha Parameter value.
#' @param steps Parameter value.
#' @param clamp Parameter value.
#' @return Return value used internally.
#' @keywords internal
.diffuse_binary_pattern <- function(x, W, alpha = 0.5, steps = 1L, clamp = TRUE) {
  v <- as.numeric(x)
  if (steps <= 0L || alpha <= 0) return(v)
  for (s in seq_len(steps)) {
    v <- v + alpha * as.numeric(W %*% v)
    if (clamp) v <- pmin(1, pmax(0, v))
  }
  v
}

#' Fuzzy Queen Diffuse Mat
#' @description
#' Internal helper for `.fuzzy_queen_diffuse_mat`.
#' @param X Parameter value.
#' @param W Parameter value.
#' @param alpha Parameter value.
#' @param steps Parameter value.
#' @return Return value used internally.
#' @keywords internal
.fuzzy_queen_diffuse_mat <- function(X, W, alpha = 0.5, steps = 1L) {
  V <- as.matrix(X)
  if (steps <= 0L || alpha <= 0) return(V)
  for (s in seq_len(steps)) {
    V <- V + alpha * as.matrix(W %*% V)
    V <- pmin(1, pmax(0, V))
  }
  V
}

#' Diffuse a signal over a fuzzy queen adjacency
#' @param X Numeric matrix to diffuse (rows = observations).
#' @param W Sparse `dgCMatrix` queen adjacency.
#' @param alpha Diffusion weight per step.
#' @param steps Number of diffusion steps.
#' @param row_normalize Whether to row-normalise `W` before diffusion.
#' @return Dense matrix of diffused values.
#' @keywords internal
fuzzy_queen_diffuse_mat <- function(X, W, alpha = 0.5, steps = 1L, row_normalize = TRUE) {
  if (!inherits(W, "dgCMatrix")) stop("W must be a dgCMatrix.")
  V <- as.matrix(X)
  if (nrow(W) != nrow(V)) stop("Row count of X must match dimension of W.")
  W_use <- if (isTRUE(row_normalize)) .row_normalize_sparse(W) else W
  .fuzzy_queen_diffuse_mat(V, W_use, alpha = alpha, steps = steps)
}

#' Fuzzy-queen Jaccard similarity between two binary patterns
#' @param x Numeric/logical vectors (0/1) of identical length.
#' @param y Numeric/logical vectors (0/1) of identical length.
#' @param W Sparse \code{dgCMatrix} queen adjacency.
#' @param alpha Diffusion weight per step (0–1).
#' @param steps Integer number of diffusion steps (≥ 0).
#' @param row_normalize Logical; row-normalise \code{W} before diffusion.
#' @return Scalar numeric fuzzy Jaccard similarity in \eqn{[0, 1]} (or \code{NA}
#'   when both patterns are empty).
#' @keywords internal
fuzzy_queen_jaccard <- function(x, y,
                                W,
                                alpha = 0.5,
                                steps = 1L,
                                row_normalize = TRUE) {
  if (length(x) != length(y)) stop("x and y must have identical length.")
  n <- length(x)
  if (!inherits(W, "dgCMatrix")) stop("W must be a dgCMatrix.")
  if (nrow(W) != n || ncol(W) != n) stop("Dimensions of W must match length of x and y.")
  if (!is.numeric(alpha) || alpha < 0) stop("alpha must be non-negative.")
  if (!is.numeric(steps) || steps < 0) stop("steps must be a non-negative integer.")

  x_vec <- as.numeric(x > 0)
  y_vec <- as.numeric(y > 0)

  W_use <- if (isTRUE(row_normalize)) .row_normalize_sparse(W) else W

  x_fuzzy <- .diffuse_binary_pattern(x_vec, W_use, alpha = alpha, steps = steps)
  y_fuzzy <- .diffuse_binary_pattern(y_vec, W_use, alpha = alpha, steps = steps)

  num <- sum(pmin(x_fuzzy, y_fuzzy))
  den <- sum(pmax(x_fuzzy, y_fuzzy))
  if (den == 0) return(NA_real_)
  as.numeric(num / den)
}

#' Spw Analyse Hex Geometry
#' @description
#' Internal helper for `.spw_analyse_hex_geometry`.
#' @param gi Parameter value.
#' @return Return value used internally.
#' @keywords internal
.spw_analyse_hex_geometry <- function(gi) {
  mk <- function(dr, dc) {
    b <- gi[, .(gx_nb = gx - dr, gy_nb = gy - dc, cx_nb = center_x, cy_nb = center_y)]
    m <- merge(gi, b, by.x = c("gx", "gy"), by.y = c("gx_nb", "gy_nb"), all = FALSE)
    if (!nrow(m)) return(NULL)
    data.table(dx = m$cx_nb - m$center_x, dy = m$cy_nb - m$center_y)
  }
  edges_xy <- rbindlist(list(
    mk(0, 1), mk(1, 0), mk(1, 1), mk(1, -1)
  ), use.names = TRUE, fill = TRUE)

  if (is.null(edges_xy) || !nrow(edges_xy)) {
    return(list(hex_like = FALSE, horiz = NA_real_, vert = NA_real_, diag_ratio = NA_real_))
  }

  theta <- abs(atan2(edges_xy$dy, edges_xy$dx) * 180 / pi)
  horiz <- mean(pmin(theta, 180 - theta) <= 10, na.rm = TRUE)
  vert <- mean(abs(theta - 90) <= 10, na.rm = TRUE)
  hex_hits <- mean(abs(theta - 60) <= 10 | abs(theta - 120) <= 10, na.rm = TRUE)
  hex_hits <- if (is.finite(hex_hits)) hex_hits else 0

  neighbor_pairs <- function(dr, dc) {
    b2 <- gi[, .(gx_nb = gx - dr, gy_nb = gy - dc, cx_nb = center_x, cy_nb = center_y)]
    m2 <- merge(gi, b2, by.x = c("gx", "gy"), by.y = c("gx_nb", "gy_nb"), all = FALSE)
    if (!nrow(m2)) return(numeric())
    sqrt((m2$cx_nb - m2$center_x)^2 + (m2$cy_nb - m2$center_y)^2)
  }
  d_ax <- c(neighbor_pairs(0, 1), neighbor_pairs(1, 0))
  d_diag <- neighbor_pairs(1, 1)
  med_ax <- if (length(d_ax)) median(d_ax, na.rm = TRUE) else NA_real_
  med_diag <- if (length(d_diag)) median(d_diag, na.rm = TRUE) else NA_real_
  diag_ratio <- if (is.finite(med_ax) && med_ax > 0 && is.finite(med_diag)) med_diag / med_ax else NA_real_
  ratio_ok <- is.finite(diag_ratio) && abs(diag_ratio - 1) <= 0.2

  list(
    hex_like = (hex_hits >= 0.4) && max(horiz, vert, na.rm = TRUE) < 0.8 && ratio_ok,
    horiz = if (is.finite(horiz)) horiz else 0,
    vert = if (is.finite(vert)) vert else 0,
    diag_ratio = if (is.finite(diag_ratio)) diag_ratio else NA_real_
  )
}

#' Spw Choose Topology
#' @description
#' Internal helper for `.spw_choose_topology`.
#' @param scope_obj A `scope_object`.
#' @param layer Parameter value.
#' @param topology Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.spw_choose_topology <- function(scope_obj, layer, topology, verbose,
                                 parent = "computeWeights", step = "S02") {
  grid_info <- as.data.table(layer$grid_info)
  has_centers <- all(c("gx", "gy", "center_x", "center_y") %in% names(grid_info))

  platform_vals <- c(
    if (!is.null(scope_obj@meta.data) && nrow(scope_obj@meta.data) && "platform" %in% names(scope_obj@meta.data))
      as.character(scope_obj@meta.data$platform) else NULL,
    if (!is.null(scope_obj@stats$platform)) as.character(scope_obj@stats$platform) else NULL
  )
  is_visium <- any(grepl("visium", platform_vals, ignore.case = TRUE))
  is_xenium <- any(grepl("xenium", platform_vals, ignore.case = TRUE))

  resolve_hex <- function(force = FALSE, msg_prefix = NULL) {
    hex_info <- .spw_infer_hex_variant(grid_info, force = force)
    if (!is.null(msg_prefix) && verbose) {
      .log_info(parent, step, paste0(msg_prefix, hex_info$variant), verbose)
    } else if (force && verbose && isFALSE(hex_info$hex_like)) {
      .log_info(parent, step, paste0(
        "Hex topology requested; defaulting to ",
        hex_info$variant,
        " due to ambiguous geometry."
      ), verbose)
    } else if (verbose && isTRUE(hex_info$hex_like)) {
      .log_info(parent, step, paste0("Auto-detected hex topology: ", hex_info$variant), verbose)
    }
    list(topology = hex_info$variant, base = "hex")
  }

  if (!identical(topology, "auto")) {
    if (identical(topology, "hex")) {
      return(resolve_hex(force = TRUE))
    }
    if (verbose) .log_info(parent, step, paste0("Topology forced by argument: ", topology), verbose)
    return(list(topology = topology, base = topology))
  }

  if (is_xenium) {
    if (verbose) .log_info(parent, step, "Platform flagged as Xenium; forcing queen adjacency.", verbose)
    return(list(topology = "queen", base = "queen"))
  }

  if (is_visium && has_centers) {
    return(resolve_hex(force = TRUE, msg_prefix = "Platform flagged as Visium; selecting hex orientation: "))
  }

  if (has_centers) {
    return(.spw_detect_topology(grid_info, verbose, parent = parent, step = step))
  }

  if (verbose) {
    .log_info(parent, step, "Defaulting to queen adjacency (insufficient geometry for auto detection).", verbose)
  }
  list(topology = "queen", base = "queen")
}

#' Spw Detect Topology
#' @description
#' Internal helper for `.spw_detect_topology`.
#' @param grid_info Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.spw_detect_topology <- function(grid_info, verbose, parent = "computeWeights", step = "S02") {
  hex_info <- .spw_infer_hex_variant(grid_info, force = FALSE)
  if (isTRUE(hex_info$hex_like)) {
    return(list(topology = hex_info$variant, base = "hex"))
  }
  if (verbose) .log_info(parent, step, "Auto-detected topology: queen", verbose)
  list(topology = "queen", base = "queen")
}

#' Spw Expand Neighbourhood
#' @description
#' Internal helper for `.spw_expand_neighbourhood`.
#' @param nb_obj Parameter value.
#' @param nb_order Parameter value.
#' @param min_neighbors Numeric threshold.
#' @param max_order Numeric threshold.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.spw_expand_neighbourhood <- function(nb_obj,
                                      nb_order = 1L,
                                      min_neighbors = NULL,
                                      max_order = NULL,
                                      verbose = FALSE,
                                      parent = "computeWeights",
                                      step = "S03") {
  order_use <- as.integer(nb_order)
  if (!is.finite(order_use) || order_use < 1L) {
    order_use <- 1L
  }

  if (is.null(max_order)) {
    max_order <- max(order_use, 4L)
  } else {
    max_order <- as.integer(max_order)
    if (!is.finite(max_order) || max_order < order_use) {
      max_order <- order_use
    }
  }

  if (!is.null(min_neighbors)) {
    min_neighbors <- as.numeric(min_neighbors)
    if (!is.finite(min_neighbors) || min_neighbors <= 0) min_neighbors <- NULL
  }

  expand_once <- function(ord) {
    if (ord <= 1L) {
      nb_obj
    } else {
      spdep::nblag_cumul(nb_obj, ord)
    }
  }

  clean_nb <- function(nb_list) {
    out <- lapply(nb_list, function(v) {
      if (length(v)) {
        unique(as.integer(v))
      } else {
        integer(0)
      }
    })
    attr(out, "class") <- "nb"
    out
  }

  nb_use <- clean_nb(expand_once(order_use))
  avg_deg <- mean(vapply(nb_use, length, integer(1)))

  if (!is.null(min_neighbors)) {
    while (avg_deg < min_neighbors && order_use < max_order) {
      order_use <- order_use + 1L
      nb_use <- clean_nb(expand_once(order_use))
      avg_deg <- mean(vapply(nb_use, length, integer(1)))
    }
  }

  attr(nb_use, "class") <- "nb"
  attr(nb_use, "region.id") <- attributes(nb_obj)[["region.id"]]
  attr(nb_use, "queen") <- attributes(nb_obj)[["queen"]]
  attr(nb_use, "topology") <- attributes(nb_obj)[["topology"]]
  attr(nb_use, "nb_order") <- order_use
  attr(nb_use, "call") <- NULL

  if (verbose && (order_use > 1L || (!is.null(min_neighbors) && order_use != nb_order))) {
    base_avg <- mean(vapply(nb_obj, length, integer(1)))
    msg <- paste0(
      "Using neighbourhood order ", order_use,
      " (avg neighbours ", round(avg_deg, 1),
      if (!is.null(min_neighbors)) paste0("; target ", min_neighbors, ", cap ", max_order) else "",
      if (order_use > 1L && is.finite(base_avg)) paste0("; base avg ", round(base_avg, 1)) else "",
      ")"
    )
    .log_info(parent, step, msg, verbose)
  }

  list(nb = nb_use, order_used = order_use, avg_neighbors = avg_deg)
}

#' Spw Hex Col Parity
#' @description
#' Internal helper for `.spw_hex_col_parity`.
#' @param gi Parameter value.
#' @return Return value used internally.
#' @keywords internal
.spw_hex_col_parity <- function(gi) {
  col_diffs <- gi[order(gx, center_y), diff(center_y), by = gx]$V1
  col_diffs <- col_diffs[is.finite(col_diffs)]
  if (!length(col_diffs)) return("hex-oddq")
  cell_height <- median(abs(col_diffs), na.rm = TRUE)
  if (!is.finite(cell_height) || cell_height == 0) return("hex-oddq")
  half_height <- cell_height / 2

  col_shift <- gi[, .(min_y = min(center_y)), by = gx]
  global_min <- min(col_shift$min_y)
  col_shift[, shift := min_y - global_min]
  col_shift[, ratio := shift / half_height]
  col_shift[, nearest := round(ratio)]
  col_shift[, diff := abs(ratio - nearest)]
  col_shift[diff > 0.35, nearest := NA_integer_]
  col_shift[, flag := nearest %% 2]

  odd_mean <- mean(col_shift[gx %% 2 == 1L, flag], na.rm = TRUE)
  even_mean <- mean(col_shift[gx %% 2 == 0L, flag], na.rm = TRUE)

  if (is.nan(odd_mean) && is.nan(even_mean)) return("hex-oddq")
  if (!is.nan(odd_mean) && (is.nan(even_mean) || odd_mean >= 0.5)) return("hex-oddq")
  if (!is.nan(even_mean) && even_mean >= 0.5) return("hex-evenq")

  if (!is.nan(odd_mean) && odd_mean > 0.25) return("hex-oddq")
  "hex-evenq"
}

#' Spw Hex Row Parity
#' @description
#' Internal helper for `.spw_hex_row_parity`.
#' @param gi Parameter value.
#' @return Return value used internally.
#' @keywords internal
.spw_hex_row_parity <- function(gi) {
  row_diffs <- gi[order(gy, center_x), diff(center_x), by = gy]$V1
  row_diffs <- row_diffs[is.finite(row_diffs)]
  if (!length(row_diffs)) return("hex-oddr")
  cell_width <- median(abs(row_diffs), na.rm = TRUE)
  if (!is.finite(cell_width) || cell_width == 0) return("hex-oddr")
  half_width <- cell_width / 2

  row_shift <- gi[, .(min_x = min(center_x)), by = gy]
  global_min <- min(row_shift$min_x)
  row_shift[, shift := min_x - global_min]
  row_shift[, ratio := shift / half_width]
  row_shift[, nearest := round(ratio)]
  row_shift[, diff := abs(ratio - nearest)]
  row_shift[diff > 0.35, nearest := NA_integer_]
  row_shift[, flag := nearest %% 2]

  odd_mean <- mean(row_shift[gy %% 2 == 1L, flag], na.rm = TRUE)
  even_mean <- mean(row_shift[gy %% 2 == 0L, flag], na.rm = TRUE)

  if (is.nan(odd_mean) && is.nan(even_mean)) return("hex-oddr")
  if (!is.nan(odd_mean) && (is.nan(even_mean) || odd_mean >= 0.5)) return("hex-oddr")
  if (!is.nan(even_mean) && even_mean >= 0.5) return("hex-evenr")

  if (!is.nan(odd_mean) && odd_mean > 0.25) return("hex-oddr")
  "hex-evenr"
}

#' Spw Infer Hex Variant
#' @description
#' Internal helper for `.spw_infer_hex_variant`.
#' @param grid_info Parameter value.
#' @param force Parameter value.
#' @return Return value used internally.
#' @keywords internal
.spw_infer_hex_variant <- function(grid_info, force = FALSE) {
  gi <- as.data.table(grid_info)[, .(gx, gy, center_x, center_y)]
  analysis <- .spw_analyse_hex_geometry(gi)
  hex_like <- isTRUE(analysis$hex_like)

  orientation <- if (hex_like || force) {
    if (is.finite(analysis$horiz) && is.finite(analysis$vert)) {
      if (analysis$horiz >= analysis$vert) "row" else "col"
    } else {
      "row"
    }
  } else {
    NA_character_
  }

  if (is.na(orientation)) {
    return(list(hex_like = hex_like, variant = "hex-oddr", orientation = orientation))
  }

  variant <- if (orientation == "row") {
    .spw_hex_row_parity(gi)
  } else {
    .spw_hex_col_parity(gi)
  }

  list(hex_like = hex_like, variant = variant, orientation = orientation)
}

#' Spw Select Layer
#' @description
#' Internal helper for `.spw_select_layer`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.spw_select_layer <- function(scope_obj, grid_name, verbose = FALSE) {
  layer <- .select_grid_layer(scope_obj, grid_name)
  if (is.null(grid_name)) {
    matches <- vapply(scope_obj@grid, identical, logical(1), layer)
    grid_name <- names(scope_obj@grid)[which(matches)][1]
  }
  if (is.null(grid_name) || !nzchar(grid_name)) {
    stop("Unable to resolve grid_name for spatial weights.")
  }
  list(layer = layer, grid_name = grid_name)
}

#' Spw Thread Guard V2
#' @description
#' Internal helper for `.spw_thread_guard_v2`.
#' @param builder_threads Parameter value.
#' @param requested_ncores Parameter value.
#' @param builder_backend Parameter value.
#' @param other_threads Parameter value.
#' @param context String label used for tracing/debug messages.
#' @return Return value used internally.
#' @keywords internal
.spw_thread_guard_v2 <- function(builder_threads = 1L,
                              requested_ncores = detectCores(logical = TRUE),
                              builder_backend = "serial",
                              other_threads = NULL,
                              context = "[geneSCOPE::.compute_weights]") {
  threads_use <- suppressWarnings(as.integer(builder_threads))
  if (!length(threads_use) || !is.finite(threads_use) || threads_use < 1L) {
    threads_use <- 1L
  }
  thread_config <- .configure_threads_for("openmp_only",
    ncores_requested = threads_use,
    restore_after = TRUE
  )
  restore_threads <- attr(thread_config, "restore_function")
  before <- attr(thread_config, "state_before")
  after <- .save_thread_state()

  if (.thread_guard_trace_enabled()) {
    .thread_guard_trace(
      context = paste0(context, "::nb(", builder_backend, ")"),
      before = if (!is.null(before)) before else .save_thread_state(),
      after = after,
      arrow_threads = thread_config[["arrow_threads"]],
      builder_threads = threads_use,
      requested_ncores = requested_ncores,
      other_threads = other_threads
    )
  }

  list(
    restore = function() {
      if (!is.null(restore_threads)) restore_threads()
    },
    config = thread_config
  )
}

#' Spw Validate Grid Metadata
#' @description
#' Internal helper for `.spw_validate_grid_metadata`.
#' @param layer Parameter value.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @return Return value used internally.
#' @keywords internal
.spw_validate_grid_metadata <- function(layer, grid_name) {
  grid_info <- layer$grid_info
  required_cols <- c("gx", "gy", "grid_id")
  if (!all(required_cols %in% names(grid_info))) {
    stop("grid_info for layer ", grid_name, " must contain columns: gx, gy, grid_id.")
  }
  if (is.null(layer$xbins_eff) || is.null(layer$ybins_eff)) {
    stop("grid layer missing xbins_eff / ybins_eff metadata")
  }
  gx <- grid_info$gx
  gy <- grid_info$gy
  xbins_eff <- layer$xbins_eff
  ybins_eff <- layer$ybins_eff
  if (any(gx < 1 | gx > xbins_eff | gy < 1 | gy > ybins_eff)) {
    stop("gx/gy exceed xbins_eff × ybins_eff range; check coordinates.")
  }
  list(
    grid_info = as.data.table(grid_info),
    grid_id = grid_info$grid_id,
    gx = gx,
    gy = gy,
    xbins_eff = xbins_eff,
    ybins_eff = ybins_eff
  )
}

# Baseline-compatible thread guard: preserve safe tracing while matching legacy signature
.spw_thread_guard_v2_internal <- .spw_thread_guard_v2

#' Guard thread configuration for spatial weights builders (baseline signature).
#' @description
#' Internal helper for `.spw_thread_guard_v2`.
#' Restores BLAS/OpenMP thread counts after .grid_nb execution; delegates to
#' -safe guard for traceability.
#' @keywords internal
#' @return Return value used internally.
.spw_thread_guard_v2 <- function() {
  guard <- .spw_thread_guard_v2_internal(
    builder_threads = 1L,
    requested_ncores = detectCores(logical = TRUE),
    builder_backend = "serial",
    other_threads = NULL,
    context = "[geneSCOPE::.compute_weights]"
  )

  blas_restore <- NULL
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA_integer_)
    RhpcBLASctl::blas_set_num_threads(1L)
    blas_restore <- function() {
      if (!is.na(old_threads)) {
        RhpcBLASctl::blas_set_num_threads(old_threads)
      }
    }
  }

  list(
    restore = function() {
      if (!is.null(blas_restore)) blas_restore()
      if (!is.null(guard$restore)) guard$restore()
    },
    config = guard$config
  )
}

# Compatibility shim for legacy thread guard name.
if (!exists(".spw_thread_guard", inherits = FALSE)) {
  .spw_thread_guard <- .spw_thread_guard_v2
}

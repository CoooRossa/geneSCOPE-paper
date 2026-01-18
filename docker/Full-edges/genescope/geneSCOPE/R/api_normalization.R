#' Normalize grid counts to Z-score expression.
#' @description
#' Builds a grid-by-gene count matrix from a grid layer, optionally performs
#' row normalization, converts to dense, and applies column-wise Z scoring
#' (sample SD) to produce an `Xz` expression layer.
#' @param scope_obj A `scope_object` that already contains a grid layer.
#' @param grid_name Target grid sub-layer; if `NULL`, the sole sub-layer is used.
#' @param keep_zero_grids Whether to keep grids whose counts are all zero.
#' @param row_norm Whether to perform row normalization before Z scoring.
#' @param zero_var_to_zero Whether to zero columns with standard deviation <= `var_eps`.
#' @param var_eps Standard deviation threshold below which columns are set to zero.
#' @param mem_dense_limit Threshold for `nrow * ncol`; values above this still
#'   coerce to dense with a warning (default `2e11`, roughly 16 GB for doubles).
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object` (invisibly) containing a new `Xz` layer.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir", grid_length = 30)
#' scope_obj <- normalizeMoleculesInGrid(scope_obj, grid_name = "grid30")
#' }
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom data.table as.data.table
#' @export
normalizeMoleculesInGrid <- function(scope_obj,
                                     grid_name = NULL,
                                     keep_zero_grids = FALSE,
                                     row_norm = TRUE,
                                     zero_var_to_zero = TRUE,
                                     var_eps = 0,
                                     mem_dense_limit = 2e11,
                                     verbose = TRUE) {
  ## ---- 0. Layer check -------------------------------------------------
  parent <- "normalizeMoleculesInGrid"
  step_s01 <- .log_step(parent, "S01", "select grid layer", verbose)
  step_s01$enter(extra = paste0(
    "grid_name=", if (is.null(grid_name)) "auto" else grid_name,
    " keep_zero_grids=", keep_zero_grids
  ))

  stopifnot(!is.null(scope_obj@grid))
  grid_layer_name <- if (is.null(grid_name)) {
    nms <- names(scope_obj@grid)
    if (length(nms) == 1L) {
      nms
    } else {
      stop("Multiple sub-layers; specify `grid_name`.")
    }
  } else {
    as.character(grid_name)
  }

  g <- scope_obj@grid[[grid_layer_name]]
  stopifnot(all(c("counts", "grid_info") %in% names(g)))
  step_s01$done(extra = paste0("grid=", grid_layer_name))

  ## ---- 1. Sparse counts → dgCMatrix ----------------------------------
  step_s02 <- .log_step(parent, "S02", "build sparse count matrix", verbose)
  step_s02$enter()
  dt <- as.data.table(g$counts)

  grids <- if (keep_zero_grids) {
    sort(unique(g$grid_info$grid_id))
  } else {
    sort(unique(dt$grid_id))
  }
  genes <- sort(unique(dt$gene))

  .log_info(parent, "S02",
    paste0("Matrix dimensions: ", length(grids), " grids x ", length(genes), " genes"),
    verbose
  )
  .log_backend(parent, "S02", "R",
    "Matrix sparse dgCMatrix build", verbose = verbose
  )

  X <- sparseMatrix(
    i = match(dt$grid_id, grids),
    j = match(dt$gene, genes),
    x = dt$count,
    dims = c(length(grids), length(genes)),
    dimnames = list(grids, genes)
  )
  step_s02$done(extra = paste0("n_grids=", length(grids), " n_genes=", length(genes)))

  ## ---- 2. Row normalization ------------------------------------------
  step_s03 <- .log_step(parent, "S03", "row normalization", verbose)
  step_s03$enter(extra = paste0("row_norm=", row_norm))
  if (row_norm) {
    .log_info(parent, "S03", "Applying row normalization...", verbose)
    lib <- rowSums(X)
    lib[lib == 0] <- 1
    X@x <- X@x / lib[X@i + 1L]
  } else {
    .log_info(parent, "S03", "Row normalization disabled; skipping.", verbose)
  }
  step_s03$done()

  ## ---- 3. Dense center + sample SD -----------------------------------
  elems <- prod(dim(X))
  step_s04 <- .log_step(parent, "S04", "dense standardization", verbose)
  step_s04$enter(extra = paste0("dense_elements=", elems))
  .log_backend(parent, "S04", "R",
    "Matrix sparse_to_dense=TRUE scale=base_scale", verbose = verbose
  )
  .log_info(parent, "S04", "Converting to dense matrix for standardization...", verbose)
  if (elems > mem_dense_limit) {
    warning(
      "Matrix size ", format(elems, scientific = FALSE),
      " > mem_dense_limit; coercing to dense may use much RAM."
    )
  }

  Xd <- as.matrix(X) # Dense copy
  .log_info(parent, "S04", "Applying column-wise standardization...", verbose)
  Xd <- scale(Xd, center = TRUE, scale = TRUE) # base::scale → sample SD

  Xd[is.na(Xd)] <- 0 # Constant columns NA → 0
  step_s04$done()

  ## ---- 4. Zero near-variance columns as requested --------------------
  step_s05 <- .log_step(parent, "S05", "store Xz output", verbose)
  step_s05$enter(extra = paste0("zero_var_to_zero=", zero_var_to_zero, " var_eps=", var_eps))
  if (zero_var_to_zero) {
    .log_info(parent, "S05", "Processing zero-variance columns...", verbose)
    csd <- attr(Xd, "scaled:scale")
    if (is.null(csd)) { # R 4.4 scale() no longer attaches attributes, can re-compute
      csd <- apply(Xd, 2, sd)
    }
    Xd[, csd <= var_eps] <- 0
  } else {
    .log_info(parent, "S05", "zero_var_to_zero=FALSE; skipping.", verbose)
  }

  ## ---- 5. Save & return ----------------------------------------------
  g$Xz <- Xd
  scope_obj@grid[[grid_layer_name]] <- g

  .log_info(parent, "S05", "Normalization completed, Xz layer stored", verbose)
  step_s05$done()

  invisible(scope_obj)
}

#' Add a log-CPM normalized layer to a scope object.
#' @description
#' Reads the specified count layer from `scope_obj@cells`, scales counts to
#' `scale_factor`, applies `log1p`, and stores the result in
#' `scope_obj@cells[[output_layer]]`.
#' @param scope_obj A `scope_object`.
#' @param input_layer Name of the existing layer in `scope_obj@cells`
#'   (default `counts`, must be a `dgCMatrix`).
#' @param output_layer Name of the layer to write to (default `logCPM`,
#'   overwrites if it exists).
#' @param scale_factor Library size each cell is scaled to before `log1p`.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object` (invisibly).
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir",　grid_length = 30)
#' scope_obj <- addSingleCells(scope_obj, data_dir = "/path/to/output_dir")
#' scope_obj <- normalizeSingleCells(scope_obj, input_layer = "counts", output_layer = "logCPM", scale_factor = 1e4)
#' }
#' @import Matrix
#' @export
normalizeSingleCells <- function(scope_obj,
                                 input_layer = "counts",
                                 output_layer = "logCPM",
                                 scale_factor = 1e4,
                                 verbose = TRUE) {
  parent <- "normalizeSingleCells"
  step_s01 <- .log_step(parent, "S01", "validate cell layer", verbose)
  step_s01$enter(extra = paste0("input_layer=", input_layer, " output_layer=", output_layer))

  stopifnot(
    inherits(scope_obj, "scope_object"),
    is.list(scope_obj@cells),
    input_layer %in% names(scope_obj@cells)
  )

  mat <- scope_obj@cells[[input_layer]]
  if (!inherits(mat, "dgCMatrix")) {
    stop("Layer '", input_layer, "' is not a dgCMatrix.")
  }
  .log_backend(parent, "S01", "R",
    "Matrix input=dgCMatrix sparse=TRUE", verbose = verbose
  )
  step_s01$done(extra = paste0("n_genes=", nrow(mat), " n_cells=", ncol(mat)))

  step_s02 <- .log_step(parent, "S02", "compute library sizes", verbose)
  step_s02$enter()
  .log_info(parent, "S02", "Computing library sizes...", verbose)
  lib <- colSums(mat) # Library size per column
  zerocell <- lib == 0
  lib[zerocell] <- 1 # Avoid division by zero
  step_s02$done(extra = paste0("zero_cells=", sum(zerocell)))

  step_s03 <- .log_step(parent, "S03", "apply log-CPM transform", verbose)
  step_s03$enter(extra = paste0("scale_factor=", scale_factor))
  .log_backend(parent, "S03", "R",
    "Matrix log1p CPM sparse=TRUE", verbose = verbose
  )
  .log_info(parent, "S03", "Applying CPM transformation and log1p...", verbose)
  ## Vectorized: match corresponding column factors for each non-zero element
  col_rep <- rep(seq_len(ncol(mat)), diff(mat@p))
  sf_vec <- (scale_factor / lib)[col_rep]
  out <- mat # Copy sparse structure
  out@x <- log1p(out@x * sf_vec)
  step_s03$done()

  if (any(zerocell)) {
    if (verbose) warning(sum(zerocell), " cells have total count = 0; they remain 0 in log-CPM.")
  }

  step_s04 <- .log_step(parent, "S04", "store output layer", verbose)
  step_s04$enter(extra = paste0("output_layer=", output_layer))
  scope_obj@cells[[output_layer]] <- out

  .log_info(parent, "S04",
    paste0("Log-CPM normalization completed, stored in @cells$", output_layer),
    verbose
  )
  step_s04$done()

  invisible(scope_obj)
}

#' Compute gene-gene correlations.
#' @description
#' Computes a gene × gene correlation matrix from either single-cell expression
#' (`level = "cell"`) or grid-level expression (`level = "grid"`), stores results
#' into the `scope_object`, and returns the updated object invisibly.
#' @param scope_obj A `scope_object` containing expression layers.
#' @param level Correlation level (`cell` or `grid`).
#' @param grid_name Optional grid name when operating on grid-level data.
#' @param layer Expression layer name to correlate (e.g., `logCPM`).
#' @param method Correlation method (`pearson`, `spearman`, `kendall`).
#' @param blocksize Block size passed to the C++ backend (tuning knob).
#' @param ncores Number of threads to use.
#' @param chunk_size Chunk size used when using file-backed computation.
#' @param memory_limit_gb In-memory size threshold (GB) above which file-backed mode is used.
#' @param use_bigmemory Whether to allow file-backed computation for large matrices.
#' @param backing_path Directory for file-backed artifacts (default `tempdir()`).
#' @param force_compute Whether to force recomputation even if cached outputs exist.
#' @param store_layer Slot name for storing the correlation matrix.
#' @param compute_fdr Whether to compute and store an FDR-adjusted p-value matrix.
#' @param fdr_method Multiple-testing adjustment method passed to `p.adjust()`.
#' @param verbose Emit progress messages when TRUE.
#' @return The updated `scope_object` (invisibly).
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir",　grid_length = 30)
#' scope_obj <- addSingleCells(scope_obj, data_dir = "/path/to/output_dir")
#' scope_obj <- normalizeSingleCells(scope_obj, input_layer = "counts", output_layer = "logCPM", scale_factor = 1e4)
#' scope_obj <- computeCorrelation(scope_obj, level = "cell", layer = "logCPM", method = "pearson", ncores = 16) 
#' scope_obj <- computeCorrelation(scope_obj, level = "grid", grid_name = "grid30", layer = "Xz", method = "pearson", ncores = 16)
#' }
#' @seealso `addSingleCells()`, `normalizeSingleCells()`, `computeL()`
#' @export
computeCorrelation <- function(
    scope_obj,
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
    .compute_correlation(
        scope_obj = scope_obj,
        level = level,
        grid_name = grid_name,
        layer = layer,
        method = method,
        blocksize = blocksize,
        ncores = ncores,
        chunk_size = chunk_size,
        memory_limit_gb = memory_limit_gb,
        use_bigmemory = use_bigmemory,
        backing_path = backing_path,
        force_compute = force_compute,
        store_layer = store_layer,
        compute_fdr = compute_fdr,
        fdr_method = fdr_method,
        verbose = verbose
    )
}

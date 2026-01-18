#' Compute Lee's L statistic.
#' @description
#' Computes Lee's L for gene pairs using a grid-level expression layer, stores
#' results under a statistics layer in `scope_obj@stats`, and returns the updated
#' object.
#' @param scope_obj A `scope_object` with at least one populated `@grid` slot.
#' @param grid_name Name of the grid layer to process (auto-selected if only one exists).
#' @param genes Optional subset of genes to include.
#' @param within If `TRUE`, restricts analysis to the selected gene set on both axes.
#' @param ncores Number of cores for parallel processing.
#' @param block_side Number of grid cells per side for block partitioning.
#' @param perms Number of permutations for Monte-Carlo p-values.
#' @param block_size Number of permutations processed per batch.
#' @param L_min Similarity threshold used when building QC similarity graphs.
#' @param norm_layer Name of the normalised expression layer (default `"Xz"`).
#' @param lee_stats_layer_name Output statistics layer name (auto-generated when NULL).
#' @param legacy_formula Use legacy denominator for compatibility.
#' @param mem_limit_GB RAM threshold that triggers streaming mode.
#' @param chunk_size Number of columns processed per chunk in streaming mode.
#' @param use_bigmemory Whether to use file-backed matrices for large computations.
#' @param backing_path Directory for temporary files (default `tempdir()`).
#' @param cache_inputs Whether to .cache preprocessed inputs for reuse across calls.
#' @param verbose Whether to emit progress messages.
#' @param ncore Deprecated alias of `ncores`.
#' @return The modified `scope_object`.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir", grid_length = 30)
#' scope_obj <- normalizeMoleculesInGrid(scope_obj, grid_name = "grid30")
#' scope_obj <- computeWeights(scope_obj, grid_name = "grid30")
#' scope_obj <- computeL(scope_obj, grid_name = "grid30", ncores = 16)
#' }
#' @seealso `computeWeights()`, `computeMH()`, `computeLvsRCurve()`, `getTopLvsR()`
#' @export
computeL <- function(
    scope_obj,
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
    .compute_l(
        scope_obj = scope_obj,
        grid_name = grid_name,
        genes = genes,
        within = within,
        ncores = ncores,
        block_side = block_side,
        perms = perms,
        block_size = block_size,
        L_min = L_min,
        norm_layer = norm_layer,
        lee_stats_layer_name = lee_stats_layer_name,
        legacy_formula = legacy_formula,
        mem_limit_GB = mem_limit_GB,
        chunk_size = chunk_size,
        use_bigmemory = use_bigmemory,
        backing_path = backing_path,
        cache_inputs = cache_inputs,
        verbose = verbose,
        ncore = ncore
    )
}

#' Compute L vs R curve.
#' @description
#' Computes a smoothed relationship between Lee's L and Pearson correlation using
#' bootstrapping and stores the fitted curve and confidence intervals into the
#' `scope_object`.
#' @param scope_obj A `scope_object` containing Lee's L and correlation matrices.
#' @param grid_name Grid layer name to operate on.
#' @param level Correlation level used for Pearson r (`grid` or `cell`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param span LOESS span parameter.
#' @param B Number of bootstrap iterations.
#' @param deg Degree for the LOESS fit.
#' @param ncores Number of threads to use.
#' @param length_out Grid size for the fitted curve.
#' @param downsample Downsampling control (ratio < 1, or target count >= 1).
#' @param n_strata Number of strata used when drawing bootstrap samples.
#' @param k_max Maximum number of points used per stratum (guards runtime).
#' @param jitter_eps Optional jitter factor applied to correlation values.
#' @param ci_method Confidence interval method (`percentile`, `basic`, `bc`).
#' @param ci_adjust Optional analytic adjustment for CI width.
#' @param min_rel_width Minimum relative CI width.
#' @param widen_span Additional span widening used when enforcing `min_rel_width`.
#' @param curve_name Name under which the curve is stored.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object` with stored curve data.
#' @examples
#' \dontrun{
#' scope_obj <- computeLvsRCurve(scope_obj, grid_name = "grid30", level = "cell", downsample = 0.1, ncores = 16)
#' }
#' @seealso `computeL()`, `computeCorrelation()`, `getTopLvsR()`
#' @export
computeLvsRCurve <- function(
    scope_obj,
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
    .compute_l_vs_r_curve(
        scope_obj = scope_obj,
        grid_name = grid_name,
        level = level,
        lee_stats_layer = lee_stats_layer,
        span = span,
        B = B,
        deg = deg,
        ncores = ncores,
        length_out = length_out,
        downsample = downsample,
        n_strata = n_strata,
        k_max = k_max,
        jitter_eps = jitter_eps,
        ci_method = ci_method,
        ci_adjust = ci_adjust,
        min_rel_width = min_rel_width,
        widen_span = widen_span,
        curve_name = curve_name,
        verbose = verbose
    )
}

#' Extract top L vs R pairs.
#' @description
#' Ranks gene pairs by the deviation between Lee's L and Pearson correlation and
#' returns a table of top pairs.
#' @param scope_obj A `scope_object` containing Lee's L and correlation matrices.
#' @param grid_name Grid layer name.
#' @param pear_level Correlation level (`cell` or `grid`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param expr_layer Optional expression layer name for gene prevalence filtering.
#' @param pear_range Range of Pearson r values to include.
#' @param L_range Range of Lee's L values to include.
#' @param top_n Number of top pairs to return.
#' @param direction Which tail to select (`largest`, `smallest`, `both`).
#' @param do_perm Whether to run permutation testing.
#' @param perms Number of permutations for p-value estimation.
#' @param block_side Block side length used when building permutation blocks.
#' @param use_blocks Whether to use spatial blocks for permutations.
#' @param ncores Number of threads to use.
#' @param clamp_mode Whether to clamp Delta using reference-only or both ends.
#' @param p_adj_mode Multiple-testing adjustment mode.
#' @param mem_limit_GB Memory threshold controlling chunking behavior.
#' @param pval_mode P-value mode used when transforming permutation counts.
#' @param curve_layer Optional curve layer from `computeLvsRCurve()`.
#' @param CI_rule Confidence-interval filtering rule (`remove_within`, `remove_outside`, `none`).
#' @param verbose Emit progress messages when TRUE.
#' @return A data.frame with ranked gene pairs and associated statistics.
#' @examples
#' \dontrun{
#' scope_obj <- computeL(scope_obj, grid_name = "grid30", ncores = 16)
#' scope_obj <- computeCorrelation(scope_obj, level = "cell", layer = "logCPM", ncores = 16)
#' top_pairs <- getTopLvsR(scope_obj, grid_name = "grid30", top_n = 200, ncores = 16, L_range = c(0.1, 1))
#' head(top_pairs)
#' }
#' @seealso `computeL()`, `computeCorrelation()`, `computeLvsRCurve()`
#' @export
getTopLvsR <- function(
    scope_obj,
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
    .get_top_l_vs_r(
        scope_obj = scope_obj,
        grid_name = grid_name,
        pear_level = pear_level,
        lee_stats_layer = lee_stats_layer,
        expr_layer = expr_layer,
        pear_range = pear_range,
        L_range = L_range,
        top_n = top_n,
        direction = direction,
        do_perm = do_perm,
        perms = perms,
        block_side = block_side,
        use_blocks = use_blocks,
        ncores = ncores,
        clamp_mode = clamp_mode,
        p_adj_mode = p_adj_mode,
        mem_limit_GB = mem_limit_GB,
        pval_mode = pval_mode,
        curve_layer = curve_layer,
        CI_rule = CI_rule,
        verbose = verbose
    )
}

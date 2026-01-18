#' Compute per-grid gene density.
#' @description
#' Aggregates counts for a gene set within each grid tile and stores density per
#' unit area under `scope_obj@density`.
#' @param scope_obj A `scope_object` with grid data.
#' @param grid_name Grid layer to use (defaults to active layer).
#' @param density_name Output column name for the density values.
#' @param genes Optional character vector of genes to aggregate.
#' @param cluster_col Optional `scope_obj@meta.data` column used to select genes by cluster.
#' @param cluster_num Optional cluster identifier used with `cluster_col`.
#' @param layer_name Grid layer slot containing counts (default `"counts"`).
#' @param normalize_method Normalization mode (`none`, `per_grid`, `global_gene`).
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object` (invisibly).
#' @examples
#' \dontrun{
#' scope_obj <- computeDensity(scope_obj, grid_name = "grid30", layer_name = "counts", genes = "CEACAM5", density_name = "CEACAM5")
#' }
#' @seealso `plotDensity()`, `plotDensityCentroids()`
#' @export
computeDensity <- function(
    scope_obj,
    grid_name = NULL,
    density_name = "density",
    genes = NULL,
    cluster_col = NULL,
    cluster_num = NULL,
    layer_name = "counts",
    normalize_method = c("none", "per_grid", "global_gene"),
    verbose = TRUE) {
    normalize_method <- match.arg(normalize_method)
    .compute_density(
        scope_obj = scope_obj,
        grid_name = grid_name,
        density_name = density_name,
        genes = genes,
        cluster_col = cluster_col,
        cluster_num = cluster_num,
        layer_name = layer_name,
        normalize_method = normalize_method,
        verbose = verbose
    )
}

#' Compute gene-gene Jaccard heatmap.
#' @description
#' Computes a gene × gene Jaccard similarity matrix from binary gene-by-bin
#' presence/absence and renders a heatmap when plotting packages are available.
#' @param scope_obj A `scope_object` with grid counts or an `Xz` matrix.
#' @param grid_name Grid layer name (default `"grid30"`).
#' @param expr_layer Expression layer to use (`Xz` or `counts`).
#' @param min_bins_per_gene Minimum bins per gene required to include.
#' @param top_n_genes Number of top genes to plot when no whitelist is given.
#' @param genes_of_interest Optional gene whitelist.
#' @param output_dir Optional directory to save outputs.
#' @param heatmap_bg Background color used for the heatmap.
#' @param cluster_col Optional `scope_obj@meta.data` column used to annotate genes by cluster.
#' @param keep_cluster_na Whether to keep genes with NA cluster labels.
#' @param split_heatmap_by_cluster Whether to split heatmap panels by cluster.
#' @return Invisibly returns a list containing the Jaccard matrix and (when available) the heatmap object.
#' @examples
#' \dontrun{
#' computeGeneGeneJaccardHeatmap(scope_obj, grid_name = "grid30", expr_layer = "Xz", cluster_col = "gene_module")
#' }
#' @seealso `computeGeneModuleWithinBetweenTest()`
#' @export
computeGeneGeneJaccardHeatmap <- function(
    scope_obj,
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
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    .compute_gene_gene_jaccard_heatmap(
        scope_obj = scope_obj,
        grid_name = grid_name,
        expr_layer = expr_layer,
        min_bins_per_gene = min_bins_per_gene,
        top_n_genes = top_n_genes,
        genes_of_interest = genes_of_interest,
        output_dir = output_dir,
        heatmap_bg = heatmap_bg,
        cluster_col = cluster_col,
        keep_cluster_na = keep_cluster_na,
        split_heatmap_by_cluster = split_heatmap_by_cluster,
        verbose = verbose
    )
}

#' Compute within/between module test.
#' @description
#' Compares within-module versus between-module Jaccard similarity distributions
#' using Wilcoxon testing or a sign-flip permutation test.
#' @param jac_mat Jaccard similarity matrix.
#' @param gene_modules Named vector assigning genes to modules.
#' @param alternative Alternative hypothesis (`greater`, `two.sided`, `less`).
#' @param method Testing method (`wilcoxon` or `permutation`).
#' @param n_perm Number of permutations when `method = "permutation"`.
#' @param seed Optional seed for reproducibility.
#' @param min_within Minimum within-module pairs required per gene.
#' @param min_between Minimum between-module pairs required per gene.
#' @return Invisibly returns a list with per-gene summaries and a test result.
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:10)
#' jac <- matrix(runif(100), nrow = 10, dimnames = list(genes, genes))
#' diag(jac) <- 1
#' mods <- setNames(rep(c("A", "B"), each = 5), genes)
#' res <- computeGeneModuleWithinBetweenTest(jac, mods, method = "wilcoxon")
#' @seealso `computeGeneGeneJaccardHeatmap()`
#' @export
computeGeneModuleWithinBetweenTest <- function(
    jac_mat,
    gene_modules,
    alternative = c("greater", "two.sided", "less"),
    method = c("wilcoxon", "permutation"),
    n_perm = 10000L,
    seed = NULL,
    min_within = 1L,
    min_between = 1L) {
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    .compute_gene_module_within_between_test(
        jac_mat = jac_mat,
        gene_modules = gene_modules,
        alternative = alternative,
        method = method,
        n_perm = n_perm,
        seed = seed,
        min_within = min_within,
        min_between = min_between,
        verbose = verbose
    )
}

#' Compute Morisita's Iδ statistic.
#' @description
#' Computes Moran's Iδ for each gene at grid or cell level and stores results
#' under `scope_obj@stats` and `scope_obj@meta.data`.
#' @param scope_obj A `scope_object` with grid or cell data.
#' @param grid_name Grid layer name when `level = "grid"`.
#' @param level Compute at `grid` or `cell` level.
#' @param ncores Number of threads to use.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object` (invisibly).
#' @examples
#' \dontrun{
#' scope_obj <- computeIDelta(scope_obj, grid_name = "grid30", level = "grid", ncores = 2)
#' }
#' @seealso `plotIDelta()`
#' @export
computeIDelta <- function(
    scope_obj,
    grid_name = NULL,
    level = c("grid", "cell"),
    ncores = 4,
    verbose = getOption("geneSCOPE.verbose", TRUE)) {
    level <- match.arg(level)
    .compute_idelta(
        scope_obj = scope_obj,
        grid_name = grid_name,
        level = level,
        ncores = ncores,
        verbose = verbose
    )
}

#' Compute Morisita-Horn statistics.
#' @description
#' Computes Morisita–Horn similarities between genes on a graph and stores the
#' resulting matrix and/or graph back into the scope object's stats layer.
#' @param scope_obj A `scope_object` containing Lee statistics.
#' @param grid_name Grid layer name.
#' @param lee_layer Lee statistics layer name.
#' @param graph_slot Graph slot holding the consensus network.
#' @param graph_slot_mh Slot name to store the Morisita–Horn graph.
#' @param matrix_slot_mh Optional slot to store the Morisita–Horn matrix.
#' @param out Output type: `matrix`, `igraph`, or `both`.
#' @param L_min Minimum Lee's L threshold for edges.
#' @param area_norm Whether to normalise by area.
#' @param ncores Number of threads to use.
#' @param use_chao Apply Chao correction when TRUE.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object` (invisibly).
#' @examples
#' \dontrun{
#' scope_obj <- computeMH(scope_obj, grid_name = "grid30", ncores = 16, graph_slot_mh = "g_morisita_30")
#' }
#' @seealso `computeL()`, `clusterGenes()`
#' @export
computeMH <- function(
    scope_obj,
    grid_name = NULL,
    lee_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    graph_slot_mh = "g_morisita",
    matrix_slot_mh = NULL,
    out = c("matrix", "igraph", "both"),
    L_min = 0,
    area_norm = TRUE,
    ncores = 8,
    use_chao = TRUE,
    verbose = TRUE) {
    .compute_mh(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_layer = lee_layer,
        graph_slot = graph_slot,
        graph_slot_mh = graph_slot_mh,
        matrix_slot_mh = matrix_slot_mh,
        out = out,
        L_min = L_min,
        area_norm = area_norm,
        ncores = ncores,
        use_chao = use_chao,
        verbose = verbose
    )
}

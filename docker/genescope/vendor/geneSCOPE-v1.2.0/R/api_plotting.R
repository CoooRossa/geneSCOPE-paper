#' Plot cluster assignments across methods.
#' @description
#' Creates a scatter-style comparison plot of cluster calls from multiple
#' methods for each gene.
#' @param scope_obj A `scope_object`.
#' @param method_cols Parameter value.
#' @param method_labels Parameter value.
#' @param point_size Parameter value.
#' @param palette Parameter value.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' p <- plotClusterComparison(scope_obj, method_cols = c("modL0.15", "modHS0.15"))
#' print(p)
#' }
#' @seealso `plotNetwork()`, `plotDendro()`
#' @export
plotClusterComparison <- function(
    scope_obj,
    method_cols,
    method_labels = method_cols,
    point_size = 3,
    palette = RColorBrewer::brewer.pal) {
    parent <- "plotClusterComparison"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve inputs and palette", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggplot2", reason = "cluster_comparison", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render comparison plot", verbose)
    step02$enter()
    plot_obj <- .plot_cluster_comparison(
        scope_obj = scope_obj,
        method_cols = method_cols,
        method_labels = method_labels,
        point_size = point_size,
        palette = palette
    )
    step02$done()
    plot_obj
}

#' Plot a dendrogram derived from the consensus graph.
#'
#' @description
#' Builds the same MST/forest skeleton as `plotDendroNetwork()` and renders a
#' dendrogram with optional cluster/annotation tip rows.
#' @param scope_obj A `scope_object` containing a LeeStats layer and a consensus graph.
#' @param grid_name Grid layer name (auto-selected when only one layer exists).
#' @param lee_stats_layer Lee statistics layer name (e.g. `"LeeStats_Xz"`).
#' @param graph_slot Slot name in the LeeStats layer holding the base graph (default `"g_consensus"`).
#' @param cluster_name `scope_obj@meta.data` column providing cluster memberships for genes.
#' @param cluster_ids Optional subset of cluster IDs to include.
#' @param IDelta_col_name Optional Iδ column used to personalize PageRank.
#' @param linkage Linkage method passed to `hclust()` (e.g. `"average"`).
#' @param plot_dend Whether to draw the dendrogram plot (TRUE) or only return objects.
#' @param weight_low_cut Edge-weight floor; edges below are treated as absent.
#' @param damping PageRank damping factor.
#' @param cluster_palette Palette used for cluster tip-dot row(s).
#' @param tip_dots Whether to draw tip-dot row for `cluster_name`.
#' @param tip_point_cex Tip dot size (graphics `cex` units).
#' @param tip_palette Optional explicit palette for tip dots (overrides `cluster_palette`).
#' @param cluster_name2 Optional second cluster column for a second tip-dot row.
#' @param tip_palette2 Optional palette for the second tip-dot row.
#' @param tip_row_offset2 Vertical offset for the second tip-dot row.
#' @param tip_label_offset Radial/vertical label offset.
#' @param tip_label_cex Tip label size (graphics `cex` units).
#' @param tip_label_adj Label adjustment (graphics `adj`).
#' @param tip_label_srt Label rotation (degrees).
#' @param tip_label_col Label color.
#' @param leaf_order Leaf ordering strategy (`"OLO"` or `"none"`).
#' @param length_mode Edge-length transform mode (`"neg_log"`, `"inverse"`, `"inverse_sqrt"`).
#' @param weight_normalize Whether to normalize weights before distance construction.
#' @param weight_clip_quantile Quantile used to clip extreme weights.
#' @param height_rescale Height rescaling mode (`"q95"`, `"max"`, `"none"`).
#' @param height_power Power applied to heights after rescaling.
#' @param height_scale Global multiplicative scale for heights.
#' @param distance_on Whether to compute distances on the `"tree"` skeleton or `"graph"`.
#' @param enforce_cluster_contiguity Whether to enforce contiguity of clusters in leaf order.
#' @param distance_smooth_power Optional smoothing power for distances (> 0).
#' @param tip_shape2 Shape used for the second tip-dot row.
#' @param tip_row1_label Legend label for tip row 1.
#' @param tip_row2_label Legend label for tip row 2.
#' @param tip_label_indent Indent fraction for labels.
#' @param legend_inline Whether to draw legend inline.
#' @param legend_files Optional file paths for legend export.
#' @param compose_outfile Optional output path for composed figure.
#' @param compose_width Output width for composed figure (pixels).
#' @param compose_height Output height for composed figure (pixels).
#' @param compose_res Output resolution (DPI) for composed figure.
#' @param legend_ncol1 Number of legend columns for tip row 1.
#' @param legend_ncol2 Number of legend columns for tip row 2.
#' @param title Plot title.
#' @return Invisibly returns a list containing dendrogram/network objects.
#' @examples
#' \dontrun{
#' out <- plotDendro(scope_obj, grid_name = "grid30", cluster_name = "modL0.15")
#' }
#' @seealso `plotDendroNetwork()`, `plotDendroNetworkMulti()`, `getDendroWalkPaths()`
#' @export
plotDendro <- function(
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
    parent <- "plotDendro"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    leaf_order_val <- if (length(leaf_order)) leaf_order[1] else leaf_order
    step01 <- .log_step(parent, "S01", "resolve inputs and ordering", verbose)
    step01$enter()
    .log_backend(parent, "S01", "leaf_order", leaf_order_val, reason = "ordering_strategy", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render dendrogram", verbose)
    step02$enter()
    plot_obj <- .plot_dendro(
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
    step02$done()
    plot_obj
}

#' Plot the dendrogram network derived from the consensus graph.
#'
#' @description
#' Builds a cluster-internal MST and cluster-level MST (optionally PageRank
#' reweighted) and renders the resulting network.
#' @param scope_obj A `scope_object` with a LeeStats layer containing a graph.
#' @param grid_name Grid layer name (auto-selected when only one layer exists).
#' @param lee_stats_layer Lee statistics layer name (e.g. `"LeeStats_Xz"`).
#' @param gene_subset Optional gene subset (character vector or selector supported by internal helpers).
#' @param L_min Minimum Lee's L threshold for keeping edges.
#' @param use_consensus_graph Whether to use the consensus graph stored in the stats layer.
#' @param graph_slot_name Graph slot name inside the stats layer (auto-detected when NULL).
#' @param cluster_vec Cluster label source; either a vector of labels or a meta.data column name.
#' @param IDelta_col_name Optional Iδ column used to personalize PageRank.
#' @param IDelta_invert Whether to invert Iδ personalization.
#' @param damping PageRank damping factor.
#' @param weight_low_cut Edge-weight floor; edges below are treated as absent.
#' @param cluster_palette Palette used for clusters.
#' @param node_size Node point size.
#' @param edge_width Edge width multiplier.
#' @param label_size Label size.
#' @param seed Random seed for layout reproducibility.
#' @param length_scale Global edge-length multiplier.
#' @param max.overlaps Maximum label overlaps for `ggrepel`.
#' @param hub_border_col Hub node border color.
#' @param hub_border_size Hub node border line width.
#' @param title Optional plot title.
#' @param k_top Number of high-weight inter-cluster edges to retain in addition to the MST.
#' @param tree_mode Tree layout mode (`"radial"` or `"rooted"`).
#' @return Invisibly returns a list containing the processed graph and plot-related artifacts.
#' @examples
#' \dontrun{
#' out <- plotDendroNetwork(scope_obj, grid_name = "grid30", cluster_vec = "modL0.15", graph_slot_name = "g_consensus")
#' }
#' @seealso `plotDendro()`, `plotDendroNetworkMulti()`
#' @export
plotDendroNetwork <- function(
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
    tree_mode = c("radial", "rooted")) {
    parent <- "plotDendroNetwork"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve inputs and graph source", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggraph", reason = "network_plot", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "build network plot", verbose)
    step02$enter()
    plot_obj <- .plot_dendro_network_impl(
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
        tree_mode = tree_mode
    )
    step02$done()
    plot_obj
}

#' Plot multiple stochastic dendrogram networks.
#'
#' @description
#' Repeats dendrogram network construction with jittered edge weights and
#' visualizes a union graph across runs.
#' @param scope_obj A `scope_object` with a LeeStats layer containing a graph.
#' @param grid_name Grid layer name (auto-selected when only one layer exists).
#' @param lee_stats_layer Lee statistics layer name (e.g. `"LeeStats_Xz"`).
#' @param gene_subset Optional gene subset (character vector or selector supported by internal helpers).
#' @param L_min Minimum Lee's L threshold for keeping edges.
#' @param use_consensus_graph Whether to use the consensus graph stored in the stats layer.
#' @param graph_slot_name Graph slot name inside the stats layer (auto-detected when NULL).
#' @param cluster_vec Cluster label source; either a vector of labels or a meta.data column name.
#' @param IDelta_col_name Optional Iδ column used to personalize PageRank.
#' @param IDelta_invert Whether to invert Iδ personalization.
#' @param damping PageRank damping factor.
#' @param weight_low_cut Edge-weight floor; edges below are treated as absent.
#' @param cluster_palette Palette used for clusters.
#' @param node_size Node point size.
#' @param edge_width Edge width multiplier.
#' @param label_size Label size.
#' @param seed Random seed for layout reproducibility.
#' @param length_scale Global edge-length multiplier.
#' @param max.overlaps Maximum label overlaps for `ggrepel`.
#' @param hub_border_col Hub node border color.
#' @param hub_border_size Hub node border line width.
#' @param title Optional plot title.
#' @param k_top Number of high-weight inter-cluster edges to retain in addition to the MST.
#' @param tree_mode Tree layout mode (`"radial"` or `"rooted"`).
#' @param n_runs Number of stochastic repeats.
#' @param noise_sd Standard deviation of multiplicative Gaussian noise applied to edge weights before MST selection.
#' @return Invisibly returns a list containing union graph and plot object.
#' @examples
#' \dontrun{
#' out <- plotDendroNetworkMulti(scope_obj, grid_name = "grid30", cluster_vec = "modL0.15", graph_slot_name = "g_consensus", n_runs = 20)
#' }
#' @seealso `plotDendroNetwork()`
#' @export
plotDendroNetworkMulti <- function(
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
    parent <- "plotDendroNetworkMulti"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve inputs and run config", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggraph", reason = paste0("runs=", n_runs), verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "aggregate networks and render", verbose)
    step02$enter()
    plot_obj <- .plot_dendro_network_multi_impl(
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
    step02$done()
    plot_obj
}

#' Plot dendrogram network with branch/subcluster annotations.
#'
#' @description
#' Extends `plotDendroNetwork()` by optionally splitting branches into
#' subcommunities and annotating the network.
#' @param scope_obj A `scope_object` with a LeeStats layer containing a graph.
#' @param grid_name Grid layer name (auto-selected when only one layer exists).
#' @param lee_stats_layer Lee statistics layer name (e.g. `"LeeStats_Xz"`).
#' @param gene_subset Optional gene subset (character vector or selector supported by internal helpers).
#' @param L_min Minimum Lee's L threshold for keeping edges.
#' @param drop_isolated Whether to drop isolated nodes.
#' @param use_consensus_graph Whether to use the consensus graph stored in the stats layer.
#' @param graph_slot_name Graph slot name inside the stats layer.
#' @param cluster_vec Cluster label source; either a vector of labels or a meta.data column name.
#' @param IDelta_col_name Optional Iδ column used to personalize PageRank.
#' @param damping PageRank damping factor.
#' @param weight_low_cut Edge-weight floor; edges below are treated as absent.
#' @param cluster_palette Palette used for clusters.
#' @param vertex_size Backwards-compatible alias controlling node size.
#' @param base_edge_mult Backwards-compatible alias controlling edge width.
#' @param label_cex Backwards-compatible alias controlling label size.
#' @param seed Random seed for layout reproducibility.
#' @param hub_factor Degree multiplier defining hub nodes.
#' @param length_scale Global edge-length multiplier.
#' @param max.overlaps Maximum label overlaps for `ggrepel`.
#' @param hub_border_col Hub node border color.
#' @param hub_border_size Hub node border line width.
#' @param show_sign Whether to distinguish negative edges by linetype.
#' @param neg_linetype Linetype used for negative edges.
#' @param neg_legend_lab Legend label for negative edges.
#' @param pos_legend_lab Legend label for positive edges.
#' @param show_qc_caption Whether to annotate QC caption when available.
#' @param title Optional plot title.
#' @param k_top Number of high-weight inter-cluster edges to retain in addition to the MST.
#' @param tree_mode Layout mode (`"rooted"`, `"radial"`, `"forest"`).
#' @param enable_subbranch Whether to enable branch/subcluster annotation.
#' @param cluster_id Optional cluster to focus on when sub-branching.
#' @param include_root Whether to include the root cluster in sub-branching.
#' @param max_subclusters Maximum number of subclusters to annotate.
#' @param fallback_community Whether to fall back to community detection when splitting fails.
#' @param min_sub_size Minimum size for subclusters.
#' @param community_method Community detection method (`"louvain"`, `"leiden"`).
#' @param subbranch_palette Palette used for subbranches.
#' @param downstream_min_size Optional minimum downstream size filter.
#' @param force_split Whether to attempt relaxed acceptance when splitting fails.
#' @param main_fraction_cap Maximum community fraction threshold.
#' @param core_periph Whether to allow core-periphery splitting.
#' @param core_degree_quantile Quantile threshold for defining the core.
#' @param core_min_fraction Minimum core fraction.
#' @param degree_gini_threshold Degree Gini threshold triggering core-periphery splitting.
#' @param verbose Emit progress messages when TRUE.
#' @return A list containing plot and intermediate network objects.
#' @examples
#' \dontrun{
#' out <- plotDendroNetworkWithBranches(scope_obj, grid_name = "grid30", cluster_vec = "modL0.15", graph_slot_name = "g_consensus")
#' }
#' @seealso `plotDendroNetwork()`, `getDendroWalkPaths()`
#' @export
plotDendroNetworkWithBranches <- function(
    scope_obj,
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
    force_split = TRUE,
    main_fraction_cap = 0.9,
    core_periph = TRUE,
    core_degree_quantile = 0.75,
    core_min_fraction = 0.05,
    degree_gini_threshold = 0.35,
    verbose = TRUE) {
    parent <- "plotDendroNetworkWithBranches"
    verbose <- if (is.null(verbose)) getOption("geneSCOPE.verbose", TRUE) else verbose
    step01 <- .log_step(parent, "S01", "build base network and subcluster plan", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggraph", reason = "branch_plot", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render branch plot", verbose)
    step02$enter()
    plot_obj <- .plot_dendro_network_with_branches_impl(
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
        tree_mode = tree_mode,
        enable_subbranch = enable_subbranch,
        cluster_id = cluster_id,
        include_root = include_root,
        max_subclusters = max_subclusters,
        fallback_community = fallback_community,
        min_sub_size = min_sub_size,
        community_method = community_method,
        subbranch_palette = subbranch_palette,
        downstream_min_size = downstream_min_size,
        force_split = force_split,
        main_fraction_cap = main_fraction_cap,
        core_periph = core_periph,
        core_degree_quantile = core_degree_quantile,
        core_min_fraction = core_min_fraction,
        degree_gini_threshold = degree_gini_threshold,
        verbose = verbose
    )
    step02$done()
    plot_obj
}

#' Plot per-grid density layers.
#'
#' @description
#' Renders one or two density layers stored under `scope_obj@density` for a given
#' grid layer.
#' @param overlay_image Whether to render an image raster as an additional overlay layer.
#' @param image_path Optional explicit image path used when `overlay_image = TRUE`.
#' @param image_alpha Alpha used when drawing the image overlay.
#' @param image_choice Which histology image to prefer (`"auto"`, `"hires"`, `"lowres"`).
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
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' p <- plotDensity(scope_obj, grid_name = "grid30", density1_name = "CEACAM5")
#' print(p)
#' }
#' @seealso `computeDensity()`
#' @export
plotDensity <- function(
    scope_obj,
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
    parent <- "plotDensity"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    bar_len_user <- if (missing(bar_len)) NULL else bar_len
    step01 <- .log_step(parent, "S01", "resolve density layers and overlays", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggplot2", reason = "density_plot", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render density plot", verbose)
    step02$enter()
    plot_obj <- .plot_density_impl(
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
        bar_len = bar_len_user,
        bar_offset = bar_offset,
        arrow_pt = arrow_pt,
        scale_legend_colour = scale_legend_colour,
        max.cutoff1 = max.cutoff1,
        max.cutoff2 = max.cutoff2,
        legend_digits = legend_digits
    )
    step02$done()
    plot_obj
}

#' Plot single-cell centroid expression overlays.
#'
#' @description
#' Overlays two genes' expression over single-cell centroids, with optional
#' segmentation overlays.
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
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' p <- plotDensityCentroids(scope_obj, gene1_name = "CEACAM5", gene2_name = "ACTA2")
#' print(p)
#' }
#' @seealso `addSingleCells()`
#' @export
plotDensityCentroids <- function(
    scope_obj,
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
    parent <- "plotDensityCentroids"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve inputs and overlays", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggplot2", reason = "centroid_plot", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render centroid plot", verbose)
    step02$enter()
    plot_obj <- .plot_density_centroids_impl(
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
    step02$done()
    plot_obj
}

#' Plot the spatial boundary of grid tiles.
#'
#' @description
#' Draws rectangles for each grid tile in a specified grid layer.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param colour Parameter value.
#' @param linewidth Parameter value.
#' @param panel_bg Parameter value.
#' @param base_size Parameter value.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' p <- plotGridBoundary(scope_obj, grid_name = "grid30")
#' print(p)
#' }
#' @seealso `addScopeHistology()`, `plotDensity()`
#' @export
plotGridBoundary <- function(
    scope_obj,
    grid_name,
    colour = "black",
    linewidth = 0.2,
    panel_bg = "#C0C0C0",
    base_size = 10) {
    parent <- "plotGridBoundary"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve grid boundary", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggplot2", reason = "boundary_plot", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render boundary plot", verbose)
    step02$enter()
    plot_obj <- .plot_grid_boundary_impl(
        scope_obj = scope_obj,
        grid_name = grid_name,
        colour = colour,
        linewidth = linewidth,
        panel_bg = panel_bg,
        base_size = base_size
    )
    step02$done()
    plot_obj
}

#' Plot Moran's Iδ per cluster.
#'
#' @description
#' Creates a faceted scatter-line plot of Iδ values for top genes in each cluster.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param cluster_col Parameter value.
#' @param top_n Parameter value.
#' @param min_genes Numeric threshold.
#' @param nrow Parameter value.
#' @param point_size Parameter value.
#' @param line_size Parameter value.
#' @param label_size Parameter value.
#' @param subCluster Parameter value.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' p <- plotIDelta(scope_obj, grid_name = "grid30", cluster_col = "modL0.15")
#' print(p)
#' }
#' @seealso `computeIDelta()`
#' @export
plotIDelta <- function(
    scope_obj,
    grid_name,
    cluster_col,
    top_n = NULL,
    min_genes = 1,
    nrow = 1,
    point_size = 3,
    line_size = 0.5,
    label_size = 2,
    subCluster = NULL) {
    parent <- "plotIDelta"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve inputs and filters", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggplot2", reason = "idelta_plot", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render I-delta plot", verbose)
    step02$enter()
    plot_obj <- .plot_idelta_impl(
        scope_obj = scope_obj,
        grid_name = grid_name,
        cluster_col = cluster_col,
        top_n = top_n,
        min_genes = min_genes,
        nrow = nrow,
        point_size = point_size,
        line_size = line_size,
        label_size = label_size,
        subCluster = subCluster
    )
    step02$done()
    plot_obj
}

#' Plot the distribution of Lee's L values.
#'
#' @description
#' Draws a histogram of Lee's L across gene pairs for a given grid layer.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param bins Parameter value.
#' @param xlim Parameter value.
#' @param title Parameter value.
#' @param use_abs Logical flag.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' p <- plotLDistribution(scope_obj, grid_name = "grid30")
#' print(p)
#' }
#' @seealso `computeL()`
#' @export
plotLDistribution <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = NULL,
    bins = 30,
    xlim = NULL,
    title = NULL,
    use_abs = FALSE) {
    parent <- "plotLDistribution"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve Lee's L matrix and bins", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggplot2", reason = "histogram", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render distribution plot", verbose)
    step02$enter()
    plot_obj <- .plot_l_distribution(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        bins = bins,
        xlim = xlim,
        title = title,
        use_abs = use_abs
    )
    step02$done()
    plot_obj
}

#' Plot Lee's L vs fold-change scatter.
#'
#' @description
#' Builds a scatter plot of Lee's L vs gene total count fold change for all gene pairs.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param title Parameter value.
#' @return A `ggplot2` object.
#' @examples
#' \dontrun{
#' p <- plotLScatter(scope_obj, grid_name = "grid30")
#' print(p)
#' }
#' @seealso `plotLDistribution()`, `plotLvsR()`
#' @export
plotLScatter <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = NULL,
    title = NULL) {
    parent <- "plotLScatter"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve inputs and gene pairs", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggplot2", reason = "scatter", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render scatter plot", verbose)
    step02$enter()
    plot_obj <- .plot_l_scatter(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        title = title
    )
    step02$done()
    plot_obj
}

#' Plot Lee's L vs Pearson correlation.
#'
#' @description
#' Builds a scatter plot comparing Lee's L with Pearson correlation for gene pairs.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param pear_level Parameter value.
#' @param lee_stats_layer Layer name.
#' @param delta_top_n Parameter value.
#' @param flip Parameter value.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' p <- plotLvsR(scope_obj, grid_name = "grid30", pear_level = "grid")
#' p2 <- p + 
#'     geom_ribbon(data = scope_obj@stats$grid30$LeeStats_Xz$LR_curve_grid30,
#'         aes(x = Pear, ymin = lo95, ymax = hi95),
#'         inherit.aes = FALSE, fill = "orange", alpha = 0.25) +
#'     geom_line(data = scope_obj@stats$grid30$LeeStats_Xz$LR_curve_grid30,
#'         aes(x = Pear, y = fit),
#'         inherit.aes = FALSE, colour = "firebrick", linewidth = 0.8)
#' print(p)
#' }
#' @seealso `computeLvsRCurve()`, `getTopLvsR()`
#' @export
plotLvsR <- function(
    scope_obj,
    grid_name,
    pear_level = c("cell", "grid"),
    lee_stats_layer = "LeeStats_Xz",
    delta_top_n = 10,
    flip = TRUE) {
    parent <- "plotLvsR"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve L/r matrices and labels", verbose)
    step01$enter()
    .log_backend(parent, "S01", "backend", "ggplot2", reason = "scatter", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render L vs R plot", verbose)
    step02$enter()
    plot_obj <- .plot_l_vs_r(
        scope_obj = scope_obj,
        grid_name = grid_name,
        pear_level = pear_level,
        lee_stats_layer = lee_stats_layer,
        delta_top_n = delta_top_n,
        flip = flip
    )
    step02$done()
    plot_obj
}

#' Plot the gene network graph.
#'
#' @description
#' Visualizes a gene network derived from Lee's L / consensus graphs, optionally
#' thresholded by similarity and FDR.
#' @param scope_obj A `scope_object` containing a LeeStats layer and a graph.
#' @param grid_name Grid layer name (auto-selected when only one layer exists).
#' @param lee_stats_layer Lee statistics layer name (e.g. `"LeeStats_Xz"`).
#' @param gene_subset Optional gene subset (character vector or selector supported by internal helpers).
#' @param L_min Minimum Lee's L threshold for positive edges.
#' @param L_min_neg Optional threshold for negative edges.
#' @param p_cut Optional p-value cutoff (when available).
#' @param use_FDR Whether to use an FDR matrix for edge filtering.
#' @param FDR_max Maximum allowed FDR.
#' @param pct_min Gene prevalence filter used by internal helper logic.
#' @param CI95_filter Deprecated flag retained for compatibility.
#' @param curve_layer Curve layer name used by CI filtering.
#' @param CI_rule Confidence-interval filtering rule (`remove_within`, `remove_outside`).
#' @param drop_isolated Whether to drop isolated nodes from the plotted graph.
#' @param weight_abs Whether to use absolute edge weights for sizing.
#' @param use_consensus_graph Whether to use an already-filtered consensus graph (skips re-filtering).
#' @param graph_slot_name Graph slot name inside the stats layer (auto-detected when NULL).
#' @param fdr_source Which FDR slot to use (e.g. `"FDR"`, `"FDR_storey"`, `"FDR_beta"`).
#' @param cluster_vec Cluster label source; either a vector of labels or a meta.data column name.
#' @param cluster_palette Palette used for cluster colors.
#' @param node_size Node point size.
#' @param edge_width Edge width multiplier.
#' @param label_size Label size.
#' @param vertex_size Deprecated alias of `node_size`.
#' @param base_edge_mult Deprecated alias of `edge_width`.
#' @param label_cex Deprecated alias of `label_size`.
#' @param layout_niter Iterations for the layout algorithm.
#' @param seed Random seed for layout reproducibility.
#' @param hub_factor Degree multiplier defining hub nodes.
#' @param length_scale Global edge-length multiplier.
#' @param max.overlaps Maximum label overlaps for `ggrepel`.
#' @param L_max Maximum Lee's L used for edge scaling.
#' @param hub_border_col Hub node border color.
#' @param hub_border_size Hub node border line width.
#' @param show_sign Whether to distinguish negative edges by linetype.
#' @param neg_linetype Linetype used for negative edges.
#' @param neg_legend_lab Legend label for negative edges.
#' @param pos_legend_lab Legend label for positive edges.
#' @param title Optional plot title.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' p <- plotNetwork(scope_obj, grid_name = "grid30", L_min = 0.15, cluster_vec = "modL0.15", graph_slot_name = "g_consensus")
#' print(p)
#' }
#' @seealso `computeL()`, `clusterGenes()`, `plotDendroNetwork()`
#' @export
plotNetwork <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL,
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
    use_consensus_graph = TRUE,
    graph_slot_name = NULL,
    fdr_source = c("FDR", "FDR_storey", "FDR_beta", "FDR_mid", "FDR_uniform", "FDR_disc"),
    cluster_vec = NULL,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
        "#A65628", "#984EA3", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
    ),
    node_size = 4,
    edge_width = 3,
    label_size = 4,
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
    parent <- "plotNetwork"
    verbose <- getOption("geneSCOPE.verbose", TRUE)
    step01 <- .log_step(parent, "S01", "resolve graph and filters", verbose)
    step01$enter()
    graph_source <- if (isTRUE(use_consensus_graph)) "consensus" else "reconstructed"
    .log_backend(parent, "S01", "graph_source", graph_source, reason = "use_consensus_graph", verbose = verbose)
    step01$done()

    step02 <- .log_step(parent, "S02", "render network plot", verbose)
    step02$enter()
    plot_obj <- .plot_network_impl(
        scope_obj = scope_obj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        gene_subset = gene_subset,
        L_min = L_min,
        L_min_neg = L_min_neg,
        p_cut = p_cut,
        use_FDR = use_FDR,
        FDR_max = FDR_max,
        pct_min = pct_min,
        CI95_filter = CI95_filter,
        curve_layer = curve_layer,
        CI_rule = CI_rule,
        drop_isolated = drop_isolated,
        weight_abs = weight_abs,
        use_consensus_graph = use_consensus_graph,
        graph_slot_name = graph_slot_name,
        fdr_source = fdr_source,
        cluster_vec = cluster_vec,
        cluster_palette = cluster_palette,
        node_size = node_size,
        edge_width = edge_width,
        label_size = label_size,
        vertex_size = vertex_size,
        base_edge_mult = base_edge_mult,
        label_cex = label_cex,
        layout_niter = layout_niter,
        seed = seed,
        hub_factor = hub_factor,
        length_scale = length_scale,
        max.overlaps = max.overlaps,
        L_max = L_max,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        show_sign = show_sign,
        neg_linetype = neg_linetype,
        neg_legend_lab = neg_legend_lab,
        pos_legend_lab = pos_legend_lab,
        title = title
    )
    step02$done()
    plot_obj
}

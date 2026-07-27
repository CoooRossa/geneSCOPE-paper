
#' Plot cluster assignments across methods
#' @description
#' Internal helper for `.plot_cluster_comparison`.
#' Creates a scatter-style comparison plot of cluster calls from multiple
#' methods for each gene.
#' @param scope_obj A `scope_object` with `meta.data` cluster columns.
#' @param method_cols Character vector of column names to compare.
#' @param method_labels Optional labels for the methods (same length as `method_cols`).
#' @param point_size Point size for plotted markers.
#' @param palette Palette function or vector passed to `scale_fill_manual`.
#' @return A `ggplot` object.
#' @keywords internal
.plot_cluster_comparison <- function(scope_obj,
                                  method_cols,
                                  method_labels = method_cols,
                                  point_size = 3,
                                  palette = RColorBrewer::brewer.pal) {
    stopifnot(length(method_cols) == length(method_labels))

    df <- scope_obj@meta.data |>
        rownames_to_column("gene") |>
        select(gene, tidyselect::all_of(method_cols))

    # keep genes with at least one non-NA assignment
    df <- df |>
        filter(if_any(-gene, ~ !is.na(.x)))

    # order genes by progressively less stringent columns
    df <- df |>
        arrange(across(all_of(method_cols), ~ desc(.x)))

    long_df <- df |>
        pivot_longer(
            cols = tidyselect::all_of(method_cols),
            names_to = "method",
            values_to = "cluster"
        ) |>
        drop_na(cluster) |>
        mutate(
            method  = factor(method, levels = method_cols, labels = method_labels),
            gene    = factor(gene, levels = df$gene),
            cluster = factor(cluster)
        )

    n_col <- nlevels(long_df$cluster)
    pal <- if (is.function(palette)) {
        if (n_col <= 12) palette(n_col, "Set3") else scales::hue_pal()(n_col)
    } else {
        if (length(palette) < n_col) {
            colorRampPalette(palette)(n_col)
        } else {
            palette[seq_len(n_col)]
        }
    }

    p <- ggplot(
        long_df,
        aes(x = gene, y = method, fill = cluster)
    ) +
        geom_point(shape = 21, size = point_size, colour = "black") +
        scale_fill_manual(values = pal) +
        labs(x = NULL, y = "Group", fill = "Cluster") +
        theme_minimal(base_size = 8) +
        theme(
            axis.text.x = element_text(
                angle = 90, vjust = .5, hjust = 1,
                size = 5
            ),
            panel.title = element_text(size = 10, colour = "black"),
            panel.background = element_rect(fill = "#c0c0c0", colour = NA),
            panel.border = element_rect(colour = "black", fill = NA, size = .5),
            axis.line = element_line(colour = "black", size = .3),
            axis.ticks = element_line(colour = "black", size = .3),
            axis.text = element_text(size = 8, colour = "black"),
            axis.title = element_text(size = 9, colour = "black"),
            legend.position = "none"
        )

    p
}

#' Plot Dendro
#' @description
#' Internal helper for `.plot_dendro`.
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
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_dendro <- function(
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
    .plot_dendro_core(
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

#' Plot L Distribution
#' @description
#' Internal helper for `.plot_l_distribution`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param bins Parameter value.
#' @param xlim Parameter value.
#' @param title Parameter value.
#' @param use_abs Logical flag.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_l_distribution <- function(scope_obj,
                              grid_name = NULL,
                              lee_stats_layer = NULL,
                              bins = 30,
                              xlim = NULL,
                              title = NULL,
                              use_abs = FALSE) {
  .plot_l_distribution_core(
    scope_obj = scope_obj,
    grid_name = grid_name,
    lee_stats_layer = lee_stats_layer,
    bins = bins,
    xlim = xlim,
    title = title,
    use_abs = use_abs
  )
}

#' Plot L Scatter
#' @description
#' Internal helper for `.plot_l_scatter`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param title Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_l_scatter <- function(scope_obj,
                         grid_name = NULL,
                         lee_stats_layer = NULL,
                         title = NULL) {
  .plot_l_scatter_core(
    scope_obj = scope_obj,
    grid_name = grid_name,
    lee_stats_layer = lee_stats_layer,
    title = title
  )
}

#' Plot L Vs R
#' @description
#' Internal helper for `.plot_l_vs_r`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param pear_level Parameter value.
#' @param lee_stats_layer Layer name.
#' @param delta_top_n Parameter value.
#' @param flip Parameter value.
#' @return A `ggplot` object (or a list of plot components) used internally.
#' @keywords internal
.plot_l_vs_r <- function(scope_obj,
                     grid_name,
                     pear_level = c("cell", "grid"),
                     lee_stats_layer = "LeeStats_Xz",
                     delta_top_n = 10,
                     flip = TRUE) {
  .plot_lvs_r_core(
    scope_obj = scope_obj,
    grid_name = grid_name,
    pear_level = pear_level,
    lee_stats_layer = lee_stats_layer,
    delta_top_n = delta_top_n,
    flip = flip
  )
}

#' Plot a gene network using Lee's L statistics.
#' @description
#' Internal helper for `.plot_network`.
#' @param grid_name Character. Grid sub-layer to plot. Defaults to the
#'                        first layer.
#' @param lee_stats_layer Name of the Lee's stats layer (default \code{"LeeStats_Xz"}).
#' @param gene_subset Optional character vector of genes to keep.
#' @param L_min Minimum positive / negative |L| thresholds.
#' @param L_min_neg Minimum positive / negative |L| thresholds.
#' @param p_cut p-value or FDR cut-offs when the corresponding
#' @param FDR_max p-value or FDR cut-offs when the corresponding
#'                        matrices are available.
#' @param pct_min Quantile string such as \code{"q80"} passed to the
#'                        internal quantile filter.
#' @param CI95_filter Logical. Remove edges according to the 95 % LR-curve
#'                        band using \code{CI_rule}.
#' @param curve_layer Name of the LR-curve layer for CI filtering.
#' @param CI_rule \code{"remove_within"} or \code{"remove_outside"}.
#' @param drop_isolated Logical. Discard isolated nodes.
#' @param weight_abs Use absolute L for edge weight (default \code{TRUE}).
#' @param use_consensus_graph Logical. Use pre-computed consensus graph stored
#'                        in \code{graph_slot_name}.
#' @param graph_slot_name Name under which the consensus graph is stored.
#' @param cluster_vec Either a named vector of cluster assignments or the
#'                        name of a column in \code{meta.data}.
#' @param cluster_palette Character vector of colours (named or unnamed).
#' @param node_size Node glyph size (replaces deprecated `vertex_size`).
#' @param edge_width Maximum edge width scaling (replaces deprecated `base_edge_mult`).
#' @param label_size Text size for node labels (replaces deprecated `label_cex`).
#' @param vertex_size Deprecated aliases kept for backwards compatibility.
#' @param base_edge_mult Deprecated aliases kept for backwards compatibility.
#' @param label_cex Deprecated aliases kept for backwards compatibility.
#' @param layout_niter Iterations for the Fruchterman-Reingold layout.
#' @param seed Random seed for layout reproducibility.
#' @param hub_factor Degree multiplier defining "hub" nodes.
#' @param length_scale Global edge-length multiplier.
#' @param show_sign Draw negative edges with distinct linetype.
#' @param neg_linetype Linetype for negative edges.
#' @param title Optional plot title.
#' @param scope_obj A `scope_object`.
#' @param use_FDR Logical flag.
#' @param fdr_source Parameter value.
#' @param max.overlaps Parameter value.
#' @param L_max Numeric threshold.
#' @param hub_border_col Parameter value.
#' @param hub_border_size Parameter value.
#' @param neg_legend_lab Parameter value.
#' @param pos_legend_lab Parameter value.
#' @return A \code{ggplot} object.
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_point geom_node_text scale_edge_width scale_edge_colour_identity scale_edge_linetype_manual
#' @importFrom ggplot2 aes labs guides guide_legend coord_fixed theme_void theme element_text unit margin
#' @importFrom igraph V E induced_subgraph graph_from_adjacency_matrix get.edgelist subgraph.edges degree
#' @importFrom grDevices col2rgb rgb rainbow colorRampPalette
#' @importFrom dplyr filter
#' @importFrom stats approxfun quantile hclust as.dendrogram
#' @importFrom Matrix drop0
#' @keywords internal
.plot_network <- function(
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
        stop("Less than two genes remain after sub‑setting.")
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

        ## (v) p‑value & FDR
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

    ## ① Derive clu vector from cluster_vec or meta.data
    clu <- setNames(rep(NA_character_, length(Vnames)), Vnames)
    is_factor_input <- FALSE
    factor_levels <- NULL

    ## helper – extract column from meta.data
    get_meta_cluster <- function(col) {
        if (is.null(scope_obj@meta.data) || !(col %in% colnames(scope_obj@meta.data))) {
            stop("Column ‘", col, "’ not found in scope_obj@meta.data.")
        }
        tmp <- scope_obj@meta.data[[col]]
        names(tmp) <- rownames(scope_obj@meta.data)
        tmp[Vnames]
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

    ## Step 3: edge colours (sign / intensity gradient)
    e_idx <- igraph::as_edgelist(g, names = FALSE)
    w_norm <- igraph::E(g)$weight / max(igraph::E(g)$weight)
    cluster_ord <- setNames(seq_along(uniq_clu), uniq_clu)
    edge_cols <- character(igraph::ecount(g))
    for (i in seq_along(edge_cols)) {
        v1 <- Vnames[e_idx[i, 1]]
        v2 <- Vnames[e_idx[i, 2]]
        if (show_sign && igraph::E(g)$sign[i] == "neg") {
            edge_cols[i] <- "gray40"
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
            edge_cols[i] <- rgb(
                1 - tval + tval * rgb_ref[1],
                1 - tval + tval * rgb_ref[2],
                1 - tval + tval * rgb_ref[3]
            )
        }
    }
    igraph::E(g)$edge_col <- edge_cols
    igraph::E(g)$linetype <- if (show_sign) igraph::E(g)$sign else "solid"

    ## ===== 6. ggraph rendering =====
    library(ggraph)
    library(ggplot2)
    library(dplyr)
    set.seed(seed)
    lay <- create_layout(g, layout = "fr", niter = layout_niter)
    lay$basecol <- basecol[lay$name]
    lay$deg <- deg_vec[lay$name]
    lay$hub <- lay$deg > hub_factor * median(lay$deg)

    ## ---- edge width / colour / linetype layers ----
    p <- ggraph(lay) +
        geom_edge_link(
            aes(
                width = weight, colour = edge_col,
                linetype = linetype
            ),
            lineend = "round",
            show.legend = FALSE
        ) +
        scale_edge_width(name = "|L|", range = c(0.3, edge_width), guide = "none") +
        scale_edge_colour_identity(guide = "none") + {
            if (show_sign) {
                scale_edge_linetype_manual(
                    name   = "Direction",
                    values = c(pos = "solid", neg = neg_linetype),
                    breaks = c("pos", "neg"),
                    labels = c(pos = pos_legend_lab, neg = neg_legend_lab),
                    guide  = "none"
                )
            }
        }

    ## ---- node points / text ----
    p <- p +
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
        labs(title = title)

    ## Optional legend order adjustment
    p <- p +
        guides(
            fill = guide_legend(order = 1),
            size = "none",
            edge_width = "none",
            linetype = "none",
            colour = "none"
        )

    return(p)
}

#' Dendrogram of the .plot_dendro_network skeleton (plotRWDendrogram8-style)
#' @description
#'   Builds the same MST/forest skeleton as in `.plot_dendro_network` (cluster-internal
#'   MST + cluster-level MST after optional PageRank reweighting) from the
#'   specified consensus graph, then renders a dendrogram in the identical
#'   visual style as `plotRWDendrogram8` with optional dual tip-dot rows,
#'   adjustable label offsets, and OLO leaf-order optimisation.
#' @param graph_slot Slot in the LeeStats layer that stores the base graph
#'                   (e.g. the consensus graph). Defaults to `"g_consensus"`.
#' @param title Custom title for the rendered dendrogram. Defaults to
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
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
#'   `"Graph-weighted Dendrogram"`.
#' @return Invisibly returns a list containing `dend`, `hclust`, `dist`,
#'   `genes`, `cluster_map`, `tree_graph` (MST skeleton), `graph` (full graph
#'   on the selected genes), and `PageRank` if applied.
#' @importFrom igraph V E as_data_frame induced_subgraph distances ecount vcount edge_attr_names page_rank ends
#' @importFrom stats hclust as.dendrogram quantile runmed
#' @importFrom graphics plot points text par
#' @keywords internal
.plot_dendro <- function(
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
    if (!exists("tip_label_indent", inherits = FALSE) || is.null(tip_label_indent)) tip_label_indent <- 0.01
    if (!exists("legend_inline", inherits = FALSE) || is.null(legend_inline)) legend_inline <- FALSE
    if (!is.null(distance_smooth_power) && distance_smooth_power <= 0) {
        stop("distance_smooth_power must be > 0 or NULL")
    }

    ## ---------- 0. Resolve grid layer & LeeStats ----------
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    if (is.null(grid_name)) {
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
    }

    ## LeeStats object
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

    ## ---------- 2. Iδ-PageRank reweighting (same as .plot_dendro_network) ----------
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

    ## ---------- 4. Distance matrix (tree-based to merge within clusters first) ----------
    genes <- igraph::V(g_tree)$name
    # keep full-graph view for return, regardless of distance_on
    g_full <- igraph::induced_subgraph(g, vids = genes)
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

#' Plot Dendrogram-style Network Layout
#' @description
#'   Creates a tree-like visualization of a gene network using minimum spanning
#'   tree algorithms within and between clusters, optionally weighted by
#'   personalized PageRank scores derived from Morisita's Iδ values.
#' @param IDelta_col_name Character. Column name in \code{meta.data} containing
#'                        Morisita's Iδ values for PageRank weighting. If \code{NULL},
#'                        uniform PageRank weighting is used.
#'                        personalisation (highest Iδ receives lowest weight).
#' @param damping Numeric. Damping factor for PageRank algorithm (default 0.85).
#' @param weight_low_cut Numeric. Minimum edge weight threshold after PageRank
#'                        weighting (default 0).
#' @param k_top Integer. Maximum number of high-weight non-tree edges
#'                        to retain between clusters (default 1).
#' @param tree_mode Character. Tree layout style: \code{"rooted"},
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param lee_stats_layer Layer name.
#' @param gene_subset Parameter value.
#' @param L_min Numeric threshold.
#' @param use_consensus_graph Logical flag.
#' @param graph_slot_name Slot name.
#' @param cluster_vec Parameter value.
#' @param IDelta_invert Parameter value.
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
#'                        \code{"radial"}, or \code{"forest"} (default "rooted").
#' @return Invisibly returns a list with components:
#'   \describe{
#'     \item{\code{graph}}{The processed \code{igraph} object.}
#'     \item{\code{pagerank}}{PageRank scores (personalised when \code{IDelta_col_name} is provided).}
#'     \item{\code{cross_edges}}{Data frame of inter-cluster edges in the final graph.}
#'   }
#' @importFrom igraph mst components subgraph.edges delete_edges delete_vertices degree E V as_data_frame ends edge_attr page_rank induced_subgraph ecount vcount
#' @importFrom stats quantile median
#' @keywords internal
.plot_dendro_network <- function(
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
    ## ---------- Iδ-PageRank options ----------
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
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else {
        grid_name
    }
    leeStat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
        scope_obj@stats[[grid_name]][[lee_stats_layer]]
    } else {
        g_layer[[lee_stats_layer]]
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

    ## ========= 2. Iδ-PageRank reweighting / random walk =========
    ## Always run a random-walk weighting step. When Iδ is provided, use it to
    ## personalise PageRank (with optional reversal); otherwise, use uniform
    ## personalization for a simple random walk.
    delta <- NULL
    if (!is.null(IDelta_col_name)) {
        delta <- scope_obj@meta.data[V(g)$name, IDelta_col_name, drop = TRUE]
        delta[is.na(delta)] <- median(delta, na.rm = TRUE)
        q <- quantile(delta, c(.1, .9))
        delta <- pmax(pmin(delta, q[2]), q[1]) # align with dendroRW behaviour
        if (isTRUE(IDelta_invert)) {
            # flip preference: higher Iδ gets lower personalization weight
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
    ## —— 5.1 Intra-cluster MST ——
    all_edges <- igraph::as_data_frame(g, "edges")
    all_edges$key <- with(
        all_edges,
        ifelse(from < to, paste(from, to, sep = "|"),
            paste(to, from, sep = "|")
        )
    )
    keep_key <- character(0)

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
        keep_key <- c(keep_key, ks$key)
    }

    ## —— 5.2 Inter-cluster MST ——
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
            ## Obtain MST for each connected component
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
                # retain the original edge with maximal weight
                for (k in ks$key) {
                    cand <- inter[inter$pair == k, ]
                    cand <- cand[order(cand$weight, decreasing = TRUE), ]
                    keep_key <- c(keep_key, cand$key[1])
                }
            }
        }
    }

    ## —— 5.3 Filter edges based on keep_key ——
    keep_eid <- which(all_edges$key %in% unique(keep_key))
    g <- igraph::subgraph.edges(g, keep_eid, delete.vertices = TRUE)
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))

    ## ============ Tree-layout helper steps ==============
    ## Retain original edge IDs for later mapping
    igraph::E(g)$eid <- seq_len(igraph::ecount(g))

    ## ---- 6.1  Intra-cluster MST ----
    keep_eid <- integer(0)
    for (cl in unique(clu)) {
        vsub <- Vnames[clu == cl]
        if (length(vsub) < 2) next
        g_sub <- igraph::induced_subgraph(g, vsub)
        if (igraph::ecount(g_sub) > 0) {
            mst_sub <- igraph::mst(g_sub, weights = igraph::E(g_sub)$weight)
            keep_eid <- c(keep_eid, igraph::edge_attr(mst_sub, "eid"))
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
        inter <- ed_tab[ed_tab$cl1 != ed_tab$cl2, ] # consider only inter-cluster edges

        if (nrow(inter)) {
            ## ---- 6.2a  Build the cluster-level graph and take its MST ----
            inter$pair <- ifelse(inter$cl1 < inter$cl2,
                paste(inter$cl1, inter$cl2, sep = "|"),
                paste(inter$cl2, inter$cl1, sep = "|")
            )
            agg <- aggregate(weight ~ pair + cl1 + cl2, data = inter, max)
            g_clu <- igraph::graph_from_data_frame(
                agg[, c("cl1", "cl2", "weight")],
                directed = FALSE, vertices = unique(clu)
            )

            ## For disconnected cases, take an MST per component
            mst_clu <- igraph::mst(
                g_clu,
                weights = 1 / (igraph::E(g_clu)$weight + 1e-9)
            )

            # Extract endpoint pairs from that MST‐forest:
            ep <- igraph::ends(mst_clu, igraph::E(mst_clu))
            keep_pairs <- ifelse(
                ep[, 1] < ep[, 2],
                paste(ep[, 1], ep[, 2], sep = "|"),
                paste(ep[, 2], ep[, 1], sep = "|")
            )
            inter_MST <- inter[inter$pair %in% keep_pairs, ]

            ## ---- 6.2b  Identify high-weight off-tree edges ----
            extra_edges <- inter[!(inter$pair %in% keep_pairs), ]
            if (nrow(extra_edges)) {
                if (!is.null(w_extra_min)) {
                    extra_edges <- extra_edges[extra_edges$weight >= w_extra_min, ]
                }
                if (!is.null(k_top) && nrow(extra_edges) > k_top) {
                    extra_edges <- extra_edges[order(extra_edges$weight, decreasing = TRUE)[seq_len(k_top)], ]
                }
            }
            ## ---- 6.2c  Collect edge IDs to retain ----
            keep_eid <- c(
                keep_eid,
                inter_MST$eid,
                if (nrow(extra_edges)) extra_edges$eid else integer(0)
            )
        }
    }
    keep_eid <- unique(keep_eid)
    g <- igraph::delete_edges(g, igraph::E(g)[!eid %in% keep_eid])
    ## Remove isolated vertices that remain
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
    ## ============ Tree layout complete ==============
    ## —— 6.3  Derive inter-cluster edges for plotting ——
    # Convert the current graph into a data.frame of edges
    edf_final <- igraph::as_data_frame(g, what = "edges")
    # Keep only edges whose endpoints belong to different clusters
    cross_edges <- edf_final[clu[edf_final$from] != clu[edf_final$to], ]

    ## ===== 7. Recompute degree / colours =====
    Vnames <- igraph::V(g)$name
    deg_vec <- igraph::degree(g)

    is_factor_input <- FALSE
    factor_levels <- NULL

    ## helper – extract column from meta.data
    get_meta_cluster <- function(col) {
        if (is.null(scope_obj@meta.data) || !(col %in% colnames(scope_obj@meta.data))) {
            stop("Column ‘", col, "’ not found in scope_obj@meta.data.")
        }
        tmp <- scope_obj@meta.data[[col]]
        names(tmp) <- rownames(scope_obj@meta.data)
        tmp[Vnames]
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
    cluster_ord <- setNames(seq_along(uniq_clu), uniq_clu)
    edge_cols <- character(igraph::ecount(g))
    for (i in seq_along(edge_cols)) {
        v1 <- Vnames[e_idx[i, 1]]
        v2 <- Vnames[e_idx[i, 2]]
        {
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
    }
    igraph::E(g)$edge_col <- edge_cols
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
        root_v <- V(g)[which.max(deg_vec)] # same approach as before
        lay <- create_layout(g, layout = "tree", root = root_v)
    } else if (tree_mode == "radial") {
        lay <- create_layout(g, layout = "tree", circular = TRUE)
    }
    lay$basecol <- basecol[lay$name]
    lay$deg <- deg_vec[lay$name]
    hub_factor <- 2
    lay$hub <- lay$deg > hub_factor * median(lay$deg)

    qc_txt <- NULL

    p <- ggraph(lay) +
        geom_edge_link(aes(width = weight, colour = edge_col, linetype = linetype),
            lineend = "round",
            show.legend = FALSE
        ) +
        scale_edge_width(name = "|L|", range = c(0.3, edge_width), guide = "none") +
        scale_edge_colour_identity(guide = "none") +
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
            repel = TRUE,
            vjust = 1.4, max.overlaps = max.overlaps
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
        labs(title = title, caption = qc_txt) +
        guides(
            fill = guide_legend(order = 1),
            size = "none",
            edge_width = "none",
            linetype = "none",
            colour = "none"
        )

    invisible(list(
        graph       = g,
        pagerank    = if (exists("pr")) pr else NULL,
        cross_edges = cross_edges,
        plot        = p
    ))
}

#' Plot Dendrogram-style Network Layout (multi-run stochastic MST union)
#' @description
#'   Repeats the PageRank-weighted MST construction used in `.plot_dendro_network`
#'   with small random perturbations, unions the resulting tree backbones, and
#'   renders a graph that can include multiple plausible inter/ intra-cluster
#'   paths instead of a single deterministic tree.
#' @param n_runs Integer. Number of stochastic repeats to perform (default 10).
#' @param noise_sd Numeric. Standard deviation of multiplicative Gaussian noise
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
#'                 applied to edge weights before MST selection (default 0.01).
#'                 Set to 0 to recover deterministic behaviour.
#' @return Invisibly returns a list with components:
#'   \describe{
#'     \item{\code{graph}}{Union graph of MST edges across runs, with weights averaged over occurrences.}
#'     \item{\code{pagerank}}{PageRank scores from the last run (for reference).}
#'     \item{\code{run_edge_counts}}{Data frame of edge occurrence counts across runs.}
#'     \item{\code{plot}}{ggraph object.}
#'   }
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

        g_layer <- .select_grid_layer(scope_obj, grid_name)
        gname <- if (is.null(grid_name)) {
            names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
        } else grid_name
        leeStat <- if (!is.null(scope_obj@stats[[gname]]) &&
            !is.null(scope_obj@stats[[gname]][[lee_stats_layer]])) {
            scope_obj@stats[[gname]][[lee_stats_layer]]
        } else {
            g_layer[[lee_stats_layer]]
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
    run_list <- vector("list", n_runs)
    for (i in seq_len(n_runs)) {
        run_list[[i]] <- build_once(i)
    }

    # Union edges and aggregate weights / counts
    all_edges <- do.call(rbind, lapply(seq_along(run_list), function(i) {
        cbind(run = i, run_list[[i]]$edges)
    }))
    agg <- aggregate(weight ~ key + from + to, data = all_edges, FUN = mean)
    counts <- aggregate(weight ~ key, data = all_edges, FUN = length)
    names(counts)[2] <- "count"
    agg <- merge(agg, counts, by = "key", all.x = TRUE)

    # Build union graph
    g_union <- igraph::graph_from_data_frame(
        agg[, c("from", "to", "weight")],
        directed = FALSE
    )
    igraph::E(g_union)$count <- agg$count

    # Use clusters from first run (assumed consistent)
    clu <- run_list[[1]]$clu
    V(g_union)$clu <- clu[V(g_union)$name]

    # Colours as in .plot_dendro_network
    Vnames <- igraph::V(g_union)$name
    deg_vec <- igraph::degree(g_union)
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
    e_idx <- igraph::as_edgelist(g_union, names = FALSE)
    w_norm <- igraph::E(g_union)$weight / max(igraph::E(g_union)$weight)
    cluster_ord <- setNames(seq_along(uniq_clu), uniq_clu)
    edge_cols <- character(igraph::ecount(g_union))
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
    igraph::E(g_union)$edge_col <- edge_cols
    igraph::E(g_union)$linetype <- "solid"

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
        root_v <- V(g_union)[which.max(deg_vec)]
        lay <- create_layout(g_union, layout = "tree", root = root_v)
    } else if (tree_mode == "radial") {
        lay <- create_layout(g_union, layout = "tree", circular = TRUE)
    }
    lay$basecol <- basecol[lay$name]
    lay$deg <- deg_vec[lay$name]
    hub_factor <- 2
    lay$hub <- lay$deg > hub_factor * median(lay$deg)

    p <- ggraph(lay) +
        geom_edge_link(aes(width = weight, colour = edge_col, linetype = linetype),
            lineend = "round",
            show.legend = FALSE
        ) +
        scale_edge_width(name = "|L|", range = c(0.3, edge_width), guide = "none") +
        scale_edge_colour_identity(guide = "none") +
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
            repel = TRUE,
            vjust = 1.4, max.overlaps = max.overlaps
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
       labs(title = title) +
       guides(
            fill = guide_legend(order = 1),
            size = "none",
            edge_width = "none",
            linetype = "none",
            colour = "none"
        )

    occ_df <- data.frame(
        edge = agg$key,
        count = agg$count,
        stringsAsFactors = FALSE
    )

    invisible(list(
        graph = g_union,
        pagerank = run_list[[length(run_list)]]$pr,
        run_edge_counts = occ_df,
        plot = p
    ))
}

#' Dendrogram-based Branch Subclustering Visualization
#' @description
#' Calls existing \code{.plot_dendro_network()} to build tree/skeleton network, then performs
#' "articulation point branching" subclustering refinement on specified clusters (or all clusters);
#' if no candidates and \code{fallback_community=TRUE} is set,
#' falls back to Louvain/Leiden community splitting. Outputs new network plot with subcluster coloring.
#' @param enable_subbranch Whether to enable branch subclustering refinement (default TRUE)
#' @param cluster_id Optional character vector: clusters to subcluster split only; NULL=all
#' @param include_root Whether articulation point method includes articulation points in each branch
#' @param max_subclusters Maximum number of subclusters to retain per cluster (greedy deduplication by size)
#' @param fallback_community Whether to fall back to community detection when articulation points yield no results
#' @param min_sub_size Minimum size threshold for subclusters in community fallback mode
#' @param community_method Priority order of community detection algorithms
#' @param subbranch_palette Subcluster color palette: 1st color = remaining/unrefined, others cycle for subclusters
#' @param downstream_min_size Downstream subtree size threshold (NULL = adaptive)
#' @param force_split Whether to attempt relaxed acceptance conditions when original community split fails
#' @param main_fraction_cap Maximum community proportion threshold (allows retention of smaller communities even when exceeded)
#' @param core_periph Allow core-periphery splitting
#' @param core_degree_quantile Core degree threshold quantile
#' @param core_min_fraction Minimum core fraction
#' @param degree_gini_threshold Degree Gini threshold that triggers core-periphery splitting
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
#' @param verbose Logical; whether to emit progress messages.
#' @return list(plot, graph, subclusters, subcluster_df, cluster_df, method_subcluster,
#'              base_network, branch_network, params)
#' @keywords internal
.plot_dendro_network_with_branches <- function(
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
                sg <- tryCatch(
                    igraph::induced_subgraph(g, vids = genes_target),
                    error = function(e) {
                        warning("[geneSCOPE::plotRWDendrogram] Subgraph construction failed for cluster ", cid, ": ", e$message)
                        return(NULL)
                    }
                )
                if (is.null(sg) || igraph::vcount(sg) == 0) {
                    return(list(method = "none", subclusters = list()))
                }

                .log_info(
                    parent,
                    step_prep,
                    paste0("Cluster ", cid, " subgraph: nodes=", igraph::vcount(sg), " edges=", igraph::ecount(sg)),
                    verbose
                )

                method_used <- "none"
                subclusters <- list()
                vcount_sg <- igraph::vcount(sg)

                # Try articulation point method
                if (vcount_sg >= 3) {
                    tryCatch(
                        {
                            arts_vs <- igraph::articulation_points(sg)
                            if (length(arts_vs)) {
                                arts_ids <- as.integer(igraph::as_ids(arts_vs))
                                arts_ids <- arts_ids[is.finite(arts_ids) & arts_ids >= 1 & arts_ids <= vcount_sg]

                                cand <- list()
                                for (a_id in arts_ids) {
                                    root_name <- igraph::V(sg)$name[a_id]
                                    sg_minus <- tryCatch(igraph::delete_vertices(sg, a_id), error = function(e) NULL)
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
                                        cand[[length(cand) + 1]] <- list(
                                            size = length(branch_genes),
                                            genes = branch_genes
                                        )
                                    }
                                }

                                if (length(cand)) {
                                    ord <- order(vapply(cand, `[[`, numeric(1), "size"), decreasing = TRUE)
                                    used <- character(0)
                                    kept <- list()
                                    for (i in ord) {
                                        gs <- cand[[i]]$genes
                                        if (!any(gs %in% used)) {
                                            kept[[length(kept) + 1]] <- cand[[i]]
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
                                louvain = igraph::cluster_louvain(sg),
                                leiden = igraph::cluster_leiden(sg),
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
                                    subclusters[[paste0(cid, "_sub", k)]] <- igraph::V(sg)$name[comm$membership == sid]
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

    # 5. Build ordered factor levels (parent cluster -> subclusters)
    numeric_sort_levels <- function(x) {
        ux <- unique(x)
        suppressWarnings({
            nx <- as.numeric(ux)
        })
        if (all(!is.na(nx))) as.character(sort(as.numeric(ux))) else sort(ux)
    }

    base_levels <- numeric_sort_levels(clu)

    if (method_final %in% c("none", "disabled")) {
        cluster_df$cluster <- factor(cluster_df$cluster, levels = base_levels)
    } else {
        sub_levels_raw <- sort(unique(na.omit(sub_attr)))
        parent_of <- sub("^([^_]+)_.*$", "\\1", sub_levels_raw)
        level_vec <- character(0)

        for (pv in base_levels) {
            level_vec <- c(level_vec, pv)
            if (pv %in% parent_of) {
                kids <- sub_levels_raw[parent_of == pv]
                level_vec <- c(level_vec, kids)
            }
        }

        orphan <- setdiff(sub_levels_raw, level_vec)
        if (length(orphan)) level_vec <- c(level_vec, orphan)

        cluster_df$cluster <- factor(cluster_df$cluster, levels = base_levels)
    }

    # 6. Build result tables and colors
    subcluster_df <- data.frame(
        gene = Vnames,
        cluster = cluster_df$cluster,
        subcluster = sub_attr,
        stringsAsFactors = FALSE
    )

    if (!(method_final %in% c("none", "disabled"))) {
        sub_levels_only <- level_vec[level_vec %in% subcluster_df$subcluster]
        subcluster_df$subcluster <- factor(subcluster_df$subcluster, levels = sub_levels_only)
    }

    if (method_final %in% c("none", "disabled")) {
        uniq_clu <- sort(unique(na.omit(clu)))
        col_map <- setNames(rep(cluster_palette, length.out = length(uniq_clu)), uniq_clu)
        node_cols <- col_map[clu]
    } else {
        sub_levels <- sort(unique(na.omit(sub_attr)))
        col_map <- setNames(rep(subbranch_palette[1], length(Vnames)), Vnames)
        if (length(sub_levels)) {
            cols_sub <- rep(subbranch_palette[-1], length.out = length(sub_levels))
            for (i in seq_along(sub_levels)) {
                gs <- names(sub_attr)[sub_attr == sub_levels[i]]
                col_map[gs] <- cols_sub[i]
            }
        }
        node_cols <- col_map[Vnames]
    }
    subcluster_df$color <- node_cols

    # 7. Generate final plot
    .log_info(parent, step_plot, "Generating final plot", verbose)
    branch_network <- NULL

    if (method_final %in% c("none", "disabled")) {
        cluster_vec_base <- factor(clu, levels = base_levels)
        new_args <- base_args
        new_args$cluster_vec <- setNames(cluster_vec_base, names(clu))
        new_args$title <- if (is.null(title)) "Base Network" else title
        branch_network <- do.call(.plot_dendro_network, new_args)
        plt <- branch_network$plot
    } else {
        cluster_vec_sub <- setNames(ifelse(is.na(sub_attr), clu, sub_attr), names(clu))
        cluster_vec_sub <- factor(cluster_vec_sub, levels = level_vec)
        new_args <- base_args
        new_args$cluster_vec <- cluster_vec_sub
        new_args$title <- if (is.null(title)) {
            paste0("Sub-branch: ", method_final)
        } else {
            paste0(title, " | Sub-branch: ", method_final)
        }
        branch_network <- do.call(.plot_dendro_network, new_args)
        plt <- branch_network$plot
    }

    # 8. Return results
    list(
        plot = plt,
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

#' Internal helper for `.plot_density` scalebar sizing.
#' @description
#' Selects a "nice" scalebar length based on the x-range and target fraction.
#' @param x_range Numeric vector of x-range in microns.
#' @param target_frac Fraction of x-range targeted for the bar length.
#' @param nice_um Vector of preferred lengths in microns.
#' @return Numeric length in microns.
#' @keywords internal
.plot_density_choose_scalebar_length <- function(x_range,
                                                target_frac = 0.18,
                                                nice_um = c(10, 20, 50, 100, 200, 400, 500, 800, 1000, 2000, 5000)) {
    x_rng <- range(x_range, na.rm = TRUE)
    dx <- diff(x_rng)
    if (!is.finite(dx) || dx <= 0) return(nice_um[1])
    target_len <- dx * target_frac
    pick <- nice_um[nice_um <= target_len]
    if (!length(pick)) return(nice_um[1])
    max(pick)
}

#' Internal helper for `.plot_density` scalebar placement.
#' @description
#' Computes scalebar coordinates in data space with adaptive padding and gaps.
#' @param x_range Numeric vector of x-range in data units.
#' @param y_range Numeric vector of y-range in data units.
#' @param scale_pos Named numeric vector with x/y fractions (0-1).
#' @param bar_len Scalebar length in data units.
#' @param flip_y Logical indicating reversed y-axis (visual flip).
#' @param bar_offset Additional y-offset as fraction of the y-span.
#' @param pad_x_frac Fractional x padding relative to span.
#' @param pad_y_frac Fractional y padding relative to span.
#' @param label_gap_frac Fractional y gap between bar and label.
#' @return List with draw flag and scalebar coordinates.
#' @keywords internal
.plot_density_compute_scalebar <- function(x_range,
                                          y_range,
                                          scale_pos,
                                          bar_len,
                                          flip_y = FALSE,
                                          bar_offset = 0,
                                          pad_x_frac = 0.04,
                                          pad_y_frac = 0.05,
                                          label_gap_frac = 0.02) {
    pad_x_frac <- if (isTRUE(is.finite(pad_x_frac))) pad_x_frac else 0.04
    pad_y_frac <- if (isTRUE(is.finite(pad_y_frac))) pad_y_frac else 0.05
    label_gap_frac <- if (isTRUE(is.finite(label_gap_frac))) label_gap_frac else 0.02
    x_min <- min(x_range, na.rm = TRUE)
    x_max <- max(x_range, na.rm = TRUE)
    y_min <- min(y_range, na.rm = TRUE)
    y_max <- max(y_range, na.rm = TRUE)
    dx <- x_max - x_min
    dy <- y_max - y_min
    if (!isTRUE(is.finite(dx)) || dx <= 0 || !isTRUE(is.finite(dy)) || dy <= 0) {
        return(list(draw = FALSE))
    }
    if (!isTRUE(is.finite(bar_len)) || bar_len <= 0) {
        return(list(draw = FALSE))
    }
    pad_x <- dx * pad_x_frac
    pad_y <- dy * pad_y_frac
    avail_dx <- dx - 2 * pad_x - bar_len
    if (!isTRUE(is.finite(avail_dx)) || avail_dx < 0) {
        return(list(draw = FALSE))
    }
    scale_pos <- if (is.null(scale_pos)) c(x = 0, y = 0) else scale_pos
    if (!is.null(names(scale_pos)) && all(c("x", "y") %in% names(scale_pos))) {
        scale_pos <- scale_pos[c("x", "y")]
    } else {
        scale_pos <- c(x = scale_pos[1], y = scale_pos[2])
    }
    scale_pos <- suppressWarnings(as.numeric(scale_pos))
    scale_pos[is.na(scale_pos) | !is.finite(scale_pos)] <- 0
    scale_pos <- pmin(1, pmax(0, scale_pos))
    if (length(scale_pos) < 2L) {
        scale_pos <- rep_len(scale_pos, 2L)
    }
    scale_pos <- setNames(scale_pos, c("x", "y"))

    x0 <- x_min + pad_x + avail_dx * scale_pos["x"]
    x1 <- x0 + bar_len
    if (!isTRUE(is.finite(x0)) || !isTRUE(is.finite(x1))) {
        return(list(draw = FALSE))
    }

    bar_offset <- if (isTRUE(is.finite(bar_offset))) bar_offset else 0
    label_gap <- dy * label_gap_frac
    if (!isTRUE(is.finite(label_gap)) || label_gap <= 0) {
        label_gap <- dy * 0.02
    }
    avail_dy <- dy - 2 * pad_y
    if (!isTRUE(is.finite(avail_dy)) || avail_dy < 0) {
        return(list(draw = FALSE))
    }
    y_bar_vis <- pad_y + avail_dy * scale_pos["y"] + bar_offset * dy
    if (!isTRUE(is.finite(y_bar_vis))) {
        y_bar_vis <- pad_y
    }
    y_bar_vis <- min(max(y_bar_vis, pad_y), dy - pad_y)
    label_y_vis <- y_bar_vis + label_gap
    if (!isTRUE(is.finite(label_y_vis))) {
        label_y_vis <- y_bar_vis + label_gap
    }
    if (!isTRUE(is.finite(label_y_vis))) {
        return(list(draw = FALSE))
    }
    if (isTRUE(label_y_vis > dy - pad_y)) label_y_vis <- y_bar_vis - label_gap
    if (isTRUE(label_y_vis < pad_y)) label_y_vis <- y_bar_vis + label_gap
    label_y_vis <- min(max(label_y_vis, pad_y), dy - pad_y)

    to_data_y <- function(y_vis) {
        if (isTRUE(flip_y)) y_max - y_vis else y_min + y_vis
    }

    list(
        draw = TRUE,
        x0 = x0,
        x1 = x1,
        y_bar = to_data_y(y_bar_vis),
        y_label = to_data_y(label_y_vis)
    )
}

#' 4.PlotDensity.r (2025-01-06 rev)
#' @description
#' Internal helper for `.plot_density`.
#' @title .plot_density - visualise grid-level density stored in @density
#' @param scope_obj A \code{scope_object}.
#' @param grid_name Name of the grid layer to plot (e.g. "grid50").
#' @param density1_name Column names to plot.
#' @param density2_name Column names to plot.
#' @param palette1 Color palettes for density1 and density2.
#' @param palette2 Color palettes for density1 and density2.
#' @param alpha1 Alpha transparency for density layers.
#' @param alpha2 Alpha transparency for density layers.
#' @param tile_shape Shape used to render each grid tile: "square", "circle", or "hex".
#' @param hex_orientation Orientation for hex tiles ("flat" = flat-top, "pointy" = pointy-top).
#' @param aspect_ratio Optional numeric width/height ratio to enforce for the grid axes;
#'        when `NULL`, the raw `grid_info` extents are used (no forced square padding).
#' @param scale_bar_pos Optional numeric vector (x,y) with fractional offsets for the scale bar
#'        relative to the plotting window (0 = left/bottom, 1 = right/top).
#' @param scale_bar_show Logical toggle controlling whether the scale bar is drawn.
#' @param scale_bar_colour Colour used for the scale bar line and text (default black).
#' @param scale_bar_corner Corner preset for the scale bar when `scale_bar_pos` is not supplied;
#'        choices are "bottom-left", "bottom-right", "top-left", "top-right".
#' @param use_histology Logical; overlay histology stored on the grid when available.
#' @param histology_level Which histology slot ("lowres" or "hires") to prioritise.
#' @param axis_mode Coordinate space for the plot: `"grid"` (microns) or `"image"` (pixels).
#' @param seg_type Segmentation overlay type: "cell", "nucleus", or "both".
#' @param colour_cell Colors for cell and nucleus segmentation.
#' @param colour_nucleus Colors for cell and nucleus segmentation.
#' @param alpha_seg Alpha transparency for segmentation overlay.
#' @param grid_gap Spacing for background grid lines.
#' @param scale_text_size Font size for scale bar text.
#' @param bar_len Length of scale bar in micrometers.
#' @param bar_offset Vertical offset of scale bar as fraction of y-range.
#' @param arrow_pt Arrow point size for scale bar.
#' @param scale_legend_colour Color of scale bar and text.
#' @param max.cutoff1 Maximum cutoff fractions for density values.
#' @param max.cutoff2 Maximum cutoff fractions for density values.
#' @param legend_digits Number of decimal places in legend.
#' @param overlay_image Parameter value.
#' @param image_path Filesystem path.
#' @param image_alpha Parameter value.
#' @param image_choice Parameter value.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot geom_tile aes scale_fill_gradient guide_colorbar geom_vline geom_hline scale_x_continuous scale_y_continuous coord_fixed theme_minimal theme element_blank element_rect element_text element_line unit margin annotate arrow
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggforce geom_shape
#' @importFrom scales alpha number_format squish
#' @importFrom data.table rbindlist
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
    accuracy_val <- 1 / (10^legend_digits)

    ## ------------------------------------------------------------------ 1
    ## validate grid layer & retrieve geometry
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
    prepare_rgba <- function(img, alpha_scale = 1) {
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
    flip_raster_vertical <- function(img) {
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

    ## ------------------------------------------------------------------ 2
    ## locate density data frame (new slot first, old slot as fallback)
    densityDF <- scope_obj@density[[grid_layer_name]]
    if (is.null(densityDF)) {
        densityDF <- g_layer$densityDF
    } # legacy
    if (is.null(densityDF)) {
        stop("No density table found for grid '", grid_layer_name, "'.")
    }

    ## make sure rownames are grid IDs
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

    ## helper to assemble heatmap data.frame
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
    heat2 <- if (!is.null(density2_name)) build_heat(density2_name, max.cutoff2)
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
    }

    ## ------------------------------------------------------------------ 3
    ## tile geometry (centre & size)
    tile_df <- function(df) {
        transform(df,
            x = (xmin + xmax) / 2,
            y = (ymin + ymax) / 2,
            w = pmax(0, xmax - xmin),
            h = pmax(0, ymax - ymin)
        )
    }

    build_shape_geom <- function(df, alpha_val) {
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

    resolve_scale_bar_pos <- function(pos, corner, default_x = 0, default_y = 0) {
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

    library(ggplot2)
    p <- ggplot()

    ## ------------------------------------------------------------------ 3.0
    ## Histology background (attached to scope_obj)
    if (histology_available) {
        img_rgba <- prepare_rgba(histology_slot$png, image_alpha)
        if (!is.null(img_rgba)) {
            if (!use_image_coords && isTRUE(flip_histology_y)) {
                img_rgba <- flip_raster_vertical(img_rgba)
            }
            if (use_image_coords) {
                p <- p + annotation_raster(
                    img_rgba,
                    xmin = 0, xmax = histology_slot$width,
                    ymin = 0, ymax = histology_slot$height,
                    interpolate = TRUE
                )
            } else {
                roi_draw <- if (!is.null(histology_roi)) histology_roi else default_roi
                p <- p + annotation_raster(
                    img_rgba,
                    xmin = roi_draw["xmin"], xmax = roi_draw["xmax"],
                    ymin = roi_draw["ymin"], ymax = roi_draw["ymax"],
                    interpolate = TRUE
                )
            }
        }
    } else if (isTRUE(overlay_image)) {
        img_info <- g_layer$image_info
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
                eff_mpp <- mpp_fullres
                is_hires <- grepl("tissue_hires_image", basename(sel_path), ignore.case = TRUE)
                is_lowres <- grepl("tissue_lowres_image", basename(sel_path), ignore.case = TRUE)
                is_aligned <- grepl("aligned_tissue_image", basename(sel_path), ignore.case = TRUE)
                if (is_hires && is.finite(hires_sf)) eff_mpp <- mpp_fullres / hires_sf
                if (is_lowres && is.finite(lowres_sf)) eff_mpp <- mpp_fullres / lowres_sf
                if (is_aligned && is.finite(aligned_sf)) eff_mpp <- mpp_fullres / aligned_sf
                if (is_aligned && !is.finite(aligned_sf)) eff_mpp <- mpp_fullres
                if (!is.finite(eff_mpp) || eff_mpp <= 0) eff_mpp <- mpp_fullres
                w_um <- as.numeric(wpx) * as.numeric(eff_mpp)
                h_um <- as.numeric(hpx) * as.numeric(eff_mpp)
                y_origin_manual <- if (!is.null(img_info$y_origin)) img_info$y_origin else "top-left"
                if (identical(y_origin_manual, "top-left")) {
                    img <- flip_raster_vertical(img)
                }
                p <- p + annotation_raster(img,
                    xmin = 0, xmax = w_um,
                    ymin = 0, ymax = h_um,
                    interpolate = TRUE
                )
            }
        }
    }

    shape_layers1 <- build_shape_geom(tile_df(heat1), alpha1)
    p <- Reduce(`+`, shape_layers1, init = p)
    p <- p +
        scale_fill_gradient(
            name = density1_name,
            low = "transparent",
            high = palette1,
            limits = c(0, unique(heat1$cut)),
            oob = scales::squish,
            labels = scales::number_format(accuracy = accuracy_val),
            na.value = "transparent",
            guide = guide_colorbar(order = 1)
        )

    ## second density with new fill scale
    if (!is.null(heat2)) {
        library(ggnewscale)
        shape_layers2 <- build_shape_geom(tile_df(heat2), alpha2)
        p <- p + new_scale_fill()
        p <- Reduce(`+`, shape_layers2, init = p)
        p <- p +
            scale_fill_gradient(
                name = density2_name,
                low = "transparent",
                high = palette2,
                limits = c(0, unique(heat2$cut)),
                oob = scales::squish,
                labels = scales::number_format(accuracy = accuracy_val),
                na.value = "transparent",
                guide = guide_colorbar(order = 2)
            )
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
        if (use_image_coords) {
            seg_df <- .coords_physical_to_level(as.data.frame(seg_dt),
                x_col = "x", y_col = "y", histo = histology_slot
            )
            seg_dt <- as.data.table(seg_df)
            set(seg_dt, j = "x", value = seg_dt$x_img)
            set(seg_dt, j = "y", value = seg_dt$y_img)
            seg_dt[, c("x_img", "y_img") := NULL]
        }
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
    ## aesthetics: square, gridlines, scale‑bar
    if (!use_image_coords) {
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
        x_rng <- ensure_extent(x_rng_phys)
        y_rng <- ensure_extent(y_rng_phys)
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
        x_rng <- c(0, histology_slot$width)
        y_rng <- if (flip_histology_y) c(histology_slot$height, 0) else c(0, histology_slot$height)
        grid_x <- phys_to_img(seq(x_rng_phys[1], x_rng_phys[2], by = grid_gap), axis = "x")
        grid_y <- phys_to_img(seq(y_rng_phys[1], y_rng_phys[2], by = grid_gap), axis = "y")
        grid_x <- grid_x[is.finite(grid_x)]
        grid_y <- grid_y[is.finite(grid_y)]
    }

    if (length(grid_x)) {
        p <- p + geom_vline(xintercept = grid_x, linewidth = 0.05, colour = "grey80")
    }
    if (length(grid_y)) {
        p <- p + geom_hline(yintercept = grid_y, linewidth = 0.05, colour = "grey80")
    }

    if (use_image_coords) {
        if (flip_histology_y) {
            p <- p +
                scale_x_continuous(limits = c(0, histology_slot$width), expand = c(0, 0)) +
                scale_y_reverse(limits = c(histology_slot$height, 0), expand = c(0, 0)) +
                coord_fixed(expand = FALSE)
        } else {
            p <- p +
                scale_x_continuous(limits = c(0, histology_slot$width), expand = c(0, 0)) +
                scale_y_continuous(limits = c(0, histology_slot$height), expand = c(0, 0)) +
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
            legend.key.width = unit(0.4, "cm"), # legend key width
            legend.key.height = unit(0.15, "cm"), # legend key height
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

    ## scale‑bar
    draw_scale_bar <- isTRUE(scale_bar_show)
    scale_pos <- resolve_scale_bar_pos(scale_bar_pos, scale_bar_corner)
    bar_len_um <- suppressWarnings(as.numeric(bar_len[1]))
    if (is.null(bar_len) || !is.finite(bar_len_um)) {
        bar_len_um <- .plot_density_choose_scalebar_length(
            if (use_image_coords) histology_roi[c("xmin", "xmax")] else x_rng_phys
        )
    }
    if (!is.finite(bar_len_um) || bar_len_um <= 0) draw_scale_bar <- FALSE

    if (isTRUE(draw_scale_bar)) {
        if (use_image_coords) {
            roi <- histology_roi
            span_um <- roi["xmax"] - roi["xmin"]
            px_per_um <- histology_slot$width / span_um
            if (!is.finite(px_per_um) || px_per_um <= 0) {
                draw_scale_bar <- FALSE
            } else {
                bar_len_display <- bar_len_um * px_per_um
                coords <- .plot_density_compute_scalebar(
                    x_range = c(0, histology_slot$width),
                    y_range = if (flip_histology_y) c(histology_slot$height, 0) else c(0, histology_slot$height),
                    scale_pos = scale_pos,
                    bar_len = bar_len_display,
                    flip_y = flip_histology_y,
                    bar_offset = bar_offset
                )
            }
        } else {
            coords <- .plot_density_compute_scalebar(
                x_range = x_rng,
                y_range = y_rng,
                scale_pos = scale_pos,
                bar_len = bar_len_um,
                flip_y = FALSE,
                bar_offset = bar_offset
            )
        }
        if (isTRUE(draw_scale_bar) && isTRUE(coords$draw)) {
            label <- paste0(bar_len_um, " \u00B5m")
            scale_df <- data.frame(
                x0 = coords$x0,
                x1 = coords$x1,
                y_bar = coords$y_bar,
                y_label = coords$y_label,
                label = label,
                stringsAsFactors = FALSE
            )
            line_width <- max(0.6, min(1, scale_text_size / 3))
            p <- p +
                geom_segment(
                    data = scale_df,
                    aes(x = x0, xend = x1, y = y_bar, yend = y_bar),
                    arrow = arrow(length = unit(arrow_pt, "pt"), ends = "both"),
                    colour = scale_bar_colour,
                    linewidth = line_width,
                    inherit.aes = FALSE
                ) +
                geom_text(
                    data = scale_df,
                    aes(x = (x0 + x1) / 2, y = y_label, label = label),
                    colour = scale_bar_colour,
                    vjust = 1,
                    size = scale_text_size,
                    inherit.aes = FALSE
                )
        }
    }
    return(p)
}

#' Plot Gene Centroid Expression with Segmentation Overlay
#' @description
#'   Plots cell centroids coloured by the normalised expression of two specified genes.
#'   Expression values are extracted from `scope_obj@cells$counts`, min-max scaled
#'   to 0-1 separately for each gene, and mapped to point transparency (`alpha`) so
#'   that higher expression appears more opaque. Centroids of cells expressing
#'   `gene1` are drawn with `palette1`; centroids of cells expressing `gene2` with
#'   `palette2`. Cells that express neither gene are not shown. Optionally overlays
#'   cell/nucleus segmentation polygons.
#' @param scope_obj A `scope_object` containing `coord$centroids`, segmentation
#'                   vertices, and `cells$counts`.
#' @param gene1_name Character. Name of the first gene.
#' @param gene2_name Character. Name of the second gene.
#' @param palette1 Colour used for cells expressing `gene1`. Default "#fc3d5d".
#' @param palette2 Colour used for cells expressing `gene2`. Default "#4753f8".
#' @param size1 Point size for `gene1` cells. Default 0.3.
#' @param size2 Point size for `gene2` cells. Default 0.3.
#' @param alpha1 Maximum alpha for `gene1` cells. Default 0.7.
#' @param alpha2 Maximum alpha for `gene2` cells. Default 0.7.
#' @param seg_type Segmentation overlay type: "none", "cell", "nucleus",
#'                   or "both". Default "none".
#' @param colour_cell Line colour for cell segmentation. Default "black".
#' @param colour_nucleus Line colour for nucleus segmentation. Default "#3182bd".
#' @param alpha_seg Alpha for segmentation polygons. Default 0.2.
#' @param grid_gap Spacing of background grid lines (µm). Default 100.
#' @param scale_text_size Font size for scale-bar text. Default 2.4.
#' @param bar_len Length of the scale bar (µm). Default 400.
#' @param bar_offset Vertical offset of the scale bar as a fraction of the y-range. Default 0.01.
#' @param arrow_pt Arrow-head size of the scale bar. Default 4.
#' @param scale_legend_colour Colour of the scale-bar text and line. Default "black".
#' @param max.cutoff1 Fraction (0-1) of gene1's maximum expression used for clipping. Default 1.
#' @param max.cutoff2 Fraction (0-1) of gene2's maximum expression used for clipping. Default 1.
#' @return           A `ggplot` object.
#' @importFrom ggplot2 ggplot geom_point aes scale_alpha_continuous scale_x_continuous scale_y_continuous coord_fixed theme_minimal theme element_blank element_rect element_line element_text unit margin annotate arrow
#' @importFrom ggforce geom_shape
#' @importFrom ggnewscale new_scale
#' @importFrom data.table rbindlist
#' @importFrom scales alpha
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

    # scale‑bar
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

#' .plot_grid_boundary
#' @description
#' Internal helper for `.plot_grid_boundary`.
#' Draw grid cell boundaries for a given grid layer in a coord object and
#' guarantee that the output image keeps the exact spatial aspect ratio
#' of the original coordinates (no stretching or squeezing), following the
#' same padding logic used in `.plot_density()`.
#' @param scope_obj A coord-style S4 object that contains a `@grid` slot,
#'                   e.g. produced by the STUtility pipeline.
#' @param grid_name Character. The name of the grid layer to plot, e.g.
#'                   "grid31". Must exist inside `scope_obj@grid`.
#' @param colour Border colour for grid rectangles. Default "black".
#' @param linewidth Border line width for rectangles. Default 0.2.
#' @param panel_bg Background colour of the panel. Default "#C0C0C0".
#' @param base_size Base font size for `theme_minimal()`. Default 10.
#' @return A `ggplot` object.
#' @importFrom ggplot2 ggplot geom_rect aes coord_fixed theme_minimal theme element_text element_rect element_blank element_line labs unit
#' @keywords internal
.plot_grid_boundary <- function(scope_obj,
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

#' Plot Morisita Iδ peRr cluster
#' @description
#'   Creates a faceted scatter-line plot showing Morisita's Iδ values for the
#'   top genes in each cluster of a chosen grid layer.  Clusters can be
#'   down-sampled, re-ordered to balance facet rows, and coloured by a custom
#'   palette.
#' @param scope_obj A \code{scope_object} containing \code{@grid} and
#'                        \code{@meta.data}.
#' @param grid_name Character. Grid sub-layer name. If \code{NULL} and
#'                        only one layer exists, that layer is used.
#' @param cluster_col Column name in \code{meta.data} that assigns each
#'                        gene to a cluster.
#'                        raw or normalised Iδ.
#' @param top_n Optional integer. Keep the top \code{n} genes
#'                        (highest Iδ) per cluster.
#' @param nrow Integer. Number of facet rows to aim for.
#' @param min_genes Integer. Minimum number of genes a cluster must keep
#'                        after filtering (default 1).
#'                        unnamed the order is taken as-is; otherwise names
#'                        are matched to cluster labels.
#' @param point_size Numeric. Aesthetics for the plot.
#' @param line_size Numeric. Aesthetics for the plot.
#' @param label_size Numeric. Aesthetics for the plot.
#' @param subCluster Optional character vector. Restrict to these cluster
#'                        labels.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_wrap labs theme_minimal theme element_blank element_rect margin coord_cartesian unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by slice_max ungroup filter arrange desc mutate
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
    label_size = 2, # default reduced to 2
    subCluster = NULL) {
    meta_col <- paste0(grid_name, "_iDelta")
    if (!meta_col %in% colnames(scope_obj@meta.data)) {
        stop("Cannot find Iδ values: meta.data column '", meta_col, "' is missing.")
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
        stop("No cluster meets min_genes ≥ ", min_genes, ".")
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
            box.padding = unit(0.15, "lines"),
            point.padding = unit(0.15, "lines"),
            segment.size = 0.3,
            segment.color = "grey50",
            force = 0.5,
            max.overlaps = Inf
        ) +
        facet_wrap(~cluster, scales = "free_x", nrow = nrow, strip.position = "bottom") +
        labs(
            x     = NULL, # remove x-axis label as well
            y     = expression(I[delta]),
            title = paste0("Iδ by Cluster (", grid_name, ")")
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

#' Plot Distribution of Lee's L Values
#' @description
#'   Extract the Lee's L matrix from a specified grid sublayer and LeeStats layer within `scope_obj`.
#'   Bin the off-diagonal elements of the Lee's L matrix into intervals and plot a histogram to
#'   visualize the distribution of Lee's L values. Optionally, take the absolute value of Lee's L
#'   before plotting. The plot adheres to publication-quality aesthetic standards.
#' @param scope_obj An object in which Lee's L has already been computed and stored under
#'                         `scope_obj@grid[[grid_name]][[lee_stats_layer]]$L`.
#' @param grid_name Character (optional). The name of the grid sublayer to use (e.g., `"grid50"`).
#'                         If NULL and there is only one sublayer under `scope_obj@grid`, that layer is used.
#'                         Otherwise, this parameter is required.
#' @param lee_stats_layer Character (optional). Name of the LeeStats layer to use (e.g., `"LeeStats_Xz"`).
#'                         If NULL, attempts to automatically detect a single layer whose name starts
#'                         with `"LeeStats_"`. If multiple matching layers are found, an error is thrown
#'                         and the user must specify `lee_stats_layer`.
#' @param bins Integer or numeric vector. Passed to `geom_histogram()` as bin settings.
#'                         If a single integer, that number of equal-width bins is used. If a numeric
#'                         vector of breakpoints, those defines the histogram bins.
#' @param xlim Numeric vector of length 2 (optional). X-axis display limits (c(min, max)).
#'                         Defaults to NULL (automatic scaling).
#' @param title Character (optional). The plot title. Defaults to `"Lee's L Distribution"` or
#'                         `"Absolute Lee's L Distribution"` if `use_abs = TRUE`.
#' @param use_abs Logical. Whether to take the absolute value of Lee's L before plotting.
#'                         Defaults to FALSE (plot original Lee's L).
#' @return A `ggplot2` object displaying the histogram of Lee's L values (or their absolute values).
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_bw theme element_rect element_line element_text coord_cartesian expansion
#' @importFrom grid unit
#' @keywords internal
.plot_l_distribution <- function(scope_obj,
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

#' Plot Lee's L vs. Gene Total Count Fold Change (All Gene Pairs)
#' @description
#'   For a given `scope_obj`, extract the Lee's L matrix from
#'   `scope_obj@grid[[grid_name]][[lee_stats_layer]]$L`.
#'   Summarize each gene's total count across all grid cells from
#'   `scope_obj@grid[[grid_name]]$counts`. For every gene pair (i, j) with i < j,
#'   compute Lee's L and the fold change of their total counts
#'   (fold change = max(total_i, total_j) / min(total_i, total_j), ≥ 1).
#'   Finally, plot a scatterplot with Lee's L on the x-axis and fold change on
#'   the y-axis, using a 25% gray background and publication-quality aesthetics.
#' @param scope_obj An object containing
#'                        `scope_obj@grid[[grid_name]]$counts` and
#'                        `scope_obj@grid[[grid_name]][[lee_stats_layer]]$L`.
#' @param grid_name Character (optional). The name of the grid sublayer to use.
#'                        If NULL and `scope_obj@grid` has exactly one sublayer,
#'                        that sublayer is used automatically. Otherwise, it must be provided.
#' @param lee_stats_layer Character (optional). The name of the LeeStats layer (e.g., `"LeeStats_Xz"`).
#'                        If NULL, the function searches for exactly one layer under
#'                        `scope_obj@grid[[grid_name]]` whose name starts with `"LeeStats_"`.
#'                        If multiple matches are found, an error is thrown and the user must
#'                        supply `lee_stats_layer` explicitly.
#' @param title Character (optional). Plot title. Defaults to
#'                        `"Lee's L vs. Fold Change (All Gene Pairs)"` if NULL.
#' @return A `ggplot2` object showing a scatterplot of Lee's L vs. fold change
#'         for all gene pairs (i < j).
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_rect element_line element_text coord_cartesian scale_x_continuous scale_y_continuous expansion
#' @importFrom dplyr group_by summarise filter left_join rename rowwise ungroup
#' @importFrom tidyr crossing
#' @importFrom grid unit
#' @keywords internal
.plot_l_scatter <- function(scope_obj,
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
      y = "Fold Change (≥1)"
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

#' Scatter plot of Lee's L versus Pearson correlation
#' @description
#'   Builds a scatter plot comparing Lee's spatial correlation statistic
#'   (L) with the conventional Pearson \(r\) for every gene pair in a
#'   selected grid layer or at the cell level.  Optionally flips axes and
#'   annotates the largest positive and negative L–rplotLvsR differences
#'   (\code{Delta}).
#' @param scope_obj A \code{scope_object} containing both Lee's L matrices
#'                    (in \code{LeeStats_Xz}) and Pearson correlation
#'                    matrices.
#' @param grid_name Character. Grid sub-layer name.
#' @param pear_level Character, \code{"cell"} or \code{"grid"}; which
#'                    Pearson matrix to use.
#' @param lee_stats_layer Character. Name of the Lee's L stats layer to use
#'                    (default "LeeStats_Xz").
#' @param delta_top_n Integer. How many extreme \code{Delta} pairs to label
#'                    on the plot.
#' @param flip Logical. If \code{TRUE}, put Pearson on the x-axis and
#'                    Lee's L on the y-axis.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{.get_top_l_vs_r}}
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline scale_x_continuous scale_y_continuous labs theme_minimal theme element_text element_rect element_line element_blank
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr mutate bind_rows slice_max slice_min distinct
#' @importFrom scales label_number
#' @keywords internal
.plot_l_vs_r <- function(scope_obj,
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

  ## ---- 2. Convert to long table + Δ ---------------------------------------------------------
  genes <- colnames(L_mat)
  ut <- upper.tri(L_mat, diag = FALSE)
  df_long <- data.frame(
    gene1 = rep(genes, each = length(genes))[ut],
    gene2 = rep(genes, length(genes))[ut],
    LeesL = L_mat[ut],
    Pear  = r_mat[ut]
  ) |>
    mutate(Delta = LeesL - Pear)

  ## ---- 3. Label points (extreme Δ) ----------------------------------------------------
  df_label <- bind_rows(
    slice_max(df_long, Delta, n = delta_top_n, with_ties = FALSE),
    slice_min(df_long, Delta, n = delta_top_n, with_ties = FALSE)
  ) |>
    distinct(gene1, gene2, .keep_all = TRUE) |>
    mutate(label = sprintf("%s–%s\nL=%.3f", gene1, gene2, LeesL))

  ## ---- 4. Plot ----------------------------------------------------------------
  if (!flip) {
    p <- ggplot(
      df_long,
      aes(x = LeesL, y = Pear)
    )
    xlab <- "Lee’s L"
    ylab <- "Pearson correlation"
    ttl <- sprintf(
      "Lee’s L vs Pearson  (%s, %s)",
      sub("grid", "Grid ", grid_name),
      ifelse(pear_level == "cell", "cell level", "grid level")
    )
  } else {
    p <- ggplot(
      df_long,
      aes(x = Pear, y = LeesL)
    )
    xlab <- "Pearson correlation"
    ylab <- "Lee’s L"
    ttl <- sprintf(
      "Pearson vs Lee’s L  (%s, %s)",
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

.plot_dendro_network_impl <- .plot_dendro_network

.plot_dendro_network_multi_impl <- .plot_dendro_network_multi

.plot_dendro_network_with_branches_impl <- .plot_dendro_network_with_branches

.plot_density_impl <- .plot_density

.plot_density_centroids_impl <- .plot_density_centroids

.plot_grid_boundary_impl <- .plot_grid_boundary

.plot_idelta_impl <- .plot_idelta

.plot_network_impl <- .plot_network

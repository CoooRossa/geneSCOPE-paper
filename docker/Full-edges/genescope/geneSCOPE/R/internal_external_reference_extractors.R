#' Extract igraph from supported inputs.
#' @keywords internal
.extract_graph_any <- function(x,
                               grid_name = "grid30",
                               stats_layer = "LeeStats_Xz",
                               graph_slot = "g_consensus") {
    if (is.list(x) && !is.null(x$graph) && igraph::is_igraph(x$graph)) {
        return(x$graph)
    }
    if (igraph::is_igraph(x)) {
        return(x)
    }
    if (inherits(x, "scope_object")) {
        return(.extract_graph_from_scope(
            scope_obj = x,
            grid_name = grid_name,
            stats_layer = stats_layer,
            graph_slot = graph_slot
        ))
    }
    stop("Unsupported input for graph extraction. Provide a scope_object, igraph, or list with $graph.")
}

#' Extract graph from plotDendroNetwork() output.
#' @keywords internal
.extract_graph_from_dendro_return <- function(x) {
    if (is.list(x) && !is.null(x$graph) && igraph::is_igraph(x$graph)) {
        return(x$graph)
    }
    stop("Input does not contain a valid $graph igraph object.")
}

#' Normalize inputs to an igraph.
#' @keywords internal
.normalize_input_to_igraph <- function(x,
                                       grid_name = "grid30",
                                       stats_layer = "LeeStats_Xz",
                                       graph_slot = "g_consensus") {
    .extract_graph_any(
        x = x,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot
    )
}

#' Extract graph from scope_object using stats path.
#' @keywords internal
.extract_graph_from_scope <- function(scope_obj,
                                      grid_name = "grid30",
                                      stats_layer = "LeeStats_Xz",
                                      graph_slot = "g_consensus") {
    if (!inherits(scope_obj, "scope_object")) {
        stop("scope_obj must be a scope_object.")
    }
    if (is.null(scope_obj@stats) || !length(scope_obj@stats)) {
        stop("scope_obj@stats is empty; run computeL()/clusterGenes() first.")
    }
    grid_names <- names(scope_obj@stats)
    if (is.null(grid_names) || !length(grid_names)) {
        stop("scope_obj@stats has no named grid layers.")
    }
    if (!grid_name %in% grid_names) {
        stop("grid_name not found in scope_obj@stats. Available: ", paste(grid_names, collapse = ", "))
    }
    stats_layer_obj <- scope_obj@stats[[grid_name]]
    stats_layers <- names(stats_layer_obj)
    if (is.null(stats_layers) || !length(stats_layers)) {
        stop("No stats layers found for grid '", grid_name, "'.")
    }
    if (!stats_layer %in% stats_layers) {
        stop("stats_layer not found for grid '", grid_name, "'. Available: ", paste(stats_layers, collapse = ", "))
    }
    layer_obj <- stats_layer_obj[[stats_layer]]
    layer_names <- names(layer_obj)
    if (is.null(layer_names) || !length(layer_names)) {
        stop("stats_layer '", stats_layer, "' is empty for grid '", grid_name, "'.")
    }
    if (!graph_slot %in% layer_names) {
        stop("graph_slot not found in stats_layer '", stats_layer, "'. Available: ", paste(layer_names, collapse = ", "))
    }
    g <- layer_obj[[graph_slot]]
    if (is.null(g) || !igraph::is_igraph(g)) {
        stop("Graph slot '", graph_slot, "' not found or not an igraph in stats layer.")
    }
    g
}

#' Extract membership vector from supported inputs.
#' @keywords internal
.extract_membership_any <- function(x,
                                    graph,
                                    membership = NULL,
                                    grid_name = "grid30",
                                    stats_layer = "LeeStats_Xz",
                                    graph_slot = "g_consensus") {
    if (!igraph::is_igraph(graph)) {
        stop("graph must be an igraph object.")
    }
    vnames <- igraph::V(graph)$name
    if (is.null(vnames)) {
        stop("Graph vertices missing names (V(graph)$name).")
    }

    memb <- NULL
    if (!is.null(membership)) {
        memb <- membership
    } else {
        for (attr in c("module", "cluster", "membership")) {
            vals <- igraph::vertex_attr(graph, attr, index = igraph::V(graph))
            if (!is.null(vals)) {
                memb <- vals
                names(memb) <- vnames
                break
            }
        }
    }

    if (is.null(memb) && inherits(x, "scope_object")) {
        stats_root <- x@stats
        if (!is.null(stats_root) && grid_name %in% names(stats_root)) {
            stats_layers <- stats_root[[grid_name]]
            if (!is.null(stats_layers) && stats_layer %in% names(stats_layers)) {
                stats_leaf <- stats_layers[[stats_layer]]
                if (!is.null(stats_leaf)) {
                    cand <- c(
                        "membership", "gene_membership", "membership_final",
                        "consensus_membership", "clusters", "partition"
                    )
                    for (nm in cand) {
                        if (!is.null(stats_leaf[[nm]])) {
                            memb <- stats_leaf[[nm]]
                            break
                        }
                    }
                }
            }
        }
    }

    if (is.null(memb)) {
        stop("No membership found. Please provide membership=...")
    }

    if (is.factor(memb)) memb <- as.character(memb)
    if (is.data.frame(memb)) {
        if (ncol(memb) != 1L) stop("membership data.frame must have a single column.")
        memb <- memb[[1]]
    }

    if (is.null(names(memb))) {
        if (length(memb) == length(vnames)) {
            names(memb) <- vnames
        } else {
            stop("membership must be a named vector or align with V(graph)$name.")
        }
    }

    keep <- intersect(names(memb), vnames)
    memb <- memb[keep]
    memb <- memb[!is.na(memb)]
    if (!length(memb)) {
        stop("membership has no overlap with graph vertices after filtering.")
    }
    memb
}

#' Extract standardized edge table from an igraph.
#' @keywords internal
.extract_edge_table_from_igraph <- function(graph, edge_weight = "weight") {
    edge_df <- .graph_edge_table(graph, edge_weight = edge_weight)
    if (!nrow(edge_df)) {
        return(data.frame(
            geneA = character(0),
            geneB = character(0),
            weight = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    if (!is.null(edge_weight) && edge_weight %in% names(edge_df)) {
        weight <- edge_df[[edge_weight]]
    } else if ("weight" %in% names(edge_df)) {
        weight <- edge_df$weight
    } else {
        weight <- edge_df$score
    }
    data.frame(
        geneA = edge_df$from,
        geneB = edge_df$to,
        weight = as.numeric(weight),
        stringsAsFactors = FALSE
    )
}

#' Extract graph for module-quality scoring.
#' @keywords internal
.mq_extract_graph <- function(x,
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

#' Extract membership vector for module-quality scoring.
#' @keywords internal
.mq_extract_membership <- function(x,
                                   graph,
                                   membership = NULL,
                                   grid_name = "grid30",
                                   stats_layer = "LeeStats_Xz",
                                   graph_slot = "g_consensus") {
    .extract_membership_any(
        x = x,
        graph = graph,
        membership = membership,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot
    )
}

#' Standardize edge weights for module-quality scoring.
#' @keywords internal
.mq_standardize_weights <- function(graph, edge_weight = "weight") {
    edge_df <- .graph_edge_table(graph, edge_weight = edge_weight)
    if (!nrow(edge_df)) {
        return(data.frame(
            from = character(0),
            to = character(0),
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
    weight <- as.numeric(weight)
    weight[!is.finite(weight)] <- NA_real_
    data.frame(
        from = edge_df$from,
        to = edge_df$to,
        weight = weight,
        stringsAsFactors = FALSE
    )
}

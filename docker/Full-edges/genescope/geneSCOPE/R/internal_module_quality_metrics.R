#' Build edge table with score column from an igraph.
#' @keywords internal
.graph_edge_table <- function(graph, edge_weight = "weight") {
    edge_df <- igraph::as_data_frame(graph, what = "edges")
    if (!nrow(edge_df)) {
        edge_df$score <- numeric(0)
        return(edge_df)
    }
    if (!is.null(edge_weight) && edge_weight %in% names(edge_df)) {
        score <- edge_df[[edge_weight]]
    } else if ("weight" %in% names(edge_df)) {
        score <- edge_df$weight
    } else {
        score <- rep(1, nrow(edge_df))
    }
    score <- as.numeric(score)
    score[!is.finite(score)] <- NA_real_
    edge_df$score <- abs(score)
    edge_df <- edge_df[edge_df$from != edge_df$to, , drop = FALSE]
    edge_df
}

#' Compute weighted modularity from a membership vector.
#' @keywords internal
.compute_weighted_modularity <- function(graph, membership, edge_weight = "weight") {
    if (!igraph::is_igraph(graph)) return(NA_real_)
    vnames <- igraph::V(graph)$name
    if (is.null(vnames)) return(NA_real_)
    memb <- membership[vnames]
    if (!is.numeric(memb)) {
        memb <- as.integer(factor(memb))
    }
    keep <- !is.na(memb)
    if (sum(keep) < 2) return(NA_real_)
    g_sub <- igraph::induced_subgraph(graph, vids = which(keep))
    if (igraph::ecount(g_sub) == 0) return(NA_real_)
    weight <- igraph::edge_attr(g_sub, edge_weight)
    if (is.null(weight)) weight <- igraph::edge_attr(g_sub, "weight")
    if (is.null(weight)) weight <- rep(1, igraph::ecount(g_sub))
    igraph::modularity(g_sub, membership = memb[keep], weights = weight)
}

#' Compute within/between STRING score summaries for modules.
#' @keywords internal
.module_within_between <- function(edge_df, membership, score_col = "string_score") {
    m_from <- membership[edge_df$from]
    m_to <- membership[edge_df$to]
    score <- edge_df[[score_col]]
    keep <- !is.na(m_from) & !is.na(m_to) & !is.na(score)
    if (!any(keep)) {
        empty <- data.frame(
            module = character(0),
            module_size = integer(0),
            within_string_mean = numeric(0),
            stringsAsFactors = FALSE
        )
        return(list(
            within_mean = NA_real_,
            between_mean = NA_real_,
            delta = NA_real_,
            per_module = empty
        ))
    }
    m_from <- m_from[keep]
    m_to <- m_to[keep]
    score <- score[keep]

    within <- m_from == m_to
    within_scores <- score[within]
    between_scores <- score[!within]
    within_mean <- if (length(within_scores)) mean(within_scores) else NA_real_
    between_mean <- if (length(between_scores)) mean(between_scores) else NA_real_

    per_module <- NULL
    mods <- sort(unique(membership[!is.na(membership)]))
    if (length(mods)) {
        per_module <- lapply(mods, function(mod) {
            in_from <- m_from == mod
            in_to <- m_to == mod
            within_idx <- in_from & in_to
            between_idx <- xor(in_from, in_to)
            within_mean_mod <- if (any(within_idx)) mean(score[within_idx]) else NA_real_
            between_mean_mod <- if (any(between_idx)) mean(score[between_idx]) else NA_real_
            data.frame(
                module = mod,
                module_size = sum(membership == mod, na.rm = TRUE),
                within_string_mean = within_mean_mod,
                between_string_mean = between_mean_mod,
                delta_within_between = within_mean_mod - between_mean_mod,
                stringsAsFactors = FALSE
            )
        })
        per_module <- do.call(rbind, per_module)
    } else {
        per_module <- data.frame(
            module = character(0),
            module_size = integer(0),
            within_string_mean = numeric(0),
            between_string_mean = numeric(0),
            delta_within_between = numeric(0),
            stringsAsFactors = FALSE
        )
    }

    list(
        within_mean = within_mean,
        between_mean = between_mean,
        delta = within_mean - between_mean,
        per_module = per_module
    )
}

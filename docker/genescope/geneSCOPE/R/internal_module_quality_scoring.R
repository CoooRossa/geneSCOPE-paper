#' Standardize membership labels for module-quality scoring.
#' @keywords internal
.mq_standardize_membership <- function(membership) {
    if (is.factor(membership)) membership <- as.character(membership)
    if (is.data.frame(membership)) {
        if (ncol(membership) != 1L) {
            stop("membership data.frame must have a single column.")
        }
        membership <- membership[[1]]
    }
    if (is.null(names(membership))) {
        stop("membership must be a named vector aligned to graph vertices.")
    }
    membership <- membership[!is.na(membership)]
    if (!length(membership)) {
        stop("membership has no non-missing entries.")
    }
    module_levels <- unique(as.character(membership))
    module_id <- match(as.character(membership), module_levels)
    names(module_id) <- names(membership)
    list(
        module_id = module_id,
        module_levels = module_levels,
        module_labels = as.character(membership)
    )
}

#' Summarize within/between edge contributions per node and module.
#' @keywords internal
.mq_edge_summaries <- function(edge_df, membership_id, module_levels) {
    if (!nrow(edge_df)) {
        empty_nodes <- data.frame(
            gene = character(0),
            module = character(0),
            module_id = integer(0),
            within_sum = numeric(0),
            within_n = integer(0),
            between_sum = numeric(0),
            between_n = integer(0),
            stringsAsFactors = FALSE
        )
        empty_modules <- data.frame(
            module = character(0),
            module_id = integer(0),
            module_size = integer(0),
            within_sum = numeric(0),
            within_n = integer(0),
            between_sum = numeric(0),
            between_n = integer(0),
            stringsAsFactors = FALSE
        )
        return(list(nodes = empty_nodes, modules = empty_modules))
    }

    m_from <- membership_id[edge_df$from]
    m_to <- membership_id[edge_df$to]
    keep <- !is.na(m_from) & !is.na(m_to) & is.finite(edge_df$weight)
    if (!any(keep)) {
        empty_nodes <- data.frame(
            gene = character(0),
            module = character(0),
            module_id = integer(0),
            within_sum = numeric(0),
            within_n = integer(0),
            between_sum = numeric(0),
            between_n = integer(0),
            stringsAsFactors = FALSE
        )
        empty_modules <- data.frame(
            module = character(0),
            module_id = integer(0),
            module_size = integer(0),
            within_sum = numeric(0),
            within_n = integer(0),
            between_sum = numeric(0),
            between_n = integer(0),
            stringsAsFactors = FALSE
        )
        return(list(nodes = empty_nodes, modules = empty_modules))
    }

    edge_df <- edge_df[keep, , drop = FALSE]
    m_from <- m_from[keep]
    m_to <- m_to[keep]
    weight <- edge_df$weight

    within <- m_from == m_to

    nodes <- names(membership_id)
    node_within_sum <- NULL
    node_within_n <- NULL
    node_between_sum <- NULL
    node_between_n <- NULL

    if (any(within)) {
        w_within <- weight[within]
        from_within <- edge_df$from[within]
        to_within <- edge_df$to[within]
        node_within_sum <- tapply(c(w_within, w_within),
            c(from_within, to_within), sum
        )
        node_within_n <- tapply(rep(1L, length(w_within) * 2L),
            c(from_within, to_within), sum
        )
    }

    if (any(!within)) {
        w_between <- weight[!within]
        from_between <- edge_df$from[!within]
        to_between <- edge_df$to[!within]
        node_between_sum <- tapply(c(w_between, w_between),
            c(from_between, to_between), sum
        )
        node_between_n <- tapply(rep(1L, length(w_between) * 2L),
            c(from_between, to_between), sum
        )
    }

    node_within_sum <- .mq_named_fill(node_within_sum, nodes, 0)
    node_within_n <- .mq_named_fill(node_within_n, nodes, 0L)
    node_between_sum <- .mq_named_fill(node_between_sum, nodes, 0)
    node_between_n <- .mq_named_fill(node_between_n, nodes, 0L)

    node_modules <- module_levels[membership_id[nodes]]

    nodes_df <- data.frame(
        gene = nodes,
        module = node_modules,
        module_id = as.integer(membership_id[nodes]),
        within_sum = as.numeric(node_within_sum),
        within_n = as.integer(node_within_n),
        between_sum = as.numeric(node_between_sum),
        between_n = as.integer(node_between_n),
        stringsAsFactors = FALSE
    )

    module_sizes <- table(membership_id)
    module_ids <- as.integer(names(module_sizes))
    module_labels <- module_levels[module_ids]

    mod_within_sum <- NULL
    mod_within_n <- NULL
    mod_between_sum <- NULL
    mod_between_n <- NULL

    if (any(within)) {
        mod_within_sum <- tapply(weight[within], m_from[within], sum)
        mod_within_n <- tapply(rep(1L, sum(within)), m_from[within], sum)
    }

    if (any(!within)) {
        mod_between_sum <- tapply(
            c(weight[!within], weight[!within]),
            c(m_from[!within], m_to[!within]),
            sum
        )
        mod_between_n <- tapply(
            rep(1L, sum(!within) * 2L),
            c(m_from[!within], m_to[!within]),
            sum
        )
    }

    mod_within_sum <- .mq_named_fill(mod_within_sum, module_ids, 0)
    mod_within_n <- .mq_named_fill(mod_within_n, module_ids, 0L)
    mod_between_sum <- .mq_named_fill(mod_between_sum, module_ids, 0)
    mod_between_n <- .mq_named_fill(mod_between_n, module_ids, 0L)

    modules_df <- data.frame(
        module = module_labels,
        module_id = module_ids,
        module_size = as.integer(module_sizes[as.character(module_ids)]),
        within_sum = as.numeric(mod_within_sum),
        within_n = as.integer(mod_within_n),
        between_sum = as.numeric(mod_between_sum),
        between_n = as.integer(mod_between_n),
        stringsAsFactors = FALSE
    )

    list(nodes = nodes_df, modules = modules_df)
}

#' Fill named vector with defaults for target names.
#' @keywords internal
.mq_named_fill <- function(vec, target_names, fill = 0) {
    out <- rep(fill, length(target_names))
    names(out) <- as.character(target_names)
    if (!is.null(vec) && length(vec)) {
        out[names(vec)] <- vec
    }
    out
}

#' Compute member-level scores.
#' @keywords internal
.mq_member_scores <- function(edge_df, membership_id, module_levels) {
    summaries <- .mq_edge_summaries(edge_df, membership_id, module_levels)
    nodes_df <- summaries$nodes
    if (!nrow(nodes_df)) return(nodes_df)

    within_mean <- ifelse(nodes_df$within_n > 0,
        nodes_df$within_sum / nodes_df$within_n,
        NA_real_
    )
    between_mean <- ifelse(nodes_df$between_n > 0,
        nodes_df$between_sum / nodes_df$between_n,
        NA_real_
    )

    nodes_df$member_score <- within_mean
    nodes_df$member_between <- between_mean
    nodes_df$member_margin <- within_mean - between_mean
    nodes_df
}

#' Compute module-level scores.
#' @keywords internal
.mq_module_scores <- function(edge_df, membership_id, module_levels, compute_conductance = TRUE) {
    summaries <- .mq_edge_summaries(edge_df, membership_id, module_levels)
    modules_df <- summaries$modules
    if (!nrow(modules_df)) return(modules_df)

    cohesion <- ifelse(modules_df$within_n > 0,
        modules_df$within_sum / modules_df$within_n,
        NA_real_
    )
    between_mean <- ifelse(modules_df$between_n > 0,
        modules_df$between_sum / modules_df$between_n,
        NA_real_
    )
    separation <- cohesion - between_mean
    conductance <- NA_real_
    if (isTRUE(compute_conductance)) {
        denom <- modules_df$within_sum + modules_df$between_sum
        conductance <- ifelse(denom > 0, modules_df$between_sum / denom, NA_real_)
    }

    modules_df$cohesion <- cohesion
    modules_df$between_mean <- between_mean
    modules_df$separation <- separation
    modules_df$conductance <- conductance
    modules_df
}

#' Compute overall within/between means.
#' @keywords internal
.mq_within_between_means <- function(edge_df, membership_id) {
    if (!nrow(edge_df)) return(list(within_mean = NA_real_, between_mean = NA_real_))
    m_from <- membership_id[edge_df$from]
    m_to <- membership_id[edge_df$to]
    keep <- !is.na(m_from) & !is.na(m_to) & is.finite(edge_df$weight)
    if (!any(keep)) return(list(within_mean = NA_real_, between_mean = NA_real_))
    m_from <- m_from[keep]
    m_to <- m_to[keep]
    weight <- edge_df$weight[keep]
    within <- m_from == m_to
    within_mean <- if (any(within)) mean(weight[within]) else NA_real_
    between_mean <- if (any(!within)) mean(weight[!within]) else NA_real_
    list(within_mean = within_mean, between_mean = between_mean)
}

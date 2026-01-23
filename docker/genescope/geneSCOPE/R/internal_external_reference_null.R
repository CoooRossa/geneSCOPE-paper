#' Degree-preserving null model for STRING evaluation.
#' @keywords internal
.external_reference_degree_null <- function(graph,
                                            edge_weight,
                                            gene_to_string,
                                            interactions,
                                            score_threshold,
                                            precision_k,
                                            precision_obs,
                                            enrichment_obs,
                                            set_defs,
                                            n_null = 50L,
                                            seed = 1,
                                            expected_positive_rate = NA_real_) {
    if (!igraph::is_igraph(graph) || n_null <= 0L) {
        return(NULL)
    }
    precision_k <- as.integer(precision_k)
    precision_k <- precision_k[is.finite(precision_k) & precision_k > 0L]
    n_prec <- length(precision_k)
    n_sets <- if (!is.null(set_defs)) nrow(set_defs) else 0L

    prec_mat <- if (n_prec) matrix(NA_real_, nrow = n_null, ncol = n_prec) else NULL
    enrich_mat <- if (n_sets) matrix(NA_real_, nrow = n_null, ncol = n_sets) else NULL
    if (n_prec) colnames(prec_mat) <- paste0("precision_at_", precision_k)
    if (n_sets) colnames(enrich_mat) <- set_defs$set_id

    set.seed(seed)
    for (i in seq_len(n_null)) {
        g_null <- igraph::rewire(
            graph,
            with = igraph::keeping_degseq(niter = max(10L, igraph::ecount(graph) * 10L))
        )
        edge_null <- .graph_edge_table(g_null, edge_weight = edge_weight)
        edge_null <- .stringdb_label_edges(
            edge_df = edge_null,
            gene_to_string = gene_to_string,
            interactions = interactions,
            score_threshold = score_threshold
        )
        if (n_prec) {
            prec_vals <- .external_reference_precision_at_k(edge_null$label, edge_null$score, k = precision_k)
            prec_mat[i, ] <- prec_vals[paste0("precision_at_", precision_k)]
        }
        if (n_sets) {
            enrich_vals <- .external_reference_top_enrichment(
                edge_null,
                score_threshold = score_threshold,
                expected_rate = expected_positive_rate,
                set_defs = set_defs
            )
            enrich_mat[i, ] <- enrich_vals$enrichment_ratio
        }
    }

    summary_rows <- list()
    if (n_prec) {
        for (k in precision_k) {
            metric <- paste0("precision_at_", k)
            obs <- precision_obs[[metric]]
            vals <- prec_mat[, metric]
            mu <- mean(vals, na.rm = TRUE)
            sd_val <- stats::sd(vals, na.rm = TRUE)
            z <- if (!is.na(obs) && is.finite(sd_val) && sd_val > 0) (obs - mu) / sd_val else NA_real_
            p <- if (!is.na(obs)) (1 + sum(vals >= obs, na.rm = TRUE)) / (length(vals) + 1) else NA_real_
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
                metric = metric,
                observed = obs,
                null_mean = mu,
                null_sd = sd_val,
                null_z = z,
                null_p = p,
                stringsAsFactors = FALSE
            )
        }
    }

    if (n_sets) {
        for (i in seq_len(n_sets)) {
            metric <- paste0("enrichment_", set_defs$set_id[i])
            obs <- NA_real_
            if (!is.null(enrichment_obs) && nrow(enrichment_obs)) {
                obs <- enrichment_obs$enrichment_ratio[match(set_defs$set_id[i], enrichment_obs$set_id)]
            }
            vals <- enrich_mat[, set_defs$set_id[i]]
            mu <- mean(vals, na.rm = TRUE)
            sd_val <- stats::sd(vals, na.rm = TRUE)
            z <- if (!is.na(obs) && is.finite(sd_val) && sd_val > 0) (obs - mu) / sd_val else NA_real_
            p <- if (!is.na(obs)) (1 + sum(vals >= obs, na.rm = TRUE)) / (length(vals) + 1) else NA_real_
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
                metric = metric,
                observed = obs,
                null_mean = mu,
                null_sd = sd_val,
                null_z = z,
                null_p = p,
                stringsAsFactors = FALSE
            )
        }
    }

    if (!length(summary_rows)) {
        return(NULL)
    }
    do.call(rbind, summary_rows)
}

#' Label-permutation null model for STRING evaluation.
#' @keywords internal
.external_reference_label_null <- function(graph,
                                           edge_weight,
                                           gene_to_string,
                                           interactions,
                                           score_threshold,
                                           precision_k,
                                           precision_obs,
                                           enrichment_obs,
                                           set_defs,
                                           n_null = 50L,
                                           seed = 1,
                                           expected_positive_rate = NA_real_) {
    if (!igraph::is_igraph(graph) || n_null <= 0L) {
        return(NULL)
    }
    precision_k <- as.integer(precision_k)
    precision_k <- precision_k[is.finite(precision_k) & precision_k > 0L]
    n_prec <- length(precision_k)
    n_sets <- if (!is.null(set_defs)) nrow(set_defs) else 0L

    prec_mat <- if (n_prec) matrix(NA_real_, nrow = n_null, ncol = n_prec) else NULL
    enrich_mat <- if (n_sets) matrix(NA_real_, nrow = n_null, ncol = n_sets) else NULL
    if (n_prec) colnames(prec_mat) <- paste0("precision_at_", precision_k)
    if (n_sets) colnames(enrich_mat) <- set_defs$set_id

    edge_df <- .graph_edge_table(graph, edge_weight = edge_weight)
    set.seed(seed)
    for (i in seq_len(n_null)) {
        mapping_perm <- gene_to_string
        mapping_perm <- setNames(sample(mapping_perm, length(mapping_perm), replace = FALSE), names(mapping_perm))
        edge_null <- .stringdb_label_edges(
            edge_df = edge_df,
            gene_to_string = mapping_perm,
            interactions = interactions,
            score_threshold = score_threshold
        )
        if (n_prec) {
            prec_vals <- .external_reference_precision_at_k(edge_null$label, edge_null$score, k = precision_k)
            prec_mat[i, ] <- prec_vals[paste0("precision_at_", precision_k)]
        }
        if (n_sets) {
            enrich_vals <- .external_reference_top_enrichment(
                edge_null,
                score_threshold = score_threshold,
                expected_rate = expected_positive_rate,
                set_defs = set_defs
            )
            enrich_mat[i, ] <- enrich_vals$enrichment_ratio
        }
    }

    summary_rows <- list()
    if (n_prec) {
        for (k in precision_k) {
            metric <- paste0("precision_at_", k)
            obs <- precision_obs[[metric]]
            vals <- prec_mat[, metric]
            mu <- mean(vals, na.rm = TRUE)
            sd_val <- stats::sd(vals, na.rm = TRUE)
            z <- if (!is.na(obs) && is.finite(sd_val) && sd_val > 0) (obs - mu) / sd_val else NA_real_
            p <- if (!is.na(obs)) (1 + sum(vals >= obs, na.rm = TRUE)) / (length(vals) + 1) else NA_real_
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
                metric = metric,
                observed = obs,
                null_mean = mu,
                null_sd = sd_val,
                null_z = z,
                null_p = p,
                stringsAsFactors = FALSE
            )
        }
    }

    if (n_sets) {
        for (i in seq_len(n_sets)) {
            metric <- paste0("enrichment_", set_defs$set_id[i])
            obs <- NA_real_
            if (!is.null(enrichment_obs) && nrow(enrichment_obs)) {
                obs <- enrichment_obs$enrichment_ratio[match(set_defs$set_id[i], enrichment_obs$set_id)]
            }
            vals <- enrich_mat[, set_defs$set_id[i]]
            mu <- mean(vals, na.rm = TRUE)
            sd_val <- stats::sd(vals, na.rm = TRUE)
            z <- if (!is.na(obs) && is.finite(sd_val) && sd_val > 0) (obs - mu) / sd_val else NA_real_
            p <- if (!is.na(obs)) (1 + sum(vals >= obs, na.rm = TRUE)) / (length(vals) + 1) else NA_real_
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
                metric = metric,
                observed = obs,
                null_mean = mu,
                null_sd = sd_val,
                null_z = z,
                null_p = p,
                stringsAsFactors = FALSE
            )
        }
    }

    if (!length(summary_rows)) {
        return(NULL)
    }
    do.call(rbind, summary_rows)
}

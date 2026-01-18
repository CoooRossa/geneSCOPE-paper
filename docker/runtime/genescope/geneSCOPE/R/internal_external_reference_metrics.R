#' Compute precision at K for a score ranking.
#' @keywords internal
.external_reference_precision_at_k <- function(labels, scores, k = c(50L, 100L, 500L)) {
    k <- as.integer(k)
    k <- k[is.finite(k) & k > 0]
    if (!length(k)) return(setNames(numeric(0), character(0)))
    vals <- vapply(k, function(kv) .precision_at_k(labels, scores, k = kv), numeric(1))
    names(vals) <- paste0("precision_at_", k)
    vals
}

#' Compute enrichment of STRING-positive edges vs background.
#' @keywords internal
.external_reference_enrichment <- function(n_pos_graph, n_graph, n_pos_bg, n_pairs_bg) {
    expected_rate <- if (!is.na(n_pairs_bg) && n_pairs_bg > 0) n_pos_bg / n_pairs_bg else NA_real_
    ratio <- NA_real_
    if (!is.na(expected_rate) && expected_rate > 0 && n_graph > 0) {
        ratio <- (n_pos_graph / n_graph) / expected_rate
    }

    pval <- NA_real_
    if (all(!is.na(c(n_pos_graph, n_graph, n_pos_bg, n_pairs_bg)))) {
        if (n_pairs_bg >= n_graph && n_pos_bg >= n_pos_graph) {
            n_non_graph <- n_pairs_bg - n_graph
            n_pos_non_graph <- n_pos_bg - n_pos_graph
            n_neg_graph <- n_graph - n_pos_graph
            n_neg_non_graph <- n_non_graph - n_pos_non_graph
            if (min(n_non_graph, n_pos_non_graph, n_neg_graph, n_neg_non_graph) >= 0) {
                tab <- matrix(c(n_pos_graph, n_neg_graph, n_pos_non_graph, n_neg_non_graph), nrow = 2)
                pval <- tryCatch(
                    stats::fisher.test(tab, alternative = "greater")$p.value,
                    error = function(e) NA_real_
                )
            }
        }
    }

    list(
        expected_positive_rate = expected_rate,
        enrichment_ratio = ratio,
        enrichment_p = pval
    )
}

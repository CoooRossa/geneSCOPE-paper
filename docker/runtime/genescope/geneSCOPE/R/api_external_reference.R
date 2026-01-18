#' Evaluate external reference concordance using STRING.
#' @description
#' Ranks graph edges by their weights and evaluates agreement with STRING
#' protein-protein associations using AUPRC (primary), AUROC, and Top-K precision.
#' @param x A `scope_object`, an `igraph`, or a list returned by `plotDendroNetwork()`.
#' @param grid_name Grid layer name for `scope_object` input (default `"grid30"`).
#' @param stats_layer Stats layer name for `scope_object` input (default `"LeeStats_Xz"`).
#' @param graph_slot Graph slot name for `scope_object` input (default `"g_consensus"`).
#' @param edge_weight Edge attribute used for ranking edges (default `"weight"`).
#' @param species NCBI taxonomy ID for STRING (default 9606).
#' @param string_score_threshold STRING combined score threshold for positives.
#' @param cache_dir Optional directory for caching STRING mappings/interactions.
#' @param keep_details Whether to retain edge-level annotations in the output.
#' @param n_null Number of degree-preserving rewires for null AUPRC (optional).
#' @param seed Random seed for null generation.
#' @param verbose Emit progress messages when TRUE.
#' @param ... Reserved for future extensions.
#' @return A list with `summary`, optional `details`, and `notes`.
#' @examples
#' \dontrun{
#' edge_eval <- evaluateExternalReference_STRING(scope_obj, grid_name = "grid30")
#' }
#' @export
evaluateExternalReference_STRING <- function(
    x,
    grid_name = "grid30",
    stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    edge_weight = "weight",
    species = 9606,
    string_score_threshold = 700,
    cache_dir = NULL,
    keep_details = FALSE,
    n_null = 0,
    seed = 1,
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...) {
    parent <- "evaluateExternalReference_STRING"
    step01 <- .log_step(parent, "S01", "resolve graph and edge table", verbose)
    step01$enter()
    g <- .extract_graph_any(
        x = x,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot
    )
    if (is.null(igraph::V(g)$name)) stop("Graph vertices missing names (V(g)$name).")
    edge_df <- .graph_edge_table(g, edge_weight = edge_weight)
    step01$done(paste0("edges=", nrow(edge_df)))

    step02 <- .log_step(parent, "S02", "map genes and load STRING interactions", verbose)
    step02$enter()
    string_db <- .stringdb_connect(
        species = species,
        score_threshold = 0,
        cache_dir = cache_dir
    )
    genes <- igraph::V(g)$name
    mapping <- .stringdb_map_genes(string_db, genes, cache_dir = cache_dir)
    interactions <- .stringdb_get_interactions(string_db, mapping, cache_dir = cache_dir)
    step02$done(paste0("mapped=", sum(!is.na(mapping)), "/", length(mapping)))

    step03 <- .log_step(parent, "S03", "score edges against STRING", verbose)
    step03$enter()
    edge_df <- .stringdb_label_edges(
        edge_df = edge_df,
        gene_to_string = mapping,
        interactions = interactions,
        score_threshold = string_score_threshold
    )
    step03$done()

    mapping_rate <- if (length(mapping)) mean(!is.na(mapping)) else NA_real_
    comparable_edge_rate <- if (nrow(edge_df)) mean(!is.na(edge_df$string_score)) else NA_real_
    positive_rate <- NA_real_
    if (nrow(edge_df)) {
        comp <- !is.na(edge_df$string_score)
        if (any(comp)) {
            positive_rate <- mean(edge_df$string_score[comp] >= string_score_threshold)
        }
    }

    auprc <- .compute_auprc(edge_df$label, edge_df$score)
    auroc <- .compute_auroc(edge_df$label, edge_df$score)
    precision_at_100 <- .precision_at_k(edge_df$label, edge_df$score, k = 100L)
    precision_at_500 <- .precision_at_k(edge_df$label, edge_df$score, k = 500L)
    precision_at_1000 <- .precision_at_k(edge_df$label, edge_df$score, k = 1000L)

    null_mean <- NA_real_
    null_sd <- NA_real_
    delta_auprc <- NA_real_
    empirical_p <- NA_real_
    if (n_null > 0 && igraph::ecount(g) > 0) {
        set.seed(seed)
        null_vals <- numeric(n_null)
        weight_vals <- igraph::edge_attr(g, edge_weight)
        if (is.null(weight_vals)) weight_vals <- igraph::edge_attr(g, "weight")
        if (is.null(weight_vals) || !length(weight_vals)) {
            weight_vals <- rep(1, igraph::ecount(g))
        }
        for (i in seq_len(n_null)) {
            g_null <- igraph::rewire(
                g,
                with = igraph::keeping_degseq(niter = max(10L, igraph::ecount(g) * 10L))
            )
            weight_sample <- sample(weight_vals, igraph::ecount(g_null), replace = FALSE)
            igraph::E(g_null)$weight <- weight_sample
            if (!identical(edge_weight, "weight")) {
                igraph::E(g_null)[[edge_weight]] <- weight_sample
            }
            edge_null <- .graph_edge_table(g_null, edge_weight = edge_weight)
            edge_null <- .stringdb_label_edges(
                edge_df = edge_null,
                gene_to_string = mapping,
                interactions = interactions,
                score_threshold = string_score_threshold
            )
            null_vals[i] <- .compute_auprc(edge_null$label, edge_null$score)
        }
        null_mean <- mean(null_vals, na.rm = TRUE)
        null_sd <- sd(null_vals, na.rm = TRUE)
        if (is.nan(null_mean)) null_mean <- NA_real_
        if (is.nan(null_sd)) null_sd <- NA_real_
        delta_auprc <- if (!is.na(auprc) && !is.na(null_mean)) auprc - null_mean else NA_real_
        if (!is.na(auprc)) {
            empirical_p <- (1 + sum(null_vals >= auprc, na.rm = TRUE)) / (n_null + 1)
        }
    }

    summary <- data.frame(
        mapping_rate = mapping_rate,
        comparable_edge_rate = comparable_edge_rate,
        positive_rate = positive_rate,
        AUPRC = auprc,
        AUROC = auroc,
        precision_at_100 = precision_at_100,
        precision_at_500 = precision_at_500,
        precision_at_1000 = precision_at_1000,
        null_mean_AUPRC = null_mean,
        null_sd_AUPRC = null_sd,
        delta_AUPRC = delta_auprc,
        empirical_p = empirical_p,
        stringsAsFactors = FALSE
    )

    notes <- paste0(
        "AUPRC provides an objective, reproducible measure of how well the graph",
        " ranks gene-gene edges that agree with STRING (score >= ", string_score_threshold, ")."
    )
    out <- list(summary = summary, notes = notes)
    if (isTRUE(keep_details)) out$details <- edge_df
    out
}

#' Evaluate external reference concordance using STRING (CamelCase API).
#' @description
#' Scores gene-gene edges against STRING and reports AUROC, AUPRC, precision@K,
#' enrichment statistics, and coverage metrics.
#' @param x A `scope_object`, an `igraph`, or a list returned by `plotDendroNetwork()`.
#' @param grid_name Grid layer name for `scope_object` input (default `"grid30"`).
#' @param stats_layer Stats layer name for `scope_object` input (default `"LeeStats_Xz"`).
#' @param graph_slot_name Graph slot name for `scope_object` input (default `"g_consensus"`).
#' @param graph_slot_mh Optional MH graph slot name (default `"g_morisita"`).
#' @param edge_weight Edge attribute used for ranking edges (default `"weight"`).
#' @param species NCBI taxonomy ID for STRING (default 9606).
#' @param string_score_threshold STRING combined score threshold for positives.
#' @param precision_k Vector of K values for precision@K.
#' @param cache_dir Optional directory for caching STRING mappings/interactions.
#' @param keep_details Whether to retain edge-level annotations in the output.
#' @param verbose Emit progress messages when TRUE.
#' @param ... Reserved for future extensions.
#' @return A list with `metrics`, `coverage`, `params`, and `notes`.
#' @examples
#' \dontrun{
#' score <- evaluateExternalReferenceSTRING(scope_obj, grid_name = "grid30")
#' }
#' @export
evaluateExternalReferenceSTRING <- function(
    x,
    grid_name = "grid30",
    stats_layer = "LeeStats_Xz",
    graph_slot_name = "g_consensus",
    graph_slot_mh = "g_morisita",
    edge_weight = "weight",
    species = 9606,
    string_score_threshold = 700,
    precision_k = c(50L, 100L, 500L),
    cache_dir = NULL,
    keep_details = FALSE,
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...) {
    parent <- "evaluateExternalReferenceSTRING"
    step01 <- .log_step(parent, "S01", "resolve graph and edge table", verbose)
    step01$enter()
    graph_slot <- graph_slot_name
    if (is.null(graph_slot) && !is.null(graph_slot_mh)) graph_slot <- graph_slot_mh
    g <- .normalize_input_to_igraph(
        x = x,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot
    )
    if (is.null(igraph::V(g)$name)) stop("Graph vertices missing names (V(g)$name).")
    edge_df <- .graph_edge_table(g, edge_weight = edge_weight)
    step01$done(paste0("edges=", nrow(edge_df)))

    step02 <- .log_step(parent, "S02", "map genes and load STRING interactions", verbose)
    step02$enter()
    string_db <- .stringdb_connect(
        species = species,
        score_threshold = 0,
        cache_dir = cache_dir
    )
    genes <- igraph::V(g)$name
    mapping <- .stringdb_map_genes(string_db, genes, cache_dir = cache_dir)
    interactions <- .stringdb_get_interactions(string_db, mapping, cache_dir = cache_dir)
    step02$done(paste0("mapped=", sum(!is.na(mapping)), "/", length(mapping)))

    step03 <- .log_step(parent, "S03", "score edges against STRING", verbose)
    step03$enter()
    edge_df <- .stringdb_label_edges(
        edge_df = edge_df,
        gene_to_string = mapping,
        interactions = interactions,
        score_threshold = string_score_threshold
    )
    step03$done()

    n_genes <- igraph::vcount(g)
    n_edges <- nrow(edge_df)
    n_mapped_genes <- sum(!is.na(mapping))
    mapping_rate <- if (length(mapping)) mean(!is.na(mapping)) else NA_real_
    n_comparable_edges <- sum(!is.na(edge_df$string_score))
    comparable_edge_rate <- if (n_edges) mean(!is.na(edge_df$string_score)) else NA_real_
    positive_rate <- NA_real_
    if (n_edges && n_comparable_edges) {
        positive_rate <- mean(edge_df$string_score[!is.na(edge_df$string_score)] >= string_score_threshold)
    }

    auprc <- .compute_auprc(edge_df$label, edge_df$score)
    auroc <- .compute_auroc(edge_df$label, edge_df$score)
    prec_vals <- .external_reference_precision_at_k(edge_df$label, edge_df$score, k = precision_k)

    mapped_ids <- unique(na.omit(as.character(mapping)))
    n_pairs_bg <- length(mapped_ids) * (length(mapped_ids) - 1) / 2
    n_pos_bg <- 0L
    if (nrow(interactions)) {
        pos_idx <- interactions$combined_score >= string_score_threshold
        if (any(pos_idx)) {
            keys <- paste(
                pmin(interactions$from[pos_idx], interactions$to[pos_idx]),
                pmax(interactions$from[pos_idx], interactions$to[pos_idx]),
                sep = "|"
            )
            n_pos_bg <- length(unique(keys))
        }
    }
    n_pos_graph <- sum(edge_df$string_score >= string_score_threshold, na.rm = TRUE)
    enrich <- .external_reference_enrichment(
        n_pos_graph = n_pos_graph,
        n_graph = n_comparable_edges,
        n_pos_bg = n_pos_bg,
        n_pairs_bg = n_pairs_bg
    )

    metrics <- data.frame(
        AUROC = auroc,
        AUPRC = auprc,
        enrichment_ratio = enrich$enrichment_ratio,
        enrichment_p = enrich$enrichment_p,
        expected_positive_rate = enrich$expected_positive_rate,
        stringsAsFactors = FALSE
    )
    if (length(prec_vals)) {
        metrics <- cbind(metrics, as.data.frame(as.list(prec_vals)))
    }

    coverage <- data.frame(
        n_genes = n_genes,
        n_edges = n_edges,
        n_mapped_genes = n_mapped_genes,
        mapping_rate = mapping_rate,
        n_comparable_edges = n_comparable_edges,
        comparable_edge_rate = comparable_edge_rate,
        positive_rate = positive_rate,
        stringsAsFactors = FALSE
    )

    params <- list(
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot,
        edge_weight = edge_weight,
        species = species,
        string_score_threshold = string_score_threshold,
        precision_k = as.integer(precision_k)
    )

    notes <- paste0(
        "STRING-based scoring uses known interactions as positives; coverage and scores can reflect",
        " study biases in public databases. Interpret AUROC/AUPRC with this bias in mind."
    )

    out <- list(metrics = metrics, coverage = coverage, params = params, notes = notes)
    if (isTRUE(keep_details)) out$details <- edge_df
    out
}

#' Evaluate module quality with STRING and graph metrics.
#' @description
#' Computes within- vs between-module STRING score differences and a classical
#' graph metric (weighted modularity) for cluster/module assignments.
#' @param x A `scope_object`, an `igraph`, or a list returned by `plotDendroNetwork()`.
#' @param membership Optional membership vector (named by gene).
#' @param grid_name Grid layer name for `scope_object` input (default `"grid30"`).
#' @param stats_layer Stats layer name for `scope_object` input (default `"LeeStats_Xz"`).
#' @param graph_slot Graph slot name for `scope_object` input (default `"g_consensus"`).
#' @param edge_weight Edge attribute used for weighting (default `"weight"`).
#' @param species NCBI taxonomy ID for STRING (default 9606).
#' @param string_score_threshold STRING combined score threshold for positives.
#' @param n_perm Number of permutations for within/between delta.
#' @param perm_mode Permutation mode (`label_shuffle` only).
#' @param seed Random seed for permutations.
#' @param verbose Emit progress messages when TRUE.
#' @param ... Reserved for future extensions.
#' @return A list with `summary`, `per_module`, and `notes`.
#' @examples
#' \dontrun{
#' mod_eval <- evaluateModuleQuality_STRING(scope_obj, grid_name = "grid30")
#' }
#' @export
evaluateModuleQuality_STRING <- function(
    x,
    membership = NULL,
    grid_name = "grid30",
    stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    edge_weight = "weight",
    species = 9606,
    string_score_threshold = 700,
    n_perm = 200,
    perm_mode = "label_shuffle",
    seed = 1,
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...) {
    parent <- "evaluateModuleQuality_STRING"
    step01 <- .log_step(parent, "S01", "resolve graph and membership", verbose)
    step01$enter()
    g <- .extract_graph_any(
        x = x,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot
    )
    m <- .extract_membership_any(
        x = x,
        graph = g,
        membership = membership,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot
    )
    step01$done(paste0("modules=", length(unique(m))))

    step02 <- .log_step(parent, "S02", "map genes and load STRING interactions", verbose)
    step02$enter()
    string_db <- .stringdb_connect(species = species, score_threshold = 0)
    genes <- igraph::V(g)$name
    mapping <- .stringdb_map_genes(string_db, genes, cache_dir = NULL)
    interactions <- .stringdb_get_interactions(string_db, mapping, cache_dir = NULL)
    step02$done(paste0("mapped=", sum(!is.na(mapping)), "/", length(mapping)))

    edge_df <- .graph_edge_table(g, edge_weight = edge_weight)
    edge_df <- .stringdb_label_edges(
        edge_df = edge_df,
        gene_to_string = mapping,
        interactions = interactions,
        score_threshold = string_score_threshold
    )

    module_stats <- .module_within_between(edge_df, m, score_col = "string_score")
    weighted_modularity <- .compute_weighted_modularity(g, m, edge_weight = edge_weight)

    perm_p <- NA_real_
    if (n_perm > 0 && !is.na(module_stats$delta)) {
        if (!identical(perm_mode, "label_shuffle")) {
            stop("perm_mode must be 'label_shuffle'.")
        }
        set.seed(seed)
        delta_perm <- numeric(n_perm)
        for (i in seq_len(n_perm)) {
            perm_m <- sample(m)
            names(perm_m) <- names(m)
            delta_perm[i] <- .module_within_between(edge_df, perm_m, score_col = "string_score")$delta
        }
        if (all(is.na(delta_perm))) {
            perm_p <- NA_real_
        } else {
            perm_p <- (1 + sum(abs(delta_perm) >= abs(module_stats$delta), na.rm = TRUE)) / (n_perm + 1)
        }
    }

    summary <- data.frame(
        weighted_modularity = weighted_modularity,
        within_string_mean = module_stats$within_mean,
        between_string_mean = module_stats$between_mean,
        delta_within_between = module_stats$delta,
        perm_p = perm_p,
        stringsAsFactors = FALSE
    )

    notes <- paste0(
        "Within-vs-between STRING score deltas provide a quantitative, external",
        " check of module coherence (STRING score >= ", string_score_threshold, ")."
    )

    list(summary = summary, per_module = module_stats$per_module, notes = notes)
}

#' Summarize edge and module evaluations.
#' @param edge_eval Output from `evaluateExternalReference_STRING()`.
#' @param module_eval Output from `evaluateModuleQuality_STRING()`.
#' @return A list with `edge_level_table` and `module_level_table`.
#' @examples
#' \dontrun{
#' summary_tables <- summarizeEvaluation(edge_eval, mod_eval)
#' }
#' @export
summarizeEvaluation <- function(edge_eval, module_eval) {
    if (is.null(edge_eval$summary)) stop("edge_eval must include $summary.")
    if (is.null(module_eval$summary)) stop("module_eval must include $summary.")
    edge_level_table <- edge_eval$summary

    mod_summary <- module_eval$summary
    per_module <- module_eval$per_module
    if (is.null(per_module)) per_module <- data.frame()
    if (nrow(per_module)) {
        if (!("between_string_mean" %in% names(per_module))) {
            per_module$between_string_mean <- mod_summary$between_string_mean
        }
        if (!("delta_within_between" %in% names(per_module))) {
            per_module$delta_within_between <- mod_summary$delta_within_between
        }
        per_module$perm_p <- mod_summary$perm_p
        per_module$weighted_modularity <- mod_summary$weighted_modularity
    }
    list(edge_level_table = edge_level_table, module_level_table = per_module)
}

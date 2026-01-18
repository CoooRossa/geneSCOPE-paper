#' Score external reference concordance using STRING and store results.
#' @description
#' Computes STRING agreement metrics for a graph and, when a `scope_object`
#' is supplied, writes a structured payload back to
#' `scope_obj@stats[[grid_name]][[stats_layer]]` under `external_reference`.
#' @param x A `scope_object`, an `igraph`, or a list returned by `plotDendroNetwork()`.
#' @param grid_name Grid layer name for `scope_object` input (default `"grid30"`).
#' @param stats_layer Stats layer name for `scope_object` input (default `"LeeStats_Xz"`).
#' @param graph_slot Graph slot name for `scope_object` input (default `"g_consensus"`).
#' @param edge_weight Edge attribute used for ranking edges (default `"weight"`).
#' @param species NCBI taxonomy ID for STRING (default 9606).
#' @param input_id_type STRINGdb input ID type (default `"gene"`).
#' @param string_score_threshold STRING combined score threshold for positives.
#' @param precision_k Vector of K values for precision@K.
#' @param top_frac Fraction of top edges to evaluate enrichment (default 0.01).
#' @param null_method Null strategy (`degree_preserving` or `label_permutation`).
#' @param run_null Whether to run null evaluation.
#' @param n_null Number of null rewires (default 50).
#' @param null_seed Random seed for null evaluation.
#' @param cache_dir Optional directory for caching STRING mappings/interactions.
#' @param keep_details Whether to retain edge-level annotations in the output.
#' @param verbose Emit progress messages when TRUE.
#' @param ... Reserved for future extensions.
#' @return If `x` is a `scope_object`, returns the updated `scope_object`.
#' Otherwise returns the scoring payload list.
#' @examples
#' \dontrun{
#' scope_obj <- scoreExternalReferenceSTRING(scope_obj, grid_name = "grid30")
#' }
#' @export
scoreExternalReferenceSTRING <- function(
    x,
    grid_name = "grid30",
    stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    edge_weight = "weight",
    species = 9606,
    input_id_type = "gene",
    string_score_threshold = 700,
    precision_k = c(50L, 100L, 500L),
    top_frac = 0.01,
    null_method = c("degree_preserving", "label_permutation"),
    run_null = FALSE,
    n_null = 50L,
    null_seed = 1,
    cache_dir = NULL,
    keep_details = FALSE,
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...) {
    res <- .external_reference_score_stringdb(
        x = x,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot,
        edge_weight = edge_weight,
        species = species,
        input_id_type = input_id_type,
        string_score_threshold = string_score_threshold,
        precision_k = precision_k,
        top_frac = top_frac,
        null_method = null_method,
        run_null = run_null,
        n_null = n_null,
        null_seed = null_seed,
        cache_dir = cache_dir,
        keep_details = keep_details,
        verbose = verbose
    )
    payload <- res$payload

    if (inherits(x, "scope_object")) {
        scope_obj <- res$cache_scope
        if (is.null(scope_obj)) scope_obj <- x
        scope_obj <- .external_ref_writeback(
            scope_obj = scope_obj,
            grid_name = grid_name,
            stats_layer = stats_layer,
            graph_slot = graph_slot,
            payload = payload
        )
        return(scope_obj)
    }

    payload
}

#' Extract stored external reference scores.
#' @description
#' Reads the STRING evaluation payload from `scope_obj@stats` and returns a
#' flattened data.frame for reporting.
#' @param x A `scope_object` or a list containing `external_reference`.
#' @param grid_name Grid layer name for `scope_object` input.
#' @param stats_layer Stats layer name for `scope_object` input.
#' @param graph_slot Graph slot name for `scope_object` input.
#' @param ref External reference name (default `"STRINGdb"`).
#' @return A data.frame with one row per graph (or empty if not found).
#' @export
getExternalReferenceScores <- function(
    x,
    grid_name = "grid30",
    stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    ref = "STRINGdb") {
    payload <- NULL
    if (inherits(x, "scope_object")) {
        payload <- .external_ref_read(
            scope_obj = x,
            grid_name = grid_name,
            stats_layer = stats_layer,
            graph_slot = graph_slot,
            ref = ref
        )
    } else if (is.list(x) && !is.null(x$external_reference) && !is.null(x$external_reference[[ref]])) {
        payload <- x$external_reference[[ref]]
    }

    if (is.null(payload)) return(data.frame())
    .external_reference_scores_to_df(payload, grid_name, stats_layer, graph_slot)
}

#' Summarize stored external reference scores.
#' @description
#' Scans a `scope_object` for stored STRING evaluation payloads and returns
#' per-graph and optional per-module summary tables.
#' @param x A `scope_object` containing stored external reference scores.
#' @param grid_name Optional grid filter.
#' @param stats_layer Optional stats layer filter.
#' @param ref External reference name (default `"STRINGdb"`).
#' @param membership Optional membership vector for per-module summaries.
#' @param wide Whether to return wide (TRUE) or long (FALSE) per-graph table.
#' @return A list with `per_graph`, `per_module`, and `summary_text`.
#' @export
summarizeExternalReferenceScores <- function(
    x,
    grid_name = NULL,
    stats_layer = NULL,
    ref = "STRINGdb",
    membership = NULL,
    wide = TRUE) {
    if (!inherits(x, "scope_object")) {
        stop("summarizeExternalReferenceScores requires a scope_object.")
    }

    stats_root <- x@stats
    entries <- list()
    per_module_rows <- list()

    for (gname in names(stats_root)) {
        if (!is.null(grid_name) && !identical(gname, grid_name)) next
        layer_root <- stats_root[[gname]]
        for (lname in names(layer_root)) {
            if (!is.null(stats_layer) && !identical(lname, stats_layer)) next
            layer_obj <- layer_root[[lname]]

            if (is.list(layer_obj$external_reference) && !is.null(layer_obj$external_reference[[ref]])) {
                ref_obj <- layer_obj$external_reference[[ref]]
                for (slot in names(ref_obj)) {
                    payload <- ref_obj[[slot]]
                    entries[[length(entries) + 1]] <- .external_reference_scores_to_df(payload, gname, lname, slot)
                }
            }

            for (slot in names(layer_obj)) {
                graph_obj <- layer_obj[[slot]]
                if (is.list(graph_obj) && !igraph::is_igraph(graph_obj)) {
                    if (!is.null(graph_obj$external_reference) && !is.null(graph_obj$external_reference[[ref]])) {
                        payload <- graph_obj$external_reference[[ref]]
                        entries[[length(entries) + 1]] <- .external_reference_scores_to_df(payload, gname, lname, slot)
                    }
                }
            }
        }
    }

    per_graph <- if (length(entries)) do.call(rbind, entries) else data.frame()

    if (!wide && nrow(per_graph)) {
        id_cols <- c("grid_name", "stats_layer", "graph_slot")
        keep <- intersect(id_cols, names(per_graph))
        value_cols <- setdiff(names(per_graph), keep)
        long_rows <- lapply(value_cols, function(col) {
            data.frame(
                per_graph[keep],
                metric = col,
                value = per_graph[[col]],
                stringsAsFactors = FALSE
            )
        })
        per_graph <- do.call(rbind, long_rows)
    }

    if (!is.null(membership) && length(membership)) {
        for (gname in names(stats_root)) {
            if (!is.null(grid_name) && !identical(gname, grid_name)) next
            layer_root <- stats_root[[gname]]
            for (lname in names(layer_root)) {
                if (!is.null(stats_layer) && !identical(lname, stats_layer)) next
                layer_obj <- layer_root[[lname]]
                for (slot in names(layer_obj)) {
                    payload <- .external_ref_read(
                        scope_obj = x,
                        grid_name = gname,
                        stats_layer = lname,
                        graph_slot = slot,
                        ref = ref
                    )
                    if (is.null(payload) || is.null(payload$details)) next
                    graph_obj <- tryCatch(.extract_graph_any(
                        x = x,
                        grid_name = gname,
                        stats_layer = lname,
                        graph_slot = slot
                    ), error = function(e) NULL)
                    if (is.null(graph_obj)) next
                    m <- tryCatch(.extract_membership_any(
                        x = x,
                        graph = graph_obj,
                        membership = membership,
                        grid_name = gname,
                        stats_layer = lname,
                        graph_slot = slot
                    ), error = function(e) NULL)
                    if (is.null(m)) next
                    module_stats <- .module_within_between(payload$details, m, score_col = "string_score")
                    per_module <- module_stats$per_module
                    if (nrow(per_module)) {
                        per_module$grid_name <- gname
                        per_module$stats_layer <- lname
                        per_module$graph_slot <- slot
                        if (!("between_string_mean" %in% names(per_module))) {
                            per_module$between_string_mean <- module_stats$between_mean
                        }
                        if (!("delta_within_between" %in% names(per_module))) {
                            per_module$delta_within_between <- module_stats$delta
                        }
                        per_module_rows[[length(per_module_rows) + 1]] <- per_module
                    }
                }
            }
        }
    }

    per_module <- if (length(per_module_rows)) do.call(rbind, per_module_rows) else data.frame()
    summary_text <- paste0(
        "Summarized STRING external reference scores for ",
        nrow(per_graph),
        " graph(s)."
    )
    list(per_graph = per_graph, per_module = per_module, summary_text = summary_text)
}

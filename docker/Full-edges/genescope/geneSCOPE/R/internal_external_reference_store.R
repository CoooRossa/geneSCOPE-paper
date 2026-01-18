# Package-level cache environment (created at load to avoid locked namespace).
.geneSCOPE_cache_env <- new.env(parent = emptyenv())

#' Hash helper for external reference caches.
#' @keywords internal
.external_reference_hash <- function(x) {
    if (!requireNamespace("digest", quietly = TRUE)) {
        return(NA_character_)
    }
    digest::digest(x, algo = "sha256")
}

#' Access package-level cache environment for STRINGdb.
#' @keywords internal
.stringdb_cache_env <- function() {
    ns <- asNamespace("geneSCOPE")
    env <- get0(".geneSCOPE_cache_env", envir = ns, inherits = FALSE)
    if (is.null(env) || !is.environment(env)) {
        env <- new.env(parent = emptyenv())
    }
    env
}

#' Build STRINGdb interaction URL and local file path.
#' @keywords internal
.stringdb_interaction_paths <- function(string_db) {
    network_type_param <- ""
    if (tolower(string_db$network_type) == "physical") {
        network_type_param <- "physical."
    }
    link_data_param <- "links.v"
    link_data <- tolower(string_db$link_data)
    if (link_data == "detailed") {
        link_data_param <- "links.detailed.v"
    } else if (link_data == "full") {
        link_data_param <- "links.full.v"
    }
    file_version <- string_db$file_version
    species <- string_db$species
    file_base <- paste0(species, ".protein.", network_type_param, link_data_param, file_version, ".txt.gz")
    url <- paste0(
        string_db$protocol,
        "://stringdb-downloads.org/download/protein.",
        network_type_param,
        link_data_param,
        file_version,
        "/",
        file_base
    )
    list(
        url = url,
        file_base = file_base,
        file_path = file.path(string_db$input_directory, file_base)
    )
}

#' Fallback loader for STRING interactions that skips NA rows.
#' @keywords internal
.stringdb_get_interactions_fallback <- function(string_db, ids) {
    ids <- unique(na.omit(as.character(ids)))
    if (!length(ids)) {
        return(data.frame(from = character(0), to = character(0), combined_score = numeric(0)))
    }

    if (!dir.exists(string_db$input_directory)) {
        dir.create(string_db$input_directory, recursive = TRUE)
    }
    paths <- .stringdb_interaction_paths(string_db)
    if (file.exists(paths$file_path)) {
        unlink(paths$file_path)
    }

    download_fun <- getFromNamespace("downloadAbsentFile", "STRINGdb")
    temp <- download_fun(paths$url, oD = string_db$input_directory)
    ppi <- utils::read.table(temp, sep = " ", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

    required <- c("protein1", "protein2", "combined_score")
    if (!all(required %in% names(ppi))) {
        stop("STRINGdb interaction file missing required columns.")
    }
    ppi <- ppi[, required, drop = FALSE]
    keep <- stats::complete.cases(ppi)
    if (any(!keep)) {
        ppi <- ppi[keep, , drop = FALSE]
    }
    ppi <- ppi[ppi$protein1 %in% ids & ppi$protein2 %in% ids, , drop = FALSE]
    if (!nrow(ppi)) {
        return(data.frame(from = character(0), to = character(0), combined_score = numeric(0)))
    }
    out <- data.frame(
        from = ppi$protein1,
        to = ppi$protein2,
        combined_score = as.numeric(ppi$combined_score),
        stringsAsFactors = FALSE
    )
    out
}

#' Build digest of scoring parameters.
#' @keywords internal
.external_reference_params_digest <- function(params) {
    .external_reference_hash(params)
}

#' Define top-edge sets for enrichment evaluation.
#' @keywords internal
.external_reference_top_sets <- function(n_edges, precision_k, top_frac = NULL, n_comparable = NULL) {
    sets <- list()
    n_comp <- n_edges
    if (!is.null(n_comparable) && is.finite(n_comparable)) {
        n_comp <- n_comparable
    }
    precision_k <- as.integer(precision_k)
    precision_k <- precision_k[is.finite(precision_k) & precision_k > 0L]
    if (length(precision_k)) {
        precision_k <- unique(precision_k)
        for (k in precision_k) {
            k_use <- min(k, n_comp)
            sets[[length(sets) + 1]] <- data.frame(
                set_id = paste0("top_k_", k_use),
                k = k_use,
                frac = NA_real_,
                stringsAsFactors = FALSE
            )
        }
    }
    if (!is.null(top_frac) && is.finite(top_frac) && top_frac > 0 && n_comp > 0) {
        k_frac <- max(1L, ceiling(top_frac * n_comp))
        sets[[length(sets) + 1]] <- data.frame(
            set_id = paste0("top_frac_", formatC(top_frac * 100, format = "f", digits = 2), "pct"),
            k = k_frac,
            frac = top_frac,
            stringsAsFactors = FALSE
        )
    }
    if (!length(sets)) {
        return(data.frame(
            set_id = character(0),
            k = integer(0),
            frac = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    out <- do.call(rbind, sets)
    out[!duplicated(out$set_id), , drop = FALSE]
}

#' Compute PR curve points.
#' @keywords internal
.external_reference_pr_curve <- function(labels, scores) {
    ok <- !is.na(labels) & !is.na(scores)
    labels <- labels[ok]
    scores <- scores[ok]
    if (!length(labels)) {
        return(data.frame(recall = numeric(0), precision = numeric(0)))
    }
    labels <- as.integer(labels) > 0
    n_pos <- sum(labels)
    if (n_pos == 0) {
        return(data.frame(recall = numeric(0), precision = numeric(0)))
    }
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    tp <- cumsum(labels)
    fp <- cumsum(!labels)
    recall <- tp / n_pos
    precision <- tp / pmax(tp + fp, 1)
    data.frame(
        recall = c(0, recall),
        precision = c(1, precision)
    )
}

#' Compute ROC curve points.
#' @keywords internal
.external_reference_roc_curve <- function(labels, scores) {
    ok <- !is.na(labels) & !is.na(scores)
    labels <- labels[ok]
    scores <- scores[ok]
    if (!length(labels)) {
        return(data.frame(fpr = numeric(0), tpr = numeric(0)))
    }
    labels <- as.integer(labels) > 0
    n_pos <- sum(labels)
    n_neg <- sum(!labels)
    if (n_pos == 0 || n_neg == 0) {
        return(data.frame(fpr = numeric(0), tpr = numeric(0)))
    }
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    tp <- cumsum(labels)
    fp <- cumsum(!labels)
    tpr <- tp / n_pos
    fpr <- fp / n_neg
    data.frame(
        fpr = c(0, fpr),
        tpr = c(0, tpr)
    )
}

#' Compute enrichment for defined top-edge sets.
#' @keywords internal
.external_reference_top_enrichment <- function(edge_df,
                                               score_threshold,
                                               expected_rate,
                                               set_defs) {
    if (is.null(set_defs) || !nrow(set_defs)) {
        return(data.frame(
            set_id = character(0),
            k = integer(0),
            frac = numeric(0),
            n_top = integer(0),
            n_comparable = integer(0),
            n_pos = integer(0),
            positive_rate = numeric(0),
            enrichment_ratio = numeric(0),
            expected_positive_rate = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    if (!nrow(edge_df)) {
        empty <- set_defs
        empty$n_top <- 0L
        empty$n_comparable <- 0L
        empty$n_pos <- 0L
        empty$positive_rate <- NA_real_
        empty$enrichment_ratio <- NA_real_
        empty$expected_positive_rate <- expected_rate
        return(empty)
    }
    comparable_idx <- !is.na(edge_df$string_score) & !is.na(edge_df$score)
    edge_df <- edge_df[comparable_idx, , drop = FALSE]
    if (!nrow(edge_df)) {
        empty <- set_defs
        empty$n_top <- 0L
        empty$n_comparable <- 0L
        empty$n_pos <- 0L
        empty$positive_rate <- NA_real_
        empty$enrichment_ratio <- NA_real_
        empty$expected_positive_rate <- expected_rate
        return(empty)
    }
    ord <- order(edge_df$score, decreasing = TRUE, na.last = NA)
    edge_df <- edge_df[ord, , drop = FALSE]

    res <- lapply(seq_len(nrow(set_defs)), function(i) {
        k <- set_defs$k[i]
        n_top <- min(k, nrow(edge_df))
        top_df <- if (n_top > 0) edge_df[seq_len(n_top), , drop = FALSE] else edge_df[integer(0), , drop = FALSE]
        n_comparable <- n_top
        n_pos <- sum(top_df$string_score >= score_threshold, na.rm = TRUE)
        positive_rate <- if (n_top > 0) n_pos / n_top else NA_real_
        enrichment_ratio <- NA_real_
        if (!is.na(expected_rate) && expected_rate > 0 && !is.na(positive_rate)) {
            enrichment_ratio <- positive_rate / expected_rate
        }
        data.frame(
            set_id = set_defs$set_id[i],
            k = set_defs$k[i],
            frac = set_defs$frac[i],
            n_top = n_top,
            n_comparable = n_comparable,
            n_pos = n_pos,
            positive_rate = positive_rate,
            enrichment_ratio = enrichment_ratio,
            expected_positive_rate = expected_rate,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, res)
}

#' Resolve writeback path for external reference scores.
#' @keywords internal
.external_ref_store_path <- function(scope_obj,
                                     grid_name,
                                     stats_layer,
                                     graph_slot) {
    if (!inherits(scope_obj, "scope_object")) {
        stop("scope_obj must be a scope_object for writeback.")
    }
    if (is.null(scope_obj@stats) || !length(scope_obj@stats)) {
        stop("scope_obj@stats is empty; cannot write external reference scores.")
    }
    if (!grid_name %in% names(scope_obj@stats)) {
        stop("grid_name not found in scope_obj@stats: ", grid_name)
    }
    stats_layer_obj <- scope_obj@stats[[grid_name]]
    if (is.null(stats_layer_obj) || !stats_layer %in% names(stats_layer_obj)) {
        stop("stats_layer not found in scope_obj@stats[[grid_name]]: ", stats_layer)
    }
    layer_obj <- stats_layer_obj[[stats_layer]]
    if (is.null(layer_obj) || !graph_slot %in% names(layer_obj)) {
        stop("graph_slot not found in stats_layer: ", graph_slot)
    }
    graph_obj <- layer_obj[[graph_slot]]
    store_in_graph <- is.list(graph_obj) && !igraph::is_igraph(graph_obj)
    list(store_in_graph = store_in_graph)
}

#' Write external reference scores back to scope_obj@stats.
#' @keywords internal
.external_ref_writeback <- function(scope_obj,
                                    grid_name,
                                    stats_layer,
                                    graph_slot,
                                    payload) {
    path <- .external_ref_store_path(scope_obj, grid_name, stats_layer, graph_slot)
    stats_layer_obj <- scope_obj@stats[[grid_name]]
    layer_obj <- stats_layer_obj[[stats_layer]]

    if (isTRUE(path$store_in_graph)) {
        graph_obj <- layer_obj[[graph_slot]]
        if (is.null(graph_obj$external_reference)) graph_obj$external_reference <- list()
        if (is.null(graph_obj$external_reference$STRINGdb)) {
            graph_obj$external_reference$STRINGdb <- list()
        }
        graph_obj$external_reference$STRINGdb <- payload
        layer_obj[[graph_slot]] <- graph_obj
    } else {
        if (is.null(layer_obj$external_reference)) layer_obj$external_reference <- list()
        if (is.null(layer_obj$external_reference$STRINGdb)) {
            layer_obj$external_reference$STRINGdb <- list()
        }
        layer_obj$external_reference$STRINGdb[[graph_slot]] <- payload
    }

    stats_layer_obj[[stats_layer]] <- layer_obj
    scope_obj@stats[[grid_name]] <- stats_layer_obj
    scope_obj
}

#' Read external reference payload from scope_obj@stats.
#' @keywords internal
.external_ref_read <- function(scope_obj,
                               grid_name,
                               stats_layer,
                               graph_slot,
                               ref = "STRINGdb") {
    if (!inherits(scope_obj, "scope_object")) {
        return(NULL)
    }
    stats_root <- scope_obj@stats
    if (is.null(stats_root) || !grid_name %in% names(stats_root)) {
        return(NULL)
    }
    layer_root <- stats_root[[grid_name]]
    if (is.null(layer_root) || !stats_layer %in% names(layer_root)) {
        return(NULL)
    }
    layer_obj <- layer_root[[stats_layer]]
    if (is.null(layer_obj)) return(NULL)

    if (!is.null(layer_obj$external_reference) && !is.null(layer_obj$external_reference[[ref]])) {
        ref_obj <- layer_obj$external_reference[[ref]]
        if (!is.null(ref_obj[[graph_slot]])) return(ref_obj[[graph_slot]])
    }
    if (!graph_slot %in% names(layer_obj)) return(NULL)
    graph_obj <- layer_obj[[graph_slot]]
    if (is.list(graph_obj) && !igraph::is_igraph(graph_obj)) {
        if (!is.null(graph_obj$external_reference) && !is.null(graph_obj$external_reference[[ref]])) {
            return(graph_obj$external_reference[[ref]])
        }
    }
    NULL
}

#' Retrieve cached STRINGdb artifacts.
#' @keywords internal
.stringdb_cache_get <- function(scope_obj = NULL,
                                cache_dir = NULL,
                                cache_type,
                                cache_key) {
    if (is.null(cache_key) || is.na(cache_key) || !nzchar(cache_key)) {
        return(list(value = NULL, cache_hit = FALSE, scope_obj = scope_obj))
    }

    env <- .stringdb_cache_env()
    if (exists(cache_key, envir = env, inherits = FALSE)) {
        cache_val <- get(cache_key, envir = env, inherits = FALSE)
        if (is.list(cache_val) && !is.null(cache_val[[cache_type]])) {
            return(list(value = cache_val[[cache_type]], cache_hit = TRUE, scope_obj = scope_obj))
        }
    }

    if (!is.null(cache_dir)) {
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
        cache_path <- file.path(cache_dir, paste0("stringdb_", cache_type, "_", cache_key, ".rds"))
        if (file.exists(cache_path)) {
            return(list(value = readRDS(cache_path), cache_hit = TRUE, scope_obj = scope_obj))
        }
    }

    list(value = NULL, cache_hit = FALSE, scope_obj = scope_obj)
}

#' Store cached STRINGdb artifacts.
#' @keywords internal
.stringdb_cache_put <- function(scope_obj = NULL,
                                cache_dir = NULL,
                                cache_type,
                                cache_key,
                                value) {
    if (is.null(cache_key) || is.na(cache_key) || !nzchar(cache_key)) {
        return(list(scope_obj = scope_obj))
    }

    env <- .stringdb_cache_env()
    cache_val <- if (exists(cache_key, envir = env, inherits = FALSE)) {
        get(cache_key, envir = env, inherits = FALSE)
    } else {
        list()
    }
    if (!is.list(cache_val)) cache_val <- list()
    cache_val[[cache_type]] <- value
    assign(cache_key, cache_val, envir = env)

    if (!is.null(cache_dir)) {
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
        cache_path <- file.path(cache_dir, paste0("stringdb_", cache_type, "_", cache_key, ".rds"))
        saveRDS(value, cache_path)
    }

    list(scope_obj = scope_obj)
}

#' Map genes to STRING IDs with caching.
#' @keywords internal
.stringdb_map_genes_cached <- function(string_db,
                                       genes,
                                       input_id_type = "gene",
                                       scope_obj = NULL,
                                       cache_dir = NULL) {
    genes <- unique(as.character(genes))
    if (!length(genes)) {
        return(list(mapping = setNames(character(0), character(0)), cache_hit = FALSE, scope_obj = scope_obj))
    }
    gene_hash <- .external_reference_hash(sort(genes))
    cache_key <- .external_reference_hash(list(
        species = string_db$species,
        input_id_type = input_id_type,
        gene_hash = gene_hash
    ))

    cache <- .stringdb_cache_get(scope_obj, cache_dir, "mapping", cache_key)
    if (isTRUE(cache$cache_hit)) {
        return(list(mapping = cache$value, cache_hit = TRUE, scope_obj = cache$scope_obj))
    }

    df <- data.frame(stats::setNames(list(genes), input_id_type), stringsAsFactors = FALSE)
    mapped <- string_db$map(df, input_id_type, removeUnmappedRows = FALSE)
    mapping <- setNames(mapped$STRING_id, mapped[[input_id_type]])

    cache_put <- .stringdb_cache_put(scope_obj, cache_dir, "mapping", cache_key, mapping)
    list(mapping = mapping, cache_hit = FALSE, scope_obj = cache_put$scope_obj)
}

#' Fetch STRING interactions with caching.
#' @keywords internal
.stringdb_get_interactions_cached <- function(string_db,
                                              string_ids,
                                              score_threshold,
                                              scope_obj = NULL,
                                              cache_dir = NULL) {
    ids <- unique(na.omit(as.character(string_ids)))
    if (!length(ids)) {
        empty <- data.frame(from = character(0), to = character(0), combined_score = numeric(0))
        return(list(interactions = empty, cache_hit = FALSE, scope_obj = scope_obj))
    }
    ids_hash <- .external_reference_hash(sort(ids))
    cache_key <- .external_reference_hash(list(
        species = string_db$species,
        score_threshold = score_threshold,
        ids_hash = ids_hash
    ))

    cache <- .stringdb_cache_get(scope_obj, cache_dir, "interactions", cache_key)
    if (isTRUE(cache$cache_hit)) {
        return(list(interactions = cache$value, cache_hit = TRUE, scope_obj = cache$scope_obj))
    }

    interactions <- tryCatch(
        string_db$get_interactions(ids),
        error = function(e) {
            msg <- conditionMessage(e)
            if (grepl("edge data frame contains NAs", msg) || grepl("graph_from_data_frame", msg)) {
                return(.stringdb_get_interactions_fallback(string_db, ids))
            }
            stop(e)
        }
    )
    if (!is.data.frame(interactions) || !nrow(interactions)) {
        interactions <- data.frame(from = character(0), to = character(0), combined_score = numeric(0))
    } else {
        interactions <- interactions[, c("from", "to", "combined_score"), drop = FALSE]
    }

    cache_put <- .stringdb_cache_put(scope_obj, cache_dir, "interactions", cache_key, interactions)
    list(interactions = interactions, cache_hit = FALSE, scope_obj = cache_put$scope_obj)
}

#' Score STRING reference agreement and return structured payload.
#' @keywords internal
.external_reference_score_stringdb <- function(x,
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
                                               verbose = getOption("geneSCOPE.verbose", TRUE)) {
    parent <- "scoreExternalReferenceSTRING"
    null_method <- match.arg(null_method)
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
    cache_dir_use <- cache_dir
    if (is.null(cache_dir_use) && !inherits(x, "scope_object")) {
        cache_dir_use <- tempdir()
    }
    string_db <- .stringdb_connect(
        species = species,
        score_threshold = 0,
        cache_dir = cache_dir_use
    )
    genes <- igraph::V(g)$name
    cache_scope <- if (inherits(x, "scope_object")) x else NULL

    mapping_res <- .stringdb_map_genes_cached(
        string_db,
        genes,
        input_id_type = input_id_type,
        scope_obj = cache_scope,
        cache_dir = cache_dir_use
    )
    mapping <- mapping_res$mapping
    cache_scope <- mapping_res$scope_obj
    interactions_res <- .stringdb_get_interactions_cached(
        string_db,
        mapping,
        score_threshold = string_score_threshold,
        scope_obj = cache_scope,
        cache_dir = cache_dir_use
    )
    interactions <- interactions_res$interactions
    cache_scope <- interactions_res$scope_obj
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
    node_coverage <- if (n_genes) n_mapped_genes / n_genes else NA_real_

    ok_idx <- !is.na(edge_df$label) & !is.na(edge_df$score)
    n_comparable_edges <- sum(ok_idx)
    labels_ok <- edge_df$label[ok_idx]
    scores_ok <- edge_df$score[ok_idx]
    n_pos_labels <- sum(labels_ok, na.rm = TRUE)
    n_neg_labels <- sum(!labels_ok, na.rm = TRUE)

    edge_coverage <- if (n_edges) n_comparable_edges / n_edges else NA_real_
    positive_rate <- if (n_comparable_edges) mean(labels_ok) else NA_real_

    auroc_reason <- "ok"
    auprc_reason <- "ok"
    auroc_defined <- TRUE
    auprc_defined <- TRUE
    if (n_comparable_edges == 0) {
        auroc_defined <- FALSE
        auprc_defined <- FALSE
        auroc_reason <- "no_comparable_edges"
        auprc_reason <- "no_comparable_edges"
    } else if (n_pos_labels == 0) {
        auroc_defined <- FALSE
        auprc_defined <- FALSE
        auroc_reason <- "no_positive_labels"
        auprc_reason <- "no_positive_labels"
    } else if (n_neg_labels == 0) {
        auroc_defined <- FALSE
        auprc_defined <- FALSE
        auroc_reason <- "no_negative_labels"
        auprc_reason <- "no_negative_labels"
    }

    auprc <- NA_real_
    auroc <- NA_real_
    if (isTRUE(auprc_defined)) {
        auprc <- .compute_auprc(labels_ok, scores_ok)
    }
    if (isTRUE(auroc_defined)) {
        auroc <- .compute_auroc(labels_ok, scores_ok)
    }

    prec_vals <- .external_reference_precision_at_k(edge_df$label, edge_df$score, k = precision_k)

    mapped_ids <- unique(na.omit(as.character(mapping)))
    bg_gene_count <- length(mapped_ids)
    n_pairs_bg <- bg_gene_count * (bg_gene_count - 1) / 2
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

    set_defs <- .external_reference_top_sets(
        n_edges,
        precision_k,
        top_frac = top_frac,
        n_comparable = n_comparable_edges
    )
    top_enrich <- .external_reference_top_enrichment(
        edge_df,
        score_threshold = string_score_threshold,
        expected_rate = enrich$expected_positive_rate,
        set_defs = set_defs
    )

    ranking_key <- "abs(weight)"
    curves <- list(
        pr = .external_reference_pr_curve(edge_df$label, edge_df$score),
        roc = .external_reference_roc_curve(edge_df$label, edge_df$score)
    )

    scores <- list(
        auprc = auprc,
        auprc_defined = auprc_defined,
        auprc_reason = auprc_reason,
        auroc = auroc,
        auroc_defined = auroc_defined,
        auroc_reason = auroc_reason,
        ranking_key = ranking_key,
        topk_precision = prec_vals
    )

    mapping_summary <- list(
        n_input_genes = length(genes),
        n_mapped = n_mapped_genes,
        coverage = node_coverage
    )

    coverage <- list(
        n_edges = n_edges,
        n_comparable_edges = n_comparable_edges,
        edge_coverage = edge_coverage,
        positive_rate = positive_rate
    )

    null_summary <- NULL
    if (isTRUE(run_null) && n_null > 0L && n_edges > 0L) {
        null_fn <- if (identical(null_method, "label_permutation")) {
            .external_reference_label_null
        } else {
            .external_reference_degree_null
        }
        null_summary <- null_fn(
            graph = g,
            edge_weight = edge_weight,
            gene_to_string = mapping,
            interactions = interactions,
            score_threshold = string_score_threshold,
            precision_k = precision_k,
            precision_obs = prec_vals,
            enrichment_obs = top_enrich,
            set_defs = set_defs,
            n_null = n_null,
            seed = null_seed,
            expected_positive_rate = enrich$expected_positive_rate
        )
    }

    params <- list(
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot,
        edge_weight = edge_weight,
        species = species,
        input_id_type = input_id_type,
        string_score_threshold = string_score_threshold,
        precision_k = as.integer(precision_k),
        top_frac = top_frac,
        ranking_key = ranking_key,
        run_null = isTRUE(run_null),
        null_method = null_method,
        n_null = as.integer(n_null),
        null_seed = null_seed
    )

    gene_hash <- .external_reference_hash(sort(genes))
    meta <- list(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
        gene_universe = list(n = n_genes, hash = gene_hash),
        species = species,
        string_version = string_db$version,
        params = params,
        params_digest = .external_reference_params_digest(params),
        cache_hit = isTRUE(mapping_res$cache_hit) && isTRUE(interactions_res$cache_hit),
        mapping_cache_hit = isTRUE(mapping_res$cache_hit),
        interactions_cache_hit = isTRUE(interactions_res$cache_hit)
    )

    enrichment <- list(
        n_pos_graph = n_pos_graph,
        n_graph = n_comparable_edges,
        n_pos_bg = n_pos_bg,
        n_pairs_bg = n_pairs_bg,
        bg_gene_count = bg_gene_count,
        expected_positive_rate = enrich$expected_positive_rate,
        enrichment_ratio = enrich$enrichment_ratio,
        enrichment_p = enrich$enrichment_p,
        top_sets = top_enrich
    )

    null_payload <- NULL
    if (!is.null(null_summary)) {
        null_payload <- list(
            method = null_method,
            B = as.integer(n_null),
            summary = null_summary
        )
    }

    notes_rebuttal <- paste0(
        "STRING-based evaluation uses known interactions as positives; AUPRC/top-K precision and coverage ",
        "offer objective measures of concordance while reflecting public database biases."
    )

    payload <- list(
        meta = meta,
        mapping = mapping_summary,
        coverage = coverage,
        scores = scores,
        enrichment = enrichment,
        null = null_payload,
        curves = curves,
        notes_rebuttal = notes_rebuttal
    )

    if (isTRUE(keep_details)) payload$details <- edge_df
    list(payload = payload, cache_scope = cache_scope)
}

#' Flatten external reference payload into a data.frame.
#' @keywords internal
.external_reference_scores_to_df <- function(payload,
                                             grid_name,
                                             stats_layer,
                                             graph_slot) {
    if (is.null(payload) || is.null(payload$scores)) {
        return(data.frame())
    }
    score_list <- list(
        auprc = payload$scores$auprc,
        auprc_defined = if (!is.null(payload$scores$auprc_defined)) payload$scores$auprc_defined else NA,
        auprc_reason = if (!is.null(payload$scores$auprc_reason)) payload$scores$auprc_reason else NA_character_,
        auroc = payload$scores$auroc,
        auroc_defined = if (!is.null(payload$scores$auroc_defined)) payload$scores$auroc_defined else NA,
        auroc_reason = if (!is.null(payload$scores$auroc_reason)) payload$scores$auroc_reason else NA_character_,
        ranking_key = if (!is.null(payload$scores$ranking_key)) payload$scores$ranking_key else NA_character_
    )
    if (!is.null(payload$scores$topk_precision)) {
        score_list <- c(score_list, as.list(payload$scores$topk_precision))
    }

    mapping_list <- payload$mapping
    if (is.null(mapping_list)) mapping_list <- list()

    coverage_list <- payload$coverage
    if (is.null(coverage_list)) coverage_list <- list()

    enrichment_list <- list()
    top_enrich_list <- list()
    if (!is.null(payload$enrichment)) {
        enrichment_list <- payload$enrichment
        top_sets <- enrichment_list$top_sets
        enrichment_list$top_sets <- NULL
        if (is.data.frame(top_sets) && nrow(top_sets)) {
            top_enrich_list <- lapply(seq_len(nrow(top_sets)), function(i) {
                row <- top_sets[i, , drop = FALSE]
                prefix <- row$set_id
                vals <- as.list(row)
                vals$set_id <- NULL
                names(vals) <- paste0(prefix, "_", names(vals))
                vals
            })
            top_enrich_list <- do.call(c, top_enrich_list)
        }
    }

    null_list <- list()
    if (!is.null(payload$null)) {
        null_list <- list(
            null_method = payload$null$method,
            null_B = payload$null$B
        )
        if (is.data.frame(payload$null$summary) && nrow(payload$null$summary)) {
            null_metrics <- lapply(seq_len(nrow(payload$null$summary)), function(i) {
                row <- payload$null$summary[i, , drop = FALSE]
                prefix <- row$metric
                vals <- as.list(row)
                vals$metric <- NULL
                names(vals) <- paste0("null_", prefix, "_", names(vals))
                vals
            })
            null_metrics <- do.call(c, null_metrics)
            null_list <- c(null_list, null_metrics)
        }
    }

    meta <- payload$meta
    mapping_id_type <- NA_character_
    if (!is.null(meta$params) && !is.null(meta$params$input_id_type)) {
        mapping_id_type <- meta$params$input_id_type
    }
    meta_list <- list(
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot,
        species = meta$species,
        string_version = meta$string_version,
        mapping_id_type = mapping_id_type,
        timestamp = meta$timestamp,
        gene_universe_n = meta$gene_universe$n,
        gene_universe_hash = meta$gene_universe$hash,
        params_digest = meta$params_digest,
        cache_hit = meta$cache_hit,
        mapping_cache_hit = meta$mapping_cache_hit,
        interactions_cache_hit = meta$interactions_cache_hit
    )

    row <- c(meta_list, mapping_list, coverage_list, score_list, enrichment_list, top_enrich_list, null_list)
    as.data.frame(row, stringsAsFactors = FALSE)
}

#' Score module quality using Hotspot-style statistics.
#' @description
#' Computes member- and module-level cohesion/separation metrics from a graph
#' and membership labels. Optionally runs a label-permutation null to derive
#' z/p/FDR statistics. Results are returned as two tables. Writeback happens
#' only when `writeback = TRUE` and input is a `scope_object`.
#' @param x A `scope_object`, `igraph`, or list with `graph` and optional `membership`.
#' @param membership Optional membership vector named by gene.
#' @param grid_name Grid layer name for `scope_object` input (default "grid30").
#' @param stats_layer Stats layer name for `scope_object` input (default "LeeStats_Xz").
#' @param graph_slot Graph slot name for `scope_object` input (default "g_consensus").
#' @param edge_weight Edge attribute used as weight (default "weight").
#' @param n_perm Number of label permutations for null (default 0).
#' @param seed Random seed for permutations.
#' @param use_cpp Use C++ permutation if available (default TRUE).
#' @param compute_conductance Whether to compute conductance (default TRUE).
#' @param writeback Write results to `scope_obj@stats[[grid_name]][[stats_layer]]$module_quality`.
#' @param verbose Emit progress messages when TRUE.
#' @param ... Reserved for future extensions.
#' @return A list with `member_scores`, `module_scores`, and `meta`.
#' @examples
#' \\dontrun{
#' out <- ScoreModuleQuality(scope_obj, grid_name = "grid30", n_perm = 50)
#' }
#' @export
ScoreModuleQuality <- function(
    x,
    membership = NULL,
    grid_name = "grid30",
    stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    edge_weight = "weight",
    n_perm = 0L,
    seed = 1,
    use_cpp = TRUE,
    compute_conductance = TRUE,
    writeback = FALSE,
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...) {
    parent <- "ScoreModuleQuality"
    step01 <- .log_step(parent, "S01", "resolve graph and membership", verbose)
    step01$enter()
    graph <- .mq_extract_graph(
        x = x,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot
    )
    memb <- .mq_extract_membership(
        x = x,
        graph = graph,
        membership = membership,
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot
    )
    step01$done(paste0("nodes=", igraph::vcount(graph)))

    step02 <- .log_step(parent, "S02", "prepare edge table", verbose)
    step02$enter()
    edge_df <- .mq_standardize_weights(graph, edge_weight = edge_weight)
    step02$done(paste0("edges=", nrow(edge_df)))

    memb_std <- .mq_standardize_membership(memb)
    membership_id <- memb_std$module_id
    module_levels <- memb_std$module_levels

    member_scores <- .mq_member_scores(edge_df, membership_id, module_levels)
    module_scores <- .mq_module_scores(edge_df, membership_id, module_levels, compute_conductance)

    null_meta <- list(
        method = "label_permutation",
        n_perm = as.integer(n_perm),
        seed = as.integer(seed),
        used_cpp = FALSE
    )
    null_summary <- NULL

    if (isTRUE(n_perm > 0)) {
        null_res <- NULL
        if (isTRUE(use_cpp)) {
            null_res <- .mq_score_null_cpp(
                edge_df = edge_df,
                membership_id = membership_id,
                n_perm = n_perm,
                seed = seed,
                compute_conductance = compute_conductance
            )
            if (!is.null(null_res)) null_meta$used_cpp <- TRUE
        }
        if (is.null(null_res)) {
            null_res <- .mq_score_null_R(
                edge_df = edge_df,
                membership_id = membership_id,
                module_levels = module_levels,
                n_perm = n_perm,
                seed = seed,
                compute_conductance = compute_conductance
            )
        }

        if (!is.null(null_res$member_margin) && nrow(member_scores)) {
            zp <- .mq_z_p_from_null(
                member_scores$member_margin,
                null_res$member_margin$mean,
                null_res$member_margin$sd
            )
            member_scores$member_z <- zp$z
            member_scores$member_p <- zp$p
            member_scores$member_fdr <- zp$fdr
        } else {
            member_scores$member_z <- NA_real_
            member_scores$member_p <- NA_real_
            member_scores$member_fdr <- NA_real_
        }

        if (!is.null(null_res$module_separation) && nrow(module_scores)) {
            zp <- .mq_z_p_from_null(
                module_scores$separation,
                null_res$module_separation$mean,
                null_res$module_separation$sd
            )
            module_scores$separation_z <- zp$z
            module_scores$separation_p <- zp$p
            module_scores$separation_fdr <- zp$fdr
        } else {
            module_scores$separation_z <- NA_real_
            module_scores$separation_p <- NA_real_
            module_scores$separation_fdr <- NA_real_
        }

        if (isTRUE(compute_conductance)) {
            if (!is.null(null_res$module_conductance) && nrow(module_scores)) {
                zp <- .mq_z_p_from_null(
                    module_scores$conductance,
                    null_res$module_conductance$mean,
                    null_res$module_conductance$sd
                )
                module_scores$conductance_z <- zp$z
                module_scores$conductance_p <- zp$p
                module_scores$conductance_fdr <- zp$fdr
            } else {
                module_scores$conductance_z <- NA_real_
                module_scores$conductance_p <- NA_real_
                module_scores$conductance_fdr <- NA_real_
            }
        } else {
            module_scores$conductance_z <- NA_real_
            module_scores$conductance_p <- NA_real_
            module_scores$conductance_fdr <- NA_real_
        }

        null_summary <- list(
            member_margin = null_res$member_margin,
            module_separation = null_res$module_separation,
            module_conductance = null_res$module_conductance
        )
    } else {
        member_scores$member_z <- NA_real_
        member_scores$member_p <- NA_real_
        member_scores$member_fdr <- NA_real_
        module_scores$separation_z <- NA_real_
        module_scores$separation_p <- NA_real_
        module_scores$separation_fdr <- NA_real_
        module_scores$conductance_z <- NA_real_
        module_scores$conductance_p <- NA_real_
        module_scores$conductance_fdr <- NA_real_
    }

    top_n <- getOption("geneSCOPE.module_quality.top_n", 5L)
    member_scores <- .mq_finalize_member_table(member_scores)
    module_scores <- .mq_finalize_module_table(
        module_scores = module_scores,
        member_scores = member_scores,
        null_summary = null_summary,
        top_n = top_n
    )

    meta <- list(
        schema_version = "1.0",
        created_at = as.character(Sys.time()),
        run_id = .mq_make_run_id(),
        n_nodes = length(membership_id),
        n_edges = nrow(edge_df),
        module_levels = module_levels,
        params = list(
            grid_name = grid_name,
            stats_layer = stats_layer,
            graph_slot = graph_slot,
            edge_weight = edge_weight,
            n_perm = as.integer(n_perm),
            seed = as.integer(seed),
            compute_conductance = isTRUE(compute_conductance)
        ),
        null = null_meta
    )

    payload <- list(
        meta = meta,
        member_scores = member_scores,
        module_scores = module_scores,
        null_summary = null_summary
    )

    if (isTRUE(writeback) && inherits(x, "scope_object")) {
        scope_obj <- .mq_writeback_scope(
            scope_obj = x,
            grid_name = grid_name,
            stats_layer = stats_layer,
            payload = payload
        )
        payload$scope_obj <- scope_obj
    }

    payload
}

#' Extract stored module-quality scores.
#' @description
#' Returns member and module tables stored under `module_quality` in a scope object.
#' @param x A `scope_object` or a list containing `module_quality`.
#' @param grid_name Optional grid filter.
#' @param stats_layer Optional stats layer filter.
#' @return A list with `member_scores` and `module_scores`.
#' @export
GetModuleQualityScores <- function(x,
                                   grid_name = NULL,
                                   stats_layer = NULL) {
    if (inherits(x, "scope_object")) {
        member_rows <- list()
        module_rows <- list()
        stats_root <- x@stats
        for (gname in names(stats_root)) {
            if (!is.null(grid_name) && !identical(gname, grid_name)) next
            layer_root <- stats_root[[gname]]
            for (lname in names(layer_root)) {
                if (!is.null(stats_layer) && !identical(lname, stats_layer)) next
                entry <- layer_root[[lname]]$module_quality
                if (is.null(entry)) next

                if (!is.null(entry$hotspot_style) && length(entry$hotspot_style)) {
                    for (run_id in names(entry$hotspot_style)) {
                        run <- entry$hotspot_style[[run_id]]
                        if (!is.null(run$member_table) && nrow(run$member_table)) {
                            df <- run$member_table
                            df$grid_name <- gname
                            df$stats_layer <- lname
                            df$run_id <- run_id
                            member_rows[[length(member_rows) + 1]] <- df
                        }
                        if (!is.null(run$module_table) && nrow(run$module_table)) {
                            df <- run$module_table
                            df$grid_name <- gname
                            df$stats_layer <- lname
                            df$run_id <- run_id
                            module_rows[[length(module_rows) + 1]] <- df
                        }
                    }
                    next
                }

                if (!is.null(entry$member_scores) && nrow(entry$member_scores)) {
                    df <- entry$member_scores
                    df$grid_name <- gname
                    df$stats_layer <- lname
                    member_rows[[length(member_rows) + 1]] <- df
                }
                if (!is.null(entry$module_scores) && nrow(entry$module_scores)) {
                    df <- entry$module_scores
                    df$grid_name <- gname
                    df$stats_layer <- lname
                    module_rows[[length(module_rows) + 1]] <- df
                }
            }
        }
        member_scores <- if (length(member_rows)) do.call(rbind, member_rows) else data.frame()
        module_scores <- if (length(module_rows)) do.call(rbind, module_rows) else data.frame()
        return(list(member_scores = member_scores, module_scores = module_scores))
    }

    if (is.list(x) && !is.null(x$module_quality)) {
        entry <- x$module_quality
        if (!is.null(entry$hotspot_style) && length(entry$hotspot_style)) {
            member_rows <- list()
            module_rows <- list()
            for (run_id in names(entry$hotspot_style)) {
                run <- entry$hotspot_style[[run_id]]
                if (!is.null(run$member_table) && nrow(run$member_table)) {
                    df <- run$member_table
                    df$run_id <- run_id
                    member_rows[[length(member_rows) + 1]] <- df
                }
                if (!is.null(run$module_table) && nrow(run$module_table)) {
                    df <- run$module_table
                    df$run_id <- run_id
                    module_rows[[length(module_rows) + 1]] <- df
                }
            }
            return(list(
                member_scores = if (length(member_rows)) do.call(rbind, member_rows) else data.frame(),
                module_scores = if (length(module_rows)) do.call(rbind, module_rows) else data.frame()
            ))
        }

        return(list(
            member_scores = if (!is.null(entry$member_scores)) entry$member_scores else data.frame(),
            module_scores = if (!is.null(entry$module_scores)) entry$module_scores else data.frame()
        ))
    }

    list(member_scores = data.frame(), module_scores = data.frame())
}

#' Summarize stored module-quality scores.
#' @description
#' Produces module- and member-level summary tables from stored results.
#' @param x A `scope_object` containing stored module-quality scores.
#' @param grid_name Optional grid filter.
#' @param stats_layer Optional stats layer filter.
#' @return A list with `module_table`, `member_table`, and `summary_text`.
#' @export
SummarizeModuleQualityScores <- function(x,
                                         grid_name = NULL,
                                         stats_layer = NULL) {
    if (!inherits(x, "scope_object")) {
        stop("SummarizeModuleQualityScores requires a scope_object.")
    }

    scores <- GetModuleQualityScores(
        x = x,
        grid_name = grid_name,
        stats_layer = stats_layer
    )
    top_n <- getOption("geneSCOPE.module_quality.top_n", 5L)
    tables <- .mq_summary_tables(
        member_scores = scores$member_scores,
        module_scores = scores$module_scores,
        top_n = top_n
    )

    per_graph <- data.frame()
    if (nrow(scores$module_scores)) {
        groups <- interaction(scores$module_scores$grid_name, scores$module_scores$stats_layer, drop = TRUE)
        rows <- lapply(split(scores$module_scores, groups), function(df) {
            data.frame(
                grid_name = df$grid_name[1],
                stats_layer = df$stats_layer[1],
                n_modules = nrow(df),
                mean_cohesion = mean(df$cohesion, na.rm = TRUE),
                mean_separation = mean(df$separation, na.rm = TRUE),
                frac_positive_sep = mean(df$separation > 0, na.rm = TRUE),
                stringsAsFactors = FALSE
            )
        })
        per_graph <- do.call(rbind, rows)
    }

    summary_text <- paste0(
        "Summarized module quality scores for ",
        nrow(per_graph),
        " graph(s)."
    )

    list(
        module_table = tables$module_table,
        member_table = tables$member_table,
        per_graph = per_graph,
        per_module = tables$module_table,
        summary_text = summary_text
    )
}

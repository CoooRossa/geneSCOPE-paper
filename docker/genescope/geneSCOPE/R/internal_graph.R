# NetworKit helpers
#' Check if the Python 'networkit' module is available through reticulate.
#' @description
#' Internal helper for `.nk_available`.
#' @keywords internal
#' @return Logical flag indicating module availability.
.nk_available <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
    ok <- FALSE
    try({ ok <- isTRUE(reticulate::py_module_available("networkit")) }, silent = TRUE)
    isTRUE(ok)
}

#' Run NetworKit ParallelLeiden across multiple seeds (igraph fallback when absent).
#' @description
#' Internal helper for `.nk_parallel_leiden_mruns`.
#' @keywords internal
#' @param vertex_names Character vector of vertex identifiers.
#' @param edge_table Data frame of edges with columns `from` and `to`.
#' @param edge_weights Numeric vector of edge weights aligned to `edge_table`.
#' @param gamma Resolution parameter.
#' @param iterations Integer number of Leiden iterations.
#' @param randomize Logical; whether to randomize the Leiden algorithm.
#' @param threads Integer thread count for NetworKit.
#' @param seeds Integer vector of seeds; one membership column is returned per seed.
#' @return Integer matrix of memberships (rows = vertices, cols = runs).
.nk_parallel_leiden_mruns <- function(vertex_names, edge_table, edge_weights, gamma = 1.0,
                                      iterations = 10L, randomize = TRUE,
                                      threads = 1L, seeds = 1L) {
    if (.nk_available()) {
        nk <- reticulate::import("networkit", delay_load = TRUE)
        threads <- as.integer(threads)
        if (!is.na(threads) && threads > 0) try(nk$setNumberOfThreads(threads), silent = TRUE)

        vcount <- length(vertex_names)
        g <- nk$Graph(as.integer(vcount), weighted = TRUE, directed = FALSE)
        if (!is.null(edge_table) && nrow(edge_table)) {
            vid <- seq_len(vcount) - 1L
            map <- structure(vid, names = as.character(vertex_names))
            u <- map[as.character(edge_table$from)]
            v <- map[as.character(edge_table$to)]
            w <- as.numeric(edge_weights)
            keep <- !is.na(u) & !is.na(v) & is.finite(w)
            if (any(keep)) {
                u <- as.integer(u[keep]); v <- as.integer(v[keep]); w <- as.numeric(w[keep])
                for (i in seq_along(u)) g$addEdge(u[[i]], v[[i]], w[[i]])
            }
        }

        seeds <- as.integer(seeds)
        memb_list <- vector("list", length(seeds))
        for (i in seq_along(seeds)) {
            if (!is.na(seeds[i])) try(nk$random$setSeed(seeds[i], TRUE), silent = TRUE)
            pl <- nk$community$ParallelLeiden(g,
                iterations = as.integer(iterations),
                randomize  = isTRUE(randomize),
                gamma      = as.numeric(gamma)
            )
            pl$run()
            part <- pl$getPartition()
            mv <- as.integer(unlist(part$getVector())) + 1L
            memb_list[[i]] <- mv
        }
        mm <- do.call(cbind, memb_list)
        rownames(mm) <- vertex_names
        colnames(mm) <- paste0("r", seq_along(seeds))
        return(mm)
    }

    g_local <- igraph::graph_from_data_frame(cbind(edge_table, weight = edge_weights), directed = FALSE, vertices = vertex_names)
    memb_list <- lapply(seq_along(seeds), function(i) {
        set.seed(as.integer(seeds[i]))
        res <- try(igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight, resolution = as.numeric(gamma)), silent = TRUE)
        if (inherits(res, "try-error")) res <- try(igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight, resolution_parameter = as.numeric(gamma)), silent = TRUE)
        if (inherits(res, "try-error")) res <- igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight)
        memb <- igraph::membership(res)
        if (is.null(names(memb))) names(memb) <- igraph::V(g_local)$name
        as.integer(memb[as.character(vertex_names)])
    })
    mm <- do.call(cbind, memb_list)
    rownames(mm) <- vertex_names
    colnames(mm) <- paste0("r", seq_along(seeds))
    mm
}

#' Run NetworKit PLM across multiple seeds (igraph Louvain fallback when absent).
#' @description
#' Internal helper for `.nk_plm_mruns`.
#' @keywords internal
#' @param vertex_names Character vector of vertex identifiers.
#' @param edge_table Data frame of edges with columns `from` and `to`.
#' @param edge_weights Numeric vector of edge weights aligned to `edge_table`.
#' @param gamma Resolution parameter.
#' @param refine Logical; whether to enable PLM refinement.
#' @param threads Integer thread count for NetworKit.
#' @param seeds Integer vector of seeds; one membership column is returned per seed.
#' @return Integer matrix of memberships (rows = vertices, cols = runs).
.nk_plm_mruns <- function(vertex_names, edge_table, edge_weights, gamma = 1.0,
                          refine = TRUE, threads = 1L, seeds = 1L) {
    if (.nk_available()) {
        nk <- reticulate::import("networkit", delay_load = TRUE)
        threads <- as.integer(threads)
        if (!is.na(threads) && threads > 0) try(nk$setNumberOfThreads(threads), silent = TRUE)

        vcount <- length(vertex_names)
        g <- nk$Graph(as.integer(vcount), weighted = TRUE, directed = FALSE)
        if (!is.null(edge_table) && nrow(edge_table)) {
            vid <- seq_len(vcount) - 1L
            map <- structure(vid, names = as.character(vertex_names))
            u <- map[as.character(edge_table$from)]
            v <- map[as.character(edge_table$to)]
            w <- as.numeric(edge_weights)
            keep <- !is.na(u) & !is.na(v) & is.finite(w)
            if (any(keep)) {
                u <- as.integer(u[keep]); v <- as.integer(v[keep]); w <- as.numeric(w[keep])
                for (i in seq_along(u)) g$addEdge(u[[i]], v[[i]], w[[i]])
            }
        }

        seeds <- as.integer(seeds)
        memb_list <- vector("list", length(seeds))
        for (i in seq_along(seeds)) {
            if (!is.na(seeds[i])) try(nk$random$setSeed(seeds[i], TRUE), silent = TRUE)
            algo <- nk$community$PLM(g, refine = isTRUE(refine), gamma = as.numeric(gamma))
            algo$run()
            part <- algo$getPartition()
            mv <- as.integer(unlist(part$getVector())) + 1L
            memb_list[[i]] <- mv
        }
        mm <- do.call(cbind, memb_list)
        rownames(mm) <- vertex_names
        colnames(mm) <- paste0("r", seq_along(seeds))
        return(mm)
    }

    g_local <- igraph::graph_from_data_frame(cbind(edge_table, weight = edge_weights), directed = FALSE, vertices = vertex_names)
    memb_list <- lapply(seq_along(seeds), function(i) {
        set.seed(as.integer(seeds[i]))
        res <- try(igraph::cluster_louvain(g_local, weights = igraph::E(g_local)$weight, resolution = as.numeric(gamma)), silent = TRUE)
        if (inherits(res, "try-error")) res <- igraph::cluster_louvain(g_local, weights = igraph::E(g_local)$weight)
        memb <- igraph::membership(res)
        if (is.null(names(memb))) names(memb) <- igraph::V(g_local)$name
        as.integer(memb[as.character(vertex_names)])
    })
    mm <- do.call(cbind, memb_list)
    rownames(mm) <- vertex_names
    colnames(mm) <- paste0("r", seq_along(seeds))
    mm
}

#' Community Detect
#' @description
#' Internal helper for `.community_detect`.
#' @param graph Parameter value.
#' @param algo Parameter value.
#' @param resolution Parameter value.
#' @param objective Parameter value.
#' @param backend Parameter value.
#' @param threads Number of threads to use.
#' @param debug Parameter value.
#' @return Return value used internally.
#' @keywords internal
.community_detect <- function(graph, algo = c("leiden", "louvain"), resolution = 1,
                               objective = c("CPM", "modularity"),
                               backend = c("igraph", "networkit"),
                               threads = 1, debug = FALSE) {
    algo <- match.arg(algo)
    objective <- match.arg(objective)
    backend <- match.arg(backend)

    if (algo == "leiden" && identical(backend, "networkit")) {
        memb <- try(.community_detect_networkit(graph, resolution = resolution,
                                               threads = threads, debug = debug), silent = TRUE)
        if (!inherits(memb, "try-error")) return(memb)
        if (debug) {
            .log_info("clusterGenes", "S05", "networkit Leiden failed; falling back to igraph implementation", debug)
        }
        backend <- "igraph"
    }

    if (identical(backend, "igraph")) {
        if (algo == "leiden") {
            comm <- try(igraph::cluster_leiden(graph, weights = igraph::E(graph)$weight,
                                               resolution = resolution, objective_function = objective), silent = TRUE)
            if (inherits(comm, "try-error")) {
                comm <- try(igraph::cluster_leiden(graph, weights = igraph::E(graph)$weight,
                                                   resolution_parameter = resolution, objective_function = objective), silent = TRUE)
            }
            if (inherits(comm, "try-error")) {
                comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight,
                                                     resolution = resolution, objective_function = objective), silent = TRUE)
                if (inherits(comm, "try-error")) {
                    comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight,
                                                         resolution_parameter = resolution, objective_function = objective), silent = TRUE)
                    if (inherits(comm, "try-error")) comm <- igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight)
                }
            }
            memb <- igraph::membership(comm)
            if (is.null(names(memb))) names(memb) <- igraph::as_ids(igraph::V(graph))
            return(memb)
        } else {
            comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight,
                                                 resolution = resolution, objective_function = objective), silent = TRUE)
            if (inherits(comm, "try-error")) {
                comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight,
                                                     resolution_parameter = resolution, objective_function = objective), silent = TRUE)
                if (inherits(comm, "try-error")) comm <- igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight)
            }
            memb <- igraph::membership(comm)
            if (is.null(names(memb))) names(memb) <- igraph::as_ids(igraph::V(graph))
            return(memb)
        }
    }

    stop("Unsupported backend", call. = FALSE)
}

#' Community Detect Networkit
#' @description
#' Internal helper for `.community_detect_networkit`.
#' @param graph Parameter value.
#' @param resolution Parameter value.
#' @param threads Number of threads to use.
#' @param debug Parameter value.
#' @return Return value used internally.
#' @keywords internal
.community_detect_networkit <- function(graph, resolution = 1, threads = 1, debug = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("reticulate package required for networkit backend", call. = FALSE)
    }

    edges <- igraph::as_data_frame(graph, what = "edges")
    node_names <- igraph::as_ids(igraph::V(graph))
    if (is.null(node_names)) node_names <- as.character(seq_len(igraph::vcount(graph)))

    if (!nrow(edges)) {
        memb <- rep(1L, igraph::vcount(graph))
        names(memb) <- node_names
        return(memb)
    }

    if (!"weight" %in% names(edges)) edges$weight <- 1
    edges$weight[is.na(edges$weight)] <- 1

    idx_map <- setNames(seq_len(igraph::vcount(graph)) - 1L, node_names)
    edges$from_idx <- idx_map[as.character(edges$from)]
    edges$to_idx   <- idx_map[as.character(edges$to)]

    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    cols <- edges[, c("from_idx", "to_idx", "weight")]
    if (requireNamespace("data.table", quietly = TRUE)) {
        fwrite(cols, tmp, sep = " ", col.names = FALSE)
    } else {
        write.table(cols, file = tmp, sep = " ", col.names = FALSE,
                    row.names = FALSE, quote = FALSE)
    }

    nk <- reticulate::import("networkit", delay_load = TRUE)
    threads_use <- suppressWarnings(as.integer(threads))
    if (length(threads_use) == 0L || is.na(threads_use) || threads_use < 1L) {
        threads_use <- 1L
    }
    try(nk$setNumberOfThreads(threads_use), silent = TRUE)
    reader <- tryCatch(
        nk$graphio$EdgeListReader(" ", 0L, "#", TRUE, FALSE),
        error = function(e_new) {
            tryCatch(
                nk$graphio$EdgeListReader(weighted = TRUE, directed = FALSE),
                error = function(e_old) stop(e_old)
            )
        }
    )
    G <- reader$read(tmp)
    if (debug) {
        .log_info("clusterGenes", "S05", sprintf("networkit ParallelLeiden threads=%d resolution=%.4f",
            threads_use, resolution), debug)
    }
    # NetworKit >=11 uses `gamma` argument; older releases still expect `resolution`.
    create_algo <- function(use_gamma = TRUE) {
        res <- as.numeric(resolution)
        if (use_gamma) {
            nk$community$ParallelLeiden(G, gamma = res)
        } else {
            nk$community$ParallelLeiden(G, resolution = res)
        }
    }
    algo <- tryCatch(
        create_algo(TRUE),
        error = function(e_gamma) {
            is_py_type_error <- inherits(e_gamma, "python.builtin.TypeError") ||
                inherits(e_gamma, "reticulate.python.builtin.type_error")
            has_gamma_kw <- grepl("unexpected keyword argument 'gamma'", conditionMessage(e_gamma), fixed = TRUE)
            has_cinit_arity <- grepl("__cinit__", conditionMessage(e_gamma), fixed = TRUE)
            if (is_py_type_error || has_gamma_kw || has_cinit_arity) {
                tryCatch(
                    create_algo(FALSE),
                    error = function(e_resolution) {
                        stop(e_resolution)
                    }
                )
            } else {
                stop(e_gamma)
            }
        }
    )
    algo$run()
    comm <- algo$getPartition()
    membership <- reticulate::py_to_r(comm$getVector())
    membership <- as.integer(membership) + 1L
    names(membership) <- node_names
    membership
}

#' Threshold Similarity Graph
#' @description
#' Internal helper for `.threshold_similarity_graph`.
#' @param similarity_matrix Parameter value.
#' @param significance_matrix Parameter value.
#' @param config Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.threshold_similarity_graph <- function(similarity_matrix,
                                  significance_matrix = NULL,
                                  config,
                                  verbose = TRUE) {
    .cluster_message(verbose, "[cluster] Step1: thresholding similarity matrix...",
        parent = "clusterGenes", step = "S05"
    )

    genes_all <- rownames(similarity_matrix)
    if (is.null(genes_all)) stop("Input similarity matrix must have rownames/colnames for node IDs.")

    similarity_filtered <- similarity_matrix
    similarity_filtered[abs(similarity_filtered) < config$min_cutoff] <- 0
    if (config$use_significance && !is.null(significance_matrix)) {
        similarity_filtered <- tryCatch(
            .align_and_filter_fdr(similarity_filtered, similarity_matrix, significance_matrix, config$significance_max),
            error = function(e) {
                warning("[cluster] significance filter disabled: ", conditionMessage(e))
                similarity_filtered
            }
        )
    }
    similarity_filtered[similarity_filtered < 0] <- 0
    similarity_filtered <- pmax(similarity_filtered, t(similarity_filtered))
    diag(similarity_filtered) <- 0
    similarity_filtered <- .filter_matrix_by_quantile(similarity_filtered, config$pct_min, "q100")
    if (!inherits(similarity_filtered, "sparseMatrix")) {
        similarity_filtered <- drop0(Matrix(similarity_filtered, sparse = TRUE))
    }

    kept_genes <- genes_all[
        if (config$drop_isolated) rowSums(similarity_filtered != 0) > 0 else rep(TRUE, nrow(similarity_filtered))
    ]
    if (length(kept_genes) < 2) stop("Fewer than 2 genes after filtering.")

    similarity_kept <- similarity_filtered[kept_genes, kept_genes, drop = FALSE]
    .cluster_message(verbose, sprintf("[cluster]   retained %d/%d nodes; |E|=%d",
                        length(kept_genes), length(genes_all), nnzero(similarity_kept) / 2),
        parent = "clusterGenes", step = "S05"
    )

    list(
        similarity = similarity_filtered,
        similarity_sub = similarity_kept,
        kept_genes = kept_genes,
        genes_all = genes_all,
        A = similarity_filtered,                 # backward compatibility
        A_sub = similarity_kept           # backward compatibility
    )
}

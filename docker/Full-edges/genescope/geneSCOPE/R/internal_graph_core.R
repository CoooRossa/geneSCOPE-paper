#' Networkit Available
#' @description
#' Internal helper for `.networkit_available`.
#' @return Return value used internally.
#' @keywords internal
.networkit_available <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
    if (!isTRUE(get0("configured", envir = .networkit_runtime_state, ifnotfound = FALSE))) {
        .autoconfigure_networkit_runtime()
        assign("configured", TRUE, envir = .networkit_runtime_state)
    }
    ok <- tryCatch({ reticulate::py_module_available("networkit") }, error = function(e) FALSE)
    isTRUE(ok)
}

#' Networkit Edge Reader
#' @description
#' Internal helper for `.networkit_edge_reader`.
#' @param nk Parameter value.
#' @return Return value used internally.
#' @keywords internal
.networkit_edge_reader <- function(nk) {
    tryCatch(
        nk$graphio$EdgeListReader(" ", 0L, "#", TRUE, FALSE),
        error = function(e_new) {
            tryCatch(
                nk$graphio$EdgeListReader(weighted = TRUE, directed = FALSE),
                error = function(e_old) stop(e_old)
            )
        }
    )
}

#' Networkit Load Graph From File
#' @description
#' Internal helper for `.networkit_load_graph_from_file`.
#' @param nk Parameter value.
#' @param edge_file Filesystem path.
#' @param vertex_count Parameter value.
#' @return Return value used internally.
#' @keywords internal
.networkit_load_graph_from_file <- function(nk, edge_file, vertex_count) {
    vertex_count <- as.integer(vertex_count)
    if (is.null(edge_file) || !nzchar(edge_file) || !file.exists(edge_file)) {
        return(nk$Graph(vertex_count, weighted = TRUE, directed = FALSE))
    }

    nk_graph <- .networkit_edge_reader(nk)$read(edge_file)
    existing_nodes <- tryCatch(nk_graph$numberOfNodes(), error = function(e) NA_integer_)
    if (is.finite(existing_nodes) && vertex_count > existing_nodes) {
        nk_graph$addNodes(as.integer(vertex_count - existing_nodes))
    }
    nk_graph
}

#' Networkit Parallel Leiden Runs
#' @description
#' Internal helper for `.networkit_parallel_leiden_runs`.
#' @param vertex_names Parameter value.
#' @param edge_table Parameter value.
#' @param edge_weights Parameter value.
#' @param gamma Parameter value.
#' @param iterations Parameter value.
#' @param randomize Parameter value.
#' @param threads Number of threads to use.
#' @param seeds Integer vector of random seeds.
#' @return Return value used internally.
#' @keywords internal
.networkit_parallel_leiden_runs <- function(vertex_names, edge_table, edge_weights,
                                            gamma = 1.0, iterations = 10L,
                                            randomize = TRUE,
                                            threads = 1L, seeds = 1L) {
    if (!.networkit_available()) {
        g_local <- igraph::graph_from_data_frame(cbind(edge_table, weight = edge_weights), directed = FALSE, vertices = vertex_names)
        membership_list <- lapply(seeds, function(seed_id) {
            res <- try(igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight, resolution_parameter = gamma), silent = TRUE)
            if (inherits(res, "try-error")) res <- try(igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight, resolution = gamma), silent = TRUE)
            if (inherits(res, "try-error")) res <- igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight)
            memberships <- igraph::membership(res)
            if (is.null(names(memberships))) names(memberships) <- igraph::V(g_local)$name
            as.integer(memberships[vertex_names])
        })
        membership_matrix <- do.call(cbind, membership_list)
        rownames(membership_matrix) <- vertex_names
        return(membership_matrix)
    }

    nk <- reticulate::import("networkit", delay_load = TRUE)
    threads <- as.integer(threads)
    if (!is.na(threads) && threads > 0) try(nk$setNumberOfThreads(threads), silent = TRUE)

    edge_ctx <- .networkit_write_edge_file(vertex_names, edge_table, edge_weights)
    on.exit(edge_ctx$cleanup(), add = TRUE)
    nk_graph <- .networkit_load_graph_from_file(nk, edge_ctx$file, edge_ctx$vertex_count)

    seeds <- as.integer(seeds)
    membership_list <- vector("list", length(seeds))
    for (i in seq_along(seeds)) {
        if (!is.na(seeds[i])) try(nk$random$setSeed(seeds[i], TRUE), silent = TRUE)
        pl <- nk$community$ParallelLeiden(nk_graph,
            iterations = as.integer(iterations),
            randomize = as.logical(randomize),
            gamma = as.numeric(gamma))
        pl$run()
        partition <- pl$getPartition()
        memberships_zero_based <- as.integer(unlist(partition$getVector()))
        membership_list[[i]] <- memberships_zero_based + 1L
    }
    membership_matrix <- do.call(cbind, membership_list)
    rownames(membership_matrix) <- vertex_names
    colnames(membership_matrix) <- paste0("r", seq_along(seeds))
    membership_matrix
}

#' Networkit Plm Runs
#' @description
#' Internal helper for `.networkit_plm_runs`.
#' @param vertex_names Parameter value.
#' @param edge_table Parameter value.
#' @param edge_weights Parameter value.
#' @param gamma Parameter value.
#' @param refine Parameter value.
#' @param threads Number of threads to use.
#' @param seeds Integer vector of random seeds.
#' @return Return value used internally.
#' @keywords internal
.networkit_plm_runs <- function(vertex_names, edge_table, edge_weights,
                                gamma = 1.0, refine = TRUE,
                                threads = 1L, seeds = 1L) {
    if (!.networkit_available()) {
        g_local <- igraph::graph_from_data_frame(cbind(edge_table, weight = edge_weights), directed = FALSE, vertices = vertex_names)
        membership_list <- lapply(seeds, function(seed_id) {
            res <- try(igraph::cluster_louvain(g_local, weights = igraph::E(g_local)$weight, resolution = gamma), silent = TRUE)
            if (inherits(res, "try-error")) res <- igraph::cluster_louvain(g_local, weights = igraph::E(g_local)$weight)
            memberships <- igraph::membership(res)
            if (is.null(names(memberships))) names(memberships) <- igraph::V(g_local)$name
            as.integer(memberships[vertex_names])
        })
        membership_matrix <- do.call(cbind, membership_list)
        rownames(membership_matrix) <- vertex_names
        return(membership_matrix)
    }

    nk <- reticulate::import("networkit", delay_load = TRUE)
    threads <- as.integer(threads)
    if (!is.na(threads) && threads > 0) try(nk$setNumberOfThreads(threads), silent = TRUE)

    edge_ctx <- .networkit_write_edge_file(vertex_names, edge_table, edge_weights)
    on.exit(edge_ctx$cleanup(), add = TRUE)
    nk_graph <- .networkit_load_graph_from_file(nk, edge_ctx$file, edge_ctx$vertex_count)

    seeds <- as.integer(seeds)
    membership_list <- vector("list", length(seeds))
    for (i in seq_along(seeds)) {
        if (!is.na(seeds[i])) try(nk$random$setSeed(seeds[i], TRUE), silent = TRUE)
        plm_algo <- nk$community$PLM(nk_graph, refine = as.logical(refine), gamma = as.numeric(gamma))
        plm_algo$run()
        partition <- plm_algo$getPartition()
        memberships_zero_based <- as.integer(unlist(partition$getVector()))
        membership_list[[i]] <- memberships_zero_based + 1L
    }
    membership_matrix <- do.call(cbind, membership_list)
    rownames(membership_matrix) <- vertex_names
    colnames(membership_matrix) <- paste0("r", seq_along(seeds))
    membership_matrix
}

#' Networkit Write Edge File
#' @description
#' Internal helper for `.networkit_write_edge_file`.
#' @param vertex_names Parameter value.
#' @param edge_table Parameter value.
#' @param edge_weights Parameter value.
#' @param tmpdir Parameter value.
#' @return Return value used internally.
#' @keywords internal
.networkit_write_edge_file <- function(vertex_names, edge_table, edge_weights, tmpdir = NULL) {
    vertex_count <- length(vertex_names)
    result <- list(
        file = NULL,
        edge_count = 0L,
        vertex_count = vertex_count,
        cleanup = function() invisible(NULL)
    )
    if (!nrow(edge_table) || !vertex_count) return(result)

    vertex_ids <- seq_len(vertex_count) - 1L
    u <- structure(vertex_ids, names = vertex_names)[as.character(edge_table$from)]
    v <- structure(vertex_ids, names = vertex_names)[as.character(edge_table$to)]
    w <- as.numeric(edge_weights)
    keep <- !is.na(u) & !is.na(v) & is.finite(w)
    if (!any(keep)) return(result)

    edge_df <- data.frame(
        from = as.integer(u[keep]),
        to = as.integer(v[keep]),
        weight = as.numeric(w[keep]),
        stringsAsFactors = FALSE
    )

    tmpdir <- if (is.null(tmpdir) || !nzchar(tmpdir)) {
        getOption(
            "networkit.edge.tmpdir",
            default = getOption("networkit.edge.tmpdir_2", default = tempdir())
        )
    } else {
        tmpdir
    }
    if (!dir.exists(tmpdir)) {
        dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
    }

    final_path <- tempfile(pattern = "nk_edges_", tmpdir = tmpdir, fileext = ".tsv")
    staging_path <- paste0(final_path, ".tmp")

    write_ok <- FALSE
    writer <- function(dest) {
        if (requireNamespace("data.table", quietly = TRUE)) {
            fwrite(edge_df, dest, sep = " ", col.names = FALSE)
        } else {
            write.table(edge_df, file = dest, sep = " ", col.names = FALSE,
                               row.names = FALSE, quote = FALSE)
        }
    }

    tryCatch({
        writer(staging_path)
        write_ok <- TRUE
    }, error = function(e) {
        write_ok <<- FALSE
        if (file.exists(staging_path)) unlink(staging_path, recursive = FALSE, force = FALSE)
        stop(e)
    })

    if (write_ok) {
        if (!file.rename(staging_path, final_path)) {
            final_path <- staging_path
        }
    } else {
        final_path <- staging_path
    }

    cleanup_paths <- unique(c(final_path, staging_path))
    result$file <- final_path
    result$edge_count <- nrow(edge_df)
    result$cleanup <- function() {
        existing <- cleanup_paths[file.exists(cleanup_paths)]
        if (length(existing)) unlink(existing, recursive = FALSE, force = FALSE)
        invisible(NULL)
    }
    result
}

#' Nk Prepare Graph Input
#' @description
#' Internal helper for `.nk_prepare_graph_input`.
#' @param graph Parameter value.
#' @return Return value used internally.
#' @keywords internal
.nk_prepare_graph_input <- function(graph) {
    .cache <- tryCatch(igraph::graph_attr(graph, "_nk_cache"), error = function(e) NULL)
    if (is.list(.cache) &&
        !is.null(.cache$edge_table) &&
        !is.null(.cache$edge_weights) &&
        !is.null(.cache$vertex_names)) {
        return(list(
            graph = graph,
            edge_table = .cache$edge_table,
            edge_weights = .cache$edge_weights,
            vertex_names = .cache$vertex_names
        ))
    }

    vertex_names <- igraph::V(graph)$name
    edge_df <- igraph::as_data_frame(graph, what = "edges")
    if (nrow(edge_df)) {
        edge_table <- edge_df[, c("from", "to"), drop = FALSE]
    } else {
        edge_table <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
    }
    weights <- igraph::E(graph)$weight
    if (is.null(weights)) {
        weights <- if (nrow(edge_table)) rep(1, nrow(edge_table)) else numeric(0)
    }
    .cache <- list(
        edge_table = edge_table,
        edge_weights = weights,
        vertex_names = vertex_names
    )
    graph <- igraph::set_graph_attr(graph, "_nk_cache", .cache)
    list(
        graph = graph,
        edge_table = edge_table,
        edge_weights = weights,
        vertex_names = vertex_names
    )
}

#' Run Graph Algorithm
#' @description
#' Internal helper for `.run_graph_algorithm`.
#' @param g Parameter value.
#' @param algo_name Parameter value.
#' @param res_param Parameter value.
#' @param objective Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_graph_algorithm <- function(g, algo_name, res_param, objective) {
    w <- igraph::E(g)$weight
    if (identical(algo_name, "leiden")) {
        res <- try(igraph::cluster_leiden(
            g,
            weights = w,
            resolution_parameter = res_param,
            objective_function = objective
        ), silent = TRUE)
        if (inherits(res, "try-error")) {
            res <- try(igraph::cluster_leiden(
                g,
                weights = w,
                resolution = res_param
            ), silent = TRUE)
            if (inherits(res, "try-error")) res <- igraph::cluster_leiden(g, weights = w)
        }
        return(res)
    }
    formn <- try(names(formals(igraph::cluster_louvain)), silent = TRUE)
    if (!inherits(formn, "try-error") && "resolution" %in% formn) {
        return(igraph::cluster_louvain(g, weights = w, resolution = res_param))
    }
    igraph::cluster_louvain(g, weights = w)
}

.nk_available <- .networkit_available
.nk_parallel_leiden_mruns <- .networkit_parallel_leiden_runs
.nk_plm_mruns <- .networkit_plm_runs

#' Consensus co-occurrence (sparse COO via C++)
#' @description
#'   Works in two modes:
#'   - Package mode: calls registered symbol `_geneSCOPE_consensus_coo_cpp`.
#'   - Sourced mode: if you ran `Rcpp::sourceCpp('src/9.ConsensusAccel.cpp')`,
#'     calls the Rcpp-exposed function `.consensus_coo_cpp()` directly.
#' @param memb Integer matrix (genes × runs), each column a run's cluster labels
#' @param thr Numeric in [0,1]. Threshold on fraction across runs (default 0)
#' @param n_threads Integer. Threads (ignored if OpenMP not available)
#' @return list(i,j,x,n) with 1-based COO
#' @keywords internal
consensus_coo <- function(memb, thr = 0, n_threads = NULL) {
    if (is.null(n_threads)) n_threads <- .safe_thread_count()
    # Try package-registered symbol
    call_sym_ok <- FALSE
    res <- try({
        .Call(`_geneSCOPE_consensus_coo_cpp`, memb, thr, as.integer(n_threads), PACKAGE = "geneSCOPE")
    }, silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    # Try direct Rcpp function if sourceCpp() was used
    if (exists(".consensus_coo_cpp", mode = "function", inherits = TRUE)) {
        return(.consensus_coo_cpp(memb, thr, as.integer(n_threads)))
    }
    stop("consensus_coo: C++ backend not found. Install/load geneSCOPE or run Rcpp::sourceCpp('genescope/src/9.ConsensusAccel.cpp').")
}

#' Consensus sparse matrix (dgCMatrix)
#' @description
#'   Convenience wrapper building a symmetric sparse matrix with edge weights
#'   equal to fraction of co-occurrence across runs for pairs meeting `thr`.
#' @param memb Parameter value.
#' @param thr Parameter value.
#' @param n_threads Number of threads to use.
#' @return A symmetric `dgCMatrix` with zero diagonal
#' @keywords internal
consensus_sparse <- function(memb, thr = 0, n_threads = NULL) {
    if (is.null(n_threads)) n_threads <- .safe_thread_count()
    coo <- consensus_coo(memb, thr = thr, n_threads = n_threads)
    # Build symmetric sparse matrix
    M <- if (length(coo$i) == 0L) {
        sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(coo$n, coo$n), symmetric = TRUE)
    } else {
        sparseMatrix(i = coo$i, j = coo$j, x = coo$x, dims = c(coo$n, coo$n), symmetric = TRUE)
    }
    # Propagate dimnames from membership matrix when available to enable name-based indexing
    rn <- try(if (is.matrix(memb) || inherits(memb, "Matrix")) rownames(memb) else rownames(as.matrix(memb)), silent = TRUE)
    if (!inherits(rn, "try-error") && !is.null(rn) && length(rn) == nrow(M)) {
        # Guard against rare cases where dimnames<- is applied to non-arrays
        try({
            if (is.matrix(M) || inherits(M, "Matrix")) dimnames(M) <- list(rn, rn)
        }, silent = TRUE)
    }
    M
}

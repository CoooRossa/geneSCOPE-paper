#' Validate STRINGdb availability.
#' @keywords internal
.stringdb_require_or_stop <- function() {
    if (identical(Sys.getenv("GENESCOPE_DISABLE_STRINGDB"), "1")) {
        stop(
            "STRINGdb is required for STRING evaluation. Install with:\n",
            "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
            "  BiocManager::install('STRINGdb')"
        )
    }
    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
        stop(
            "STRINGdb is required for STRING evaluation. Install with:\n",
            "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
            "  BiocManager::install('STRINGdb')"
        )
    }
    invisible(TRUE)
}

#' Read update flag from env/option.
#' @keywords internal
.stringdb_update_flag <- function() {
    opt <- getOption("geneSCOPE.stringdb_update", FALSE)
    if (isTRUE(opt)) return(TRUE)
    env <- Sys.getenv("GENESCOPE_STRINGDB_UPDATE", "")
    if (nzchar(env) && !identical(env, "0")) return(TRUE)
    FALSE
}

#' Download STRINGdb file with md5 tracking.
#' @keywords internal
.stringdb_download_with_md5 <- function(url, dest, force_update = FALSE) {
    dest_dir <- dirname(dest)
    if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
    md5_path <- paste0(dest, ".md5")

    if (!force_update && file.exists(dest) && file.info(dest)$size > 0) {
        if (file.exists(md5_path)) return(invisible(TRUE))
        current <- tryCatch(tools::md5sum(dest), error = function(e) NA_character_)
        current <- unname(current)
        if (!is.na(current)) writeLines(current, md5_path)
        return(invisible(TRUE))
    }

    tmp <- paste0(dest, ".tmp")
    if (file.exists(tmp)) unlink(tmp)

    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    if (!is.finite(old_timeout) || old_timeout < 600) options(timeout = 600)

    last_err <- NULL
    download_ok <- FALSE
    for (attempt in seq_len(3L)) {
        warn_msgs <- character()
        ok <- tryCatch({
            withCallingHandlers(
                utils::download.file(url, tmp, mode = "wb", quiet = TRUE),
                warning = function(w) {
                    warn_msgs <<- c(warn_msgs, conditionMessage(w))
                    invokeRestart("muffleWarning")
                }
            )
            TRUE
        }, error = function(e) {
            last_err <<- e
            FALSE
        })

        bad_warn <- any(grepl("downloaded length|timeout", warn_msgs, ignore.case = TRUE))
        if (isTRUE(ok) && !bad_warn && file.exists(tmp) && file.info(tmp)$size > 0) {
            download_ok <- TRUE
            break
        }
        if (file.exists(tmp)) unlink(tmp)
        Sys.sleep(attempt)
    }

    if (!download_ok || !file.exists(tmp) || file.info(tmp)$size == 0) {
        if (file.exists(tmp)) unlink(tmp)
        msg <- if (!is.null(last_err)) conditionMessage(last_err) else "download failed"
        stop("Failed to download STRINGdb data: ", basename(dest), " (", msg, ")")
    }

    new_md5 <- unname(tryCatch(tools::md5sum(tmp), error = function(e) NA_character_))

    if (isTRUE(force_update) && file.exists(dest)) {
        old_md5 <- NA_character_
        if (file.exists(md5_path)) {
            old_md5 <- readLines(md5_path, warn = FALSE)
            old_md5 <- if (length(old_md5)) trimws(old_md5[1]) else NA_character_
        } else {
            old_md5 <- unname(tryCatch(tools::md5sum(dest), error = function(e) NA_character_))
        }
        if (!is.na(old_md5) && !is.na(new_md5) && nzchar(old_md5) && identical(old_md5, new_md5)) {
            unlink(tmp)
            if (!file.exists(md5_path) && !is.na(new_md5)) writeLines(new_md5, md5_path)
            return(invisible(TRUE))
        }
    }

    ok_rename <- file.rename(tmp, dest)
    if (!isTRUE(ok_rename)) {
        file.copy(tmp, dest, overwrite = TRUE)
        unlink(tmp)
    }
    if (!is.na(new_md5)) writeLines(new_md5, md5_path)
    invisible(TRUE)
}

#' Prepare STRINGdb download cache with md5 checks.
#' @keywords internal
.stringdb_prepare_cache_files <- function(string_db, cache_dir, force_update = NULL) {
    if (is.null(cache_dir) || !nzchar(cache_dir)) return(invisible(FALSE))
    if (is.null(force_update)) force_update <- .stringdb_update_flag()

    protocol <- string_db$protocol
    species <- string_db$species
    file_version <- string_db$file_version

    aliases_base <- paste0(species, ".protein.aliases.v", file_version, ".txt.gz")
    info_base <- paste0(species, ".protein.info.v", file_version, ".txt.gz")

    aliases_url <- paste0(
        protocol,
        "://stringdb-downloads.org/download/protein.aliases.v",
        file_version,
        "/",
        aliases_base
    )
    info_url <- paste0(
        protocol,
        "://stringdb-downloads.org/download/protein.info.v",
        file_version,
        "/",
        info_base
    )

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
    links_base <- paste0(
        species,
        ".protein.",
        network_type_param,
        link_data_param,
        file_version,
        ".txt.gz"
    )
    links_url <- paste0(
        protocol,
        "://stringdb-downloads.org/download/protein.",
        network_type_param,
        link_data_param,
        file_version,
        "/",
        links_base
    )

    .stringdb_download_with_md5(aliases_url, file.path(cache_dir, aliases_base), force_update = force_update)
    .stringdb_download_with_md5(info_url, file.path(cache_dir, info_base), force_update = force_update)
    .stringdb_download_with_md5(links_url, file.path(cache_dir, links_base), force_update = force_update)
    invisible(TRUE)
}

#' Create STRINGdb instance.
#' @keywords internal
.stringdb_connect <- function(species = 9606,
                              score_threshold = 0,
                              version = "11.5",
                              cache_dir = NULL) {
    .stringdb_require_or_stop()
    cache_dir_use <- cache_dir
    if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
        cache_dir_use <- Sys.getenv("GENESCOPE_STRINGDB_DIR", "")
    }
    if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
        cache_dir_use <- getOption("geneSCOPE.stringdb_dir", "")
    }
    if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
        cache_dir_use <- tryCatch(
            tools::R_user_dir("geneSCOPE", "cache"),
            error = function(e) ""
        )
    }
    if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
        return(STRINGdb::STRINGdb$new(
            version = version,
            species = species,
            score_threshold = score_threshold
        ))
    }
    if (!dir.exists(cache_dir_use)) {
        ok <- tryCatch({
            dir.create(cache_dir_use, recursive = TRUE)
            TRUE
        }, error = function(e) FALSE)
        if (!isTRUE(ok)) {
            return(STRINGdb::STRINGdb$new(
                version = version,
                species = species,
                score_threshold = score_threshold
            ))
        }
    }
    string_db <- STRINGdb::STRINGdb$new(
        version = version,
        species = species,
        score_threshold = score_threshold,
        input_directory = cache_dir_use
    )
    .stringdb_prepare_cache_files(string_db, cache_dir_use)
    string_db
}

#' Map genes to STRING IDs with optional caching.
#' @keywords internal
.stringdb_map_genes <- function(string_db, genes, cache_dir = NULL) {
    genes <- unique(as.character(genes))
    if (!length(genes)) {
        return(setNames(character(0), character(0)))
    }
    cache_path <- NULL
    if (!is.null(cache_dir) && requireNamespace("digest", quietly = TRUE)) {
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
        key <- digest::digest(list(genes = sort(genes), species = string_db$species))
        cache_path <- file.path(cache_dir, paste0("string_mapping_", key, ".rds"))
        if (file.exists(cache_path)) {
            return(readRDS(cache_path))
        }
    }
    df <- data.frame(gene = genes, stringsAsFactors = FALSE)
    mapped <- string_db$map(df, "gene", removeUnmappedRows = FALSE)
    mapping <- setNames(mapped$STRING_id, mapped$gene)
    if (!is.null(cache_path)) saveRDS(mapping, cache_path)
    mapping
}

#' Fetch STRING interactions for provided STRING IDs with optional caching.
#' @keywords internal
.stringdb_get_interactions <- function(string_db, string_ids, cache_dir = NULL) {
    ids <- unique(na.omit(as.character(string_ids)))
    if (!length(ids)) {
        return(data.frame(from = character(0), to = character(0), combined_score = numeric(0)))
    }
    cache_path <- NULL
    if (!is.null(cache_dir) && requireNamespace("digest", quietly = TRUE)) {
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
        key <- digest::digest(list(ids = sort(ids), species = string_db$species))
        cache_path <- file.path(cache_dir, paste0("string_interactions_", key, ".rds"))
        if (file.exists(cache_path)) {
            return(readRDS(cache_path))
        }
    }
    interactions <- string_db$get_interactions(ids)
    if (!is.data.frame(interactions) || !nrow(interactions)) {
        interactions <- data.frame(from = character(0), to = character(0), combined_score = numeric(0))
    } else {
        interactions <- interactions[, c("from", "to", "combined_score"), drop = FALSE]
    }
    if (!is.null(cache_path)) saveRDS(interactions, cache_path)
    interactions
}

#' Annotate edges with STRING scores and labels.
#' @keywords internal
.stringdb_label_edges <- function(edge_df, gene_to_string, interactions, score_threshold = 700) {
    edge_df$string_from <- gene_to_string[edge_df$from]
    edge_df$string_to <- gene_to_string[edge_df$to]
    comparable <- !is.na(edge_df$string_from) & !is.na(edge_df$string_to)

    if (nrow(interactions)) {
        key_int <- paste(
            pmin(interactions$from, interactions$to),
            pmax(interactions$from, interactions$to),
            sep = "|"
        )
        score_map <- setNames(interactions$combined_score, key_int)
    } else {
        score_map <- setNames(numeric(0), character(0))
    }

    key_edges <- paste(
        pmin(edge_df$string_from, edge_df$string_to),
        pmax(edge_df$string_from, edge_df$string_to),
        sep = "|"
    )
    string_score <- score_map[key_edges]
    string_score[is.na(string_score)] <- 0
    string_score[!comparable] <- NA_real_

    edge_df$string_score <- as.numeric(string_score)
    edge_df$label <- ifelse(!is.na(edge_df$string_score), edge_df$string_score >= score_threshold, NA)
    edge_df
}

#' Compute AUPRC for a score ranking.
#' @keywords internal
.compute_auprc <- function(labels, scores) {
    ok <- !is.na(labels) & !is.na(scores)
    labels <- labels[ok]
    scores <- scores[ok]
    if (!length(labels)) return(NA_real_)
    labels <- as.integer(labels) > 0
    n_pos <- sum(labels)
    n_neg <- sum(!labels)
    if (n_pos == 0 || n_neg == 0) return(NA_real_)
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    tp <- cumsum(labels)
    fp <- cumsum(!labels)
    recall <- tp / n_pos
    precision <- tp / pmax(tp + fp, 1)
    recall <- c(0, recall)
    precision <- c(1, precision)
    sum(diff(recall) * precision[-1])
}

#' Compute AUROC for a score ranking.
#' @keywords internal
.compute_auroc <- function(labels, scores) {
    ok <- !is.na(labels) & !is.na(scores)
    labels <- labels[ok]
    scores <- scores[ok]
    if (!length(labels)) return(NA_real_)
    labels <- as.integer(labels) > 0
    n_pos <- sum(labels)
    n_neg <- sum(!labels)
    if (n_pos == 0 || n_neg == 0) return(NA_real_)
    ranks <- rank(scores, ties.method = "average")
    sum_ranks_pos <- sum(ranks[labels])
    (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' Precision at K for a score ranking.
#' @keywords internal
.precision_at_k <- function(labels, scores, k = 100L) {
    ok <- !is.na(labels) & !is.na(scores)
    labels <- labels[ok]
    scores <- scores[ok]
    if (!length(labels)) return(NA_real_)
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    k <- min(length(labels), as.integer(k))
    if (k <= 0L) return(NA_real_)
    mean(as.integer(labels[seq_len(k)]))
}

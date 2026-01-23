#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
  library(STRINGdb)
})

options(stringsAsFactors = FALSE)

split_csv <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(character(0))
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  trimws(parts[nzchar(parts)])
}

parse_int_list <- function(x) {
  vals <- suppressWarnings(as.integer(split_csv(x)))
  vals <- vals[is.finite(vals)]
  unique(vals)
}

normalize_string_version <- function(x) {
  x <- trimws(as.character(x))
  if (!nzchar(x) || is.na(x)) return(NA_character_)
  x
}

sha256_file <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return("")
  cmd <- Sys.which("sha256sum")
  if (nzchar(cmd)) {
    out <- tryCatch(system2(cmd, path, stdout = TRUE, stderr = TRUE), error = function(e) "")
    if (length(out) > 0) return(strsplit(out[1], "\\s+")[[1]][1])
  }
  cmd2 <- Sys.which("shasum")
  if (nzchar(cmd2)) {
    out <- tryCatch(system2(cmd2, c("-a", "256", path), stdout = TRUE, stderr = TRUE), error = function(e) "")
    if (length(out) > 0) return(strsplit(out[1], "\\s+")[[1]][1])
  }
  as.character(tools::md5sum(path))
}

get_first_repeat_dir <- function(bench_root, method) {
  method_dir <- file.path(bench_root, method)
  reps <- list.files(method_dir, pattern = "^repeat_", full.names = TRUE)
  if (length(reps) == 0) {
    stop("No repeat_* directories found for method: ", method)
  }
  reps <- sort(reps)
  reps[1]
}

file_ends_with_newline <- function(path) {
  info <- suppressWarnings(file.info(path))
  size <- suppressWarnings(as.numeric(info$size))
  if (!is.finite(size) || size < 1) return(TRUE)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = size - 1, origin = "start")
  last <- readBin(con, what = "raw", n = 1)
  if (!length(last)) return(TRUE)
  identical(as.integer(last), 10L)
}

safe_fread_edges_tsv <- function(path) {
  info <- suppressWarnings(file.info(path))
  size <- suppressWarnings(as.numeric(info$size))
  needs_pipe <- !file_ends_with_newline(path) && is.finite(size) && size > 0 && (size %% 4096 == 0)
  if (needs_pipe) {
    cmd <- paste("cat", shQuote(path), "; echo")
    return(fread(cmd = cmd, sep = "\t", fill = TRUE, showProgress = FALSE))
  }
  fread(path, sep = "\t", fill = TRUE, showProgress = FALSE)
}

read_edges_all <- function(path) {
  edges <- safe_fread_edges_tsv(path)
  required <- c("gene_a", "gene_b", "weight")
  missing <- setdiff(required, names(edges))
  if (length(missing)) {
    stop("Missing columns in edges file ", path, ": ", paste(missing, collapse = ", "))
  }
  if ("edge_type" %in% names(edges)) {
    edge_type_val <- as.character(edges[["edge_type"]])
    edges <- edges[!is.na(edge_type_val) & nzchar(edge_type_val)]
  }
  edges <- edges[!is.na(gene_a) & !is.na(gene_b)]
  edges <- edges[nzchar(as.character(gene_a)) & nzchar(as.character(gene_b))]
  has_fdr <- "fdr" %in% names(edges)
  edges <- edges[, .(
    gene_a = as.character(gene_a),
    gene_b = as.character(gene_b),
    weight = suppressWarnings(as.numeric(weight)),
    fdr = if (has_fdr) suppressWarnings(as.numeric(fdr)) else NA_real_
  )]
  edges <- edges[gene_a != gene_b]
  edges <- edges[is.finite(weight)]
  edges
}

.stringdb_update_flag <- function() {
  opt <- getOption("stringdb_update", FALSE)
  if (isTRUE(opt)) return(TRUE)
  env <- Sys.getenv("STRINGDB_UPDATE", "")
  if (nzchar(env) && !identical(env, "0")) return(TRUE)
  FALSE
}

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

  ok_rename <- file.rename(tmp, dest)
  if (!isTRUE(ok_rename)) {
    file.copy(tmp, dest, overwrite = TRUE)
    unlink(tmp)
  }
  if (!is.na(new_md5)) writeLines(new_md5, md5_path)
  invisible(TRUE)
}

.stringdb_safe_link_data <- function(string_db) {
  link_data <- tryCatch(string_db$link_data, error = function(e) NULL)
  if (is.null(link_data) || length(link_data) == 0 || is.na(link_data[[1]])) {
    return("links")
  }
  link_data <- as.character(link_data[[1]])
  if (!nzchar(link_data)) return("links")
  tolower(link_data)
}

.stringdb_interaction_paths <- function(string_db) {
  network_type_param <- ""
  if (tolower(string_db$network_type) == "physical") {
    network_type_param <- "physical."
  }
  link_data_param <- "links.v"
  link_data <- .stringdb_safe_link_data(string_db)
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

.stringdb_prepare_cache_files <- function(string_db, cache_dir, force_update = NULL) {
  if (is.null(cache_dir) || !nzchar(cache_dir)) return(invisible(FALSE))
  if (is.null(force_update)) force_update <- .stringdb_update_flag()

  protocol <- string_db$protocol
  species <- string_db$species
  file_version <- string_db$file_version

  aliases_base <- paste0(species, ".protein.aliases.v", file_version, ".txt.gz")
  info_base <- paste0(species, ".protein.info.v", file_version, ".txt.gz")

  aliases_url <- paste0(protocol, "://stringdb-downloads.org/download/protein.aliases.v", file_version, "/", aliases_base)
  info_url <- paste0(protocol, "://stringdb-downloads.org/download/protein.info.v", file_version, "/", info_base)

  network_type_param <- ""
  if (tolower(string_db$network_type) == "physical") {
    network_type_param <- "physical."
  }
  link_data_param <- "links.v"
  link_data <- .stringdb_safe_link_data(string_db)
  if (link_data == "detailed") {
    link_data_param <- "links.detailed.v"
  } else if (link_data == "full") {
    link_data_param <- "links.full.v"
  }
  links_base <- paste0(species, ".protein.", network_type_param, link_data_param, file_version, ".txt.gz")
  links_url <- paste0(protocol, "://stringdb-downloads.org/download/protein.", network_type_param, link_data_param, file_version, "/", links_base)

  .stringdb_download_with_md5(aliases_url, file.path(cache_dir, aliases_base), force_update = force_update)
  .stringdb_download_with_md5(info_url, file.path(cache_dir, info_base), force_update = force_update)
  .stringdb_download_with_md5(links_url, file.path(cache_dir, links_base), force_update = force_update)
  invisible(TRUE)
}

.stringdb_connect <- function(species = 9606, score_threshold = 0, version = "12.0", cache_dir = NULL) {
  cache_dir_use <- cache_dir
  if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
    cache_dir_use <- tempdir()
  }
  dir.create(cache_dir_use, recursive = TRUE, showWarnings = FALSE)

  args <- list(species = as.integer(species), score_threshold = score_threshold, version = version)
  if (!is.null(cache_dir_use) && nzchar(cache_dir_use)) args$input_directory <- cache_dir_use
  string_db <- tryCatch(
    do.call(STRINGdb$new, args),
    error = function(e) STRINGdb$new(species = as.integer(species), score_threshold = score_threshold, version = version)
  )
  .stringdb_prepare_cache_files(string_db, cache_dir_use)
  string_db
}

.external_reference_hash <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    return(NA_character_)
  }
  digest::digest(x, algo = "sha256")
}

.standalone_stringdb_cache_env <- new.env(parent = emptyenv())

.stringdb_cache_env <- function() {
  .standalone_stringdb_cache_env
}

.stringdb_cache_get <- function(scope_obj = NULL, cache_dir = NULL, cache_type, cache_key) {
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

  if (!is.null(cache_dir) && nzchar(cache_dir)) {
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    cache_path <- file.path(cache_dir, paste0("stringdb_", cache_type, "_", cache_key, ".rds"))
    if (file.exists(cache_path)) {
      return(list(value = readRDS(cache_path), cache_hit = TRUE, scope_obj = scope_obj))
    }
  }

  list(value = NULL, cache_hit = FALSE, scope_obj = scope_obj)
}

.stringdb_cache_put <- function(scope_obj = NULL, cache_dir = NULL, cache_type, cache_key, value) {
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

  if (!is.null(cache_dir) && nzchar(cache_dir)) {
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    cache_path <- file.path(cache_dir, paste0("stringdb_", cache_type, "_", cache_key, ".rds"))
    saveRDS(value, cache_path)
  }

  list(scope_obj = scope_obj)
}

.stringdb_map_genes_cached <- function(string_db, genes, input_id_type = "gene", scope_obj = NULL, cache_dir = NULL) {
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
  data.frame(
    from = ppi$protein1,
    to = ppi$protein2,
    combined_score = as.numeric(ppi$combined_score),
    stringsAsFactors = FALSE
  )
}

.stringdb_get_interactions_cached <- function(string_db, string_ids, score_threshold, scope_obj = NULL, cache_dir = NULL) {
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

stringdb_subscore_columns <- function() {
  c("nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore")
}

stringdb_obj_cache <- new.env(parent = emptyenv())

is_valid_string_version <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
}

get_stringdb_obj_for_details <- function(string_version, cache_dir) {
  key <- if (is_valid_string_version(string_version)) paste0("v", string_version) else "default"
  if (exists(key, envir = stringdb_obj_cache, inherits = FALSE)) {
    return(get(key, envir = stringdb_obj_cache))
  }
  args <- list(species = 9606, score_threshold = 0)
  if (is_valid_string_version(string_version)) args$version <- string_version
  if (!is.null(cache_dir) && nzchar(cache_dir)) args$input_directory <- cache_dir
  obj <- tryCatch(
    do.call(STRINGdb$new, args),
    error = function(e) STRINGdb$new(species = 9606, score_threshold = 0)
  )
  tryCatch({
    obj$link_data <- "detailed"
  }, error = function(e) NULL)
  assign(key, obj, envir = stringdb_obj_cache)
  obj
}

stringdb_detailed_links_path <- function(string_db, cache_dir) {
  network_type_param <- ""
  if (tolower(string_db$network_type) == "physical") {
    network_type_param <- "physical."
  }
  file_version <- string_db$file_version
  species <- string_db$species
  file_base <- paste0(species, ".protein.", network_type_param, "links.detailed.v", file_version, ".txt.gz")
  url <- paste0(
    string_db$protocol,
    "://stringdb-downloads.org/download/protein.",
    network_type_param,
    "links.detailed.v",
    file_version,
    "/",
    file_base
  )
  base_dir <- if (!is.null(cache_dir) && nzchar(cache_dir)) cache_dir else string_db$input_directory
  list(url = url, file_path = file.path(base_dir, file_base))
}

download_if_missing <- function(url, path) {
  if (file.exists(path) && file.info(path)$size > 0) return(invisible(TRUE))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  tmp <- paste0(path, ".tmp")
  if (file.exists(tmp)) unlink(tmp)

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  if (!is.finite(old_timeout) || old_timeout < 600) options(timeout = 600)

  ok <- FALSE
  for (attempt in seq_len(3L)) {
    warn_msgs <- character()
    ok <- tryCatch({
      withCallingHandlers(
        utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      TRUE
    }, error = function(e) {
      FALSE
    })
    bad_warn <- any(grepl("downloaded length|timeout", warn_msgs, ignore.case = TRUE))
    if (isTRUE(ok) && !bad_warn && file.exists(tmp) && file.info(tmp)$size > 0) break
    if (file.exists(tmp)) unlink(tmp)
    Sys.sleep(attempt)
  }

  if (isTRUE(ok) && file.exists(tmp) && file.info(tmp)$size > 0) {
    ok_rename <- file.rename(tmp, path)
    if (!isTRUE(ok_rename)) {
      file.copy(tmp, path, overwrite = TRUE)
      unlink(tmp)
    }
  }
  invisible(file.exists(path) && file.info(path)$size > 0)
}

read_stringdb_detailed_subset <- function(file_path, ids, cols) {
  ids <- unique(na.omit(as.character(ids)))
  if (!length(ids)) return(data.table())
  zcat_cmd <- Sys.which("gzcat")
  if (!nzchar(zcat_cmd)) zcat_cmd <- Sys.which("zcat")
  gzip_cmd <- Sys.which("gzip")
  if (!nzchar(zcat_cmd) && nzchar(gzip_cmd)) zcat_cmd <- gzip_cmd
  awk_cmd <- Sys.which("awk")

  dt <- data.table()
  if (nzchar(zcat_cmd) && nzchar(awk_cmd)) {
    ids_file <- tempfile(fileext = ".txt")
    on.exit(unlink(ids_file), add = TRUE)
    writeLines(ids, ids_file)
    decompress <- if (basename(zcat_cmd) == "gzip") {
      paste(zcat_cmd, "-dc", shQuote(file_path))
    } else {
      paste(zcat_cmd, shQuote(file_path))
    }
    cmd <- paste(
      decompress,
      "|",
      awk_cmd,
      shQuote("NR==FNR{a[$1]=1;next} FNR==1{print;next} ($1 in a) && ($2 in a){print}"),
      shQuote(ids_file),
      "-"
    )
    dt <- tryCatch(fread(cmd = cmd, showProgress = FALSE), error = function(e) data.table())
  } else {
    dt <- tryCatch(fread(file_path, showProgress = FALSE), error = function(e) data.table())
  }
  if (!nrow(dt)) return(data.table())
  if (!all(cols %in% names(dt))) return(data.table())
  dt[, ..cols]
}

append_stringdb_subscores_bulk <- function(detail_dt, string_version, cache_dir) {
  required <- c("string_from", "string_to")
  if (!all(required %in% names(detail_dt))) {
    return(detail_dt)
  }

  sub_cols <- stringdb_subscore_columns()
  existing <- intersect(sub_cols, names(detail_dt))
  if (length(existing)) {
    has_any <- any(vapply(existing, function(col) any(!is.na(detail_dt[[col]])), logical(1)))
    if (has_any) {
      return(detail_dt)
    }
    detail_dt[, (existing) := NULL]
  }
  ids <- unique(na.omit(c(detail_dt$string_from, detail_dt$string_to)))
  if (!length(ids)) {
    detail_dt[, (stringdb_subscore_columns()) := NA_real_]
    return(detail_dt)
  }

  string_db <- get_stringdb_obj_for_details(string_version, cache_dir)
  links_info <- stringdb_detailed_links_path(string_db, cache_dir)
  if (!file.exists(links_info$file_path)) {
    download_if_missing(links_info$url, links_info$file_path)
  }
  if (!file.exists(links_info$file_path)) {
    detail_dt[, (stringdb_subscore_columns()) := NA_real_]
    return(detail_dt)
  }

  cols <- c(
    "protein1",
    "protein2",
    "neighborhood",
    "fusion",
    "cooccurence",
    "coexpression",
    "experimental",
    "database",
    "textmining",
    "combined_score"
  )
  interactions <- read_stringdb_detailed_subset(links_info$file_path, ids, cols)
  if (!nrow(interactions)) {
    detail_dt[, (stringdb_subscore_columns()) := NA_real_]
    return(detail_dt)
  }

  interactions <- interactions[protein1 %chin% ids & protein2 %chin% ids]
  if (!nrow(interactions)) {
    detail_dt[, (stringdb_subscore_columns()) := NA_real_]
    return(detail_dt)
  }

  interactions[, key := paste(pmin(protein1, protein2), pmax(protein1, protein2), sep = "|")]
  interactions <- unique(interactions, by = "key")
  interactions[, `:=`(
    nscore = as.numeric(neighborhood),
    fscore = as.numeric(fusion),
    pscore = as.numeric(cooccurence),
    ascore = as.numeric(coexpression),
    escore = as.numeric(experimental),
    dscore = as.numeric(database),
    tscore = as.numeric(textmining)
  )]
  sub <- interactions[, .(key, nscore, fscore, pscore, ascore, escore, dscore, tscore)]

  detail_dt[, .row := .I]
  detail_dt[, key := paste(pmin(string_from, string_to), pmax(string_from, string_to), sep = "|")]
  detail_dt <- merge(detail_dt, sub, by = "key", all.x = TRUE, sort = FALSE)
  setorder(detail_dt, .row)
  detail_dt[, c("key", ".row") := NULL]
  detail_dt
}

option_list <- list(
  make_option(c("--bench_root"), type = "character", help = "Benchmark output root (contains <method>/repeat_*/edges_all.tsv)"),
  make_option(c("--outdir"), type = "character", help = "Output directory for mapped edges"),
  make_option(c("--methods"), type = "character", default = "genescope,giotto,hotspot"),
  make_option(c("--species"), type = "integer", default = 9606),
  make_option(c("--input_id_type"), type = "character", default = "gene"),
  make_option(c("--string_version"), type = "character", default = "12.0",
              help = "STRING dataset version (passed to STRINGdb$new; default 12.0)."),
  make_option(c("--string_score_thresholds"), type = "character", default = "400,700",
              help = "Comma-separated STRING score thresholds to precompute background positive counts (e.g., 400,700)."),
  make_option(c("--keep_subscores"), type = "integer", default = 1,
              help = "If 1, append STRING subscores (nscore..tscore) into mapped output (requires links.detailed).")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$bench_root) || is.null(opt$outdir)) {
  stop("Missing --bench_root or --outdir")
}

methods <- split_csv(opt$methods)
if (!length(methods)) stop("No methods provided")

thresholds <- parse_int_list(opt$string_score_thresholds)
thresholds <- thresholds[is.finite(thresholds) & thresholds >= 0L & thresholds <= 1000L]
if (!length(thresholds)) thresholds <- c(700L)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
string_version <- normalize_string_version(opt$string_version)
if (is.na(string_version)) stop("Invalid --string_version: ", opt$string_version)
cache_dir <- file.path(opt$outdir, paste0("STRINGdb_cache_v", string_version))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

score_source <- "standalone"
connect_fn <- .stringdb_connect
map_fn <- .stringdb_map_genes_cached
interactions_fn <- .stringdb_get_interactions_cached
label_fn <- .stringdb_label_edges

mapped_rows <- list()
per_method <- list()

for (method in methods) {
  rep_dir <- get_first_repeat_dir(opt$bench_root, method)
  edges_all_path <- file.path(rep_dir, "edges_all.tsv")
  meta_path <- file.path(rep_dir, "meta.json")
  if (!file.exists(edges_all_path)) {
    stop("Missing edges file for method ", method, ": ", edges_all_path)
  }
  edges_all <- read_edges_all(edges_all_path)
  genes <- unique(c(edges_all$gene_a, edges_all$gene_b))

  string_db <- connect_fn(
    species = as.integer(opt$species),
    score_threshold = 0,
    version = string_version,
    cache_dir = cache_dir
  )
  string_version_actual <- tryCatch(as.character(string_db$version), error = function(e) NA_character_)
  mapping_res <- map_fn(string_db, genes, input_id_type = opt$input_id_type, scope_obj = NULL, cache_dir = cache_dir)
  mapping <- mapping_res$mapping
  interactions_res <- interactions_fn(string_db, mapping, score_threshold = 0, scope_obj = NULL, cache_dir = cache_dir)
  interactions <- interactions_res$interactions

  edge_df <- edges_all[, .(from = gene_a, to = gene_b, weight = weight, fdr = fdr)]
  edge_df <- label_fn(edge_df, gene_to_string = mapping, interactions = interactions, score_threshold = max(thresholds))
  edge_df[, method := method]
  if (opt$keep_subscores != 0) {
    edge_df <- append_stringdb_subscores_bulk(edge_df, string_version_actual, cache_dir)
  }

  mapped_rows[[method]] <- edge_df

  bench_meta <- NULL
  if (file.exists(meta_path)) {
    bench_meta <- tryCatch(read_json(meta_path, simplifyVector = TRUE), error = function(e) NULL)
  }

  mapped_ids <- unique(na.omit(as.character(mapping)))
  bg_gene_count <- length(mapped_ids)
  n_pairs_bg <- bg_gene_count * (bg_gene_count - 1) / 2
  n_pos_bg_by_threshold <- list()
  if (is.data.frame(interactions) && nrow(interactions)) {
    keys_all <- paste(pmin(interactions$from, interactions$to), pmax(interactions$from, interactions$to), sep = "|")
    for (cut in thresholds) {
      pos_idx <- interactions$combined_score >= cut
      n_pos_bg_by_threshold[[as.character(cut)]] <- as.integer(length(unique(keys_all[pos_idx])))
    }
  } else {
    for (cut in thresholds) {
      n_pos_bg_by_threshold[[as.character(cut)]] <- 0L
    }
  }

  per_method[[method]] <- list(
    method = method,
    rep_dir = rep_dir,
    edges_all_path = edges_all_path,
    edges_all_sha256 = sha256_file(edges_all_path),
    meta_path = meta_path,
    meta_sha256 = sha256_file(meta_path),
    bench_meta = bench_meta,
    string_version = string_version_actual,
    string_version_requested = string_version,
    string_version_actual = string_version_actual,
    mapping = list(
      n_input_genes = length(genes),
      n_mapped = sum(!is.na(mapping)),
      coverage = if (length(genes)) sum(!is.na(mapping)) / length(genes) else NA_real_
    ),
    bg = list(
      bg_gene_count = as.integer(bg_gene_count),
      n_pairs_bg = as.numeric(n_pairs_bg),
      n_pos_bg_by_threshold = n_pos_bg_by_threshold
    ),
    mapped_edges = list(
      n_edges = nrow(edge_df),
      n_comparable_edges = sum(!is.na(edge_df$string_score))
    )
  )
}

mapped_dt <- rbindlist(mapped_rows, use.names = TRUE, fill = TRUE)
out_tsv <- file.path(opt$outdir, "stringdb_edge_mapped_all.tsv")
fwrite(mapped_dt, out_tsv, sep = "\t")

meta <- list(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  bench_root = opt$bench_root,
  outdir = opt$outdir,
  methods = methods,
  species = as.integer(opt$species),
  input_id_type = opt$input_id_type,
  string_version = string_version,
  string_score_thresholds = thresholds,
  keep_subscores = opt$keep_subscores != 0,
  score_function_source = score_source,
  outputs = list(
    mapped_edges_tsv = out_tsv,
    cache_dir = cache_dir
  ),
  per_method = per_method
)
write_json(meta, file.path(opt$outdir, "stringdb_edge_mapped_meta.json"), auto_unbox = TRUE, pretty = TRUE)

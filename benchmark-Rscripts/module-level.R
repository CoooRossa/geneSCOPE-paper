#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(Rcpp)
  library(RcppParallel)
})

options(stringsAsFactors = FALSE)

mean_score_label <- "Mean score (missing=0)"

make_percent_progress_logger <- function(total, label = "Progress", step_percent = 5L) {
  total <- suppressWarnings(as.integer(total))
  step_percent <- suppressWarnings(as.integer(step_percent))
  if (!is.finite(step_percent) || step_percent < 1L) step_percent <- 5L
  if (!is.finite(total) || total <= 0L) {
    return(function(i, detail = NULL) invisible(FALSE))
  }
  last_bucket <- -1L
  function(i, detail = NULL) {
    i <- suppressWarnings(as.integer(i))
    if (!is.finite(i) || i < 0L) return(invisible(FALSE))
    if (i > total) i <- total
    pct <- 100 * i / total
    bucket <- as.integer(floor(pct / step_percent))
    if (bucket <= last_bucket && i < total) return(invisible(FALSE))
    last_bucket <<- bucket
    msg <- sprintf("%s: %d/%d (%.1f%%)", label, i, total, pct)
    if (!is.null(detail) && nzchar(as.character(detail))) {
      msg <- paste0(msg, " - ", as.character(detail))
    }
    message(msg)
    invisible(TRUE)
  }
}

split_csv <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || !nzchar(x)) return(character(0))
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  trimws(parts[nzchar(parts)])
}

display_method_name <- function(x) {
  x_chr <- as.character(x)
  key <- tolower(trimws(x_chr))
  if (key == "genescope" || startsWith(key, "genescope")) return("geneSCOPE")
  if (key == "giotto" || key == "giotto_r" || startsWith(key, "giotto")) return("Giotto Suite")
  if (key == "hotspot" || key == "hotspot_py" || startsWith(key, "hotspot")) return("Hotspot")
  if (key == "seagal" || startsWith(key, "seagal")) return("Seagal")
  if (key == "smoothie") return("Smoothie")
  if (key == "spatialcorr") return("SpatialCorr")
  x_chr
}

bench_base_color_map <- function(method_labels) {
  method_labels <- unique(as.character(method_labels))
  if (!length(method_labels)) return(character(0))
  base <- c(
    "geneSCOPE" = "#F8766D",
    "Giotto Suite" = "#00BA38",
    "Hotspot" = "#619CFF",
    "Seagal" = "#9B59B6"
  )
  out <- base[method_labels]
  missing <- is.na(out)
  if (any(missing)) {
    out[missing] <- "#7F7F7F"
  }
  stats::setNames(as.character(out), method_labels)
}

boost_color_vividness <- function(color_hex, vividness = 0.25) {
  vividness <- suppressWarnings(as.numeric(vividness))
  if (!is.finite(vividness)) vividness <- 0.25
  vividness <- max(0, vividness)
  if (is.na(color_hex) || !nzchar(color_hex)) return("#000000")
  rgb <- tryCatch(grDevices::col2rgb(color_hex), error = function(e) NULL)
  if (is.null(rgb)) return(color_hex)
  hsv <- grDevices::rgb2hsv(rgb[1, ] / 255, rgb[2, ] / 255, rgb[3, ] / 255)
  h <- hsv["h", ]
  s <- hsv["s", ]
  v <- hsv["v", ]
  s2 <- pmin(1, s + (1 - s) * vividness)
  grDevices::hsv(h = h, s = s2, v = v)
}

theme_bench_like_stringdb <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

merge_title_subtitle <- function(title, subtitle) {
  title <- as.character(title)
  subtitle <- as.character(subtitle)
  if (is.na(subtitle) || !nzchar(subtitle)) return(title)
  paste0(title, " — ", subtitle)
}

percent_label_0_100 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- ifelse(!is.finite(x), NA_character_, paste0(formatC(x, format = "f", digits = 1), "%"))
  sub("\\.0%$", "%", out)
}

format_p_value <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  if (!is.finite(p)) return("p=NA")
  if (p < 1e-4) return("p<1e-4")
  paste0("p=", signif(p, 3))
}

as_int1 <- function(x, default = NA_integer_) {
  x <- suppressWarnings(as.integer(x))
  if (!is.finite(x) || is.na(x)) return(default)
  x
}

safe_choose2 <- function(n) {
  n <- suppressWarnings(as.numeric(n))
  if (!is.finite(n) || n < 2) return(0)
  n * (n - 1) / 2
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

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

load_cpp_bindings <- function() {
  script_dir <- get_script_dir()
  cpp_path <- file.path(script_dir, "module-level.cpp")
  if (!file.exists(cpp_path)) {
    stop("Missing C++ source file: ", cpp_path)
  }
  Rcpp::sourceCpp(cpp_path, verbose = FALSE)
}

map_gene_indices <- function(genes, gene_index_map) {
  idx <- suppressWarnings(as.integer(gene_index_map[genes]))
  idx <- idx[is.finite(idx)]
  unique(idx)
}

read_modules <- function(path) {
  dt <- fread(path)
  required <- c("gene", "module_id")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop("Missing columns in modules file ", path, ": ", paste(missing, collapse = ", "))
  }
  dt <- dt[, .(
    gene = trimws(as.character(gene)),
    module_id = trimws(as.character(module_id))
  )]
  dt <- dt[nzchar(gene)]
  dt <- dt[!is.na(module_id) & nzchar(module_id)]
  dt <- dt[module_id != "-1"]
  dt <- unique(dt, by = c("module_id", "gene"))
  dt
}

read_background_genes <- function(path) {
  dt <- fread(path)
  if (!"gene" %in% names(dt)) {
    stop("Missing column in modules file ", path, ": gene")
  }
  dt <- dt[, .(gene = trimws(as.character(gene)))]
  dt <- dt[nzchar(gene)]
  unique(dt$gene)
}

get_git_commit <- function() {
  cmd <- Sys.which("git")
  if (!nzchar(cmd)) return("")
  is_git <- suppressWarnings(tryCatch(
    system2(cmd, c("rev-parse", "--is-inside-work-tree"), stdout = TRUE, stderr = TRUE),
    error = function(e) ""
  ))
  if (!length(is_git) || !identical(trimws(is_git[1]), "true")) return("")
  out <- suppressWarnings(tryCatch(
    system2(cmd, c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE),
    error = function(e) ""
  ))
  if (!length(out)) return("")
  commit <- trimws(out[1])
  if (!nzchar(commit) || grepl("fatal", commit, ignore.case = TRUE)) return("")
  commit
}

save_empty_plot <- function(path, title, width, height, dpi) {
  p <- ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available", hjust = 0.5, vjust = 0.5)
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
}

sanitize_filename <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x) | !nzchar(x), "NA", x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "NA")
}

extract_first_number <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  direct <- suppressWarnings(as.numeric(x))
  if (is.finite(direct)) return(direct)
  m <- regmatches(x, regexpr("-?[0-9]+", x, perl = TRUE))
  if (length(m) && nzchar(m)) {
    out <- suppressWarnings(as.numeric(m))
    if (is.finite(out)) return(out)
  }
  NA_real_
}

write_tsv_gz <- function(dt, path) {
  dt <- as.data.table(dt)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    tmp <- tempfile(fileext = ".tsv")
    fwrite(dt, tmp, sep = "\t")
    gzip_cmd <- Sys.which("gzip")
    if (nzchar(gzip_cmd)) {
      err <- tempfile(fileext = ".log")
      status <- suppressWarnings(system2(gzip_cmd, c("-c", tmp), stdout = path, stderr = err))
      unlink(tmp)
      if (!identical(status, 0L)) {
        err_msg <- ""
        if (file.exists(err)) {
          err_msg <- paste(readLines(err, warn = FALSE), collapse = "\n")
        }
        unlink(err)
        if (nzchar(err_msg)) {
          stop("gzip failed writing: ", path, "\n", err_msg)
        }
        stop("gzip failed writing: ", path)
      }
      unlink(err)
      return(invisible(path))
    }
    con <- gzfile(path, open = "wb")
    on.exit(close(con), add = TRUE)
    writeLines(paste(names(dt), collapse = "\t"), con)
    if (nrow(dt)) {
      write.table(dt, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    unlink(tmp)
    return(invisible(path))
  }
  fwrite(dt, path, sep = "\t")
  invisible(path)
}

plot_null_hist <- function(values, observed, title, xlab, subtitle, out_path,
                           width = 6, height = 4, dpi = 150, value_is_percent = FALSE) {
  values <- suppressWarnings(as.numeric(values))
  observed <- suppressWarnings(as.numeric(observed))
  if (!length(values) || !any(is.finite(values)) || !is.finite(observed)) {
    save_empty_plot(out_path, title, width, height, dpi)
    return(invisible(FALSE))
  }
  plot_dt <- data.table(value = values[is.finite(values)])
  p <- ggplot(plot_dt, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "grey80", color = "grey40", linewidth = 0.3) +
    geom_density(color = "grey20", linewidth = 0.6, na.rm = TRUE) +
    geom_vline(xintercept = observed, color = "#B2182B", linewidth = 0.7) +
    theme_bench_like_stringdb() +
    labs(title = merge_title_subtitle(title, subtitle), x = xlab, y = "Density")
  if (isTRUE(value_is_percent)) {
    p <- p + scale_x_continuous(labels = percent_label_0_100)
  }
  ggsave(filename = out_path, plot = p, width = width, height = height, dpi = dpi)
  invisible(TRUE)
}

plot_null_boxplot <- function(values, observed, title, ylab, subtitle, out_path,
                              width = 6, height = 4, dpi = 150, max_points = 3000L, seed = 1L, value_is_percent = FALSE) {
  values <- suppressWarnings(as.numeric(values))
  observed <- suppressWarnings(as.numeric(observed))
  if (!length(values) || !any(is.finite(values)) || !is.finite(observed)) {
    save_empty_plot(out_path, title, width, height, dpi)
    return(invisible(FALSE))
  }
  vals <- values[is.finite(values)]
  plot_dt <- data.table(group = "null", value = vals)

  point_dt <- plot_dt
  if (is.finite(max_points) && max_points > 0L && nrow(point_dt) > max_points) {
    sel <- sample_int_deterministic(nrow(point_dt), max_points, seed = seed)
    point_dt <- point_dt[sel]
  }

  p <- ggplot(plot_dt, aes(x = group, y = value)) +
    geom_boxplot(fill = "grey80", color = "grey40", linewidth = 0.3, outlier.size = 0.3, na.rm = TRUE) +
    geom_point(
      data = point_dt,
      mapping = aes(x = group, y = value),
      position = position_jitter(width = 0.18, height = 0),
      size = 0.2,
      alpha = 0.2,
      color = "grey20",
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    geom_hline(yintercept = observed, color = "#B2182B", linewidth = 0.7) +
    theme_bench_like_stringdb() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(title = merge_title_subtitle(title, subtitle), x = NULL, y = ylab)
  if (isTRUE(value_is_percent)) {
    p <- p + scale_y_continuous(labels = percent_label_0_100)
  }

  ggsave(filename = out_path, plot = p, width = width, height = height, dpi = dpi)
  invisible(TRUE)
}

sample_int_deterministic <- function(n, k, seed) {
  n <- suppressWarnings(as.integer(n))
  k <- suppressWarnings(as.integer(k))
  seed <- suppressWarnings(as.integer(seed))
  if (!is.finite(n) || n <= 0L || !is.finite(k) || k <= 0L) return(integer(0))
  if (k >= n) return(seq_len(n))
  if (!is.finite(seed)) seed <- 1L
  old_seed <- get0(".Random.seed", ifnotfound = NULL, inherits = TRUE)
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  sample.int(n, size = k, replace = FALSE)
}

option_list <- list(
  make_option(c("--map_dir"), type = "character", help = "Directory containing stringdb_edge_mapped_all.tsv"),
  make_option(c("--outdir"), type = "character", help = "Output directory"),
  make_option(c("--methods"), type = "character", default = "genescope,giotto,hotspot"),
  make_option(c("--modules_tsv_by_method"), type = "character", help = "Comma-separated modules.tsv paths"),
  make_option(c("--string_score_threshold"), type = "integer", default = 700,
              help = "STRING score threshold used for the positive-edge fraction metric (score >= threshold)"),
  make_option(c("--n_random"), type = "integer", default = 100000),
  make_option(c("--seed"), type = "integer", default = 1),
  make_option(c("--min_module_genes"), type = "integer", default = 2),
  make_option(c("--min_valid_null_draws"), type = "integer", default = NA_integer_),
  make_option(c("--max_resample_attempts"), type = "integer", default = 20),
  make_option(c("--plot_max_null_points_per_method"), type = "integer", default = 20000,
              help = "Max number of null draw dots shown per method in plots (default: 20000; values are sampled)"),
  make_option(c("--plot_max_modules_per_method"), type = "integer", default = 200)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$map_dir) || !nzchar(opt$map_dir)) stop("Missing --map_dir")
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("Missing --outdir")
if (is.null(opt$modules_tsv_by_method) || !nzchar(opt$modules_tsv_by_method)) {
  stop("Missing --modules_tsv_by_method")
}

methods <- split_csv(opt$methods)
module_paths <- split_csv(opt$modules_tsv_by_method)
if (!length(methods)) stop("No methods provided via --methods")
if (length(module_paths) != length(methods)) {
  stop("Length mismatch: --methods has ", length(methods),
       " entries but --modules_tsv_by_method has ", length(module_paths))
}

string_score_threshold <- as_int1(opt$string_score_threshold, default = 700L)
n_random <- as_int1(opt$n_random, default = 100000L)
seed <- as_int1(opt$seed, default = 1L)
min_module_genes <- as_int1(opt$min_module_genes, default = 2L)
min_valid_null_draws <- as_int1(opt$min_valid_null_draws, default = NA_integer_)
max_resample_attempts <- as_int1(opt$max_resample_attempts, default = 20L)
plot_max_null_points_per_method <- as_int1(opt$plot_max_null_points_per_method, default = 20000L)
plot_max_modules_per_method <- as_int1(opt$plot_max_modules_per_method, default = 200L)

if (!is.finite(string_score_threshold)) stop("Invalid --string_score_threshold")
if (!is.finite(n_random) || n_random < 0L) stop("Invalid --n_random")
if (!is.finite(seed)) stop("Invalid --seed")
if (!is.finite(min_module_genes) || min_module_genes < 1L) stop("Invalid --min_module_genes")
if (!is.finite(max_resample_attempts) || max_resample_attempts < 0L) {
  stop("Invalid --max_resample_attempts")
}
if (!is.finite(plot_max_null_points_per_method) || plot_max_null_points_per_method < 0L) {
  stop("Invalid --plot_max_null_points_per_method")
}
if (!is.finite(plot_max_modules_per_method) || plot_max_modules_per_method < 1L) {
  stop("Invalid --plot_max_modules_per_method")
}
if (!is.finite(min_valid_null_draws)) {
  min_valid_null_draws <- floor(0.8 * n_random)
}
if (!is.finite(min_valid_null_draws) || min_valid_null_draws < 0L) {
  stop("Invalid --min_valid_null_draws")
}
if (is.finite(min_valid_null_draws) && is.finite(n_random) && min_valid_null_draws > n_random) {
  stop("Invalid --min_valid_null_draws (", min_valid_null_draws, ") > --n_random (", n_random, ")")
}

t_start <- Sys.time()
message(sprintf(
  "[INFO] Start anchor-cross enrichment: methods=%s, n_random=%d, string_score_threshold=%d",
  paste(methods, collapse = ","),
  n_random,
  string_score_threshold
))
message("[INFO] Loading C++ bindings...")
load_cpp_bindings()

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

draw_dir <- file.path(opt$outdir, "null_draws_anchor_cross")
plot_dir <- file.path(opt$outdir, "null_plots_anchor_cross")
dir.create(draw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

map_path <- file.path(opt$map_dir, "stringdb_edge_mapped_all.tsv")
if (!file.exists(map_path)) stop("Missing mapping file: ", map_path)

map_meta_path <- file.path(opt$map_dir, "stringdb_edge_mapped_meta.json")
map_meta <- NULL
if (file.exists(map_meta_path)) {
  map_meta <- tryCatch(jsonlite::read_json(map_meta_path, simplifyVector = TRUE), error = function(e) NULL)
}

set.seed(seed)

modules_by_method <- list()
module_hashes <- list()
background_gene_lists <- list()
method_progress <- make_percent_progress_logger(length(methods), label = "Load modules", step_percent = 25L)
for (i in seq_along(methods)) {
  method <- methods[[i]]
  method_progress(i, detail = paste0("method=", method))
  path <- module_paths[[i]]
  if (!file.exists(path)) stop("Missing modules file for method ", method, ": ", path)
  module_hashes[[method]] <- sha256_file(path)
  background_gene_lists[[length(background_gene_lists) + 1L]] <- read_background_genes(path)
  mod_dt <- read_modules(path)
  if (nrow(mod_dt)) {
    mod_counts <- mod_dt[, .N, by = .(module_id)]
    keep_ids <- mod_counts[N >= min_module_genes, module_id]
    mod_dt <- mod_dt[module_id %in% keep_ids]
  }
  modules_by_method[[method]] <- mod_dt
}

background_genes <- sort(unique(unlist(background_gene_lists, use.names = FALSE)))
background_size <- length(background_genes)
gene_index_map <- stats::setNames(seq_along(background_genes) - 1L, background_genes)

map_dt <- fread(map_path, select = c("from", "to", "string_score"))
map_dt <- map_dt[!is.na(from) & !is.na(to)]
map_dt[, from := trimws(as.character(from))]
map_dt[, to := trimws(as.character(to))]
map_dt <- map_dt[nzchar(from) & nzchar(to)]
map_dt[, string_score := suppressWarnings(as.numeric(string_score))]
map_dt <- map_dt[is.finite(string_score)]
map_dt <- map_dt[from != to]
map_dt[, gene1 := ifelse(from <= to, from, to)]
map_dt[, gene2 := ifelse(from <= to, to, from)]
edges <- map_dt[, .(string_score = max(string_score, na.rm = TRUE)), by = .(gene1, gene2)]
setkey(edges, gene1, gene2)

edge_i <- suppressWarnings(as.integer(gene_index_map[edges$gene1]))
edge_j <- suppressWarnings(as.integer(gene_index_map[edges$gene2]))
edge_score <- suppressWarnings(as.numeric(edges$string_score))
keep_edges <- is.finite(edge_i) & is.finite(edge_j)
edge_map_ptr <- build_edge_map_cpp(edge_i[keep_edges], edge_j[keep_edges], edge_score[keep_edges])

module_rows <- list()
module_counter <- 0L
null_draw_samples <- list()
null_box_stats_eavg <- list()
null_box_stats_pos <- list()

module_counts <- vapply(methods, function(method) {
  mod_dt <- modules_by_method[[method]]
  if (is.null(mod_dt) || !nrow(mod_dt)) return(0L)
  as.integer(length(unique(mod_dt$module_id)))
}, integer(1))
names(module_counts) <- methods
total_modules <- as.integer(sum(module_counts))
step_pct <- if (total_modules <= 200L) 1L else 5L
module_progress <- make_percent_progress_logger(total_modules, label = "Score modules", step_percent = step_pct)
message(sprintf("[INFO] Scoring modules: total=%d", total_modules))

plot_width <- 6
plot_height <- 4
plot_dpi <- 150

for (method in methods) {
  mod_dt <- modules_by_method[[method]]
  if (is.null(mod_dt) || !nrow(mod_dt)) next
  module_ids <- sort(unique(mod_dt$module_id))

  for (mid in module_ids) {
    module_counter <- module_counter + 1L
    module_progress(module_counter, detail = paste0("method=", method, ", module_id=", mid))
    genes <- unique(mod_dt[module_id == mid, gene])
    module_size <- length(genes)
    n_pairs_total <- safe_choose2(module_size)
    n_edges_target <- as_int1(n_pairs_total, default = NA_integer_)

    module_idx <- map_gene_indices(genes, gene_index_map)
    edge_stats <- module_pair_stats_threshold_cpp(module_idx, edge_map_ptr, threshold = string_score_threshold)
    n_pairs_mapped_obs <- as.integer(edge_stats$n_pairs_mapped)
    e_sum_obs <- as.numeric(edge_stats$sum_mapped)
    n_pairs_pos_obs <- as.integer(edge_stats$n_pairs_ge_threshold)
    e_avg_obs <- if (n_pairs_total > 0) e_sum_obs / n_pairs_total else NA_real_
    mapped_pair_frac_obs <- if (n_pairs_total > 0) n_pairs_mapped_obs / n_pairs_total else NA_real_
    pos_pair_frac_obs <- if (n_pairs_total > 0) n_pairs_pos_obs / n_pairs_total else NA_real_

    status <- "ok"
    status_reason <- NA_character_
    null_mean_eavg <- NA_real_
    null_sd_eavg <- NA_real_
    delta_eavg <- NA_real_
    empirical_p_eavg <- NA_real_
    null_mean_pos_pair_frac <- NA_real_
    null_sd_pos_pair_frac <- NA_real_
    delta_pos_pair_frac <- NA_real_
    empirical_p_pos_pair_frac <- NA_real_
    n_valid_null_draws <- 0L
    cross_edge_space <- NA_real_

    null_eavg <- numeric(0)
    null_pos_pair_frac <- numeric(0)

    if (n_pairs_total <= 0 || is.na(n_edges_target) || !is.finite(n_edges_target) || n_edges_target <= 0L) {
      status <- "skipped"
      status_reason <- "n_pairs_total_zero"
    } else if (n_random <= 0L) {
      status <- "null_insufficient"
      status_reason <- "n_random_nonpositive"
    } else {
      bg_genes <- setdiff(background_genes, genes)
      bg_idx <- map_gene_indices(bg_genes, gene_index_map)
      m <- length(module_idx)
      b <- length(bg_idx)
      cross_edge_space <- as.numeric(m) * as.numeric(b)
      if (m <= 0L || b <= 0L || !is.finite(cross_edge_space) || cross_edge_space < n_pairs_total) {
        status <- "null_insufficient"
        status_reason <- "cross_edge_space_too_small"
      } else {
        module_seed <- seed + as.integer(module_counter) * 104729L
        null_payload <- draw_null_anchor_cross_eavg_posfrac_cpp(
          n_random = n_random,
          module_indices = module_idx,
          bg_indices = bg_idx,
          n_edges_target = n_edges_target,
          edge_map_ptr = edge_map_ptr,
          threshold = string_score_threshold,
          seed = module_seed
        )
        null_eavg <- null_payload$null_eavg
        null_pos_pair_frac <- null_payload$null_pos_pair_frac

        if (!length(null_eavg) || !length(null_pos_pair_frac)) {
          status <- "null_insufficient"
          status_reason <- "cross_edge_space_too_small"
        } else {
          valid_mask <- is.finite(null_eavg) & is.finite(null_pos_pair_frac)
          valid_null_eavg <- null_eavg[valid_mask]
          valid_null_pos_pair_frac <- null_pos_pair_frac[valid_mask]
          n_valid_null_draws <- length(valid_null_eavg)
          if (n_valid_null_draws < min_valid_null_draws) {
            status <- "null_insufficient"
            status_reason <- "valid_null_draws_lt_min"
          } else if (!is.finite(e_avg_obs)) {
            status <- "null_insufficient"
            status_reason <- "E_avg_obs_nonfinite"
          } else if (!is.finite(pos_pair_frac_obs)) {
            status <- "null_insufficient"
            status_reason <- "pos_pair_frac_obs_nonfinite"
          } else {
            null_mean_eavg <- mean(valid_null_eavg)
            null_sd_eavg <- stats::sd(valid_null_eavg)
            delta_eavg <- e_avg_obs - null_mean_eavg
            empirical_p_eavg <- (1 + sum(valid_null_eavg >= e_avg_obs)) / (1 + n_valid_null_draws)
            null_mean_pos_pair_frac <- mean(valid_null_pos_pair_frac)
            null_sd_pos_pair_frac <- stats::sd(valid_null_pos_pair_frac)
            delta_pos_pair_frac <- pos_pair_frac_obs - null_mean_pos_pair_frac
            empirical_p_pos_pair_frac <- (1 + sum(valid_null_pos_pair_frac >= pos_pair_frac_obs)) / (1 + n_valid_null_draws)
            status <- "ok"

            stats_eavg <- grDevices::boxplot.stats(valid_null_eavg, coef = 1.5)$stats
            stats_pos <- grDevices::boxplot.stats(valid_null_pos_pair_frac, coef = 1.5)$stats
            null_box_stats_eavg[[length(null_box_stats_eavg) + 1L]] <- data.table(
              method = as.character(method),
              method_label = display_method_name(method),
              module_id = as.character(mid),
              module_id_num = extract_first_number(mid),
              ymin = as.numeric(stats_eavg[1]),
              lower = as.numeric(stats_eavg[2]),
              middle = as.numeric(stats_eavg[3]),
              upper = as.numeric(stats_eavg[4]),
              ymax = as.numeric(stats_eavg[5]),
              observed = as.numeric(e_avg_obs),
              n_valid = as.integer(n_valid_null_draws)
            )
            null_box_stats_pos[[length(null_box_stats_pos) + 1L]] <- data.table(
              method = as.character(method),
              method_label = display_method_name(method),
              module_id = as.character(mid),
              module_id_num = extract_first_number(mid),
              ymin = as.numeric(stats_pos[1]),
              lower = as.numeric(stats_pos[2]),
              middle = as.numeric(stats_pos[3]),
              upper = as.numeric(stats_pos[4]),
              ymax = as.numeric(stats_pos[5]),
              observed = as.numeric(pos_pair_frac_obs),
              n_valid = as.integer(n_valid_null_draws)
            )

            if (plot_max_null_points_per_method > 0L) {
              denom_modules <- module_counts[[method]]
              if (!is.finite(denom_modules) || denom_modules < 1L) denom_modules <- 1L
              points_per_module <- as.integer(floor(plot_max_null_points_per_method / denom_modules))
              if (!is.finite(points_per_module) || points_per_module < 1L) points_per_module <- 1L
              sample_n <- min(points_per_module, n_valid_null_draws)
              if (sample_n > 0L) {
                sel <- sample_int_deterministic(n_valid_null_draws, sample_n, seed = module_seed + 11L)
                null_draw_samples[[length(null_draw_samples) + 1L]] <- data.table(
                  method = as.character(method),
                  method_label = display_method_name(method),
                  module_id = as.character(mid),
                  obs_eavg = as.numeric(e_avg_obs),
                  obs_pos_pair_frac = as.numeric(pos_pair_frac_obs),
                  null_eavg = as.numeric(valid_null_eavg[sel]),
                  null_pos_pair_frac = as.numeric(valid_null_pos_pair_frac[sel])
                )
              }
            }
          }
        }
      }
    }

    file_id <- sanitize_filename(paste0(method, "__", mid))

    draws_path <- file.path(draw_dir, paste0("module_", file_id, ".tsv.gz"))
    if (length(null_eavg) && length(null_pos_pair_frac)) {
      valid_mask <- is.finite(null_eavg) & is.finite(null_pos_pair_frac)
      draws_dt <- data.table(
        draw_idx = which(valid_mask),
        null_eavg = as.numeric(null_eavg[valid_mask]),
        null_pos_pair_frac = as.numeric(null_pos_pair_frac[valid_mask]),
        module_id = as.character(mid)
      )
    } else {
      draws_dt <- data.table(
        draw_idx = integer(),
        null_eavg = numeric(),
        null_pos_pair_frac = numeric(),
        module_id = character()
      )
    }
    write_tsv_gz(draws_dt, draws_path)

    title_base <- paste0(display_method_name(method), " ", mid)
    if (identical(status, "ok")) {
      subtitle_eavg <- sprintf(
        "k=%d, n_random=%d, n_valid=%d, %s",
        n_edges_target,
        n_random,
        n_valid_null_draws,
        format_p_value(empirical_p_eavg)
      )
      plot_null_hist(
        values = null_eavg,
        observed = e_avg_obs,
        title = paste0(title_base, ": Null score"),
        xlab = paste0("Null ", mean_score_label, " (cross edges)"),
        subtitle = subtitle_eavg,
        out_path = file.path(plot_dir, paste0("module_", file_id, "__Eavg.png")),
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi
      )
      plot_null_boxplot(
        values = null_eavg,
        observed = e_avg_obs,
        title = paste0(title_base, ": Null score (box)"),
        ylab = paste0("Null ", mean_score_label, " (cross edges)"),
        subtitle = subtitle_eavg,
        out_path = file.path(plot_dir, paste0("module_", file_id, "__Eavg_boxplot.png")),
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi,
        max_points = 3000L,
        seed = module_seed + 101L
      )

      subtitle_pos <- sprintf(
        "k=%d, n_random=%d, n_valid=%d, %s",
        n_edges_target,
        n_random,
        n_valid_null_draws,
        format_p_value(empirical_p_pos_pair_frac)
      )
      plot_null_hist(
        values = 100 * null_pos_pair_frac,
        observed = 100 * pos_pair_frac_obs,
        title = paste0(title_base, ": Null positive edge %"),
        xlab = paste0("Null positive edge % (score \u2265 ", string_score_threshold, ")"),
        subtitle = subtitle_pos,
        out_path = file.path(plot_dir, paste0("module_", file_id, "__PosFrac.png")),
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi,
        value_is_percent = TRUE
      )
      plot_null_boxplot(
        values = 100 * null_pos_pair_frac,
        observed = 100 * pos_pair_frac_obs,
        title = paste0(title_base, ": Null positive edge % (box)"),
        ylab = paste0("Null positive edge % (score \u2265 ", string_score_threshold, ")"),
        subtitle = subtitle_pos,
        out_path = file.path(plot_dir, paste0("module_", file_id, "__PosFrac_boxplot.png")),
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi,
        max_points = 3000L,
        seed = module_seed + 103L,
        value_is_percent = TRUE
      )
    } else {
      save_empty_plot(
        path = file.path(plot_dir, paste0("module_", file_id, "__Eavg.png")),
        title = paste0(title_base, ": E_avg (", status, ")"),
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi
      )
      save_empty_plot(
        path = file.path(plot_dir, paste0("module_", file_id, "__PosFrac.png")),
        title = paste0(title_base, ": positive edge % (", status, ")"),
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi
      )
      save_empty_plot(
        path = file.path(plot_dir, paste0("module_", file_id, "__Eavg_boxplot.png")),
        title = paste0(title_base, ": E_avg box (", status, ")"),
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi
      )
      save_empty_plot(
        path = file.path(plot_dir, paste0("module_", file_id, "__PosFrac_boxplot.png")),
        title = paste0(title_base, ": positive edge % box (", status, ")"),
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi
      )
    }

    module_rows[[length(module_rows) + 1L]] <- data.table(
      method = as.character(method),
      module_id = as.character(mid),
      module_size = as.integer(module_size),
      background_size = as.integer(background_size),
      n_pairs_total = n_pairs_total,
      n_edges_target = as.integer(n_edges_target),
      cross_edge_space = as.numeric(cross_edge_space),
      E_sum_obs = e_sum_obs,
      E_avg_obs = e_avg_obs,
      n_pairs_mapped_obs = n_pairs_mapped_obs,
      mapped_pair_frac_obs = mapped_pair_frac_obs,
      n_pairs_pos_obs = n_pairs_pos_obs,
      pos_pair_frac_obs = pos_pair_frac_obs,
      null_mean_Eavg = null_mean_eavg,
      null_sd_Eavg = null_sd_eavg,
      delta_Eavg = delta_eavg,
      empirical_p_Eavg = empirical_p_eavg,
      q_Eavg = NA_real_,
      null_mean_pos_pair_frac = null_mean_pos_pair_frac,
      null_sd_pos_pair_frac = null_sd_pos_pair_frac,
      delta_pos_pair_frac = delta_pos_pair_frac,
      empirical_p_pos_pair_frac = empirical_p_pos_pair_frac,
      q_pos_pair_frac = NA_real_,
      n_valid_null_draws = as.integer(n_valid_null_draws),
      status = status,
      status_reason = status_reason
    )
  }
}

out_dt <- rbindlist(module_rows, use.names = TRUE, fill = TRUE)
if (!nrow(out_dt)) {
  warning("No modules scored. Outputs will be empty.")
  out_dt <- data.table(
    method = character(),
    module_id = character(),
    module_size = integer(),
    background_size = integer(),
    n_pairs_total = numeric(),
    n_edges_target = integer(),
    cross_edge_space = numeric(),
    E_sum_obs = numeric(),
    E_avg_obs = numeric(),
    n_pairs_mapped_obs = integer(),
    mapped_pair_frac_obs = numeric(),
    n_pairs_pos_obs = integer(),
    pos_pair_frac_obs = numeric(),
    null_mean_Eavg = numeric(),
    null_sd_Eavg = numeric(),
    delta_Eavg = numeric(),
    empirical_p_Eavg = numeric(),
    q_Eavg = numeric(),
    null_mean_pos_pair_frac = numeric(),
    null_sd_pos_pair_frac = numeric(),
    delta_pos_pair_frac = numeric(),
    empirical_p_pos_pair_frac = numeric(),
    q_pos_pair_frac = numeric(),
    n_valid_null_draws = integer(),
    status = character(),
    status_reason = character()
  )
}

out_dt[, q_Eavg := NA_real_]
for (method in unique(out_dt$method)) {
  idx <- which(out_dt$method == method & is.finite(out_dt$empirical_p_Eavg))
  if (!length(idx)) next
  out_dt$q_Eavg[idx] <- stats::p.adjust(out_dt$empirical_p_Eavg[idx], method = "BH")
}

out_dt[, q_pos_pair_frac := NA_real_]
for (method in unique(out_dt$method)) {
  idx <- which(out_dt$method == method & is.finite(out_dt$empirical_p_pos_pair_frac))
  if (!length(idx)) next
  out_dt$q_pos_pair_frac[idx] <- stats::p.adjust(out_dt$empirical_p_pos_pair_frac[idx], method = "BH")
}

out_path_summary <- file.path(opt$outdir, "module_enrichment_anchor_cross_null_summary.tsv")
fwrite(out_dt, out_path_summary, sep = "\t")

plot_dir <- file.path(opt$outdir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

method_levels <- unique(vapply(methods, display_method_name, character(1)))
color_map <- bench_base_color_map(method_levels)
obs_line_color_map <- stats::setNames(
  vapply(color_map, boost_color_vividness, character(1), vividness = 0.25),
  names(color_map)
)

null_draw_plot_dt <- rbindlist(null_draw_samples, use.names = TRUE, fill = TRUE)
if (nrow(null_draw_plot_dt) && plot_max_null_points_per_method > 0L) {
  null_draw_plot_dt[, method_label := factor(method_label, levels = method_levels)]
  if (is.finite(plot_max_null_points_per_method) && plot_max_null_points_per_method > 0L) {
    downsampled <- list()
    for (ml in method_levels) {
      rows <- which(null_draw_plot_dt$method_label == ml)
      if (!length(rows)) next
      if (length(rows) <= plot_max_null_points_per_method) {
        downsampled[[length(downsampled) + 1L]] <- null_draw_plot_dt[rows]
      } else {
        sel <- sample_int_deterministic(length(rows), plot_max_null_points_per_method, seed = seed + 97L + match(ml, method_levels))
        downsampled[[length(downsampled) + 1L]] <- null_draw_plot_dt[rows[sel]]
      }
    }
    null_draw_plot_dt <- rbindlist(downsampled, use.names = TRUE, fill = TRUE)
  }
}

plot_dt_eavg_obs <- out_dt[is.finite(E_avg_obs)]
plot_dt_null_eavg <- null_draw_plot_dt[is.finite(null_eavg)]
plot_dt_null_mean_eavg <- out_dt[is.finite(null_mean_Eavg)]

score_ylim <- NULL
score_vals <- c(
  plot_dt_eavg_obs$E_avg_obs,
  plot_dt_null_eavg$null_eavg,
  plot_dt_null_mean_eavg$null_mean_Eavg,
  out_dt$delta_Eavg
)
score_vals <- suppressWarnings(as.numeric(score_vals))
score_vals <- score_vals[is.finite(score_vals)]
if (length(score_vals)) {
  score_ylim <- range(score_vals)
  if (isTRUE(all.equal(score_ylim[1], score_ylim[2]))) {
    pad <- if (score_ylim[2] == 0) 1 else abs(score_ylim[2]) * 0.05
    score_ylim <- c(score_ylim[1] - pad, score_ylim[2] + pad)
  }
}

plot_dt_delta <- out_dt[is.finite(delta_Eavg)]
if (nrow(plot_dt_delta)) {
  plot_dt_delta[, method_label := vapply(method, display_method_name, character(1))]
  plot_dt_delta[, method_label := factor(method_label, levels = method_levels)]
  p_delta <- ggplot(plot_dt_delta, aes(x = method_label, y = delta_Eavg, fill = method_label)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    guides(fill = "none") +
    theme_bench_like_stringdb() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Delta score (obs - null mean)",
      x = "Method",
      y = "Observed - null mean"
    )
  if (!is.null(score_ylim)) {
    p_delta <- p_delta + coord_cartesian(ylim = score_ylim)
  }
  ggsave(
    filename = file.path(plot_dir, "delta_Eavg_boxplot_by_method.png"),
    plot = p_delta,
    width = 6,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "delta_Eavg_boxplot_by_method.png"),
                  "Delta score (obs - null mean)",
    6,
    4,
    150
  )
}

if (nrow(plot_dt_eavg_obs)) {
  plot_dt_eavg_obs[, method_label := vapply(method, display_method_name, character(1))]
  plot_dt_eavg_obs[, method_label := factor(method_label, levels = method_levels)]
  p_eavg_obs <- ggplot(plot_dt_eavg_obs, aes(x = method_label, y = E_avg_obs, fill = method_label)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    guides(fill = "none") +
    theme_bench_like_stringdb() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Observed score",
      x = "Method",
      y = mean_score_label
    )
  if (!is.null(score_ylim)) {
    p_eavg_obs <- p_eavg_obs + coord_cartesian(ylim = score_ylim)
  }
  ggsave(
    filename = file.path(plot_dir, "Eavg_obs_boxplot_by_method.png"),
    plot = p_eavg_obs,
    width = 6,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "Eavg_obs_boxplot_by_method.png"),
    "Observed score",
    6,
    4,
    150
  )
}

if (nrow(plot_dt_null_eavg) && plot_max_null_points_per_method > 0L) {
  p_null_eavg <- ggplot(plot_dt_null_eavg, aes(x = method_label, y = null_eavg, fill = method_label)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.7, na.rm = TRUE) +
    geom_point(
      aes(color = method_label),
      position = position_jitter(width = 0.15, height = 0),
      size = 0.25,
      alpha = 0.25,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_color_manual(values = color_map, drop = FALSE) +
    guides(fill = "none", color = "none") +
    theme_bench_like_stringdb() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Null score (permutation)",
      x = "Method",
      y = paste0("Null ", mean_score_label)
    )
  if (!is.null(score_ylim)) {
    p_null_eavg <- p_null_eavg + coord_cartesian(ylim = score_ylim)
  }
  ggsave(
    filename = file.path(plot_dir, "null_Eavg_boxplot_by_method.png"),
    plot = p_null_eavg,
    width = 6,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "null_Eavg_boxplot_by_method.png"),
    "Null score (permutation)",
    6,
    4,
    150
  )
}

plot_dt_null_pos <- null_draw_plot_dt[is.finite(null_pos_pair_frac)]
box_pos_dt <- rbindlist(null_box_stats_pos, use.names = TRUE, fill = TRUE)
if (nrow(plot_dt_null_pos)) {
  plot_dt_null_pos[, null_pos_pair_pct := 100 * null_pos_pair_frac]
}
if (nrow(box_pos_dt)) {
  box_pos_dt[, c("ymin", "lower", "middle", "upper", "ymax", "observed") := lapply(
    .SD,
    function(x) 100 * suppressWarnings(as.numeric(x))
  ), .SDcols = c("ymin", "lower", "middle", "upper", "ymax", "observed")]
}

pos_ylim <- NULL
pos_vals <- c(
  100 * out_dt$pos_pair_frac_obs,
  plot_dt_null_pos$null_pos_pair_pct,
  100 * out_dt$null_mean_pos_pair_frac
)
if (nrow(box_pos_dt)) {
  pos_vals <- c(pos_vals, box_pos_dt$ymin, box_pos_dt$ymax, box_pos_dt$observed)
}
pos_vals <- suppressWarnings(as.numeric(pos_vals))
pos_vals <- pos_vals[is.finite(pos_vals)]
if (length(pos_vals)) {
  pos_ylim <- range(pos_vals)
  pos_ylim[1] <- max(0, pos_ylim[1])
  pos_ylim[2] <- min(100, pos_ylim[2])
  if (isTRUE(all.equal(pos_ylim[1], pos_ylim[2]))) {
    pad <- if (pos_ylim[2] == 0) 1 else abs(pos_ylim[2]) * 0.05
    pos_ylim <- c(max(0, pos_ylim[1] - pad), min(100, pos_ylim[2] + pad))
  }
}
if (nrow(plot_dt_null_pos) && plot_max_null_points_per_method > 0L) {
  p_null_pos <- ggplot(plot_dt_null_pos, aes(x = method_label, y = null_pos_pair_pct, fill = method_label)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.7, na.rm = TRUE) +
    geom_point(
      aes(color = method_label),
      position = position_jitter(width = 0.15, height = 0),
      size = 0.25,
      alpha = 0.25,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_y_continuous(labels = percent_label_0_100) +
    guides(fill = "none", color = "none") +
    theme_bench_like_stringdb() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Null positive edge % (permutation)",
      x = "Method",
      y = "Positive edge %"
    )
  if (!is.null(pos_ylim)) {
    p_null_pos <- p_null_pos + coord_cartesian(ylim = pos_ylim)
  }
  ggsave(
    filename = file.path(plot_dir, "null_pos_pair_frac_boxplot_by_method.png"),
    plot = p_null_pos,
    width = 6,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "null_pos_pair_frac_boxplot_by_method.png"),
    "Null positive edge % (permutation)",
    6,
    4,
    150
  )
}

if (nrow(plot_dt_null_mean_eavg)) {
  plot_dt_null_mean_eavg[, method_label := vapply(method, display_method_name, character(1))]
  plot_dt_null_mean_eavg[, method_label := factor(method_label, levels = method_levels)]
  p_null_mean_eavg <- ggplot(plot_dt_null_mean_eavg, aes(x = method_label, y = null_mean_Eavg, fill = method_label)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    guides(fill = "none") +
    theme_bench_like_stringdb() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Null mean score",
      x = "Method",
      y = paste0("Null mean ", mean_score_label)
    )
  if (!is.null(score_ylim)) {
    p_null_mean_eavg <- p_null_mean_eavg + coord_cartesian(ylim = score_ylim)
  }
  ggsave(
    filename = file.path(plot_dir, "null_mean_Eavg_boxplot_by_method.png"),
    plot = p_null_mean_eavg,
    width = 6,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "null_mean_Eavg_boxplot_by_method.png"),
    "Null mean score",
    6,
    4,
    150
  )
}

plot_dt_null_mean_pos <- out_dt[is.finite(null_mean_pos_pair_frac)]
if (nrow(plot_dt_null_mean_pos)) {
  plot_dt_null_mean_pos[, method_label := vapply(method, display_method_name, character(1))]
  plot_dt_null_mean_pos[, method_label := factor(method_label, levels = method_levels)]
  plot_dt_null_mean_pos[, null_mean_pos_pair_pct := 100 * null_mean_pos_pair_frac]
  p_null_mean_pos <- ggplot(plot_dt_null_mean_pos, aes(x = method_label, y = null_mean_pos_pair_pct, fill = method_label)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_y_continuous(labels = percent_label_0_100) +
    guides(fill = "none") +
    theme_bench_like_stringdb() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Null mean positive edge %",
      x = "Method",
      y = "Null mean positive edge %"
    )
  if (!is.null(pos_ylim)) {
    p_null_mean_pos <- p_null_mean_pos + coord_cartesian(ylim = pos_ylim)
  }
  ggsave(
    filename = file.path(plot_dir, "null_mean_pos_pair_frac_boxplot_by_method.png"),
    plot = p_null_mean_pos,
    width = 6,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "null_mean_pos_pair_frac_boxplot_by_method.png"),
    "Null mean positive edge %",
    6,
    4,
    150
  )
}

box_eavg_dt <- rbindlist(null_box_stats_eavg, use.names = TRUE, fill = TRUE)
if (nrow(box_eavg_dt)) {
  box_eavg_dt[, method_label := factor(method_label, levels = method_levels)]
  setorder(box_eavg_dt, method_label, module_id_num, module_id)
  box_eavg_dt[, module_label := paste0(as.character(method_label), "-", module_id)]
  box_eavg_dt[, module_label := factor(module_label, levels = unique(module_label))]
  n_boxes <- nlevels(box_eavg_dt$module_label)
  plot_w <- max(6, min(24, 4 + 0.18 * n_boxes))

  p_mod_eavg <- ggplot(box_eavg_dt, aes(x = module_label, fill = method_label)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity",
      alpha = 0.7,
      linewidth = 0.3,
      color = "grey40",
      na.rm = TRUE
    ) +
    geom_errorbar(
      data = box_eavg_dt,
      mapping = aes(x = module_label, ymin = observed, ymax = observed, color = "Observed"),
      inherit.aes = FALSE,
      width = 0.6,
      linewidth = 0.6,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_color_manual(values = c("Observed" = "black"), drop = FALSE) +
    guides(fill = "none", color = guide_legend(title = NULL)) +
    theme_bench_like_stringdb() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key = element_blank()
    ) +
    labs(
      title = "Null score by module",
      x = "Method / module",
      y = mean_score_label
    )
  ggsave(
    filename = file.path(plot_dir, "null_Eavg_boxplot_by_method_module.png"),
    plot = p_mod_eavg,
    width = plot_w,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "null_Eavg_boxplot_by_method_module.png"),
    "Null score by module",
    10,
    4,
    150
  )
}

if (nrow(box_pos_dt)) {
  box_pos_dt[, method_label := factor(method_label, levels = method_levels)]
  setorder(box_pos_dt, method_label, module_id_num, module_id)
  box_pos_dt[, module_label := paste0(as.character(method_label), "-", module_id)]
  box_pos_dt[, module_label := factor(module_label, levels = unique(module_label))]
  n_boxes <- nlevels(box_pos_dt$module_label)
  plot_w <- max(6, min(24, 4 + 0.18 * n_boxes))

  p_mod_pos <- ggplot(box_pos_dt, aes(x = module_label, fill = method_label)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity",
      alpha = 0.7,
      linewidth = 0.3,
      color = "grey40",
      na.rm = TRUE
    ) +
    geom_errorbar(
      data = box_pos_dt,
      mapping = aes(x = module_label, ymin = observed, ymax = observed, color = "Observed"),
      inherit.aes = FALSE,
      width = 0.6,
      linewidth = 0.6,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_color_manual(values = c("Observed" = "black"), drop = FALSE) +
    scale_y_continuous(labels = percent_label_0_100) +
    guides(fill = "none", color = guide_legend(title = NULL)) +
    theme_bench_like_stringdb() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key = element_blank()
    ) +
    labs(
      title = "Null positive edge % by module",
      x = "Method / module",
      y = "Positive edge %"
    )
  if (!is.null(pos_ylim)) {
    p_mod_pos <- p_mod_pos + coord_cartesian(ylim = pos_ylim)
  }
  ggsave(
    filename = file.path(plot_dir, "null_pos_pair_frac_boxplot_by_method_module.png"),
    plot = p_mod_pos,
    width = plot_w,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "null_pos_pair_frac_boxplot_by_method_module.png"),
    "Null positive edge % by module",
    10,
    4,
    150
  )
}

plot_dt_pos_obs <- out_dt[is.finite(pos_pair_frac_obs)]
if (nrow(plot_dt_pos_obs)) {
  plot_dt_pos_obs[, method_label := vapply(method, display_method_name, character(1))]
  plot_dt_pos_obs[, method_label := factor(method_label, levels = method_levels)]
  plot_dt_pos_obs[, pos_pair_pct_obs := 100 * pos_pair_frac_obs]
  p_pos_obs <- ggplot(plot_dt_pos_obs, aes(x = method_label, y = pos_pair_pct_obs, fill = method_label)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    scale_y_continuous(labels = percent_label_0_100) +
    guides(fill = "none") +
    theme_bench_like_stringdb() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Observed positive edge %",
      x = "Method",
      y = "Positive edge %"
    )
  if (!is.null(pos_ylim)) {
    p_pos_obs <- p_pos_obs + coord_cartesian(ylim = pos_ylim)
  }
  ggsave(
    filename = file.path(plot_dir, "pos_pair_frac_obs_boxplot_by_method.png"),
    plot = p_pos_obs,
    width = 6,
    height = 4,
    dpi = 150
  )
} else {
  save_empty_plot(
    file.path(plot_dir, "pos_pair_frac_obs_boxplot_by_method.png"),
    "Observed positive edge %",
    6,
    4,
    150
  )
}

run_meta <- list(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  map_dir = opt$map_dir,
  map_path = map_path,
  map_sha256 = sha256_file(map_path),
  map_meta_path = if (file.exists(map_meta_path)) map_meta_path else "",
  map_meta = map_meta,
  outdir = opt$outdir,
  methods = methods,
  modules_tsv_by_method = module_paths,
  modules_sha256_by_method = module_hashes,
  params = list(
    n_random = n_random,
    seed = seed,
    min_module_genes = min_module_genes,
    min_valid_null_draws = min_valid_null_draws,
    max_resample_attempts = max_resample_attempts,
    plot_max_modules_per_method = plot_max_modules_per_method,
    plot_max_null_points_per_method = plot_max_null_points_per_method,
    string_score_threshold = string_score_threshold
  ),
  background_size = background_size,
  null_definition = paste0(
    "anchor_cross_edges: sample k=choose(m,2) edges with one endpoint in module and the other in background(outside), ",
    "without replacement; missing edges treated as 0"
  ),
  notes = list(
    observed_definition = "module_internal_pairs: fixed denom choose(m,2); missing edges treated as 0",
    null_k = "k = choose(module_size,2)",
    cross_edge_space = "m*b where m=module_size and b=outside background size",
    module_seed_strategy = "module_seed = seed + module_counter * 104729"
  ),
  git_commit = get_git_commit()
)

out_path_meta <- file.path(opt$outdir, "meta.json")
write_json(run_meta, out_path_meta, auto_unbox = TRUE, pretty = TRUE, null = "null")

wall_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
message(sprintf("[DONE] Anchor-cross enrichment completed. wall_sec=%.2f", wall_sec))

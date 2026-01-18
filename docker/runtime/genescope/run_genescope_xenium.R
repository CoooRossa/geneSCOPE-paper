args <- commandArgs(trailingOnly = TRUE)

data_dir  <- args[[1]]
outdir    <- args[[2]]
grid_um   <- as.integer(args[[3]])
ncores    <- as.integer(args[[4]])
seed      <- as.integer(args[[5]])
parallel_backend <- as.character(args[[6]])
coord_csv <- ifelse(length(args) >= 7 && nzchar(args[[7]]), args[[7]], NA_character_)
dataset_id <- ifelse(length(args) >= 8 && nzchar(args[[8]]), args[[8]], basename(data_dir))
roi_id <- ifelse(
  length(args) >= 9 && nzchar(args[[9]]),
  args[[9]],
  ifelse(is.na(coord_csv), "full", tools::file_path_sans_ext(basename(coord_csv)))
)
stats_tsv <- ifelse(
  length(args) >= 10 && nzchar(args[[10]]),
  args[[10]],
  file.path(outdir, "stats.tsv")
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(geneSCOPE)
  library(arrow)
  library(sf)
  library(jsonlite)
  library(Matrix)
  library(data.table)
})

set.seed(seed)

start_time <- Sys.time()

cat("[run] data_dir=", data_dir, "\n", sep = "")
cat("[run] grid_um=", grid_um, " ncores=", ncores, " seed=", seed, " backend=", parallel_backend, "\n", sep = "")

grid_name <- paste0("grid", grid_um)

read_cells_parquet <- function(data_dir) {
  cells_pq <- file.path(data_dir, "cells.parquet")
  if (!file.exists(cells_pq)) stop("Missing: ", cells_pq)

  cells_dt <- as.data.frame(arrow::read_parquet(cells_pq))
  lc <- setNames(names(cells_dt), tolower(names(cells_dt)))
  idcol <- lc[["cell_id"]]
  if (is.null(idcol)) idcol <- lc[["barcode"]]
  if (is.null(idcol)) idcol <- lc[["id"]]
  if (is.null(idcol)) stop("cells.parquet: cannot find cell id column.")

  xcol <- lc[["x_centroid"]]
  if (is.null(xcol)) xcol <- lc[["x"]]
  ycol <- lc[["y_centroid"]]
  if (is.null(ycol)) ycol <- lc[["y"]]
  if (is.null(xcol) || is.null(ycol)) stop("cells.parquet: cannot find x/y columns.")

  data.frame(
    cell_id = as.character(cells_dt[[idcol]]),
    x = as.numeric(cells_dt[[xcol]]),
    y = as.numeric(cells_dt[[ycol]]),
    stringsAsFactors = FALSE
  )
}

read_roi_polygon <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) > 0) lines[1] <- sub("^\\ufeff", "", lines[1])

  lines2 <- lines[!grepl("^\\s*#", lines)]
  lines2 <- lines2[nzchar(trimws(lines2))]
  if (length(lines2) < 4) stop("ROI CSV after dropping comment lines is too short: ", path)

  df <- utils::read.csv(text = paste(lines2, collapse = "\n"), stringsAsFactors = FALSE)
  cn <- tolower(names(df))
  xcol <- which(cn %in% c("x"))
  ycol <- which(cn %in% c("y"))
  if (length(xcol) != 1 || length(ycol) != 1) {
    stop(sprintf("ROI CSV must have columns X,Y (case-insensitive). cols=%s",
                 paste(names(df), collapse = ",")))
  }

  xx <- suppressWarnings(as.numeric(df[[xcol]]))
  yy <- suppressWarnings(as.numeric(df[[ycol]]))
  keep <- is.finite(xx) & is.finite(yy)
  pts <- cbind(x = xx[keep], y = yy[keep])

  if (nrow(pts) >= 2) {
    dxy <- abs(pts[-1, , drop = FALSE] - pts[-nrow(pts), , drop = FALSE])
    same_prev <- c(FALSE, rowSums(dxy) < 1e-9)
    pts <- pts[!same_prev, , drop = FALSE]
  }

  if (nrow(pts) < 3) stop("ROI polygon needs >=3 vertices after cleaning.")

  if (sum(abs(pts[1, ] - pts[nrow(pts), ])) > 1e-9) {
    pts <- rbind(pts, pts[1, , drop = FALSE])
  }

  pts
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

sha256_text <- function(text) {
  tf <- tempfile()
  writeLines(text, tf, useBytes = TRUE)
  on.exit(unlink(tf), add = TRUE)
  sha256_file(tf)
}

read_resource_stats <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(list())
  df <- tryCatch(utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(list(stats_path = path))

  cols <- names(df)
  pick_col <- function(candidates) {
    hit <- intersect(candidates, cols)
    if (length(hit) == 0) return(NULL)
    hit[[1]]
  }

  mem_cur_col <- pick_col(c("mem_current_bytes", "memory_current_bytes", "mem_bytes"))
  mem_peak_col <- pick_col(c("mem_peak_bytes", "memory_peak_bytes", "mem_max_bytes", "memory_max_bytes"))
  cpu_col <- pick_col(c("cpu_usage_usec"))

  res <- list(stats_path = path)
  if (!is.null(mem_cur_col)) {
    v <- suppressWarnings(as.numeric(df[[mem_cur_col]]))
    if (any(is.finite(v))) res$max_mem_current_bytes <- max(v, na.rm = TRUE)
  }
  if (!is.null(mem_peak_col)) {
    v <- suppressWarnings(as.numeric(df[[mem_peak_col]]))
    if (any(is.finite(v))) res$max_mem_peak_bytes <- max(v, na.rm = TRUE)
  } else if (!is.null(mem_cur_col) && !is.null(res$max_mem_current_bytes)) {
    res$max_mem_peak_bytes <- res$max_mem_current_bytes
  }
  if (!is.null(cpu_col)) {
    v <- suppressWarnings(as.numeric(df[[cpu_col]]))
    v <- v[is.finite(v)]
    if (length(v) > 0) res$cpu_usage_usec <- v[[length(v)]]
  }
  res
}

clip_cells_to_roi <- function(cells_df, roi_vertices) {
  poly <- sf::st_polygon(list(as.matrix(roi_vertices)))
  roi <- sf::st_sfc(poly)
  sf::st_crs(roi) <- NA

  pts <- sf::st_as_sf(cells_df, coords = c("x", "y"), remove = FALSE, crs = NA)
  inside <- sf::st_within(pts, roi, sparse = FALSE)[, 1]
  cells_df[inside, , drop = FALSE]
}

roi_hash <- ""
roi_source_sha256 <- ""
obs_id_hash <- ""
if (!is.na(coord_csv) && file.exists(coord_csv)) {
  roi_source_sha256 <- sha256_file(coord_csv)
}

n_obs_input <- NA_integer_
n_obs_raw <- NA_integer_
n_obs_roi <- NA_integer_
roi_applied <- 0L
if (file.exists(file.path(data_dir, "cells.parquet"))) {
  cells_df <- read_cells_parquet(data_dir)
  n_obs_raw <- nrow(cells_df)
  n_obs_input <- n_obs_raw
  if (!is.na(coord_csv) && file.exists(coord_csv)) {
    roi <- read_roi_polygon(coord_csv)
    cells_df <- clip_cells_to_roi(cells_df, roi)
    n_obs_roi <- nrow(cells_df)
    n_obs_input <- n_obs_roi
    roi_applied <- 1L
  }
  obs_id_hash <- sha256_text(paste(sort(as.character(cells_df$cell_id)), collapse = "\n"))
}

if (nzchar(obs_id_hash)) {
  roi_hash <- sha256_text(paste0(roi_source_sha256, "\n", obs_id_hash))
}

if (is.na(n_obs_raw)) n_obs_raw <- n_obs_input
if (is.na(n_obs_roi)) n_obs_roi <- n_obs_input

# NOTE:
# - parallel_backend is used only for ROI clipping stages (per your patch).
# - In Docker, default serial is the safest. You can override via --parallel_backend.
scope_obj <- createSCOPE(
  data_dir = data_dir,
  grid_length = grid_um,
  seg_type = "cell",
  coord_file = ifelse(is.na(coord_csv), NULL, coord_csv),
  ncores = ncores,
  parallel_backend = parallel_backend
)

scope_obj <- normalizeMoleculesInGrid(scope_obj = scope_obj, grid_name = grid_name)

scope_obj <- computeWeights(scope_obj = scope_obj, grid_name = grid_name, ncores = ncores)

scope_obj <- computeL(scope_obj = scope_obj, grid_name = grid_name, use_bigmemory = FALSE, ncores = ncores)

scope_obj <- clusterGenes(
  scope_obj = scope_obj,
  grid_name = grid_name,
  pct_min = "q95",
  algo = "leiden",
  resolution = 0.1,
  ncores = ncores,
  consensus_thr = 0.95,
  n_restart = 200
)

membership_col <- "modL0.00"
genes <- rownames(scope_obj@meta.data)
membership <- scope_obj@meta.data[[membership_col]]
if (is.null(membership)) {
  membership <- rep(NA_integer_, length(genes))
}
if (is.factor(membership)) {
  membership <- as.character(membership)
}
module_id <- suppressWarnings(as.integer(membership))
module_id[is.na(module_id)] <- -1L

modules_df <- data.frame(
  gene = genes,
  module_id = module_id,
  stringsAsFactors = FALSE
)
utils::write.table(
  modules_df,
  file.path(outdir, "modules.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

extract_lee_stats <- function(scope_obj, grid_name) {
  if (is.null(scope_obj@stats) || is.null(scope_obj@stats[[grid_name]])) {
    stop("Missing scope_obj@stats[[", grid_name, "]] after computeL()")
  }
  stats_grid <- scope_obj@stats[[grid_name]]
  layer_names <- names(stats_grid)
  if (is.null(layer_names) || length(layer_names) == 0) {
    stop("Empty scope_obj@stats[[", grid_name, "]]")
  }

  layer_names <- layer_names[!grepl("^_", layer_names)]
  if (length(layer_names) == 0) {
    stop("No stats layers found in scope_obj@stats[[", grid_name, "]]")
  }

  lee_layers <- layer_names[grepl("^LeeStats_", layer_names)]
  layer_name <- if (length(lee_layers) > 0) lee_layers[[1]] else layer_names[[1]]

  lee <- stats_grid[[layer_name]]
  if (is.null(lee) || is.null(lee$L) || is.null(lee$FDR)) {
    stop("Lee stats layer missing L/FDR: ", layer_name)
  }
  list(
    layer_name = layer_name,
    L = lee$L,
    FDR = lee$FDR,
    P = lee$P,
    Z = lee$Z,
    meta = if (!is.null(lee$meta)) lee$meta else list()
  )
}

infer_grid_n_obs <- function(scope_obj, grid_name) {
  if (is.null(scope_obj@grid) || is.null(scope_obj@grid[[grid_name]])) return(NA_integer_)
  g <- scope_obj@grid[[grid_name]]
  for (nm in c("Xz", "counts", "raw_counts", "expr", "X", "data", "logCPM")) {
    x <- g[[nm]]
    if (is.null(x)) next
    d <- tryCatch(dim(x), error = function(e) NULL)
    if (!is.null(d) && length(d) == 2) return(as.integer(d[[1]]))
  }
  if (!is.null(g$grid_info)) {
    d2 <- tryCatch(dim(g$grid_info), error = function(e) NULL)
    if (!is.null(d2) && length(d2) >= 1) return(as.integer(d2[[1]]))
  }
  NA_integer_
}

subset_square_matrix <- function(mat, idx) {
  if (inherits(mat, "big.matrix")) {
    return(as.matrix(mat[idx, idx, drop = FALSE]))
  }
  if (inherits(mat, "Matrix")) {
    return(as.matrix(mat[idx, idx, drop = FALSE]))
  }
  if (is.matrix(mat)) return(mat[idx, idx, drop = FALSE])
  as.matrix(mat)[idx, idx, drop = FALSE]
}

emit_gene_graph_edges <- function(lee, genes, module_genes, n_obs, outdir, edge_type = "lee_L") {
  out_path <- file.path(outdir, "gene_graph_edges.tsv")
  empty <- data.table::data.table(
    gene_a = character(),
    gene_b = character(),
    weight = numeric(),
    p_value = numeric(),
    fdr = numeric(),
    n_obs = integer(),
    edge_type = character()
  )

  genes <- as.character(genes)
  module_genes <- as.character(module_genes)

  if (is.null(genes) || length(genes) < 2 || is.null(module_genes) || length(module_genes) < 2) {
    data.table::fwrite(empty, out_path, sep = "\t")
    return(list(
      edges_n_total = 0L,
      edges_n_obs = as.integer(n_obs),
      edge_type = edge_type,
      edges_gene_set = "module_genes",
      edges_genes_n = 0L
    ))
  }

  idx_genes <- match(module_genes, genes)
  idx_genes <- idx_genes[!is.na(idx_genes)]
  idx_genes <- sort(unique(as.integer(idx_genes)))
  if (length(idx_genes) < 2) {
    data.table::fwrite(empty, out_path, sep = "\t")
    return(list(
      edges_n_total = 0L,
      edges_n_obs = as.integer(n_obs),
      edge_type = edge_type,
      edges_gene_set = "module_genes",
      edges_genes_n = 0L
    ))
  }

  gene_sub <- genes[idx_genes]

  Lm <- subset_square_matrix(lee$L, idx_genes)
  Fm <- subset_square_matrix(lee$FDR, idx_genes)
  if (is.null(dim(Lm)) || nrow(Lm) != length(gene_sub) || ncol(Lm) != length(gene_sub)) {
    stop("gene_graph_edges: L dims mismatch with module genes")
  }
  if (is.null(dim(Fm)) || nrow(Fm) != length(gene_sub) || ncol(Fm) != length(gene_sub)) {
    stop("gene_graph_edges: FDR dims mismatch with module genes")
  }

  idx <- which(upper.tri(Lm, diag = FALSE), arr.ind = TRUE)
  if (nrow(idx) == 0) {
    data.table::fwrite(empty, out_path, sep = "\t")
    return(list(
      edges_n_total = 0L,
      edges_n_obs = as.integer(n_obs),
      edge_type = edge_type,
      edges_gene_set = "module_genes",
      edges_genes_n = as.integer(length(gene_sub))
    ))
  }

  pairs <- cbind(idx[, 1], idx[, 2])
  weight <- as.numeric(Lm[pairs])
  fdr <- as.numeric(Fm[pairs])

  p_value <- NULL
  if (!is.null(lee$P)) {
    Pm <- subset_square_matrix(lee$P, idx_genes)
    p_value <- as.numeric(Pm[pairs])
  } else if (!is.null(lee$Z)) {
    Zm <- subset_square_matrix(lee$Z, idx_genes)
    z <- as.numeric(Zm[pairs])
    p_value <- 2 * stats::pnorm(-abs(z))
  } else {
    p_value <- rep(NA_real_, length(weight))
  }

  for (v in list(weight = weight, p_value = p_value, fdr = fdr)) {
    # no-op (forces eager eval)
  }
  p_value[!is.finite(p_value)] <- NA_real_
  fdr[!is.finite(fdr)] <- NA_real_
  weight[!is.finite(weight)] <- NA_real_
  p_value[p_value < 0] <- 0
  p_value[p_value > 1] <- 1
  fdr[fdr < 0] <- 0
  fdr[fdr > 1] <- 1

  edges <- data.table::data.table(
    gene_a = gene_sub[idx[, 1]],
    gene_b = gene_sub[idx[, 2]],
    weight = weight,
    p_value = p_value,
    fdr = fdr,
    n_obs = as.integer(n_obs),
    edge_type = edge_type
  )
  data.table::fwrite(edges, out_path, sep = "\t")
  list(
    edges_n_total = nrow(edges),
    edges_n_obs = as.integer(n_obs),
    edge_type = edge_type,
    edges_gene_set = "module_genes",
    edges_genes_n = as.integer(length(gene_sub))
  )
}

lee <- extract_lee_stats(scope_obj, grid_name)
gg_genes <- if (!is.null(rownames(lee$L))) rownames(lee$L) else rownames(scope_obj@meta.data)
gg_genes_n <- length(gg_genes)

n_obs_grid <- infer_grid_n_obs(scope_obj, grid_name)
module_genes <- modules_df$gene[modules_df$module_id > 0]
edges_out <- emit_gene_graph_edges(lee, gg_genes, module_genes, n_obs_grid, outdir, edge_type = "lee_L")

method_version <- tryCatch(as.character(utils::packageVersion("geneSCOPE")), error = function(e) "unknown")
timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
run_id <- paste0("genescope_", format(Sys.time(), "%Y%m%dT%H%M%S"))

n_vars_used <- nrow(scope_obj@meta.data)
n_vars_input <- n_vars_used
n_obs_used <- n_obs_input

resources <- list()
if (!is.na(stats_tsv) && nzchar(stats_tsv)) {
  resources <- read_resource_stats(stats_tsv)
}

wall_time_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

meta <- list(
  method = "genescope",
  method_version = method_version,
  run_id = run_id,
  timestamp = timestamp,
  gg_weight_type = "Lee_L",
  gg_fdr_source = "native_api",
  gg_dim = list(n = gg_genes_n),
  gg_genes_n = gg_genes_n,
  gg_outputs_emitted = FALSE,
  edges_all_fullgene_disabled = TRUE,
  edges_file = "gene_graph_edges.tsv",
  edges_schema = c("gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"),
  edges_weight_type = "Lee_L",
  edges_stats_layer = lee$layer_name,
  edges_p_source = if (!is.null(lee$meta$p_source)) lee$meta$p_source else "unknown",
  edges_p_resolution = if (!is.null(lee$meta$p_resolution)) lee$meta$p_resolution else "",
  edges_fdr_method = if (!is.null(lee$meta$FDR_main_method)) lee$meta$FDR_main_method else "unknown",
  edges_edge_type = edges_out$edge_type,
  edges_gene_set = edges_out$edges_gene_set,
  edges_genes_n = edges_out$edges_genes_n,
  edges_n_total = edges_out$edges_n_total,
  edges_n_obs = edges_out$edges_n_obs,
  input_dataset_id = dataset_id,
  roi_id = roi_id,
  roi_hash = roi_hash,
  roi_source_path = ifelse(is.na(coord_csv), "", coord_csv),
  roi_source_sha256 = roi_source_sha256,
  obs_id_hash = obs_id_hash,
  n_obs_raw = n_obs_raw,
  n_obs_roi = n_obs_roi,
  n_obs_input = n_obs_input,
  n_obs_used = n_obs_used,
  n_vars_input = n_vars_input,
  n_vars_used = n_vars_used,
  gene_filtering = list(pct_min = "q95"),
  random_seed = seed,
  stochastic = TRUE,
  params = list(
    grid_um = grid_um,
    ncores = ncores,
    parallel_backend = parallel_backend,
    membership_col = membership_col,
    roi_applied = roi_applied
  ),
  runtime = list(
    wall_time_sec = wall_time_sec,
    n_threads = ncores
  ),
  resources = resources,
  membership_type = "hard"
)

jsonlite::write_json(meta, file.path(outdir, "meta.json"), auto_unbox = TRUE, pretty = TRUE)

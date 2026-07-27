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
cluster_pct <- ifelse(length(args) >= 11 && nzchar(args[[11]]), args[[11]], "q95")
n_restart <- ifelse(length(args) >= 12 && nzchar(args[[12]]), as.integer(args[[12]]), 1000L)
perms <- ifelse(length(args) >= 13 && nzchar(args[[13]]), as.integer(args[[13]]), 1000L)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(geneSCOPE)
  library(arrow)
  library(sf)
  library(jsonlite)
  library(Matrix)
  library(data.table)
})

if (!identical(as.character(utils::packageVersion("geneSCOPE")), "1.0.2")) {
  stop("This frozen runner requires geneSCOPE 1.0.2")
}
if (!cluster_pct %in% c("q95", "q99.9")) stop("Unsupported cluster_pct: ", cluster_pct)
if (!is.finite(n_restart) || n_restart < 1L) stop("n_restart must be positive")
if (!is.finite(perms) || perms < 1L) stop("perms must be positive")

reset_rng <- function() {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
}
reset_rng()

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
    out <- tryCatch(system2(cmd, shQuote(path), stdout = TRUE, stderr = TRUE),
                    error = function(e) "")
  } else {
    cmd <- Sys.which("shasum")
    if (!nzchar(cmd)) stop("sha256sum or shasum is required for input provenance.")
    out <- tryCatch(
      system2(cmd, c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE),
      error = function(e) ""
    )
  }
  status <- attr(out, "status", exact = TRUE)
  digest <- if (length(out)) strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]] else ""
  if ((!is.null(status) && status != 0L) || !grepl("^[0-9a-fA-F]{64}$", digest)) {
    stop("Could not hash input: ", path)
  }
  tolower(digest)
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

input_files <- c(
  cell_feature_matrix = file.path(data_dir, "cell_feature_matrix.h5"),
  cells = file.path(data_dir, "cells.parquet"),
  transcripts = file.path(data_dir, "transcripts.parquet")
)
input_files <- input_files[file.exists(input_files)]
input_sha256 <- as.list(vapply(input_files, sha256_file, character(1L)))

if (is.na(n_obs_raw)) n_obs_raw <- n_obs_input
if (is.na(n_obs_roi)) n_obs_roi <- n_obs_input

# NOTE:
# - parallel_backend is used only for ROI clipping stages (per your patch).
# - In Docker, default serial is the safest. You can override via --parallel_backend.
coord_file_arg <- if (is.na(coord_csv)) NULL else coord_csv
scope_obj <- createSCOPE(
  data_dir = data_dir,
  grid_length = grid_um,
  seg_type = "cell",
  coord_file = coord_file_arg,
  ncores = ncores,
  parallel_backend = parallel_backend
)

scope_obj <- normalizeMoleculesInGrid(scope_obj = scope_obj, grid_name = grid_name)

scope_obj <- computeWeights(
  scope_obj = scope_obj,
  grid_name = grid_name,
  style = "B",
  topology = "auto",
  store_mat = TRUE,
  store_listw = TRUE,
  ncores = ncores
)

reset_rng()
scope_obj <- computeL(
  scope_obj = scope_obj,
  grid_name = grid_name,
  use_bigmemory = FALSE,
  ncores = ncores,
  perms = perms,
  use_blocks = FALSE,
  norm_layer = "Xz"
)

cluster_name <- paste0("correction_", cluster_pct, "_res0.1_", grid_name, "_freq0.95")
reset_rng()
scope_obj <- clusterGenes(
  scope_obj = scope_obj,
  grid_name = grid_name,
  pct_min = cluster_pct,
  algo = "leiden",
  resolution = 0.1,
  ncores = ncores,
  consensus_thr = 0.95,
  n_restart = n_restart,
  use_log1p_weight = TRUE,
  use_consensus = TRUE,
  cluster_name = cluster_name,
  graph_slot_name = cluster_name
)

membership_col <- cluster_name
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

write_mm <- function(mat, path) {
  if (inherits(mat, "big.matrix")) {
    mat <- as.matrix(mat[, ])
  }
  if (is.null(dim(mat)) || length(dim(mat)) != 2) {
    stop("Not a 2D matrix: ", path)
  }
  if (!inherits(mat, "Matrix")) mat <- Matrix::Matrix(mat, sparse = FALSE)

  # Matrix::writeMM() does not support dense dgeMatrix for full g×g outputs.
  # Emit MatrixMarket array format (column-major) instead.
  con <- file(path, open = "wt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  dims <- dim(mat)
  writeLines("%%MatrixMarket matrix array real general", con)
  writeLines("%", con)
  writeLines(paste(dims[[1]], dims[[2]]), con)

  vals <- if (inherits(mat, "dgeMatrix")) {
    mat@x
  } else {
    as.numeric(mat)
  }
  vals[is.na(vals)] <- NaN

  chunk_size <- 1000000L
  n <- length(vals)
  for (i in seq.int(1L, n, by = chunk_size)) {
    j <- min(n, i + chunk_size - 1L)
    cat(sprintf("%.17g", vals[i:j]), sep = "\n", file = con, append = TRUE)
    cat("\n", file = con, append = TRUE)
  }
}

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
  valid_layers <- lee_layers[vapply(lee_layers, function(nm) {
    x <- stats_grid[[nm]]
    meta <- x$meta
    is.list(meta) &&
      meta$formula_id %in% c("Lee2009_S2_v1", "Lee_S2_v1") &&
      identical(meta$use_blocks, FALSE) &&
      identical(as.integer(meta$perms), as.integer(perms)) &&
      identical(meta$permutation_scheme, "global_joint_shuffle")
  }, logical(1))]
  if (length(valid_layers) != 1L) {
    stop("Expected exactly one canonical global-shuffle Lee layer; found: ",
         paste(valid_layers, collapse = ", "))
  }
  layer_name <- valid_layers[[1L]]

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

emit_gg_matrices <- function(lee, genes, outdir) {
  if (is.null(lee) || is.null(lee$L) || is.null(lee$FDR)) stop("emit_gg_matrices: lee is missing L/FDR")

  L <- lee$L
  FDR <- lee$FDR

  if (!is.null(rownames(L))) {
    genes <- rownames(L)
  }
  if (is.null(genes) || length(genes) == 0) {
    stop("Cannot determine gene order for gg matrices")
  }

  if (is.null(dimnames(L))) dimnames(L) <- list(genes, genes)
  if (is.null(dimnames(FDR))) dimnames(FDR) <- list(genes, genes)

  gg_genes_df <- data.frame(idx = seq_along(genes), gene = genes, stringsAsFactors = FALSE)
  utils::write.table(
    gg_genes_df,
    file.path(outdir, "gg_genes.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  write_mm(L, file.path(outdir, "gg_weight_matrix.mtx"))
  write_mm(FDR, file.path(outdir, "gg_fdr_matrix.mtx"))

  list(
    gg_genes_n = length(genes),
    gg_layer_name = lee$layer_name,
    gg_genes = genes
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

as_dense_matrix <- function(mat) {
  if (inherits(mat, "big.matrix")) {
    return(as.matrix(mat[, ]))
  }
  if (inherits(mat, "Matrix")) {
    return(as.matrix(mat))
  }
  if (is.matrix(mat)) return(mat)
  as.matrix(mat)
}

emit_edges_all <- function(lee, genes, n_obs, outdir, edge_type = "lee_L") {
  out_path <- file.path(outdir, "all_edges.tsv")
  empty <- data.table::data.table(
    gene_a = character(),
    gene_b = character(),
    weight = numeric(),
    p_value = numeric(),
    fdr = numeric(),
    n_obs = integer(),
    edge_type = character()
  )

  if (is.null(genes) || length(genes) < 2) {
    data.table::fwrite(empty, out_path, sep = "\t")
    return(list(edges_n_total = 0L, edges_n_obs = as.integer(n_obs), edge_type = edge_type))
  }

  Lm <- as_dense_matrix(lee$L)
  Fm <- as_dense_matrix(lee$FDR)
  if (is.null(dim(Lm)) || nrow(Lm) != length(genes) || ncol(Lm) != length(genes)) {
    stop("edges_all: L dims mismatch with genes")
  }
  if (is.null(dim(Fm)) || nrow(Fm) != length(genes) || ncol(Fm) != length(genes)) {
    stop("edges_all: FDR dims mismatch with genes")
  }

  idx <- which(upper.tri(Lm, diag = FALSE), arr.ind = TRUE)
  if (nrow(idx) == 0) {
    data.table::fwrite(empty, out_path, sep = "\t")
    return(list(edges_n_total = 0L, edges_n_obs = as.integer(n_obs), edge_type = edge_type))
  }

  pairs <- cbind(idx[, 1], idx[, 2])
  weight <- as.numeric(Lm[pairs])
  fdr <- as.numeric(Fm[pairs])

  p_value <- NULL
  if (!is.null(lee$P)) {
    Pm <- as_dense_matrix(lee$P)
    p_value <- as.numeric(Pm[pairs])
  } else if (!is.null(lee$Z)) {
    Zm <- as_dense_matrix(lee$Z)
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
    gene_a = genes[idx[, 1]],
    gene_b = genes[idx[, 2]],
    weight = weight,
    p_value = p_value,
    fdr = fdr,
    n_obs = as.integer(n_obs),
    edge_type = edge_type
  )
  data.table::fwrite(edges, out_path, sep = "\t")
  list(edges_n_total = nrow(edges), edges_n_obs = as.integer(n_obs), edge_type = edge_type)
}

lee <- extract_lee_stats(scope_obj, grid_name)
gg_out <- emit_gg_matrices(lee, rownames(scope_obj@meta.data), outdir)

n_obs_grid <- infer_grid_n_obs(scope_obj, grid_name)
edges_out <- emit_edges_all(lee, gg_out$gg_genes, n_obs_grid, outdir, edge_type = "lee_L")

method_version <- tryCatch(as.character(utils::packageVersion("geneSCOPE")), error = function(e) "unknown")
read_single_line <- function(path) {
  if (!file.exists(path)) return("")
  x <- readLines(path, n = 1L, warn = FALSE)
  if (length(x)) trimws(x[[1L]]) else ""
}
package_source_commit <- read_single_line("/tmp/geneSCOPE_provenance/source_commit")
package_vendor_tree_sha256 <- read_single_line("/tmp/geneSCOPE_provenance/vendor_tree_sha256")
conda_explicit_sha256 <- sha256_file("/tmp/geneSCOPE_provenance/conda-explicit.txt")
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
  package_source_commit = package_source_commit,
  package_vendor_tree_sha256 = package_vendor_tree_sha256,
  conda_explicit_sha256 = conda_explicit_sha256,
  run_id = run_id,
  timestamp = timestamp,
  gg_weight_type = "Lee_L",
  gg_fdr_source = "native_api",
  gg_dim = list(n = gg_out$gg_genes_n),
  gg_genes_n = gg_out$gg_genes_n,
  edges_file = "all_edges.tsv",
  edges_schema = c("gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"),
  edges_weight_type = "Lee_L",
  edges_stats_layer = gg_out$gg_layer_name,
  edges_p_source = if (!is.null(lee$meta$p_source)) lee$meta$p_source else "unknown",
  edges_p_resolution = if (!is.null(lee$meta$p_resolution)) lee$meta$p_resolution else "",
  edges_fdr_method = if (!is.null(lee$meta$FDR_main_method)) lee$meta$FDR_main_method else "unknown",
  edges_edge_type = edges_out$edge_type,
  edges_n_total = edges_out$edges_n_total,
  edges_n_obs = edges_out$edges_n_obs,
  input_dataset_id = dataset_id,
  input_sha256 = input_sha256,
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
  gene_filtering = list(pct_min = cluster_pct),
  random_seed = seed,
  stochastic = TRUE,
  params = list(
    grid_um = grid_um,
    ncores = ncores,
    parallel_backend = parallel_backend,
    membership_col = membership_col,
    lee_formula = "n/S2 * crossprod(W %*% Xz) / sqrt(crossprod norms)",
    lee_formula_id = if (!is.null(lee$meta$formula_id)) lee$meta$formula_id else "",
    permutation_scheme = "global_joint_shuffle",
    use_blocks = FALSE,
    permutations = perms,
    rng_kind = "L'Ecuyer-CMRG",
    seed = seed,
    cluster_pct = cluster_pct,
    cluster_resolution = 0.1,
    cluster_consensus_threshold = 0.95,
    cluster_n_restart = n_restart,
    weight_provenance = scope_obj@grid[[grid_name]]$weight_provenance,
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

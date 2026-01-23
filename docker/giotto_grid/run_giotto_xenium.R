#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(getopt)
  library(Giotto)
})

options(giotto.update_param = FALSE)

start_time <- Sys.time()

spec <- matrix(
  c(
    "data_dir", "d", 1, "character",
    "outdir", "o", 1, "character",
    "grid_stepsize", "Z", 1, "double",
    "spatial_grid_name", "Y", 1, "character",
    "min_cells_per_grid", "P", 1, "integer",
    "n_corr_genes", "G", 1, "integer",
    "k_modules", "k", 1, "integer",
    "expression_values", "V", 1, "character",
    "seed", "s", 1, "integer",
    "coord_file", "c", 1, "character",
    "roi_flip_y", "F", 0, "logical",

    # legacy flags (accepted but ignored)
    "max_cells", "m", 1, "integer",
    "repeat", "R", 1, "integer",
    "dataset_id", "D", 1, "character",
    "roi_id", "r", 1, "character",
    "threads", "t", 1, "integer",
    "ncores", "n", 1, "integer"
  ),
  byrow = TRUE,
  ncol = 4
)

opt <- getopt(spec)

if (is.null(opt$data_dir) || is.null(opt$outdir) || is.null(opt$grid_stepsize)) {
  cat(
    paste0(
      "Usage:\n",
      "  run_giotto_xenium.R --data_dir <xenium_outs> --outdir <outdir> --grid_stepsize <num>\n",
      "    [--expression_values normalized|scaled|custom]\n",
      "    [--n_corr_genes 500] [--k_modules 5] [--min_cells_per_grid 4]\n",
      "    [--spatial_grid_name spatial_grid] [--seed 1]\n",
      "    [--coord_file roi.csv --roi_flip_y]\n"
    ),
    file = stderr()
  )
  quit(status = 2)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

grid_stepsize <- suppressWarnings(as.numeric(opt$grid_stepsize))
if (!is.finite(grid_stepsize) || grid_stepsize <= 0) stop("--grid_stepsize must be > 0")

spatial_grid_name <- if (!is.null(opt$spatial_grid_name) && nzchar(opt$spatial_grid_name)) {
  as.character(opt$spatial_grid_name)
} else {
  "spatial_grid"
}

min_cells_per_grid <- if (!is.null(opt$min_cells_per_grid)) suppressWarnings(as.integer(opt$min_cells_per_grid)) else 4L
if (!is.finite(min_cells_per_grid) || min_cells_per_grid < 1L) min_cells_per_grid <- 4L

n_corr_genes <- if (!is.null(opt$n_corr_genes)) suppressWarnings(as.integer(opt$n_corr_genes)) else 500L
if (!is.finite(n_corr_genes) || n_corr_genes < 2L) n_corr_genes <- 500L

k_modules <- if (!is.null(opt$k_modules)) suppressWarnings(as.integer(opt$k_modules)) else 5L
if (!is.finite(k_modules) || k_modules < 1L) k_modules <- 5L

expression_values <- if (!is.null(opt$expression_values) && nzchar(opt$expression_values)) {
  tolower(as.character(opt$expression_values))
} else {
  "normalized"
}
if (!expression_values %in% c("normalized", "scaled", "custom")) {
  stop("--expression_values must be one of: normalized|scaled|custom")
}

if (!is.null(opt$seed)) {
  seed <- suppressWarnings(as.integer(opt$seed))
  if (is.finite(seed)) set.seed(seed)
}

cat(
  sprintf(
    "[INFO] data_dir=%s outdir=%s grid_stepsize=%s spatial_grid_name=%s expression_values=%s n_corr_genes=%d k_modules=%d min_cells_per_grid=%d\n",
    as.character(opt$data_dir),
    as.character(opt$outdir),
    as.character(grid_stepsize),
    as.character(spatial_grid_name),
    as.character(expression_values),
    as.integer(n_corr_genes),
    as.integer(k_modules),
    as.integer(min_cells_per_grid)
  ),
  file = stderr()
)
flush(stderr())

# ---------------------------
# 1) createGiottoXeniumObject (10x aggregated expression)
# ---------------------------
cat("[STEP] createGiottoXeniumObject\n", file = stderr()); flush(stderr())

if (!is.character(opt$data_dir) || length(opt$data_dir) != 1L || !nzchar(opt$data_dir)) {
  stop("--data_dir must be a single non-empty path string")
}
if (!dir.exists(opt$data_dir)) {
  stop("--data_dir does not exist or is not a directory: ", opt$data_dir)
}

expr_candidates <- c(
  file.path(opt$data_dir, "cell_feature_matrix.h5"),
  file.path(opt$data_dir, "cell_feature_matrix", "cell_feature_matrix.h5")
)
expr_h5_src <- expr_candidates[file.exists(expr_candidates)][1]
if (is.na(expr_h5_src) || !nzchar(expr_h5_src)) {
  found <- list.files(opt$data_dir, pattern = "cell_feature_matrix.*\\.h5$", recursive = TRUE, full.names = TRUE)
  if (length(found) > 0) {
    expr_h5_src <- found[1]
  } else {
    listing <- paste(utils::head(list.files(opt$data_dir), 50), collapse = ", ")
    stop(
      "Could not find Xenium expression H5 under --data_dir. Expected 'cell_feature_matrix.h5'. ",
      "Top-level entries: ", listing
    )
  }
}

expr_h5_use <- tempfile(pattern = "cell_feature_matrix_", fileext = ".h5")
if (file.exists(expr_h5_src)) {
  expr_h5_use <- tempfile(pattern = "cell_feature_matrix_", fileext = ".h5")
  ok_copy <- try(file.copy(expr_h5_src, expr_h5_use, overwrite = TRUE), silent = TRUE)
  if (inherits(ok_copy, "try-error") || !isTRUE(ok_copy)) stop("Failed to copy expression H5: ", expr_h5_src)
  cat(sprintf("[INFO] expression_h5_tmp=%s\n", expr_h5_use), file = stderr()); flush(stderr())
} else {
  stop("Expression H5 path resolved but does not exist: ", expr_h5_src)
}

gobject <- Giotto::createGiottoXeniumObject(
  xenium_dir = opt$data_dir,
  expression_path = as.character(expr_h5_use),
  load_transcripts = FALSE,
  load_expression = TRUE,
  load_cellmeta = TRUE,
  load_images = NULL,
  load_aligned_images = NULL,
  bounds_path = list(cell = "cell"),
  verbose = FALSE
)

cat("[STEP] createGiottoXeniumObject done\n", file = stderr()); flush(stderr())

# ---------------------------
# 2) ROI subset (optional)
# ---------------------------
coord_file <- if (!is.null(opt$coord_file) && nzchar(opt$coord_file)) as.character(opt$coord_file) else ""
if (nzchar(coord_file)) {
  cat(sprintf("[STEP] ROI subset coord_file=%s flip_y=%s\n", coord_file, ifelse(isTRUE(opt$roi_flip_y), "TRUE", "FALSE")),
      file = stderr()); flush(stderr())

  if (!file.exists(coord_file)) stop("Missing --coord_file: ", coord_file)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required for ROI selection.")

  roi_dt <- data.table::fread(coord_file)
  names(roi_dt) <- tolower(names(roi_dt))
  if (!all(c("x", "y") %in% names(roi_dt))) stop("--coord_file must have X and Y columns (case-insensitive).")

  pts <- as.matrix(roi_dt[, c("x", "y"), with = FALSE])
  if (!is.matrix(pts) || nrow(pts) < 3) stop("ROI polygon requires >= 3 vertices.")
  pts[, 1] <- suppressWarnings(as.numeric(pts[, 1]))
  pts[, 2] <- suppressWarnings(as.numeric(pts[, 2]))
  if (any(!is.finite(pts))) stop("ROI polygon has non-finite coordinates.")
  if (isTRUE(opt$roi_flip_y)) pts[, 2] <- -pts[, 2]
  if (sum(abs(pts[1, ] - pts[nrow(pts), ])) > 1e-9) pts <- rbind(pts, pts[1, ])

  roi_df <- data.frame(
    poly_ID = rep("roi", nrow(pts)),
    x = pts[, 1],
    y = pts[, 2],
    stringsAsFactors = FALSE
  )

  gpoly <- GiottoClass::createGiottoPolygon(roi_df, name = "roi", calc_centroids = FALSE, verbose = FALSE)
  gobject <- GiottoClass::addGiottoPolygons(gobject, gpolygons = list(gpoly))

  cells_sv <- Giotto::getCellsFromPolygon(gobject, polygon_name = "roi", spat_unit = "cell", spat_loc_name = "raw")
  cells_df <- terra::as.data.frame(cells_sv)
  if (is.null(cells_df) || !is.data.frame(cells_df) || nrow(cells_df) == 0) stop("ROI polygon selected 0 cells.")

  cn <- tolower(names(cells_df))
  col_i <- NA_integer_
  if (any(cn == "cell_id")) col_i <- which(cn == "cell_id")[1]
  else if (any(grepl("cell_id", cn))) col_i <- which(grepl("cell_id", cn))[1]
  else if (any(cn == "cell")) col_i <- which(cn == "cell")[1]
  else if (any(grepl("cell", cn))) col_i <- which(grepl("cell", cn))[1]
  if (!is.finite(col_i)) stop("Cannot infer cell id column from ROI result. Columns: ", paste(names(cells_df), collapse = ", "))

  cell_ids <- unique(as.character(cells_df[[col_i]]))
  cell_ids <- cell_ids[nzchar(cell_ids)]
  if (length(cell_ids) == 0) stop("ROI polygon selected 0 cells.")

  cat(sprintf("[INFO] ROI selected cells: %d\n", length(cell_ids)), file = stderr()); flush(stderr())
  gobject <- Giotto::subsetGiotto(gobject, cell_ids = cell_ids, verbose = FALSE)
  cat("[STEP] ROI subset done\n", file = stderr()); flush(stderr())
}

# ---------------------------
# 3) processGiotto (no filtering)
# ---------------------------
cat("[STEP] processGiotto\n", file = stderr()); flush(stderr())
gobject <- Giotto::processGiotto(
  gobject = gobject,
  filter_params = list(
    expression_threshold   = 0,
    feat_det_in_min_cells  = 0,
    min_det_feats_per_cell = 0
  ),
  norm_params = list(),
  stat_params = list(),
  adjust_params = NULL,
  verbose = FALSE
)
cat("[STEP] processGiotto done\n", file = stderr()); flush(stderr())

# Fix defaults explicitly (avoid multi-value default resolution)
GiottoClass::activeSpatUnit(gobject) <- "cell"
GiottoClass::activeFeatType(gobject) <- "rna"

# ---------------------------
# 4) createSpatialGrid (explicit default grid creator)
# ---------------------------
cat("[STEP] createSpatialGrid\n", file = stderr()); flush(stderr())
gobject <- GiottoClass::createSpatialDefaultGrid(
  gobject         = gobject,
  spat_unit       = "cell",
  feat_type       = "rna",
  spat_loc_name   = "raw",
  sdimx_stepsize  = grid_stepsize,
  sdimy_stepsize  = grid_stepsize,
  minimum_padding = 1,
  name            = spatial_grid_name,
  return_gobject  = TRUE
)
cat("[STEP] createSpatialGrid done\n", file = stderr()); flush(stderr())

# ---------------------------
# 5) getExpression + gene selection (mean)
# ---------------------------
if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")
cat("[STEP] getExpression + gene selection\n", file = stderr()); flush(stderr())

expr_mat <- Giotto::getExpression(
  gobject = gobject,
  values = expression_values,
  spat_unit = "cell",
  feat_type = "rna",
  output = "matrix"
)
if (is.null(expr_mat) || is.null(rownames(expr_mat)) || is.null(colnames(expr_mat))) stop("getExpression() returned invalid matrix.")

scores <- if (requireNamespace("Matrix", quietly = TRUE) && inherits(expr_mat, "Matrix")) {
  Matrix::rowMeans(expr_mat)
} else {
  base::rowMeans(expr_mat)
}

genes <- as.character(rownames(expr_mat))
scores <- suppressWarnings(as.numeric(scores))
gene_dt <- data.table::data.table(gene = genes, score = scores)
gene_dt <- gene_dt[nzchar(gene) & is.finite(score)]
if (nrow(gene_dt) < 2) stop("Too few genes available after scoring.")
data.table::setorder(gene_dt, -score, gene)

n_take <- min(n_corr_genes, nrow(gene_dt))
selected_genes <- unique(gene_dt$gene[seq_len(n_take)])
selected_genes <- selected_genes[nzchar(selected_genes)]
if (length(selected_genes) < 2) stop("Selected < 2 genes after filtering.")
if (k_modules > length(selected_genes)) stop("--k_modules must be <= selected_genes count.")
cat(sprintf("[INFO] selected_genes=%d\n", length(selected_genes)), file = stderr()); flush(stderr())

expr_mat_sub <- expr_mat[selected_genes, , drop = FALSE]

# ---------------------------
# 6) detectSpatialCorFeatsMatrix (grid-averaging)
#    Use spatial_locs as data.table, spatial_grid as gridDT with compat class tag.
# ---------------------------
cat("[STEP] detectSpatialCorFeatsMatrix (grid)\n", file = stderr()); flush(stderr())

spat_locs_dt <- GiottoClass::getSpatialLocations(
  gobject      = gobject,
  spat_unit    = "cell",
  name         = "raw",
  output       = "data.table",
  copy_obj     = TRUE,
  verbose      = FALSE,
  set_defaults = TRUE
)

grid_dt <- GiottoClass::getSpatialGrid(
  gobject         = gobject,
  spat_unit       = "cell",
  feat_type       = "rna",
  name            = spatial_grid_name,
  return_grid_Obj = FALSE,
  set_defaults    = TRUE
)
if (is.null(grid_dt) || !is.data.frame(grid_dt)) stop("getSpatialGrid(return_grid_Obj=FALSE) failed.")

# Minimal compat tag to satisfy internal type checks while keeping $ column access
class(grid_dt) <- unique(c("spatialGridObj", class(grid_dt)))

spat_cor <- Giotto::detectSpatialCorFeatsMatrix(
  expression_matrix   = expr_mat_sub,
  method              = "grid",
  spatial_network     = NULL,
  spatial_grid        = grid_dt,
  spatial_locs        = spat_locs_dt,
  subset_feats        = NULL,
  min_cells_per_grid  = min_cells_per_grid,
  cor_method          = "pearson"
)

cat("[STEP] detectSpatialCorFeatsMatrix done\n", file = stderr()); flush(stderr())

# ---------------------------
# 7) clusterSpatialCorFeats
# ---------------------------
cat("[STEP] clusterSpatialCorFeats\n", file = stderr()); flush(stderr())

spat_cor <- Giotto::clusterSpatialCorFeats(
  spatCorObject = spat_cor,
  name = "spat_grid_clus",
  hclust_method = "ward.D",
  k = k_modules,
  return_obj = TRUE
)

cat("[STEP] clusterSpatialCorFeats done\n", file = stderr()); flush(stderr())

# ---------------------------
# 8) outputs: modules.tsv + all_edges.tsv only
# ---------------------------
clusters <- NULL
if (!is.null(spat_cor[["cor_clusters"]])) {
  clusters <- spat_cor[["cor_clusters"]][["spat_grid_clus"]]
}
if (is.null(clusters) || length(clusters) == 0 || is.null(names(clusters))) {
  stop("No clusters found in spatCorObject for name 'spat_grid_clus'.")
}

modules_dt <- data.table::data.table(
  gene = as.character(names(clusters)),
  module_id = suppressWarnings(as.integer(clusters))
)
modules_dt <- modules_dt[nzchar(gene)]
data.table::setorder(modules_dt, module_id, gene)

cor_dt <- data.table::as.data.table(spat_cor[["cor_DT"]])
if (is.null(cor_dt) || !all(c("feat_ID", "variable", "spat_cor") %in% names(cor_dt))) {
  stop("spatCorObject is missing expected cor_DT columns: feat_ID, variable, spat_cor")
}

edges_dt <- cor_dt[, list(
  gene_a = as.character(feat_ID),
  gene_b = as.character(variable),
  weight = suppressWarnings(as.numeric(spat_cor))
)]
edges_dt <- edges_dt[nzchar(gene_a) & nzchar(gene_b) & gene_a != gene_b]
edges_dt[, `:=`(ga = pmin(gene_a, gene_b), gb = pmax(gene_a, gene_b))]
edges_dt <- edges_dt[ga < gb]
edges_dt <- edges_dt[, list(
  weight = if (all(is.na(weight))) NA_real_ else mean(weight, na.rm = TRUE)
), by = list(gene_a = ga, gene_b = gb)]
data.table::setorder(edges_dt, gene_a, gene_b)

cat(sprintf("[INFO] edges=%d modules=%d\n", nrow(edges_dt), nrow(modules_dt)), file = stderr()); flush(stderr())

data.table::fwrite(edges_dt, file.path(opt$outdir, "all_edges.tsv"), sep = "\t")
data.table::fwrite(modules_dt, file.path(opt$outdir, "modules.tsv"), sep = "\t")

wall_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("[DONE] wall_sec=%.2f\n", wall_sec), file = stderr()); flush(stderr())

#!/usr/bin/env Rscript
# Giotto grid-only runner.
# References:
# - createSpatialGrid: https://giottosuite.com/reference/createSpatialGrid.html
# - detectSpatialCorFeats: https://giottosuite.com/reference/detectSpatialCorFeats.html
# - clusterSpatialCorFeats: https://giottosuite.com/reference/clusterSpatialCorFeats.html
#
# This script intentionally:
# - always runs the grid workflow (no network branch)
# - never computes all-genes all-to-all edges
# - leaves correlation gene count to user (HPC)

suppressPackageStartupMessages({
  library(getopt)
  library(data.table)
  library(Matrix)
  library(jsonlite)
  library(Giotto)
})

options(giotto.update_param = FALSE)

N_CORR_GENES_DEFAULT <- 500L

eprint <- function(...) {
  cat(sprintf(...), "\n", file = stderr())
  flush(stderr())
}

step_start <- function(label) {
  eprint("[STEP] %s start", label)
  Sys.time()
}

step_end <- function(label, t0) {
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  eprint("[STEP] %s done (wall_sec=%.2f)", label, wall)
}

as_int <- function(x, default = 0L) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) return(default)
  v
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) return(default)
  v
}

call_giotto_name <- function(fn_name, args) {
  fn <- get(fn_name, mode = "function")
  fml <- tryCatch(names(formals(fn)), error = function(e) NULL)
  if (is.null(fml) || "..." %in% fml) {
    return(do.call(fn_name, args))
  }
  args2 <- args[names(args) %in% fml]
  do.call(fn_name, args2)
}

find_col <- function(df, candidates) {
  lc <- setNames(names(df), tolower(names(df)))
  for (cand in candidates) {
    hit <- lc[[tolower(cand)]]
    if (!is.null(hit)) return(hit)
  }
  NULL
}

read_roi_polygon <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) > 0) lines[1] <- sub("^\\ufeff", "", lines[1])

  lines2 <- lines[!grepl("^\\s*#", lines)]
  lines2 <- lines2[nzchar(trimws(lines2))]
  if (length(lines2) < 4) stop("ROI CSV after dropping comment lines is too short: ", path)

  dt <- data.table::fread(text = paste(lines2, collapse = "\n"))
  cn <- tolower(names(dt))
  xcol <- which(cn %in% c("x"))
  ycol <- which(cn %in% c("y"))
  if (length(xcol) != 1 || length(ycol) != 1) {
    stop(sprintf("ROI CSV must have columns X,Y (case-insensitive). cols=%s",
                 paste(names(dt), collapse = ",")))
  }

  xx <- suppressWarnings(as.numeric(dt[[xcol]]))
  yy <- suppressWarnings(as.numeric(dt[[ycol]]))
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

bbox_from_xy <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NULL)
  list(
    xmin = min(x),
    xmax = max(x),
    ymin = min(y),
    ymax = max(y)
  )
}

format_bbox <- function(bbox) {
  if (is.null(bbox)) return("NA")
  sprintf("x=[%.3f, %.3f] y=[%.3f, %.3f]", bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax)
}

select_roi_cells <- function(gobject, roi_vertices, name) {
  roi_df <- data.frame(
    poly_ID = rep(name, nrow(roi_vertices)),
    x = roi_vertices[, 1],
    y = roi_vertices[, 2],
    stringsAsFactors = FALSE
  )
  gpoly <- createGiottoPolygon(roi_df, name = name, calc_centroids = FALSE, verbose = FALSE)
  gobject <- addGiottoPolygons(gobject, gpolygons = list(gpoly))
  polygon_cells <- getCellsFromPolygon(gobject, polygon_name = name, spat_unit = "cell")
  poly_df <- terra::as.data.frame(polygon_cells)
  poly_cell_col <- find_col(poly_df, c("cell_ID", "cell_id", "cell", "barcode", "obs_id"))
  if (is.null(poly_cell_col)) stop("Cannot find cell id column from ROI polygon selection.")
  cell_ids <- unique(as.character(poly_df[[poly_cell_col]]))
  cell_ids <- cell_ids[nzchar(cell_ids)]
  list(gobject = gobject, cell_ids = cell_ids)
}

sha256_file <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return("")
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

resolve_xenium_expression_h5 <- function(xenium_dir) {
  candidates <- c(
    file.path(xenium_dir, "cell_feature_matrix.h5"),
    file.path(xenium_dir, "cell_feature_matrix", "cell_feature_matrix.h5"),
    file.path(xenium_dir, "filtered_feature_bc_matrix.h5"),
    file.path(xenium_dir, "filtered_feature_bc_matrix", "filtered_feature_bc_matrix.h5")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

copy_h5_if_readonly <- function(h5_path) {
  if (is.null(h5_path) || !nzchar(h5_path) || !file.exists(h5_path)) return(h5_path)
  writable <- suppressWarnings(file.access(h5_path, 2) == 0)
  if (isTRUE(writable)) return(h5_path)

  tmp_dir <- file.path(tempdir(), "xenium_shadow")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  tmp_path <- file.path(tmp_dir, basename(h5_path))

  ok <- tryCatch(file.copy(h5_path, tmp_path, overwrite = TRUE), error = function(e) FALSE)
  if (!isTRUE(ok) || !file.exists(tmp_path)) {
    stop("Expression H5 is not writable and copy failed. h5_path=", h5_path)
  }
  tmp_path
}

compute_grid_n_obs <- function(gobject,
                               spat_unit,
                               feat_type,
                               spat_loc_name,
                               spatial_grid_name,
                               min_cells_per_grid) {
  grid_obj <- getSpatialGrid(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type,
    name = spatial_grid_name,
    return_grid_Obj = FALSE,
    set_defaults = TRUE
  )

  spat_locs <- getSpatialLocations(
    gobject = gobject,
    spat_unit = spat_unit,
    name = spat_loc_name,
    output = "data.table",
    copy_obj = TRUE,
    verbose = FALSE
  )
  cell_col <- find_col(spat_locs, c("cell_ID", "cell_id", "cell", "barcode", "obs_id"))
  x_col <- find_col(spat_locs, c("sdimx", "x", "x_centroid", "x_coord"))
  y_col <- find_col(spat_locs, c("sdimy", "y", "y_centroid", "y_coord"))
  if (is.null(cell_col) || is.null(x_col) || is.null(y_col)) {
    stop("Cannot resolve cell_ID/sdimx/sdimy columns for grid n_obs computation.")
  }

  spat_locs2 <- data.table::data.table(
    cell_ID = as.character(spat_locs[[cell_col]]),
    sdimx = as.numeric(spat_locs[[x_col]]),
    sdimy = as.numeric(spat_locs[[y_col]])
  )
  spat_locs2 <- spat_locs2[nzchar(cell_ID) & is.finite(sdimx) & is.finite(sdimy)]

  spat_locs_annot <- GiottoClass:::annotate_spatlocs_with_spatgrid_2D(
    spatloc = spat_locs2,
    spatgrid = grid_obj
  )
  if (!"gr_loc" %in% names(spat_locs_annot)) stop("annotate_spatlocs_with_spatgrid_2D missing gr_loc.")

  gr <- as.character(spat_locs_annot$gr_loc)
  gr <- gr[nzchar(gr)]
  if (!length(gr)) return(0L)

  min_cells_per_grid <- suppressWarnings(as.integer(min_cells_per_grid))
  if (!is.finite(min_cells_per_grid) || min_cells_per_grid < 1L) min_cells_per_grid <- 1L

  cells_per_grid <- table(gr)
  sum(cells_per_grid >= min_cells_per_grid)
}

compute_auto_grid_stepsize <- function(gobject, spat_unit, spat_loc_name, target_bins) {
  spat_locs <- getSpatialLocations(
    gobject = gobject,
    spat_unit = spat_unit,
    name = spat_loc_name,
    output = "data.table",
    copy_obj = TRUE,
    verbose = FALSE
  )
  cell_col <- find_col(spat_locs, c("cell_ID", "cell_id", "cell", "barcode", "obs_id"))
  x_col <- find_col(spat_locs, c("sdimx", "x", "x_centroid", "x_coord"))
  y_col <- find_col(spat_locs, c("sdimy", "y", "y_centroid", "y_coord"))
  if (is.null(cell_col) || is.null(x_col) || is.null(y_col)) {
    stop("Cannot resolve cell_ID/sdimx/sdimy columns for auto grid stepsize computation.")
  }

  bb <- bbox_from_xy(spat_locs[[x_col]], spat_locs[[y_col]])
  if (is.null(bb)) stop("Cannot compute bbox for auto grid stepsize computation.")

  dx <- suppressWarnings(as.numeric(bb$xmax - bb$xmin))
  dy <- suppressWarnings(as.numeric(bb$ymax - bb$ymin))
  if (!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0) {
    stop("Invalid bbox for auto grid stepsize computation.")
  }
  area <- dx * dy

  target_bins <- suppressWarnings(as.integer(target_bins))
  if (!is.finite(target_bins) || target_bins < 1L) target_bins <- 1L

  step <- sqrt(area / target_bins)
  if (!is.finite(step) || step <= 0) stop("Auto grid stepsize computation failed.")
  step
}

do_spatial_grid_averaging_gridDT <- function(expression_matrix,
                                             spatial_grid_dt,
                                             spatial_locs,
                                             min_cells_per_grid = 4) {
  expr_values <- expression_matrix

  if (!is.data.frame(spatial_grid_dt)) {
    stop("do_spatial_grid_averaging_gridDT: spatial_grid_dt must be a data.frame/data.table.")
  }
  needed_cols <- c("x_start", "x_end", "y_start", "y_end", "gr_name", "gr_x_name", "gr_y_name")
  if (!all(needed_cols %in% colnames(spatial_grid_dt))) {
    stop("do_spatial_grid_averaging_gridDT: spatial_grid_dt missing required columns.")
  }

  Giotto:::evaluate_provided_spatial_locations(spatial_locs)

  if (all(c("sdimx", "sdimy", "sdimz") %in% colnames(spatial_locs))) {
    spatial_locs <- GiottoClass:::annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid_dt)
  } else if (all(c("sdimx", "sdimy") %in% colnames(spatial_locs))) {
    spatial_locs <- GiottoClass:::annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid_dt)
  } else {
    stop("do_spatial_grid_averaging_gridDT: spatial_locs missing sdimx/sdimy columns.")
  }

  gr_loc <- NULL
  min_cells_per_grid <- suppressWarnings(as.integer(min_cells_per_grid))
  if (!is.finite(min_cells_per_grid) || min_cells_per_grid < 1L) min_cells_per_grid <- 1L

  cells_per_grid <- sort(table(spatial_locs$gr_loc))
  cells_per_grid <- cells_per_grid[cells_per_grid >= min_cells_per_grid]
  loc_names <- names(cells_per_grid)
  if (length(loc_names) == 0) {
    stop("do_spatial_grid_averaging_gridDT: no grid bins pass min_cells_per_grid=", min_cells_per_grid)
  }

  loc_av_expr_list <- vector("list", length(loc_names))
  names(loc_av_expr_list) <- loc_names
  for (loc_name in loc_names) {
    loc_cell_IDs <- spatial_locs[gr_loc == loc_name]$cell_ID
    subset_expr <- expr_values[, colnames(expr_values) %in% loc_cell_IDs]
    if (is.vector(subset_expr) == TRUE) {
      loc_av_expr <- subset_expr
    } else {
      loc_av_expr <- rowMeans(subset_expr)
    }
    loc_av_expr_list[[loc_name]] <- loc_av_expr
  }

  loc_av_expr_matrix <- do.call("cbind", loc_av_expr_list)
  as.matrix(loc_av_expr_matrix)
}

detectSpatialCorFeats_grid_compat <- function(gobject,
                                             spat_unit = NULL,
                                             feat_type = NULL,
                                             spat_loc_name = "raw",
                                             expression_values = c("raw", "counts"),
                                             subset_feats = NULL,
                                             spatial_grid_name = "spatial_grid",
                                             min_cells_per_grid = 4,
                                             cor_method = c("pearson", "kendall", "spearman")) {
  # Giotto 4.2.2 grid workflow bug:
  # - detectSpatialCorFeats(method='grid') passes a gridDT (return_grid_Obj=FALSE)
  # - but do_spatial_grid_averaging() incorrectly requires spatialGridObj and stops.
  # This compat re-implements the documented computation using the gridDT.

  spat_unit <- set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
  feat_type <- set_default_feat_type(gobject = gobject, spat_unit = spat_unit, feat_type = feat_type)

  cor_method <- match.arg(cor_method, choices = c("pearson", "kendall", "spearman"))
  values <- match.arg(as.character(expression_values)[[1]], unique(c("raw", "counts", expression_values)))

  sanitize_cor_matrix <- function(mat, label = "cor") {
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    if (!is.numeric(mat)) storage.mode(mat) <- "double"
    bad <- !is.finite(mat)
    n_bad <- sum(bad)
    if (n_bad > 0) {
      eprint("[WARN] %s has %d non-finite entries (NA/NaN/Inf); setting to 0", label, n_bad)
      mat[bad] <- 0
    }
    if (!is.null(dim(mat)) && nrow(mat) == ncol(mat)) {
      diag(mat) <- 1
    }
    mat[mat > 1] <- 1
    mat[mat < -1] <- -1
    mat
  }

  cor_rows_pearson_sparse <- function(m) {
    m <- Matrix::Matrix(m, sparse = TRUE)
    genes <- rownames(m)
    if (is.null(genes) || !length(genes)) {
      genes <- as.character(seq_len(nrow(m)))
    } else {
      genes <- as.character(genes)
    }
    n <- as.integer(ncol(m))
    if (!is.finite(n) || n < 2L) stop("Need >=2 observations to compute correlation.")

    s <- Matrix::rowSums(m)
    s2 <- Matrix::rowSums(m ^ 2)
    s <- as.numeric(s)
    s2 <- as.numeric(s2)

    ss <- s2 - (s * s) / n
    ss[!is.finite(ss)] <- 0
    ss[ss < 0] <- 0
    denom <- sqrt(ss)
    denom[denom == 0] <- NA_real_

    xtx <- Matrix::tcrossprod(m)
    xtx <- as.matrix(xtx)
    numer <- xtx - tcrossprod(s) / n

    cor_mat <- numer / denom
    cor_mat <- t(t(cor_mat) / denom)
    cor_mat[!is.finite(cor_mat)] <- 0
    diag(cor_mat) <- 1
    cor_mat[cor_mat > 1] <- 1
    cor_mat[cor_mat < -1] <- -1
    rownames(cor_mat) <- genes
    colnames(cor_mat) <- genes
    cor_mat
  }

  expr_values <- getExpression(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type,
    values = values,
    output = "matrix"
  )
  if (!is.null(subset_feats)) {
    expr_values <- expr_values[rownames(expr_values) %in% subset_feats, , drop = FALSE]
  }

  spatial_locs <- getSpatialLocations(
    gobject = gobject,
    spat_unit = spat_unit,
    name = spat_loc_name,
    output = "data.table",
    copy_obj = TRUE
  )

  spatial_grid_dt <- getSpatialGrid(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type,
    name = spatial_grid_name,
    return_grid_Obj = FALSE
  )

  loc_av_expr_matrix <- do_spatial_grid_averaging_gridDT(
    expression_matrix = expr_values,
    spatial_grid_dt = spatial_grid_dt,
    spatial_locs = spatial_locs,
    min_cells_per_grid = min_cells_per_grid
  )

  cor_spat_matrix <- suppressWarnings(stats::cor(t(loc_av_expr_matrix), method = cor_method))
  cor_spat_matrix <- sanitize_cor_matrix(cor_spat_matrix, "cor_spat_matrix")
  cor_spat_matrixDT <- data.table::as.data.table(cor_spat_matrix)
  cor_spat_matrixDT[, `:=`(feat_ID, rownames(cor_spat_matrix))]
  cor_spat_DT <- data.table::melt.data.table(data = cor_spat_matrixDT, id.vars = "feat_ID", value.name = "spat_cor")

  cordiff <- spat_cor <- expr_cor <- spatrank <- exprrank <- rankdiff <- NULL
  if (cor_method == "pearson" && inherits(expr_values, "Matrix")) {
    cor_matrix <- cor_rows_pearson_sparse(expr_values)
  } else {
    cor_matrix <- suppressWarnings(stats::cor(t(as.matrix(expr_values)), method = cor_method))
  }
  cor_matrix <- sanitize_cor_matrix(cor_matrix, "cor_matrix")
  cor_matrixDT <- data.table::as.data.table(cor_matrix)
  cor_matrixDT[, `:=`(feat_ID, rownames(cor_matrix))]
  cor_DT <- data.table::melt.data.table(data = cor_matrixDT, id.vars = "feat_ID", value.name = "expr_cor")

  data.table::setorder(cor_spat_DT, feat_ID, variable)
  data.table::setorder(cor_DT, feat_ID, variable)
  doubleDT <- cbind(cor_spat_DT, expr_cor = cor_DT[["expr_cor"]])
  doubleDT[, `:=`(cordiff, spat_cor - expr_cor)]
  doubleDT[, `:=`(spatrank, data.table::frank(-spat_cor, ties.method = "first")), by = feat_ID]
  doubleDT[, `:=`(exprrank, data.table::frank(-expr_cor, ties.method = "first")), by = feat_ID]
  doubleDT[, `:=`(rankdiff, spatrank - exprrank)]
  data.table::setorder(doubleDT, feat_ID, -spat_cor)

  spatCorObject <- list(
    cor_DT = doubleDT,
    feat_order = rownames(cor_spat_matrix),
    cor_hclust = list(),
    cor_clusters = list()
  )
  class(spatCorObject) <- append("spatCorObject", class(spatCorObject))
  spatCorObject
}

extract_complete_edges <- function(spat_cor_obj, gene_set) {
  cor_dt <- spat_cor_obj[["cor_DT"]]
  if (is.null(cor_dt)) stop("spatCorObject missing cor_DT")
  cor_dt <- data.table::as.data.table(cor_dt)
  if (!all(c("feat_ID", "variable", "spat_cor") %in% names(cor_dt))) {
    stop("spatCorObject cor_DT missing feat_ID/variable/spat_cor columns")
  }

  genes <- as.character(gene_set)
  genes <- genes[nzchar(genes)]
  genes <- sort(unique(genes))
  if (length(genes) < 2) {
    return(data.table(gene_a = character(), gene_b = character(), weight = numeric()))
  }

  weights_dt <- cor_dt[, .(
    gene_a = as.character(feat_ID),
    gene_b = as.character(variable),
    weight = suppressWarnings(as.numeric(spat_cor))
  )]
  weights_dt <- weights_dt[gene_a %in% genes & gene_b %in% genes]
  weights_dt <- weights_dt[gene_a != gene_b]
  weights_dt[, `:=`(ga = pmin(gene_a, gene_b), gb = pmax(gene_a, gene_b))]
  weights_dt <- weights_dt[ga < gb]
  weights_dt <- weights_dt[, .(
    weight = if (all(is.na(weight))) NA_real_ else mean(weight, na.rm = TRUE)
  ), by = .(gene_a = ga, gene_b = gb)]

  idx <- data.table::CJ(i = seq_along(genes), j = seq_along(genes))
  idx <- idx[i < j]
  edges <- data.table::data.table(
    gene_a = genes[idx$i],
    gene_b = genes[idx$j]
  )
  data.table::setkey(weights_dt, gene_a, gene_b)
  edges[weights_dt, weight := i.weight, on = .(gene_a, gene_b)]
  edges[, weight := as.numeric(weight)]
  edges
}

extract_modules <- function(spat_cor, cluster_name) {
  clusters <- NULL
  if (!is.null(spat_cor[["cor_clusters"]])) {
    clusters <- spat_cor[["cor_clusters"]][[cluster_name]]
  }
  if (is.null(clusters) || length(clusters) == 0) {
    stop("No clusters found in spatCorObject for name: ", cluster_name)
  }
  data.frame(
    gene = names(clusters),
    module_id = as.integer(clusters),
    stringsAsFactors = FALSE
  )
}

select_corr_genes <- function(expr_mat, n_corr_genes, method = c("mean", "variance")) {
  method <- match.arg(method)
  n_corr_genes <- suppressWarnings(as.integer(n_corr_genes))
  if (!is.finite(n_corr_genes) || n_corr_genes < 2L) stop("--n_corr_genes must be >= 2.")

  genes <- rownames(expr_mat)
  if (is.null(genes) || !length(genes)) stop("Expression matrix missing gene rownames.")
  genes <- as.character(genes)

  if (inherits(expr_mat, "dgCMatrix")) {
    mu <- Matrix::rowMeans(expr_mat)
    if (method == "mean") {
      score <- mu
    } else {
      # Avoid allocating expr_mat^2 (can be large): sum squares via dgCMatrix slots.
      n_obs <- as.integer(ncol(expr_mat))
      if (!is.finite(n_obs) || n_obs < 1L) stop("Expression matrix has no columns (cells).")
      rsq <- Matrix::rowSums(expr_mat ^ 2)
      m2 <- rsq / n_obs
      score <- m2 - mu ^ 2
      score[score < 0] <- 0
    }
  } else if (inherits(expr_mat, "Matrix")) {
    mu <- Matrix::rowMeans(expr_mat)
    if (method == "mean") {
      score <- mu
    } else {
      m2 <- Matrix::rowMeans(expr_mat ^ 2)
      score <- m2 - mu ^ 2
      score[score < 0] <- 0
    }
  } else {
    mu <- base::rowMeans(expr_mat)
    if (method == "mean") {
      score <- mu
    } else {
      m2 <- base::rowMeans(expr_mat ^ 2)
      score <- m2 - mu ^ 2
      score[score < 0] <- 0
    }
  }

  dt <- data.table(gene = genes, score = suppressWarnings(as.numeric(score)))
  dt <- dt[nzchar(gene) & is.finite(score)]
  if (nrow(dt) < 2) stop("Too few genes available after scoring for gene selection.")
  data.table::setorder(dt, -score, gene)
  n_take <- min(n_corr_genes, nrow(dt))
  selected <- dt$gene[seq_len(n_take)]
  selected <- selected[nzchar(selected)]
  selected <- unique(selected)
  if (length(selected) < 2) stop("Too few selected genes (n<2) after filtering.")
  selected
}

spec <- matrix(c(
  "data_dir",   "d", 1, "character",
  "outdir",     "o", 1, "character",
  "coord_file", "c", 1, "character",
  "roi_flip_y", "F", 0, "logical",
  "ncores",     "n", 1, "integer",
  "seed",       "s", 1, "integer",
  "dataset_id", "D", 1, "character",
  "roi_id",     "r", 1, "character",
  "max_cells",  "m", 1, "integer",
  "load_expression", "e", 0, "integer",
  "load_cellmeta", "l", 0, "integer",
  "load_transcripts", "t", 0, "integer",
  # Compatibility flags (ignored or validated; grid workflow only)
  "coexpr_mode", "Q", 1, "character",
  "detect_method", "q", 1, "character",
  "grid_stepsize", "Z", 1, "double",
  "grid_um", "u", 1, "double",
  "spatial_grid_um", "U", 1, "double",
  "spatial_grid_name", "Y", 1, "character",
  "min_cells_per_grid", "P", 1, "integer",
  "n_corr_genes", "G", 1, "integer",
  "n_spatial_genes", "g", 1, "integer",
  "gene_select", "M", 1, "character",
  "expression_values", "V", 1, "character",
  "k_modules",  "k", 1, "integer",
  "cor_method", "C", 1, "character",
  "emit_edge_stats", "E", 1, "integer",
  "pattern_dimensions", "J", 1, "character",
  "top_pos_genes", "A", 1, "integer",
  "top_neg_genes", "T", 1, "integer",
  "min_pos_cor", "p", 1, "double",
  "min_neg_cor", "N", 1, "double",
  "n_pattern_genes_max", "W", 1, "integer",
  "k_neighbors", "K", 1, "integer",
  "maximum_distance_knn", "X", 1, "double",
  "bin_method", "B", 1, "character",
  "calc_hub", "H", 1, "integer",
  "hub_min_int", "I", 1, "integer",
  "stats_file", "S", 1, "character"
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

if (is.null(opt$data_dir) || is.null(opt$outdir)) {
  cat("Usage: run_giotto_xenium.R --data_dir <Xenium outs> --outdir <out>\n", file = stderr())
  quit(status = 2)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

coord_path <- if (!is.null(opt$coord_file) && nzchar(opt$coord_file)) opt$coord_file else ""
roi_applied <- if (nzchar(coord_path)) 1L else 0L
roi_flip_y <- !is.null(opt$roi_flip_y) && !identical(opt$roi_flip_y, FALSE) && opt$roi_flip_y != 0

dataset_id <- if (!is.null(opt$dataset_id) && nzchar(opt$dataset_id)) opt$dataset_id else basename(opt$data_dir)
roi_id <- if (!is.null(opt$roi_id) && nzchar(opt$roi_id)) {
  opt$roi_id
} else if (nzchar(coord_path)) {
  tools::file_path_sans_ext(basename(coord_path))
} else {
  "full"
}

ncores <- as_int(opt$ncores, 1L)
if (!is.finite(ncores) || ncores < 1L) ncores <- 1L
data.table::setDTthreads(threads = ncores)
eprint("[INFO] data.table threads set to %d", ncores)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(ncores)
  RhpcBLASctl::omp_set_num_threads(ncores)
  eprint("[INFO] RhpcBLASctl threads set to %d", ncores)
}

seed <- as_int(opt$seed, 1L)
set.seed(seed)

max_cells <- as_int(opt$max_cells, 0L)
if (max_cells < 0L) max_cells <- 0L

load_expression <- if (is.null(opt$load_expression)) TRUE else (as_int(opt$load_expression, 1L) != 0L)
load_cellmeta <- if (is.null(opt$load_cellmeta)) FALSE else (as_int(opt$load_cellmeta, 0L) != 0L)
load_transcripts <- if (is.null(opt$load_transcripts)) FALSE else (as_int(opt$load_transcripts, 1L) != 0L)

requested_modes <- character()
if (!is.null(opt$coexpr_mode) && nzchar(opt$coexpr_mode)) {
  requested_modes <- c(requested_modes, tolower(as.character(opt$coexpr_mode)))
}
if (!is.null(opt$detect_method) && nzchar(opt$detect_method)) {
  requested_modes <- c(requested_modes, tolower(as.character(opt$detect_method)))
}
requested_modes <- unique(requested_modes[nzchar(requested_modes)])
if (length(requested_modes) > 0 && any(requested_modes != "grid")) {
  stop("This runner is grid-only; --coexpr_mode/--detect_method must be 'grid' (requested: ",
       paste(requested_modes, collapse = ","), ").")
}

grid_stepsize_user <- as_num(opt$grid_stepsize, NA_real_)
if (!is.finite(grid_stepsize_user)) grid_stepsize_user <- NA_real_
if (is.finite(grid_stepsize_user) && grid_stepsize_user <= 0) {
  stop("--grid_stepsize must be > 0 (or omitted to use auto bbox logic).")
}
if (!is.finite(grid_stepsize_user)) {
  spatial_grid_um <- as_num(opt$spatial_grid_um, NA_real_)
  if (is.finite(spatial_grid_um) && spatial_grid_um > 0) {
    grid_stepsize_user <- spatial_grid_um
  }
}
if (!is.finite(grid_stepsize_user)) {
  grid_um <- as_num(opt$grid_um, NA_real_)
  if (is.finite(grid_um) && grid_um > 0) {
    grid_stepsize_user <- grid_um
  }
}

spatial_grid_name <- if (!is.null(opt$spatial_grid_name) && nzchar(opt$spatial_grid_name)) {
  opt$spatial_grid_name
} else {
  "spatial_grid"
}

min_cells_per_grid <- as_int(opt$min_cells_per_grid, 4L)
if (!is.finite(min_cells_per_grid) || min_cells_per_grid < 1L) min_cells_per_grid <- 4L

n_corr_genes <- as_int(opt$n_corr_genes, N_CORR_GENES_DEFAULT)
if (!is.finite(n_corr_genes) || n_corr_genes < 2L) n_corr_genes <- N_CORR_GENES_DEFAULT
if (is.null(opt$n_corr_genes) && !is.null(opt$n_spatial_genes)) {
  n_corr_genes <- as_int(opt$n_spatial_genes, n_corr_genes)
}

gene_select <- if (!is.null(opt$gene_select) && nzchar(opt$gene_select)) tolower(opt$gene_select) else "mean"
if (!gene_select %in% c("mean", "variance")) stop("Invalid --gene_select: ", gene_select, " (expected mean|variance)")

expression_values_for_corr <- if (!is.null(opt$expression_values) && nzchar(opt$expression_values)) tolower(opt$expression_values) else "raw"
if (!expression_values_for_corr %in% c("raw", "counts")) {
  stop("Invalid --expression_values: ", expression_values_for_corr, " (expected raw|counts)")
}

k_modules <- as_int(opt$k_modules, 5L)
if (!is.finite(k_modules) || k_modules < 1L) k_modules <- 5L

cor_method <- if (!is.null(opt$cor_method) && nzchar(opt$cor_method)) tolower(opt$cor_method) else "pearson"
if (cor_method != "pearson") stop("This runner supports only --cor_method pearson.")

emit_edge_stats <- as_int(opt$emit_edge_stats, 0L)
if (!is.finite(emit_edge_stats) || emit_edge_stats < 0L) emit_edge_stats <- 0L
if (emit_edge_stats != 0L) {
  eprint("[WARN] --emit_edge_stats ignored (edges_all.tsv uses constant p_value/fdr=1.0).")
}

stats_file <- if (!is.null(opt$stats_file) && nzchar(opt$stats_file)) opt$stats_file else ""

roi_source_sha256 <- ""
if (nzchar(coord_path) && file.exists(coord_path)) {
  roi_source_sha256 <- sha256_file(coord_path)
}

start_time <- Sys.time()

if (!isTRUE(load_expression)) {
  stop("This grid-only PC runner requires a precomputed cell-level expression matrix. Set --load_expression 1.")
}

eprint("[INFO] Loading Xenium outs: %s", opt$data_dir)
eprint("[INFO] Xenium load_expression=%s load_transcripts=%s load_cellmeta=%s",
       ifelse(isTRUE(load_expression), "TRUE", "FALSE"),
       ifelse(isTRUE(load_transcripts), "TRUE", "FALSE"),
       ifelse(isTRUE(load_cellmeta), "TRUE", "FALSE"))

expr_h5_use <- NULL
expr_h5 <- resolve_xenium_expression_h5(opt$data_dir)
if (!is.null(expr_h5) && nzchar(expr_h5)) {
  expr_h5_use <- copy_h5_if_readonly(expr_h5)
  if (!identical(expr_h5_use, expr_h5)) {
    eprint("[INFO] Copied expression H5 to writable path: %s", expr_h5_use)
  }
}

create_opts <- list(
  xenium_dir = opt$data_dir,
  load_expression = isTRUE(load_expression),
  load_transcripts = isTRUE(load_transcripts),
  load_cellmeta = isTRUE(load_cellmeta),
  load_images = NULL,
  verbose = FALSE
)
if (!is.null(expr_h5_use) && nzchar(expr_h5_use)) {
  create_opts$expression_path <- expr_h5_use
}

gobject <- call_giotto_name("createGiottoXeniumObject", create_opts)

spat_locs_all <- getSpatialLocations(
  gobject,
  spat_unit = "cell",
  output = "data.table",
  copy_obj = TRUE,
  verbose = FALSE
)
cell_col <- find_col(spat_locs_all, c("cell_ID", "cell_id", "cell", "barcode", "obs_id"))
if (is.null(cell_col)) stop("Cannot find cell id column in spatial locations.")

x_col_all <- find_col(spat_locs_all, c("sdimx", "x", "x_centroid", "x_coord"))
y_col_all <- find_col(spat_locs_all, c("sdimy", "y", "y_centroid", "y_coord"))
if (!is.null(x_col_all) && !is.null(y_col_all)) {
  cell_bbox <- bbox_from_xy(spat_locs_all[[x_col_all]], spat_locs_all[[y_col_all]])
  eprint("[INFO] Cell bbox: %s", format_bbox(cell_bbox))
}

all_cell_ids <- as.character(spat_locs_all[[cell_col]])
all_cell_ids <- all_cell_ids[nzchar(all_cell_ids)]
n_obs_raw <- length(all_cell_ids)
if (n_obs_raw == 0) stop("No cells found in Giotto object.")
eprint("[INFO] cells loaded: %d", n_obs_raw)

roi_cell_ids <- all_cell_ids
roi_hash <- ""
obs_id_hash <- ""
n_obs_roi <- n_obs_raw

if (nzchar(coord_path)) {
  if (!file.exists(coord_path)) stop("Missing: ", coord_path)
  roi_vertices <- read_roi_polygon(coord_path)
  eprint("[INFO] ROI bbox (raw): %s", format_bbox(bbox_from_xy(roi_vertices[, 1], roi_vertices[, 2])))

  if (roi_flip_y) {
    roi_vertices[, 2] <- -roi_vertices[, 2]
    eprint("[INFO] ROI bbox (flip_y): %s", format_bbox(bbox_from_xy(roi_vertices[, 1], roi_vertices[, 2])))
  }

  res <- select_roi_cells(gobject, roi_vertices, ifelse(roi_flip_y, "roi_flip", "roi"))
  gobject <- res$gobject
  roi_cell_ids <- res$cell_ids
  n_obs_roi <- length(roi_cell_ids)
  eprint("[INFO] ROI clip on cells: %d -> %d", n_obs_raw, n_obs_roi)
  if (n_obs_roi == 0) {
    stop("No cells found within ROI polygon. If coordinates are inverted, retry with --roi_flip_y.")
  }
}

obs_ids_used <- roi_cell_ids
downsample_applied <- 0L
if (max_cells > 0L && length(obs_ids_used) > max_cells) {
  set.seed(seed)
  obs_ids_used <- sample(obs_ids_used, max_cells)
  downsample_applied <- 1L
}

obs_ids_used <- as.character(obs_ids_used)
obs_ids_used <- obs_ids_used[nzchar(obs_ids_used)]
obs_ids_used <- unique(obs_ids_used)

n_obs_selected <- length(obs_ids_used)
if (n_obs_selected == 0) stop("No cells selected after ROI/downsample.")
eprint("[INFO] n_obs_selected: %d", n_obs_selected)

obs_id_hash <- sha256_text(paste(sort(obs_ids_used), collapse = "\n"))
if (nzchar(obs_id_hash)) {
  roi_hash <- sha256_text(paste0(roi_source_sha256, "\n", obs_id_hash))
}

if (roi_applied == 1L || downsample_applied == 1L) {
  gobject <- subsetGiotto(gobject, cell_ids = obs_ids_used)
}

spat_unit_use <- tryCatch(set_default_spat_unit(gobject = gobject, spat_unit = NULL), error = function(e) "cell")
if (is.null(spat_unit_use) || !nzchar(spat_unit_use)) spat_unit_use <- "cell"

feat_type_use <- tryCatch(
  set_default_feat_type(gobject = gobject, spat_unit = spat_unit_use, feat_type = NULL),
  error = function(e) "rna"
)
if (is.null(feat_type_use) || !nzchar(feat_type_use)) feat_type_use <- "rna"
eprint("[INFO] default spat_unit=%s feat_type=%s", spat_unit_use, feat_type_use)

expr_raw <- tryCatch(
  getExpression(
    gobject = gobject,
    spat_unit = spat_unit_use,
    feat_type = feat_type_use,
    values = expression_values_for_corr,
    output = "matrix"
  ),
  error = function(e) NULL
)
if (is.null(expr_raw)) {
  stop("No expression matrix found for values=", expression_values_for_corr, " (expected raw|counts).")
}
if (is.null(rownames(expr_raw)) || is.null(colnames(expr_raw))) {
  stop("Expression matrix missing rownames (genes) or colnames (cell_IDs).")
}

n_obs_used <- as.integer(ncol(expr_raw))
n_vars_used <- as.integer(nrow(expr_raw))
eprint("[INFO] expression_values=%s n_obs_used=%d n_vars_used=%d", expression_values_for_corr, n_obs_used, n_vars_used)

selected_genes <- select_corr_genes(expr_raw, n_corr_genes = n_corr_genes, method = gene_select)
eprint("[INFO] selected_genes=%d (gene_select=%s n_corr_genes=%d)", length(selected_genes), gene_select, n_corr_genes)

if (k_modules > (length(selected_genes) - 1L)) {
  stop("--k_modules=", k_modules, " must be <= selected_genes-1=", length(selected_genes) - 1L)
}

spat_loc_name_use <- "raw"
grid_stepsize_effective <- NA_real_
grid_stepsize_effective_source <- "auto_bbox_target_bins"

t_step <- step_start("createSpatialGrid")
if (is.finite(grid_stepsize_user)) {
  grid_stepsize_effective <- grid_stepsize_user
  grid_stepsize_effective_source <- "user"
} else {
  target_bins <- max(1L, as.integer(floor(n_obs_selected / max(1L, min_cells_per_grid))))
  grid_stepsize_effective <- compute_auto_grid_stepsize(
    gobject = gobject,
    spat_unit = spat_unit_use,
    spat_loc_name = spat_loc_name_use,
    target_bins = target_bins
  )
}
eprint(
  "[INFO] spatial_grid stepsize_user=%s stepsize_effective=%.6g source=%s",
  ifelse(is.finite(grid_stepsize_user), format(grid_stepsize_user, digits = 8), "NA"),
  grid_stepsize_effective,
  grid_stepsize_effective_source
)
grid_opts <- list(
  gobject = gobject,
  spat_unit = spat_unit_use,
  spat_loc_name = spat_loc_name_use,
  name = spatial_grid_name,
  sdimx_stepsize = grid_stepsize_effective,
  sdimy_stepsize = grid_stepsize_effective,
  minimum_padding = 1,
  return_gobject = TRUE,
  verbose = FALSE
)
gobject <- call_giotto_name("createSpatialGrid", grid_opts)
step_end("createSpatialGrid", t_step)

t_step <- step_start("detectSpatialCorFeats_grid")
spat_cor_obj <- detectSpatialCorFeats_grid_compat(
  gobject = gobject,
  spat_unit = spat_unit_use,
  feat_type = feat_type_use,
  spat_loc_name = spat_loc_name_use,
  expression_values = expression_values_for_corr,
  subset_feats = selected_genes,
  spatial_grid_name = spatial_grid_name,
  min_cells_per_grid = min_cells_per_grid,
  cor_method = cor_method
)
step_end("detectSpatialCorFeats_grid", t_step)

cluster_name <- "spat_grid_clus"
t_step <- step_start("clusterSpatialCorFeats")
spat_cor_obj <- call_giotto_name(
  "clusterSpatialCorFeats",
  list(
    spatCorObject = spat_cor_obj,
    name = cluster_name,
    k = k_modules,
    verbose = FALSE
  )
)
step_end("clusterSpatialCorFeats", t_step)

edges_n_obs_source <- "grid_bins>=min_cells_per_grid"
edges_n_obs_try <- tryCatch(
  compute_grid_n_obs(
    gobject = gobject,
    spat_unit = spat_unit_use,
    feat_type = feat_type_use,
    spat_loc_name = spat_loc_name_use,
    spatial_grid_name = spatial_grid_name,
    min_cells_per_grid = min_cells_per_grid
  ),
  error = function(e) {
    eprint("[WARN] compute_grid_n_obs failed: %s", conditionMessage(e))
    NA_integer_
  }
)
edges_n_obs <- suppressWarnings(as.integer(edges_n_obs_try))
if (!is.finite(edges_n_obs) || edges_n_obs < 0L) {
  edges_n_obs <- as.integer(n_obs_used)
  edges_n_obs_source <- "fallback_n_obs_used"
}

modules_df_core <- extract_modules(spat_cor_obj, cluster_name)
expr_genes <- as.character(rownames(expr_raw))
modules_dt <- data.table(
  gene = expr_genes,
  module_id = rep(-1L, length(expr_genes))
)
idx_all <- match(modules_df_core$gene, modules_dt$gene)
keep_idx <- !is.na(idx_all)
modules_dt$module_id[idx_all[keep_idx]] <- as.integer(modules_df_core$module_id[keep_idx])

edges_all_base <- extract_complete_edges(spat_cor_obj, gene_set = selected_genes)
edges_all_dt <- data.table::as.data.table(edges_all_base)
if (!all(c("gene_a", "gene_b", "weight") %in% names(edges_all_dt))) {
  stop("extract_complete_edges returned unexpected schema.")
}
edges_all_dt[, `:=`(
  weight = suppressWarnings(as.numeric(weight)),
  p_value = 1.0,
  fdr = 1.0,
  n_obs = as.integer(edges_n_obs),
  edge_type = "spat_cor_grid"
)]
setcolorder(edges_all_dt, c("gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"))
if (!identical(names(edges_all_dt), c("gene_a", "gene_b", "weight", "p_value", "fdr", "n_obs", "edge_type"))) {
  stop("edges_all.tsv schema mismatch after setcolorder.")
}

fwrite(edges_all_dt, file.path(opt$outdir, "edges_all.tsv"), sep = "\t")
fwrite(modules_dt, file.path(opt$outdir, "modules.tsv"), sep = "\t")

spat_locs_use <- getSpatialLocations(
  gobject,
  spat_unit = "cell",
  output = "data.table",
  copy_obj = TRUE,
  verbose = FALSE
)
cell_col_use <- find_col(spat_locs_use, c("cell_ID", "cell_id", "cell", "barcode", "obs_id"))
x_col <- find_col(spat_locs_use, c("sdimx", "x", "x_centroid", "x_coord"))
y_col <- find_col(spat_locs_use, c("sdimy", "y", "y_centroid", "y_coord"))
if (is.null(cell_col_use) || is.null(x_col) || is.null(y_col)) {
  stop("Cannot resolve obs_id/x/y columns in spatial locations.")
}
coords_dt <- data.table(
  obs_id = as.character(spat_locs_use[[cell_col_use]]),
  x = as.numeric(spat_locs_use[[x_col]]),
  y = as.numeric(spat_locs_use[[y_col]])
)
obs_ids_expr <- colnames(expr_raw)
ord <- match(obs_ids_expr, coords_dt$obs_id)
if (any(is.na(ord))) stop("obs_id mismatch between expression matrix and coordinates.")
coords_dt <- coords_dt[ord, ]
fwrite(coords_dt, file.path(opt$outdir, "obs_coords.tsv"), sep = "\t")

method_version <- as.character(packageVersion("Giotto"))
wall <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

meta <- list(
  workflow = "giotto_grid_pc",
  method = "giotto",
  method_version = method_version,
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  input_dataset_id = dataset_id,
  roi_id = roi_id,
  roi_applied = roi_applied,
  roi_flip_y = roi_flip_y,
  downsample_applied = downsample_applied,
  max_cells = as.integer(max_cells),
  load_expression = isTRUE(load_expression),
  load_transcripts = isTRUE(load_transcripts),
  load_cellmeta = isTRUE(load_cellmeta),
  expression_values_for_corr = expression_values_for_corr,
  gene_select = gene_select,
  n_corr_genes = as.integer(n_corr_genes),
  selected_genes_n = as.integer(length(selected_genes)),
  selected_genes_sha256 = sha256_text(paste(selected_genes, collapse = "\n")),
  k_modules = as.integer(k_modules),
  n_obs_used = as.integer(n_obs_used),
  edges_n_obs = as.integer(edges_n_obs),
  edges_n_obs_source = edges_n_obs_source,
  edges_gene_set = "selected_genes",
  grid_stepsize_user = if (is.finite(grid_stepsize_user)) grid_stepsize_user else NA_real_,
  grid_stepsize_effective = if (is.finite(grid_stepsize_effective)) grid_stepsize_effective else NA_real_,
  grid_stepsize_effective_source = grid_stepsize_effective_source,
  min_cells_per_grid = as.integer(min_cells_per_grid),
  edges_all_outfile = "edges_all.tsv",
  modules_outfile = "modules.tsv",
  random_seed = as.integer(seed),
  cor_method = cor_method,
  runtime = list(
    wall_time_sec = wall,
    n_threads = as.integer(ncores)
  ),
  resources = list(stats_file = stats_file),
  roi_source_path = coord_path,
  roi_source_sha256 = roi_source_sha256,
  obs_id_hash = obs_id_hash,
  roi_hash = roi_hash
)

jsonlite::write_json(meta, file.path(opt$outdir, "meta.json"), auto_unbox = TRUE, pretty = TRUE)

eprint("[DONE] Giotto grid workflow completed. wall_sec=%.2f", wall)

#' Loess Residual Bootstrap Compat
#' @description
#' Internal helper for `.loess_residual_bootstrap_compat`.
#' @param x Parameter value.
#' @param y Parameter value.
#' @param strat Parameter value.
#' @param grid Parameter value.
#' @param B Parameter value.
#' @param span Parameter value.
#' @param degree Parameter value.
#' @param nthreads Parameter value.
#' @param ... Additional arguments (currently unused).
#' @return Return value used internally.
#' @keywords internal
.loess_residual_bootstrap_compat <- function(x, y, strat, grid,
                                             B = 1000L, span = 0.45,
                                             degree = 1L, nthreads = 1L,
                                             ...) {
    .loess_residual_bootstrap(x, y, strat, grid,
        B = B, span = span,
        deg = degree,
        n_threads = nthreads,
        ...
    )
}

#' Loess Residual Bootstrap Compat V2
#' @description
#' Internal helper for `.loess_residual_bootstrap_compat_v2`.
#' @param ... Additional arguments (currently unused).
#' @keywords internal
#' @return Return value used internally.
.loess_residual_bootstrap_compat_v2 <- function(...) {
    .loess_residual_bootstrap_compat(...)
}

#' Compute per-grid gene density
#' @description
#' Internal helper for `.compute_density`.
#' Aggregates counts for a gene set within each grid tile and stores density
#' per unit area under `scope_obj@density`.
#' @param scope_obj A `scope_object` with grid data.
#' @param grid_name Grid layer to use (defaults to active layer).
#' @param density_name Output column name for the density values.
#' @param genes Optional genes to aggregate.
#' @param cluster_col Optional meta column to pick genes by cluster.
#' @param cluster_num Optional cluster identifier to select genes.
#' @param layer_name Layer containing counts (default `counts`).
#' @param normalize_method Normalization mode (`none`, `per_grid`, `global_gene`).
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @keywords internal
.compute_density <- function(scope_obj,
                           grid_name = NULL,
                           density_name = "density",
                           genes = NULL,
                           cluster_col = NULL,
                           cluster_num = NULL,
                           layer_name = "counts",
                           normalize_method = c("none", "per_grid", "global_gene"),
                           verbose = TRUE) {
    normalize_method <- match.arg(normalize_method)

    step_s01 <- .log_step("computeDensity", "S01", "select grid layer", verbose)
    step_s01$enter(extra = paste0(
        "grid_name=", if (is.null(grid_name)) "auto" else grid_name,
        " layer_name=", layer_name,
        " normalize_method=", normalize_method
    ))

    ## --------------------------------------------------------------------- 2
    ## pick grid layer & sanity check
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_layer_name <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]

    .check_grid_content(scope_obj, grid_layer_name)

    grid_info <- g_layer$grid_info
    counts_dt <- g_layer[[layer_name]]
    if (is.null(counts_dt)) {
        stop("Layer '", layer_name, "' not found in grid '", grid_layer_name, "'.")
    }

    if (!all(c("grid_id", "gene", "count") %in% colnames(counts_dt))) {
        stop("Layer '", layer_name, "' must contain columns grid_id, gene, count.")
    }

    counts_dt <- as.data.table(counts_dt)

    ## --------------------------------------------------------------------- 3
    ## decide gene set
    sel_genes <- .get_gene_subset(scope_obj,
        genes        = genes,
        cluster_col  = cluster_col,
        cluster_num  = cluster_num
    )

    .log_info("computeDensity", "S01",
        paste0("Processing ", length(sel_genes), " genes"),
        verbose
    )
    step_s01$done(extra = paste0("grid=", grid_layer_name, " n_genes=", length(sel_genes)))

    ## --------------------------------------------------------------------- 4
    ## optional normalisation
    step_s02 <- .log_step("computeDensity", "S02", "apply normalization", verbose)
    step_s02$enter(extra = paste0("method=", normalize_method))
    .log_backend("computeDensity", "S02", "R",
        paste0("data.table normalize_method=", normalize_method),
        verbose = verbose
    )
    if (normalize_method != "none") {
        .log_info("computeDensity", "S02",
            paste0("Applying ", normalize_method, " normalization"),
            verbose
        )
    } else {
        .log_info("computeDensity", "S02", "Normalization disabled; keeping raw counts.", verbose)
    }

    counts_proc <- copy(counts_dt)

    if (normalize_method == "per_grid") {
        tot_grid <- counts_proc[, .(total = sum(count)), by = grid_id]
        counts_proc <- counts_proc[tot_grid, on = "grid_id"]
        counts_proc[, count := fifelse(total == 0, 0, count / total)][, total := NULL]
    } else if (normalize_method == "global_gene") {
        tot_gene <- counts_proc[, .(g_tot = sum(count)), by = gene]
        counts_proc <- counts_proc[tot_gene, on = "gene"]
        counts_proc[, count := fifelse(g_tot == 0, 0, count / g_tot)][, g_tot := NULL]
    }
    ## else "none": keep raw counts
    step_s02$done()

    ## --------------------------------------------------------------------- 5
    ## sum across selected genes -> per-grid counts
    step_s03 <- .log_step("computeDensity", "S03", "aggregate counts", verbose)
    step_s03$enter(extra = paste0("n_genes=", length(sel_genes)))
    sel_counts <- counts_proc[gene %in% sel_genes,
        .(count = sum(count)),
        by = grid_id
    ]
    step_s03$done(extra = paste0("n_grids=", length(unique(sel_counts$grid_id))))

    ## --------------------------------------------------------------------- 6
    ## attach cell area & compute density
    step_s04 <- .log_step("computeDensity", "S04", "compute density", verbose)
    step_s04$enter(extra = paste0("n_grids=", nrow(grid_info)))
    grid_dt <- as.data.table(grid_info)[
        ,
        .(grid_id, xmin, xmax, ymin, ymax,
            area = (xmax - xmin) * (ymax - ymin)
        )
    ]

    dens_dt <- merge(grid_dt, sel_counts, by = "grid_id", all.x = TRUE)
    dens_dt[is.na(count), count := 0]
    dens_dt[, density := count / area]
    step_s04$done()

    ## --------------------------------------------------------------------- 7
    ## create / update @density entry for this grid layer
    step_s05 <- .log_step("computeDensity", "S05", "store density outputs", verbose)
    step_s05$enter(extra = paste0("density_name=", density_name))
    all_ids <- grid_dt$grid_id
    dens_df <- if (!is.null(scope_obj@density[[grid_layer_name]])) {
        scope_obj@density[[grid_layer_name]]
    } else {
        data.frame(row.names = all_ids)
    }

    # ensure all grids represented
    if (nrow(dens_df) == 0) {
        dens_df <- data.frame(row.names = all_ids)
    }

    dens_vec <- dens_dt$density
    names(dens_vec) <- dens_dt$grid_id
    dens_df[[density_name]] <- dens_vec[rownames(dens_df)]
    dens_df[[density_name]][is.na(dens_df[[density_name]])] <- 0

    scope_obj@density[[grid_layer_name]] <- dens_df

    .log_info("computeDensity", "S05", "Density computation completed", verbose)
    step_s05$done()

    invisible(scope_obj)
}

#' Compute gene-gene Jaccard heatmap
#' @description
#' Internal helper for `.compute_gene_gene_jaccard_heatmap`.
#' Builds a binary gene x bin matrix and renders a Jaccard similarity heatmap,
#' optionally stratified by cluster annotations.
#' @param scope_obj A `scope_object` with grid counts or Xz matrix.
#' @param grid_name Grid layer name (default `grid30`).
#' @param expr_layer Expression layer to use (`Xz` or `counts`).
#' @param min_bins_per_gene Minimum bins per gene required to include.
#' @param top_n_genes Number of top genes to plot when no whitelist is given.
#' @param genes_of_interest Optional gene whitelist.
#' @param output_dir Optional directory to save outputs.
#' @param cluster_col Optional meta column for cluster labels.
#' @param split_heatmap_by_cluster Whether to split heatmap by cluster.
#' @param heatmap_bg Parameter value.
#' @param keep_cluster_na Logical flag.
#' @return Invisibly returns the heatmap object (if created).
#' @keywords internal
.compute_gene_gene_jaccard_heatmap <- function(scope_obj,
                                          grid_name = "grid30",
                                          expr_layer = c("Xz", "counts"),
                                          min_bins_per_gene = 8L,
                                          top_n_genes = 60L,
                                          genes_of_interest = NULL,
                                          output_dir = NULL,
                                          heatmap_bg = "#c0c0c0",
                                          cluster_col = NULL,
                                          keep_cluster_na = TRUE,
                                          split_heatmap_by_cluster = FALSE,
                                          verbose = getOption("geneSCOPE.verbose", TRUE)) {
  expr_layer <- match.arg(expr_layer)
  step_s01 <- .log_step("computeGeneGeneJaccardHeatmap", "S01", "validate inputs and dependencies", verbose)
  step_s01$enter(extra = paste0("grid_name=", grid_name, " expr_layer=", expr_layer))
  dep_heatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)
  .log_backend("computeGeneGeneJaccardHeatmap", "S01", "dependency",
    paste0("ComplexHeatmap available=", dep_heatmap),
    verbose = verbose
  )
  if (!dep_heatmap) {
    stop("Package 'ComplexHeatmap' is required to draw the Jaccard heatmap.")
  }
  dep_circlize <- requireNamespace("circlize", quietly = TRUE)
  .log_backend("computeGeneGeneJaccardHeatmap", "S01", "dependency",
    paste0("circlize available=", dep_circlize),
    verbose = verbose
  )
  if (!dep_circlize) {
    stop("Package 'circlize' is required to build the heatmap palette.")
  }
  step_s01$done()

  cluster_lookup <- NULL
  if (!is.null(cluster_col)) {
    meta_df <- scope_obj@meta.data
    if (is.null(meta_df) || !is.data.frame(meta_df)) {
      warning("scope_obj@meta.data is missing or not a data.frame; ignore cluster_col.")
    } else if (!(cluster_col %in% colnames(meta_df))) {
      warning("cluster column ", cluster_col, " is not present in meta.data; ignore this parameter.")
    } else if (is.null(rownames(meta_df))) {
      warning("meta.data has no rownames, cannot map gene -> cluster; ignore cluster_col.")
    } else {
      cluster_lookup <- as.character(meta_df[[cluster_col]])
      names(cluster_lookup) <- rownames(meta_df)
    }
  }

  grid_layer <- scope_obj@grid[[grid_name]]
  if (is.null(grid_layer)) {
    warning("Grid layer ", grid_name, " is missing; cannot compute Jaccard.")
    return(invisible(NULL))
  }

  step_s02 <- .log_step("computeGeneGeneJaccardHeatmap", "S02", "build gene-by-bin matrix", verbose)
  step_s02$enter(extra = paste0("expr_layer=", expr_layer))
  gene_bin_mat <- NULL
  if (identical(expr_layer, "counts")) {
    if (is.null(grid_layer$counts)) {
      warning("Grid layer ", grid_name, " has no counts; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    counts_df <- as.data.frame(grid_layer$counts)
    if (!nrow(counts_df)) {
      warning("Grid layer ", grid_name, " has no gene x bin records.")
      return(invisible(NULL))
    }

    cols_expected <- c("gene", "grid_id", "count")
    if (!all(cols_expected %in% colnames(counts_df))) {
      warning("counts is missing required columns: ", paste(setdiff(cols_expected, colnames(counts_df)), collapse = ", "))
      return(invisible(NULL))
    }

    counts_df$gene <- as.character(counts_df$gene)
    counts_df$grid_id <- as.character(counts_df$grid_id)
    counts_df$count <- as.numeric(as.character(counts_df$count))
    counts_df <- counts_df[!is.na(counts_df$gene) & !is.na(counts_df$grid_id) & counts_df$count > 0, , drop = FALSE]
    if (!nrow(counts_df)) {
      warning("No gene x bin records remain after filtering NA/zero.")
      return(invisible(NULL))
    }

    counts_df <- aggregate(count ~ gene + grid_id, data = counts_df, FUN = sum)
    genes <- sort(unique(counts_df$gene))
    bins <- sort(unique(counts_df$grid_id))

    gene_idx <- match(counts_df$gene, genes)
    bin_idx <- match(counts_df$grid_id, bins)
    gene_bin_mat <- sparseMatrix(
      i = gene_idx,
      j = bin_idx,
      x = as.numeric(counts_df$count > 0),
      dims = c(length(genes), length(bins)),
      dimnames = list(genes, bins)
    )
  } else {
    if (is.null(grid_layer$Xz)) {
      warning("Grid layer ", grid_name, " has no Xz; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    xz_mat <- grid_layer$Xz
    xz_mat <- as.matrix(xz_mat)
    if (!nrow(xz_mat) || !ncol(xz_mat)) {
      warning(grid_name, " Xz layer is empty; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    if (is.null(colnames(xz_mat)) || is.null(rownames(xz_mat))) {
      warning("Xz layer lacks row/col names; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    bin_mask <- xz_mat > 0
    bin_mask[is.na(bin_mask)] <- FALSE
    gene_bin_mat <- Matrix(t(bin_mask) * 1, sparse = TRUE)
    rownames(gene_bin_mat) <- colnames(xz_mat)
    colnames(gene_bin_mat) <- rownames(xz_mat)
  }

  if (is.null(gene_bin_mat) || !nrow(gene_bin_mat)) {
    warning("No usable gene x bin information; cannot compute Jaccard.")
    return(invisible(NULL))
  }

  if (!is.null(genes_of_interest)) {
    genes_keep <- intersect(unique(genes_of_interest), rownames(gene_bin_mat))
    if (!length(genes_keep)) {
      warning("genes_of_interest not found in current grid.")
      return(invisible(NULL))
    }
    gene_bin_mat <- gene_bin_mat[genes_keep, , drop = FALSE]
  }

  gene_clusters <- NULL
  if (!is.null(cluster_lookup)) {
    gene_names_all <- rownames(gene_bin_mat)
    idx <- match(gene_names_all, names(cluster_lookup))
    found <- !is.na(idx)
    if (!any(found)) {
      warning("cluster column ", cluster_col, " has no overlapping genes with current grid.")
      return(invisible(NULL))
    }
    if (!keep_cluster_na) {
      found <- found & !is.na(cluster_lookup[idx])
      if (!any(found)) {
        warning("All overlapping genes have NA cluster labels; cannot continue Jaccard.")
        return(invisible(NULL))
      }
    }
    gene_bin_mat <- gene_bin_mat[found, , drop = FALSE]
    idx <- idx[found]
    gene_clusters <- cluster_lookup[idx]
    names(gene_clusters) <- rownames(gene_bin_mat)
  }

  if (!nrow(gene_bin_mat)) {
    warning("No genes remain after filtering.")
    return(invisible(NULL))
  }

  bin_presence <- rowSums(gene_bin_mat)
  keep_idx <- which(bin_presence >= min_bins_per_gene)
  if (!length(keep_idx)) {
    warning("No genes meet the minimum bin count of ", min_bins_per_gene, ".")
    return(invisible(NULL))
  }

  gene_bin_mat <- gene_bin_mat[keep_idx, , drop = FALSE]
  bin_presence <- bin_presence[keep_idx]
  if (!is.null(gene_clusters)) {
    gene_clusters <- gene_clusters[rownames(gene_bin_mat)]
  }

  if (!is.null(top_n_genes) && is.finite(top_n_genes) && top_n_genes > 0 &&
      nrow(gene_bin_mat) > top_n_genes) {
    ord <- order(bin_presence, decreasing = TRUE)
    keep_ord <- ord[seq_len(top_n_genes)]
    gene_bin_mat <- gene_bin_mat[keep_ord, , drop = FALSE]
    bin_presence <- bin_presence[keep_ord]
    if (!is.null(gene_clusters)) {
      gene_clusters <- gene_clusters[keep_ord]
    }
  }

  if (!is.null(gene_clusters) && split_heatmap_by_cluster) {
    cluster_order <- order(gene_clusters, rownames(gene_bin_mat))
    gene_bin_mat <- gene_bin_mat[cluster_order, , drop = FALSE]
    bin_presence <- bin_presence[cluster_order]
    gene_clusters <- gene_clusters[cluster_order]
  }

  .log_backend("computeGeneGeneJaccardHeatmap", "S02", "matrix",
    paste0("dgCMatrix n_genes=", nrow(gene_bin_mat), " n_bins=", ncol(gene_bin_mat)),
    verbose = verbose
  )
  step_s02$done()

  step_s03 <- .log_step("computeGeneGeneJaccardHeatmap", "S03", "compute Jaccard heatmap", verbose)
  step_s03$enter(extra = paste0("n_genes=", nrow(gene_bin_mat)))
  gene_names <- rownames(gene_bin_mat)
  inter_mat <- tcrossprod(gene_bin_mat)
  inter_dense <- as.matrix(inter_mat)
  union_dense <- outer(bin_presence, bin_presence, "+") - inter_dense
  jac_mat <- inter_dense
  positive_union <- union_dense > 0
  jac_mat[positive_union] <- jac_mat[positive_union] / union_dense[positive_union]
  jac_mat[!positive_union] <- NA_real_
  diag(jac_mat) <- 1
  rownames(jac_mat) <- gene_names
  colnames(jac_mat) <- gene_names

  col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#113A70", "#f7f7f7", "#B11226"))
  row_split_factor <- column_split_factor <- NULL
  if (!is.null(gene_clusters) && split_heatmap_by_cluster) {
    cluster_labels_plot <- ifelse(is.na(gene_clusters), "NA", as.character(gene_clusters))
    cluster_levels <- unique(cluster_labels_plot)
    row_split_factor <- factor(cluster_labels_plot, levels = cluster_levels)
    column_split_factor <- row_split_factor
  }
  column_title_txt <- NULL
  .log_backend("computeGeneGeneJaccardHeatmap", "S03", "R",
    "tcrossprod + ComplexHeatmap", verbose = verbose
  )
  ht <- ComplexHeatmap::Heatmap(
    jac_mat,
    name = "Jaccard",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    column_title = column_title_txt,
    column_title_gp = grid::gpar(fontsize = 10),
    row_split = row_split_factor,
    column_split = column_split_factor
  )
  step_s03$done()

  invisible(list(
    matrix = jac_mat,
    genes = gene_names,
    clusters = gene_clusters,
    grid = grid_name,
    min_bins = min_bins_per_gene,
    top_n = top_n_genes,
    heatmap = ht,
    palette = col_fun
  ))
}

#' Test within- vs between-module Jaccard similarity
#' @description
#' Internal helper for `.compute_gene_module_within_between_test`.
#' Performs Wilcoxon or permutation testing on Jaccard similarity matrices to
#' compare within-module versus between-module distributions.
#' @param jac_mat Jaccard similarity matrix.
#' @param gene_modules Named vector assigning genes to modules.
#' @param alternative Alternative hypothesis for the test.
#' @param method Testing method (`wilcoxon` or `permutation`).
#' @param n_perm Number of permutations when `method = \"permutation\"`.
#' @param seed Optional seed for reproducibility.
#' @param min_within Minimum within-module pairs required.
#' @param min_between Minimum between-module pairs required.
#' @return List with test statistics and p-values.
#' @keywords internal
.compute_gene_module_within_between_test <- function(jac_mat,
                                               gene_modules,
                                               alternative = c("greater", "two.sided", "less"),
                                               method = c("wilcoxon", "permutation"),
                                               n_perm = 10000L,
                                               seed = NULL,
                                               min_within = 1L,
                                               min_between = 1L,
                                               verbose = getOption("geneSCOPE.verbose", TRUE)) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  step_s01 <- .log_step("computeGeneModuleWithinBetweenTest", "S01", "validate inputs", verbose)
  step_s01$enter(extra = paste0("method=", method, " alternative=", alternative))

  if (is.null(jac_mat) || !is.matrix(jac_mat)) {
    stop("jac_mat must be a numeric matrix with gene names as dimnames.")
  }
  if (is.null(rownames(jac_mat)) || is.null(colnames(jac_mat))) {
    stop("jac_mat must have rownames/colnames for genes.")
  }
  if (is.null(gene_modules)) {
    stop("gene_modules (named vector of module labels) is required.")
  }
  if (is.null(names(gene_modules))) {
    stop("gene_modules must be named by gene.")
  }

  gene_names <- intersect(rownames(jac_mat), names(gene_modules))
  if (!length(gene_names)) {
    stop("No overlap between jac_mat genes and gene_modules names.")
  }
  step_s01$done(extra = paste0("n_genes=", length(gene_names)))

  jac_mat <- jac_mat[gene_names, gene_names, drop = FALSE]
  gene_modules <- gene_modules[gene_names]

  # Drop genes with NA module labels
  jac_mat <- jac_mat[!is.na(gene_modules), !is.na(gene_modules), drop = FALSE]
  gene_modules <- gene_modules[!is.na(gene_modules)]
  gene_names <- names(gene_modules)

  if (!length(gene_names)) {
    stop("All genes have NA module labels after filtering.")
  }

  step_s02 <- .log_step("computeGeneModuleWithinBetweenTest", "S02", "compute within/between summaries", verbose)
  step_s02$enter(extra = paste0("min_within=", min_within, " min_between=", min_between))
  per_gene <- lapply(seq_along(gene_names), function(i) {
    g <- gene_names[[i]]
    mod <- gene_modules[[i]]
    row_vals <- jac_mat[g, ]

    within_mask <- gene_modules == mod & names(gene_modules) != g
    between_mask <- gene_modules != mod

    within_vals <- row_vals[within_mask]
    between_vals <- row_vals[between_mask]

    within_vals <- within_vals[!is.na(within_vals)]
    between_vals <- between_vals[!is.na(between_vals)]

    if (length(within_vals) < min_within || length(between_vals) < min_between) {
      return(NULL)
    }

    data.frame(
      gene = g,
      module = mod,
      within_mean = mean(within_vals),
      between_mean = mean(between_vals),
      within_n = length(within_vals),
      between_n = length(between_vals),
      diff = mean(within_vals) - mean(between_vals),
      stringsAsFactors = FALSE
    )
  })

  per_gene_df <- do.call(rbind, per_gene)
  if (is.null(per_gene_df) || !nrow(per_gene_df)) {
    warning("No genes retained for within/between summary (check min_within/min_between).")
    return(invisible(NULL))
  }
  step_s02$done(extra = paste0("n_genes=", nrow(per_gene_df)))

  delta <- per_gene_df$diff
  delta <- delta[is.finite(delta)]

  if (!length(delta)) {
    warning("No finite within-between differences to test.")
    return(invisible(list(gene_summary = per_gene_df, test = NULL)))
  }

  step_s03 <- .log_step("computeGeneModuleWithinBetweenTest", "S03", "run within/between test", verbose)
  step_s03$enter(extra = paste0("method=", method, " n_perm=", n_perm))
  .log_backend("computeGeneModuleWithinBetweenTest", "S03", "R",
    paste0("test_method=", method),
    verbose = verbose
  )
  test_res <- switch(method,
    wilcoxon = {
      wilcox.test(delta, mu = 0, alternative = alternative, paired = FALSE, exact = FALSE)
    },
    permutation = {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      if (!is.finite(n_perm) || n_perm < 1) {
        stop("n_perm must be a positive integer for permutation test.")
      }
      obs <- mean(delta)
      perm_stats <- replicate(n_perm, {
        mean(delta * sample(c(-1, 1), length(delta), replace = TRUE))
      })
      p_val <- switch(alternative,
        greater = (sum(perm_stats >= obs) + 1) / (n_perm + 1),
        less = (sum(perm_stats <= obs) + 1) / (n_perm + 1),
        two.sided = (sum(abs(perm_stats) >= abs(obs)) + 1) / (n_perm + 1)
      )
      list(
        statistic = obs,
        p.value = p_val,
        method = "paired sign-flip permutation on per-gene within-minus-between means",
        alternative = alternative,
        n_perm = n_perm
      )
    }
  )
  step_s03$done()

  invisible(list(
    gene_summary = per_gene_df,
    test = test_res,
    alternative = alternative,
    method = method
  ))
}

#' Compute the Idelta statistic for spatial genes
#' @description
#' Internal helper for `.compute_idelta`.
#' Calculates Moran's Idelta for each gene at grid or cell level and stores results
#' under `scope_obj@stats`.
#' @param scope_obj A `scope_object` with grid or cell data.
#' @param grid_name Grid layer name when `level = \"grid\"`.
#' @param level Compute at `grid` or `cell` level.
#' @param ncores Number of threads to use.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @keywords internal
.compute_idelta <- function(scope_obj,
                          grid_name = NULL,
                          level = c("grid", "cell"),
                          ncores = 4,
                          verbose = getOption("geneSCOPE.verbose", TRUE)) {
    level <- match.arg(level)
    ncores_requested <- suppressWarnings(as.integer(ncores))
    cores_detected <- suppressWarnings(detectCores(logical = TRUE))
    ncores_safe <- max(1L, min(ncores_requested, cores_detected))
    ncores <- ncores_safe
    step_s01 <- .log_step("computeIDelta", "S01", "prepare sparse matrix", verbose)
    step_s01$enter(extra = paste0(
        "level=", level,
        " grid_name=", if (is.null(grid_name)) "auto" else grid_name
    ))

    if (level == "grid") {
        ## ---- 1. Select grid layer ---------------------------------------------------
        g_layer <- .select_grid_layer(scope_obj, grid_name)
        # Parse actual grid_name (selectGridLayer might auto-select)
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
        .log_info("computeIDelta", "S01",
            paste0("Selected grid layer=", grid_name),
            verbose
        )

        ## ---- 2. Build gene x grid sparse matrix ------------------------------------
        dt <- g_layer$counts[g_layer$counts$count > 0, ] # Keep only positive counts
        genes <- sort(unique(dt$gene))
        grids <- g_layer$grid_info$grid_id

        Gsp <- sparseMatrix(
            i = match(dt$gene, genes),
            j = match(dt$grid_id, grids),
            x = dt$count,
            dims = c(length(genes), length(grids)),
            dimnames = list(genes, grids)
        )

        storage_key <- grid_name
        col_suffix <- paste0(grid_name, "_iDelta")
    } else { # level == "cell"
        ## ---- 1. Check for cell count matrix ----------------------------------------
        if (is.null(scope_obj@cells) || is.null(scope_obj@cells$counts)) {
            stop("No cell count matrix found. Please run .add_single_cells() first.")
        }

        ## ---- 2. Use existing gene x cell sparse matrix -----------------------------
        Gsp <- scope_obj@cells$counts # Already genes x cells
        if (!inherits(Gsp, "dgCMatrix")) {
            stop("Cell count matrix must be a dgCMatrix.")
        }

        genes <- rownames(Gsp)

        storage_key <- "cell"
        col_suffix <- "cell_iDelta"
    }

    ## ---- 3. Call C++ kernel to compute Idelta -----------------------------------------
    .log_backend("computeIDelta", "S01", "matrix", "dgCMatrix sparse", verbose = verbose)
    step_s01$done(extra = paste0("n_genes=", length(genes), " n_cols=", ncol(Gsp)))

    step_s02 <- .log_step("computeIDelta", "S02", "memory guard and thread clamp", verbose)
    step_s02$enter(extra = paste0("threads=", ncores_safe))
    if (is.finite(ncores_safe) && is.finite(ncores_requested) && ncores_safe < ncores_requested) {
        .log_backend("computeIDelta", "S02", "threads",
            paste0(ncores_safe, " requested=", ncores_requested),
            reason = "clamped", verbose = verbose
        )
    } else {
        .log_backend("computeIDelta", "S02", "threads",
            paste0(ncores_safe, " requested=", ncores_requested),
            verbose = verbose
        )
    }
    # Memory guard: approximate sparse footprint (values + indices) per thread
    nnz <- length(Gsp@x)
    approx_gb <- (nnz * 16) / (1024^3) # rough estimate of numeric + index bytes
    sys_mem_gb <- .get_system_memory_gb()
    est_total_gb <- approx_gb * ncores
    .log_info("computeIDelta", "S02",
        paste0("Estimated sparse footprint: nnz=", nnz, " approx_gb=", round(approx_gb, 4)),
        verbose
    )
    if (est_total_gb > sys_mem_gb) {
        stop(
            "[geneSCOPE::.compute_idelta] Estimated memory requirement (",
            round(est_total_gb, 1), " GB) exceeds system capacity (",
            round(sys_mem_gb, 1), " GB). Reduce ncores or gene/cell counts."
        )
    }
    step_s02$done(extra = paste0("est_total_gb=", round(est_total_gb, 4)))

    step_s03 <- .log_step("computeIDelta", "S03", "compute I-delta", verbose)
    step_s03$enter(extra = paste0("n_genes=", length(genes)))
    .log_backend("computeIDelta", "S03", "C++",
        paste0("idelta_sparse_cpp threads=", ncores),
        verbose = verbose
    )
    delta_raw <- .idelta_sparse_cpp(Gsp, n_threads = ncores)
    names(delta_raw) <- genes
    step_s03$done()

    ## ---- 4. Ensure delta_raw is a numeric vector, not matrix ----------------------
    delta_raw <- as.numeric(delta_raw)
    names(delta_raw) <- genes

    ## ---- 5. Write to scope_obj -----------------------------------------------------
    step_s04 <- .log_step("computeIDelta", "S04", "store outputs", verbose)
    step_s04$enter(extra = paste0("storage_key=", storage_key))
    # 5a. @stats new location
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[storage_key]])) {
        scope_obj@stats[[storage_key]] <- list()
    }
    scope_obj@stats[[storage_key]]$iDeltaStats <- list(
        genes      = genes,
        delta_raw  = delta_raw, # Ensure this is a numeric vector
        level      = level
    )

    # 5b. @meta.data maintain old behavior (column name based on level)
    if (is.null(scope_obj@meta.data)) {
        scope_obj@meta.data <- data.frame(row.names = genes)
    }

    # Ensure all genes have rows in meta.data
    missing_genes <- setdiff(genes, rownames(scope_obj@meta.data))
    if (length(missing_genes) > 0) {
        new_rows <- data.frame(row.names = missing_genes)
        scope_obj@meta.data <- rbind(scope_obj@meta.data, new_rows)
    }

    # Assign scalar values (not matrix) to each gene
    scope_obj@meta.data[genes, col_suffix] <- delta_raw[genes]

    .log_info("computeIDelta", "S04", "I-delta computation completed", verbose)
    step_s04$done(extra = paste0("meta_col=", col_suffix))
    invisible(scope_obj)
}

#' Compute Morisita-Horn similarity between genes
#' @description
#' Internal helper for `.compute_mh`.
#' Builds Morisita-Horn similarity across genes using a consensus graph and
#' stores the result as matrix and/or igraph.
#' @param scope_obj A `scope_object` containing Lee statistics.
#' @param grid_name Grid layer name.
#' @param lee_layer Lee statistics layer name.
#' @param graph_slot Graph slot holding the consensus network.
#' @param graph_slot_mh Slot name to store the Morisita-Horn graph.
#' @param matrix_slot_mh Optional slot to store the Morisita-Horn matrix.
#' @param out Output type: matrix, igraph, or both.
#' @param L_min Minimum Lee's L threshold for edges.
#' @param area_norm Whether to normalise by area.
#' @param ncores Number of threads to use.
#' @param use_chao Apply Chao correction when TRUE.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @keywords internal
.compute_mh <- function(scope_obj,
                      grid_name = NULL,
                      lee_layer = "LeeStats_Xz",
                      graph_slot = "g_consensus",
                      graph_slot_mh = "g_morisita",
                      matrix_slot_mh = NULL,
                      out = c("matrix", "igraph", "both"),
                      L_min = 0,
                      area_norm = TRUE,
                      ncores = 8,
                      use_chao = TRUE,
                      verbose = TRUE) {
    out <- match.arg(out)
    step_s01 <- .log_step("computeMH", "S01", "select grid and Lee layer", verbose)
    step_s01$enter(extra = paste0(
        "grid_name=", if (is.null(grid_name)) "auto" else grid_name,
        " lee_layer=", lee_layer
    ))
    .log_info("computeMH", "S01",
        paste0(
            "L_min=", L_min,
            " area_norm=", area_norm,
            " use_chao=", use_chao,
            " ncores=", ncores,
            " out=", out
        ),
        verbose
    )

    ## ---- 1. Select and check grid layer ---- ##
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    gname <- names(scope_obj@grid)[
        vapply(scope_obj@grid, identical, logical(1), g_layer)
    ]
    .check_grid_content(scope_obj, gname)

    .log_info("computeMH", "S01", paste0("Selected grid layer: ", gname), verbose)

    ## ---- 2. Get Lee's L ---- ##
    Lmat <- .get_lee_matrix(scope_obj, grid_name = gname, lee_layer = lee_layer)

    .log_info("computeMH", "S01",
        paste0("Lee's L matrix dimensions: ", nrow(Lmat), "x", ncol(Lmat)),
        verbose
    )
    if (verbose) {
        lee_stats <- summary(as.vector(Lmat[upper.tri(Lmat)]))
        .log_info("computeMH", "S01",
            paste0(
                "Lee's L range: [",
                round(lee_stats[1], 4), ", ", round(lee_stats[6], 4), "]"
            ),
            verbose
        )
    }
    step_s01$done(extra = paste0("n_genes=", nrow(Lmat)))

    ## ---- 3. Try to read existing graph (only when graph_slot is not NULL) ---- ##
    step_s02 <- .log_step("computeMH", "S02", "resolve or build graph", verbose)
    step_s02$enter(extra = paste0("graph_slot=", if (is.null(graph_slot)) "NULL" else graph_slot))
    gnet <- NULL
    if (!is.null(graph_slot)) {
        # (1) First check @stats
        if (!is.null(scope_obj@stats) &&
            !is.null(scope_obj@stats[[gname]]) &&
            !is.null(scope_obj@stats[[gname]][[lee_layer]]) &&
            !is.null(scope_obj@stats[[gname]][[lee_layer]][[graph_slot]])) {
            gnet <- scope_obj@stats[[gname]][[lee_layer]][[graph_slot]]
            .log_backend("computeMH", "S02", "igraph", "reuse_stats", verbose = verbose)
            .log_info("computeMH", "S02", "Found existing graph in @stats", verbose)
        }
        # (2) Then check legacy @grid
        else if (!is.null(g_layer[[lee_layer]]) &&
            !is.null(g_layer[[lee_layer]][[graph_slot]])) {
            gnet <- g_layer[[lee_layer]][[graph_slot]]
            .log_backend("computeMH", "S02", "igraph", "reuse_grid_legacy", verbose = verbose)
            .log_info("computeMH", "S02", "Found existing graph in @grid (legacy)", verbose)
        }
    } else {
        .log_backend("computeMH", "S02", "igraph", "build_temporary", verbose = verbose)
        .log_info("computeMH", "S02", "graph_slot is NULL, will build temporary graph", verbose)
    }

    ## ---- 4. If still no network, build based on L_min ---- ##
    if (is.null(gnet)) {
        .log_backend("computeMH", "S02", "igraph",
            paste0("build_from_L L_min=", L_min),
            verbose = verbose
        )
        .log_info("computeMH", "S02", "Building new network from Lee's L matrix", verbose)
        keep <- (Lmat >= L_min)
        diag(keep) <- FALSE
        idx <- which(keep, arr.ind = TRUE)
        if (!nrow(idx)) {
            stop("`L_min` too high, no edges in Lee's L matrix meet criteria.")
        }

        if (verbose) {
            .log_info("computeMH", "S02",
                paste0("Edges meeting L_min threshold: ", nrow(idx)),
                verbose
            )
            .log_info("computeMH", "S02",
                paste0(
                    "Edge weight range: [",
                    round(min(Lmat[idx]), 4), ", ", round(max(Lmat[idx]), 4), "]"
                ),
                verbose
            )
        }

        edges_df <- data.frame(
            from = rownames(Lmat)[idx[, 1]],
            to = colnames(Lmat)[idx[, 2]],
            weight = Lmat[idx],
            stringsAsFactors = FALSE
        )
        gnet <- igraph::graph_from_data_frame(edges_df, directed = FALSE)
    }

    if (verbose) {
        n_vertices <- igraph::vcount(gnet)
        n_edges <- igraph::ecount(gnet)
        .log_info("computeMH", "S02",
            paste0("Final network: ", n_vertices, " vertices, ", n_edges, " edges"),
            verbose
        )
    }
    step_s02$done()

    ## ---- 5. gene x grid sparse count matrix (same as legacy) ---- ##
    step_s03 <- .log_step("computeMH", "S03", "build gene-by-grid matrix", verbose)
    step_s03$enter()
    counts_dt <- g_layer$counts[count > 0]
    genes <- sort(unique(counts_dt$gene))
    grids <- g_layer$grid_info$grid_id

    if (verbose) {
        .log_info("computeMH", "S03",
            paste0("Unique genes: ", length(genes)),
            verbose
        )
        .log_info("computeMH", "S03",
            paste0("Grid cells: ", length(grids)),
            verbose
        )
        .log_info("computeMH", "S03",
            paste0("Non-zero counts: ", nrow(counts_dt)),
            verbose
        )
    }

    Gsp <- sparseMatrix(
        i = match(counts_dt$gene, genes),
        j = match(counts_dt$grid_id, grids),
        x = counts_dt$count,
        dims = c(length(genes), length(grids)),
        dimnames = list(genes, grids)
    )

    .log_backend("computeMH", "S03", "R", "Matrix dgCMatrix sparse", verbose = verbose)
    if (verbose) {
        sparsity <- (1 - length(Gsp@x) / (nrow(Gsp) * ncol(Gsp))) * 100
        .log_info("computeMH", "S03",
            paste0(
                "Sparse matrix: ", nrow(Gsp), "x", ncol(Gsp),
                " (", round(sparsity, 2), "% sparse)"
            ),
            verbose
        )
    }
    step_s03$done(extra = paste0("n_genes=", length(genes)))

    ## ---- 5a. Optional area normalization ---- ##
    step_s04 <- .log_step("computeMH", "S04", "apply area normalization", verbose)
    step_s04$enter(extra = paste0("area_norm=", area_norm))
    if (area_norm) {
        .log_info("computeMH", "S04", "Applying area normalization to count matrix", verbose)
        gi <- g_layer$grid_info
        if (all(c("width", "height") %in% names(gi))) {
            gi <- gi[match(colnames(Gsp), gi$grid_id), ]
            area_vec <- gi$width * gi$height
            if (length(unique(area_vec)) == 1L) {
                .log_info("computeMH", "S04",
                    paste0("Uniform grid area: ", unique(area_vec)),
                    verbose
                )
                Gsp@x <- Gsp@x / unique(area_vec)
            } else {
                if (verbose) {
                    area_stats <- summary(area_vec)
                    .log_info("computeMH", "S04",
                        paste0(
                            "Variable grid areas: [",
                            round(area_stats[1], 3), ", ", round(area_stats[6], 3), "]"
                        ),
                        verbose
                    )
                }
                Gsp <- Gsp %*% Diagonal(x = 1 / area_vec)
            }
        } else {
            .log_info("computeMH", "S04",
                "Warning: grid_info lacks width/height; skipping area_norm",
                verbose
            )
        }
    } else {
        .log_info("computeMH", "S04", "Skipping area normalization (area_norm = FALSE)", verbose)
    }
    step_s04$done()

    ## ---- 6. Morisita-Horn (C++) ---- ##
    step_s05 <- .log_step("computeMH", "S05", "compute Morisita-Horn", verbose)
    step_s05$enter()
    ed <- t(igraph::as_edgelist(gnet, names = FALSE) - 1L)
    storage.mode(ed) <- "integer"

    if (verbose) {
        .log_info("computeMH", "S05",
            paste0("Processing ", ncol(ed), " edges"),
            verbose
        )
    }
    .log_backend("computeMH", "S05", "C++",
        paste0("morisita_horn_sparse threads=", ncores, " use_chao=", use_chao),
        verbose = verbose
    )

    mh_vec <- .morisita_horn_sparse(
        G = Gsp, edges = ed,
        use_chao = use_chao, nthreads = ncores
    )

    if (verbose) {
        # Check if mh_vec contains valid numeric values
        if (is.numeric(mh_vec) && length(mh_vec) > 0) {
            finite_vals <- mh_vec[is.finite(mh_vec)]
            if (length(finite_vals) > 0) {
                mh_min <- min(finite_vals, na.rm = TRUE)
                mh_max <- max(finite_vals, na.rm = TRUE)
                mh_mean <- mean(finite_vals, na.rm = TRUE)
                .log_info("computeMH", "S05",
                    paste0(
                        "Morisita-Horn similarity range: [",
                        round(mh_min, 4), ", ", round(mh_max, 4), "]"
                    ),
                    verbose
                )
                .log_info("computeMH", "S05",
                    paste0("Mean similarity: ", round(mh_mean, 4)),
                    verbose
                )
                if (length(finite_vals) < length(mh_vec)) {
                    n_invalid <- length(mh_vec) - length(finite_vals)
                    .log_info("computeMH", "S05",
                        paste0("Warning: ", n_invalid, " non-finite values found"),
                        verbose
                    )
                }
            } else {
                .log_info("computeMH", "S05",
                    "Warning: All similarity values are non-finite",
                    verbose
                )
            }
        } else {
            .log_info("computeMH", "S05",
                "Warning: Similarity computation returned non-numeric results",
                verbose
            )
        }
    }
    step_s05$done()

    if (out %in% c("igraph", "both")) {
        igraph::edge_attr(gnet, "CMH") <- mh_vec
    }

    ## ---- 7. Build sparse CMH matrix aligned to Lee's L gene set ---- ##
    step_s06 <- .log_step("computeMH", "S06", "build MH matrix", verbose)
    step_s06$enter()
    gene_all <- rownames(Lmat)
    if (is.null(gene_all)) {
        # Fallback to vertex names if Lmat lacks dimnames
        gene_all <- igraph::V(gnet)$name
    }
    ed_names <- igraph::as_edgelist(gnet, names = TRUE)
    ii <- match(ed_names[,1], gene_all)
    jj <- match(ed_names[,2], gene_all)
    ok <- !is.na(ii) & !is.na(jj)
    if (!all(ok)) {
        .log_info("computeMH", "S06",
            paste0("Note: ", sum(!ok), " edges not in L gene set; skipped in matrix build."),
            verbose
        )
    }
    ii <- as.integer(ii[ok]); jj <- as.integer(jj[ok]); xv <- as.numeric(mh_vec[ok])
    MH <- sparseMatrix(i = c(ii, jj), j = c(jj, ii), x = c(xv, xv),
                               dims = c(length(gene_all), length(gene_all)),
                               dimnames = list(gene_all, gene_all))
    MH <- drop0(as(MH, "dgCMatrix"))
    if (verbose) {
        .log_info("computeMH", "S06",
            paste0("CMH matrix built: ", nrow(MH), "x", ncol(MH), "; nnz=", nnzero(MH)),
            verbose
        )
    }
    step_s06$done(extra = paste0("nnz=", nnzero(MH)))

    ## ---- 8. Write to @stats ---- ##
    step_s07 <- .log_step("computeMH", "S07", "store outputs", verbose)
    step_s07$enter(extra = paste0("out=", out))
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()
    if (is.null(scope_obj@stats[[gname]][[lee_layer]])) scope_obj@stats[[gname]][[lee_layer]] <- list()

    # Determine matrix slot name if not provided
    if (is.null(matrix_slot_mh) || !nzchar(matrix_slot_mh)) {
        if (startsWith(graph_slot_mh, "g_")) {
            matrix_slot_mh <- sub("^g_", "m_", graph_slot_mh)
        } else {
            matrix_slot_mh <- paste0("m_", graph_slot_mh)
        }
    }

    if (out %in% c("igraph", "both")) {
        scope_obj@stats[[gname]][[lee_layer]][[graph_slot_mh]] <- gnet
    }
    if (out %in% c("matrix", "both")) {
        scope_obj@stats[[gname]][[lee_layer]][[matrix_slot_mh]] <- MH
    }

    .log_info("computeMH", "S07", "Analysis completed successfully", verbose)
    if (out %in% c("igraph","both")) {
        .log_info("computeMH", "S07",
            paste0("Graph stored in slot: ", graph_slot_mh, " (edge attr 'CMH')"),
            verbose
        )
    }
    if (out %in% c("matrix","both")) {
        .log_info("computeMH", "S07",
            paste0("Matrix stored in slot: ", matrix_slot_mh),
            verbose
        )
    }
    step_s07$done()

    invisible(scope_obj)
}

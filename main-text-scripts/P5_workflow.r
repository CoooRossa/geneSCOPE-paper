#!/usr/bin/env Rscript

library(geneSCOPE)
library(ggplot2)

gray_bg_theme <- ggplot2::theme(
  text = ggplot2::element_text(size = 8, face = "plain"),
  plot.background = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
  panel.background = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
  legend.background = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
  legend.box.background = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
  legend.key = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
  plot.title = ggplot2::element_text(size = 10, face = "plain"),
  axis.title = ggplot2::element_text(size = 9, face = "plain"),
  legend.title = ggplot2::element_text(size = 8, face = "plain"),
  legend.text = ggplot2::element_text(size = 8),
  strip.text = ggplot2::element_text(size = 9, face = "plain"),
  axis.text = ggplot2::element_text(size = 8)
)

P5.path <- "/path/to/GSE280314_Xenium_V1_Human_Colon_Cancer_P5_CRC_Add_on_FFPE_outs"
P5.coord_file <- "/path/to/P5_roi.csv"
grid_um <- 30
grid_name <- paste0("grid", grid_um)

P5.coord <- createSCOPE(
  data_dir = P5.path,
  grid_length = c(30),
  seg_type = "cell",
  coord_file = P5.coord_file,
  ncores = 96
)

P5.coord <- addSingleCells(
  scope_obj = P5.coord,
  xenium_dir = P5.path
)

P5.coord <- normalizeSingleCells(
  scope_obj = P5.coord,
  input_layer = "counts",
  output_layer = "logCPM",
  scale_factor = 1e4
)

P5.coord <- normalizeMoleculesInGrid(
  scope_obj = P5.coord,
  grid_name = grid_name
)

P5.coord <- computeWeights(
  scope_obj = P5.coord,
  grid_name = grid_name
)

P5.coord <- computeL(
  scope_obj = P5.coord,
  use_bigmemory = FALSE,
  grid_name = grid_name,
  ncores = 64
)

P5.coord <- computeCorrelation(
  scope_obj = P5.coord,
  level = "cell",
  layer = "logCPM",
  method = "pearson",
  blocksize = 2000,
  ncores = 64
)

curve_name <- paste0("LR_curve_", grid_um)
P5.coord <- computeLvsRCurve(
  scope_obj = P5.coord,
  level = "cell",
  grid_name = grid_name,
  ncores = 64,
  downsample = 0.05,
  k_max = 2000,
  n_strata = 1000,
  min_rel_width = 0.15,
  widen_span = 0.1,
  curve_name = curve_name
)

p_lvsr <- plotLvsR(
  scope_obj = P5.coord,
  grid_name = grid_name,
  pear_level = "cell",
  delta_top_n = 0,
  flip = TRUE
)

p_lvsr <- p_lvsr +
  ggplot2::geom_ribbon(
    data = P5.coord@stats[[grid_name]]$LeeStats_Xz[[curve_name]],
    ggplot2::aes(x = Pear, ymin = lo95, ymax = hi95),
    inherit.aes = FALSE,
    fill = "orange",
    alpha = 0.25
  ) +
  ggplot2::geom_line(
    data = P5.coord@stats[[grid_name]]$LeeStats_Xz[[curve_name]],
    ggplot2::aes(x = Pear, y = fit),
    inherit.aes = FALSE,
    colour = "firebrick",
    linewidth = 0.8
  ) +
  gray_bg_theme

dir.create("./LvsR", showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = file.path("./LvsR", paste0("LvsR_grid", grid_um, ".png")),
  plot = p_lvsr,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)

options(future.globals.maxSize = 500000 * 1024^2)

cluster_col <- paste0("q95_res0.1_grid", grid_um, "_log1p_freq0.95")
P5.coord <- clusterGenes(
  scope_obj = P5.coord,
  grid_name = grid_name,
  L_min = 0,
  algo = "leiden",
  resolution = 0.10,
  pct_min = "q95",
  cluster_name = cluster_col,
  graph_slot_name = cluster_col,
  use_log1p_weight = TRUE,
  use_consensus = TRUE,
  consensus_thr = 0.95,
  n_restart = 1000
)

P5.coord@meta.data[[cluster_col]] <- factor(
  P5.coord@meta.data[[cluster_col]],
  levels = as.character(sort(unique(na.omit(P5.coord@meta.data[[cluster_col]]))))
)

network_dir <- file.path(".", paste0("grid", grid_um), "network")
dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)

p_network <- plotNetwork(
  scope_obj = P5.coord,
  lee_stats_layer = "LeeStats_Xz",
  grid_name = grid_name,
  use_consensus_graph = TRUE,
  graph_slot_name = cluster_col,
  cluster_vec = cluster_col,
  show_sign = TRUE,
  drop_isolated = TRUE,
  neg_linetype = "dashed",
  vertex_size = 4,
  max.overlaps = 10,
  base_edge_mult = 3,
  label_cex = 4,
  layout_niter = 1000,
  seed = 1,
  hub_factor = 1.5,
  L_min = 0.11,
  L_min_neg = 0.11,
  title = " "
) +
  gray_bg_theme

p_dendro_network <- plotDendroNetwork(
  scope_obj = P5.coord,
  lee_stats_layer = "LeeStats_Xz",
  grid_name = grid_name,
  use_consensus_graph = TRUE,
  graph_slot_name = cluster_col,
  cluster_vec = cluster_col,
  IDelta_col_name = NULL,
  node_size = 4,
  edge_width = 3,
  label_size = 4,
  seed = 1,
  max.overlaps = 10,
  title = " ",
  tree_mode = "radial"
) +
  gray_bg_theme

ggsave(
  filename = file.path(network_dir, paste0("network_", cluster_col, ".png")),
  plot = p_network,
  width = 15,
  height = 15,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(network_dir, paste0("dendro_network_", cluster_col, ".png")),
  plot = p_dendro_network,
  width = 15,
  height = 15,
  units = "in",
  dpi = 600
)

P5.coord <- computeDensity(
  scope_obj = P5.coord,
  grid_name = grid_name,
  layer_name = "counts",
  normalize_method = "none",
  density_name = "CEACAM5",
  genes = "CEACAM5"
)

P5.coord <- computeDensity(
  scope_obj = P5.coord,
  grid_name = grid_name,
  layer_name = "counts",
  normalize_method = "none",
  density_name = "CEACAM6",
  genes = "CEACAM6"
)

P5.coord <- computeDensity(
  scope_obj = P5.coord,
  grid_name = grid_name,
  layer_name = "counts",
  normalize_method = "none",
  density_name = "ACTA2",
  genes = "ACTA2"
)

grid_dir <- file.path(".", paste0("grid", grid_um))
density_dir <- file.path(grid_dir, "density")
dir.create(density_dir, recursive = TRUE, showWarnings = FALSE)

p_density_CEACAM5_ACTA2 <- plotDensity(
  scope_obj = P5.coord,
  density1_name = "CEACAM5",
  density2_name = "ACTA2",
  max.cutoff1 = 0.3,
  max.cutoff2 = 0.3,
  seg_type = "cell",
  grid_name = grid_name,
  alpha_seg = 0.5,
  alpha1 = 1,
  alpha2 = 0.5,
  legend_digits = 3,
  bar_offset = 0.02,
  arrow_pt = 2,
  scale_text_size = 1.5
) +
  ggplot2::ggtitle(paste0("P5 Density Overlay CEACAM5 vs ACTA2\nGrid Size ", grid_um)) +
  gray_bg_theme

ggsave(
  filename = file.path(density_dir, paste0("CEACAM5_ACTA2_grid", grid_um, ".png")),
  plot = p_density_CEACAM5_ACTA2,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

p_density_CEACAM5_CEACAM6 <- plotDensity(
  scope_obj = P5.coord,
  density1_name = "CEACAM5",
  density2_name = "CEACAM6",
  max.cutoff1 = 0.3,
  max.cutoff2 = 0.3,
  seg_type = "cell",
  grid_name = grid_name,
  alpha_seg = 0.5,
  alpha1 = 1,
  alpha2 = 0.5,
  legend_digits = 3,
  bar_offset = 0.02,
  arrow_pt = 2,
  scale_text_size = 1.5
) +
  ggplot2::ggtitle(paste0("P5 Density Overlay CEACAM5 vs CEACAM6\nGrid Size ", grid_um)) +
  gray_bg_theme

ggsave(
  filename = file.path(density_dir, paste0("CEACAM5_CEACAM6_grid", grid_um, ".png")),
  plot = p_density_CEACAM5_CEACAM6,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

cells_density_dir <- file.path(".", "cells", "density")
dir.create(cells_density_dir, recursive = TRUE, showWarnings = FALSE)

p_centroids_CEACAM5_CEACAM6 <- plotDensityCentroids(
  scope_obj = P5.coord,
  gene1_name = "CEACAM5",
  gene2_name = "CEACAM6",
  seg_type = "cell",
  max.cutoff1 = 0.3,
  max.cutoff2 = 0.3,
  alpha_seg = 0.50,
  alpha1 = 1,
  alpha2 = 0.5,
  bar_offset = 0.02,
  arrow_pt = 2,
  scale_text_size = 1.5
) +
  ggplot2::ggtitle("P5 CEACAM5 vs CEACAM6 Centroids") +
  gray_bg_theme +
  ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 3))

ggsave(
  filename = file.path(cells_density_dir, "CEACAM5_CEACAM6_centroids_grid_cells.png"),
  plot = p_centroids_CEACAM5_CEACAM6,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

p_centroids_CEACAM5_ACTA2 <- plotDensityCentroids(
  scope_obj = P5.coord,
  gene1_name = "CEACAM5",
  gene2_name = "ACTA2",
  seg_type = "cell",
  max.cutoff1 = 0.3,
  max.cutoff2 = 0.3,
  alpha_seg = 0.50,
  alpha1 = 1,
  alpha2 = 0.5,
  bar_offset = 0.02,
  arrow_pt = 2,
  scale_text_size = 1.5
) +
  ggplot2::ggtitle("P5 CEACAM5 vs ACTA2 Centroids") +
  gray_bg_theme +
  ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 3))

ggsave(
  filename = file.path(cells_density_dir, "CEACAM5_ACTA2_centroids_grid_cells.png"),
  plot = p_centroids_CEACAM5_ACTA2,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

p_grid_boundary <- plotGridBoundary(
  scope_obj = P5.coord,
  grid_name = grid_name
) +
  gray_bg_theme

ggsave(
  filename = file.path(grid_dir, paste0("grid", grid_um, "_boundary.png")),
  plot = p_grid_boundary,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

# ---- Density plots for selected clusters ----
cluster_ids <- c(1:24)
if (!cluster_col %in% colnames(P5.coord@meta.data)) {
  stop("Cluster column not found in meta.data: ", cluster_col)
}

cluster_genes <- rownames(P5.coord@meta.data)[P5.coord@meta.data[[cluster_col]] %in% cluster_ids]
cluster_genes <- cluster_genes[order(P5.coord@meta.data[cluster_genes, cluster_col], cluster_genes)]

cluster_density_dir <- file.path(grid_dir, "density", paste0("clusters_", paste(cluster_ids, collapse = "_")))
dir.create(cluster_density_dir, recursive = TRUE, showWarnings = FALSE)

for (gene in cluster_genes) {
  P5.coord <- computeDensity(
    scope_obj = P5.coord,
    grid_name = grid_name,
    layer_name = "counts",
    normalize_method = "none",
    density_name = gene,
    genes = gene
  )

  p_gene_density <- plotDensity(
    scope_obj = P5.coord,
    density1_name = gene,
    max.cutoff1 = 0.6,
    seg_type = "cell",
    grid_name = grid_name,
    alpha_seg = 0.15,
    alpha1 = 1,
    legend_digits = 3,
    bar_offset = 0.02,
    arrow_pt = 2,
    scale_text_size = 1.5
  ) +
    ggplot2::ggtitle(paste0("P5 Density Overlay for ", gene, " (Cluster ", P5.coord@meta.data[gene, cluster_col], ")")) +
    gray_bg_theme

  ggsave(
    filename = file.path(cluster_density_dir, paste0(gene, "_density.png")),
    plot = p_gene_density,
    width = 5,
    height = 5,
    units = "in",
    dpi = 600
  )
}

top.delta.l <- getTopLvsR(
  scope_obj = P5.coord,
  grid_name = grid_name,
  pear_level = "cell",
  L_range = c(0.0, 1),
  top_n = 40000,
  ncores = 64,
  direction = "largest",
  pval_mode = "uniform",
  curve_layer = curve_name,
  CI_rule = "remove_within"
)
top.delta.l <- top.delta.l[top.delta.l$fdr < 0.05 & top.delta.l$pct1 >= 5 & top.delta.l$pct2 >= 5, , drop = FALSE]

if (nrow(top.delta.l) < 20) {
  stop("top.delta.l has fewer than 20 rows after filtering (n = ", nrow(top.delta.l), ").")
}

topdelta_dir <- file.path(".", "TopDelta")
dir.create(topdelta_dir, recursive = TRUE, showWarnings = FALSE)

for (j in 1:20) {
  gene1 <- top.delta.l$gene1[j]
  gene2 <- top.delta.l$gene2[j]

  P5.coord <- computeDensity(
    scope_obj = P5.coord,
    grid_name = grid_name,
    layer_name = "counts",
    normalize_method = "none",
    density_name = gene1,
    genes = gene1
  )
  P5.coord <- computeDensity(
    scope_obj = P5.coord,
    grid_name = grid_name,
    layer_name = "counts",
    normalize_method = "none",
    density_name = gene2,
    genes = gene2
  )

  p_td <- plotDensity(
    scope_obj = P5.coord,
    density1_name = gene1,
    density2_name = gene2,
    max.cutoff1 = 0.25,
    max.cutoff2 = 0.25,
    seg_type = "cell",
    grid_name = grid_name,
    alpha_seg = 0.5,
    alpha1 = 1,
    alpha2 = 0.5,
    legend_digits = 3,
    bar_offset = 0.02,
    arrow_pt = 2,
    scale_text_size = 1.5
  ) +
    ggplot2::ggtitle(paste0("P5 Density Overlay ", gene1, " vs ", gene2, " - Grid Size ", grid_um)) +
    gray_bg_theme

  ggsave(
    filename = file.path(topdelta_dir, paste0(gene1, "_", gene2, "_grid", grid_um, "_L_0.1.1.png")),
    plot = p_td,
    width = 5,
    height = 5,
    units = "in",
    dpi = 600
  )
}

layer_name <- grid_name
mat <- P5.coord@stats[[layer_name]]$LeeStats_Xz$L

heatmap_dir <- file.path(grid_dir, "heatmap")
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(heatmap_dir, paste0("LeeL_heatmap_grid", grid_um, "_", cluster_col, ".png"))
png(
  filename = out_file,
  width = 20,
  height = 20,
  units = "in",
  res = 600,
  bg = "#c0c0c0"
)

ht <- ComplexHeatmap::Heatmap(
  mat,
  name = "Lee's L",
  col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "#c0c0c0", "red")),
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_side = "left",
  column_dend_side = "top",
  row_title = NULL,
  column_title = NULL
)
ComplexHeatmap::draw(ht, background = "#c0c0c0")
ComplexHeatmap::decorate_heatmap_body(
  "Lee's L",
  grid::grid.rect(gp = grid::gpar(fill = NA, col = "#1a1a1a", lwd = 2))
)
dev.off()
message("Saved: ", out_file)

# ---- I_delta by cluster ----
P5.coord <- computeIDelta(
  scope_obj = P5.coord,
  grid_name = grid_name,
  level = "grid",
  ncores = 64
)

idelta_dir <- file.path(grid_dir, "idelta")
dir.create(idelta_dir, recursive = TRUE, showWarnings = FALSE)

p_idelta <- plotIDelta(
  scope_obj = P5.coord,
  grid_name = grid_name,
  cluster_col = cluster_col,
  nrow = 2
) +
  ggplot2::theme(legend.position = "none")

ggsave(
  filename = file.path(idelta_dir, paste0(cluster_col, "_idelta_by_cluster.png")),
  plot = p_idelta,
  width = 15,
  height = 7,
  units = "in",
  dpi = 600
)

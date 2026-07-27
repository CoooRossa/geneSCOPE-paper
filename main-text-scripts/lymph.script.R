#!/usr/bin/env Rscript

library(geneSCOPE)
library(ggplot2)

script_dir <- local({
  x <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(x)) dirname(normalizePath(sub("^--file=", "", x[[1L]]))) else normalizePath(getwd())
})
source(file.path(script_dir, "freeze_helpers.R"))

if (!identical(as.character(utils::packageVersion("geneSCOPE")), "1.0.2")) {
  stop("This workflow requires geneSCOPE 1.0.2.")
}
freeze_source <- require_freeze_source_metadata()
gate_max_abs_L_diff <- canonical_lee_s2_gate()
seed <- integer_env("GENESCOPE_SEED", 1L)
ncores <- integer_env("GENESCOPE_THREADS", 64L)
configure_freeze_runtime()

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

Lymph.path <- required_directory_env("GENESCOPE_LN_OUTS")
Lymph.scope_file <- required_file_env("GENESCOPE_LN_SCOPE_RDS")
Lymph.top_pairs_file <- required_file_env("GENESCOPE_LN_TOP_PAIRS")
Lymph.analysis_manifest_file <- normalizePath(
  file.path(dirname(Lymph.scope_file), "manifest.json"), mustWork = TRUE
)
Lymph.delta_manifest_file <- normalizePath(
  file.path(dirname(Lymph.top_pairs_file), "manifest.json"), mustWork = TRUE
)
Lymph.generator_file <- normalizePath(
  file.path(dirname(Lymph.top_pairs_file), "recompute_ln_complete_delta.R"),
  mustWork = TRUE
)
Lymph.coord_file <- normalizePath(
  Sys.getenv("GENESCOPE_LN_ROI", file.path(script_dir, "..", "ROI-coordinate-files", "lymph_roi.csv")),
  mustWork = TRUE
)
output_root <- Sys.getenv("GENESCOPE_LN_OUTPUT", file.path(getwd(), "LN_correction_output"))
output_root <- assert_fresh_output_dir(output_root)
setwd(output_root)

grid_um <- 30
grid_name <- paste0("grid", grid_um)
raw_input_gate <- assert_reference_raw_inputs(
  script_dir, "LN", Lymph.path, Lymph.coord_file
)
curve_name <- paste0("LR_curve_", grid_um, "_shuffle")
authoritative_cluster_col <- "shuffle_q99.9_res0.1_grid30"
analysis_sources <- list(
  scope = assert_reference_artifact(
    script_dir, "MAIN_RESULTS", "LN/LN_scope_shuffle_v102.rds",
    Lymph.scope_file
  ),
  top_pairs = assert_reference_artifact(
    script_dir, "LN_COMPLETE_RESULTS", "LN_top_pairs_complete_delta_v102.tsv",
    Lymph.top_pairs_file
  ),
  analysis_manifest = assert_reference_artifact(
    script_dir, "MAIN_RESULTS", "LN/manifest.json", Lymph.analysis_manifest_file
  ),
  delta_manifest = assert_reference_artifact(
    script_dir, "LN_COMPLETE_RESULTS", "manifest.json", Lymph.delta_manifest_file
  ),
  generator = assert_reference_artifact(
    script_dir, "LN_COMPLETE_RESULTS", "recompute_ln_complete_delta.R",
    Lymph.generator_file
  )
)
analysis_provenance_gate <- assert_analysis_provenance(
  "LN", Lymph.analysis_manifest_file, top_pair_rows = 75405L,
  delta_manifest_path = Lymph.delta_manifest_file,
  generator_path = Lymph.generator_file
)
Lymph.coord <- readRDS(Lymph.scope_file)
analysis_scope_gate <- assert_authoritative_scope(
  Lymph.coord, "LN", grid_name, curve_name, authoritative_cluster_col
)
raw_identity_gate <- assert_scope_xenium_identity(
  Lymph.coord, Lymph.path, Lymph.coord_file, "LN"
)

p_lvsr <- plotLvsR(
  scope_obj = Lymph.coord,
  grid_name = grid_name,
  pear_level = "cell",
  delta_top_n = 0,
  flip = TRUE
)

p_lvsr <- p_lvsr +
  ggplot2::geom_ribbon(
    data = Lymph.coord@stats[[grid_name]]$LeeStats_Xz[[curve_name]],
    ggplot2::aes(x = Pear, ymin = lo95, ymax = hi95),
    inherit.aes = FALSE,
    fill = "orange",
    alpha = 0.25
  ) +
  ggplot2::geom_line(
    data = Lymph.coord@stats[[grid_name]]$LeeStats_Xz[[curve_name]],
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

pct_mins <- "q99.9"
cluster_cols <- paste0(pct_mins, "_res0.1_grid", grid_um, "_log1p_freq0.95")

network_dir <- file.path(".", paste0("grid", grid_um), "network")
dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)

for (idx in seq_along(pct_mins)) {
  pct_min <- pct_mins[[idx]]
  cluster_col <- cluster_cols[[idx]]

  display_mapping_path <- file.path(
    script_dir, "..", "correction-analysis", "display-mappings", "LN_display_mapping.tsv"
  )
  display_mapping <- read_display_mapping(script_dir, "LN")
  Lymph.coord@meta.data[[paste0(cluster_col, "_raw")]] <- as.character(
    Lymph.coord@meta.data[[authoritative_cluster_col]]
  )
  membership_gate <- assert_reference_membership(
    script_dir, "LN", rownames(Lymph.coord@meta.data),
    Lymph.coord@meta.data[[paste0(cluster_col, "_raw")]]
  )
  Lymph.coord@meta.data[[cluster_col]] <- apply_display_mapping(
    Lymph.coord@meta.data[[paste0(cluster_col, "_raw")]], display_mapping
  )
  cluster_palette <- display_palette(display_mapping)
  Lymph.coord@meta.data[[cluster_col]] <- factor(
    Lymph.coord@meta.data[[cluster_col]],
    levels = as.character(sort(unique(na.omit(Lymph.coord@meta.data[[cluster_col]]))))
  )

  p_network <- plotNetwork(
    scope_obj = Lymph.coord,
    lee_stats_layer = "LeeStats_Xz",
    grid_name = grid_name,
    use_consensus_graph = TRUE,
    graph_slot_name = authoritative_cluster_col,
    cluster_vec = cluster_col,
    cluster_palette = cluster_palette,
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

  dendro_out <- plotDendroNetwork(
    scope_obj = Lymph.coord,
    lee_stats_layer = "LeeStats_Xz",
    grid_name = grid_name,
    use_consensus_graph = TRUE,
    graph_slot_name = authoritative_cluster_col,
    cluster_vec = cluster_col,
    cluster_palette = cluster_palette,
    IDelta_col_name = NULL,
    node_size = 4,
    edge_width = 3,
    label_size = 4,
    seed = 1,
    max.overlaps = 10,
    title = " ",
    tree_mode = "radial"
  )

  dendro_plot <- if (inherits(dendro_out, "ggplot")) {
    dendro_out
  } else if (is.list(dendro_out) && inherits(dendro_out$plot, "ggplot")) {
    dendro_out$plot
  } else {
    stop("plotDendroNetwork returned an unexpected type.")
  }
  dendro_plot <- dendro_plot + gray_bg_theme

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
    plot = dendro_plot,
    width = 15,
    height = 15,
    units = "in",
    dpi = 600
  )
}

top.delta.all <- read_authoritative_top_pairs(
  Lymph.top_pairs_file, "LN", expected_rows = 75405L
)
pair_scope_gate <- assert_scope_pair_table(Lymph.coord, top.delta.all, "LN")
top.delta.l <- filter_display_pairs(top.delta.all)
assert_reference_top6(
  file.path(script_dir, "..", "correction-analysis"), "LN", top.delta.l
)

if (!file.copy(Lymph.top_pairs_file, "LN_top_pairs_all.tsv", overwrite = FALSE)) {
  stop("Could not copy the authoritative LN pair table into the figure bundle.")
}
utils::write.table(top.delta.l, "LN_top_pairs_display_filter.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(utils::head(top.delta.l, 6L), "LN_Top6.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

density_genes <- c("ITGB2", "PDGFRA", "PTPN6")
for (gene in density_genes) {
  Lymph.coord <- computeDensity(
    scope_obj = Lymph.coord,
    grid_name = grid_name,
    layer_name = "counts",
    normalize_method = "none",
    density_name = gene,
    genes = gene
  )
}

grid_dir <- file.path(".", paste0("grid", grid_um))
density_dir <- file.path(grid_dir, "density")
dir.create(density_dir, recursive = TRUE, showWarnings = FALSE)

p_density_ITGB2_PDGFRA <- plotDensity(
  scope_obj = Lymph.coord,
  density1_name = "ITGB2",
  density2_name = "PDGFRA",
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
  ggplot2::ggtitle(paste0("Lymph Density Overlay ITGB2 vs PDGFRA\nGrid Size ", grid_um)) +
  gray_bg_theme

ggsave(
  filename = file.path(density_dir, paste0("ITGB2_PDGFRA_grid", grid_um, ".png")),
  plot = p_density_ITGB2_PDGFRA,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

p_density_ITGB2_PTPN6 <- plotDensity(
  scope_obj = Lymph.coord,
  density1_name = "ITGB2",
  density2_name = "PTPN6",
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
  ggplot2::ggtitle(paste0("Lymph Density Overlay ITGB2 vs PTPN6\nGrid Size ", grid_um)) +
  gray_bg_theme

ggsave(
  filename = file.path(density_dir, paste0("ITGB2_PTPN6_grid", grid_um, ".png")),
  plot = p_density_ITGB2_PTPN6,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

cells_density_dir <- file.path(".", "cells", "density")
dir.create(cells_density_dir, recursive = TRUE, showWarnings = FALSE)

p_centroids_ITGB2_PDGFRA <- plotDensityCentroids(
  scope_obj = Lymph.coord,
  gene1_name = "ITGB2",
  gene2_name = "PDGFRA",
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
  ggplot2::ggtitle("Lymph ITGB2 vs PDGFRA Centroids") +
  gray_bg_theme +
  ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 3))

ggsave(
  filename = file.path(cells_density_dir, "ITGB2_PDGFRA_centroids_grid_cells.png"),
  plot = p_centroids_ITGB2_PDGFRA,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

p_centroids_ITGB2_PTPN6 <- plotDensityCentroids(
  scope_obj = Lymph.coord,
  gene1_name = "ITGB2",
  gene2_name = "PTPN6",
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
  ggplot2::ggtitle("Lymph ITGB2 vs PTPN6 Centroids") +
  gray_bg_theme +
  ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 3))

ggsave(
  filename = file.path(cells_density_dir, "ITGB2_PTPN6_centroids_grid_cells.png"),
  plot = p_centroids_ITGB2_PTPN6,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600
)

p_grid_boundary <- plotGridBoundary(
  scope_obj = Lymph.coord,
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

output_gate <- assert_required_figure_outputs(
  output_root, "LN",
  required_relative = c(
    "LvsR/LvsR_grid30.png",
    file.path(network_dir, paste0("network_", cluster_col, ".png")),
    file.path(network_dir, paste0("dendro_network_", cluster_col, ".png")),
    file.path(grid_dir, paste0("grid", grid_um, "_boundary.png")),
    "LN_top_pairs_all.tsv", "LN_top_pairs_display_filter.tsv", "LN_Top6.tsv"
  ),
  expected_png_count = 8L
)

write_freeze_output_manifest(
  output_root, "LN", freeze_source, gate_max_abs_L_diff, membership_gate,
  input_dir = Lymph.path,
  roi_file = Lymph.coord_file,
  display_mapping_path = display_mapping_path,
  workflow_path = file.path(script_dir, "lymph.script.R"),
  analysis_sources = analysis_sources,
  parameters = list(
    grid_um = grid_um, seed = seed, ncores = ncores,
    analysis_mode = "render_from_hash_pinned_authoritative_results",
    display_derivations = "computeDensity_only",
    authoritative_scope_gate = analysis_scope_gate,
    raw_identity_gate = raw_identity_gate,
    analysis_provenance_gate = analysis_provenance_gate,
    pair_scope_gate = pair_scope_gate,
    L_permutations = 1000L, delta_permutations = 1000L,
    delta_adjustment = "BH_universe",
    display_filter = list(q_Delta = "<0.05", L = ">0", r = "<0.05",
                          pct1 = ">20", pct2 = ">20")
  ),
  output_gate = output_gate
)

library(geneSCOPE)

# Paths (replace with your local paths)
scope.path <- "path/to/xenium_scope_outs"
scope.coord_file <- "path/to/scope_tumor_region_coordinates.csv"
scope.idelta_rds <- "path/to/scope_coord_with_iDelta_all_grids.rds"

figures_dir <- "path/to/output/figures"
figures_cosmx_dir <- "path/to/output/figures-cosmx"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_cosmx_dir, recursive = TRUE, showWarnings = FALSE)

scope.coord <- createSCOPE(
  data_dir = scope.path,
  grid_length = seq(5, 150, 1),
  seg_type = "cell",
 filtergenes = T,
  max_dist_mol_nuc = 25,
  filterqv = 20,
  max_gene_types = 8000,
  min_gene_types = 1,
  min_seg_points = 1,
 coord_file = scope.coord_file,
  ncores = 128,
  verbose = TRUE,
  flip_y = T
)

for (i in seq(5, 120, 1)) {
  scope.coord <- computeIDelta(
    scope_obj = scope.coord,
    grid_name = paste0("grid", i)
  )
}
# saveRDS(scope.coord, file = scope.idelta_rds)
scope.coord <- readRDS(scope.idelta_rds)
scope.coord <- addSingleCells(
  scope_obj = scope.coord,
  xenium_dir = scope.path
)

scope.coord <- computeIDelta(
  scope_obj = scope.coord,
  level = "cell"
)

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(inflection)


# ---- Stable knee detection helpers ----
get_smooth_xy <- function(df, span = 0.3, by = 1, method = c("loess", "approx")) {
  method <- match.arg(method)
  df <- df %>% arrange(grid)
  if (n_distinct(df$grid) < 6 || all(is.na(df$iDelta))) {
    return(tibble(grid = numeric(0), iDelta = numeric(0)))
  }
  xx <- seq(min(df$grid, na.rm = TRUE), max(df$grid, na.rm = TRUE), by = by)
  if (length(xx) < 2) {
    return(tibble(grid = df$grid, iDelta = df$iDelta))
  }
  if (method == "loess") {
    fit <- try(loess(iDelta ~ grid, data = df, span = span), silent = TRUE)
    if (inherits(fit, "try-error")) {
      yy <- approx(df$grid, df$iDelta, xout = xx, rule = 2)$y
      return(tibble(grid = xx, iDelta = yy))
    } else {
      preds <- predict(fit, newdata = data.frame(grid = xx))
      preds <- pmax(preds, 0)
      return(tibble(grid = xx, iDelta = preds))
    }
  } else {
    yy <- approx(df$grid, df$iDelta, xout = xx, rule = 2)$y
    yy <- pmax(yy, 0)
    return(tibble(grid = xx, iDelta = yy))
  }
}

safe_uik_idx <- function(x, y) {
  idx <- NA_integer_
  try({
    res <- inflection::uik(x, y)
    if (!is.null(res) && length(res) > 0 && is.finite(res[1])) idx <- res[1]
  }, silent = TRUE)
  idx
}

curvature_idx <- function(y) {
  n <- length(y)
  if (n < 5 || all(is.na(y))) return(NA_integer_)
  y1 <- c(NA, (y[3:n] - y[1:(n - 2)]) / 2, NA)
  y2 <- c(NA, y[3:n] - 2 * y[2:(n - 1)] + y[1:(n - 2)], NA)
  curv <- abs(y2) / (1 + y1^2)^(3/2)
  curv[!is.finite(curv)] <- NA
  k <- which.max(curv[2:(n - 1)])
  if (length(k) == 0 || is.infinite(k) || is.na(k)) return(NA_integer_)
  k + 1
}

detect_knee_general <- function(x, y, method = c("uik", "curvature"), return_index = FALSE) {
  method <- match.arg(method)
  ord <- order(x)
  x_sorted <- x[ord]
  y_sorted <- y[ord]
  df <- tibble(x = x_sorted, y = y_sorted) %>%
    group_by(x) %>% summarise(y = mean(y, na.rm = TRUE), .groups = "drop")
  x_sorted <- df$x; y_sorted <- df$y
  if (length(x_sorted) < 4 || all(is.na(y_sorted))) {
    return(if (return_index) NA_integer_ else NA_real_)
  }
  idx <- if (method == "uik") safe_uik_idx(x_sorted, y_sorted) else NA_integer_
  if (is.na(idx)) idx <- curvature_idx(y_sorted)
  if (is.na(idx)) return(if (return_index) NA_integer_ else NA_real_)
  if (return_index) return(idx)
  x_sorted[idx]
}

compute_knee_stable <- function(
  delta_long,
  do_outlier_filter = TRUE,
  q_lower = 0.1, q_upper = 0.9,
  smooth_method = c("loess", "approx"),
  loess_span = 0.3,
  resample_by = 1,
  log_x = FALSE,
  knee_method = c("uik", "curvature"),
  min_points = 6,
  verbose = TRUE
) {
  smooth_method <- match.arg(smooth_method)
  knee_method <- match.arg(knee_method)

  dl <- delta_long %>%
    select(gene, grid, iDelta) %>%
    mutate(grid = as.numeric(grid), iDelta = pmax(as.numeric(iDelta), 0)) %>%
    filter(!is.na(gene) & !is.na(grid) & !is.na(iDelta))

  if (do_outlier_filter) {
    iqr_tbl <- dl %>% group_by(gene) %>% summarise(
      ql = quantile(iDelta, q_lower, na.rm = TRUE),
      qu = quantile(iDelta, q_upper, na.rm = TRUE),
      iqr = qu - ql,
      .groups = "drop"
    ) %>% mutate(lower = ql - 1.5 * iqr, upper = qu + 1.5 * iqr)
    dl <- dl %>% inner_join(iqr_tbl, by = "gene") %>%
      filter(iDelta >= lower & iDelta <= upper) %>%
      select(gene, grid, iDelta)
  }

  dl_avg <- dl %>% group_by(gene, grid) %>%
    summarise(iDelta = mean(iDelta, na.rm = TRUE), .groups = "drop")

  delta_smooth <- dl_avg %>% group_by(gene) %>%
    group_modify(~ get_smooth_xy(.x, span = loess_span, by = resample_by, method = smooth_method)) %>%
    ungroup() %>% filter(!is.na(iDelta))

  gene_knee_results <- delta_smooth %>% group_by(gene) %>%
    group_modify(~ {
      g <- .x$grid; y <- .x$iDelta
      if (n_distinct(g) < min_points) return(tibble(uik_knee = NA_real_))
      xk <- if (log_x) log(g) else g
      knee_x <- detect_knee_general(xk, y, method = knee_method, return_index = FALSE)
      tibble(uik_knee = knee_x)
    }) %>% ungroup()

  valid_knees <- gene_knee_results %>% filter(!is.na(uik_knee))
  individual_knee_mean <- mean(valid_knees$uik_knee, na.rm = TRUE)
  individual_knee_sd   <- sd(valid_knees$uik_knee, na.rm = TRUE)

  avg_curve_mean <- delta_smooth %>% group_by(grid) %>%
    summarise(iDelta_mean = mean(iDelta, na.rm = TRUE), .groups = "drop") %>%
    arrange(grid)
  x_overall <- if (log_x) log(avg_curve_mean$grid) else avg_curve_mean$grid
  idx_overall <- detect_knee_general(x_overall, avg_curve_mean$iDelta_mean, method = knee_method, return_index = TRUE)
  overall_knee_mean <- if (is.na(idx_overall)) NA_real_ else avg_curve_mean$grid[idx_overall]

  gene_knee_with_y <- delta_smooth %>% inner_join(valid_knees, by = "gene") %>%
    group_by(gene) %>%
    slice_min(abs(grid - uik_knee), with_ties = FALSE) %>%
    ungroup() %>% select(gene, uik_knee, knee_y = iDelta)

  list(
    delta_smooth = delta_smooth,
    gene_knee_results = gene_knee_results,
    valid_knees = valid_knees,
    individual_knee_mean = individual_knee_mean,
    individual_knee_sd = individual_knee_sd,
    avg_curve_mean = avg_curve_mean,
    overall_knee_mean = overall_knee_mean,
    gene_knee_with_y = gene_knee_with_y
  )
}

meta <- scope.coord@meta.data %>%
  tibble::rownames_to_column(var = "gene")
meta <- meta[meta$gene %in% unique(scope.coord@grid$grid5$counts$gene), ]

cell_idelta <- meta %>%
  select(gene, cell_iDelta) %>%
  filter(!is.na(cell_iDelta))

grid_cols <- grep("^grid.*_iDelta([.][0-9]+)?$", colnames(meta), value = TRUE)[1:116] 

delta_long_tmp <- meta %>%
  select(gene, all_of(grid_cols)) %>%
  pivot_longer(
    -gene,
    names_to = c("grid", "rep"),
    names_pattern = "^grid(?:_lenGrid)?([0-9]+(?:[.][0-9]+)?)_iDelta(?:[.]([0-9]+))?$",
    values_to = "iDelta",
    values_drop_na = FALSE
  ) %>%
  mutate(grid = as.numeric(grid), rep = dplyr::coalesce(suppressWarnings(as.integer(rep)), 0L)) %>%
  arrange(gene, desc(grid), rep) %>%
  group_by(gene, grid) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  select(gene, grid, iDelta)

iqr_tbl <- delta_long_tmp %>%
  group_by(gene) %>%
  summarise(
    q10  = quantile(iDelta, 0.05, na.rm = TRUE),
    q90  = quantile(iDelta, 0.95, na.rm = TRUE),
    iqr  = q90 - q10,
    .groups = "drop"
  ) %>%
  mutate(
    lower = q10 - 1.5 * iqr,
    upper = q90 + 1.5 * iqr
  )

delta_long <- delta_long_tmp %>%
  inner_join(iqr_tbl, by = "gene") %>%
  filter(iDelta >= lower & iDelta <= upper) %>%
  select(gene, grid, iDelta) %>%
  arrange(gene, desc(grid))

# delta_long <- delta_only %>%
#   pivot_longer(-gene, names_to = "grid", values_to = "iDelta") %>%
#   mutate(grid = as.integer(grid)) %>%
#   arrange(gene, desc(grid))

delta_long[is.na(delta_long$grid), "grid"] <- -1

res <- compute_knee_stable(
  delta_long,
  do_outlier_filter = FALSE,
  smooth_method = "loess",
  loess_span = 0.1,
  resample_by = 1,
  log_x = FALSE,
  knee_method = "uik"
)

gene_knee_results    <- res$gene_knee_results
valid_knees          <- res$valid_knees
individual_knee_mean <- res$individual_knee_mean
individual_knee_sd   <- res$individual_knee_sd
avg_curve_mean       <- res$avg_curve_mean
overall_knee_mean    <- res$overall_knee_mean
gene_knee_with_y     <- res$gene_knee_with_y

cat("\n=== Individual gene UIK knee summary (stable) ===\n")
cat("Mean knee:", round(individual_knee_mean, 1), "μm\n")
cat("Standard deviation:", round(individual_knee_sd, 1), "μm\n")
if (nrow(valid_knees) > 0) {
  cat("Range:", round(min(valid_knees$uik_knee), 1), "-", round(max(valid_knees$uik_knee), 1), "μm\n")
} else {
  cat("Range: NA\n")
}

cat("\n=== Overall curve UIK knee (stable) ===\n")
cat("UIK knee of the mean curve:", round(overall_knee_mean, 1), "μm\n")

cat("\n=== UIK method comparison ===\n")
cat("Mean of per-gene knees:", round(individual_knee_mean, 1), "μm\n")
cat("Knee of the mean curve:", round(overall_knee_mean, 1), "μm\n")
cat("Difference:", round(abs(individual_knee_mean - overall_knee_mean), 1), "μm\n")

red_color <- "#E41A1C"
blue_color <- "#377EB8"

knee_range <- if (nrow(valid_knees) > 0) {
  data.frame(x = c(min(valid_knees$uik_knee), max(valid_knees$uik_knee)), ymin = 0, ymax = 100)
} else {
  data.frame(x = range(avg_curve_mean$grid, na.rm = TRUE), ymin = 0, ymax = 100)
}

p_main <- ggplot(avg_curve_mean, aes(x = grid, y = iDelta_mean)) +
  geom_ribbon(
    data = knee_range,
    aes(x = x, ymin = ymin - 5, ymax = ymax + 5),
    fill = blue_color, alpha = 0.5,
    linewidth = 0, color = blue_color,
    inherit.aes = FALSE
  ) +
  geom_line(color = "black", linewidth = 1) +
  {
    if (!is.na(overall_knee_mean)) {
      geom_vline(
        xintercept = overall_knee_mean,
        color = red_color, linewidth = 1, linetype = "solid"
      )
    }
  } + 
  scale_x_reverse(
    breaks = seq(
      min(avg_curve_mean$grid, na.rm = TRUE),
      max(avg_curve_mean$grid, na.rm = TRUE),
      by = 5
    )
  ) +
  labs(
    title = "scope - Knee Point Analysis",
    x = "Grid size (μm)",
    y = "Mean Curve Iδ"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_rect(fill = "#C0C0C0", color = NA),
    plot.background = element_rect(fill = "#ffffff", color = NA),
    axis.title = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 10),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  scale_y_continuous(
    limits = c(0, min(100, max(avg_curve_mean$iDelta_mean, na.rm = TRUE))),
    oob = scales::squish,
    breaks = scales::breaks_width(5)
  )
p_main 

ggsave(
  filename = file.path(figures_cosmx_dir, "scope_iDelta_mean_curve_knee_stable.png"),
  plot = p_main,
  width = 6,
  height = 5,
  dpi = 600
)
cell_curve_points <- res$delta_smooth %>%
inner_join(cell_idelta, by = "gene") %>%
group_by(gene) %>%
mutate(diff = abs(iDelta - cell_iDelta)) %>%
slice_min(diff, n = 1, with_ties = FALSE) %>%
ungroup()

vline_df <- data.frame(type = "Knee of Mean Curve", x = overall_knee_mean)

p_all_genes_annotated <- ggplot(res$delta_smooth, aes(x = grid, y = iDelta, group = gene)) +
geom_line(alpha = 0.5, linewidth = 0.5, color = "gray60") +
geom_ribbon(
data = knee_range,
aes(x = x, ymin = ymin - 5, ymax = ymax + 5),
fill = blue_color, alpha = 0.3,
linewidth = 0, color = blue_color,
inherit.aes = FALSE
) +
geom_line(
data = avg_curve_mean,
aes(x = grid, y = iDelta_mean, group = 1),
color = "black", linewidth = 1
) +
geom_point(
data = gene_knee_with_y,
    aes(x = uik_knee, y = knee_y, color = "Gene Wise Knee"),
size = 1, alpha = 0.5, show.legend = TRUE
) +
{
if (!is.na(overall_knee_mean)) {
geom_vline(
data = vline_df,
aes(xintercept = x, color = type),
linetype = "solid", linewidth = 1, show.legend = TRUE
)
}
} +
  scale_color_manual(
    name = NULL,
    breaks = c("Gene Wise Knee", "Knee of Mean Curve"),
    values = c(
      "Gene Wise Knee" = blue_color,
      "Knee of Mean Curve" = red_color
    )
  ) +
guides(color = guide_legend(override.aes = list(
linetype = c(0, 1), # dot-only, line-only
shape = c(16, NA) # dot for gene knee; no dot for mean knee
)))  +
scale_x_reverse(
breaks = seq(
min(res$delta_smooth$grid, na.rm = TRUE),
max(res$delta_smooth$grid, na.rm = TRUE),
by = 10
)
) +
  labs(
    title = "Iδ vs Width (scope)",
    x = "Grid width (μm)",
    y = "Iδ"
  ) +
theme_minimal(base_size = 8) +
theme(
panel.border = element_rect(color = "black", fill = NA),
panel.background = element_rect(fill = "#c0c0c0", color = NA),
plot.background = element_rect(fill = "#ffffff", color = NA),
axis.title = element_text(size = 9),
axis.text.y = element_text(size = 8),
axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
axis.ticks = element_line(color = "black"),
plot.title = element_text(hjust = 0.5, size = 10),
legend.position = c(0.02, 0.98),
legend.justification = c(0, 1),
legend.background = element_blank(),
legend.key = element_blank(),
legend.text = element_text(size = 7),
plot.margin = unit(c(1, 1, 1, 1), "cm")
)+
  scale_y_continuous(
    limits = c(0, min(100, max(avg_curve_mean$iDelta_mean, na.rm = TRUE))),
    oob = scales::squish,
    breaks = scales::breaks_width(5)
  )

ggsave(
filename = file.path(figures_cosmx_dir, "scope_iDelta_all_genes_annotated_knee_stable.png"),
plot = p_all_genes_annotated,
width = 6,
height = 5,
dpi = 600
)

# ==============================
# scope: Additional plots (LQC style)
# 1) Gene-wise knee distribution: violin + histogram
# 2) iDelta at overall knee vs cell-level iDelta: violin + scatter
# 3) iDelta at each-gene knee vs cell-level iDelta: violin + scatter
# Aesthetics aligned with genescope/R/7.LeeLQCPlot.r
# ==============================

theme_lqc <- function() {
  theme_bw(base_size = 8) +
    theme(
      plot.title         = element_text(size = 12, colour = "black"),
      panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.background = element_rect(fill = "#ffffff", color = NA),
      panel.grid.major.y = element_line(colour = "#c0c0c0", size = .2),
      panel.border       = element_rect(colour = "black", fill = NA, size = .5),
      axis.line          = element_line(colour = "black", size = .3),
      axis.ticks         = element_line(colour = "black", size = .3),
      axis.text          = element_text(size = 10, colour = "black"),
      axis.title         = element_text(size = 12, colour = "black"),
      plot.margin        = unit(c(0, 0, 0, 0), "cm")
    )
}

# ---------- 1) Gene-wise knee distribution ----------
if (nrow(valid_knees) > 0) {
  knee_df <- valid_knees %>% select(gene, uik_knee) %>% filter(is.finite(uik_knee))

  p_scope_knee_violin <- ggplot(knee_df, aes(x = "Genes", y = uik_knee)) +
    geom_violin(fill = "white", colour = "black", size = .3, trim = FALSE) +
    geom_boxplot(width = .15, outlier.size = .4, outlier.stroke = .2,
                 fill = "white", colour = "black", size = .25) +
    labs(title = "Gene Wise Knee (UIK, scope2)", x = NULL, y = "Knee (μm)") +
    theme_lqc()

  p_scope_knee_hist <- ggplot(knee_df, aes(x = uik_knee)) +
    geom_histogram(bins = 30, fill = "white", colour = "black", size = .3) +
    labs(title = "Gene Wise Knee (UIK) Histogram (scope)", x = "Knee (μm)", y = "Frequency") +
    theme_lqc()

  ggsave(
    filename = file.path(figures_cosmx_dir, "scope_geneWiseKnee_UIK_violin.png"),
    plot = p_scope_knee_violin,
    width = 3,
    height = 2.5,
    dpi = 600
  )
  ggsave(
    filename = file.path(figures_cosmx_dir, "scope_geneWiseKnee_UIK_histogram.png"),
    plot = p_scope_knee_hist,
    width = 6,
    height = 5,
    dpi = 600
  )
}

# ---------- 2) iDelta at overall knee vs cell-level iDelta ----------
if (is.finite(overall_knee_mean)) {
  curve_overall <- res$delta_smooth %>%
    group_by(gene) %>%
    slice_min(abs(grid - overall_knee_mean), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(gene, grid, iDelta_curve = iDelta)

  curve_grid30 <- res$delta_smooth %>%
    group_by(gene) %>%
    slice_min(abs(grid - 30), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(gene, iDelta_grid30 = iDelta)

  df_overall_vs_cell <- curve_overall %>%
    inner_join(cell_idelta, by = "gene") %>%
    filter(is.finite(iDelta_curve), is.finite(cell_iDelta))

  df_overall_long <- df_overall_vs_cell %>%
    select(gene, `Knee of Mean Curve` = iDelta_curve, `Cell` = cell_iDelta) %>%
    pivot_longer(-gene, names_to = "type", values_to = "iDelta")

  p_scope_overall_violin <- ggplot(df_overall_long, aes(x = type, y = iDelta)) +
    geom_violin(fill = "white", colour = "black", size = .3, trim = FALSE) +
    geom_boxplot(width = .25, outlier.size = .4, outlier.stroke = .2,
                 fill = "white", colour = "black", size = .25) +
    labs(title = "Iδ: Knee of Mean Curve vs Cell (scope)", x = NULL, y = "Iδ") +
    theme_lqc()

  df_grid30_vs_cell <- curve_grid30 %>%
    inner_join(cell_idelta, by = "gene") %>%
    filter(is.finite(iDelta_grid30), is.finite(cell_iDelta))

  iqr_x <- IQR(df_grid30_vs_cell$iDelta_grid30, na.rm = TRUE)
  qx1  <- quantile(df_grid30_vs_cell$iDelta_grid30, 0.25, na.rm = TRUE)
  qx3  <- quantile(df_grid30_vs_cell$iDelta_grid30, 0.75, na.rm = TRUE)
  x_low <- qx1 - 1.5 * iqr_x
  x_high <- qx3 + 1.5 * iqr_x

  iqr_y <- IQR(df_grid30_vs_cell$cell_iDelta, na.rm = TRUE)
  qy1  <- quantile(df_grid30_vs_cell$cell_iDelta, 0.25, na.rm = TRUE)
  qy3  <- quantile(df_grid30_vs_cell$cell_iDelta, 0.75, na.rm = TRUE)
  y_low <- qy1 - 1.5 * iqr_y
  y_high <- qy3 + 1.5 * iqr_y

  df_grid30_vs_cell <- df_grid30_vs_cell %>%
    mutate(
      x_out = iDelta_grid30 < x_low | iDelta_grid30 > x_high,
      y_out = cell_iDelta    < y_low | cell_iDelta    > y_high,
      out_flag = case_when(
        x_out & y_out ~ "Both",
        x_out ~ "X only",
        y_out ~ "Y only",
        TRUE ~ "In range"
      )
    )

  axis_range <- range(
    c(df_grid30_vs_cell$iDelta_grid30, df_grid30_vs_cell$cell_iDelta),
    na.rm = TRUE
  )

  p_scope_overall_scatter <- ggplot(df_grid30_vs_cell, aes(x = iDelta_grid30, y = cell_iDelta)) +
    geom_point(aes(fill = out_flag),
               shape = 21, size = 1.2, stroke = .35, colour = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = .3) +
    labs(title = "Iδ: Grid 30 μm vs Cell (scope)",
         x = "Iδ at Grid 30 μm", y = "Cell-level Iδ") +
    scale_x_continuous(limits = axis_range, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = axis_range, expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(
      values = c(
        "X only"  = "#1f78b4",
        "Y only"  = "#e31a1c",
        "Both"    = "#6a3d9a",
        "In range"= "white"
      ),
      name = "Outlier (1.5×IQR)"
    ) +
    theme_lqc() +
    theme(
      legend.background = element_blank(),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10),
      legend.key        = element_blank()
    )

  ggsave(
    filename = file.path(figures_dir, "scope_overallKnee_iDelta_vs_cell_violin.png"),
    plot = p_scope_overall_violin,
    width = 6,
    height = 5,
    dpi = 600
  )
  ggsave(
    filename = file.path(figures_dir, "scope_grid30_iDelta_vs_cell_scatter.png"),
    plot = p_scope_overall_scatter,
    width = 6,
    height = 5,
    dpi = 600
  )
}
# ---------- 3) iDelta at each-gene knee vs cell-level iDelta ----------
if (nrow(res$gene_knee_with_y) > 0) {
  df_geneK_vs_cell <- res$gene_knee_with_y %>%
    inner_join(cell_idelta, by = "gene") %>%
    rename(iDelta_curve_knee = knee_y) %>%
    filter(is.finite(iDelta_curve_knee), is.finite(cell_iDelta))

  df_geneK_long <- df_geneK_vs_cell %>%
    select(gene, `Gene Wise Knee` = iDelta_curve_knee, `Cell` = cell_iDelta) %>%
    pivot_longer(-gene, names_to = "type", values_to = "iDelta")

  p_scope_geneK_violin <- ggplot(df_geneK_long, aes(x = type, y = iDelta)) +
    geom_violin(fill = "white", colour = "black", size = .3, trim = FALSE) +
    geom_boxplot(width = .25, outlier.size = .4, outlier.stroke = .2,
                 fill = "white", colour = "black", size = .25) +
    labs(title = "Iδ: Gene Wise Knee vs Cell (scope)", x = NULL, y = "Iδ") +
    theme_lqc()

  p_scope_geneK_scatter <- ggplot(df_geneK_vs_cell, aes(x = iDelta_curve_knee, y = cell_iDelta)) +
    geom_point(shape = 21, size = 1.2, fill = "white", colour = "black", stroke = .2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = .3) +
    labs(title = "Iδ: Gene Wise Knee vs Cell (scope)",
         x = "Curve Iδ at Gene Wise Knee", y = "Cell-level Iδ") +
    theme_lqc()

  ggsave(
    filename = file.path(figures_dir, "scope_geneKnee_iDelta_vs_cell_violin.png"),
    plot = p_scope_geneK_violin,
    width = 6,
    height = 5,
    dpi = 600
  )
  ggsave(
    filename = file.path(figures_dir, "scope_geneKnee_iDelta_vs_cell_scatter.png"),
    plot = p_scope_geneK_scatter,
    width = 6,
    height = 5,
    dpi = 600
  )
}

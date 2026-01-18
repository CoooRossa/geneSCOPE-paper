#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

display_method_name <- function(x) {
  x_chr <- as.character(x)
  key <- tolower(trimws(x_chr))
  if (key == "genescope" || startsWith(key, "genescope")) return("geneSCOPE")
  if (key == "giotto" || key == "giotto_r" || startsWith(key, "giotto")) return("Giotto Suite")
  if (key == "hotspot" || key == "hotspot_py" || startsWith(key, "hotspot")) return("Hotspot")
  if (key == "seagal" || startsWith(key, "seagal")) return("Seagal")
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
  if (any(missing)) out[missing] <- "#7F7F7F"
  stats::setNames(as.character(out), method_labels)
}

theme_paper_like <- function() {
  theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

safe_fread <- function(path) {
  tryCatch(
    data.table::fread(path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
}

stats_detect_schema <- function(nms) {
  nms <- as.character(nms)
  time_col <- NA_character_
  if ("t_sec" %in% nms) time_col <- "t_sec"
  if ("epoch" %in% nms) time_col <- "epoch"
  if ("ts_epoch" %in% nms) time_col <- "ts_epoch"
  if ("epoch_sec" %in% nms) time_col <- "epoch_sec"
  if (!is.na(time_col)) return(list(time_col = time_col))
  NULL
}

parse_stats_tsv <- function(path) {
  dt <- safe_fread(path)
  if (is.null(dt) || !nrow(dt)) return(NULL)

  schema <- stats_detect_schema(names(dt))
  if (is.null(schema)) return(NULL)

  time_abs <- suppressWarnings(as.numeric(dt[[schema$time_col]]))
  if (!any(is.finite(time_abs))) return(NULL)
  t0 <- min(time_abs[is.finite(time_abs)], na.rm = TRUE)
  dt[, t_sec := time_abs - t0]

  mem_candidates <- c("mem_current_bytes", "memory_current_bytes", "mem_bytes", "memory_current_bytes")
  mem_col <- intersect(mem_candidates, names(dt))[1]
  if (is.na(mem_col)) return(NULL)
  mem_bytes <- suppressWarnings(as.numeric(dt[[mem_col]]))
  mem_bytes[!is.finite(mem_bytes)] <- NA_real_
  dt[, mem_gb := mem_bytes / (1024^3)]

  cpu_cum_col <- intersect(c("cpu_usage_usec", "cpu_usec"), names(dt))[1]
  cpu_cum <- if (!is.na(cpu_cum_col)) suppressWarnings(as.numeric(dt[[cpu_cum_col]])) else rep(NA_real_, nrow(dt))
  cpu_cum[!is.finite(cpu_cum)] <- NA_real_

  dt_time <- c(NA_real_, diff(time_abs))
  dt_time[!is.finite(dt_time) | dt_time <= 0] <- NA_real_

  cpu_delta <- NULL
  if ("cpu_usage_delta_usec" %in% names(dt)) {
    cpu_delta <- suppressWarnings(as.numeric(dt[["cpu_usage_delta_usec"]]))
  } else if (any(is.finite(cpu_cum))) {
    cpu_delta <- c(NA_real_, diff(cpu_cum))
  } else {
    cpu_delta <- rep(NA_real_, nrow(dt))
  }
  cpu_delta[!is.finite(cpu_delta) | cpu_delta < 0] <- NA_real_

  cpu_cores <- cpu_delta / (dt_time * 1e6)
  cpu_cores[!is.finite(cpu_cores)] <- NA_real_
  dt[, cpu := cpu_cores]

  dt <- dt[is.finite(t_sec) & t_sec >= 0]
  if (!nrow(dt)) return(NULL)
  setorder(dt, t_sec)
  dt[, .(t_sec, cpu, mem_gb)]
}

parse_iso8601_to_epoch <- function(x) {
  x <- as.character(x)
  x <- gsub("([+-]\\d{2}):(\\d{2})$", "\\1\\2", x, perl = TRUE)
  ts <- suppressWarnings(as.POSIXct(x, format = "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"))
  as.numeric(ts)
}

parse_seagal_repeat_starts <- function(path) {
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (!length(lines)) return(data.table())

  pat <- paste0(
    "^\\[(\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}[+-]\\d{2}:\\d{2})\\]\\s+",
    "\\[INFO\\]\\s+",
    "(?:repeat_idx|replicate_idx)=(\\d+)\\s*$"
  )
  m <- regexec(pat, lines, perl = TRUE)
  hits <- regmatches(lines, m)
  hits <- hits[lengths(hits) == 3]
  if (!length(hits)) return(data.table())

  iso <- vapply(hits, `[[`, character(1), 2)
  idx <- suppressWarnings(as.integer(vapply(hits, `[[`, character(1), 3)))
  epoch <- parse_iso8601_to_epoch(iso)

  dt <- data.table(repeat_idx = idx, start_iso = iso, start_epoch = epoch)
  dt <- dt[is.finite(start_epoch) & is.finite(repeat_idx)]
  if (!nrow(dt)) return(data.table())
  dt <- unique(dt)
  setorder(dt, repeat_idx, start_epoch)
  dt
}

parse_stats_epoch_cpu_mem <- function(path) {
  dt <- safe_fread(path)
  if (is.null(dt) || !nrow(dt)) return(NULL)

  schema <- stats_detect_schema(names(dt))
  if (is.null(schema)) return(NULL)
  epoch <- suppressWarnings(as.numeric(dt[[schema$time_col]]))
  if (!any(is.finite(epoch))) return(NULL)

  mem_candidates <- c("mem_current_bytes", "memory_current_bytes", "mem_bytes")
  mem_col <- intersect(mem_candidates, names(dt))[1]
  if (is.na(mem_col)) return(NULL)
  mem_bytes <- suppressWarnings(as.numeric(dt[[mem_col]]))
  mem_bytes[!is.finite(mem_bytes) | mem_bytes < 0] <- NA_real_

  cpu_cum_col <- intersect(c("cpu_usage_usec", "cpu_usec"), names(dt))[1]
  if (is.na(cpu_cum_col)) return(NULL)
  cpu_cum_usec <- suppressWarnings(as.numeric(dt[[cpu_cum_col]]))
  cpu_cum_usec[!is.finite(cpu_cum_usec) | cpu_cum_usec < 0] <- NA_real_

  out <- data.table(epoch = epoch, cpu_cum_usec = cpu_cum_usec, mem_bytes = mem_bytes)
  out <- out[is.finite(epoch)]
  if (!nrow(out)) return(NULL)
  setorder(out, epoch)
  out
}

segment_epoch_timeseries <- function(dt, start_epoch, end_epoch = NA_real_) {
  if (is.null(dt) || !nrow(dt)) return(NULL)
  seg <- dt[epoch >= start_epoch]
  if (is.finite(end_epoch)) seg <- seg[epoch < end_epoch]
  if (nrow(seg) < 2) return(NULL)

  t0 <- min(seg$epoch[is.finite(seg$epoch)], na.rm = TRUE)
  seg[, t_sec := epoch - t0]

  dt_time <- c(NA_real_, diff(seg$epoch))
  dt_time[!is.finite(dt_time) | dt_time <= 0] <- NA_real_
  cpu_delta <- c(NA_real_, diff(seg$cpu_cum_usec))
  cpu_delta[!is.finite(cpu_delta) | cpu_delta < 0] <- NA_real_

  seg[, cpu := cpu_delta / (dt_time * 1e6)]
  seg[, mem_gb := mem_bytes / (1024^3)]

  seg <- seg[is.finite(t_sec) & t_sec >= 0]
  if (!nrow(seg)) return(NULL)
  seg[, .(t_sec, cpu, mem_gb)]
}

list_repeat_dirs <- function(method_root) {
  if (!dir.exists(method_root)) return(character(0))
  entries <- list.files(method_root, full.names = TRUE, recursive = FALSE)
  entries <- entries[file.info(entries)$isdir %in% TRUE]
  keep <- grepl("^repeat_[0-9]+$", basename(entries))
  entries <- entries[keep]
  idx <- suppressWarnings(as.integer(sub("^repeat_([0-9]+)$", "\\1", basename(entries))))
  ord <- order(idx, basename(entries), na.last = TRUE)
  entries[ord]
}

detect_seagal_segmented_timeseries <- function(method_root, rep_dirs) {
  stats_path <- file.path(method_root, "stats.tsv")
  runlog_path <- file.path(method_root, "run.log")
  if (!file.exists(stats_path) || !file.exists(runlog_path)) return(list())

  starts <- parse_seagal_repeat_starts(runlog_path)
  if (!nrow(starts)) return(list())

  stats_dt <- parse_stats_epoch_cpu_mem(stats_path)
  if (is.null(stats_dt) || !nrow(stats_dt)) return(list())

  rep_ids <- basename(rep_dirs)
  rep_idx <- suppressWarnings(as.integer(sub("^repeat_([0-9]+)$", "\\1", rep_ids)))
  ord <- order(rep_idx, rep_ids, na.last = TRUE)
  rep_ids <- rep_ids[ord]
  rep_idx <- rep_idx[ord]

  out <- list()
  for (i in seq_along(rep_ids)) {
    rid <- rep_ids[[i]]
    idx <- rep_idx[[i]]
    if (!is.finite(idx)) next

    s <- starts[repeat_idx == idx][1]
    if (!nrow(s) || !is.finite(s$start_epoch)) next
    n <- starts[repeat_idx == (idx + 1)][1]
    end_epoch <- if (nrow(n) && is.finite(n$start_epoch)) n$start_epoch else NA_real_

    seg <- segment_epoch_timeseries(stats_dt, start_epoch = s$start_epoch, end_epoch = end_epoch)
    if (is.null(seg) || nrow(seg) < 2) next
    out[[rid]] <- seg
  }
  out
}

resample_to_grid <- function(dt, t_grid, y_col) {
  dt <- as.data.table(dt)
  if (!nrow(dt)) return(rep(NA_real_, length(t_grid)))
  sub <- dt[is.finite(t_sec) & is.finite(get(y_col))]
  if (nrow(sub) < 2) return(rep(NA_real_, length(t_grid)))
  sub <- sub[order(t_sec)]
  stats::approx(x = sub$t_sec, y = sub[[y_col]], xout = t_grid, rule = 1, ties = "ordered")$y
}

compute_mean_curves <- function(rep_long, t_step, xmax) {
  rep_long <- as.data.table(rep_long)
  if (!nrow(rep_long)) return(data.table())

  xmax <- suppressWarnings(as.numeric(xmax))
  if (!is.finite(xmax) || xmax <= 0) {
    xmax <- suppressWarnings(max(rep_long$t_sec, na.rm = TRUE))
  }
  if (!is.finite(xmax) || xmax <= 0) return(data.table())

  t_step <- suppressWarnings(as.numeric(t_step))
  if (!is.finite(t_step) || t_step <= 0) t_step <- 1
  t_grid <- seq(0, floor(xmax / t_step) * t_step, by = t_step)

  out <- list()
  for (m in unique(rep_long$method)) {
    sub_m <- rep_long[method == m]
    reps <- unique(sub_m$replicate_id)
    if (!length(reps)) next

    cpu_mat <- vapply(reps, function(rid) resample_to_grid(sub_m[replicate_id == rid], t_grid, "cpu_pct"), numeric(length(t_grid)))
    mem_mat <- vapply(reps, function(rid) resample_to_grid(sub_m[replicate_id == rid], t_grid, "mem_gb"), numeric(length(t_grid)))

    cpu_mean <- rowMeans(cpu_mat, na.rm = TRUE)
    mem_mean <- rowMeans(mem_mat, na.rm = TRUE)
    cpu_n <- rowSums(is.finite(cpu_mat))
    mem_n <- rowSums(is.finite(mem_mat))

    out[[m]] <- data.table(
      method = m,
      t_sec = t_grid,
      cpu_pct = ifelse(cpu_n > 0, cpu_mean, NA_real_),
      mem_gb = ifelse(mem_n > 0, mem_mean, NA_real_),
      n_cpu = cpu_n,
      n_mem = mem_n
    )
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

time_weighted_mean <- function(t_sec, y) {
  t_sec <- suppressWarnings(as.numeric(t_sec))
  y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(t_sec) & is.finite(y)
  if (sum(ok) < 2) return(NA_real_)

  dt <- data.table(t_sec = t_sec[ok], y = y[ok])
  setorder(dt, t_sec)
  dt <- dt[!duplicated(t_sec)]
  if (nrow(dt) < 2) return(NA_real_)

  dt_t <- diff(dt$t_sec)
  dt_t[!is.finite(dt_t) | dt_t <= 0] <- NA_real_
  if (!any(is.finite(dt_t))) return(NA_real_)

  y0 <- dt$y[-nrow(dt)]
  y1 <- dt$y[-1]
  area <- sum(((y0 + y1) / 2) * dt_t, na.rm = TRUE)
  duration <- sum(dt_t, na.rm = TRUE)
  if (!is.finite(duration) || duration <= 0) return(NA_real_)
  area / duration
}

plot_panels <- function(rep_long, mean_long, metric_col, y_label, out_path,
                        method_levels, color_map, xmax, ymax, title, y_breaks = waiver()) {
  rep_long <- as.data.table(rep_long)
  mean_long <- as.data.table(mean_long)

  rep_plot <- rep_long[, .(
    method = factor(method, levels = method_levels),
    replicate_id = replicate_id,
    t_sec = t_sec,
    value = get(metric_col)
  )]
  mean_plot <- mean_long[, .(
    method = factor(method, levels = method_levels),
    t_sec = t_sec,
    value = get(metric_col)
  )]

  p <- ggplot() +
    geom_line(
      data = rep_plot,
      aes(x = t_sec, y = value, group = replicate_id, color = method),
      alpha = 0.30,
      linewidth = 1.1,
      na.rm = TRUE
    ) +
    geom_line(
      data = mean_plot,
      aes(x = t_sec, y = value, color = method),
      alpha = 1.0,
      linewidth = 1.6,
      na.rm = TRUE
    ) +
    facet_grid(method ~ ., scales = "fixed") +
    scale_color_manual(values = color_map, guide = "none") +
    scale_x_continuous(
      limits = c(0, xmax),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = y_breaks,
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    theme_paper_like() +
    theme(
      strip.background = element_blank(),
      strip.text.y = element_text(face = "bold", size = 12, angle = 0),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.spacing = grid::unit(0.7, "lines"),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    ) +
    labs(title = title, x = "Time (s)", y = y_label)

  ggsave(out_path, p, width = 8.2 * 1.5, height = 10.5, dpi = 180)
  invisible(p)
}

option_list <- list(
  make_option(c("--run_root"), type = "character", default = "",
              help = "Run root containing genescope/giotto/hotspot/seagal directories"),
  make_option(c("--outdir"), type = "character", default = "",
              help = "Output directory for plots (default: <run_root>/plots_runtime)"),
  make_option(c("--n_runs"), type = "integer", default = 5L,
              help = "Number of repeats per method to plot (default: 5)"),
  make_option(c("--t_step"), type = "double", default = 1.0,
              help = "Resampling step (sec) for mean curves (default: 1)"),
  make_option(c("--xmax_sec"), type = "double", default = NA_real_,
              help = "Fixed x-axis max seconds (default: auto)"),
  make_option(c("--cpu_ymax"), type = "double", default = 2000.0,
              help = "CPU y-axis max percent (default: 2000)"),
  make_option(c("--mem_ymax"), type = "double", default = NA_real_,
              help = "Memory y-axis max GB (default: auto)")
)

opt <- parse_args(OptionParser(option_list = option_list))

run_root <- as.character(opt$run_root)
if (!nzchar(run_root)) stop("--run_root is required")
run_root <- normalizePath(run_root, mustWork = TRUE)

outdir <- as.character(opt$outdir)
if (!nzchar(outdir)) outdir <- file.path(run_root, "plots_runtime")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_runs <- suppressWarnings(as.integer(opt$n_runs))
if (!is.finite(n_runs) || n_runs < 1L) n_runs <- 5L

method_levels <- c("geneSCOPE", "Giotto Suite", "Hotspot", "Seagal")
color_map <- bench_base_color_map(method_levels)

method_cfg <- data.table(
  method_key = c("genescope", "giotto", "hotspot", "seagal"),
  method_label = vapply(c("genescope", "giotto", "hotspot", "seagal"), display_method_name, character(1)),
  dir_name = c("genescope", "giotto", "hotspot", "seagal")
)

rep_series <- list()
for (i in seq_len(nrow(method_cfg))) {
  m <- method_cfg[i]
  method_root <- file.path(run_root, m$dir_name)
  rep_dirs <- list_repeat_dirs(method_root)
  if (!length(rep_dirs)) next

  if (m$method_key == "seagal") {
    # Seagal: stats.tsv at method_root; segment into per-repeat curves using run.log timestamps.
    rep_dirs <- rep_dirs[file.exists(file.path(rep_dirs, "modules.tsv"))]
    rep_dirs <- head(rep_dirs, n_runs)
    segs <- detect_seagal_segmented_timeseries(method_root, rep_dirs)
    if (!length(segs)) next
    for (rid in names(segs)) {
      dt <- segs[[rid]]
      if (is.null(dt) || nrow(dt) < 2) next
      dt[, `:=`(
        method = m$method_label,
        replicate_id = rid,
        cpu_pct = cpu * 100.0
      )]
      rep_series[[paste0(m$method_key, "::", rid)]] <- dt[, .(method, replicate_id, t_sec, cpu_pct, mem_gb)]
    }
  } else {
    rep_dirs <- head(rep_dirs, n_runs)
    for (repdir in rep_dirs) {
      rid <- basename(repdir)
      stats_path <- file.path(repdir, "stats.tsv")
      if (!file.exists(stats_path)) next
      dt <- parse_stats_tsv(stats_path)
      if (is.null(dt) || nrow(dt) < 2) next
      dt[, `:=`(
        method = m$method_label,
        replicate_id = rid,
        cpu_pct = cpu * 100.0
      )]
      rep_series[[paste0(m$method_key, "::", rid)]] <- dt[, .(method, replicate_id, t_sec, cpu_pct, mem_gb)]
    }
  }
}

rep_long <- rbindlist(rep_series, use.names = TRUE, fill = TRUE)
if (!nrow(rep_long)) stop("No stats time series detected under: ", run_root)

rep_long <- rep_long[method %in% method_levels]
rep_long[, method := factor(method, levels = method_levels)]

data.table::fwrite(rep_long, file.path(outdir, "runtime_timeseries_long.tsv"), sep = "\t")

xmax <- suppressWarnings(as.numeric(opt$xmax_sec))
if (!is.finite(xmax) || xmax <= 0) {
  xmax <- suppressWarnings(max(rep_long$t_sec, na.rm = TRUE))
  if (!is.finite(xmax) || xmax <= 0) xmax <- 1
  xmax <- ceiling(xmax)
}

mean_long <- compute_mean_curves(rep_long, t_step = opt$t_step, xmax = xmax)
mean_long <- mean_long[method %in% method_levels]
mean_long[, method := factor(method, levels = method_levels)]

rep_stats <- rep_long[
  ,
  .(
    n_points = .N,
    duration_sec = suppressWarnings(max(t_sec, na.rm = TRUE)),
    cpu_pct_mean = time_weighted_mean(t_sec, cpu_pct),
    cpu_pct_max = suppressWarnings(max(cpu_pct, na.rm = TRUE)),
    mem_gb_mean = time_weighted_mean(t_sec, mem_gb),
    mem_gb_max = suppressWarnings(max(mem_gb, na.rm = TRUE))
  ),
  by = .(method, replicate_id)
]
rep_stats[, method := factor(method, levels = method_levels)]
setorder(rep_stats, method, replicate_id)
data.table::fwrite(rep_stats, file.path(outdir, "runtime_run_stats.tsv"), sep = "\t")

method_stats <- rep_stats[
  ,
  .(
    n_runs = .N,
    duration_sec_mean = mean(duration_sec, na.rm = TRUE),
    duration_sec_sd = stats::sd(duration_sec, na.rm = TRUE),
    cpu_pct_mean_mean = mean(cpu_pct_mean, na.rm = TRUE),
    cpu_pct_mean_sd = stats::sd(cpu_pct_mean, na.rm = TRUE),
    cpu_pct_max_mean = mean(cpu_pct_max, na.rm = TRUE),
    cpu_pct_max_sd = stats::sd(cpu_pct_max, na.rm = TRUE),
    mem_gb_mean_mean = mean(mem_gb_mean, na.rm = TRUE),
    mem_gb_mean_sd = stats::sd(mem_gb_mean, na.rm = TRUE),
    mem_gb_max_mean = mean(mem_gb_max, na.rm = TRUE),
    mem_gb_max_sd = stats::sd(mem_gb_max, na.rm = TRUE)
  ),
  by = .(method)
]
method_stats[, method := factor(method, levels = method_levels)]
setorder(method_stats, method)
data.table::fwrite(method_stats, file.path(outdir, "runtime_method_stats.tsv"), sep = "\t")
data.table::fwrite(mean_long, file.path(outdir, "runtime_mean_curves.tsv"), sep = "\t")

cpu_ymax <- suppressWarnings(as.numeric(opt$cpu_ymax))
if (!is.finite(cpu_ymax) || cpu_ymax <= 0) cpu_ymax <- 2000

mem_ymax <- suppressWarnings(as.numeric(opt$mem_ymax))
if (!is.finite(mem_ymax) || mem_ymax <= 0) {
  mem_ymax <- suppressWarnings(max(rep_long$mem_gb, na.rm = TRUE))
  if (!is.finite(mem_ymax) || mem_ymax <= 0) mem_ymax <- 1
  mem_ymax <- ceiling(mem_ymax)
}

plot_panels(
  rep_long, mean_long,
  metric_col = "cpu_pct",
  y_label = "CPU Usage (%)",
  out_path = file.path(outdir, "runtime_cpu_panels.png"),
  method_levels = method_levels,
  color_map = color_map,
  xmax = xmax,
  ymax = cpu_ymax,
  title = sprintf("Runtime CPU Usage (%d runs)", n_runs),
  y_breaks = seq(0, cpu_ymax, by = 500)
)

plot_panels(
  rep_long, mean_long,
  metric_col = "mem_gb",
  y_label = "Memory (GB)",
  out_path = file.path(outdir, "runtime_mem_panels.png"),
  method_levels = method_levels,
  color_map = color_map,
  xmax = xmax,
  ymax = mem_ymax,
  title = sprintf("Runtime Memory Usage (%d runs)", n_runs)
)

message("Wrote: ", file.path(outdir, "runtime_cpu_panels.png"))
message("Wrote: ", file.path(outdir, "runtime_mem_panels.png"))
message("Wrote: ", file.path(outdir, "runtime_timeseries_long.tsv"))
message("Wrote: ", file.path(outdir, "runtime_run_stats.tsv"))
message("Wrote: ", file.path(outdir, "runtime_method_stats.tsv"))
message("Wrote: ", file.path(outdir, "runtime_mean_curves.tsv"))

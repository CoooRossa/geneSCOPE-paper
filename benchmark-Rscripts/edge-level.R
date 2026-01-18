#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

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
  x <- trimws(x)
  if (!nzchar(x)) return(character(0))
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  trimws(parts[nzchar(parts)])
}

parse_top_n_list <- function(x) {
  vals <- suppressWarnings(as.integer(split_csv(x)))
  vals <- vals[is.finite(vals) & vals > 0L]
  if (!length(vals)) return(integer(0))
  unique(vals)
}

parse_double_list_allow_na <- function(x) {
  parts <- split_csv(x)
  if (!length(parts)) return(numeric(0))
  out <- vapply(parts, function(tok) {
    tok <- trimws(tolower(tok))
    if (!nzchar(tok) || tok %in% c("na", "nan", "none", "null")) return(NA_real_)
    suppressWarnings(as.numeric(tok))
  }, numeric(1))
  out
}

ranking_key_desc <- "weight"

normalize_ranking_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("weight", "w")) return("weight")
  if (x %in% c("abs", "abs_weight", "absweight", "abs(weight)", "abs_w")) return("abs(weight)")
  stop("Invalid ranking key: ", x, " (expected: weight or abs_weight)")
}

parse_ranking_keys <- function(keys_csv) {
  keys <- split_csv(keys_csv)
  if (!length(keys)) return(character(0))
  keys <- vapply(keys, normalize_ranking_key, character(1))
  unique(keys)
}

ranking_key_slug <- function(ranking_key) {
  ranking_key <- normalize_ranking_key(ranking_key)
  if (identical(ranking_key, "weight")) return("weight")
  if (identical(ranking_key, "abs(weight)")) return("absweight")
  stop("Unsupported ranking_key: ", ranking_key)
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

display_score_name <- function(x) {
  x <- as.character(x)
  if (!nzchar(x) || is.na(x)) return(x)
  if (x == "weight") return("Edge Weight")
  if (x == "string_score") return("String Score")
  if (x == "nscore") return("Neighborhood Score")
  if (x == "fscore") return("Fusion Score")
  if (x == "pscore") return("Co-occurrence Score")
  if (x == "ascore") return("Coexpression Score")
  if (x == "escore") return("Experimental Score")
  if (x == "dscore") return("Database Score")
  if (x == "tscore") return("Text-mining Score")
  x
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

format_pred_true_definitions <- function(fdr_positive_threshold, positive_only, string_score_threshold) {
  fdr_positive_threshold <- suppressWarnings(as.numeric(fdr_positive_threshold))
  if (!is.finite(fdr_positive_threshold)) fdr_positive_threshold <- NA_real_
  string_score_threshold <- suppressWarnings(as.integer(string_score_threshold))
  if (!is.finite(string_score_threshold)) string_score_threshold <- NA_integer_

  pred <- if (is.finite(fdr_positive_threshold)) {
    paste0("P: fdr < ", signif(fdr_positive_threshold, 3))
  } else {
    "P: fdr < <unset>"
  }
  if (isTRUE(positive_only)) pred <- paste0(pred, " & weight > 0")

  truth <- if (is.finite(string_score_threshold)) {
    paste0("T: STRING score >= ", string_score_threshold)
  } else {
    "T: STRING score >= <unset>"
  }
  list(predicted = pred, truth = truth)
}

format_predicted_definition_by_method <- function(methods, fdr_thresholds, positive_only) {
  methods <- as.character(methods)
  fdr_thresholds <- suppressWarnings(as.numeric(fdr_thresholds))
  if (length(fdr_thresholds) != length(methods)) {
    return(format_pred_true_definitions(NA_real_, positive_only = positive_only, string_score_threshold = NA_integer_)$predicted)
  }
  parts <- vapply(seq_along(methods), function(i) {
    label <- display_method_name(methods[[i]])
    thr <- fdr_thresholds[[i]]
    if (!is.finite(thr) || thr >= 1) {
      paste0(label, ": no FDR filter")
    } else if (thr <= 0) {
      paste0(label, ": invalid FDR threshold")
    } else {
      paste0(label, ": fdr < ", signif(thr, 3))
    }
  }, character(1))
  pred <- paste0("P: ", paste(parts, collapse = "; "))
  if (isTRUE(positive_only)) pred <- paste0(pred, " & weight > 0")
  pred
}

paper_rename_summary_cols <- function(dt) {
  dt <- data.table::copy(as.data.table(dt))
  if ("method" %in% names(dt)) {
    dt[, method := vapply(method, display_method_name, character(1))]
  }
  rename_map <- c(
    method = "Method",
    top_n = "Top N",
    n_edges_used = "Edges Used",
    n_comparable_edges = "Comparable Edges",
    edge_coverage = "Edge Coverage",
    auprc = "AUPRC",
    auprc_defined = "AUPRC Defined",
    auprc_reason = "AUPRC Reason",
    auroc = "AUROC",
    auroc_defined = "AUROC Defined",
    auroc_reason = "AUROC Reason",
    ranking_key = "Ranking Key",
    refill_to_top_n_comparable = "Refill To Top N Comparable",
    refill_max_multiplier = "Refill Max Multiplier",
    refill_status = "Refill Status",
    refill_shortfall = "Refill Shortfall",
    n_candidates_pre_refill = "Candidates Pre Refill",
    n_mappable_candidates = "Mappable Candidates",
    fdr_threshold = "FDR Threshold",
    fdr_pool_size_requested = "FDR Pool Size Requested",
    fdr_pool_size_used = "FDR Pool Size Used",
    edges_all_n = "Edges All",
    edges_fdr_n = "Edges FDR",
    edges_pool_n = "Edges Pool",
    "precision@50" = "Precision at 50",
    "precision@100" = "Precision at 100",
    "precision@500" = "Precision at 500",
    "top_frac_1%_enrichment" = "Top 1 Percent Enrichment",
    top1p_n_edges = "Top 1 Percent Edges",
    top1p_pos_count = "Top 1 Percent Positive Count",
    top_pos_rate = "Top Positive Rate",
    bg_gene_count = "Background Gene Count",
    bg_pair_count = "Background Pair Count",
    bg_pos_count = "Background Positive Count",
    bg_pos_rate = "Background Positive Rate",
    n_obs_used = "Observations Used",
    n_vars_used = "Variables Used",
    edges_total_n = "Edges Total",
    edges_sig_n = "Edges Significant",
    edges_top_n = "Edges Top N",
    edges_weight_type = "Edge Weight Type",
    edges_pvalue_method = "P-value Method"
  )
  old <- intersect(names(rename_map), names(dt))
  if (length(old)) setnames(dt, old = old, new = unname(rename_map[old]))
  dt
}

paper_rename_details_cols <- function(dt) {
  dt <- data.table::copy(as.data.table(dt))
  if ("method" %in% names(dt)) {
    dt[, method := vapply(method, display_method_name, character(1))]
  }
  rename_map <- c(
    from = "From",
    to = "To",
    weight = "Weight",
    score = "Score",
    string_from = "STRING From",
    string_to = "STRING To",
    string_score = "String Score",
    label = "Label",
    method = "Method",
    top_n = "Top N",
    nscore = "Neighborhood Score",
    fscore = "Fusion Score",
    pscore = "Co-occurrence Score",
    ascore = "Coexpression Score",
    escore = "Experimental Score",
    dscore = "Database Score",
    tscore = "Text-mining Score"
  )
  old <- intersect(names(rename_map), names(dt))
  if (length(old)) setnames(dt, old = old, new = unname(rename_map[old]))
  dt
}

filter_edges_for_scoring <- function(edges_all, fdr_threshold, positive_only) {
  edges_all <- as.data.table(edges_all)
  if (!nrow(edges_all)) return(edges_all)
  if (isTRUE(positive_only)) {
    edges_all[is.finite(fdr) & fdr < fdr_threshold & is.finite(weight) & weight > 0]
  } else {
    edges_all[is.finite(fdr) & fdr < fdr_threshold & is.finite(weight)]
  }
}

subset_edges_fdr_pool <- function(edges_sig, pool_size) {
  pool_size <- as.integer(pool_size)
  if (!is.finite(pool_size) || pool_size <= 0L) return(edges_sig)
  if (!nrow(edges_sig) || nrow(edges_sig) <= pool_size) return(edges_sig)
  dt <- as.data.table(edges_sig)
  setorder(dt, fdr, from, to)
  dt[seq_len(pool_size)]
}

build_quantile_reference <- function(weight_list, n_bins = 1000L) {
  if (!length(weight_list)) return(NULL)
  probs <- seq(0, 1, length.out = n_bins)
  qmat <- vapply(
    weight_list,
    function(w) stats::quantile(w, probs = probs, na.rm = TRUE, type = 8, names = FALSE),
    numeric(length(probs))
  )
  ref <- rowMeans(qmat, na.rm = TRUE)
  list(probs = probs, ref = ref)
}

apply_quantile_map <- function(weights, ref) {
  if (is.null(ref) || !length(weights)) return(weights)
  ranks <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::frank(weights, ties.method = "average", na.last = "keep")
  } else {
    rank(weights, ties.method = "average", na.last = "keep")
  }
  p <- (ranks - 0.5) / length(weights)
  stats::approx(ref$probs, ref$ref, xout = p, rule = 2)$y
}

add_ranking_score <- function(dt) {
  if (!nrow(dt)) return(dt)
  if (identical(ranking_key_desc, "abs(weight)")) {
    dt[, score := abs(weight)]
  } else {
    dt[, score := as.numeric(weight)]
  }
  dt
}

sort_edges_for_ranking <- function(dt) {
  dt <- as.data.table(dt)
  if (!nrow(dt)) return(dt)
  if (identical(ranking_key_desc, "abs(weight)")) {
    dt[, score := abs(weight)]
    setorder(dt, -score)
    dt[, score := NULL]
  } else {
    setorder(dt, -weight)
  }
  dt
}

select_edges_with_refill <- function(edges_sorted, top_n, refill_to_top_n_comparable, refill_max_multiplier) {
  top_n <- as.integer(top_n)
  if (!is.finite(top_n) || top_n < 1L) top_n <- 1L
  refill_max_multiplier <- as.integer(refill_max_multiplier)
  if (!is.finite(refill_max_multiplier) || refill_max_multiplier < 1L) refill_max_multiplier <- 1L

  edges_sorted <- as.data.table(edges_sorted)
  if (!nrow(edges_sorted)) {
    return(list(
      edges_eval = edges_sorted[0],
      refill_status = if (isTRUE(refill_to_top_n_comparable)) "insufficient_mappable_edges" else "disabled",
      refill_shortfall = top_n,
      n_candidates_pre_refill = 0L,
      n_mappable_candidates = 0L
    ))
  }

  max_rows <- nrow(edges_sorted)
  if (!isTRUE(refill_to_top_n_comparable)) {
    candidate_n <- min(max_rows, top_n)
    candidates <- edges_sorted[seq_len(candidate_n)]
    mappable_mask <- is.finite(candidates$string_score)
    n_mappable <- sum(mappable_mask, na.rm = TRUE)
    return(list(
      edges_eval = candidates,
      refill_status = "disabled",
      refill_shortfall = 0L,
      n_candidates_pre_refill = as.integer(candidate_n),
      n_mappable_candidates = as.integer(n_mappable)
    ))
  }

  # Avoid integer overflow when top_n is very large.
  if (top_n >= max_rows) {
    candidate_n <- max_rows
  } else {
    candidate_n_raw <- suppressWarnings(as.double(top_n) * as.double(refill_max_multiplier))
    if (!is.finite(candidate_n_raw) || candidate_n_raw >= as.double(max_rows)) {
      candidate_n <- max_rows
    } else {
      candidate_n <- as.integer(candidate_n_raw)
      if (!is.finite(candidate_n) || candidate_n < 1L) candidate_n <- max_rows
    }
  }
  candidates <- edges_sorted[seq_len(candidate_n)]
  mappable_mask <- is.finite(candidates$string_score)
  n_mappable <- sum(mappable_mask, na.rm = TRUE)
  mapped <- candidates[mappable_mask]
  if (nrow(mapped) > top_n) mapped <- mapped[seq_len(top_n)]

  status <- if (n_mappable >= top_n) "ok" else "insufficient_mappable_edges"
  shortfall <- max(0L, top_n - nrow(mapped))
  list(
    edges_eval = mapped,
    refill_status = status,
    refill_shortfall = as.integer(shortfall),
    n_candidates_pre_refill = as.integer(candidate_n),
    n_mappable_candidates = as.integer(n_mappable)
  )
}

# --- Metrics (copied from geneSCOPE internal code) ---
.compute_auprc <- function(labels, scores) {
  ok <- !is.na(labels) & !is.na(scores)
  labels <- labels[ok]
  scores <- scores[ok]
  if (!length(labels)) return(NA_real_)
  labels <- as.integer(labels) > 0
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  tp <- cumsum(labels)
  fp <- cumsum(!labels)
  recall <- tp / n_pos
  precision <- tp / pmax(tp + fp, 1)
  recall <- c(0, recall)
  precision <- c(1, precision)
  sum(diff(recall) * precision[-1])
}

.compute_auroc <- function(labels, scores) {
  ok <- !is.na(labels) & !is.na(scores)
  labels <- labels[ok]
  scores <- scores[ok]
  if (!length(labels)) return(NA_real_)
  labels <- as.integer(labels) > 0
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  sum_ranks_pos <- sum(ranks[labels])
  (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

.precision_at_k <- function(labels, scores, k = 100L) {
  ok <- !is.na(labels) & !is.na(scores)
  labels <- labels[ok]
  scores <- scores[ok]
  if (!length(labels)) return(NA_real_)
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  k <- min(length(labels), as.integer(k))
  if (k <= 0L) return(NA_real_)
  mean(as.integer(labels[seq_len(k)]))
}

.external_reference_precision_at_k <- function(labels, scores, k = c(50L, 100L, 500L)) {
  k <- as.integer(k)
  k <- k[is.finite(k) & k > 0]
  if (!length(k)) return(setNames(numeric(0), character(0)))
  vals <- vapply(k, function(kv) .precision_at_k(labels, scores, k = kv), numeric(1))
  names(vals) <- paste0("precision_at_", k)
  vals
}

.external_reference_enrichment <- function(n_pos_graph, n_graph, n_pos_bg, n_pairs_bg) {
  expected_rate <- if (!is.na(n_pairs_bg) && n_pairs_bg > 0) n_pos_bg / n_pairs_bg else NA_real_
  ratio <- NA_real_
  if (!is.na(expected_rate) && expected_rate > 0 && n_graph > 0) {
    ratio <- (n_pos_graph / n_graph) / expected_rate
  }

  pval <- NA_real_
  if (all(!is.na(c(n_pos_graph, n_graph, n_pos_bg, n_pairs_bg)))) {
    if (n_pairs_bg >= n_graph && n_pos_bg >= n_pos_graph) {
      n_non_graph <- n_pairs_bg - n_graph
      n_pos_non_graph <- n_pos_bg - n_pos_graph
      n_neg_graph <- n_graph - n_pos_graph
      n_neg_non_graph <- n_non_graph - n_pos_non_graph
      if (min(n_non_graph, n_pos_non_graph, n_neg_graph, n_neg_non_graph) >= 0) {
        tab <- matrix(c(n_pos_graph, n_neg_graph, n_pos_non_graph, n_neg_non_graph), nrow = 2)
        pval <- tryCatch(stats::fisher.test(tab, alternative = "greater")$p.value, error = function(e) NA_real_)
      }
    }
  }

  list(
    expected_positive_rate = expected_rate,
    enrichment_ratio = ratio,
    enrichment_p = pval
  )
}

.external_reference_top_sets <- function(n_edges, precision_k, top_frac = NULL, n_comparable = NULL) {
  sets <- list()
  n_comp <- n_edges
  if (!is.null(n_comparable) && is.finite(n_comparable)) n_comp <- n_comparable
  precision_k <- as.integer(precision_k)
  precision_k <- precision_k[is.finite(precision_k) & precision_k > 0L]
  if (length(precision_k)) {
    precision_k <- unique(precision_k)
    for (k in precision_k) {
      k_use <- min(k, n_comp)
      sets[[length(sets) + 1]] <- data.frame(
        set_id = paste0("top_k_", k_use),
        k = k_use,
        frac = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!is.null(top_frac) && is.finite(top_frac) && top_frac > 0 && n_comp > 0) {
    k_frac <- max(1L, ceiling(top_frac * n_comp))
    sets[[length(sets) + 1]] <- data.frame(
      set_id = paste0("top_frac_", formatC(top_frac * 100, format = "f", digits = 2), "pct"),
      k = k_frac,
      frac = top_frac,
      stringsAsFactors = FALSE
    )
  }
  if (!length(sets)) {
    return(data.frame(set_id = character(0), k = integer(0), frac = numeric(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, sets)
  out[!duplicated(out$set_id), , drop = FALSE]
}

.external_reference_top_enrichment <- function(edge_df, score_threshold, expected_rate, set_defs) {
  if (is.null(set_defs) || !nrow(set_defs)) {
    return(data.frame(
      set_id = character(0),
      k = integer(0),
      frac = numeric(0),
      n_top = integer(0),
      n_comparable = integer(0),
      n_pos = integer(0),
      positive_rate = numeric(0),
      enrichment_ratio = numeric(0),
      expected_positive_rate = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  if (!nrow(edge_df)) {
    empty <- set_defs
    empty$n_top <- 0L
    empty$n_comparable <- 0L
    empty$n_pos <- 0L
    empty$positive_rate <- NA_real_
    empty$enrichment_ratio <- NA_real_
    empty$expected_positive_rate <- expected_rate
    return(empty)
  }
  comparable_idx <- !is.na(edge_df$string_score) & !is.na(edge_df$score)
  edge_df <- edge_df[comparable_idx, , drop = FALSE]
  if (!nrow(edge_df)) {
    empty <- set_defs
    empty$n_top <- 0L
    empty$n_comparable <- 0L
    empty$n_pos <- 0L
    empty$positive_rate <- NA_real_
    empty$enrichment_ratio <- NA_real_
    empty$expected_positive_rate <- expected_rate
    return(empty)
  }
  ord <- order(edge_df$score, decreasing = TRUE, na.last = NA)
  edge_df <- edge_df[ord, , drop = FALSE]

  res <- lapply(seq_len(nrow(set_defs)), function(i) {
    k <- set_defs$k[i]
    n_top <- min(k, nrow(edge_df))
    top_df <- if (n_top > 0) edge_df[seq_len(n_top), , drop = FALSE] else edge_df[integer(0), , drop = FALSE]
    n_comparable <- n_top
    n_pos <- sum(top_df$string_score >= score_threshold, na.rm = TRUE)
    positive_rate <- if (n_top > 0) n_pos / n_top else NA_real_
    enrichment_ratio <- NA_real_
    if (!is.na(expected_rate) && expected_rate > 0 && !is.na(positive_rate)) {
      enrichment_ratio <- positive_rate / expected_rate
    }
    data.frame(
      set_id = set_defs$set_id[i],
      k = set_defs$k[i],
      frac = set_defs$frac[i],
      n_top = n_top,
      n_comparable = n_comparable,
      n_pos = n_pos,
      positive_rate = positive_rate,
      enrichment_ratio = enrichment_ratio,
      expected_positive_rate = expected_rate,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

score_edges_from_mapped <- function(edge_df, string_score_threshold, precision_k, top_frac, bg_gene_count, n_pairs_bg, n_pos_bg, keep_details) {
  edge_df <- as.data.table(edge_df)
  if (!"string_score" %in% names(edge_df)) edge_df[, string_score := NA_real_]
  if (!"score" %in% names(edge_df)) {
    if ("weight" %in% names(edge_df)) {
      if (identical(ranking_key_desc, "abs(weight)")) {
        edge_df[, score := abs(weight)]
      } else {
        edge_df[, score := as.numeric(weight)]
      }
    } else {
      edge_df[, score := NA_real_]
    }
  }
  edge_df[, string_score := suppressWarnings(as.numeric(string_score))]
  edge_df[, score := suppressWarnings(as.numeric(score))]
  edge_df[, label := ifelse(!is.na(string_score), string_score >= as.numeric(string_score_threshold), NA)]

  n_edges <- nrow(edge_df)
  ok_idx <- !is.na(edge_df$label) & !is.na(edge_df$score)
  n_comparable_edges <- sum(ok_idx)
  edge_coverage <- if (n_edges) n_comparable_edges / n_edges else NA_real_
  labels_ok <- edge_df$label[ok_idx]
  scores_ok <- edge_df$score[ok_idx]
  n_pos_labels <- sum(labels_ok, na.rm = TRUE)
  n_neg_labels <- sum(!labels_ok, na.rm = TRUE)

  auroc_reason <- "ok"
  auprc_reason <- "ok"
  auroc_defined <- TRUE
  auprc_defined <- TRUE
  if (n_comparable_edges == 0) {
    auroc_defined <- FALSE
    auprc_defined <- FALSE
    auroc_reason <- "no_comparable_edges"
    auprc_reason <- "no_comparable_edges"
  } else if (n_pos_labels == 0) {
    auroc_defined <- FALSE
    auprc_defined <- FALSE
    auroc_reason <- "no_positive_labels"
    auprc_reason <- "no_positive_labels"
  } else if (n_neg_labels == 0) {
    auroc_defined <- FALSE
    auprc_defined <- FALSE
    auroc_reason <- "no_negative_labels"
    auprc_reason <- "no_negative_labels"
  }

  auprc <- if (isTRUE(auprc_defined)) .compute_auprc(labels_ok, scores_ok) else NA_real_
  auroc <- if (isTRUE(auroc_defined)) .compute_auroc(labels_ok, scores_ok) else NA_real_
  prec_vals <- .external_reference_precision_at_k(edge_df$label, edge_df$score, k = precision_k)

  n_pos_graph <- sum(edge_df$string_score >= as.numeric(string_score_threshold), na.rm = TRUE)
  enrich <- .external_reference_enrichment(
    n_pos_graph = n_pos_graph,
    n_graph = n_comparable_edges,
    n_pos_bg = n_pos_bg,
    n_pairs_bg = n_pairs_bg
  )
  set_defs <- .external_reference_top_sets(n_edges, precision_k, top_frac = top_frac, n_comparable = n_comparable_edges)
  top_enrich <- .external_reference_top_enrichment(edge_df, score_threshold = as.numeric(string_score_threshold), expected_rate = enrich$expected_positive_rate, set_defs = set_defs)

  payload <- list(
    mapping = list(bg_gene_count = as.integer(bg_gene_count)),
    coverage = list(
      n_edges = n_edges,
      n_comparable_edges = n_comparable_edges,
      edge_coverage = edge_coverage,
      positive_rate = if (n_comparable_edges) mean(labels_ok) else NA_real_
    ),
    scores = list(
      auprc = auprc,
      auprc_defined = auprc_defined,
      auprc_reason = auprc_reason,
      auroc = auroc,
      auroc_defined = auroc_defined,
      auroc_reason = auroc_reason,
      ranking_key = ranking_key_desc,
      topk_precision = prec_vals
    ),
    enrichment = list(
      bg_gene_count = as.integer(bg_gene_count),
      n_pairs_bg = as.numeric(n_pairs_bg),
      n_pos_bg = as.integer(n_pos_bg),
      expected_positive_rate = enrich$expected_positive_rate,
      enrichment_ratio = enrich$enrichment_ratio,
      enrichment_p = enrich$enrichment_p,
      top_sets = top_enrich
    )
  )
  if (isTRUE(keep_details)) payload$details <- edge_df
  payload
}

score_row_from_payload <- function(method,
                                   payload,
                                   precision_k,
                                   top_frac,
                                   bench_meta,
                                   top_n,
                                   refill_meta = NULL,
                                   refill_to_top_n_comparable = FALSE,
                                   refill_max_multiplier = NA_integer_,
                                   fdr_pool_size_requested = NA_integer_,
                                   fdr_pool_size_used = NA_integer_,
                                   fdr_threshold = NA_real_,
                                   edges_all_n = NA_integer_,
                                   edges_fdr_n = NA_integer_,
                                   edges_pool_n = NA_integer_) {
  precision_k <- as.integer(precision_k)
  score <- payload$scores
  coverage <- payload$coverage
  enrich <- payload$enrichment

  out <- list(
    method = method,
    top_n = as.integer(top_n),
    n_edges_used = coverage$n_edges,
    n_comparable_edges = coverage$n_comparable_edges,
    edge_coverage = coverage$edge_coverage,
    auprc = score$auprc,
    auprc_defined = score$auprc_defined,
    auprc_reason = score$auprc_reason,
    auroc = score$auroc,
    auroc_defined = score$auroc_defined,
    auroc_reason = score$auroc_reason,
    ranking_key = ranking_key_desc,
    refill_to_top_n_comparable = as.integer(isTRUE(refill_to_top_n_comparable)),
    refill_max_multiplier = as.integer(refill_max_multiplier),
    refill_status = if (!is.null(refill_meta$refill_status)) refill_meta$refill_status else NA_character_,
    refill_shortfall = if (!is.null(refill_meta$refill_shortfall)) as.integer(refill_meta$refill_shortfall) else NA_integer_,
    n_candidates_pre_refill = if (!is.null(refill_meta$n_candidates_pre_refill)) as.integer(refill_meta$n_candidates_pre_refill) else NA_integer_,
    n_mappable_candidates = if (!is.null(refill_meta$n_mappable_candidates)) as.integer(refill_meta$n_mappable_candidates) else NA_integer_,
    fdr_threshold = as.numeric(fdr_threshold),
    fdr_pool_size_requested = as.integer(fdr_pool_size_requested),
    fdr_pool_size_used = as.integer(fdr_pool_size_used),
    edges_all_n = as.integer(edges_all_n),
    edges_fdr_n = as.integer(edges_fdr_n),
    edges_pool_n = as.integer(edges_pool_n)
  )

  for (k in precision_k) {
    key <- paste0("precision_at_", k)
    val <- NA_real_
    if (!is.null(score$topk_precision) && key %in% names(score$topk_precision)) {
      val <- as.numeric(score$topk_precision[[key]])
    }
    out[[paste0("precision@", k)]] <- val
  }

  top_enrich <- NA_real_
  top_n_edges <- NA_integer_
  top_pos <- NA_integer_
  top_pos_rate <- NA_real_
  if (!is.null(enrich) && is.data.frame(enrich$top_sets) && nrow(enrich$top_sets)) {
    target <- paste0("top_frac_", formatC(top_frac * 100, format = "f", digits = 2), "pct")
    row <- enrich$top_sets[enrich$top_sets$set_id == target, , drop = FALSE]
    if (!nrow(row)) {
      row <- enrich$top_sets[grepl("^top_frac_", enrich$top_sets$set_id), , drop = FALSE]
    }
    if (nrow(row)) {
      top_enrich <- row$enrichment_ratio[[1]]
      top_n_edges <- row$n_top[[1]]
      top_pos <- row$n_pos[[1]]
      top_pos_rate <- row$positive_rate[[1]]
    }
  }
  out[["top_frac_1%_enrichment"]] <- top_enrich
  out[["top1p_n_edges"]] <- top_n_edges
  out[["top1p_pos_count"]] <- top_pos
  out[["top_pos_rate"]] <- top_pos_rate

  out[["bg_gene_count"]] <- if (!is.null(enrich$bg_gene_count)) enrich$bg_gene_count else NA_real_
  out[["bg_pair_count"]] <- if (!is.null(enrich$n_pairs_bg)) enrich$n_pairs_bg else NA_real_
  out[["bg_pos_count"]] <- if (!is.null(enrich$n_pos_bg)) enrich$n_pos_bg else NA_real_
  out[["bg_pos_rate"]] <- if (!is.null(enrich$expected_positive_rate)) enrich$expected_positive_rate else NA_real_

  if (!is.null(bench_meta)) {
    out[["n_obs_used"]] <- bench_meta$n_obs_used
    out[["n_vars_used"]] <- bench_meta$n_vars_used
    out[["edges_total_n"]] <- bench_meta$edges_total_n
    out[["edges_sig_n"]] <- bench_meta$edges_sig_n
    out[["edges_top_n"]] <- coverage$n_edges
    out[["edges_weight_type"]] <- bench_meta$edges_weight_type
    out[["edges_pvalue_method"]] <- bench_meta$edges_pvalue_method
  }

  as.data.table(out)
}

plot_positive_pool_string_score_violins <- function(edges_sig_list,
                                                    outdir,
                                                    max_edges_per_method = 200000L,
                                                    fdr_positive_threshold,
                                                    positive_only,
                                                    string_score_threshold,
                                                    predicted_definition = NULL) {
  if (is.null(edges_sig_list) || !length(edges_sig_list)) return(invisible(FALSE))
  rows <- list()
  rows_matched <- list()
  for (method in names(edges_sig_list)) {
    dt <- edges_sig_list[[method]]
    if (is.null(dt) || !nrow(dt) || !"string_score" %in% names(dt)) next
    ss <- suppressWarnings(as.numeric(dt$string_score))
    ss <- ss[is.finite(ss)]
    if (!length(ss)) next
    if (length(ss) > max_edges_per_method) {
      set.seed(1)
      ss <- sample(ss, max_edges_per_method)
    }
    rows[[method]] <- data.table(method = method, string_score = ss)
    ss_m <- ss[ss >= as.numeric(string_score_threshold)]
    if (length(ss_m)) {
      rows_matched[[method]] <- data.table(method = method, string_score = ss_m)
    }
  }
  plot_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(plot_dt)) return(invisible(FALSE))

  plot_dt[, method_label := vapply(method, display_method_name, character(1))]
  color_map <- bench_base_color_map(plot_dt$method_label)
  pred_def <- format_pred_true_definitions(
    fdr_positive_threshold = fdr_positive_threshold,
    positive_only = positive_only,
    string_score_threshold = NA_integer_
  )$predicted

  plot_dir <- file.path(outdir, "plots_positive_pool")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  pred_def <- if (!is.null(predicted_definition) && nzchar(as.character(predicted_definition))) {
    as.character(predicted_definition)
  } else {
    format_pred_true_definitions(
      fdr_positive_threshold = fdr_positive_threshold,
      positive_only = positive_only,
      string_score_threshold = NA_integer_
    )$predicted
  }

  p <- ggplot(plot_dt, aes(x = method_label, y = string_score, fill = method_label)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.6, na.rm = TRUE) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5, na.rm = TRUE) +
    scale_fill_manual(values = color_map) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(
      title = paste0("Positive pool STRING score distribution\n(", pred_def, ")"),
      x = "Method",
      y = "STRING score"
    )
  ggsave(
    filename = file.path(plot_dir, "stringdb_positive_weight_violin.png"),
    plot = p,
    width = 6, height = 4, dpi = 150
  )

  plot_dt_matched <- rbindlist(rows_matched, use.names = TRUE, fill = TRUE)
  if (nrow(plot_dt_matched)) {
    plot_dt_matched[, method_label := vapply(method, display_method_name, character(1))]
    p2 <- ggplot(plot_dt_matched, aes(x = method_label, y = string_score, fill = method_label)) +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.6, na.rm = TRUE) +
      geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5, na.rm = TRUE) +
      scale_fill_manual(values = color_map) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(
        title = paste0(
          "Positive + matched (STRING score \u2265 ", string_score_threshold, ") STRING score distribution\n(",
          pred_def, ")"
        ),
        x = "Method",
        y = "STRING score"
      )
    ggsave(
      filename = file.path(plot_dir, paste0("stringdb_positive_string_score_ge", string_score_threshold, "_violin.png")),
      plot = p2,
      width = 6, height = 4, dpi = 150
    )
  }
  invisible(TRUE)
}

plot_metrics_lines <- function(details_dt, outdir, cut2, fdr_positive_threshold, positive_only, predicted_definition = NULL) {
  if (is.null(details_dt) || !nrow(details_dt)) return(invisible(FALSE))
  if (!all(c("method", "top_n", "string_score") %in% names(details_dt))) return(invisible(FALSE))

  dt <- as.data.table(details_dt)
  dt[, top_n := as.integer(top_n)]
  dt[, string_score := suppressWarnings(as.numeric(string_score))]
  dt <- dt[is.finite(top_n)]
  if (!nrow(dt)) return(invisible(FALSE))

  # Matched to STRING = string_score >= cut2; then bin within matched.
  bin1_lo <- as.numeric(cut2)
  bin1_hi <- min(1000, bin1_lo + 100)
  bin2_lo <- bin1_hi
  bin2_hi <- min(1000, bin2_lo + 100)
  bin3_lo <- bin2_hi
  bin3_hi <- 1000

  agg <- dt[, {
    ss <- string_score
    matched_vec <- is.finite(ss) & ss >= bin1_lo
    matched_n <- sum(matched_vec, na.rm = TRUE)
    if (matched_n <= 0) {
      list(
        matched_edges_n = 0L,
        share_score_1000 = NA_real_,
        share_score_bin3 = NA_real_,
        share_score_bin2 = NA_real_,
        share_score_bin1 = NA_real_
      )
    } else {
      list(
        matched_edges_n = as.integer(matched_n),
        share_score_1000 = sum(ss == 1000 & matched_vec, na.rm = TRUE) / matched_n,
        share_score_bin3 = sum(ss >= bin3_lo & ss < bin3_hi, na.rm = TRUE) / matched_n,
        share_score_bin2 = sum(ss >= bin2_lo & ss < bin2_hi, na.rm = TRUE) / matched_n,
        share_score_bin1 = sum(ss >= bin1_lo & ss < bin1_hi, na.rm = TRUE) / matched_n
      )
    }
  }, by = .(method, top_n)]
  if (!nrow(agg)) return(invisible(FALSE))

  metric_map <- c(
    matched_edges_n = paste0("Matched edges (STRING score \u2265 ", cut2, ")"),
    share_score_1000 = "Share (STRING score = 1000)",
    share_score_bin3 = paste0("Share (", format(bin3_lo, trim = TRUE), "\u2013", format(bin3_hi, trim = TRUE), ")"),
    share_score_bin2 = paste0("Share (", format(bin2_lo, trim = TRUE), "\u2013", format(bin2_hi, trim = TRUE), ")"),
    share_score_bin1 = paste0("Share (", format(bin1_lo, trim = TRUE), "\u2013", format(bin1_hi, trim = TRUE), ")")
  )

  line_long <- melt(agg, id.vars = c("method", "top_n"), variable.name = "metric", value.name = "value")
  if (!nrow(line_long)) return(invisible(FALSE))

  line_long[, method_label := vapply(method, display_method_name, character(1))]
  line_long[, metric := unname(metric_map[metric])]
  color_map <- bench_base_color_map(line_long$method_label)

  x_ticks <- c(10, 100, 1000, 5000)
  x_ticks <- x_ticks[x_ticks >= min(line_long$top_n) & x_ticks <= max(line_long$top_n)]
  x_labels <- as.character(x_ticks)

  defs <- format_pred_true_definitions(
    fdr_positive_threshold = fdr_positive_threshold,
    positive_only = positive_only,
    string_score_threshold = cut2
  )
  if (!is.null(predicted_definition) && nzchar(as.character(predicted_definition))) {
    defs$predicted <- as.character(predicted_definition)
  }

  plot_dir <- file.path(outdir, paste0("plots_recalc_cut", cut2))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  p_lines <- ggplot(line_long, aes(x = top_n, y = value, color = method_label, group = method_label)) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_color_manual(values = color_map) +
    scale_x_log10(breaks = x_ticks, labels = x_labels) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(
      title = paste0("STRING matched-edge count and score bins across Top N\n(", defs$predicted, "; ", defs$truth, ")"),
      x = "Top N (log10)",
      y = "Value",
      color = "Method"
    )
  ggsave(
    filename = file.path(plot_dir, paste0("stringdb_metrics_lines_cut", cut2, ".png")),
    plot = p_lines,
    width = 9, height = 6, dpi = 150
  )
  invisible(TRUE)
}

compute_precision_recall_at_n <- function(edges_sorted, string_score_threshold, n_points, total_true = NA_integer_) {
  edges_sorted <- as.data.table(edges_sorted)
  if (!nrow(edges_sorted)) {
    return(data.table(
      n = integer(0),
      n_positive_total = integer(0),
      total_true = integer(0),
      tp = integer(0),
      fp = integer(0),
      fn = integer(0),
      precision = numeric(0),
      recall = numeric(0)
    ))
  }
  n_positive_total <- nrow(edges_sorted)
  if (n_positive_total == 0) {
    return(data.table(
      n = integer(0),
      n_positive_total = integer(0),
      total_true = integer(0),
      tp = integer(0),
      fp = integer(0),
      fn = integer(0),
      precision = numeric(0),
      recall = numeric(0)
    ))
  }

  edges_sorted[, string_score := suppressWarnings(as.numeric(string_score))]
  true_mask <- is.finite(edges_sorted$string_score) & edges_sorted$string_score >= as.numeric(string_score_threshold)
  true_mask[is.na(true_mask)] <- FALSE
  tp_cum <- cumsum(true_mask)

  total_true <- as.integer(total_true)
  if (!is.finite(total_true) || total_true < 0L) {
    total_true <- sum(true_mask)
  }

  n_points <- as.integer(n_points)
  n_points <- n_points[is.finite(n_points) & n_points >= 1L & n_points <= n_positive_total]
  if (!length(n_points)) {
    n_points <- n_positive_total
  }
  n_points <- unique(n_points)
  n_points <- sort(n_points)

  tp <- tp_cum[n_points]
  fp <- n_points - tp
  fn <- pmax(0L, as.integer(total_true) - tp)
  precision <- tp / n_points
  recall <- if (total_true > 0) tp / total_true else rep(NA_real_, length(tp))

  data.table(
    n = as.integer(n_points),
    n_positive_total = as.integer(n_positive_total),
    total_true = as.integer(total_true),
    tp = as.integer(tp),
    fp = as.integer(fp),
    fn = as.integer(fn),
    precision = as.numeric(precision),
    recall = as.numeric(recall)
  )
}

plot_precision_recall_over_n <- function(edges_sorted_list,
                                        outdir,
                                        string_score_threshold,
                                        method_order,
                                        positive_only,
                                        predicted_definition = NULL,
                                        bin_size = 100L,
                                        total_true_by_method = NULL) {
  if (is.null(edges_sorted_list) || !length(edges_sorted_list)) return(invisible(FALSE))

  method_order <- as.character(method_order)
  bin_size <- as.integer(bin_size)
  if (!is.finite(bin_size) || bin_size <= 0L) bin_size <- 100L

  rows <- list()
  for (method_key in names(edges_sorted_list)) {
    edges_sorted <- edges_sorted_list[[method_key]]
    if (is.null(edges_sorted) || !nrow(edges_sorted)) next
    # Decide per-method N grid on the "positive" pool size (not limited by STRING mapping).
    n_positive_total <- nrow(edges_sorted)
    if (n_positive_total <= 0) next
    n_points <- seq(bin_size, n_positive_total, by = bin_size)
    if (!length(n_points)) n_points <- as.integer(n_positive_total)
    if (tail(n_points, 1) != n_positive_total) n_points <- c(n_points, as.integer(n_positive_total))

    total_true <- NA_integer_
    if (!is.null(total_true_by_method) && method_key %in% names(total_true_by_method)) {
      total_true <- as.integer(total_true_by_method[[method_key]])
    }
    dt <- compute_precision_recall_at_n(edges_sorted, string_score_threshold, n_points = n_points, total_true = total_true)
    if (!nrow(dt)) next
    dt[, method := method_key]
    rows[[method_key]] <- dt
  }
  pr_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(pr_dt)) return(invisible(FALSE))

  pr_dt[, method_label := vapply(method, display_method_name, character(1))]
  if (length(method_order)) {
    method_levels <- vapply(method_order, display_method_name, character(1))
    method_levels <- unique(method_levels) # avoid duplicated display names from different methods
    pr_dt[, method_label := factor(method_label, levels = method_levels)]
  }
  pr_dt[, n_left := data.table::shift(n, fill = 0L), by = method]
  pr_dt[, `:=`(
    precision = as.numeric(precision),
    recall = as.numeric(recall),
    fn = as.numeric(fn)
  )]

  long <- melt(
    pr_dt,
    id.vars = c("method", "method_label", "n", "n_left", "n_positive_total", "total_true", "tp", "fp", "fn"),
    measure.vars = c("precision", "recall", "fn"),
    variable.name = "metric",
    value.name = "value"
  )
  long[, metric := fifelse(
    metric == "precision",
    "Precision (TP/(TP+FP))",
    fifelse(metric == "recall", "Recall (TP/(TP+FN))", "FN (outside Top-N; STRING>=threshold; +/- weight)")
  )]
  long[, metric := factor(metric, levels = c(
    "Precision (TP/(TP+FP))",
    "Recall (TP/(TP+FN))",
    "FN (outside Top-N; STRING>=threshold; +/- weight)"
  ))]

  plot_dir <- file.path(outdir, paste0("plots_recalc_cut", as.integer(string_score_threshold)))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    pr_dt[order(method_label, n)],
    file.path(plot_dir, paste0("stringdb_precision_recall_summary_cut", as.integer(string_score_threshold), ".tsv")),
    sep = "\t"
  )

  color_map <- bench_base_color_map(long$method_label)
  pred_def <- if (!is.null(predicted_definition) && nzchar(as.character(predicted_definition))) {
    as.character(predicted_definition)
  } else {
    format_pred_true_definitions(
      fdr_positive_threshold = NA_real_,
      positive_only = positive_only,
      string_score_threshold = NA_integer_
    )$predicted
  }
  truth_def <- paste0("T: STRING score >= ", as.integer(string_score_threshold))

  p <- ggplot(long, aes(xmin = n_left, xmax = n, ymin = 0, ymax = value, fill = method_label)) +
    geom_rect(alpha = 0.75, na.rm = TRUE) +
    facet_grid(metric ~ method_label, scales = "free") +
    scale_fill_manual(values = color_map) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = paste0("Precision / Recall / FN over Top-N (bin_size=", bin_size, "; All=positive-weight pool)\\n(", pred_def, "; ", truth_def, ")"),
      x = "Top N",
      y = "Value",
      fill = "Method"
    )

  # Keep legacy filename so downstream paths stay stable.
  ggsave(
    filename = file.path(plot_dir, paste0("stringdb_metrics_lines_cut", as.integer(string_score_threshold), ".png")),
    plot = p,
    width = 10, height = 7, dpi = 150
  )
  invisible(TRUE)
}

compute_precision_recall_at_single_n <- function(edges_sorted, string_score_threshold, n_point, total_true = NA_integer_) {
  edges_sorted <- as.data.table(edges_sorted)
  n_point <- as.integer(n_point)
  if (!is.finite(n_point) || n_point < 1L) n_point <- 1L
  if (!nrow(edges_sorted)) {
    return(list(
      n_used = 0L,
      n_positive_total = 0L,
      total_true = as.integer(if (is.finite(total_true)) total_true else 0L),
      tp = 0L,
      fp = 0L,
      fn = as.integer(if (is.finite(total_true)) total_true else 0L),
      precision = NA_real_,
      recall = NA_real_
    ))
  }

  edges_sorted[, string_score := suppressWarnings(as.numeric(string_score))]
  true_mask <- is.finite(edges_sorted$string_score) & edges_sorted$string_score >= as.numeric(string_score_threshold)
  true_mask[is.na(true_mask)] <- FALSE
  tp_cum <- cumsum(true_mask)

  n_positive_total <- nrow(edges_sorted)
  n_used <- min(n_point, n_positive_total)
  tp <- as.integer(tp_cum[[n_used]])

  total_true <- as.integer(total_true)
  if (!is.finite(total_true) || total_true < 0L) {
    total_true <- sum(true_mask)
  }

  fp <- as.integer(n_used - tp)
  fn <- pmax(0L, as.integer(total_true) - tp)
  precision <- if (n_used > 0) tp / n_used else NA_real_
  recall <- if (total_true > 0) tp / total_true else NA_real_

  list(
    n_used = as.integer(n_used),
    n_positive_total = as.integer(n_positive_total),
    total_true = as.integer(total_true),
    tp = as.integer(tp),
    fp = as.integer(fp),
    fn = as.integer(fn),
    precision = as.numeric(precision),
    recall = as.numeric(recall)
  )
}

plot_precision_recall_bars_at_topn <- function(edges_sorted_list,
                                               outdir,
                                               string_score_threshold,
                                               method_order,
                                               positive_only,
                                               predicted_definition = NULL,
                                               top_n_list = c(10L, 100L, 1000L),
                                               write_separate_panels = FALSE,
                                               total_true_by_method = NULL) {
  if (is.null(edges_sorted_list) || !length(edges_sorted_list)) return(invisible(FALSE))

  method_order <- as.character(method_order)
  top_n_list <- as.integer(top_n_list)
  top_n_list <- top_n_list[is.finite(top_n_list) & top_n_list > 0L]
  if (!length(top_n_list)) top_n_list <- c(10L, 100L, 1000L)
  top_n_list <- unique(top_n_list)
  top_n_list <- sort(top_n_list)

  pred_def <- if (!is.null(predicted_definition) && nzchar(as.character(predicted_definition))) {
    as.character(predicted_definition)
  } else {
    format_pred_true_definitions(
      fdr_positive_threshold = NA_real_,
      positive_only = positive_only,
      string_score_threshold = NA_integer_
    )$predicted
  }
  truth_def <- paste0("T: STRING score >= ", as.integer(string_score_threshold))

  plot_dir <- file.path(outdir, paste0("plots_recalc_cut", as.integer(string_score_threshold)))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (top_n in top_n_list) {
    for (method_key in names(edges_sorted_list)) {
      edges_sorted <- edges_sorted_list[[method_key]]
      total_true <- NA_integer_
      if (!is.null(total_true_by_method) && method_key %in% names(total_true_by_method)) {
        total_true <- as.integer(total_true_by_method[[method_key]])
      }
      stats <- compute_precision_recall_at_single_n(edges_sorted, string_score_threshold, n_point = top_n, total_true = total_true)
      rows[[paste(method_key, top_n, sep = "_")]] <- data.table(
        method = method_key,
        top_n = as.integer(top_n),
        n_used = stats$n_used,
        n_positive_total = stats$n_positive_total,
        total_true = stats$total_true,
        tp = stats$tp,
        fp = stats$fp,
        fn = stats$fn,
        precision = stats$precision,
        recall = stats$recall
      )
    }
  }
  pr_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(pr_dt)) return(invisible(FALSE))

  pr_dt[, method_label := vapply(method, display_method_name, character(1))]
  if (length(method_order)) {
    method_levels <- vapply(method_order, display_method_name, character(1))
    method_levels <- unique(method_levels) # avoid duplicated display names from different methods
    pr_dt[, method_label := factor(method_label, levels = method_levels)]
  }

  fwrite(
    pr_dt[order(top_n, method_label)],
    file.path(plot_dir, paste0("stringdb_precision_recall_topN_summary_cut", as.integer(string_score_threshold), ".tsv")),
    sep = "\t"
  )

  color_map <- bench_base_color_map(pr_dt$method_label)

  if (isTRUE(write_separate_panels)) {
    for (top_n_val in top_n_list) {
      # Use a different variable name than the column ("top_n") to avoid data.table
      # treating it as a self-comparison (which would select all rows).
      dt_n <- pr_dt[top_n == as.integer(top_n_val)]
      if (!nrow(dt_n)) next

      p_prec <- ggplot(dt_n, aes(x = method_label, y = precision, fill = method_label)) +
        geom_col(width = 0.75, na.rm = TRUE) +
        scale_fill_manual(values = color_map) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        ) +
        labs(title = paste0("Precision at Top", top_n_val), x = "Method", y = "Precision")

      ggsave(
        filename = file.path(plot_dir, paste0("stringdb_precision_topN", top_n_val, "_cut", as.integer(string_score_threshold), ".png")),
        plot = p_prec,
        width = 6.5, height = 4.5, dpi = 150
      )

      p_rec <- ggplot(dt_n, aes(x = method_label, y = recall, fill = method_label)) +
        geom_col(width = 0.75, na.rm = TRUE) +
        scale_fill_manual(values = color_map) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        ) +
        labs(title = paste0("Recall at Top", top_n_val), x = "Method", y = "Recall")

      ggsave(
        filename = file.path(plot_dir, paste0("stringdb_recall_topN", top_n_val, "_cut", as.integer(string_score_threshold), ".png")),
        plot = p_rec,
        width = 6.5, height = 4.5, dpi = 150
      )
    }
  }

  # Combined 2xK grid (row 1: Precision, row 2: Recall), with fixed y-scale within each row.
  prec_dt <- pr_dt[, .(method_label, top_n, value = precision)]
  prec_dt[, facet_label := paste0("Precision at Top", top_n)]
  prec_dt[, facet_label := factor(facet_label, levels = paste0("Precision at Top", top_n_list))]
  max_prec <- suppressWarnings(max(prec_dt$value, na.rm = TRUE))
  if (!is.finite(max_prec) || max_prec <= 0) max_prec <- 1

  p_prec_row <- ggplot(prec_dt, aes(x = method_label, y = value, fill = method_label)) +
    geom_col(width = 0.75, na.rm = TRUE) +
    scale_fill_manual(values = color_map) +
    facet_wrap(~ facet_label, nrow = 1, scales = "fixed") +
    scale_y_continuous(
      limits = c(0, max_prec * 1.05),
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = NULL, y = "Precision")

  rec_dt <- pr_dt[, .(method_label, top_n, value = recall)]
  rec_dt[, facet_label := paste0("Recall at Top", top_n)]
  rec_dt[, facet_label := factor(facet_label, levels = paste0("Recall at Top", top_n_list))]
  max_rec <- suppressWarnings(max(rec_dt$value, na.rm = TRUE))
  if (!is.finite(max_rec) || max_rec <= 0) max_rec <- 1

  p_rec_row <- ggplot(rec_dt, aes(x = method_label, y = value, fill = method_label)) +
    geom_col(width = 0.75, na.rm = TRUE) +
    scale_fill_manual(values = color_map) +
    facet_wrap(~ facet_label, nrow = 1, scales = "fixed") +
    scale_y_continuous(
      limits = c(0, max_rec * 1.05),
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = "Method", y = "Recall")

  combined_path <- file.path(plot_dir, paste0("stringdb_precision_recall_topN_grid_cut", as.integer(string_score_threshold), ".png"))
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    grob <- gridExtra::arrangeGrob(p_prec_row, p_rec_row, ncol = 1)
    ggsave(filename = combined_path, plot = grob, width = 12, height = 7, dpi = 150)
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    p_grid <- p_prec_row / p_rec_row
    ggsave(filename = combined_path, plot = p_grid, width = 12, height = 7, dpi = 150)
  } else {
    # Fallback: write the two rows separately.
    ggsave(filename = file.path(plot_dir, paste0("stringdb_precision_topN_row_cut", as.integer(string_score_threshold), ".png")),
           plot = p_prec_row, width = 12, height = 3.5, dpi = 150)
    ggsave(filename = file.path(plot_dir, paste0("stringdb_recall_topN_row_cut", as.integer(string_score_threshold), ".png")),
           plot = p_rec_row, width = 12, height = 3.5, dpi = 150)
  }

  invisible(TRUE)
}

plot_precision_recall_curves_over_topn_list <- function(edges_sorted_list,
                                                        outdir,
                                                        string_score_threshold,
                                                        method_order,
                                                        positive_only,
                                                        predicted_definition = NULL,
                                                        top_n_list,
                                                        total_true_by_method = NULL) {
  if (is.null(edges_sorted_list) || !length(edges_sorted_list)) return(invisible(FALSE))

  top_n_list <- as.integer(top_n_list)
  top_n_list <- top_n_list[is.finite(top_n_list) & top_n_list > 0L]
  if (!length(top_n_list)) return(invisible(FALSE))
  top_n_list <- unique(top_n_list)
  top_n_list <- sort(top_n_list)

  pred_def <- if (!is.null(predicted_definition) && nzchar(as.character(predicted_definition))) {
    as.character(predicted_definition)
  } else {
    format_pred_true_definitions(
      fdr_positive_threshold = NA_real_,
      positive_only = positive_only,
      string_score_threshold = NA_integer_
    )$predicted
  }
  truth_def <- paste0("T: STRING score >= ", as.integer(string_score_threshold))

  plot_dir <- file.path(outdir, paste0("plots_recalc_cut", as.integer(string_score_threshold)))
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  method_order <- as.character(method_order)
  rows <- list()
  for (method_key in names(edges_sorted_list)) {
    edges_sorted <- edges_sorted_list[[method_key]]
    if (is.null(edges_sorted) || !nrow(edges_sorted)) next

    total_true <- NA_integer_
    if (!is.null(total_true_by_method) && method_key %in% names(total_true_by_method)) {
      total_true <- as.integer(total_true_by_method[[method_key]])
    }

    dt <- compute_precision_recall_at_n(
      edges_sorted = edges_sorted,
      string_score_threshold = string_score_threshold,
      n_points = top_n_list,
      total_true = total_true
    )
    if (!nrow(dt)) next
    dt[, method := method_key]
    rows[[method_key]] <- dt
  }
  pr_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(pr_dt)) return(invisible(FALSE))

  pr_dt[, top_n := as.integer(n)]
  pr_dt[, n := NULL]
  pr_dt[, method_label := vapply(method, display_method_name, character(1))]
  if (length(method_order)) {
    method_levels <- vapply(method_order, display_method_name, character(1))
    method_levels <- unique(method_levels)
    pr_dt[, method_label := factor(method_label, levels = method_levels)]
  }

  fwrite(
    pr_dt[order(method_label, top_n)],
    file.path(plot_dir, paste0("stringdb_precision_recall_topN_curve_summary_cut", as.integer(string_score_threshold), ".tsv")),
    sep = "\t"
  )

  long <- melt(
    pr_dt,
    id.vars = c("method", "method_label", "top_n", "n_positive_total", "total_true", "tp", "fp", "fn"),
    measure.vars = c("precision", "recall"),
    variable.name = "metric",
    value.name = "value"
  )
  long[, metric := fifelse(metric == "precision", "Precision", "Recall")]
  long[, metric := factor(metric, levels = c("Precision", "Recall"))]

  x_ticks <- c(10, 100, 1000, 5000)
  x_ticks <- x_ticks[x_ticks >= min(long$top_n) & x_ticks <= max(long$top_n)]
  x_labels <- as.character(x_ticks)

  color_map <- bench_base_color_map(long$method_label)

  p <- ggplot(long, aes(x = top_n, y = value, color = method_label, group = method_label)) +
    geom_line(linewidth = 0.65, na.rm = TRUE) +
    facet_wrap(~ metric, ncol = 1, scales = "fixed") +
    scale_color_manual(values = color_map) +
    scale_x_log10(breaks = x_ticks, labels = x_labels) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(
      title = paste0("Precision / Recall across Top N\n(", pred_def, "; ", truth_def, ")"),
      x = "Top N (log10)",
      y = "Value",
      color = "Method"
    )

  ggsave(
    filename = file.path(plot_dir, paste0("stringdb_precision_recall_topN_curve_cut", as.integer(string_score_threshold), ".png")),
    plot = p,
    width = 8, height = 6, dpi = 150
  )

  invisible(TRUE)
}

option_list <- list(
  make_option(c("--map_dir"), type = "character", help = "Directory produced by stringdb_edge_map_task83_style.R"),
  make_option(c("--outdir"), type = "character", help = "Output directory (will contain rank_*/ like the one-step script)"),
  make_option(c("--methods"), type = "character", default = "genescope,giotto,hotspot"),
  make_option(c("--ranking_keys"), type = "character", default = "weight,abs_weight"),
  make_option(c("--string_score_threshold"), type = "integer", default = 700),
  make_option(c("--precision_k"), type = "character", default = "50,100,500"),
  make_option(c("--top_frac"), type = "double", default = 0.01),
  make_option(c("--top_n_list"), type = "character", default = "10,20,30,40,50,100,200,300,400,500,1000,2000,3000,4000,5000"),
  make_option(c("--edge_fdr"), type = "double", default = NA_real_),
  make_option(c("--edge_fdr_by_method"), type = "character", default = NA_character_,
              help = "Comma-separated FDR thresholds aligned with --methods order (e.g. '0.05,0.05,1'). Values >=1 or NA disable FDR filtering for that method. Overrides --edge_fdr and meta.json edges_fdr_threshold."),
  make_option(c("--fdr_pool_size"), type = "integer", default = NA_integer_),
  make_option(c("--positive_weights_only"), type = "integer", default = 1),
  make_option(c("--weight_scale"), type = "character", default = "quantile"),
  make_option(c("--weight_scale_bins"), type = "integer", default = 1000),
  make_option(c("--refill_to_top_n_comparable"), type = "integer", default = 1),
  make_option(c("--refill_max_multiplier"), type = "integer", default = 5),
  make_option(c("--keep_details"), type = "integer", default = 1),
  make_option(c("--pr_bin_size"), type = "integer", default = 100L,
              help = "Bin size for Top-N grid when plotting Precision/Recall/FN curve (e.g., 100). Script always includes the final All point per method."),
  make_option(c("--pr_top_n_list"), type = "character", default = "10,100,1000",
              help = "Comma-separated Top-N values for bar plots (Precision/Recall), used in the combined grid and (optionally) separate single-panel images."),
  make_option(c("--pr_write_separate_panels"), type = "integer", default = 0L,
              help = "If 1, also write separate single-panel bar plots per Top-N (stringdb_precision_topN*_cut*.png and stringdb_recall_topN*_cut*.png). Default: 0 (off).")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$map_dir) || is.null(opt$outdir)) stop("Missing --map_dir or --outdir")

mapped_tsv <- file.path(opt$map_dir, "stringdb_edge_mapped_all.tsv")
meta_json <- file.path(opt$map_dir, "stringdb_edge_mapped_meta.json")
if (!file.exists(mapped_tsv)) stop("Missing mapped edges TSV: ", mapped_tsv)
if (!file.exists(meta_json)) stop("Missing mapped meta JSON: ", meta_json)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

mapped_dt <- fread(mapped_tsv)
meta <- read_json(meta_json, simplifyVector = TRUE)
methods <- split_csv(opt$methods)
if (!length(methods)) stop("No methods provided")

edge_fdr_by_method_vec <- NULL
if (!is.null(opt$edge_fdr_by_method) && !is.na(opt$edge_fdr_by_method) && nzchar(trimws(opt$edge_fdr_by_method))) {
  vals <- parse_double_list_allow_na(opt$edge_fdr_by_method)
  if (length(vals) != length(methods)) {
    stop("--edge_fdr_by_method must have the same number of values as --methods (methods=", length(methods), ", fdr=", length(vals), ")")
  }
  bad <- which(is.finite(vals) & vals < 0)
  if (length(bad)) {
    stop("Invalid negative FDR threshold(s) in --edge_fdr_by_method: ", paste(bad, collapse = ", "))
  }
  edge_fdr_by_method_vec <- stats::setNames(vals, methods)
}

ranking_keys <- parse_ranking_keys(opt$ranking_keys)
if (!length(ranking_keys)) stop("No valid ranking keys provided")

top_n_list <- parse_top_n_list(opt$top_n_list)
if (!length(top_n_list)) stop("No valid --top_n_list provided")

precision_k <- as.integer(split_csv(opt$precision_k))
precision_k <- precision_k[is.finite(precision_k) & precision_k > 0L]
if (!length(precision_k)) precision_k <- c(50L, 100L, 500L)

positive_only <- opt$positive_weights_only != 0
weight_scale <- tolower(trimws(opt$weight_scale))
if (!nzchar(weight_scale)) weight_scale <- "none"
weight_scale_bins <- as.integer(opt$weight_scale_bins)
if (!is.finite(weight_scale_bins) || weight_scale_bins < 10L) weight_scale_bins <- 1000L

refill_to_top_n_comparable <- opt$refill_to_top_n_comparable != 0
refill_max_multiplier <- as.integer(opt$refill_max_multiplier)
if (!is.finite(refill_max_multiplier) || refill_max_multiplier < 1L) refill_max_multiplier <- 5L

keep_details_flag <- opt$keep_details != 0

per_method_meta <- meta$per_method
if (is.null(per_method_meta) || !length(per_method_meta)) stop("Invalid meta JSON: missing per_method")

pr_bin_size <- as.integer(opt$pr_bin_size)
if (!is.finite(pr_bin_size) || pr_bin_size <= 0L) pr_bin_size <- 100L
pr_top_n_list <- parse_top_n_list(opt$pr_top_n_list)
if (!length(pr_top_n_list)) pr_top_n_list <- c(10L, 100L, 1000L)
pr_write_separate_panels <- opt$pr_write_separate_panels != 0

run_one_ranking_key <- function(ranking_key_desc_use) {
  ranking_key_desc <<- ranking_key_desc_use
  rk_slug <- ranking_key_slug(ranking_key_desc_use)
  branch_outdir <- file.path(opt$outdir, paste0("rank_", rk_slug))
  dir.create(branch_outdir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("[rank_%s] Start: methods=%d, top_n=%d", rk_slug, length(methods), length(top_n_list)))

  edges_sig_list <- list()
  edges_pool_raw_list <- list()
  edges_sorted_list <- list()
  rows <- list()
  details_list <- list()
  meta_list <- list()

  method_progress <- make_percent_progress_logger(length(methods), label = sprintf("[rank_%s] Load methods", rk_slug), step_percent = 25L)
  method_i <- 0L
  for (method_key in methods) {
    method_i <- method_i + 1L
    method_progress(method_i, detail = paste0("method=", method_key))
    required <- c("from", "to", "weight", "fdr", "string_score")
    missing <- setdiff(required, names(mapped_dt))
    if (length(missing)) stop("Mapped edges file missing columns: ", paste(missing, collapse = ", "))
    subscore_cols <- c("nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore")
    cols_use <- c("from", "to", "weight", "fdr")
    if ("string_from" %in% names(mapped_dt)) cols_use <- c(cols_use, "string_from")
    if ("string_to" %in% names(mapped_dt)) cols_use <- c(cols_use, "string_to")
    cols_use <- c(cols_use, "string_score", intersect(subscore_cols, names(mapped_dt)))
    dt <- mapped_dt[method == method_key, cols_use, with = FALSE]
    if (!nrow(dt)) stop("No mapped edges for method=", method_key, " in ", mapped_tsv)
    if (!"string_from" %in% names(dt)) dt[, string_from := NA_character_]
    if (!"string_to" %in% names(dt)) dt[, string_to := NA_character_]
    data.table::setcolorder(dt, c("from", "to", "weight", "fdr", "string_from", "string_to", "string_score", intersect(subscore_cols, names(dt))))
    if (!is.character(dt$from)) dt[, from := as.character(from)]
    if (!is.character(dt$to)) dt[, to := as.character(to)]
    if (!is.character(dt$string_from)) dt[, string_from := as.character(string_from)]
    if (!is.character(dt$string_to)) dt[, string_to := as.character(string_to)]
    if (!is.numeric(dt$weight)) dt[, weight := suppressWarnings(as.numeric(weight))]
    if (!is.numeric(dt$fdr)) dt[, fdr := suppressWarnings(as.numeric(fdr))]
    if (!is.numeric(dt$string_score)) dt[, string_score := suppressWarnings(as.numeric(string_score))]

    bench_meta <- NULL
    bg <- NULL
    string_version <- NA_character_
    if (!is.null(per_method_meta[[method_key]])) {
      bench_meta <- per_method_meta[[method_key]]$bench_meta
      bg <- per_method_meta[[method_key]]$bg
      string_version <- per_method_meta[[method_key]]$string_version
    }

    fdr_threshold <- NA_real_
    if (!is.null(edge_fdr_by_method_vec) && method_key %in% names(edge_fdr_by_method_vec)) {
      fdr_threshold <- edge_fdr_by_method_vec[[method_key]]
    } else {
      fdr_threshold <- opt$edge_fdr
      if (!is.finite(fdr_threshold) && !is.null(bench_meta) && !is.null(bench_meta$edges_fdr_threshold)) {
        fdr_threshold <- suppressWarnings(as.numeric(bench_meta$edges_fdr_threshold))
      }
      if (!is.finite(fdr_threshold) || fdr_threshold <= 0) fdr_threshold <- 0.05
    }
    if (is.finite(fdr_threshold) && fdr_threshold < 0) {
      stop("Invalid FDR threshold for method=", method_key, ": ", fdr_threshold)
    }

    fdr_filtering_enabled <- is.finite(fdr_threshold) && fdr_threshold > 0 && fdr_threshold < 1
    # Truth universe for FN/recall: all edges (positive or negative weights), regardless of FDR.
    # Predicted positives are still restricted by the positive pool (positive weights + optional FDR filter).
    total_true_universe <- dt[
      is.finite(weight) & is.finite(string_score) & string_score >= as.numeric(opt$string_score_threshold),
      .N
    ]

    edges_sig <- if (isTRUE(fdr_filtering_enabled)) {
      filter_edges_for_scoring(dt, fdr_threshold, positive_only)
    } else {
      out <- dt[is.finite(weight)]
      if (isTRUE(positive_only)) out <- out[weight > 0]
      out
    }
    edges_sig_list[[method_key]] <- edges_sig
    meta_list[[method_key]] <- list(
      bench_meta = bench_meta,
      bg = bg,
      string_version = string_version,
      fdr_threshold = fdr_threshold,
      fdr_filtering_enabled = isTRUE(fdr_filtering_enabled),
      total_true_universe = as.integer(total_true_universe),
      n_edges_all = nrow(dt),
      n_edges_sig = nrow(edges_sig),
      n_edges_pool = nrow(edges_sig)
    )
  }

  fdr_thresholds_used <- vapply(methods, function(m) meta_list[[m]]$fdr_threshold, numeric(1))
  predicted_def_str <- format_predicted_definition_by_method(methods, fdr_thresholds_used, positive_only = positive_only)
  plot_positive_pool_string_score_violins(
    edges_sig_list = edges_sig_list,
    outdir = branch_outdir,
    fdr_positive_threshold = opt$edge_fdr,
    positive_only = positive_only,
    string_score_threshold = as.integer(opt$string_score_threshold),
    predicted_definition = predicted_def_str
  )

  fdr_pool_size <- as.integer(opt$fdr_pool_size)
  if (!is.finite(fdr_pool_size) || fdr_pool_size <= 0L) fdr_pool_size <- NA_integer_
  fdr_pool_size_used <- NA_integer_
  if (is.finite(fdr_pool_size)) {
    sig_sizes <- vapply(edges_sig_list, nrow, integer(1))
    if (!length(sig_sizes)) stop("No edges after FDR filtering.")
    fdr_pool_size_used <- min(fdr_pool_size, min(sig_sizes))
  }

  for (method_key in methods) {
    pool_dt <- edges_sig_list[[method_key]]
    if (is.finite(fdr_pool_size_used)) {
      pool_dt <- subset_edges_fdr_pool(pool_dt, fdr_pool_size_used)
    }
    edges_pool_raw_list[[method_key]] <- pool_dt
    meta_list[[method_key]]$n_edges_pool <- nrow(pool_dt)
  }

  if (weight_scale != "none") {
    message(sprintf("[rank_%s] Weight scaling: %s (bins=%d)", rk_slug, weight_scale, as.integer(weight_scale_bins)))
    weight_list <- lapply(edges_pool_raw_list, function(df) df$weight)
    ref <- build_quantile_reference(weight_list, n_bins = weight_scale_bins)
    ws_progress <- make_percent_progress_logger(length(methods), label = sprintf("[rank_%s] Weight scaling", rk_slug), step_percent = 25L)
    ws_i <- 0L
    for (method_key in methods) {
      ws_i <- ws_i + 1L
      ws_progress(ws_i, detail = paste0("method=", method_key))
      df <- edges_pool_raw_list[[method_key]]
      if (!is.null(df) && nrow(df)) {
        df$weight <- apply_quantile_map(df$weight, ref)
      }
      edges_pool_raw_list[[method_key]] <- df
    }
  }

  sort_progress <- make_percent_progress_logger(length(methods), label = sprintf("[rank_%s] Sort edges", rk_slug), step_percent = 25L)
  sort_i <- 0L
  for (method_key in methods) {
    sort_i <- sort_i + 1L
    sort_progress(sort_i, detail = paste0("method=", method_key))
    edges_sorted_list[[method_key]] <- sort_edges_for_ranking(edges_pool_raw_list[[method_key]])
  }

  total_true_by_method <- vapply(methods, function(m) meta_list[[m]]$total_true_universe, integer(1))
  names(total_true_by_method) <- methods

  n_score_tasks <- length(methods) * length(top_n_list)
  step_pct <- if (n_score_tasks <= 200) 1L else 5L
  score_progress <- make_percent_progress_logger(
    n_score_tasks,
    label = sprintf("[rank_%s] Scoring grid", rk_slug),
    step_percent = step_pct
  )
  score_i <- 0L

  # Bar plots at requested Top-N points (default: 10/100/1000).
  plot_precision_recall_bars_at_topn(
    edges_sorted_list = edges_sorted_list,
    outdir = branch_outdir,
    string_score_threshold = as.integer(opt$string_score_threshold),
    method_order = methods,
    positive_only = positive_only,
    predicted_definition = predicted_def_str,
    top_n_list = pr_top_n_list,
    write_separate_panels = isTRUE(pr_write_separate_panels),
    total_true_by_method = total_true_by_method
  )

  # Curves across the full benchmarking Top-N list.
  plot_precision_recall_curves_over_topn_list(
    edges_sorted_list = edges_sorted_list,
    outdir = branch_outdir,
    string_score_threshold = as.integer(opt$string_score_threshold),
    method_order = methods,
    positive_only = positive_only,
    predicted_definition = predicted_def_str,
    top_n_list = top_n_list,
    total_true_by_method = total_true_by_method
  )

  for (method_key in methods) {
    bench_meta <- meta_list[[method_key]]$bench_meta
    bg <- meta_list[[method_key]]$bg
    bg_gene_count <- if (!is.null(bg) && !is.null(bg$bg_gene_count)) as.integer(bg$bg_gene_count) else NA_integer_
    n_pairs_bg <- if (!is.null(bg) && !is.null(bg$n_pairs_bg)) as.numeric(bg$n_pairs_bg) else NA_real_
    n_pos_bg <- NA_integer_
    if (!is.null(bg) && !is.null(bg$n_pos_bg_by_threshold)) {
      key <- as.character(as.integer(opt$string_score_threshold))
      if (!is.null(bg$n_pos_bg_by_threshold[[key]])) {
        n_pos_bg <- as.integer(bg$n_pos_bg_by_threshold[[key]])
      }
    }

    edges_sorted <- edges_sorted_list[[method_key]]
    for (top_n in top_n_list) {
      score_i <- score_i + 1L
      score_progress(score_i, detail = paste0("method=", method_key, ", top_n=", top_n))
      refill_meta <- select_edges_with_refill(
        edges_sorted = edges_sorted,
        top_n = top_n,
        refill_to_top_n_comparable = refill_to_top_n_comparable,
        refill_max_multiplier = refill_max_multiplier
      )
      edges_eval <- refill_meta$edges_eval
      payload <- score_edges_from_mapped(
        edge_df = edges_eval,
        string_score_threshold = as.integer(opt$string_score_threshold),
        precision_k = precision_k,
        top_frac = opt$top_frac,
        bg_gene_count = bg_gene_count,
        n_pairs_bg = n_pairs_bg,
        n_pos_bg = n_pos_bg,
        keep_details = keep_details_flag
      )

      rows[[paste(method_key, top_n, sep = "_")]] <- score_row_from_payload(
        method = method_key,
        payload = payload,
        precision_k = precision_k,
        top_frac = opt$top_frac,
        bench_meta = bench_meta,
        top_n = top_n,
        refill_meta = refill_meta,
        refill_to_top_n_comparable = refill_to_top_n_comparable,
        refill_max_multiplier = refill_max_multiplier,
        fdr_pool_size_requested = fdr_pool_size,
        fdr_pool_size_used = fdr_pool_size_used,
        fdr_threshold = meta_list[[method_key]]$fdr_threshold,
        edges_all_n = meta_list[[method_key]]$n_edges_all,
        edges_fdr_n = meta_list[[method_key]]$n_edges_sig,
        edges_pool_n = meta_list[[method_key]]$n_edges_pool
      )

      if (keep_details_flag && !is.null(payload$details) && nrow(payload$details)) {
        detail_dt <- as.data.table(payload$details)
        detail_dt[, method := method_key]
        detail_dt[, top_n := as.integer(top_n)]
        details_list[[paste(method_key, top_n, sep = "_")]] <- detail_dt
      }
    }
  }

  out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  fwrite(paper_rename_summary_cols(out_dt), file.path(branch_outdir, "stringdb_edge_scores.tsv"), sep = "\t")

  details_dt_out <- NULL
  if (length(details_list)) {
    details_dt <- rbindlist(details_list, use.names = TRUE, fill = TRUE)
    fwrite(paper_rename_details_cols(details_dt), file.path(branch_outdir, "stringdb_edge_scores_details.tsv"), sep = "\t")
    details_dt_out <- details_dt
  }

  if (!is.null(details_dt_out)) {
    for (cut in unique(c(400L, as.integer(opt$string_score_threshold)))) {
      if (!is.finite(cut) || cut <= 0L) next
      dt_cut <- details_dt_out[is.finite(string_score) & string_score >= cut]
      fwrite(
        paper_rename_details_cols(dt_cut),
        file.path(branch_outdir, paste0("stringdb_edge_scores_details_ge", cut, ".tsv")),
        sep = "\t"
      )
    }
  }

  run_meta <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    map_dir = opt$map_dir,
    outdir = opt$outdir,
    methods = methods,
    params = list(
      string_score_threshold = as.integer(opt$string_score_threshold),
      precision_k = precision_k,
      top_frac = opt$top_frac,
      top_n_list = top_n_list,
      edge_fdr = opt$edge_fdr,
      edge_fdr_by_method = opt$edge_fdr_by_method,
      fdr_threshold_by_method_used = stats::setNames(as.list(fdr_thresholds_used), methods),
      fdr_pool_size_requested = fdr_pool_size,
      fdr_pool_size_used = fdr_pool_size_used,
      positive_weights_only = positive_only,
      weight_scale = weight_scale,
      weight_scale_bins = weight_scale_bins,
      ranking_key = ranking_key_desc_use,
      ranking_key_slug = rk_slug,
      refill_to_top_n_comparable = isTRUE(refill_to_top_n_comparable),
      refill_max_multiplier = as.integer(refill_max_multiplier),
      keep_details = keep_details_flag,
      pr_bin_size = as.integer(pr_bin_size),
      pr_top_n_list = pr_top_n_list,
      pr_write_separate_panels = isTRUE(pr_write_separate_panels)
    )
  )
  write_json(run_meta, file.path(branch_outdir, "stringdb_edge_run_meta.json"), auto_unbox = TRUE, pretty = TRUE)
  message(sprintf("[rank_%s] Done", rk_slug))
  invisible(TRUE)
}

for (rk in ranking_keys) {
  run_one_ranking_key(rk)
}

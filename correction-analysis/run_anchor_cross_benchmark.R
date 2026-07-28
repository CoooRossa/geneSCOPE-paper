#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(optparse)
})

options(stringsAsFactors = FALSE)

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (!length(hit)) return(normalizePath(getwd()))
  dirname(normalizePath(sub("^--file=", "", hit[[1]])))
}

sha256_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  sha <- Sys.which("sha256sum")
  if (nzchar(sha)) {
    value <- system2(sha, shQuote(path), stdout = TRUE, stderr = TRUE)
  } else {
    sha <- Sys.which("shasum")
    if (!nzchar(sha)) stop("Neither sha256sum nor shasum is available")
    value <- system2(sha, c("-a", "256", shQuote(path)),
                     stdout = TRUE, stderr = TRUE)
  }
  status <- attr(value, "status", exact = TRUE)
  digest <- if (length(value)) strsplit(value[[1]], "[[:space:]]+")[[1]][[1]] else ""
  if ((!is.null(status) && status != 0L) || !grepl("^[0-9a-fA-F]{64}$", digest)) {
    stop("Could not hash: ", path)
  }
  tolower(digest)
}

package_version_or_na <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}

assert_identical <- function(observed, expected, label) {
  if (!identical(as.character(observed), as.character(expected))) {
    stop(label, " mismatch: observed=", observed, ", expected=", expected)
  }
}

assert_map_meta_semantics <- function(current_meta_path, frozen_map_meta, sample) {
  if (!file.exists(current_meta_path) || is.null(frozen_map_meta)) {
    stop(sample, " map SHA-256 drifted and semantic map metadata are unavailable")
  }
  current <- fromJSON(current_meta_path, simplifyVector = FALSE)
  assert_identical(current$string_version, frozen_map_meta$string_version,
                   paste(sample, "STRING version"))
  # Giotto was later restaged with a larger input universe in the LN tree. The
  # other three maps span the published union and therefore provide the stable
  # STRING projection used by the module benchmark.
  fields <- list(
    c("mapping", "n_input_genes"), c("mapping", "n_mapped"),
    c("bg", "bg_gene_count"), c("bg", "n_pairs_bg"),
    c("bg", "n_pos_bg_by_threshold", "400"),
    c("bg", "n_pos_bg_by_threshold", "700"),
    c("mapped_edges", "n_edges"), c("mapped_edges", "n_comparable_edges")
  )
  descend <- function(x, path) {
    for (name in path) x <- x[[name]]
    x
  }
  for (method in c("genescope", "hotspot", "seagal")) {
    for (field in fields) {
      assert_identical(
        descend(current$per_method[[method]], field),
        descend(frozen_map_meta$per_method[[method]], field),
        paste(sample, method, paste(field, collapse = "."))
      )
    }
  }
  invisible(TRUE)
}

clean_gene <- function(x) trimws(as.character(x))

read_modules_checked <- function(path) {
  x <- fread(path)
  missing <- setdiff(c("gene", "module_id"), names(x))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  x <- x[, .(gene = clean_gene(gene), module_id = as.character(module_id))]
  x <- x[nzchar(gene)]
  unique(x, by = "gene")
}

stage_genescope_modules <- function(old_modules_path, new_membership_path, output_path) {
  old <- read_modules_checked(old_modules_path)
  old[, module_id := "-1"]

  new <- fread(new_membership_path)
  missing <- setdiff(c("gene", "shuffle_module"), names(new))
  if (length(missing)) {
    stop("Missing columns in ", new_membership_path, ": ", paste(missing, collapse = ", "))
  }
  new <- new[, .(
    gene = clean_gene(gene),
    module_id = trimws(as.character(shuffle_module))
  )]
  new <- new[nzchar(gene)]
  new[is.na(module_id) | !nzchar(module_id) | module_id == "NA", module_id := "-1"]
  if (new[, anyDuplicated(gene)]) stop("Duplicate genes in ", new_membership_path)

  assigned_new <- new[module_id != "-1"]
  missing_assigned <- assigned_new[!gene %in% old$gene]
  if (nrow(missing_assigned)) old <- rbind(old, missing_assigned, use.names = TRUE)
  old[new, on = "gene", module_id := i.module_id]

  if (old[module_id != "-1", .N] != assigned_new[, .N]) {
    stop("Assigned-gene count changed while staging ", new_membership_path)
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(old, output_path, sep = "\t")
  invisible(old)
}

stage_background_only <- function(input_path, output_path) {
  x <- read_modules_checked(input_path)
  x[, module_id := "-1"]
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, output_path, sep = "\t")
  invisible(x)
}

ensure_symlink <- function(target, link) {
  target <- normalizePath(target)
  dir.create(dirname(link), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(link)) {
    if (!identical(normalizePath(link), target)) {
      stop("Existing staged link has the wrong target: ", link)
    }
    return(invisible(link))
  }
  if (!isTRUE(file.symlink(target, link))) stop("Failed to link ", link, " -> ", target)
  invisible(link)
}

rebuild_ln_map_if_needed <- function(tree_root, compare_root, sample_out, old_meta,
                                     mapping_script, rscript, current_map_path) {
  expected_map_sha256 <- as.character(old_meta$map_sha256)
  current_map_sha256 <- sha256_file(current_map_path)
  if (identical(current_map_sha256, expected_map_sha256)) {
    return(list(
      map_dir = dirname(current_map_path),
      rebuilt = FALSE,
      original_map_sha256 = current_map_sha256,
      rebuilt_map_sha256 = current_map_sha256
    ))
  }

  methods <- c("genescope", "giotto", "hotspot", "seagal")
  edge_sources <- setNames(
    file.path(compare_root, methods, "repeat_001", "edges_all.tsv"), methods
  )
  edge_sources[["giotto"]] <- file.path(tree_root, "giotto-modules", "edges_all.tsv")
  missing <- edge_sources[!file.exists(edge_sources)]
  if (length(missing)) stop("Missing LN frozen map inputs: ", paste(missing, collapse = ", "))
  for (method in methods) {
    assert_identical(
      sha256_file(edge_sources[[method]]),
      old_meta$map_meta$per_method[[method]]$edges_all_sha256,
      paste("LN", method, "frozen edges hash")
    )
  }

  rebuild_root <- file.path(sample_out, "rebuilt_frozen_string_map")
  staged_bench <- file.path(rebuild_root, "benchmark_inputs")
  rebuilt_map_dir <- file.path(rebuild_root, "map")
  dir.create(rebuilt_map_dir, recursive = TRUE, showWarnings = FALSE)
  for (method in methods) {
    rep_dir <- file.path(staged_bench, method, "repeat_001")
    ensure_symlink(edge_sources[[method]], file.path(rep_dir, "edges_all.tsv"))
    source_meta <- file.path(compare_root, method, "repeat_001", "meta.json")
    if (file.exists(source_meta)) ensure_symlink(source_meta, file.path(rep_dir, "meta.json"))
  }

  cache_source <- file.path(tree_root, "string1step", "STRINGdb_cache_v12.0")
  if (!dir.exists(cache_source)) stop("Missing frozen STRING v12 cache: ", cache_source)
  cache_source_files <- list.files(cache_source, full.names = TRUE)
  cache_source_files <- cache_source_files[file.info(cache_source_files)$isdir %in% FALSE]
  if (!length(cache_source_files) || any(grepl("[.]rds$", cache_source_files, ignore.case = TRUE))) {
    stop("Frozen STRING source cache is missing or contains generated RDS files: ", cache_source)
  }
  local_cache <- file.path(rebuilt_map_dir, "STRINGdb_cache_v12.0")
  legacy_cache_link <- Sys.readlink(local_cache)
  if (nzchar(legacy_cache_link)) {
    if (!identical(normalizePath(local_cache), normalizePath(cache_source))) {
      stop("Existing cache-directory link has the wrong target: ", local_cache)
    }
    unlink(local_cache)
  }
  dir.create(local_cache, recursive = TRUE, showWarnings = FALSE)
  for (source_file in cache_source_files) {
    ensure_symlink(source_file, file.path(local_cache, basename(source_file)))
  }

  rebuilt_map_path <- file.path(rebuilt_map_dir, "stringdb_edge_mapped_all.tsv")
  rebuild_log <- file.path(rebuild_root, "mapping.log")
  if (!file.exists(rebuilt_map_path) ||
      !identical(sha256_file(rebuilt_map_path), expected_map_sha256)) {
    message("[LN] rebuilding frozen STRING map with giotto-modules override")
    args <- c(
      mapping_script,
      "--bench_root", staged_bench,
      "--outdir", rebuilt_map_dir,
      "--methods", paste(methods, collapse = ","),
      "--species", "9606",
      "--input_id_type", "gene",
      "--string_version", "12.0",
      "--string_score_thresholds", "400,700",
      "--keep_subscores", "1"
    )
    status <- system2(rscript, args, stdout = rebuild_log, stderr = rebuild_log)
    if (!identical(status, 0L)) stop("LN map rebuild failed; see ", rebuild_log)
  }
  rebuilt_map_sha256 <- sha256_file(rebuilt_map_path)
  assert_identical(rebuilt_map_sha256, expected_map_sha256, "LN rebuilt STRING map hash")

  generated_cache_files <- list.files(local_cache, full.names = TRUE)
  generated_cache_files <- generated_cache_files[
    file.info(generated_cache_files)$isdir %in% FALSE &
      !nzchar(Sys.readlink(generated_cache_files))
  ]
  generated_cache_archive <- file.path(rebuild_root, "generated_rds_cache")
  if (length(generated_cache_files)) {
    dir.create(generated_cache_archive, recursive = TRUE, showWarnings = FALSE)
    for (generated_file in generated_cache_files) {
      archive_file <- file.path(generated_cache_archive, basename(generated_file))
      if (file.exists(archive_file)) {
        assert_identical(
          sha256_file(generated_file), sha256_file(archive_file),
          paste("LN generated STRING cache", basename(generated_file))
        )
        unlink(generated_file)
      } else if (!file.rename(generated_file, archive_file)) {
        stop("Could not archive generated STRING cache file: ", generated_file)
      }
    }
  }
  archived_cache_files <- if (dir.exists(generated_cache_archive)) {
    files <- list.files(generated_cache_archive, full.names = TRUE)
    files[file.info(files)$isdir %in% FALSE]
  } else {
    character()
  }
  list(
    map_dir = rebuilt_map_dir,
    rebuilt = TRUE,
    original_map_path = current_map_path,
    original_map_sha256 = current_map_sha256,
    rebuilt_map_path = rebuilt_map_path,
    rebuilt_map_sha256 = rebuilt_map_sha256,
    expected_map_sha256 = expected_map_sha256,
    giotto_override_path = edge_sources[["giotto"]],
    giotto_override_sha256 = sha256_file(edge_sources[["giotto"]]),
    edge_source_sha256 = lapply(edge_sources, sha256_file),
    cache_source = cache_source,
    cache_source_file_sha256 = setNames(
      lapply(cache_source_files, sha256_file), basename(cache_source_files)
    ),
    generated_rds_cache = list(
      path = generated_cache_archive,
      file_sha256 = setNames(
        lapply(archived_cache_files, sha256_file), basename(archived_cache_files)
      )
    ),
    source_cache_restored = !any(grepl("[.]rds$", cache_source_files, ignore.case = TRUE)),
    mapping_script = mapping_script,
    mapping_script_sha256 = sha256_file(mapping_script),
    runtime = list(
      R = R.version.string,
      package_versions = setNames(
        lapply(c("data.table", "jsonlite", "optparse", "Rcpp", "RcppParallel", "STRINGdb"),
               package_version_or_na),
        c("data.table", "jsonlite", "optparse", "Rcpp", "RcppParallel", "STRINGdb")
      )
    ),
    log = rebuild_log
  )
}

sample_tree <- c(P1 = "outP1", P2 = "outP2", P5 = "final_out", LN = "out")

option_list <- list(
  make_option("--benchmark-root", dest = "benchmark_root", type = "character",
              help = "Root containing outP1, outP2, final_out, and out"),
  make_option("--reanalysis-root", dest = "reanalysis_root", type = "character",
              help = "Root containing <sample>/<sample>_module_membership_shuffle.tsv"),
  make_option("--output-root", dest = "output_root", type = "character",
              help = "New output directory"),
  make_option("--samples", type = "character", default = "P1,P2,P5,LN",
              help = "Comma-separated subset (default: P1,P2,P5,LN)"),
  make_option("--resume", action = "store_true", default = FALSE,
              help = "Reuse completed sample outputs in a non-empty output directory"),
  make_option("--n-random", dest = "n_random", type = "integer", default = 200000L),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--min-module-genes", dest = "min_module_genes", type = "integer", default = 3L)
)

opt <- parse_args(OptionParser(option_list = option_list))
for (name in c("benchmark_root", "reanalysis_root", "output_root")) {
  value <- opt[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", gsub("_", "-", name))
}
if (opt$n_random != 200000L) stop("Correction freeze requires --n-random 200000")
if (opt$seed != 1L) stop("Correction freeze requires --seed 1")
if (opt$min_module_genes != 3L) stop("Correction freeze requires --min-module-genes 3")

samples <- trimws(strsplit(opt$samples, ",", fixed = TRUE)[[1]])
samples <- samples[nzchar(samples)]
unknown <- setdiff(samples, names(sample_tree))
if (length(unknown)) stop("Unknown samples: ", paste(unknown, collapse = ", "))
if (!length(samples)) stop("No samples selected")

benchmark_root <- normalizePath(opt$benchmark_root)
reanalysis_root <- normalizePath(opt$reanalysis_root)
output_root <- normalizePath(opt$output_root, mustWork = FALSE)
if (!isTRUE(opt$resume) && dir.exists(output_root) &&
    length(list.files(output_root, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty: ", output_root)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

repo_root <- normalizePath(file.path(script_dir(), ".."))
module_script <- file.path(repo_root, "benchmark-Rscripts", "module-level.R")
module_cpp <- file.path(repo_root, "benchmark-Rscripts", "module-level.cpp")
mapping_script <- file.path(repo_root, "benchmark-Rscripts", "mapping.R")
if (!file.exists(module_script) || !file.exists(module_cpp) || !file.exists(mapping_script)) {
  stop("Published benchmark implementation is missing")
}
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript is unavailable")

all_modules <- list()
sample_manifests <- list()

for (sample in samples) {
  completed_path <- file.path(
    output_root, sample, "module_enrichment_anchor_cross_null_summary.tsv"
  )
  completed_manifest <- file.path(output_root, sample, "manifest.json")
  if (isTRUE(opt$resume) && file.exists(completed_path) && file.exists(completed_manifest)) {
    message("[", sample, "] reusing completed output")
    all_modules[[sample]] <- fread(completed_path)
    sample_manifests[[sample]] <- completed_manifest
    next
  }
  message("[", sample, "] validating frozen benchmark inputs")
  sample_out <- file.path(output_root, sample)
  staged_dir <- file.path(sample_out, "staged_inputs")
  run_out <- file.path(sample_out, "genescope_recomputed")
  tree_root <- file.path(benchmark_root, sample_tree[[sample]])
  compare_root <- file.path(tree_root, "for_compare")
  map_dir <- file.path(tree_root, "string1step")
  old_anchor_dir <- file.path(tree_root, "anchor-null")
  old_anchor_summary <- file.path(old_anchor_dir, "module_enrichment_anchor_cross_null_summary.tsv")
  old_anchor_meta <- file.path(old_anchor_dir, "meta.json")
  new_membership <- file.path(
    reanalysis_root, sample, paste0(sample, "_module_membership_shuffle.tsv")
  )
  map_path <- file.path(map_dir, "stringdb_edge_mapped_all.tsv")

  required <- c(old_anchor_summary, old_anchor_meta, new_membership, map_path)
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Missing inputs for ", sample, ": ", paste(missing, collapse = ", "))

  source_modules <- setNames(
    file.path(compare_root, c("genescope", "giotto", "hotspot", "seagal"),
              "repeat_001", "modules.tsv"),
    c("genescope", "giotto", "hotspot", "seagal")
  )
  missing <- source_modules[!file.exists(source_modules)]
  if (length(missing)) stop("Missing frozen modules for ", sample, ": ", paste(missing, collapse = ", "))

  old_meta <- fromJSON(old_anchor_meta, simplifyVector = FALSE)
  if (!grepl("anchor_cross_edges", old_meta$null_definition, fixed = TRUE)) {
    stop("Frozen comparator summary does not declare anchor-cross null: ", old_anchor_meta)
  }
  assert_identical(old_meta$params$n_random, 200000L, paste(sample, "frozen n_random"))
  assert_identical(old_meta$params$seed, 1L, paste(sample, "frozen seed"))
  assert_identical(old_meta$params$min_module_genes, 3L, paste(sample, "frozen min_module_genes"))
  assert_identical(old_meta$params$string_score_threshold, 700L,
                   paste(sample, "frozen STRING threshold"))
  observed_map_sha256 <- sha256_file(map_path)
  expected_map_sha256 <- old_meta$map_sha256
  map_rebuild <- NULL
  if (sample == "LN" && !identical(observed_map_sha256, expected_map_sha256)) {
    map_rebuild <- rebuild_ln_map_if_needed(
      tree_root = tree_root,
      compare_root = compare_root,
      sample_out = sample_out,
      old_meta = old_meta,
      mapping_script = mapping_script,
      rscript = rscript,
      current_map_path = map_path
    )
    map_dir <- map_rebuild$map_dir
    map_path <- file.path(map_dir, "stringdb_edge_mapped_all.tsv")
    observed_map_sha256 <- sha256_file(map_path)
  }
  map_validation <- "exact_sha256"
  if (!identical(observed_map_sha256, expected_map_sha256)) {
    assert_map_meta_semantics(
      file.path(map_dir, "stringdb_edge_mapped_meta.json"), old_meta$map_meta, sample
    )
    map_validation <- "semantic_metadata_plus_observed_module_gate"
  }
  for (method in names(source_modules)) {
    assert_identical(
      sha256_file(source_modules[[method]]),
      old_meta$modules_sha256_by_method[[method]],
      paste(sample, method, "modules hash")
    )
  }

  dir.create(staged_dir, recursive = TRUE, showWarnings = FALSE)

  if (map_validation != "exact_sha256") {
    message("[", sample, "] full map hash drifted; running observed-statistics semantic gate")
    semantic_out <- file.path(sample_out, "map_semantic_gate")
    semantic_log <- file.path(sample_out, "map_semantic_gate.log")
    semantic_args <- c(
      module_script,
      "--map_dir", map_dir,
      "--outdir", semantic_out,
      "--methods", "genescope,giotto,hotspot,seagal",
      "--modules_tsv_by_method", paste(unname(source_modules), collapse = ","),
      "--string_score_threshold", "700",
      "--n_random", "0",
      "--seed", "1",
      "--min_module_genes", "3",
      "--min_valid_null_draws", "0",
      "--plot_max_null_points_per_method", "0",
      "--plot_max_modules_per_method", "200"
    )
    semantic_status <- system2(rscript, semantic_args, stdout = semantic_log, stderr = semantic_log)
    if (!identical(semantic_status, 0L)) {
      stop("STRING map semantic gate failed to run for ", sample, "; see ", semantic_log)
    }
    semantic_summary <- fread(file.path(
      semantic_out, "module_enrichment_anchor_cross_null_summary.tsv"
    ))
    frozen_summary_for_gate <- fread(old_anchor_summary)
    gate_columns <- c(
      "module_size", "background_size", "n_pairs_total", "E_sum_obs", "E_avg_obs",
      "n_pairs_mapped_obs", "mapped_pair_frac_obs", "n_pairs_pos_obs", "pos_pair_frac_obs"
    )
    gate <- merge(
      semantic_summary[, c("method", "module_id", gate_columns), with = FALSE],
      frozen_summary_for_gate[, c("method", "module_id", gate_columns), with = FALSE],
      by = c("method", "module_id"), suffixes = c("_current", "_frozen"), all = TRUE
    )
    if (nrow(gate) != nrow(frozen_summary_for_gate)) {
      stop(sample, " semantic map gate row count mismatch")
    }
    for (column in gate_columns) {
      current <- suppressWarnings(as.numeric(gate[[paste0(column, "_current")]]))
      frozen <- suppressWarnings(as.numeric(gate[[paste0(column, "_frozen")]]))
      if (any(!is.finite(current)) || any(!is.finite(frozen)) ||
          any(abs(current - frozen) > 1e-10)) {
        differences <- gate[
          abs(get(paste0(column, "_current")) - get(paste0(column, "_frozen"))) > 1e-10,
          c("method", "module_id", paste0(column, c("_current", "_frozen"))), with = FALSE
        ]
        fwrite(differences, file.path(sample_out, "map_semantic_differences.tsv"), sep = "\t")
        stop(sample, " semantic map gate mismatch in ", column)
      }
    }

    message("[", sample, "] observed semantics match; running full 200,000-draw equivalence gate")
    full_gate_out <- file.path(sample_out, "map_semantic_full_200k_gate")
    full_gate_log <- file.path(sample_out, "map_semantic_full_200k_gate.log")
    full_gate_args <- c(
      module_script,
      "--map_dir", map_dir,
      "--outdir", full_gate_out,
      "--methods", "genescope,giotto,hotspot,seagal",
      "--modules_tsv_by_method", paste(unname(source_modules), collapse = ","),
      "--string_score_threshold", "700",
      "--n_random", "200000",
      "--seed", "1",
      "--min_module_genes", "3",
      "--min_valid_null_draws", "80000",
      "--max_resample_attempts", "20",
      "--plot_max_null_points_per_method", "20000",
      "--plot_max_modules_per_method", "200"
    )
    full_gate_status <- system2(
      rscript, full_gate_args, stdout = full_gate_log, stderr = full_gate_log
    )
    if (!identical(full_gate_status, 0L)) {
      stop("Full 200,000-draw STRING map gate failed for ", sample, "; see ", full_gate_log)
    }
    full_gate <- fread(file.path(
      full_gate_out, "module_enrichment_anchor_cross_null_summary.tsv"
    ))
    null_columns <- c(
      "null_mean_Eavg", "null_sd_Eavg", "delta_Eavg", "empirical_p_Eavg",
      "null_mean_pos_pair_frac", "null_sd_pos_pair_frac", "delta_pos_pair_frac",
      "empirical_p_pos_pair_frac", "n_valid_null_draws"
    )
    full_compare <- merge(
      full_gate[, c("method", "module_id", null_columns), with = FALSE],
      frozen_summary_for_gate[, c("method", "module_id", null_columns), with = FALSE],
      by = c("method", "module_id"), suffixes = c("_current", "_frozen"), all = TRUE
    )
    if (nrow(full_compare) != nrow(frozen_summary_for_gate)) {
      stop(sample, " full 200,000-draw map gate row count mismatch")
    }
    for (column in null_columns) {
      current <- suppressWarnings(as.numeric(full_compare[[paste0(column, "_current")]]))
      frozen <- suppressWarnings(as.numeric(full_compare[[paste0(column, "_frozen")]]))
      if (any(!is.finite(current)) || any(!is.finite(frozen)) ||
          any(abs(current - frozen) > 1e-12)) {
        stop(sample, " full 200,000-draw map gate mismatch in ", column)
      }
    }
    map_validation <- "semantic_equivalence_observed_and_full_200k"
  }

  staged_genescope <- file.path(staged_dir, "genescope_modules_v102.tsv")
  staged <- stage_genescope_modules(source_modules[["genescope"]], new_membership, staged_genescope)

  background_paths <- setNames(character(3), c("giotto", "hotspot", "seagal"))
  for (method in names(background_paths)) {
    background_paths[[method]] <- file.path(staged_dir, paste0("background_", method, ".tsv"))
    stage_background_only(source_modules[[method]], background_paths[[method]])
  }

  # The original geneSCOPE gene list is retained as an unassigned benchmark
  # universe, while its module IDs are replaced by the historical v1.0.2
  # candidate memberships that underlie the audited result series.
  # This keeps comparator inputs and their published Monte Carlo baselines frozen.
  background_union <- unique(c(
    staged$gene,
    unlist(lapply(background_paths, function(path) fread(path, select = "gene")$gene), use.names = FALSE)
  ))
  assert_identical(length(background_union), old_meta$background_size,
                   paste(sample, "frozen benchmark background size"))

  message("[", sample, "] recomputing geneSCOPE (200,000 anchor-cross draws per module)")
  run_args <- c(
    module_script,
    "--map_dir", map_dir,
    "--outdir", run_out,
    "--methods", "genescope,bg_giotto,bg_hotspot,bg_seagal",
    "--modules_tsv_by_method", paste(
      c(staged_genescope, unname(background_paths)), collapse = ","
    ),
    "--string_score_threshold", "700",
    "--n_random", "200000",
    "--seed", "1",
    "--min_module_genes", "3",
    "--min_valid_null_draws", "80000",
    "--max_resample_attempts", "20",
    "--plot_max_null_points_per_method", "20000",
    "--plot_max_modules_per_method", "200"
  )
  log_path <- file.path(sample_out, "module-level.log")
  status <- system2(rscript, run_args, stdout = log_path, stderr = log_path)
  if (!identical(status, 0L)) stop("module-level.R failed for ", sample, "; see ", log_path)

  recomputed_path <- file.path(run_out, "module_enrichment_anchor_cross_null_summary.tsv")
  if (!file.exists(recomputed_path)) stop("Missing recomputed summary for ", sample)
  recomputed <- fread(recomputed_path)
  recomputed <- recomputed[method == "genescope"]
  if (!nrow(recomputed)) stop("No recomputed geneSCOPE modules for ", sample)
  if (recomputed[, any(status != "ok")]) stop("Non-ok geneSCOPE module in ", recomputed_path)
  if (recomputed[, any(n_valid_null_draws != 200000L)]) {
    stop("Not every geneSCOPE module has 200,000 valid null draws for ", sample)
  }

  frozen <- fread(old_anchor_summary)
  comparators <- frozen[method %in% c("giotto", "hotspot", "seagal")]
  if (!nrow(comparators) || comparators[, any(status != "ok")]) {
    stop("Frozen comparator rows are missing or invalid for ", sample)
  }
  if (comparators[, any(n_valid_null_draws != 200000L)]) {
    stop("Frozen comparator rows do not all have 200,000 valid draws for ", sample)
  }

  combined <- rbindlist(list(recomputed, comparators), use.names = TRUE, fill = TRUE)
  combined[, sample := sample]
  setcolorder(combined, c("sample", setdiff(names(combined), "sample")))
  combined_path <- file.path(sample_out, "module_enrichment_anchor_cross_null_summary.tsv")
  fwrite(combined, combined_path, sep = "\t")
  all_modules[[sample]] <- combined

  sample_manifest <- list(
    sample = sample,
    definition = paste(
      "For each module, sample k=choose(m,2) pairs without replacement from",
      "the module-by-outside-background cross-pair space; STRING-absent pairs are zero."
    ),
    parameters = list(n_random = 200000L, seed = 1L, min_module_genes = 3L,
                      string_score_threshold = 700L),
    correction_scope = paste(
      "geneSCOPE modules recomputed from the historical v1.0.2 candidate membership;",
      "comparator module rows and",
      "Monte Carlo summaries reused only after exact input-hash and parameter validation."
    ),
    background_strategy = paste(
      "Frozen published benchmark gene universe retained; old geneSCOPE module IDs reset",
      "to -1 and replaced with the historical v1.0.2 candidate assignments."
    ),
    background_size = length(background_union),
    files = list(
      new_membership = list(path = new_membership, sha256 = sha256_file(new_membership)),
      staged_genescope = list(path = staged_genescope, sha256 = sha256_file(staged_genescope)),
      string_map = list(
        path = map_path,
        sha256 = observed_map_sha256,
        frozen_expected_sha256 = expected_map_sha256,
        validation = map_validation
      ),
      frozen_anchor_summary = list(path = old_anchor_summary,
                                   sha256 = sha256_file(old_anchor_summary)),
      combined_summary = list(path = combined_path, sha256 = sha256_file(combined_path))
    ),
    comparator_module_sha256 = lapply(source_modules[c("giotto", "hotspot", "seagal")], sha256_file),
    implementation = list(
      module_level_R = list(path = module_script, sha256 = sha256_file(module_script)),
      module_level_cpp = list(path = module_cpp, sha256 = sha256_file(module_cpp))
    ),
    map_rebuild = map_rebuild
  )
  manifest_path <- file.path(sample_out, "manifest.json")
  write_json(sample_manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = NA)
  sample_manifests[[sample]] <- manifest_path
}

all_dt <- rbindlist(all_modules, use.names = TRUE, fill = TRUE)
all_path <- file.path(output_root, "ALL_module_enrichment_anchor_cross_null_summary.tsv")
fwrite(all_dt, all_path, sep = "\t")

ok <- all_dt[
  status == "ok" & is.finite(E_avg_obs) & is.finite(delta_Eavg)
]
summary_dt <- ok[, .(
  n_modules = .N,
  median_E_avg_obs = median(E_avg_obs),
  mean_E_avg_obs = mean(E_avg_obs),
  q25_E_avg_obs = as.numeric(quantile(E_avg_obs, 0.25, names = FALSE)),
  q75_E_avg_obs = as.numeric(quantile(E_avg_obs, 0.75, names = FALSE)),
  median_delta_Eavg = median(delta_Eavg),
  mean_delta_Eavg = mean(delta_Eavg),
  q25_delta_Eavg = as.numeric(quantile(delta_Eavg, 0.25, names = FALSE)),
  q75_delta_Eavg = as.numeric(quantile(delta_Eavg, 0.75, names = FALSE)),
  positive_delta_modules = sum(delta_Eavg > 0),
  positive_delta_fraction = mean(delta_Eavg > 0)
), by = .(sample, method)]

summary_dt[, rank_median_E_avg_obs := frank(-median_E_avg_obs, ties.method = "min"), by = sample]
summary_dt[, rank_median_delta_Eavg := frank(-median_delta_Eavg, ties.method = "min"), by = sample]
summary_dt[, rank_mean_E_avg_obs := frank(-mean_E_avg_obs, ties.method = "min"), by = sample]
summary_dt[, rank_mean_delta_Eavg := frank(-mean_delta_Eavg, ties.method = "min"), by = sample]
setorder(summary_dt, sample, rank_median_delta_Eavg, rank_median_E_avg_obs, method)

summary_path <- file.path(output_root, "summary_by_sample_method.tsv")
fwrite(summary_dt, summary_path, sep = "\t")

gene_scope <- summary_dt[method == "genescope"]
if (nrow(gene_scope) != length(samples)) stop("Missing geneSCOPE method summary")
strength <- gene_scope[, .(
  sample,
  top_median_E_avg_obs = rank_median_E_avg_obs == 1L,
  top_median_delta_Eavg = rank_median_delta_Eavg == 1L,
  top_mean_E_avg_obs = rank_mean_E_avg_obs == 1L,
  top_mean_delta_Eavg = rank_mean_delta_Eavg == 1L,
  module_level_strong = rank_median_E_avg_obs == 1L & rank_median_delta_Eavg == 1L
)]
strength_path <- file.path(output_root, "genescope_module_level_strength.tsv")
fwrite(strength, strength_path, sep = "\t")

conclusion <- if (strength[, all(module_level_strong)]) {
  "geneSCOPE remains strongest at module level: it ranks first by the median observed and null-adjusted module STRING scores in every sample."
} else {
  paste0(
    "geneSCOPE remains competitive at module level, but it is not first by both median module metrics in every sample. ",
    "See genescope_module_level_strength.tsv."
  )
}
writeLines(conclusion, file.path(output_root, "CONCLUSION.txt"), useBytes = TRUE)

root_manifest <- list(
  definition_source = list(
    article = "bbag302.pdf, Methods p.6 and Figure 6 caption p.13",
    null = paste(
      "k=choose(m,2) pairs sampled without replacement from the module-background",
      "cross-pair space; one endpoint in the focal module and one outside; missing STRING pairs score 0."
    )
  ),
  parameters = list(n_random = 200000L, seed = 1L, min_module_genes = 3L,
                    string_score_threshold = 700L),
  samples = samples,
  conclusion = conclusion,
  outputs = list(
    all_modules = list(path = all_path, sha256 = sha256_file(all_path)),
    summary = list(path = summary_path, sha256 = sha256_file(summary_path)),
    strength = list(path = strength_path, sha256 = sha256_file(strength_path))
  ),
  sample_manifests = sample_manifests
)
write_json(root_manifest, file.path(output_root, "manifest.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

message("Wrote: ", output_root)
message(conclusion)

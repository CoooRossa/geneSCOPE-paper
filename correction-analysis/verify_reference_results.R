#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
result_root <- if (length(args)) args[[1L]] else Sys.getenv("RESULT_ROOT", unset = "")
if (!nzchar(result_root)) stop("Usage: verify_reference_results.R RESULT_ROOT")
result_root <- normalizePath(result_root, mustWork = TRUE)

script_dir <- local({
  x <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", x[[1L]])))
})
source(file.path(script_dir, "..", "main-text-scripts", "freeze_helpers.R"))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.")

sample_ref <- utils::read.delim(file.path(script_dir, "samples.tsv"), stringsAsFactors = FALSE)
top6_ref <- utils::read.delim(file.path(script_dir, "reference_top6.tsv"), stringsAsFactors = FALSE)
wilcox_ref <- utils::read.delim(file.path(script_dir, "reference_wilcoxon.tsv"), stringsAsFactors = FALSE)
scale_ref <- utils::read.delim(file.path(script_dir, "reference_p5_multiscale.tsv"), stringsAsFactors = FALSE)
dendro_ref <- utils::read.delim(file.path(script_dir, "reference_p5_dendro_path.tsv"), stringsAsFactors = FALSE)
membership_ref <- utils::read.delim(
  file.path(script_dir, "reference_membership_digests.tsv"), stringsAsFactors = FALSE
)
artifact_hashes <- utils::read.delim(file.path(script_dir, "reference_artifact_sha256.tsv"),
                                     stringsAsFactors = FALSE)
external_hashes <- utils::read.delim(file.path(script_dir, "reference_external_bundle_hashes.tsv"),
                                     stringsAsFactors = FALSE)
artifact_hashes <- rbind(artifact_hashes, external_hashes)

failures <- character()
checks <- character()
record <- function(ok, label, detail = "") {
  line <- paste0(if (isTRUE(ok)) "PASS  " else "FAIL  ", label,
                 if (nzchar(detail)) paste0(" - ", detail) else "")
  message(line)
  checks <<- c(checks, line)
  if (!isTRUE(ok)) failures <<- c(failures, line)
  invisible(ok)
}

find_first <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1L]] else ""
}

sha256_program <- if (nzchar(Sys.which("sha256sum"))) {
  "sha256sum"
} else if (nzchar(Sys.which("shasum"))) {
  "shasum"
} else {
  ""
}
sha256_executable <- if (nzchar(sha256_program)) Sys.which(sha256_program) else ""

sha256_file <- function(path) {
  if (!nzchar(sha256_executable) || !file.exists(path)) return("")
  command_args <- if (identical(sha256_program, "shasum")) {
    c("-a", "256", shQuote(path))
  } else {
    shQuote(path)
  }
  out <- system2(sha256_executable, command_args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status", exact = TRUE)
  digest <- if (length(out)) strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]] else ""
  if ((!is.null(status) && status != 0L) || !grepl("^[0-9a-fA-F]{64}$", digest)) return("")
  tolower(digest)
}

resolve_result_dir <- function(dir_name, markers) {
  direct_match <- all(file.exists(file.path(result_root, markers)))
  candidates <- unique(c(
    if (direct_match) result_root else character(),
    file.path(result_root, dir_name),
    file.path(dirname(result_root), dir_name)
  ))
  hit <- candidates[dir.exists(candidates)]
  if (length(hit)) normalizePath(hit[[1L]], mustWork = TRUE) else ""
}

main_root <- resolve_result_dir("full_reanalysis_shuffle_v102_20260726",
                                c("P1", "P2", "P5", "LN"))
extra_root <- resolve_result_dir("shuffle_extra_figs_20260726", c("S15", "S6"))
benchmark_root <- resolve_result_dir("shuffle_jaccard_benchmark_20260726",
                                    c("ALL_four_method_jaccard.tsv",
                                      "edge_benchmark_manifest.json"))

ln_complete_dir <- Sys.getenv("LN_COMPLETE_DELTA_DIR", unset = "")
if (!nzchar(ln_complete_dir)) {
  ln_complete_dir <- find_first(c(
    file.path(result_root, "LN_delta_complete_v102"),
    file.path(result_root, "audit_results", "LN_delta_complete_v102"),
    file.path(dirname(result_root), "audit_results", "LN_delta_complete_v102"),
    file.path(script_dir, "..", "..", "..", "audit_results", "LN_delta_complete_v102")
  ))
}
if (nzchar(ln_complete_dir) && dir.exists(ln_complete_dir)) {
  ln_complete_dir <- normalizePath(ln_complete_dir, mustWork = TRUE)
}

anchor_cross_dir <- Sys.getenv("ANCHOR_CROSS_RESULTS", unset = "")
if (!nzchar(anchor_cross_dir)) {
  anchor_cross_dir <- find_first(c(
    file.path(result_root, "anchor_cross_benchmark_v102"),
    file.path(result_root, "audit_results", "anchor_cross_benchmark_v102"),
    file.path(dirname(result_root), "audit_results", "anchor_cross_benchmark_v102"),
    file.path(script_dir, "..", "..", "..", "audit_results",
              "anchor_cross_benchmark_v102")
  ))
}
if (nzchar(anchor_cross_dir) && dir.exists(anchor_cross_dir)) {
  anchor_cross_dir <- normalizePath(anchor_cross_dir, mustWork = TRUE)
}

expected_hash_columns <- c("root_token", "relative_path", "sha256")
if (!identical(names(artifact_hashes), expected_hash_columns)) {
  stop("Artifact SHA256 references must have exactly these columns: ",
       paste(expected_hash_columns, collapse = ", "))
}
known_root_tokens <- c("MAIN_RESULTS", "EXTRA_RESULTS", "BENCHMARK_RESULTS",
                       "LN_COMPLETE_RESULTS", "ANCHOR_CROSS_RESULTS")
if (any(!artifact_hashes$root_token %in% known_root_tokens)) {
  stop("Unknown artifact root token(s): ",
       paste(unique(artifact_hashes$root_token[
         !artifact_hashes$root_token %in% known_root_tokens]), collapse = ", "))
}
unsafe_relative_path <- vapply(strsplit(artifact_hashes$relative_path, "/", fixed = TRUE),
                               function(parts) any(parts %in% c("", ".", "..")), logical(1L)) |
  grepl("^/", artifact_hashes$relative_path)
if (any(unsafe_relative_path)) {
  stop("Artifact hash list contains an unsafe relative path: ",
       artifact_hashes$relative_path[[which(unsafe_relative_path)[[1L]]]])
}
if (any(!grepl("^[[:xdigit:]]{64}$", artifact_hashes$sha256))) {
  stop("Artifact hash list contains a malformed SHA256 value.")
}
artifact_keys <- paste(artifact_hashes$root_token, artifact_hashes$relative_path, sep = ":")
if (anyDuplicated(artifact_keys)) stop("Artifact hash list contains duplicate paths.")

artifact_roots <- c(
  MAIN_RESULTS = main_root,
  EXTRA_RESULTS = extra_root,
  BENCHMARK_RESULTS = benchmark_root,
  LN_COMPLETE_RESULTS = ln_complete_dir,
  ANCHOR_CROSS_RESULTS = anchor_cross_dir
)
record(nzchar(sha256_executable), "SHA256 verifier available", sha256_executable)
for (root_token in known_root_tokens) {
  token_root <- artifact_roots[[root_token]]
  record(nzchar(token_root) && dir.exists(token_root),
         paste0("artifact root: ", root_token), token_root)
  token_rows <- artifact_hashes[artifact_hashes$root_token == root_token, , drop = FALSE]
  for (i in seq_len(nrow(token_rows))) {
    artifact_path <- if (nzchar(token_root)) {
      file.path(token_root, token_rows$relative_path[[i]])
    } else {
      ""
    }
    observed_sha <- sha256_file(artifact_path)
    expected_sha <- tolower(token_rows$sha256[[i]])
    record(identical(tolower(observed_sha), expected_sha),
           paste0("artifact SHA256: ", root_token, "/", token_rows$relative_path[[i]]),
           if (nzchar(observed_sha)) observed_sha else "missing or unreadable")
  }
}

if (nzchar(anchor_cross_dir) && dir.exists(anchor_cross_dir)) {
  anchor_summary_path <- file.path(anchor_cross_dir, "summary_by_sample_method.tsv")
  anchor_all_path <- file.path(
    anchor_cross_dir, "ALL_module_enrichment_anchor_cross_null_summary.tsv"
  )
  anchor_manifest_path <- file.path(anchor_cross_dir, "manifest.json")
  anchor_summary <- utils::read.delim(anchor_summary_path, stringsAsFactors = FALSE)
  gene_scope_summary <- anchor_summary[anchor_summary$method == "genescope", , drop = FALSE]
  expected_samples <- c("P1", "P2", "P5", "LN")
  rank_columns <- c(
    "rank_median_E_avg_obs", "rank_median_delta_Eavg",
    "rank_mean_E_avg_obs", "rank_mean_delta_Eavg"
  )
  rank_ok <- nrow(gene_scope_summary) == length(expected_samples) &&
    setequal(gene_scope_summary$sample, expected_samples) &&
    all(rank_columns %in% names(gene_scope_summary)) &&
    all(as.matrix(gene_scope_summary[, rank_columns, drop = FALSE]) == 1)
  record(rank_ok, "published anchor-cross module benchmark: geneSCOPE ranks first",
         "P1/P2/P5/LN; observed and enrichment medians/means")

  anchor_all <- utils::read.delim(anchor_all_path, stringsAsFactors = FALSE)
  gene_scope_modules <- anchor_all[anchor_all$method == "genescope", , drop = FALSE]
  null_ok <- nrow(gene_scope_modules) > 0L &&
    setequal(gene_scope_modules$sample, expected_samples) &&
    all(gene_scope_modules$status == "ok") &&
    all(as.integer(gene_scope_modules$n_valid_null_draws) == 200000L)
  record(null_ok, "published anchor-cross module benchmark: complete null draws",
         paste0(nrow(gene_scope_modules), " geneSCOPE modules; 200000 valid draws each"))

  anchor_manifest <- jsonlite::read_json(anchor_manifest_path, simplifyVector = TRUE)
  parameter_ok <- identical(as.integer(anchor_manifest$parameters$n_random), 200000L) &&
    identical(as.integer(anchor_manifest$parameters$seed), 1L) &&
    identical(as.integer(anchor_manifest$parameters$min_module_genes), 3L)
  record(parameter_ok, "published anchor-cross module benchmark: frozen parameters",
         "n_random=200000; seed=1; minimum module size=3")

  definition_ok <- TRUE
  for (sample_id in expected_samples) {
    sample_manifest_path <- file.path(anchor_cross_dir, sample_id, "manifest.json")
    sample_manifest <- jsonlite::read_json(sample_manifest_path, simplifyVector = TRUE)
    definition_ok <- definition_ok &&
      grepl("module-by-outside-background cross-pair space", sample_manifest$definition,
            fixed = TRUE) &&
      grepl("STRING-absent pairs are zero", sample_manifest$definition, fixed = TRUE)
  }
  record(definition_ok, "published anchor-cross module benchmark: null definition",
         "module-to-background cross-pairs; missing STRING pairs scored zero")

  ln_anchor_manifest <- jsonlite::read_json(
    file.path(anchor_cross_dir, "LN", "manifest.json"), simplifyVector = TRUE
  )
  ln_map <- ln_anchor_manifest$map_rebuild
  ln_map_ok <- is.list(ln_map) &&
    identical(ln_anchor_manifest$files$string_map$validation, "exact_sha256") &&
    identical(ln_map$rebuilt_map_sha256, ln_map$expected_map_sha256)
  record(ln_map_ok, "LN published STRING mapping reconstructed exactly",
         if (ln_map_ok) ln_map$rebuilt_map_sha256 else "missing exact mapping provenance")

  source_hashes <- unlist(ln_map$cache_source_file_sha256, use.names = TRUE)
  source_paths <- file.path(ln_map$cache_source, names(source_hashes))
  source_cache_ok <- identical(ln_map$source_cache_restored, TRUE) &&
    length(source_hashes) > 0L &&
    !any(grepl("[.]rds$", names(source_hashes), ignore.case = TRUE)) &&
    all(vapply(seq_along(source_paths), function(i) {
      identical(tolower(sha256_file(source_paths[[i]])), tolower(source_hashes[[i]]))
    }, logical(1L)))
  record(source_cache_ok, "LN frozen STRING source cache unchanged",
         paste0(length(source_hashes), " raw cache files; no generated RDS"))

  generated_cache <- ln_map$generated_rds_cache
  generated_hashes <- unlist(generated_cache$file_sha256, use.names = TRUE)
  generated_paths <- file.path(generated_cache$path, names(generated_hashes))
  generated_cache_ok <- length(generated_hashes) > 0L &&
    all(grepl("[.]rds$", names(generated_hashes), ignore.case = TRUE)) &&
    all(vapply(seq_along(generated_paths), function(i) {
      identical(tolower(sha256_file(generated_paths[[i]])), tolower(generated_hashes[[i]]))
    }, logical(1L)))
  record(generated_cache_ok, "LN generated STRING cache archived separately",
         paste0(length(generated_hashes), " generated RDS files"))
}

ln_complete_manifest <- NULL
if (nzchar(ln_complete_dir) && dir.exists(ln_complete_dir)) {
  ln_manifest_path <- file.path(ln_complete_dir, "manifest.json")
  if (file.exists(ln_manifest_path)) {
    ln_complete_manifest <- jsonlite::read_json(ln_manifest_path, simplifyVector = TRUE)
  }
}

for (sample_id in sample_ref$sample_id) {
  cfg <- sample_ref[sample_ref$sample_id == sample_id, , drop = FALSE]
  sample_dir <- file.path(main_root, sample_id)
  manifest_path <- find_first(c(file.path(sample_dir, "manifest.json"),
                                file.path(sample_dir, paste0(sample_id, "_manifest.json"))))
  summary_path <- find_first(c(file.path(sample_dir, paste0(sample_id, "_summary.tsv"))))
  manifest <- if (nzchar(manifest_path)) {
    jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  } else NULL
  delta_complete <- FALSE
  delta_detail <- "missing delta_permutation provenance"

  record(!is.null(manifest), paste0(sample_id, " manifest present"), manifest_path)
  if (!is.null(manifest)) {
    pkg <- if (!is.null(manifest$package_version)) as.character(manifest$package_version) else ""
    record(identical(pkg, "1.0.2"), paste0(sample_id, " package version"), pkg)

    b <- if (!is.null(manifest$permutations)) manifest$permutations else manifest$permutation$permutations
    seed <- if (!is.null(manifest$seed)) manifest$seed else manifest$permutation$seed
    scheme <- if (is.character(manifest$permutation)) manifest$permutation else manifest$permutation$scheme
    record(identical(as.integer(b), 1000L), paste0(sample_id, " B=1000"), as.character(b))
    record(identical(as.integer(seed), 1L), paste0(sample_id, " seed=1"), as.character(seed))
    record(grepl("global_joint_shuffle|all-grid joint shuffle", scheme),
           paste0(sample_id, " global joint shuffle"), scheme)

    gate <- manifest$gate_max_abs_L_diff
    gate_ok <- is.numeric(gate) && length(gate) == 1L && is.finite(gate) && gate <= 1e-12
    record(gate_ok, paste0(sample_id, " canonical formula gate"),
           if (gate_ok) format(gate, scientific = TRUE) else "missing or non-numeric")

    dp <- if (identical(sample_id, "LN") && !is.null(ln_complete_manifest)) {
      ln_complete_manifest$delta_permutation
    } else {
      manifest$delta_permutation
    }
    delta_complete <- !is.null(dp) && identical(dp$p_adj_mode, "BH_universe") &&
      identical(dp$use_blocks, FALSE) &&
      is.finite(dp$total_universe) && is.finite(dp$selected_pairs) &&
      dp$total_universe == dp$selected_pairs
    delta_detail <- if (is.null(dp)) {
      "missing delta_permutation provenance"
    } else {
      paste0("eligible=", dp$total_universe, ", selected=", dp$selected_pairs,
             ", adjustment=", dp$p_adj_mode)
    }
  }

  if (nzchar(summary_path)) {
    sm <- utils::read.delim(summary_path, stringsAsFactors = FALSE)
    mods <- sm$shuffle_modules[[1L]]
    assigned <- sm$shuffle_assigned[[1L]]
  } else if (!is.null(manifest$clustering)) {
    mods <- manifest$clustering$modules
    assigned <- manifest$clustering$assigned
  } else {
    mods <- assigned <- NA_integer_
  }
  record(identical(as.integer(mods), as.integer(cfg$expected_modules)) &&
           identical(as.integer(assigned), as.integer(cfg$expected_assigned)),
         paste0(sample_id, " module anchor"),
         paste0(mods, " modules / ", assigned, " assigned"))

  membership_path <- find_first(c(
    file.path(sample_dir, paste0(sample_id, "_module_membership.tsv")),
    file.path(sample_dir, paste0(sample_id, "_module_membership_shuffle.tsv"))
  ))
  record(nzchar(membership_path), paste0(sample_id, " membership table present"), membership_path)
  if (nzchar(membership_path)) {
    membership_table <- utils::read.delim(
      membership_path, stringsAsFactors = FALSE, check.names = FALSE
    )
    module_col <- intersect(c("raw_module_id", "shuffle_module"), names(membership_table))
    membership_expected <- membership_ref[
      membership_ref$sample_id == sample_id, , drop = FALSE
    ]
    membership_ok <- "gene" %in% names(membership_table) &&
      length(module_col) == 1L && nrow(membership_expected) == 1L
    if (membership_ok) {
      membership_observed <- membership_digests(
        membership_table$gene, membership_table[[module_col[[1L]]]]
      )
      exact_ok <- identical(
        membership_observed$exact_gene_module_sha256,
        membership_expected$exact_gene_module_sha256[[1L]]
      )
      partition_ok <- identical(
        membership_observed$partition_sha256,
        membership_expected$partition_sha256[[1L]]
      )
      membership_ok <- exact_ok && partition_ok
      membership_detail <- paste0(
        "exact=", membership_observed$exact_gene_module_sha256,
        "; partition=", membership_observed$partition_sha256
      )
    } else {
      membership_detail <- "missing required gene/module columns or reference row"
    }
    record(membership_ok, paste0(sample_id, " exact membership/partition gate"),
           membership_detail)
  }

  top_path <- find_first(c(
    if (identical(sample_id, "LN") && nzchar(ln_complete_dir)) {
      file.path(ln_complete_dir, "LN_top_pairs_complete_delta_v102.tsv")
    } else "",
    file.path(sample_dir, paste0(sample_id, "_top_pairs_all.tsv")),
    file.path(sample_dir, paste0(sample_id, "_toplvsr_all_shuffleFDR.tsv"))
  ))
  record(nzchar(top_path), paste0(sample_id, " all-pair table present"), top_path)
  if (nzchar(top_path)) {
    top <- utils::read.delim(top_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (delta_complete && identical(sample_id, "LN")) {
      generator_gate <- try(
        assert_generator_provenance(
          "LN", file.path(ln_complete_dir, "recompute_ln_complete_delta.R"),
          top_pair_rows = nrow(top)
        ),
        silent = TRUE
      )
      delta_complete <- !inherits(generator_gate, "try-error")
      delta_detail <- if (delta_complete) {
        paste0(
          delta_detail, "; generator Top-N=", generator_gate$top_n,
          ", seed=", generator_gate$seed, ", RNG=", generator_gate$rng
        )
      } else {
        paste0("generator provenance failed: ", as.character(generator_gate))
      }
    }
    if (!delta_complete && sample_id %in% c("P1", "P2", "P5")) {
      generator_path <- file.path(main_root, "run_shuffle_reanalysis.R")
      generator_gate <- try(
        assert_generator_provenance(
          sample_id, generator_path, top_pair_rows = nrow(top)
        ),
        silent = TRUE
      )
      delta_complete <- !inherits(generator_gate, "try-error")
      delta_detail <- if (delta_complete) {
        paste0(
          "generator-pinned complete universe: selected=", nrow(top),
          ", Top-N=", generator_gate$top_n,
          ", seed=", generator_gate$seed,
          ", RNG=", generator_gate$rng
        )
      } else {
        paste0("generator provenance failed: ", as.character(generator_gate))
      }
    }
    record(delta_complete, paste0(sample_id, " complete eligible Delta universe"), delta_detail)
    filtered <- filter_display_pairs(top)
    got <- utils::head(filtered[, c("gene1", "gene2"), drop = FALSE], 6L)
    expected <- top6_ref[top6_ref$sample_id == sample_id, c("gene1", "gene2"), drop = FALSE]
    ok <- nrow(got) == 6L && identical(unname(as.matrix(got)), unname(as.matrix(expected)))
    record(ok, paste0(sample_id, " final-filter Top6"),
           "qDelta<.05; L>0; r<.05; pct1,pct2>20; full-precision Delta")
  }
}

wilcox_path <- file.path(extra_root, "S15", "S15_within_between_wilcoxon.tsv")
record(file.exists(wilcox_path), "Wilcoxon table present", wilcox_path)
if (file.exists(wilcox_path)) {
  w <- utils::read.delim(wilcox_path, stringsAsFactors = FALSE)
  for (sample_id in wilcox_ref$sample_id) {
    got <- w$p_value[w$sample == sample_id]
    expected <- wilcox_ref$p_value[wilcox_ref$sample_id == sample_id]
    ok <- length(got) == 1L && isTRUE(all.equal(as.numeric(got), as.numeric(expected),
                                               tolerance = 1e-12))
    record(ok, paste0(sample_id, " Wilcoxon p"), if (length(got)) format(got) else "missing")
  }
}

validate_p5_scale_manifest <- function(m, manifest_path, grid_um) {
  label <- paste0("P5 ", grid_um, " um manifest")
  gate <- m$gate_max_abs_L_diff
  gate_ok <- is.numeric(gate) && length(gate) == 1L && is.finite(gate) && gate <= 1e-12
  record(gate_ok, paste0(label, " canonical formula gate"),
         if (gate_ok) format(gate, scientific = TRUE) else "missing or non-numeric")

  permutation <- m$permutation
  nested_permutation <- is.list(permutation)
  b <- if (nested_permutation) permutation$permutations else m$permutations
  seed <- if (nested_permutation) permutation$seed else m$seed
  scheme <- if (nested_permutation) permutation$scheme else permutation
  use_blocks <- if (nested_permutation) permutation$use_blocks else NULL
  record(is.numeric(b) && length(b) == 1L && identical(as.integer(b), 1000L),
         paste0(label, " B=1000"), if (length(b)) as.character(b) else "missing")
  record(is.numeric(seed) && length(seed) == 1L && identical(as.integer(seed), 1L),
         paste0(label, " seed=1"), if (length(seed)) as.character(seed) else "missing")

  legacy_relative_path <- paste0("S6/grid", grid_um, "/grid", grid_um, "_manifest.json")
  legacy_reference <- artifact_hashes[
    artifact_hashes$root_token == "EXTRA_RESULTS" &
      artifact_hashes$relative_path == legacy_relative_path, , drop = FALSE
  ]
  legacy_global_shuffle <- nrow(legacy_reference) == 1L &&
    identical(tolower(sha256_file(manifest_path)), tolower(legacy_reference$sha256[[1L]])) &&
    is.numeric(m$frac_FDR_lt05_shuffle)
  explicit_global_shuffle <- is.character(scheme) && length(scheme) == 1L &&
    grepl("global_joint_shuffle|all-grid joint shuffle", scheme) &&
    (!nested_permutation || identical(use_blocks, FALSE))
  record(explicit_global_shuffle || legacy_global_shuffle,
         paste0(label, " global joint shuffle"),
         if (explicit_global_shuffle) as.character(scheme) else if (legacy_global_shuffle) {
           "legacy authoritative manifest; audited global-shuffle SHA256 is pinned"
         } else {
           "missing or incompatible permutation provenance"
         })

  new_schema <- !is.null(m$package_source_commit) ||
    !is.null(m$package_vendor_tree_sha256) || is.list(m$inputs) || is.list(m$outputs)
  if (!new_schema) return(invisible(NULL))

  package_version <- if (is.null(m$package_version)) "" else as.character(m$package_version)
  record(identical(package_version, "1.0.2"), paste0(label, " package version"),
         package_version)
  source_commit <- if (is.null(m$package_source_commit)) "" else m$package_source_commit
  vendor_sha <- if (is.null(m$package_vendor_tree_sha256)) "" else m$package_vendor_tree_sha256
  record(is.character(source_commit) && length(source_commit) == 1L &&
           grepl("^[0-9a-f]{40}$", source_commit),
         paste0(label, " package source commit"), source_commit)
  record(is.character(vendor_sha) && length(vendor_sha) == 1L &&
           grepl("^[0-9a-f]{64}$", vendor_sha),
         paste0(label, " package vendor-tree SHA256"), vendor_sha)
  expected_commit <- Sys.getenv("GENESCOPE_SOURCE_COMMIT", unset = "")
  if (nzchar(expected_commit)) {
    record(identical(source_commit, expected_commit),
           paste0(label, " package source commit matches freeze environment"), source_commit)
  }
  expected_vendor_sha <- Sys.getenv("GENESCOPE_VENDOR_TREE_SHA256", unset = "")
  if (nzchar(expected_vendor_sha)) {
    record(identical(vendor_sha, expected_vendor_sha),
           paste0(label, " vendor-tree SHA256 matches freeze environment"), vendor_sha)
  }

  input_hashes <- unlist(m$inputs$sha256, use.names = TRUE)
  input_names <- names(input_hashes)
  record(length(input_hashes) > 0L && !is.null(input_names) &&
           all(grepl("^[0-9a-f]{64}$", input_hashes)),
         paste0(label, " input SHA256 declarations"),
         paste(input_names, collapse = ","))
  input_paths <- c(
    roi = m$inputs$roi_file,
    roi_file = m$inputs$roi_file,
    cell_feature_matrix = file.path(m$inputs$xenium_outs, "cell_feature_matrix.h5"),
    cells_parquet = file.path(m$inputs$xenium_outs, "cells.parquet"),
    transcripts_parquet = file.path(m$inputs$xenium_outs, "transcripts.parquet")
  )
  for (input_name in input_names) {
    input_path <- unname(input_paths[[input_name]])
    observed_sha <- if (length(input_path) && nzchar(input_path)) sha256_file(input_path) else ""
    record(identical(tolower(observed_sha), tolower(input_hashes[[input_name]])),
           paste0(label, " input SHA256: ", input_name),
           if (nzchar(observed_sha)) observed_sha else "missing or unreadable")
  }

  output_hashes <- unlist(m$outputs$sha256, use.names = TRUE)
  output_names <- names(output_hashes)
  required_outputs <- c("membership", "scope")
  record(!is.null(output_names) && all(required_outputs %in% output_names) &&
           all(grepl("^[0-9a-f]{64}$", output_hashes)),
         paste0(label, " output SHA256 declarations"),
         paste(output_names, collapse = ","))
  output_paths <- c(
    membership = file.path(dirname(manifest_path), paste0("grid", grid_um, "_membership.tsv")),
    scope = file.path(dirname(manifest_path), paste0("grid", grid_um, "_scope_v102.rds"))
  )
  for (output_name in intersect(required_outputs, output_names)) {
    observed_sha <- sha256_file(unname(output_paths[[output_name]]))
    record(identical(tolower(observed_sha), tolower(output_hashes[[output_name]])),
           paste0(label, " output SHA256: ", output_name),
           if (nzchar(observed_sha)) observed_sha else "missing or unreadable")
  }
  invisible(NULL)
}

for (i in seq_len(nrow(scale_ref))) {
  grid_um <- scale_ref$grid_um[[i]]
  if (grid_um == 30L) {
    mp <- file.path(main_root, "P5", "manifest.json")
    summary_path <- file.path(main_root, "P5", "P5_summary.tsv")
    if (file.exists(summary_path)) {
      sm <- utils::read.delim(summary_path, stringsAsFactors = FALSE)
      mods <- sm$shuffle_modules[[1L]]
      assigned <- sm$shuffle_assigned[[1L]]
    } else {
      manifest <- jsonlite::read_json(mp, simplifyVector = TRUE)
      mods <- manifest$clustering$modules
      assigned <- manifest$clustering$assigned
    }
  } else {
    mp <- find_first(c(
      file.path(extra_root, "S6", paste0("grid", grid_um), paste0("grid", grid_um, "_manifest.json")),
      file.path(result_root, "P5_multiscale", paste0("grid", grid_um), "manifest.json"),
      file.path(result_root, paste0("grid", grid_um), "manifest.json")
    ))
    if (file.exists(mp)) {
      m <- jsonlite::read_json(mp, simplifyVector = TRUE)
      mods <- if (!is.null(m$modules)) m$modules else m$clustering$modules
      assigned <- if (!is.null(m$assigned)) m$assigned else m$clustering$assigned
    } else {
      mods <- assigned <- NA_integer_
    }
  }
  record(file.exists(mp), paste0("P5 ", grid_um, " um manifest present"), mp)
  if (file.exists(mp)) {
    scale_manifest <- jsonlite::read_json(mp, simplifyVector = TRUE)
    validate_p5_scale_manifest(scale_manifest, mp, grid_um)
  }
  ok <- identical(as.integer(mods), as.integer(scale_ref$modules[[i]])) &&
    identical(as.integer(assigned), as.integer(scale_ref$assigned[[i]]))
  record(ok, paste0("P5 ", grid_um, " um multiscale anchor"),
         paste0(mods, " modules / ", assigned, " assigned"))
}

dendro_path <- Sys.getenv("P5_DENDRO_AUDIT_TSV", unset = "")
if (!nzchar(dendro_path)) {
  dendro_path <- file.path(main_root, "P5", "P5_dendro_path_audit.tsv")
}
record(file.exists(dendro_path), "P5 dendrogram-path audit present", dendro_path)
if (file.exists(dendro_path)) {
  d <- utils::read.delim(dendro_path, stringsAsFactors = FALSE)
  if ("version" %in% names(d)) d <- d[d$version == "shuffle_v1.0.2", , drop = FALSE]
  ok <- nrow(d) == 1L &&
    identical(d$raw_module_path[[1L]], dendro_ref$raw_module_path[[1L]]) &&
    identical(d$display_module_path[[1L]], dendro_ref$display_module_path[[1L]])
  record(ok, "P5 dendrogram module route",
         if (nrow(d)) paste0("raw=", d$raw_module_path[[1L]],
                             "; display=", d$display_module_path[[1L]]) else "missing row")
}

if (length(failures)) {
  stop(length(failures), " freeze gate(s) failed. See messages above.")
}
message("All correction reference-result freeze gates passed.")

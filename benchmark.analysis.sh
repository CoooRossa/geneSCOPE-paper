#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

BENCH_ROOT="${SCRIPT_DIR}"
RSCRIPTS_DIR="$(cd "${SCRIPT_DIR}/../../benchmark-Rscripts" && pwd)"
GENESCOPE_ROOT="${REPO_ROOT}/docker/genescope/geneSCOPE"

METHODS="genescope,giotto,hotspot,seagal"
MAP_DIR="${BENCH_ROOT}/stringdb-annotation"
EDGE_OUTDIR="${BENCH_ROOT}/edge-level"
MODULE_OUTDIR="${BENCH_ROOT}/module-level"
MODULES_TSV_BY_METHOD="${BENCH_ROOT}/genescope/repeat_001/modules.tsv,${BENCH_ROOT}/giotto/repeat_001/modules.tsv,${BENCH_ROOT}/hotspot/repeat_001/modules.tsv,${BENCH_ROOT}/seagal/repeat_001/modules.tsv"

Rscript "${RSCRIPTS_DIR}/mapping.R" \
  --bench_root "${BENCH_ROOT}" \
  --outdir "${MAP_DIR}" \
  --methods "${METHODS}" \
  --genescope_root "${GENESCOPE_ROOT}" \
  --keep_subscores 1

Rscript "${RSCRIPTS_DIR}/edge-level.R" \
  --map_dir "${MAP_DIR}" \
  --outdir "${EDGE_OUTDIR}" \
  --methods "${METHODS}" \
  --refill_to_top_n_comparable 0 \
  --ranking_keys "weight" \
  --pr_top_n_list "10,30,50,100,1000" \
  --edge_fdr_by_method "1,1,1,1" \
  --string_score_threshold 700

Rscript "${RSCRIPTS_DIR}/module-level.R" \
  --map_dir "${MAP_DIR}" \
  --outdir "${MODULE_OUTDIR}" \
  --methods "${METHODS}" \
  --modules_tsv_by_method "${MODULES_TSV_BY_METHOD}" \
  --min_module_genes 3 \
  --n_random 200000 \
  --seed 1 \
  --min_valid_null_draws 80000 \
  --max_resample_attempts 20 \
  --plot_max_modules_per_method 200 \
  --string_score_threshold 700

#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: run.sh --data_dir DIR --outdir DIR [options]

Required:
  --data_dir       Xenium outs directory (container path)
  --outdir         Output directory (container path)

Optional:
  --coord_file     ROI CSV (optional)
  --roi_flip_y     Flip ROI Y coordinates (flag; default off)
  --ncores         Parallel cores (default 1)
  --threads        Alias of --ncores
  --seed           Random seed (default 1)
  --max_cells      Downsample to N cells (0 disables; default 0)

Xenium loader (PC-friendly defaults):
  --load_expression   Load 10x expression matrix (0/1; default 1)
  --load_cellmeta     Load 10x cell metadata (0/1; default 0)
  --load_transcripts  Load transcripts (0/1; default 0)

Workflow (grid-only):
  --coexpr_mode       Must be 'grid' (legacy; other values will error)

Grid-binning coexpression (defaults):
  --grid_stepsize         Grid stepsize (positive number); if omitted, runner auto-computes
  --spatial_grid_name     Spatial grid name (default spatial_grid)
  --min_cells_per_grid    Minimum cells per grid bin (default 4)
  --n_corr_genes          Genes used for correlation (default 500)
  --n_spatial_genes       Legacy alias of --n_corr_genes
  --k_modules             clusterSpatialCorFeats k (default 5)

Other:
  --cor_method     Correlation method: pearson (default pearson)
  --emit_edge_stats  Emit edge-level p_value/fdr (0/1; default 0)
  --sample_sec     Resource sampling seconds (default 1)
  --repeat         Repeat count (default 1)
  --dataset_id     Dataset identifier (optional)
  --roi_id         ROI identifier (optional)
USAGE
}

DATA_DIR=""
OUTDIR=""
COORD_FILE=""
ROI_FLIP_Y="0"
NCORES="1"
SEED="1"
MAX_CELLS="0"

LOAD_EXPRESSION="1"
LOAD_CELLMETA="0"
LOAD_TRANSCRIPTS="0"

COEXPR_MODE="grid"
GRID_STEPSIZE=""
SPATIAL_GRID_NAME="spatial_grid"
MIN_CELLS_PER_GRID="4"
PATTERN_DIMENSIONS="1:5"
TOP_POS_GENES="10"
TOP_NEG_GENES="10"
MIN_POS_COR="0.5"
MIN_NEG_COR="-0.5"
N_PATTERN_GENES_MAX="500"

K_NEIGHBORS="6"
MAXIMUM_DISTANCE_KNN="400"
BIN_METHOD="rank"
CALC_HUB="1"
HUB_MIN_INT="5"
N_CORR_GENES="500"
K_MODULES="5"
COR_METHOD="pearson"
EMIT_EDGE_STATS="0"

SAMPLE_SEC="1"
REPEAT="1"
DATASET_ID=""
ROI_ID=""
MON_PID=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data_dir)
      DATA_DIR="$2"; shift 2;;
    --data_dir=*)
      DATA_DIR="${1#*=}"; shift;;
    --outdir)
      OUTDIR="$2"; shift 2;;
    --outdir=*)
      OUTDIR="${1#*=}"; shift;;
    --coord_file)
      COORD_FILE="$2"; shift 2;;
    --coord_file=*)
      COORD_FILE="${1#*=}"; shift;;
    --roi_flip_y)
      ROI_FLIP_Y="1"; shift;;
    --roi_flip_y=*)
      ROI_FLIP_Y="${1#*=}"; shift;;
    --ncores)
      NCORES="$2"; shift 2;;
    --ncores=*)
      NCORES="${1#*=}"; shift;;
    --threads)
      NCORES="$2"; shift 2;;
    --threads=*)
      NCORES="${1#*=}"; shift;;
    --seed)
      SEED="$2"; shift 2;;
    --seed=*)
      SEED="${1#*=}"; shift;;
    --max_cells)
      MAX_CELLS="$2"; shift 2;;
    --max_cells=*)
      MAX_CELLS="${1#*=}"; shift;;
    --load_expression)
      LOAD_EXPRESSION="$2"; shift 2;;
    --load_expression=*)
      LOAD_EXPRESSION="${1#*=}"; shift;;
    --load_cellmeta)
      LOAD_CELLMETA="$2"; shift 2;;
    --load_cellmeta=*)
      LOAD_CELLMETA="${1#*=}"; shift;;
    --load_transcripts)
      LOAD_TRANSCRIPTS="$2"; shift 2;;
    --load_transcripts=*)
      LOAD_TRANSCRIPTS="${1#*=}"; shift;;
    --coexpr_mode)
      COEXPR_MODE="$2"; shift 2;;
    --coexpr_mode=*)
      COEXPR_MODE="${1#*=}"; shift;;
    --detect_method)
      COEXPR_MODE="$2"; shift 2;;
    --detect_method=*)
      COEXPR_MODE="${1#*=}"; shift;;
    --grid_stepsize)
      GRID_STEPSIZE="$2"; shift 2;;
    --grid_stepsize=*)
      GRID_STEPSIZE="${1#*=}"; shift;;
    --grid_um)
      GRID_STEPSIZE="$2"; shift 2;;
    --grid_um=*)
      GRID_STEPSIZE="${1#*=}"; shift;;
    --spatial_grid_um)
      GRID_STEPSIZE="$2"; shift 2;;
    --spatial_grid_um=*)
      GRID_STEPSIZE="${1#*=}"; shift;;
    --spatial_grid_name)
      SPATIAL_GRID_NAME="$2"; shift 2;;
    --spatial_grid_name=*)
      SPATIAL_GRID_NAME="${1#*=}"; shift;;
    --min_cells_per_grid)
      MIN_CELLS_PER_GRID="$2"; shift 2;;
    --min_cells_per_grid=*)
      MIN_CELLS_PER_GRID="${1#*=}"; shift;;
    --pattern_dimensions)
      PATTERN_DIMENSIONS="$2"; shift 2;;
    --pattern_dimensions=*)
      PATTERN_DIMENSIONS="${1#*=}"; shift;;
    --top_pos_genes)
      TOP_POS_GENES="$2"; shift 2;;
    --top_pos_genes=*)
      TOP_POS_GENES="${1#*=}"; shift;;
    --top_neg_genes)
      TOP_NEG_GENES="$2"; shift 2;;
    --top_neg_genes=*)
      TOP_NEG_GENES="${1#*=}"; shift;;
    --min_pos_cor)
      MIN_POS_COR="$2"; shift 2;;
    --min_pos_cor=*)
      MIN_POS_COR="${1#*=}"; shift;;
    --min_neg_cor)
      MIN_NEG_COR="$2"; shift 2;;
    --min_neg_cor=*)
      MIN_NEG_COR="${1#*=}"; shift;;
    --n_pattern_genes_max)
      N_PATTERN_GENES_MAX="$2"; shift 2;;
    --n_pattern_genes_max=*)
      N_PATTERN_GENES_MAX="${1#*=}"; shift;;
    --k_neighbors)
      K_NEIGHBORS="$2"; shift 2;;
    --k_neighbors=*)
      K_NEIGHBORS="${1#*=}"; shift;;
    --maximum_distance_knn)
      MAXIMUM_DISTANCE_KNN="$2"; shift 2;;
    --maximum_distance_knn=*)
      MAXIMUM_DISTANCE_KNN="${1#*=}"; shift;;
    --bin_method)
      BIN_METHOD="$2"; shift 2;;
    --bin_method=*)
      BIN_METHOD="${1#*=}"; shift;;
    --calc_hub)
      CALC_HUB="$2"; shift 2;;
    --calc_hub=*)
      CALC_HUB="${1#*=}"; shift;;
    --hub_min_int)
      HUB_MIN_INT="$2"; shift 2;;
	    --hub_min_int=*)
	      HUB_MIN_INT="${1#*=}"; shift;;
	    --n_corr_genes)
	      N_CORR_GENES="$2"; shift 2;;
	    --n_corr_genes=*)
	      N_CORR_GENES="${1#*=}"; shift;;
	    --n_spatial_genes)
	      N_CORR_GENES="$2"; shift 2;;
	    --n_spatial_genes=*)
	      N_CORR_GENES="${1#*=}"; shift;;
	    --k_modules)
	      K_MODULES="$2"; shift 2;;
	    --k_modules=*)
	      K_MODULES="${1#*=}"; shift;;
    --cor_method)
      COR_METHOD="$2"; shift 2;;
    --cor_method=*)
      COR_METHOD="${1#*=}"; shift;;
    --emit_edge_stats)
      EMIT_EDGE_STATS="$2"; shift 2;;
    --emit_edge_stats=*)
      EMIT_EDGE_STATS="${1#*=}"; shift;;
    --sample_sec)
      SAMPLE_SEC="$2"; shift 2;;
    --sample_sec=*)
      SAMPLE_SEC="${1#*=}"; shift;;
    --repeat)
      REPEAT="$2"; shift 2;;
    --repeat=*)
      REPEAT="${1#*=}"; shift;;
    --dataset_id)
      DATASET_ID="$2"; shift 2;;
    --dataset_id=*)
      DATASET_ID="${1#*=}"; shift;;
    --roi_id)
      ROI_ID="$2"; shift 2;;
    --roi_id=*)
      ROI_ID="${1#*=}"; shift;;
    -h|--help)
      usage; exit 0;;
    --*)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
    *)
      shift;;
  esac
done

if [[ -z "${DATA_DIR}" || -z "${OUTDIR}" ]]; then
  usage
  exit 2
fi

is_int() {
  [[ "$1" =~ ^-?[0-9]+$ ]]
}

is_num() {
  [[ "$1" =~ ^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]
}

sanitize_int() {
  local v="$1" default="$2"
  if is_int "${v}"; then
    echo "${v}"
  else
    echo "${default}"
  fi
}

sanitize_num() {
  local v="$1" default="$2"
  if is_num "${v}"; then
    echo "${v}"
  else
    echo "${default}"
  fi
}

sanitize_pos_num() {
  local v default
  default="$2"
  v="$(sanitize_num "$1" "$default")"
  awk -v v="$v" -v d="$default" 'BEGIN { if (v+0>0) print v; else print d }'
}

sanitize_nonneg_int() {
  local v
  v="$(sanitize_int "$1" "$2")"
  if [[ "${v}" -lt 0 ]]; then
    echo "$2"
  else
    echo "${v}"
  fi
}

sanitize_positive_int() {
  local v
  v="$(sanitize_int "$1" "$2")"
  if [[ "${v}" -lt 1 ]]; then
    echo "$2"
  else
    echo "${v}"
  fi
}

normalize_bool() {
  local v="${1:-}" default="${2:-0}"
  case "${v}" in
    1|true|TRUE|yes|YES|y|Y) echo "1";;
    0|false|FALSE|no|NO|n|N|"") echo "0";;
    *) echo "${default}";;
  esac
}

BIN_METHOD="${BIN_METHOD,,}"
case "${BIN_METHOD}" in
  rank|kmeans) ;;
  *) BIN_METHOD="rank" ;;
esac

COR_METHOD="${COR_METHOD,,}"
if [[ "${COR_METHOD}" != "pearson" ]]; then
  echo "This runner supports only --cor_method pearson." >&2
  exit 2
fi

COEXPR_MODE="${COEXPR_MODE,,}"
if [[ "${COEXPR_MODE}" != "grid" ]]; then
  echo "This runner is grid-only; --coexpr_mode must be 'grid' (got: ${COEXPR_MODE})." >&2
  exit 2
fi

if [[ -n "${GRID_STEPSIZE}" && "${GRID_STEPSIZE^^}" != "NA" ]]; then
  GRID_STEPSIZE="$(sanitize_pos_num "${GRID_STEPSIZE}" "INVALID")"
  if [[ "${GRID_STEPSIZE}" == "INVALID" ]]; then
    echo "Invalid --grid_stepsize (expected positive number or NA): ${GRID_STEPSIZE}" >&2
    exit 2
  fi
fi

if [[ -z "${SPATIAL_GRID_NAME}" ]]; then
  SPATIAL_GRID_NAME="spatial_grid"
fi

MIN_CELLS_PER_GRID="$(sanitize_positive_int "${MIN_CELLS_PER_GRID}" 4)"
TOP_POS_GENES="$(sanitize_positive_int "${TOP_POS_GENES}" 10)"
TOP_NEG_GENES="$(sanitize_positive_int "${TOP_NEG_GENES}" 10)"
N_PATTERN_GENES_MAX="$(sanitize_positive_int "${N_PATTERN_GENES_MAX}" 500)"
MIN_POS_COR="$(sanitize_num "${MIN_POS_COR}" 0.5)"
MIN_NEG_COR="$(sanitize_num "${MIN_NEG_COR}" -0.5)"
if [[ -z "${PATTERN_DIMENSIONS}" ]]; then
  PATTERN_DIMENSIONS="1:5"
fi

GRID_STEPSIZE_ARG=()
if [[ -n "${GRID_STEPSIZE}" ]]; then
  GRID_STEPSIZE_ARG=(--grid_stepsize "${GRID_STEPSIZE}")
fi

NCORES="$(sanitize_positive_int "${NCORES}" 1)"
SEED="$(sanitize_int "${SEED}" 1)"
MAX_CELLS="$(sanitize_nonneg_int "${MAX_CELLS}" 0)"

LOAD_EXPRESSION="$(normalize_bool "${LOAD_EXPRESSION}" 1)"
LOAD_CELLMETA="$(normalize_bool "${LOAD_CELLMETA}" 0)"
LOAD_TRANSCRIPTS="$(normalize_bool "${LOAD_TRANSCRIPTS}" 0)"
ROI_FLIP_Y="$(normalize_bool "${ROI_FLIP_Y}" 0)"

K_NEIGHBORS="$(sanitize_positive_int "${K_NEIGHBORS}" 6)"
MAXIMUM_DISTANCE_KNN="$(sanitize_pos_num "${MAXIMUM_DISTANCE_KNN}" 400)"
CALC_HUB="$(normalize_bool "${CALC_HUB}" 1)"
HUB_MIN_INT="$(sanitize_nonneg_int "${HUB_MIN_INT}" 5)"
N_CORR_GENES="$(sanitize_positive_int "${N_CORR_GENES}" 500)"
K_MODULES="$(sanitize_positive_int "${K_MODULES}" 5)"
EMIT_EDGE_STATS="$(normalize_bool "${EMIT_EDGE_STATS}" 0)"

SAMPLE_SEC="$(sanitize_nonneg_int "${SAMPLE_SEC}" 1)"
REPEAT="$(sanitize_positive_int "${REPEAT}" 1)"

export OMP_NUM_THREADS="${NCORES}"
export OPENBLAS_NUM_THREADS="${NCORES}"
export MKL_NUM_THREADS="${NCORES}"
export VECLIB_MAXIMUM_THREADS="${NCORES}"

mkdir -p "${OUTDIR}"

ROI_FLIP_Y_ARG=()
if [[ "${ROI_FLIP_Y}" != "0" ]]; then
  ROI_FLIP_Y_ARG=(--roi_flip_y)
fi

for rep_idx in $(seq 1 "${REPEAT}"); do
  repdir="${OUTDIR}/repeat_$(printf '%03d' "${rep_idx}")"
  mkdir -p "${repdir}"

  : > "${repdir}/run.log"
  STATS_TSV="${repdir}/stats.tsv"
  if [[ "${SAMPLE_SEC}" != "0" ]]; then
    /bin/bash /opt/app/monitor_cgroup.sh --out "${STATS_TSV}" --sample_sec "${SAMPLE_SEC}" &
    MON_PID=$!
  else
    : > "${STATS_TSV}"
  fi

  set +e
  micromamba run -n tool Rscript /opt/app/run_giotto_xenium.R \
    --data_dir "${DATA_DIR}" \
    --outdir "${repdir}" \
    --coord_file "${COORD_FILE}" \
    "${ROI_FLIP_Y_ARG[@]}" \
    --ncores "${NCORES}" \
    --seed "${SEED}" \
    --max_cells "${MAX_CELLS}" \
    --load_expression "${LOAD_EXPRESSION}" \
    --load_cellmeta "${LOAD_CELLMETA}" \
    --load_transcripts "${LOAD_TRANSCRIPTS}" \
    --coexpr_mode "${COEXPR_MODE}" \
    ${GRID_STEPSIZE_ARG[@]} \
    --spatial_grid_name "${SPATIAL_GRID_NAME}" \
    --min_cells_per_grid "${MIN_CELLS_PER_GRID}" \
    --pattern_dimensions "${PATTERN_DIMENSIONS}" \
    --top_pos_genes "${TOP_POS_GENES}" \
    --top_neg_genes "${TOP_NEG_GENES}" \
    --min_pos_cor "${MIN_POS_COR}" \
    --min_neg_cor "${MIN_NEG_COR}" \
    --n_pattern_genes_max "${N_PATTERN_GENES_MAX}" \
    --k_neighbors "${K_NEIGHBORS}" \
    --maximum_distance_knn "${MAXIMUM_DISTANCE_KNN}" \
	    --bin_method "${BIN_METHOD}" \
	    --calc_hub "${CALC_HUB}" \
	    --hub_min_int "${HUB_MIN_INT}" \
	    --n_corr_genes "${N_CORR_GENES}" \
	    --k_modules "${K_MODULES}" \
	    --cor_method "${COR_METHOD}" \
	    --emit_edge_stats "${EMIT_EDGE_STATS}" \
	    --dataset_id "${DATASET_ID}" \
    --roi_id "${ROI_ID}" \
    --stats_file "${STATS_TSV}" \
    >> "${repdir}/run.log" 2>&1
  RC="$?"
  set -e
  echo "[DONE] repeat_idx=${rep_idx}" >> "${repdir}/run.log"

  if [[ -n "${MON_PID:-}" ]]; then
    kill "${MON_PID}" >/dev/null 2>&1 || true
    wait "${MON_PID}" >/dev/null 2>&1 || true
    MON_PID=""
  fi

  if [[ "${RC}" -ne 0 ]]; then
    if [[ "${RC}" -eq 137 || "${RC}" -eq 143 ]]; then
      echo "[ERROR] Rscript exited with ${RC} (possible OOM/kill)" >> "${repdir}/run.log"
    else
      echo "[ERROR] Rscript exited with ${RC}" >> "${repdir}/run.log"
    fi
    exit "${RC}"
  fi
done

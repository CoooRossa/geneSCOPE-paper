#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: run.sh --data_dir DIR --outdir DIR --grid_stepsize NUM [options]

Required:
  --data_dir            Xenium outs directory (container path)
  --outdir              Output directory (container path)
  --grid_stepsize       Grid stepsize (positive number)

Optional (standard grid workflow):
  --coord_file          ROI CSV (optional; columns X/Y)
  --roi_flip_y          Flip ROI Y coordinates (flag)
  --threads, --ncores   Thread hint for BLAS/OMP (default 1)
  --seed                Random seed (optional)
  --sample_sec          Sample cgroup CPU/mem every N sec (default 0=disable)
  --spatial_grid_name    (default spatial_grid)
  --min_cells_per_grid   (default 4)
  --n_corr_genes         (default 500)
  --k_modules            (default 5)
  --expression_values    normalized|scaled|custom (default normalized)

Legacy benchmark flags (accepted but ignored; must be defaults):
  --max_cells 0
  --repeat 1
  --dataset_id <id>
  --roi_id <id>
USAGE
}

DATA_DIR=""
OUTDIR=""
GRID_STEPSIZE=""

COORD_FILE=""
ROI_FLIP_Y="0"
NCORES="1"
SEED=""
SAMPLE_SEC="0"

SPATIAL_GRID_NAME="spatial_grid"
MIN_CELLS_PER_GRID="4"
N_CORR_GENES="500"
K_MODULES="5"
EXPRESSION_VALUES="normalized"

MAX_CELLS="0"
REPEAT="1"
DATASET_ID=""
ROI_ID=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data_dir) DATA_DIR="$2"; shift 2;;
    --data_dir=*) DATA_DIR="${1#*=}"; shift;;
    --outdir) OUTDIR="$2"; shift 2;;
    --outdir=*) OUTDIR="${1#*=}"; shift;;
    --grid_stepsize) GRID_STEPSIZE="$2"; shift 2;;
    --grid_stepsize=*) GRID_STEPSIZE="${1#*=}"; shift;;
    --coord_file) COORD_FILE="$2"; shift 2;;
    --coord_file=*) COORD_FILE="${1#*=}"; shift;;
    --roi_flip_y) ROI_FLIP_Y="1"; shift;;
    --roi_flip_y=*) ROI_FLIP_Y="${1#*=}"; shift;;
    --threads|--ncores) NCORES="$2"; shift 2;;
    --threads=*|--ncores=*) NCORES="${1#*=}"; shift;;
    --seed) SEED="$2"; shift 2;;
    --seed=*) SEED="${1#*=}"; shift;;
    --sample_sec) SAMPLE_SEC="$2"; shift 2;;
    --sample_sec=*) SAMPLE_SEC="${1#*=}"; shift;;
    --spatial_grid_name) SPATIAL_GRID_NAME="$2"; shift 2;;
    --spatial_grid_name=*) SPATIAL_GRID_NAME="${1#*=}"; shift;;
    --min_cells_per_grid) MIN_CELLS_PER_GRID="$2"; shift 2;;
    --min_cells_per_grid=*) MIN_CELLS_PER_GRID="${1#*=}"; shift;;
    --n_corr_genes|--n_spatial_genes) N_CORR_GENES="$2"; shift 2;;
    --n_corr_genes=*|--n_spatial_genes=*) N_CORR_GENES="${1#*=}"; shift;;
    --k_modules) K_MODULES="$2"; shift 2;;
    --k_modules=*) K_MODULES="${1#*=}"; shift;;
    --expression_values) EXPRESSION_VALUES="$2"; shift 2;;
    --expression_values=*) EXPRESSION_VALUES="${1#*=}"; shift;;
    --max_cells) MAX_CELLS="$2"; shift 2;;
    --max_cells=*) MAX_CELLS="${1#*=}"; shift;;
    --repeat) REPEAT="$2"; shift 2;;
    --repeat=*) REPEAT="${1#*=}"; shift;;
    --dataset_id) DATASET_ID="$2"; shift 2;;
    --dataset_id=*) DATASET_ID="${1#*=}"; shift;;
    --roi_id) ROI_ID="$2"; shift 2;;
    --roi_id=*) ROI_ID="${1#*=}"; shift;;
    -h|--help) usage; exit 0;;
    --*) echo "Unknown option: $1" >&2; usage; exit 2;;
    *) shift;;
  esac
done

if [[ -z "${DATA_DIR}" || -z "${OUTDIR}" || -z "${GRID_STEPSIZE}" ]]; then
  usage
  exit 2
fi

if ! [[ "${NCORES}" =~ ^[0-9]+$ ]] || [[ "${NCORES}" -lt 1 ]]; then
  NCORES="1"
fi

export OMP_NUM_THREADS="${NCORES}"
export OPENBLAS_NUM_THREADS="${NCORES}"
export MKL_NUM_THREADS="${NCORES}"
export VECLIB_MAXIMUM_THREADS="${NCORES}"

if ! [[ "${MAX_CELLS}" =~ ^[0-9]+$ ]] || [[ "${MAX_CELLS}" -ne 0 ]]; then
  echo "ERROR: --max_cells is not supported by this minimal runner (must be 0)." >&2
  exit 2
fi

if ! [[ "${REPEAT}" =~ ^[0-9]+$ ]] || [[ "${REPEAT}" -ne 1 ]]; then
  echo "ERROR: --repeat is not supported by this minimal runner (must be 1)." >&2
  exit 2
fi

_unused="${DATASET_ID}${ROI_ID}"

if ! [[ "${GRID_STEPSIZE}" =~ ^[+-]?[0-9]*[.]?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
  echo "ERROR: --grid_stepsize must be a positive number." >&2
  exit 2
fi

if awk -v v="${GRID_STEPSIZE}" 'BEGIN{exit !(v+0>0)}'; then :; else
  echo "ERROR: --grid_stepsize must be > 0." >&2
  exit 2
fi

if [[ -z "${SPATIAL_GRID_NAME}" ]]; then
  SPATIAL_GRID_NAME="spatial_grid"
fi

if ! [[ "${MIN_CELLS_PER_GRID}" =~ ^[0-9]+$ ]] || [[ "${MIN_CELLS_PER_GRID}" -lt 1 ]]; then
  MIN_CELLS_PER_GRID="4"
fi

if ! [[ "${N_CORR_GENES}" =~ ^[0-9]+$ ]] || [[ "${N_CORR_GENES}" -lt 2 ]]; then
  N_CORR_GENES="500"
fi

if ! [[ "${K_MODULES}" =~ ^[0-9]+$ ]] || [[ "${K_MODULES}" -lt 1 ]]; then
  K_MODULES="5"
fi

EXPRESSION_VALUES="$(echo "${EXPRESSION_VALUES}" | tr '[:upper:]' '[:lower:]')"
case "${EXPRESSION_VALUES}" in
  normalized|scaled|custom) ;;
  *) echo "ERROR: --expression_values must be normalized|scaled|custom." >&2; exit 2;;
esac

ROI_FLIP_Y="$(echo "${ROI_FLIP_Y}" | tr '[:upper:]' '[:lower:]')"
case "${ROI_FLIP_Y}" in
  1|true|yes|y) ROI_FLIP_Y="1" ;;
  0|false|no|n|"") ROI_FLIP_Y="0" ;;
  *) ROI_FLIP_Y="0" ;;
esac

mkdir -p "${OUTDIR}"
REPDIR="${OUTDIR}/repeat_001"
mkdir -p "${REPDIR}"

MON_PID=""
cleanup_monitor() {
  if [[ -n "${MON_PID:-}" ]]; then
    kill "${MON_PID}" >/dev/null 2>&1 || true
    wait "${MON_PID}" >/dev/null 2>&1 || true
    MON_PID=""
  fi
}
trap cleanup_monitor EXIT

SAMPLE_SEC="$(echo "${SAMPLE_SEC}" | tr -d '[:space:]')"
if [[ -z "${SAMPLE_SEC}" ]]; then SAMPLE_SEC="0"; fi
if ! [[ "${SAMPLE_SEC}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "ERROR: --sample_sec must be a non-negative number." >&2
  exit 2
fi

STATS_TSV="${OUTDIR}/stats.tsv"
if [[ "${SAMPLE_SEC}" != "0" ]]; then
  /opt/app/monitor_cgroup.sh --out "${STATS_TSV}" --sample_sec "${SAMPLE_SEC}" &
  MON_PID=$!
fi

ARGS=(
  --data_dir "${DATA_DIR}"
  --outdir "${REPDIR}"
  --grid_stepsize "${GRID_STEPSIZE}"
  --spatial_grid_name "${SPATIAL_GRID_NAME}"
  --min_cells_per_grid "${MIN_CELLS_PER_GRID}"
  --n_corr_genes "${N_CORR_GENES}"
  --k_modules "${K_MODULES}"
  --expression_values "${EXPRESSION_VALUES}"
)

if [[ -n "${SEED}" ]]; then
  ARGS+=(--seed "${SEED}")
fi

if [[ -n "${COORD_FILE}" ]]; then
  ARGS+=(--coord_file "${COORD_FILE}")
  if [[ "${ROI_FLIP_Y}" == "1" ]]; then
    ARGS+=(--roi_flip_y)
  fi
fi

exec micromamba run -n tool Rscript /opt/app/run_giotto_xenium.R "${ARGS[@]}"

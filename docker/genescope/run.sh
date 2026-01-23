#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=""
COORD_FILE=""
OUTDIR=""
GRID_UM="30"
NCORES="8"
SEED="1"
SAMPLE_SEC="0"
PARALLEL_BACKEND="serial"
DATASET_ID=""
ROI_ID=""
MON_PID=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data_dir) DATA_DIR="$2"; shift 2;;
    --coord_file) COORD_FILE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --grid_um) GRID_UM="$2"; shift 2;;
    --ncores) NCORES="$2"; shift 2;;
    --threads) NCORES="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --sample_sec) SAMPLE_SEC="$2"; shift 2;;
    --parallel_backend) PARALLEL_BACKEND="$2"; shift 2;;
    --dataset_id) DATASET_ID="$2"; shift 2;;
    --roi_id) ROI_ID="$2"; shift 2;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

is_int() {
  [[ "${1:-}" =~ ^-?[0-9]+$ ]]
}

sanitize_int() {
  local v="${1:-}" default="$2"
  if is_int "${v}"; then
    echo "${v}"
  else
    echo "${default}"
  fi
}

sanitize_nonneg_int() {
  local v
  v="$(sanitize_int "${1:-}" "$2")"
  if [[ "${v}" -lt 0 ]]; then
    echo "$2"
  else
    echo "${v}"
  fi
}

sanitize_positive_int() {
  local v
  v="$(sanitize_int "${1:-}" "$2")"
  if [[ "${v}" -lt 1 ]]; then
    echo "$2"
  else
    echo "${v}"
  fi
}

if [[ -z "${DATA_DIR}" || -z "${OUTDIR}" ]]; then
  echo "Usage: --data_dir <XENIUM_OUTS_DIR> --outdir <OUTDIR> [--coord_file <ROI.csv>] [--grid_um 30] [--ncores 8|--threads 8] [--seed 1] [--sample_sec 1] [--parallel_backend serial] [--dataset_id <id>] [--roi_id <id>]" >&2
  exit 2
fi

GRID_UM="$(sanitize_positive_int "${GRID_UM}" 30)"
NCORES="$(sanitize_positive_int "${NCORES}" 8)"
SEED="$(sanitize_nonneg_int "${SEED}" 1)"
SAMPLE_SEC="$(sanitize_nonneg_int "${SAMPLE_SEC}" 0)"
if [[ -z "${PARALLEL_BACKEND}" ]]; then PARALLEL_BACKEND="serial"; fi

mkdir -p "${OUTDIR}"
REPDIR="${OUTDIR}/repeat_001"
mkdir -p "${REPDIR}"
LOG="${OUTDIR}/run.log"
: > "${LOG}"

cleanup_monitor() {
  if [[ -n "${MON_PID:-}" ]]; then
    kill "${MON_PID}" >/dev/null 2>&1 || true
    wait "${MON_PID}" >/dev/null 2>&1 || true
    MON_PID=""
  fi
}
trap cleanup_monitor EXIT

STATS_TSV="${OUTDIR}/stats.tsv"
if [[ "${SAMPLE_SEC}" != "0" ]]; then
  /opt/app/monitor_cgroup.sh --out "${STATS_TSV}" --sample_sec "${SAMPLE_SEC}" &
  MON_PID=$!
fi

micromamba run -n tool Rscript /opt/app/run_genescope_xenium.R \
  "${DATA_DIR}" "${REPDIR}" "${GRID_UM}" "${NCORES}" "${SEED}" "${PARALLEL_BACKEND}" "${COORD_FILE}" \
  "${DATASET_ID}" "${ROI_ID}" "${OUTDIR}/stats.tsv" \
  > "${LOG}" 2>&1

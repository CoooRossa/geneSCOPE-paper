#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: run.sh --data_dir DIR --outdir DIR [options]

Required:
  --data_dir             Xenium outs directory (container path) or .h5ad
  --outdir               Output directory (container path)

ROI (choose one):
  --coord_file           ROI CSV (QuPath format, optional)
  --roi_csv              Alias of --coord_file

Optional:
  --max_cells            Downsample cap (0 disables, default 0)
  --seed                 Random seed (default 1)
  --n_neighbors          Hotspot KNN neighbors (default 300)
  --model                Hotspot model (default bernoulli)
  --fdr_autocorr         FDR cutoff for autocorr (default 0.05)
  --min_gene_threshold   Min genes per module (default 20)
  --k_modules            Force module count via linkage cut (0 disables, default 0)
  --emit_gg_matrices     Emit full gg matrices (0/1, default 0)
  --core_only            Core-only modules (flag or 0/1, default false)
  --write_module_scores  Write module score matrix (flag or 0/1, default 1)
  --repeat               Repeat count (default 1)
  --sample_sec           Cgroup sample seconds (default 1)
  --ncores               Parallel cores (default 8)
  --threads              Alias of --ncores
  --dataset_id           Dataset identifier (optional)
  --roi_id               ROI identifier (optional)
USAGE
}

DATA_DIR=""
OUTDIR=""
COORD_FILE=""
MAX_CELLS="0"
SEED="1"
N_NEIGHBORS="300"
MODEL="bernoulli"
FDR_AUTOCORR="0.05"
MIN_GENE_THRESHOLD="20"
K_MODULES="0"
EMIT_GG_MATRICES="0"
CORE_ONLY="false"
WRITE_MODULE_SCORES="1"
REPEAT="1"
SAMPLE_SEC="1"
NCORES="8"
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
    --roi_csv)
      COORD_FILE="$2"; shift 2;;
    --roi_csv=*)
      COORD_FILE="${1#*=}"; shift;;
    --max_cells)
      MAX_CELLS="$2"; shift 2;;
    --max_cells=*)
      MAX_CELLS="${1#*=}"; shift;;
    --seed)
      SEED="$2"; shift 2;;
    --seed=*)
      SEED="${1#*=}"; shift;;
    --n_neighbors)
      N_NEIGHBORS="$2"; shift 2;;
    --n_neighbors=*)
      N_NEIGHBORS="${1#*=}"; shift;;
    --model)
      MODEL="$2"; shift 2;;
    --model=*)
      MODEL="${1#*=}"; shift;;
    --fdr_autocorr)
      FDR_AUTOCORR="$2"; shift 2;;
    --fdr_autocorr=*)
      FDR_AUTOCORR="${1#*=}"; shift;;
    --min_gene_threshold)
      MIN_GENE_THRESHOLD="$2"; shift 2;;
    --min_gene_threshold=*)
      MIN_GENE_THRESHOLD="${1#*=}"; shift;;
    --k_modules)
      K_MODULES="$2"; shift 2;;
    --k_modules=*)
      K_MODULES="${1#*=}"; shift;;
    --emit_gg_matrices)
      EMIT_GG_MATRICES="$2"; shift 2;;
    --emit_gg_matrices=*)
      EMIT_GG_MATRICES="${1#*=}"; shift;;
    --core_only)
      if [[ $# -ge 2 && ! "$2" =~ ^-- ]]; then
        CORE_ONLY="$2"; shift 2
      else
        CORE_ONLY="true"; shift
      fi
      ;;
    --core_only=*)
      CORE_ONLY="${1#*=}"; shift;;
    --write_module_scores)
      if [[ $# -ge 2 && ! "$2" =~ ^-- ]]; then
        WRITE_MODULE_SCORES="$2"; shift 2
      else
        WRITE_MODULE_SCORES="1"; shift
      fi
      ;;
    --write_module_scores=*)
      WRITE_MODULE_SCORES="${1#*=}"; shift;;
    --repeat)
      REPEAT="$2"; shift 2;;
    --repeat=*)
      REPEAT="${1#*=}"; shift;;
    --sample_sec)
      SAMPLE_SEC="$2"; shift 2;;
    --sample_sec=*)
      SAMPLE_SEC="${1#*=}"; shift;;
    --ncores)
      NCORES="$2"; shift 2;;
    --ncores=*)
      NCORES="${1#*=}"; shift;;
    --threads)
      NCORES="$2"; shift 2;;
    --threads=*)
      NCORES="${1#*=}"; shift;;
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
      if [[ "$1" == *=* ]]; then
        shift
      else
        if [[ $# -ge 2 && ! "$2" =~ ^-- ]]; then
          shift 2
        else
          shift
        fi
      fi
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

sanitize_int() {
  local v="$1" default="$2"
  if is_int "${v}"; then
    echo "${v}"
  else
    echo "${default}"
  fi
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

NCORES="$(sanitize_positive_int "${NCORES}" 8)"
SEED="$(sanitize_int "${SEED}" 1)"
SAMPLE_SEC="$(sanitize_nonneg_int "${SAMPLE_SEC}" 1)"
REPEAT="$(sanitize_positive_int "${REPEAT}" 1)"
MAX_CELLS="$(sanitize_nonneg_int "${MAX_CELLS}" 0)"
N_NEIGHBORS="$(sanitize_positive_int "${N_NEIGHBORS}" 300)"
MIN_GENE_THRESHOLD="$(sanitize_nonneg_int "${MIN_GENE_THRESHOLD}" 20)"
K_MODULES="$(sanitize_nonneg_int "${K_MODULES}" 0)"
EMIT_GG_MATRICES="$(sanitize_int "${EMIT_GG_MATRICES}" 0)"
if [[ "${EMIT_GG_MATRICES}" != "0" ]]; then
  EMIT_GG_MATRICES="1"
fi
CORE_ONLY="$(normalize_bool "${CORE_ONLY}" 0)"
WRITE_MODULE_SCORES="$(normalize_bool "${WRITE_MODULE_SCORES}" 1)"

mkdir -p "${OUTDIR}"

write_versions() {
  local out="$1" log_file="$2"
  if ! micromamba run -n tool python - <<'PY' > "${out}" 2>> "${log_file}"; then
import sys
from importlib.metadata import version

def v(pkg):
    try:
        return version(pkg)
    except Exception:
        return "unknown"

try:
    import numpy, pandas, scipy
    npv = numpy.__version__
    pdv = pandas.__version__
    spv = scipy.__version__
except Exception:
    npv = pdv = spv = "unknown"

try:
    import anndata
    adv = anndata.__version__
except Exception:
    adv = "unknown"

try:
    import pyarrow
    pav = pyarrow.__version__
except Exception:
    pav = "unknown"

try:
    import shapely
    shv = shapely.__version__
except Exception:
    shv = "unknown"

try:
    import sklearn
    skv = sklearn.__version__
except Exception:
    skv = "unknown"

try:
    import statsmodels
    stv = statsmodels.__version__
except Exception:
    stv = "unknown"

try:
    import h5py
    h5v = h5py.__version__
except Exception:
    h5v = "unknown"

lines = [
    f"python={sys.version.split()[0]}",
    f"hotspot={v('hotspot')}",
    f"numpy={npv}",
    f"pandas={pdv}",
    f"scipy={spv}",
    f"anndata={adv}",
    f"pyarrow={pav}",
    f"shapely={shv}",
    f"scikit-learn={skv}",
    f"statsmodels={stv}",
    f"h5py={h5v}",
]
print("\n".join(lines))
PY
    echo "versions_failed=1" > "${out}"
  fi
}

write_run_args() {
  local out="$1" rep_idx="$2" log_file="$3"
  if ! DATA_DIR="${DATA_DIR}" OUTDIR="${OUTDIR}" COORD_FILE="${COORD_FILE}" \
      NCORES="${NCORES}" SEED="${SEED}" SAMPLE_SEC="${SAMPLE_SEC}" REPEAT="${REPEAT}" \
      REPEAT_IDX="${rep_idx}" MAX_CELLS="${MAX_CELLS}" N_NEIGHBORS="${N_NEIGHBORS}" \
      MODEL="${MODEL}" FDR_AUTOCORR="${FDR_AUTOCORR}" MIN_GENE_THRESHOLD="${MIN_GENE_THRESHOLD}" \
      K_MODULES="${K_MODULES}" EMIT_GG_MATRICES="${EMIT_GG_MATRICES}" \
      CORE_ONLY="${CORE_ONLY}" WRITE_MODULE_SCORES="${WRITE_MODULE_SCORES}" RUN_ARGS_JSON="${out}" \
      micromamba run -n tool python - <<'PY' 2>> "${log_file}"; then
import json
import os

out = os.environ.get("RUN_ARGS_JSON")

payload = {
    "data_dir": os.environ.get("DATA_DIR"),
    "outdir": os.environ.get("OUTDIR"),
    "coord_file": os.environ.get("COORD_FILE"),
    "ncores": int(os.environ.get("NCORES", "1")),
    "seed": int(os.environ.get("SEED", "1")),
    "sample_sec": int(os.environ.get("SAMPLE_SEC", "1")),
    "repeat": int(os.environ.get("REPEAT", "1")),
    "repeat_idx": int(os.environ.get("REPEAT_IDX", "1")),
    "max_cells": int(os.environ.get("MAX_CELLS", "0")),
    "n_neighbors": int(os.environ.get("N_NEIGHBORS", "300")),
    "model": os.environ.get("MODEL"),
    "fdr_autocorr": float(os.environ.get("FDR_AUTOCORR", "0.05")),
    "min_gene_threshold": int(os.environ.get("MIN_GENE_THRESHOLD", "20")),
    "k_modules": int(os.environ.get("K_MODULES", "0")),
    "emit_gg_matrices": bool(int(os.environ.get("EMIT_GG_MATRICES", "0"))),
    "core_only": bool(int(os.environ.get("CORE_ONLY", "0"))),
    "write_module_scores": bool(int(os.environ.get("WRITE_MODULE_SCORES", "1"))),
}

with open(out, "w", encoding="utf-8") as f:
    json.dump(payload, f, indent=2, ensure_ascii=False)
PY
    printf '{"data_dir":"%s","outdir":"%s","coord_file":"%s","ncores":%s,"seed":%s,"sample_sec":%s,"repeat":%s,"repeat_idx":%s,"max_cells":%s,"n_neighbors":%s,"model":"%s","fdr_autocorr":%s,"min_gene_threshold":%s,"k_modules":%s,"emit_gg_matrices":%s,"core_only":%s,"write_module_scores":%s}\n' \
      "${DATA_DIR}" "${OUTDIR}" "${COORD_FILE}" "${NCORES}" "${SEED}" "${SAMPLE_SEC}" "${REPEAT}" "${rep_idx}" \
      "${MAX_CELLS}" "${N_NEIGHBORS}" "${MODEL}" "${FDR_AUTOCORR}" "${MIN_GENE_THRESHOLD}" "${K_MODULES}" "${EMIT_GG_MATRICES}" "${CORE_ONLY}" "${WRITE_MODULE_SCORES}" \
      > "${out}"
    echo "[WARN] run_args.json written without python json" >> "${log_file}"
  fi
}

TIME_BIN="$(command -v time || true)"
if [[ "${TIME_BIN}" == "time" && -x "/usr/bin/time" ]]; then
  TIME_BIN="/usr/bin/time"
fi

for rep_idx in $(seq 1 "${REPEAT}"); do
  repdir="${OUTDIR}/repeat_$(printf '%03d' "${rep_idx}")"
  mkdir -p "${repdir}"

  : > "${repdir}/run.log"
  STATS_TSV="${repdir}/stats.tsv"
  if [[ "${SAMPLE_SEC}" != "0" ]]; then
    /bin/bash /opt/app/monitor_cgroup.sh --out "${STATS_TSV}" --sample_sec "${SAMPLE_SEC}" &
    MON_PID=$!
  fi

  set +e
  micromamba run -n tool python /opt/app/run_hotspot_xenium.py \
    --data_dir "${DATA_DIR}" \
    --outdir "${repdir}" \
    --coord_file "${COORD_FILE}" \
    --max_cells "${MAX_CELLS}" \
    --seed "${SEED}" \
    --n_neighbors "${N_NEIGHBORS}" \
    --model "${MODEL}" \
    --fdr_autocorr "${FDR_AUTOCORR}" \
    --min_gene_threshold "${MIN_GENE_THRESHOLD}" \
    --k_modules "${K_MODULES}" \
    --emit_gg_matrices "${EMIT_GG_MATRICES}" \
    --core_only "${CORE_ONLY}" \
    --write_module_scores "${WRITE_MODULE_SCORES}" \
    --ncores "${NCORES}" \
    --dataset_id "${DATASET_ID}" \
    --roi_id "${ROI_ID}" \
    --stats_file "${STATS_TSV}" \
    >> "${repdir}/run.log" 2>&1
  RC="$?"
  set -e

  if [[ -n "${MON_PID:-}" ]]; then
    kill "${MON_PID}" >/dev/null 2>&1 || true
    wait "${MON_PID}" >/dev/null 2>&1 || true
    MON_PID=""
  fi

  if [[ "${RC}" -ne 0 ]]; then
    exit "${RC}"
  fi

done

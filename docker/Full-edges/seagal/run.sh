#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: run.sh --data_dir DIR --outdir DIR [options]

Required:
  --data_dir       Xenium outs directory (container path)
  --outdir         Output directory (container path)

Optional:
  --coord_file       ROI CSV (optional)
  --ncores           Parallel cores (default 1)
  --threads          Alias of --ncores
  --seed             Random seed (default 1)
  --sample_sec       Resource sampling seconds (default 1)
  --repeat           Repeat count (default 1)
  --max_cells        Downsample cap (0 disables)
  --top_genes        Keep top N genes by mean (0 disables)
  --grid_um          Grid binning step size (um; 0 disables)
  --modules_nmax     Max candidate module count (pattern-gene modules; optional)
  --extra_args_json  Extra JSON string or file path
USAGE
}

DATA_DIR=""
OUTDIR=""
COORD_FILE=""
NCORES="1"
SEED="1"
SAMPLE_SEC="1"
REPEAT="1"
MAX_CELLS="0"
TOP_GENES="0"
GRID_UM="0"
MODULES_NMAX=""
EXTRA_ARGS_JSON=""

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
    --sample_sec)
      SAMPLE_SEC="$2"; shift 2;;
    --sample_sec=*)
      SAMPLE_SEC="${1#*=}"; shift;;
    --repeat)
      REPEAT="$2"; shift 2;;
    --repeat=*)
      REPEAT="${1#*=}"; shift;;
    --max_cells)
      MAX_CELLS="$2"; shift 2;;
    --max_cells=*)
      MAX_CELLS="${1#*=}"; shift;;
    --top_genes)
      TOP_GENES="$2"; shift 2;;
    --top_genes=*)
      TOP_GENES="${1#*=}"; shift;;
    --grid_um)
      GRID_UM="$2"; shift 2;;
    --grid_um=*)
      GRID_UM="${1#*=}"; shift;;
    --modules_nmax)
      MODULES_NMAX="$2"; shift 2;;
    --modules_nmax=*)
      MODULES_NMAX="${1#*=}"; shift;;
    --extra_args_json)
      EXTRA_ARGS_JSON="$2"; shift 2;;
    --extra_args_json=*)
      EXTRA_ARGS_JSON="${1#*=}"; shift;;
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

NCORES="$(sanitize_positive_int "${NCORES}" 1)"
SEED="$(sanitize_int "${SEED}" 1)"
SAMPLE_SEC="$(sanitize_nonneg_int "${SAMPLE_SEC}" 1)"
REPEAT="$(sanitize_positive_int "${REPEAT}" 1)"
MAX_CELLS="$(sanitize_nonneg_int "${MAX_CELLS}" 0)"
TOP_GENES="$(sanitize_nonneg_int "${TOP_GENES}" 0)"
if [[ -n "${MODULES_NMAX}" ]]; then
  MODULES_NMAX="$(sanitize_nonneg_int "${MODULES_NMAX}" 0)"
fi

mkdir -p "${OUTDIR}"
LOG="${OUTDIR}/run.log"
: > "${LOG}"

log() {
  echo "[$(date -Iseconds)] $*" >> "${LOG}"
}

MON_PID=""
cleanup() {
  if [[ -n "${MON_PID}" ]]; then
    kill "${MON_PID}" >/dev/null 2>&1 || true
    wait "${MON_PID}" >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT

write_resource_limits() {
  local out="$1"
  {
    echo "date=$(date -Iseconds)"
    echo "cgroup_root=/sys/fs/cgroup"
    if [[ -f /sys/fs/cgroup/memory.max ]]; then
      echo "memory.max=$(cat /sys/fs/cgroup/memory.max 2>/dev/null || echo NA)"
      echo "memory.current=$(cat /sys/fs/cgroup/memory.current 2>/dev/null || echo NA)"
    elif [[ -f /sys/fs/cgroup/memory/memory.limit_in_bytes ]]; then
      echo "memory.limit_in_bytes=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo NA)"
      echo "memory.usage_in_bytes=$(cat /sys/fs/cgroup/memory/memory.usage_in_bytes 2>/dev/null || echo NA)"
    fi
    if [[ -f /sys/fs/cgroup/cpu.max ]]; then
      echo "cpu.max=$(cat /sys/fs/cgroup/cpu.max 2>/dev/null || echo NA)"
      if [[ -f /sys/fs/cgroup/cpu.stat ]]; then
        echo "cpu.stat=$(tr '\n' ' ' < /sys/fs/cgroup/cpu.stat 2>/dev/null || echo NA)"
      fi
    elif [[ -f /sys/fs/cgroup/cpu/cpu.cfs_quota_us ]]; then
      echo "cpu.cfs_quota_us=$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us 2>/dev/null || echo NA)"
      echo "cpu.cfs_period_us=$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us 2>/dev/null || echo NA)"
    fi
  } > "${out}" 2>&1 || true
}

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

lines = [
    f"python={sys.version.split()[0]}",
    f"seagal={v('seagal')}",
    f"numpy={npv}",
    f"pandas={pdv}",
    f"scipy={spv}",
]
print("\n".join(lines))

try:
    import subprocess
    res = subprocess.run([sys.executable, "-m", "pip", "show", "seagal"], capture_output=True, text=True, check=False)
    if res.stdout.strip():
        print("pip_show_seagal:")
        print(res.stdout.strip())
except Exception:
    pass
PY
    echo "versions_failed=1" > "${out}"
  fi
}

write_run_args() {
  local out="$1" rep_idx="$2" log_file="$3"
  if ! DATA_DIR="${DATA_DIR}" OUTDIR="${OUTDIR}" COORD_FILE="${COORD_FILE}" \
      NCORES="${NCORES}" SEED="${SEED}" SAMPLE_SEC="${SAMPLE_SEC}" REPEAT="${REPEAT}" \
      REPEAT_IDX="${rep_idx}" MAX_CELLS="${MAX_CELLS}" TOP_GENES="${TOP_GENES}" \
      GRID_UM="${GRID_UM}" MODULES_NMAX="${MODULES_NMAX}" EXTRA_ARGS_JSON="${EXTRA_ARGS_JSON}" RUN_ARGS_JSON="${out}" \
      micromamba run -n tool python - <<'PY' 2>> "${log_file}"; then
import json
import os

out = os.environ.get("RUN_ARGS_JSON")
extra_raw = os.environ.get("EXTRA_ARGS_JSON", "")
extra_val = None
extra_err = None
if extra_raw:
    try:
        if os.path.exists(extra_raw):
            with open(extra_raw, "r", encoding="utf-8") as f:
                extra_val = json.load(f)
        else:
            extra_val = json.loads(extra_raw)
    except Exception as exc:
        extra_err = str(exc)
        extra_val = extra_raw

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
    "top_genes": int(os.environ.get("TOP_GENES", "0")),
    "grid_um": os.environ.get("GRID_UM"),
    "modules_nmax": (
        int(os.environ.get("MODULES_NMAX", "0")) if os.environ.get("MODULES_NMAX") else None
    ),
    "extra_args_json": extra_val,
}
if extra_err:
    payload["extra_args_error"] = extra_err

	with open(out, "w", encoding="utf-8") as f:
	    json.dump(payload, f, indent=2, ensure_ascii=False)
PY
    modules_nmax_json="null"
    if [[ -n "${MODULES_NMAX}" ]]; then
      modules_nmax_json="${MODULES_NMAX}"
    fi
    printf '{"data_dir":"%s","outdir":"%s","coord_file":"%s","ncores":%s,"seed":%s,"sample_sec":%s,"repeat":%s,"repeat_idx":%s,"max_cells":%s,"top_genes":%s,"grid_um":"%s","modules_nmax":%s,"extra_args_json":"%s"}\n' \
      "${DATA_DIR}" "${OUTDIR}" "${COORD_FILE}" "${NCORES}" "${SEED}" "${SAMPLE_SEC}" "${REPEAT}" "${rep_idx}" "${MAX_CELLS}" "${TOP_GENES}" "${GRID_UM}" "${modules_nmax_json}" "${EXTRA_ARGS_JSON}" \
      > "${out}"
    echo "[WARN] run_args.json written without python json" >> "${log_file}"
  fi
}

write_resource_limits "${OUTDIR}/resource_limits.txt"
write_versions "${OUTDIR}/version.txt" "${LOG}"

if [[ ! -f "${OUTDIR}/stats.tsv" ]]; then
  echo -e "epoch\tiso_time\tmem_current_bytes\tmem_peak_bytes\tcpu_usage_usec\tcpu_usage_delta_usec" > "${OUTDIR}/stats.tsv"
fi

if [[ "${SAMPLE_SEC}" != "0" && -f /opt/app/monitor_cgroup.sh ]]; then
  /bin/bash /opt/app/monitor_cgroup.sh --out "${OUTDIR}/stats.tsv" --sample_sec "${SAMPLE_SEC}" \
    >> "${LOG}" 2>&1 &
  MON_PID="$!"
  sleep 0.1
  if ! kill -0 "${MON_PID}" 2>/dev/null; then
    log "[WARN] monitor_cgroup.sh failed to start"
    MON_PID=""
  fi
else
  log "[WARN] monitor_cgroup.sh not found or disabled"
fi

TIME_BIN="$(command -v time || true)"
if [[ "${TIME_BIN}" == "time" && -x "/usr/bin/time" ]]; then
  TIME_BIN="/usr/bin/time"
fi
TIME_CMD=()
if [[ -n "${TIME_BIN}" && -x "${TIME_BIN}" ]]; then
  TIME_CMD=("${TIME_BIN}" -v -o "${OUTDIR}/time.txt" -a)
fi

for rep_idx in $(seq 1 "${REPEAT}"); do
  repdir="${OUTDIR}/repeat_$(printf '%03d' "${rep_idx}")"
  mkdir -p "${repdir}"

  : > "${repdir}/run.log"
  log "[INFO] repeat_idx=${rep_idx}"

  write_run_args "${repdir}/run_args.json" "${rep_idx}" "${repdir}/run.log"

  START_EPOCH="$(date +%s)"

  MODULES_NMAX_ARGS=()
  if [[ -n "${MODULES_NMAX}" ]]; then
    MODULES_NMAX_ARGS=(--modules_nmax "${MODULES_NMAX}")
  fi

  set +e
  "${TIME_CMD[@]}" \
    micromamba run -n tool python /opt/app/run_seagal_xenium.py \
      --data_dir "${DATA_DIR}" \
      --outdir "${repdir}" \
      --coord_file "${COORD_FILE}" \
      --ncores "${NCORES}" \
      --seed "${SEED}" \
      --max_cells "${MAX_CELLS}" \
      --top_genes "${TOP_GENES}" \
      --grid_um "${GRID_UM}" \
      "${MODULES_NMAX_ARGS[@]}" \
      --extra_args_json "${EXTRA_ARGS_JSON}" \
      >> "${repdir}/run.log" 2>&1
  RC="$?"
  set -e

  END_EPOCH="$(date +%s)"
  WALL_SEC="$((END_EPOCH - START_EPOCH))"

  if [[ "${RC}" -ne 0 ]]; then
    if [[ "${RC}" -eq 137 ]]; then
      msg="[ERROR] SEAGAL runner killed (exit 137, likely OOM). Reduce n_obs via --max_cells and/or use --grid_um; see ${OUTDIR}/time.txt and ${repdir}/run.log"
      echo "${msg}" >&2
      echo "${msg}" >> "${repdir}/run.log"
      log "${msg}"
    fi
    log "[FAIL] data_dir=${DATA_DIR} outdir=${OUTDIR} repeat_idx=${rep_idx} seed=${SEED}"
    exit "${RC}"
  fi

  n_cells="NA"
  n_genes="NA"
  n_sig="NA"
  roi_applied="NA"
  modules_enabled="1"
  n_modules="NA"
  module_k_best="NA"

  if [[ ! -s "${repdir}/run_stats.tsv" ]]; then
    echo "[ERROR] Missing or empty run_stats.tsv: ${repdir}/run_stats.tsv" >&2
    log "[ERROR] Missing or empty run_stats.tsv: ${repdir}/run_stats.tsv"
    exit 1
  fi

  if [[ -f "${repdir}/run_stats.tsv" ]]; then
    while IFS=$'\t' read -r key val; do
      case "${key}" in
        n_cells) n_cells="${val}";;
        n_genes) n_genes="${val}";;
        n_sig) n_sig="${val}";;
        roi_applied) roi_applied="${val}";;
        modules_enabled) modules_enabled="${val}";;
        n_modules) n_modules="${val}";;
        module_k_best) module_k_best="${val}";;
      esac
    done < "${repdir}/run_stats.tsv"
  fi

  if [[ ! -s "${repdir}/edges_all.tsv" ]]; then
    echo "[ERROR] Missing or empty edges_all.tsv: ${repdir}/edges_all.tsv" >&2
    log "[ERROR] Missing or empty edges_all.tsv: ${repdir}/edges_all.tsv"
    exit 1
  fi

  if [[ "${modules_enabled}" != "0" ]]; then
    if [[ ! -s "${repdir}/modules.tsv" ]]; then
      echo "[ERROR] Missing or empty modules.tsv: ${repdir}/modules.tsv (modules_enabled=${modules_enabled})" >&2
      log "[ERROR] Missing or empty modules.tsv: ${repdir}/modules.tsv (modules_enabled=${modules_enabled})"
      exit 1
    fi
  fi

  {
    echo -e "wall_time_sec\t${WALL_SEC}"
    echo -e "n_cells\t${n_cells}"
    echo -e "n_genes\t${n_genes}"
    echo -e "n_sig\t${n_sig}"
    echo -e "seed\t${SEED}"
    echo -e "roi_applied\t${roi_applied}"
    echo -e "max_cells\t${MAX_CELLS}"
    echo -e "top_genes\t${TOP_GENES}"
    echo -e "repeat_idx\t${rep_idx}"
    echo -e "modules_enabled\t${modules_enabled}"
    echo -e "n_modules\t${n_modules}"
    echo -e "module_k_best\t${module_k_best}"
  } > "${repdir}/metrics.tsv"

  {
    echo -e "n_cells\t${n_cells}"
    echo -e "n_genes\t${n_genes}"
    echo -e "n_sig\t${n_sig}"
    echo -e "roi_applied\t${roi_applied}"
    echo -e "modules_enabled\t${modules_enabled}"
    echo -e "n_modules\t${n_modules}"
    echo -e "module_k_best\t${module_k_best}"
  } > "${repdir}/summary.tsv"
done

if [[ -n "${MON_PID}" ]]; then
  kill "${MON_PID}" >/dev/null 2>&1 || true
  wait "${MON_PID}" >/dev/null 2>&1 || true
  MON_PID=""
fi

if [[ ! -f "${OUTDIR}/stats.tsv" ]]; then
  echo -e "epoch\tiso_time\tmem_current_bytes\tmem_peak_bytes\tcpu_usage_usec\tcpu_usage_delta_usec" > "${OUTDIR}/stats.tsv"
fi

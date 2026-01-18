#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNTIME_DOCKER_DIR="${SCRIPT_DIR}/docker/runtime"

# ---- build runtime images ----
docker build -t genescope:final-rt "${RUNTIME_DOCKER_DIR}/genescope"
docker build -t giotto:final10x-rt  "${RUNTIME_DOCKER_DIR}/giotto_10x"
docker build -t hotspot:final-rt   "${RUNTIME_DOCKER_DIR}/hotspot"
docker build -t seagal:final-rt    "${RUNTIME_DOCKER_DIR}/seagal"

# ---- EDIT PATHS (host) ----
OUTS="/path/to/xenium_outs"
ROI_CSV="/path/to/roi.csv"

# Output root (host). Default is this folder: genescope_docker_bench/最终论文使用版
BENCH_ROOT="${SCRIPT_DIR}"

THREADS=16
MEM=64g
SEED=1
DATASET_ID="dataset_id"
ROI_ID="roi_id"

ROI_FILENAME="$(basename "${ROI_CSV}")"

if [[ "${OUTS}" == "/path/to/"* || "${ROI_CSV}" == "/path/to/"* ]]; then
  echo "[ERROR] Please edit OUTS/ROI_CSV paths before running." >&2
  exit 2
fi
if [[ ! -d "${OUTS}" ]]; then
  echo "[ERROR] OUTS directory not found: ${OUTS}" >&2
  exit 2
fi
if [[ ! -f "${ROI_CSV}" ]]; then
  echo "[ERROR] ROI CSV not found: ${ROI_CSV}" >&2
  exit 2
fi

LIMIT_ARGS=(--cpus "${THREADS}" --memory "${MEM}")
USER_ARGS=()
if [[ "$(uname -s)" == "Linux" ]]; then
  USER_ARGS=(--user "$(id -u)":"$(id -g)")
  LIMIT_ARGS+=(--memory-swap "${MEM}")
fi
ENV_ARGS=(
  -e "OMP_NUM_THREADS=${THREADS}"
  -e "OPENBLAS_NUM_THREADS=${THREADS}"
  -e "MKL_NUM_THREADS=${THREADS}"
  -e "NUMEXPR_NUM_THREADS=${THREADS}"
)

mkdir -p "${BENCH_ROOT}/${DATASET_ID}__${ROI_ID}"/{seagal,giotto,genescope,hotspot}

docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
  -v "${OUTS}:/data/outs:ro" \
  -v "${ROI_CSV}:/coord/${ROI_FILENAME}:ro" \
  -v "${BENCH_ROOT}/${DATASET_ID}__${ROI_ID}/seagal:/out" \
  seagal:final-rt \
  --data_dir /data/outs \
  --coord_file "/coord/${ROI_FILENAME}" \
  --outdir /out \
  --ncores "${THREADS}" \
  --seed "${SEED}" \
  --repeat 10 \
  --sample_sec 1 \
  --max_cells 0 \
  --grid_um 30 \
  --top_genes 0 \
  --extra_args_json '{"svg_topk":200,"use_pattern_genes":1,"min_counts":150,"min_cells":10,"n_permutation":99,"permute_ratio":0.2,"fdr_cutoff":0.05,"l_cutoff":0.1,"indep":1,"modules_enabled":1,"modules_nmax":30}'

docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
  -v "${OUTS}:/data/outs:ro" \
  -v "${ROI_CSV}:/coord/${ROI_FILENAME}:ro" \
  -v "${BENCH_ROOT}/${DATASET_ID}__${ROI_ID}/giotto:/out" \
  giotto:final10x-rt \
  --data_dir /data/outs \
  --coord_file "/coord/${ROI_FILENAME}" \
  --outdir /out \
  --max_cells 0 \
  --seed "${SEED}" \
  --threads "${THREADS}" \
  --repeat 10 \
  --dataset_id "${DATASET_ID}" \
  --roi_id "${ROI_ID}" \
  --roi_flip_y 1

docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
  -v "${OUTS}:/data/outs:ro" \
  -v "${ROI_CSV}:/coord/${ROI_FILENAME}:ro" \
  -v "${BENCH_ROOT}/${DATASET_ID}__${ROI_ID}/genescope:/out" \
  genescope:final-rt \
  --data_dir /data/outs \
  --coord_file "/coord/${ROI_FILENAME}" \
  --outdir /out \
  --seed "${SEED}" \
  --repeat 10 \
  --threads "${THREADS}" \
  --sample_sec 1 \
  --dataset_id "${DATASET_ID}" \
  --roi_id "${ROI_ID}"

docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
  -v "${OUTS}:/data/outs:ro" \
  -v "${ROI_CSV}:/coord/${ROI_FILENAME}:ro" \
  -v "${BENCH_ROOT}/${DATASET_ID}__${ROI_ID}/hotspot:/out" \
  hotspot:final-rt \
  --data_dir /data/outs \
  --coord_file "/coord/${ROI_FILENAME}" \
  --outdir /out \
  --threads "${THREADS}" \
  --repeat 10 \
  --fdr_autocorr 0.05 \
  --dataset_id "${DATASET_ID}" \
  --roi_id "${ROI_ID}" \
  --core_only 1

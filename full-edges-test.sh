#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FULL_DOCKER_DIR="${SCRIPT_DIR}/docker/Full-edges"
BENCH_ROOT="${BENCH_ROOT:-${SCRIPT_DIR}}"

# ---- Docker images ----
GENESCOPE_IMAGE="genescope:full-edges"
GIOTTO_IMAGE="giotto_grid:full-edges"
HOTSPOT_IMAGE="hotspot:full-edges"
SEAGAL_IMAGE="seagal:full-edges"

# ---- build full-edges images ----
docker build -t "${GENESCOPE_IMAGE}" "${FULL_DOCKER_DIR}/genescope"
docker build -t "${GIOTTO_IMAGE}"    "${FULL_DOCKER_DIR}/giotto_grid"
docker build -t "${HOTSPOT_IMAGE}"   "${FULL_DOCKER_DIR}/hotspot"
docker build -t "${SEAGAL_IMAGE}"    "${FULL_DOCKER_DIR}/seagal"

# ---- resources ----
THREADS="${THREADS:-96}"
MEM="${MEM:-}"
SEED="${SEED:-1}"

LIMIT_ARGS=(--cpus "${THREADS}")
if [[ -n "${MEM}" ]]; then
  LIMIT_ARGS+=(--memory "${MEM}")
fi

USER_ARGS=()
if [[ "$(uname -s)" == "Linux" ]]; then
  USER_ARGS=(--user "$(id -u)":"$(id -g)")
  if [[ -n "${MEM}" ]]; then
    LIMIT_ARGS+=(--memory-swap "${MEM}")
  fi
fi

ENV_ARGS=(
  -e "OMP_NUM_THREADS=${THREADS}"
  -e "OPENBLAS_NUM_THREADS=${THREADS}"
  -e "MKL_NUM_THREADS=${THREADS}"
  -e "NUMEXPR_NUM_THREADS=${THREADS}"
)

WORK_ROOT="${WORK_ROOT:-${BENCH_ROOT}/_work_full_edges}"
mkdir -p "${WORK_ROOT}"

mk_workdir() {
  local dataset="$1" method="$2" variant="$3"
  local d="${WORK_ROOT}/${dataset}/${method}_${variant}"
  if [[ -z "${WORK_ROOT}" || "${WORK_ROOT}" == "/" ]]; then
    echo "[ERROR] WORK_ROOT is unsafe: ${WORK_ROOT}" >&2
    exit 2
  fi
  rm -rf "${d}"
  mkdir -p "${d}"
  echo "${d}"
}

stage_two_files() {
  local src_repdir="$1" dst_repdir="$2"
  mkdir -p "${dst_repdir}"
  mv -f "${src_repdir}/edges_all.tsv" "${dst_repdir}/edges_all.tsv"
  mv -f "${src_repdir}/modules.tsv" "${dst_repdir}/modules.tsv"
}

run_dataset() {
  local dataset_name="$1"
  local outs="$2"
  local roi_csv="$3"
  local dataset_id="$4"
  local roi_id="$5"
  local roi_flip_y="$6"

  if [[ "${outs}" == "/path/to/"* || "${roi_csv}" == "/path/to/"* ]]; then
    echo "[ERROR] Please edit OUTS/ROI_CSV paths before running. dataset=${dataset_name}" >&2
    exit 2
  fi
  if [[ ! -d "${outs}" ]]; then
    echo "[ERROR] OUTS directory not found: ${outs}" >&2
    exit 2
  fi
  if [[ ! -f "${roi_csv}" ]]; then
    echo "[ERROR] ROI CSV not found: ${roi_csv}" >&2
    exit 2
  fi

  local roi_filename
  roi_filename="$(basename "${roi_csv}")"

  local bench_ds="${BENCH_ROOT}/${dataset_name}"
  mkdir -p "${bench_ds}"

  # ---- geneSCOPE (one run; stage repeat_001) ----
  {
    local work
    work="$(mk_workdir "${dataset_name}" "genescope" "run")"
    docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
      -v "${outs}:/data/outs:ro" \
      -v "${roi_csv}:/coord/${roi_filename}:ro" \
      -v "${work}:/out" \
      "${GENESCOPE_IMAGE}" \
      --data_dir /data/outs \
      --coord_file "/coord/${roi_filename}" \
      --outdir /out/repeat_001 \
      --grid_um 30 \
      --threads "${THREADS}" \
      --seed "${SEED}" \
      --sample_sec 0 \
      --parallel_backend serial \
      --dataset_id "${dataset_id}" \
      --roi_id "${roi_id}"

    stage_two_files "${work}/repeat_001" "${bench_ds}/genescope/repeat_001"
    rm -rf "${work}"
  }

  # ---- SEAGAL (one run; stage repeat_001) ----
  {
    local work
    work="$(mk_workdir "${dataset_name}" "seagal" "run")"
    docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
      -v "${outs}:/data/outs:ro" \
      -v "${roi_csv}:/coord/${roi_filename}:ro" \
      -v "${work}:/out" \
      "${SEAGAL_IMAGE}" \
      --data_dir /data/outs \
      --coord_file "/coord/${roi_filename}" \
      --outdir /out \
      --ncores "${THREADS}" \
      --seed "${SEED}" \
      --repeat 1 \
      --sample_sec 0 \
      --max_cells 0 \
      --grid_um 30 \
      --top_genes 0 \
      --modules_nmax 30 \
      --extra_args_json '{"seagal_dense_connectivities_max_obs":0,"svg_topk":1000,"use_pattern_genes":1,"min_counts":150,"min_cells":10,"n_permutation":99,"permute_ratio":0.2,"fdr_cutoff":0.05,"l_cutoff":0.1,"indep":1,"modules_enabled":1,"modules_nmax":30}' \
      --dataset_id "${dataset_id}" \
      --roi_id "${roi_id}"

    stage_two_files "${work}/repeat_001" "${bench_ds}/seagal/repeat_001"
    rm -rf "${work}"
  }

  # ---- Hotspot (two runs: modules + full edges) ----
  {
    local work_modules work_edges
    work_modules="$(mk_workdir "${dataset_name}" "hotspot" "modules")"
    work_edges="$(mk_workdir "${dataset_name}" "hotspot" "edges")"

    docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
      -v "${outs}:/data/outs:ro" \
      -v "${roi_csv}:/coord/${roi_filename}:ro" \
      -v "${work_modules}:/out" \
      "${HOTSPOT_IMAGE}" \
      --data_dir /data/outs \
      --coord_file "/coord/${roi_filename}" \
      --outdir /out \
      --threads "${THREADS}" \
      --seed "${SEED}" \
      --repeat 1 \
      --sample_sec 0 \
      --fdr_autocorr 0.05 \
      --core_only 1 \
      --emit_edges 0 \
      --emit_gg_matrices 0 \
      --write_module_scores 0 \
      --dataset_id "${dataset_id}" \
      --roi_id "${roi_id}"

    docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
      -v "${outs}:/data/outs:ro" \
      -v "${roi_csv}:/coord/${roi_filename}:ro" \
      -v "${work_edges}:/out" \
      "${HOTSPOT_IMAGE}" \
      --data_dir /data/outs \
      --coord_file "/coord/${roi_filename}" \
      --outdir /out \
      --threads "${THREADS}" \
      --seed "${SEED}" \
      --repeat 1 \
      --sample_sec 0 \
      --fdr_autocorr 1 \
      --core_only 0 \
      --emit_edges 1 \
      --emit_gg_matrices 0 \
      --write_module_scores 0 \
      --dataset_id "${dataset_id}" \
      --roi_id "${roi_id}"

    mkdir -p "${bench_ds}/hotspot/repeat_001"
    mv -f "${work_modules}/repeat_001/modules.tsv" "${bench_ds}/hotspot/repeat_001/modules.tsv"
    mv -f "${work_edges}/repeat_001/edges_all.tsv" "${bench_ds}/hotspot/repeat_001/edges_all.tsv"

    rm -rf "${work_modules}" "${work_edges}"
  }

  # ---- Giotto (two runs: modules + full edges) ----
  {
    local roi_flip_arg=()
    if [[ "${roi_flip_y}" == "1" ]]; then
      roi_flip_arg=(--roi_flip_y)
    fi

    local work_modules work_edges
    work_modules="$(mk_workdir "${dataset_name}" "giotto" "modules")"
    work_edges="$(mk_workdir "${dataset_name}" "giotto" "edges")"

    docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
      -v "${outs}:/data/outs:ro" \
      -v "${roi_csv}:/coord/${roi_filename}:ro" \
      -v "${work_modules}:/out" \
      "${GIOTTO_IMAGE}" \
      --data_dir /data/outs \
      --coord_file "/coord/${roi_filename}" \
      --outdir /out \
      "${roi_flip_arg[@]}" \
      --max_cells 0 \
      --seed "${SEED}" \
      --threads "${THREADS}" \
      --repeat 1 \
      --n_corr_genes 500 \
      --dataset_id "${dataset_id}" \
      --roi_id "${roi_id}"

    docker run --rm "${LIMIT_ARGS[@]}" "${USER_ARGS[@]}" "${ENV_ARGS[@]}" \
      -v "${outs}:/data/outs:ro" \
      -v "${roi_csv}:/coord/${roi_filename}:ro" \
      -v "${work_edges}:/out" \
      "${GIOTTO_IMAGE}" \
      --data_dir /data/outs \
      --coord_file "/coord/${roi_filename}" \
      --outdir /out \
      "${roi_flip_arg[@]}" \
      --max_cells 0 \
      --seed "${SEED}" \
      --threads "${THREADS}" \
      --repeat 1 \
      --n_corr_genes 5000 \
      --dataset_id "${dataset_id}" \
      --roi_id "${roi_id}"

    mkdir -p "${bench_ds}/giotto/repeat_001"
    mv -f "${work_modules}/repeat_001/modules.tsv" "${bench_ds}/giotto/repeat_001/modules.tsv"
    mv -f "${work_edges}/repeat_001/edges_all.tsv" "${bench_ds}/giotto/repeat_001/edges_all.tsv"

    rm -rf "${work_modules}" "${work_edges}"
  }
}

# ---- EDIT PATHS (host) ----
# Each dataset is written to:
#   ${BENCH_ROOT}/{P1,P2,P5,Lymph}/{genescope,giotto,hotspot,seagal}/repeat_001/{edges_all.tsv,modules.tsv}

LYMPH_OUTS="/path/to/Xenium_Lymph_node_outs"
LYMPH_ROI_CSV="/path/to/lymph_select2.csv"
run_dataset "Lymph" "${LYMPH_OUTS}" "${LYMPH_ROI_CSV}" "Xenium_Lymph_node" "Xenium_Lymph_node_ROI" 1

P1_OUTS="/path/to/GSE280314_Xenium_V1_Human_Colon_Cancer_P1_outs"
P1_ROI_CSV="/path/to/P1_tumor_region_coord.csv"
run_dataset "P1" "${P1_OUTS}" "${P1_ROI_CSV}" "Xenium_P1" "Xenium_P1" 0

P2_OUTS="/path/to/GSE280314_Xenium_V1_Human_Colon_Cancer_P2_outs"
P2_ROI_CSV="/path/to/P2_tumor_region_coord.csv"
run_dataset "P2" "${P2_OUTS}" "${P2_ROI_CSV}" "Xenium_P2" "Xenium_P2" 0

P5_OUTS="/path/to/GSE280314_Xenium_V1_Human_Colon_Cancer_P5_outs"
P5_ROI_CSV="/path/to/P5_tumor_region_coord.csv"
run_dataset "P5" "${P5_OUTS}" "${P5_ROI_CSV}" "Xenium_P5" "Xenium_P5" 0

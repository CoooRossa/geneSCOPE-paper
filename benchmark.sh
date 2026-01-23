#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKER_DIR="${SCRIPT_DIR}/docker"
RSCRIPTS_DIR="${SCRIPT_DIR}/benchmark-Rscripts"

OUT_ROOT="${SCRIPT_DIR}/final_out"
BENCH_ROOT_DEFAULT="${OUT_ROOT}/for_compare"
BENCH_ROOT="${BENCH_ROOT:-${BENCH_ROOT_DEFAULT}}" # mapping.R expects: <bench_root>/<method>/repeat_001/edges_all.tsv (staged from tool outputs all_edges.tsv)

OUTS="${OUTS:-/path/to/GSE280314_Xenium_V1_Human_Colon_Cancer_P5_CRC_Add_on_FFPE_outs}"
ROI_CSV="${ROI_CSV:-/path/to/P5_roi.csv}"
ROI_FILENAME="$(basename "${ROI_CSV}")"
ROI_IN_CONTAINER="/roi/${ROI_FILENAME}"

THREADS="${THREADS:-16}"
SEED="${SEED:-${seed:-1}}"
DATASET_ID="${DATASET_ID:-GSE280314_P5}"
ROI_ID="${ROI_ID:-P5_tumor_region}"
N_RUNS="${N_RUNS:-5}"

die() { echo "[ERROR] $*" >&2; exit 1; }
info() { echo "[INFO] $*" >&2; }

if [[ ! -d "${OUTS}" ]]; then
  die "OUTS directory not found: ${OUTS}"
fi
if [[ ! -f "${ROI_CSV}" ]]; then
  die "ROI CSV not found: ${ROI_CSV}"
fi

ROI_HOST_DIR="$(cd "$(dirname "${ROI_CSV}")" && pwd)"

if ! [[ "${N_RUNS}" =~ ^[0-9]+$ ]] || [[ "${N_RUNS}" -lt 1 ]]; then
  die "N_RUNS must be a positive integer (got: ${N_RUNS})"
fi

mkdir -p "${OUT_ROOT}" "${BENCH_ROOT}"

# ---- build docker images ----
docker build -t giotto:bench "${DOCKER_DIR}/giotto_grid"
docker build -t hotspot:bench "${DOCKER_DIR}/hotspot"
docker build -t seagal:bench "${DOCKER_DIR}/seagal"
docker build -t genescope:bench "${DOCKER_DIR}/genescope"

# ---- output dirs (non-alledges = module clustering; alledges = full edges) ----
GENESCOPE_OUT="${OUT_ROOT}/genescope"
GIOTTO_OUT="${OUT_ROOT}/giotto"
GIOTTO_ALLEDGES_OUT="${OUT_ROOT}/giotto-alledges"
HOTSPOT_OUT="${OUT_ROOT}/hotspot"
HOTSPOT_ALLEDGES_OUT="${OUT_ROOT}/hotspot-alledges"
SEAGAL_OUT="${OUT_ROOT}/seagal"
SEAGAL_ALLEDGES_OUT="${OUT_ROOT}/seagal-alledges"

mkdir -p "${HOTSPOT_OUT}" "${HOTSPOT_ALLEDGES_OUT}" "${GENESCOPE_OUT}" "${GIOTTO_OUT}" "${GIOTTO_ALLEDGES_OUT}" "${SEAGAL_OUT}" "${SEAGAL_ALLEDGES_OUT}"

# ---- Hotspot: non-alledges for modules five runs (modules + runtime) ----
for ((r=1; r<=N_RUNS; r++)); do
  REP_NAME="$(printf 'repeat_%03d' "${r}")"
  TMP_OUT="${HOTSPOT_OUT}/.tmp_${REP_NAME}"
  rm -rf "${TMP_OUT}"
  mkdir -p "${TMP_OUT}"

  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -e XDG_CACHE_HOME=/tmp/.cache \
    -v "${OUTS}:/data:ro" \
    -v "${ROI_HOST_DIR}:/roi:ro" \
    -v "${TMP_OUT}:/out" \
    hotspot:bench \
    --data_dir /data \
    --coord_file "${ROI_IN_CONTAINER}" \
    --outdir /out \
    --threads "${THREADS}" \
    --repeat 1 \
    --sample_sec 1 \
    --seed "${SEED}" \
    --dataset_id "${DATASET_ID}" \
    --roi_id "${ROI_ID}" \
    --core_only 1

  rm -rf "${HOTSPOT_OUT:?}/${REP_NAME}"
  mv "${TMP_OUT}/repeat_001" "${HOTSPOT_OUT}/${REP_NAME}"
  rm -rf "${TMP_OUT}"
done

# ---- Hotspot: alledges for full edges ----
docker run --rm -u "$(id -u):$(id -g)" \
  -v "${OUTS}:/data:ro" \
  -v "${ROI_HOST_DIR}:/roi:ro" \
  -v "${HOTSPOT_ALLEDGES_OUT}:/out" \
  hotspot:bench \
  --data_dir /data \
  --coord_file "${ROI_IN_CONTAINER}" \
  --outdir /out \
  --threads "${THREADS}" \
  --repeat 1 \
  --sample_sec 0 \
  --seed "${SEED}" \
  --dataset_id "${DATASET_ID}" \
  --roi_id "${ROI_ID}" \
  --fdr_autocorr 1 \
  --core_only 0

# ---- geneSCOPE: five runs (edges + modules + runtime) ----
for ((r=1; r<=N_RUNS; r++)); do
  REP_NAME="$(printf 'repeat_%03d' "${r}")"
  TMP_OUT="${GENESCOPE_OUT}/.tmp_${REP_NAME}"
  rm -rf "${TMP_OUT}"
  mkdir -p "${TMP_OUT}"

  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -e XDG_CACHE_HOME=/tmp/.cache \
    -v "${OUTS}:/data:ro" \
    -v "${ROI_HOST_DIR}:/roi:ro" \
    -v "${TMP_OUT}:/out" \
    genescope:bench \
    --data_dir /data \
    --coord_file "${ROI_IN_CONTAINER}" \
    --outdir /out \
    --seed "${SEED}" \
    --sample_sec 1 \
    --threads "${THREADS}" \
    --dataset_id "${DATASET_ID}" \
    --roi_id "${ROI_ID}"

  rm -rf "${GENESCOPE_OUT:?}/${REP_NAME}"
  mv "${TMP_OUT}/repeat_001" "${GENESCOPE_OUT}/${REP_NAME}"
  if [[ -f "${TMP_OUT}/stats.tsv" ]]; then
    mv "${TMP_OUT}/stats.tsv" "${GENESCOPE_OUT}/${REP_NAME}/stats.tsv"
  fi
  if [[ -f "${TMP_OUT}/run.log" ]]; then
    mv "${TMP_OUT}/run.log" "${GENESCOPE_OUT}/${REP_NAME}/run.log"
  fi
  rm -rf "${TMP_OUT}"

done

# ---- Giotto: alledges for full edges ----
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -e XDG_CACHE_HOME=/tmp/.cache \
  -v "${OUTS}:/data:ro" \
  -v "${ROI_HOST_DIR}:/roi:ro" \
  -v "${GIOTTO_ALLEDGES_OUT}:/out" \
  giotto:bench \
  --data_dir /data \
  --coord_file "${ROI_IN_CONTAINER}" \
  --outdir /out \
  --max_cells 0 \
  --repeat 1 \
  --grid_stepsize 30 \
  --threads "${THREADS}" \
  --seed "${SEED}" \
  --roi_flip_y --n_corr_genes 422

# ---- Giotto: non-alledges for modules five runs (modules + runtime)----
for ((r=1; r<=N_RUNS; r++)); do
  REP_NAME="$(printf 'repeat_%03d' "${r}")"
  TMP_OUT="${GIOTTO_OUT}/.tmp_${REP_NAME}"
  rm -rf "${TMP_OUT}"
  mkdir -p "${TMP_OUT}"

  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -e XDG_CACHE_HOME=/tmp/.cache \
    -v "${OUTS}:/data:ro" \
    -v "${ROI_HOST_DIR}:/roi:ro" \
    -v "${TMP_OUT}:/out:rw" \
    giotto:bench \
    --data_dir /data \
    --coord_file "${ROI_IN_CONTAINER}" \
    --outdir /out \
    --max_cells 0 \
    --repeat 1 \
    --grid_stepsize 30 \
    --threads "${THREADS}" \
    --seed "${SEED}" \
    --sample_sec 1 \
    --roi_flip_y --n_corr_genes 200

  rm -rf "${GIOTTO_OUT:?}/${REP_NAME}"
  mv "${TMP_OUT}/repeat_001" "${GIOTTO_OUT}/${REP_NAME}"
  if [[ -f "${TMP_OUT}/stats.tsv" ]]; then
    mv "${TMP_OUT}/stats.tsv" "${GIOTTO_OUT}/${REP_NAME}/stats.tsv"
  fi
  rm -rf "${TMP_OUT}"
done

# ---- SEAGAL: non-alledges for modules five runs (modules + runtime)----
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -e XDG_CACHE_HOME=/tmp/.cache \
  -v "${OUTS}:/data:ro" \
  -v "${ROI_HOST_DIR}:/roi:ro" \
  -v "${SEAGAL_OUT}:/out" \
  seagal:bench \
  --data_dir /data \
  --coord_file "${ROI_IN_CONTAINER}" \
  --outdir /out \
  --ncores "${THREADS}" \
  --seed "${SEED}" \
  --repeat "${N_RUNS}" \
  --sample_sec 1 \
  --max_cells 0 \
  --grid_um 30 \
  --top_genes 0 \
  --svg_topk 1000 --n_permutation 99 --permute_ratio 0.2 --fdr_cutoff 0.05 --l_cutoff 0.1 --indep 1 --modules_nmax 20

# ---- SEAGAL: alledges for full edges ----
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -e XDG_CACHE_HOME=/tmp/.cache \
  -v "${OUTS}:/data:ro" \
  -v "${ROI_HOST_DIR}:/roi:ro" \
  -v "${SEAGAL_ALLEDGES_OUT}:/out" \
  seagal:bench \
  --data_dir /data \
  --coord_file "${ROI_IN_CONTAINER}" \
  --outdir /out \
  --ncores "${THREADS}" \
  --seed "${SEED}" \
  --repeat 1 \
  --sample_sec 1 \
  --max_cells 0 \
  --grid_um 30 \
  --top_genes 0 \
  --svg_topk 1000 --n_permutation 99 --permute_ratio 0.2 --fdr_cutoff 1 --l_cutoff 0.1 --indep 1 --modules_nmax 20

# ---- Data format integration (stage outputs into mapping.R expected layout) ----
find_first_file() {
  local outdir="$1" fname="$2"
  if [[ -f "${outdir}/${fname}" ]]; then
    echo "${outdir}/${fname}"
    return 0
  fi
  if [[ -f "${outdir}/repeat_001/${fname}" ]]; then
    echo "${outdir}/repeat_001/${fname}"
    return 0
  fi
  local rep
  rep="$(ls -d "${outdir}"/repeat_* 2>/dev/null | sort | head -n 1 || true)"
  if [[ -n "${rep}" && -f "${rep}/${fname}" ]]; then
    echo "${rep}/${fname}"
    return 0
  fi
  return 1
}

validate_edges_all_tsv() {
  local path="$1"
  local header
  header="$(head -n 1 "${path}" | tr -d '\r')"
  [[ "${header}" == $'gene_a\tgene_b\tweight'* ]] || die "edges_all.tsv header must start with: gene_a<TAB>gene_b<TAB>weight (got: ${header}) path=${path}"
}

validate_modules_tsv() {
  local path="$1"
  local header
  header="$(head -n 1 "${path}" | tr -d '\r')"
  [[ "${header}" == $'gene\tmodule_id'* ]] || die "modules.tsv header must start with: gene<TAB>module_id (got: ${header}) path=${path}"
}

stage_method() {
  local method="$1" edges_outdir="$2" modules_outdir="$3"
  local repdir="${BENCH_ROOT}/${method}/repeat_001"
  mkdir -p "${repdir}"

  local edges_src
  edges_src="$(find_first_file "${edges_outdir}" "all_edges.tsv")" || die "Missing all_edges.tsv under: ${edges_outdir}"
  cp -f "${edges_src}" "${repdir}/edges_all.tsv"
  validate_edges_all_tsv "${repdir}/edges_all.tsv"

  local modules_src
  modules_src="$(find_first_file "${modules_outdir}" "modules.tsv")" || die "Missing modules.tsv under: ${modules_outdir}"
  cp -f "${modules_src}" "${repdir}/modules.tsv"
  validate_modules_tsv "${repdir}/modules.tsv"

  local meta_src=""
  if meta_src="$(find_first_file "${edges_outdir}" "meta.json" 2>/dev/null)"; then
    cp -f "${meta_src}" "${repdir}/meta.json"
  fi

  info "Staged ${method}: edges=$(basename "${edges_outdir}") modules=$(basename "${modules_outdir}") -> ${repdir}"
}

# geneSCOPE: edges+modules from one output dir
stage_method "genescope" "${GENESCOPE_OUT}" "${GENESCOPE_OUT}"

# other 3 tools: edges from alledges output dir; modules from non-alledges output dir
stage_method "giotto" "${GIOTTO_ALLEDGES_OUT}" "${GIOTTO_OUT}"
stage_method "hotspot" "${HOTSPOT_ALLEDGES_OUT}" "${HOTSPOT_OUT}"
stage_method "seagal" "${SEAGAL_ALLEDGES_OUT}" "${SEAGAL_OUT}"

# ---- Downstream: mapping + edge-level + module-level ----
Rscript "${RSCRIPTS_DIR}/mapping.R" \
  --bench_root "${BENCH_ROOT}" \
  --outdir "${OUT_ROOT}/string1step" \
  --methods "genescope,giotto,hotspot,seagal" \
  --keep_subscores 1

Rscript "${RSCRIPTS_DIR}/edge-level.R" \
  --map_dir "${OUT_ROOT}/string1step" \
  --outdir "${OUT_ROOT}/edge-level" \
  --methods "genescope,giotto,hotspot,seagal" \
  --pr_top_n_list "10,30,50,100,1000" \
  --edge_fdr_by_method "1,1,1,1" \
  --string_score_threshold 700

Rscript "${RSCRIPTS_DIR}/module-level.R" \
  --map_dir "${OUT_ROOT}/string1step" \
  --outdir "${OUT_ROOT}/module-level" \
  --methods "genescope,giotto,hotspot,seagal" \
  --modules_tsv_by_method "${BENCH_ROOT}/genescope/repeat_001/modules.tsv,${BENCH_ROOT}/giotto/repeat_001/modules.tsv,${BENCH_ROOT}/hotspot/repeat_001/modules.tsv,${BENCH_ROOT}/seagal/repeat_001/modules.tsv" \
  --min_module_genes 3 \
  --n_random 200000 \
  --seed 1 \
  --min_valid_null_draws 80000 \
  --max_resample_attempts 20 \
  --plot_max_modules_per_method 200 \
  --string_score_threshold 700

# ---- Runtime panels (collect runtime logs/stats -> runtime-panels.R) ----
RT_SRC_GENESCOPE="${OUT_ROOT}/genescope"
RT_SRC_GIOTTO="${OUT_ROOT}/giotto"
RT_SRC_HOTSPOT="${OUT_ROOT}/hotspot"
RT_SRC_SEAGAL="${OUT_ROOT}/seagal"

RUNTIME_RUN_ROOT="${OUT_ROOT}/runtime_panels_input"
RUNTIME_PLOTS_OUT="${RUNTIME_RUN_ROOT}/plots_runtime"

safe_rm_dir() {
  local d="$1"
  [[ -n "${d}" && "${d}" != "/" ]] || die "Refusing to remove unsafe dir: ${d}"
  [[ "${d}" == "${OUT_ROOT}/"* ]] || die "Refusing to remove outside OUT_ROOT: ${d}"
  rm -rf "${d}"
}

copy_repeat_stats_only() {
  local src_root="$1" dst_root="$2"
  [[ -d "${src_root}" ]] || die "Missing runtime directory: ${src_root}"
  mkdir -p "${dst_root}"

  local rep base num num_int dst_rep
  for rep in "${src_root}"/repeat_*; do
    [[ -d "${rep}" ]] || continue
    base="$(basename "${rep}")"
    num="${base#repeat_}"
    [[ "${num}" =~ ^[0-9]+$ ]] || continue
    num_int=$((10#${num}))
    dst_rep="${dst_root}/repeat_$(printf '%03d' "${num_int}")"
    mkdir -p "${dst_rep}"
    [[ -f "${rep}/stats.tsv" ]] || die "Missing stats.tsv: ${rep}/stats.tsv"
    cp -f "${rep}/stats.tsv" "${dst_rep}/stats.tsv"
  done
}

copy_seagal_runtime() {
  local src_root="$1" dst_root="$2"
  [[ -d "${src_root}" ]] || die "Missing runtime directory: ${src_root}"
  [[ -f "${src_root}/stats.tsv" ]] || die "Missing seagal stats.tsv: ${src_root}/stats.tsv"
  [[ -f "${src_root}/run.log" ]] || die "Missing seagal run.log: ${src_root}/run.log"

  mkdir -p "${dst_root}"
  cp -f "${src_root}/stats.tsv" "${dst_root}/stats.tsv"
  cp -f "${src_root}/run.log" "${dst_root}/run.log"

  local rep base num num_int dst_rep
  for rep in "${src_root}"/repeat_*; do
    [[ -d "${rep}" ]] || continue
    base="$(basename "${rep}")"
    num="${base#repeat_}"
    [[ "${num}" =~ ^[0-9]+$ ]] || continue
    num_int=$((10#${num}))
    dst_rep="${dst_root}/repeat_$(printf '%03d' "${num_int}")"
    mkdir -p "${dst_rep}"
    if [[ -f "${rep}/modules.tsv" ]]; then
      cp -f "${rep}/modules.tsv" "${dst_rep}/modules.tsv"
    elif [[ -f "${rep}/module.tsv" ]]; then
      cp -f "${rep}/module.tsv" "${dst_rep}/modules.tsv"
    else
      die "Missing seagal module file (modules.tsv or module.tsv): ${rep}"
    fi
  done
}

mkdir -p "${RUNTIME_RUN_ROOT}"
for sub in genescope giotto hotspot seagal plots_runtime; do
  if [[ -d "${RUNTIME_RUN_ROOT}/${sub}" ]]; then
    safe_rm_dir "${RUNTIME_RUN_ROOT}/${sub}"
  fi
done

copy_repeat_stats_only "${RT_SRC_GENESCOPE}" "${RUNTIME_RUN_ROOT}/genescope"
copy_repeat_stats_only "${RT_SRC_GIOTTO}" "${RUNTIME_RUN_ROOT}/giotto"
copy_repeat_stats_only "${RT_SRC_HOTSPOT}" "${RUNTIME_RUN_ROOT}/hotspot"
copy_seagal_runtime "${RT_SRC_SEAGAL}" "${RUNTIME_RUN_ROOT}/seagal"

Rscript "${RSCRIPTS_DIR}/runtime-panels.R" \
  --run_root "${RUNTIME_RUN_ROOT}" \
  --outdir "${RUNTIME_PLOTS_OUT}" \
  --n_runs "${N_RUNS}"

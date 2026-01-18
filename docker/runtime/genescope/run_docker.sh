#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: run_docker.sh --data_dir PATH --outdir DIR [options passed through]

This is a host-side Docker wrapper:
- Mounts --data_dir to /data (read-only)
- Mounts --outdir  to /out
- Translates --data_dir/--outdir/--coord_file to container paths (/data, /out, /data/...)
- Builds the image from this directory if missing

Options:
  --image_ref IMAGE   Docker image ref (or env IMAGE_REF)
USAGE
}

METHOD="genescope"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE_REF="${IMAGE_REF:-genescope/runtime-${METHOD}:latest}"

abs_path_any() {
  local p="$1"
  if [[ -d "${p}" ]]; then
    (cd "${p}" && pwd)
  else
    local d b
    d="$(cd "$(dirname "${p}")" && pwd)"
    b="$(basename "${p}")"
    printf '%s/%s\n' "${d}" "${b}"
  fi
}

ARGS=("$@")
DATA_DIR=""
OUTDIR=""
COORD_FILE=""
NCORES=""

for ((i=0; i<${#ARGS[@]}; i++)); do
  a="${ARGS[$i]}"
  case "${a}" in
    --data_dir) DATA_DIR="${ARGS[$((i+1))]}"; i=$((i+1));;
    --data_dir=*) DATA_DIR="${a#*=}";;
    --outdir) OUTDIR="${ARGS[$((i+1))]}"; i=$((i+1));;
    --outdir=*) OUTDIR="${a#*=}";;
    --coord_file) COORD_FILE="${ARGS[$((i+1))]}"; i=$((i+1));;
    --coord_file=*) COORD_FILE="${a#*=}";;
    --ncores|--threads) NCORES="${ARGS[$((i+1))]}"; i=$((i+1));;
    --ncores=*|--threads=*) NCORES="${a#*=}";;
    --image_ref) IMAGE_REF="${ARGS[$((i+1))]}"; i=$((i+1));;
    --image_ref=*) IMAGE_REF="${a#*=}";;
    -h|--help) usage; exit 0;;
  esac
done

if [[ -z "${DATA_DIR}" || -z "${OUTDIR}" ]]; then
  usage
  exit 2
fi

DATA_DIR_ABS="$(abs_path_any "${DATA_DIR}")"
mkdir -p "${OUTDIR}"
OUTDIR_ABS="$(cd "${OUTDIR}" && pwd)"

coord_in_container=""
if [[ -n "${COORD_FILE}" ]]; then
  if [[ ! -d "${DATA_DIR_ABS}" ]]; then
    echo "[ERROR] --coord_file requires --data_dir to be a directory (Xenium outs), not a file: ${DATA_DIR_ABS}" >&2
    exit 2
  fi
  COORD_ABS="$(abs_path_any "${COORD_FILE}")"
  if [[ "${COORD_ABS}" != "${DATA_DIR_ABS}/"* ]]; then
    echo "[ERROR] --coord_file must be under --data_dir so it is visible in-container at /data: ${COORD_ABS}" >&2
    exit 2
  fi
  coord_in_container="/data/${COORD_ABS#${DATA_DIR_ABS}/}"
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "[ERROR] docker not found in PATH" >&2
  exit 127
fi

if ! docker image inspect "${IMAGE_REF}" >/dev/null 2>&1; then
  echo "[INFO] Building image: ${IMAGE_REF} (context: ${SCRIPT_DIR})" >&2
  docker build -t "${IMAGE_REF}" "${SCRIPT_DIR}"
fi

USER_ARGS=()
if [[ "$(uname -s)" == "Linux" ]]; then
  USER_ARGS=(--user "$(id -u)":"$(id -g)")
fi

ENV_ARGS=()
if [[ -n "${NCORES}" && "${NCORES}" =~ ^[0-9]+$ && "${NCORES}" -ge 1 ]]; then
  ENV_ARGS+=(
    -e "OMP_NUM_THREADS=${NCORES}"
    -e "OPENBLAS_NUM_THREADS=${NCORES}"
    -e "MKL_NUM_THREADS=${NCORES}"
    -e "NUMEXPR_NUM_THREADS=${NCORES}"
  )
fi

CONTAINER_ARGS=()
i=0
while [[ "${i}" -lt "${#ARGS[@]}" ]]; do
  a="${ARGS[$i]}"
  case "${a}" in
    --data_dir) i=$((i+1)); CONTAINER_ARGS+=("--data_dir" "/data");;
    --data_dir=*) CONTAINER_ARGS+=("--data_dir=/data");;
    --outdir) i=$((i+1)); CONTAINER_ARGS+=("--outdir" "/out");;
    --outdir=*) CONTAINER_ARGS+=("--outdir=/out");;
    --coord_file) i=$((i+1)); CONTAINER_ARGS+=("--coord_file" "${coord_in_container}");;
    --coord_file=*) CONTAINER_ARGS+=("--coord_file=${coord_in_container}");;
    --image_ref) i=$((i+1));;
    --image_ref=*) ;;
    -h|--help) usage; exit 0;;
    *) CONTAINER_ARGS+=("${a}");;
  esac
  i=$((i+1))
done

exec docker run --rm \
  "${USER_ARGS[@]}" \
  -v "${DATA_DIR_ABS}:/data:ro" \
  -v "${OUTDIR_ABS}:/out" \
  "${ENV_ARGS[@]}" \
  "${IMAGE_REF}" \
  "${CONTAINER_ARGS[@]}"


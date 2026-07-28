#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPER_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${PAPER_ROOT}/correction-analysis/frozen_package.sh"
require_clean_paper_commit

SAMPLE_ID="$(printf '%s' "${1:-}" | tr '[:lower:]' '[:upper:]')"
case "${SAMPLE_ID}" in
  P5)
    FINAL_WORKFLOW_OUTPUT="${GENESCOPE_P5_OUTPUT:-${PAPER_ROOT}/correction-figures/P5}"
    WORKFLOW_SCRIPT="${SCRIPT_DIR}/P5_workflow.r"
    ;;
  LN)
    FINAL_WORKFLOW_OUTPUT="${GENESCOPE_LN_OUTPUT:-${PAPER_ROOT}/correction-figures/LN}"
    WORKFLOW_SCRIPT="${SCRIPT_DIR}/lymph.script.R"
    ;;
  *)
    echo "Usage: $0 P5|LN (set matching OUTS, SCOPE_RDS, and TOP_PAIRS variables)" >&2
    exit 2
    ;;
esac

output_parent="$(dirname "${FINAL_WORKFLOW_OUTPUT}")"
output_name="$(basename "${FINAL_WORKFLOW_OUTPUT}")"
if [[ -z "${output_name}" || "${output_name}" == "." || "${output_name}" == "/" ]]; then
  echo "Invalid figure output directory: ${FINAL_WORKFLOW_OUTPUT}" >&2
  exit 2
fi
mkdir -p "${output_parent}"
output_parent="$(cd "${output_parent}" && pwd)"
FINAL_WORKFLOW_OUTPUT="${output_parent}/${output_name}"
if [[ -e "${FINAL_WORKFLOW_OUTPUT}" || -L "${FINAL_WORKFLOW_OUTPUT}" ]]; then
  echo "Refusing an existing figure output path: ${FINAL_WORKFLOW_OUTPUT}" >&2
  exit 2
fi

OUTPUT_LOCK="${FINAL_WORKFLOW_OUTPUT}.publish-lock"
if ! mkdir "${OUTPUT_LOCK}" 2>/dev/null; then
  echo "Another frozen workflow owns the output lock: ${OUTPUT_LOCK}" >&2
  exit 2
fi

WORKFLOW_OUTPUT=""
cleanup_stage() {
  if [[ -n "${WORKFLOW_OUTPUT:-}" && -d "${WORKFLOW_OUTPUT}" ]]; then
    case "${WORKFLOW_OUTPUT}" in
      "${output_parent}/.${output_name}.staging."*) rm -rf -- "${WORKFLOW_OUTPUT}" ;;
      *) echo "Refusing to clean unexpected staging path: ${WORKFLOW_OUTPUT}" >&2 ;;
    esac
  fi
  if [[ -n "${OUTPUT_LOCK:-}" && -d "${OUTPUT_LOCK}" ]]; then
    rmdir "${OUTPUT_LOCK}" 2>/dev/null ||
      echo "Could not release frozen-workflow output lock: ${OUTPUT_LOCK}" >&2
  fi
}
trap cleanup_stage EXIT
WORKFLOW_OUTPUT="$(mktemp -d "${output_parent}/.${output_name}.staging.XXXXXX")"

if [[ "${SAMPLE_ID}" == "P5" ]]; then
  export GENESCOPE_P5_OUTPUT="${WORKFLOW_OUTPUT}"
else
  export GENESCOPE_LN_OUTPUT="${WORKFLOW_OUTPUT}"
fi

prepare_frozen_genescope "${WORKFLOW_OUTPUT}"
Rscript --vanilla "${PAPER_ROOT}/docker/genescope/freeze_assertions.R" 1.2.0
Rscript --vanilla "${WORKFLOW_SCRIPT}"
Rscript --vanilla "${PAPER_ROOT}/correction-analysis/verify_figure_bundle.R" \
  "${WORKFLOW_OUTPUT}" "${SAMPLE_ID}"
rm -rf -- "${WORKFLOW_OUTPUT}/.geneSCOPE-v1.2.0-library"
if [[ -e "${FINAL_WORKFLOW_OUTPUT}" || -L "${FINAL_WORKFLOW_OUTPUT}" ]]; then
  echo "Refusing to publish over a newly created output path: ${FINAL_WORKFLOW_OUTPUT}" >&2
  exit 2
fi
mv "${WORKFLOW_OUTPUT}" "${FINAL_WORKFLOW_OUTPUT}"
if [[ ! -d "${FINAL_WORKFLOW_OUTPUT}" ||
      ! -f "${FINAL_WORKFLOW_OUTPUT}/${SAMPLE_ID}_figure_manifest.json" ]]; then
  echo "Frozen figure output was not published as the expected directory." >&2
  exit 2
fi
WORKFLOW_OUTPUT=""
rmdir "${OUTPUT_LOCK}"
OUTPUT_LOCK=""
trap - EXIT
echo "Published frozen figure output: ${FINAL_WORKFLOW_OUTPUT}"

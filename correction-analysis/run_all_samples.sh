#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULT_ROOT="${RESULT_ROOT:-${SCRIPT_DIR}/results}"
THREADS="${THREADS:-8}"
source "${SCRIPT_DIR}/frozen_package.sh"
require_clean_paper_commit
prepare_frozen_genescope "${RESULT_ROOT}"
Rscript --vanilla "${SCRIPT_DIR}/../docker/genescope/freeze_assertions.R" 1.0.2

for sample_id in P1 P2 P5 LN; do
  input_var="GENESCOPE_${sample_id}_OUTS"
  input_dir="${!input_var:-}"
  if [[ -z "${input_dir}" ]]; then
    echo "Set ${input_var} to the ${sample_id} Xenium outs directory." >&2
    exit 2
  fi
  Rscript --vanilla "${SCRIPT_DIR}/run_v102_reanalysis.R" \
    "${sample_id}" "${input_dir}" "${RESULT_ROOT}" "${THREADS}"
done

if [[ -n "${VERIFY_REFERENCE_ROOT:-}" ]]; then
  Rscript --vanilla "${SCRIPT_DIR}/verify_reference_results.R" "${VERIFY_REFERENCE_ROOT}"
else
  echo "Four-sample rerun completed. Set VERIFY_REFERENCE_ROOT to a complete result bundle to run all freeze gates."
fi

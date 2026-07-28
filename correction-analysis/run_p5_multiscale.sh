#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
P5_OUTS="${GENESCOPE_P5_OUTS:-}"
RESULT_ROOT="${RESULT_ROOT:-${SCRIPT_DIR}/results/P5_multiscale}"
THREADS="${THREADS:-8}"
source "${SCRIPT_DIR}/frozen_package.sh"
require_clean_paper_commit
prepare_frozen_genescope "${RESULT_ROOT}"
Rscript --vanilla "${SCRIPT_DIR}/../docker/genescope/freeze_assertions.R" 1.2.0

if [[ -z "${P5_OUTS}" ]]; then
  echo "Set GENESCOPE_P5_OUTS to the P5 Xenium outs directory." >&2
  exit 2
fi

Rscript --vanilla "${SCRIPT_DIR}/run_p5_multiscale.R" "${P5_OUTS}" "${RESULT_ROOT}" "${THREADS}"

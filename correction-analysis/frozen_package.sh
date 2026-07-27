#!/usr/bin/env bash

FROZEN_PACKAGE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

require_clean_paper_commit() {
  local paper_root paper_commit
  paper_root="$(cd "${FROZEN_PACKAGE_SCRIPT_DIR}/.." && pwd)"
  if ! git -C "${paper_root}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "The frozen workflow requires git provenance for geneSCOPE-paper." >&2
    return 2
  fi
  if [[ -n "$(git -C "${paper_root}" status --porcelain --untracked-files=all)" ]]; then
    echo "Refusing to run from a dirty geneSCOPE-paper worktree." >&2
    return 2
  fi
  paper_commit="$(git -C "${paper_root}" rev-parse HEAD)"
  if [[ ! "${paper_commit}" =~ ^[0-9a-f]{40}$ ]]; then
    echo "Invalid geneSCOPE-paper commit: ${paper_commit}" >&2
    return 2
  fi
  export GENESCOPE_PAPER_COMMIT="${paper_commit}"
}

freeze_hash_tree() {
  local source_dir="$1"
  (
    cd "${source_dir}"
    find . -type f ! -name '.freeze-source-commit' ! -name '.freeze-tree-sha256' -print |
      LC_ALL=C sort |
      while IFS= read -r source_file; do
        if command -v sha256sum >/dev/null 2>&1; then
          sha256sum "${source_file}"
        elif command -v shasum >/dev/null 2>&1; then
          shasum -a 256 "${source_file}"
        else
          echo "sha256sum or shasum is required." >&2
          return 2
        fi
      done
  ) | if command -v sha256sum >/dev/null 2>&1; then
        sha256sum | awk '{print $1}'
      else
        shasum -a 256 | awk '{print $1}'
      fi
}

prepare_frozen_genescope() {
  local result_root="$1"
  local correction_dir paper_root vendor_dir source_commit expected_tree observed_tree
  local freeze_lib install_tmp install_source

  correction_dir="${FROZEN_PACKAGE_SCRIPT_DIR}"
  paper_root="$(cd "${correction_dir}/.." && pwd)"
  vendor_dir="${paper_root}/docker/genescope/vendor/geneSCOPE-v1.0.2"
  if [[ ! -f "${vendor_dir}/DESCRIPTION" ]]; then
    echo "Frozen geneSCOPE vendor is missing: ${vendor_dir}" >&2
    return 2
  fi

  source_commit="$(tr -d '[:space:]' < "${vendor_dir}/.freeze-source-commit")"
  expected_tree="$(tr -d '[:space:]' < "${vendor_dir}/.freeze-tree-sha256")"
  if [[ ! "${source_commit}" =~ ^[0-9a-f]{40}$ ]]; then
    echo "Invalid frozen source commit: ${source_commit}" >&2
    return 2
  fi
  if [[ ! "${expected_tree}" =~ ^[0-9a-f]{64}$ ]]; then
    echo "Invalid frozen vendor tree hash: ${expected_tree}" >&2
    return 2
  fi
  observed_tree="$(freeze_hash_tree "${vendor_dir}")"
  if [[ "${observed_tree}" != "${expected_tree}" ]]; then
    echo "Frozen vendor tree hash mismatch: ${observed_tree} != ${expected_tree}" >&2
    return 2
  fi

  mkdir -p "${result_root}"
  freeze_lib="${result_root}/.geneSCOPE-v1.0.2-library"
  mkdir -p "${freeze_lib}"
  # Reinstall on every invocation. R CMD INSTALL stages and replaces the
  # package directory, so a stale or locally modified cached install is never
  # trusted merely because a marker file is present.
  install_tmp="$(mktemp -d "${TMPDIR:-/tmp}/genescope-v102-install.XXXXXX")"
  install_source="${install_tmp}/geneSCOPE"
  mkdir -p "${install_source}"
  rsync -a "${vendor_dir}/" "${install_source}/"
  if ! R CMD INSTALL --preclean --library="${freeze_lib}" "${install_source}"; then
    rm -rf "${install_tmp}"
    return 2
  fi
  rm -rf "${install_tmp}"

  export R_LIBS_USER="${freeze_lib}"
  export GENESCOPE_SOURCE_COMMIT="${source_commit}"
  export GENESCOPE_VENDOR_TREE_SHA256="${expected_tree}"
}

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${1:-}"
TARGET_DIR="${SCRIPT_DIR}/vendor/geneSCOPE-v1.0.2"

if [[ -z "${SOURCE_DIR}" || ! -f "${SOURCE_DIR}/DESCRIPTION" ]]; then
  echo "Usage: $0 /absolute/path/to/geneSCOPE-v1.0.2" >&2
  exit 2
fi

SOURCE_DIR="$(cd "${SOURCE_DIR}" && pwd)"
SOURCE_VERSION="$(sed -n 's/^Version:[[:space:]]*//p' "${SOURCE_DIR}/DESCRIPTION")"
if [[ "${SOURCE_VERSION}" != "1.0.2" ]]; then
  echo "Refusing source with Version=${SOURCE_VERSION}; expected 1.0.2" >&2
  exit 2
fi
if ! grep -Fq 'return arma::dot(row_sums, row_sums);' "${SOURCE_DIR}/src/2.LeeL.cpp"; then
  echo "Refusing source without the canonical Lee S2 implementation" >&2
  exit 2
fi

if ! git -C "${SOURCE_DIR}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Refusing a source directory without git provenance." >&2
  exit 2
fi
if [[ -n "$(git -C "${SOURCE_DIR}" status --porcelain --untracked-files=all)" ]]; then
  echo "Refusing a source tree with uncommitted or untracked files." >&2
  exit 2
fi
SOURCE_COMMIT="$(git -C "${SOURCE_DIR}" rev-parse HEAD)"
if [[ ! "${SOURCE_COMMIT}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "Invalid source commit: ${SOURCE_COMMIT}" >&2
  exit 2
fi

# Materialize the recorded commit, rather than copying the working directory.
# This excludes ignored build products and editor metadata by construction.
ARCHIVE_DIR="$(mktemp -d)"
trap 'rm -rf "${ARCHIVE_DIR}"' EXIT
git -C "${SOURCE_DIR}" archive --format=tar "${SOURCE_COMMIT}" | tar -xf - -C "${ARCHIVE_DIR}"

mkdir -p "${TARGET_DIR}"
rsync -a --delete "${ARCHIVE_DIR}/" "${TARGET_DIR}/"

TREE_SHA256="$({
  cd "${TARGET_DIR}"
  find . -type f ! -name '.freeze-source-commit' ! -name '.freeze-tree-sha256' -print | LC_ALL=C sort |
    while IFS= read -r source_file; do
      if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "${source_file}"
      else
        shasum -a 256 "${source_file}"
      fi
    done
} | if command -v sha256sum >/dev/null 2>&1; then sha256sum; else shasum -a 256; fi | awk '{print $1}')"
printf '%s\n' "${SOURCE_COMMIT}" > "${TARGET_DIR}/.freeze-source-commit"
printf '%s\n' "${TREE_SHA256}" > "${TARGET_DIR}/.freeze-tree-sha256"

echo "Vendored geneSCOPE ${SOURCE_VERSION} at ${TARGET_DIR}"
echo "Source commit: ${SOURCE_COMMIT}"
echo "Vendor tree SHA-256: ${TREE_SHA256}"

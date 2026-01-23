#!/usr/bin/env bash
set -euo pipefail

IMAGE_REF="${IMAGE_REF:-seagal:full-edges}"

if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: docker not found in PATH" >&2
  exit 127
fi

EXTRA_ARGS="${DOCKER_BUILD_EXTRA_ARGS:-}"

# shellcheck disable=SC2086
exec docker build ${EXTRA_ARGS} -t "${IMAGE_REF}" .

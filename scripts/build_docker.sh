#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-build-docker}"
IMAGE="${DOLFINX_IMAGE:-dolfinx/dolfinx:stable}"

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "${ROOT_DIR}:/work" \
  -w /work \
  "${IMAGE}" \
  bash -lc "cmake -S . -B ${BUILD_DIR} && cmake --build ${BUILD_DIR} -j && ./${BUILD_DIR}/mie_scattering_2d"

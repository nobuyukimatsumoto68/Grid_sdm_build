#!/bin/bash
# build_dweofa_v5_gmult_claude.sh
#
# Build dweofa_mobius_HSDM_v5_gmult_claude (copy of v5; gauge-level MD multiplier from
# env GAUGE_MULT, default 4) via Makefile_claude, then stage the binary into ./bin/ where
# the flux submit scripts look for it. Output tee'd to build_dweofa_v5_gmult_claude.log.
#
# Run ON THE CLUSTER from anywhere (the script cds to its own directory, i.e. the
# src/gauge_gen_Nc4 it lives in):
#   source /usr/workspace/lsd/matsumoto5/su4_32c/env.sh
#   ./build_dweofa_v5_gmult_claude.sh

set -eu

SRC_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG="${SRC_DIR}/build_dweofa_v5_gmult_claude.log"
NAME="dweofa_mobius_HSDM_v5_gmult_claude"

cd "${SRC_DIR}"
: > "${LOG}"

echo "================ BUILD ${NAME} in ${SRC_DIR} ================" | tee -a "${LOG}"
make -f Makefile_claude "${NAME}" 2>&1 | tee -a "${LOG}"

if [ ! -x "${SRC_DIR}/${NAME}" ]; then
    echo "BUILD FAILED: ${SRC_DIR}/${NAME} not found" | tee -a "${LOG}"
    exit 1
fi

echo "================ STAGE INTO bin/ ================" | tee -a "${LOG}"
mkdir -p "${SRC_DIR}/bin"
cp "${SRC_DIR}/${NAME}" "${SRC_DIR}/bin/${NAME}"
ls -l "${SRC_DIR}/bin/${NAME}" | tee -a "${LOG}"

echo "================ DONE ================" | tee -a "${LOG}"

#!/bin/bash
# build_dweofa_v5_claude.sh
#
# Build the MX_inner=10000 test variant dweofa_mobius_HSDM_v5_claude (a copy of
# v4 with the single inner-CG iteration cap raised 1000 -> 10000). Builds in
# place via Makefile_claude, then copies the binary into ./bin/ where the flux
# run scripts look for it. All output tee'd to build_dweofa_v4_claude.log.
#
# Run from src/gauge_gen_Nc4/ ON THE CLUSTER:   ./build_dweofa_v4_claude.sh
#
# After it succeeds, point the run at the new binary:
#   in test/run_32c/flux.sh set   APP=../bin/dweofa_mobius_HSDM_v5_claude
#   then submit flux.sh as usual (same ip_hmc_mobius_test.xml, same grid/mpi).
# Compare per-trajectory wall time and the CG iteration counts in the log
# against the v4 baseline.

set -eu

SRC_DIR="/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4"
LOG="${SRC_DIR}/build_dweofa_v4_claude.log"
NAME="dweofa_mobius_HSDM_v5_claude"

cd "${SRC_DIR}"
: > "${LOG}"

echo "================ BUILD ${NAME} ================" | tee -a "${LOG}"
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
echo "Next: set APP=../bin/${NAME} in test/run_32c/flux.sh and submit." | tee -a "${LOG}"

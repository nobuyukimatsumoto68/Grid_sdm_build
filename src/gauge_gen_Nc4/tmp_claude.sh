#!/bin/bash
# tmp_claude.sh
#
# Chunk 1 handoff: build the momentum-projected connected meson driver
# (baryons_0000_dirac_claude.cc) and run a cold-config smoke test, then verify
# that the p=0 momentum dataset reproduces the original zero-momentum "_t"
# dataset. All output is tee'd to baryons_momproj_cold_claude.log.
#
# Run from src/gauge_gen_Nc4/:   ./tmp_claude.sh

set -eu

SRC_DIR="/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4"
ROOT="/mnt/baracuda_14/grid_claude/Grid_sdm_build"
LOG="${SRC_DIR}/baryons_momproj_cold_claude.log"
RUN_DIR="${SRC_DIR}/test/run_baryon_momproj_claude"

# fresh log
: > "${LOG}"

echo "================ BUILD ================" | tee -a "${LOG}"
# Reuse the proven compile path (grid-config from the local Grid build tree).
bash "${ROOT}/compile_two_baryon_claude.sh" baryons_0000_dirac_claude.cc 2>&1 | tee -a "${LOG}"

BIN="${SRC_DIR}/bin/baryons_0000_dirac_claude"
if [ ! -x "${BIN}" ]; then
    echo "BUILD FAILED: ${BIN} not found" | tee -a "${LOG}"
    exit 1
fi

echo "================ RUN (cold config) ================" | tee -a "${LOG}"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"
# No config args -> cold configuration; writes ColdConfig.h5 in this dir.
# Small grid (8.8.8.16): the p0 == _t verification is volume-independent, so a
# tiny lattice fits easily on one GPU and runs in seconds. Production 24.24.24.48
# stays on the cluster (multi-GPU).
"${BIN}" --grid 8.8.8.16 --mpi 1.1.1.1 --threads 8 --accelerator-threads 8 \
    --decomposition --comms-concurrent --comms-overlap --shm 2048 --shm-mpi 1 \
    2>&1 | tee -a "${LOG}"

H5="${RUN_DIR}/ColdConfig.h5"
if [ ! -f "${H5}" ]; then
    echo "RUN FAILED: ${H5} not produced" | tee -a "${LOG}"
    exit 1
fi

echo "================ VERIFY p0 == _t ================" | tee -a "${LOG}"
python3 "${SRC_DIR}/verify_momproj_p0_claude.py" "${H5}" 2>&1 | tee -a "${LOG}"

echo "================ DONE ================" | tee -a "${LOG}"

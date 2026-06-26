#!/bin/bash
# run_local_8c_v5_claude.sh
#
# LOCAL debug build+run of dweofa_mobius_HSDM_v5_claude on an 8^4 lattice using
# the LOCAL CUDA Grid build (NVIDIA TITAN V, sm_70) -- NOT the cluster hipcc
# Makefile. Builds with flags from build/grid-config, then runs a short 8^4
# ColdStart HMC for code debugging. All output tee'd to run_local_8c_v5_claude.log.
#
# Run from src/gauge_gen_Nc4/ on the LOCAL box:   ./run_local_8c_v5_claude.sh
#
# NOTE: 8^4 is a CORRECTNESS/debug check only -- it will NOT show the
# MX_inner=10000 speedup (that needs large volume + small mass).

set -eu

SRC_DIR="/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4"
GC="/mnt/baracuda_14/grid_claude/build/grid-config"
SRC="${SRC_DIR}/dweofa_mobius_HSDM_v5_claude.cc"
NAME="dweofa_mobius_HSDM_v5_claude"
RUN_DIR="${SRC_DIR}/test/run_8c_local_claude"
XML="${SRC_DIR}/ip_hmc_mobius_8c_claude.xml"
LOG="${SRC_DIR}/run_local_8c_v5_claude.log"
BIN="${RUN_DIR}/${NAME}"

# polite local resource cap
export OMP_NUM_THREADS=4

mkdir -p "${RUN_DIR}"
: > "${LOG}"

echo "================ BUILD (local CUDA via grid-config) ================" | tee -a "${LOG}"
# The install include tree (build/include) is incomplete (missing threads/), so
# compile against the SOURCE Grid at the project root, exactly as the local
# build/examples targets do: -I<builddir>/Grid (generated Config.h) and
# -I<srcdir>/Grid (the headers), with -DHAVE_CONFIG_H.
ROOT="/mnt/baracuda_14/grid_claude"
GRID_INC="-DHAVE_CONFIG_H -I${ROOT}/build/Grid -I${ROOT}/Grid"
CXX="$(${GC} --cxx)"
CXXFLAGS="$(${GC} --cxxflags)"
LDFLAGS="$(${GC} --ldflags)"
LIBS="$(${GC} --libs)"
echo "CXX: ${CXX}" | tee -a "${LOG}"
set -x
${CXX} ${GRID_INC} ${CXXFLAGS} "${SRC}" -o "${BIN}" ${LDFLAGS} ${LIBS} 2>&1 | tee -a "${LOG}"
set +x

if [ ! -x "${BIN}" ]; then
    echo "BUILD FAILED: ${BIN} not found" | tee -a "${LOG}"
    exit 1
fi

echo "================ RUN (8^4 ColdStart, 1 GPU) ================" | tee -a "${LOG}"
cd "${RUN_DIR}"
"${BIN}" --grid 8.8.8.8 --mpi 1.1.1.1 --threads 4 --accelerator-threads 8 \
    --shm 2048 --ParameterFile "${XML}" 2>&1 | tee -a "${LOG}"

echo "================ DONE ================" | tee -a "${LOG}"

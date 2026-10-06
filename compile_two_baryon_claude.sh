#!/bin/bash
# compile_two_baryon_claude.sh
#
# Compile a single gauge_gen_Nc4 source file against the local Grid build and
# place the binary in the source dir's ./bin, mirroring src/gauge_gen_Nc4/Makefile.
#
# Flags are taken from Grid's grid-config (cxx / cxxflags / ldflags / libs),
# exactly as the Makefile does. The only difference is GRID points at the Grid
# build tree that actually exists on this machine (../build), instead of the
# Makefile's ../../install/Grid_omp_Nc4.
#
# Usage:
#   ./compile_two_baryon_claude.sh                       # builds two_baryon_corr_claude.cc
#   ./compile_two_baryon_claude.sh some_other.cc         # builds another source
#   GRID=/path/to/grid/build ./compile_two_baryon_claude.sh   # override Grid build

set -euo pipefail

# Directory of this script = SDM repo root (Grid_sdm_build).
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Grid build tree (provides grid-config, include/, lib/libGrid.a).
GRID="${GRID:-${ROOT}/../build}"
# grid-config may live directly in the build dir or under bin/ (install layout).
if [ -f "${GRID}/bin/grid-config" ]; then
    CONFIG="${GRID}/bin/grid-config"
else
    CONFIG="${GRID}/grid-config"
fi

# Source file to compile (default: the two-baryon correlator driver).
SRC_NAME="${1:-two_baryon_corr_claude.cc}"
SRC_DIR="${ROOT}/src/gauge_gen_Nc4"
SRC="${SRC_DIR}/${SRC_NAME}"
BIN_DIR="${SRC_DIR}/bin"
BIN_NAME="$(basename "${SRC_NAME}" .cc)"

if [ ! -x "${CONFIG}" ] && [ ! -f "${CONFIG}" ]; then
    echo "ERROR: grid-config not found at ${CONFIG}" >&2
    echo "       Set GRID=/path/to/grid/build and retry." >&2
    exit 1
fi
if [ ! -f "${SRC}" ]; then
    echo "ERROR: source file not found: ${SRC}" >&2
    exit 1
fi

# Pull the same four variables the Makefile uses.
CXX="$(${CONFIG} --cxx)"
# The local Grid build was configured with /usr/local/cuda-12.6, which was removed
# (2026-09-29). If the configured compiler binary is gone, swap in the nvcc from
# NVCC_DIR (default /usr/local/cuda/bin, currently 12.9; same CUDA major as the
# libGrid.a build), keeping the remaining arguments (-std=c++17 -x cu).
CXX_BIN="${CXX%% *}"
if [ ! -x "${CXX_BIN}" ]; then
    NVCC_DIR="${NVCC_DIR:-/usr/local/cuda/bin}"
    echo "WARNING: configured compiler ${CXX_BIN} not found; using ${NVCC_DIR}/nvcc" >&2
    CXX="${NVCC_DIR}/nvcc ${CXX#* }"
fi
# nvcc 12.x accepts host gcc <= 14, but the system g++ is now 15 (OS upgrade). The
# -ccbin is Open MPI's mpic++, which calls ${OMPI_CXX:-g++}; point it at g++-12 unless
# the caller set OMPI_CXX.
GXX_MAJOR="$(g++ -dumpversion | cut -d. -f1)"
if [ -z "${OMPI_CXX:-}" ] && [ "${GXX_MAJOR}" -gt 14 ] && command -v g++-12 >/dev/null 2>&1; then
    echo "WARNING: system g++ is version ${GXX_MAJOR} (> 14, unsupported by nvcc 12); using OMPI_CXX=g++-12" >&2
    export OMPI_CXX=g++-12
fi
CXXFLAGS="$(${CONFIG} --cxxflags)"
LDFLAGS="$(${CONFIG} --ldflags)"
LIBS="$(${CONFIG} --libs)"

# This Grid build was not fully "make install"-ed, so its build/include header
# tree is incomplete (e.g. Grid/threads/Pragmas.h is missing). The complete
# headers live in the Grid source tree; add it as the first include path so
# all headers resolve there, while the generated Grid/Config.h is still picked
# up from the build dir's include path already in CXXFLAGS.
GRID_SRC="${GRID_SRC:-${ROOT}/../Grid}"
if [ -f "${GRID_SRC}/Grid/Grid.h" ]; then
    CXXFLAGS="-I${GRID_SRC} ${CXXFLAGS}"
else
    echo "WARNING: Grid source tree not found at ${GRID_SRC}; relying on ${GRID}/include only." >&2
fi

# GridStd.h does a bare #include "Config.h"; the generated Config.h lives in the
# build dir under Grid/, so add that directory directly to the include path.
if [ -f "${GRID}/Grid/Config.h" ]; then
    CXXFLAGS="-I${GRID}/Grid ${CXXFLAGS}"
fi

mkdir -p "${BIN_DIR}"

echo "Grid build : ${GRID}"
echo "Compiler   : ${CXX}"
echo "Source     : ${SRC}"
echo "Output     : ${BIN_DIR}/${BIN_NAME}"
echo "Compiling ..."

# Mirror Grid's native Makefile link rule (CXXLINK): the object comes BEFORE the
# libraries -- $(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(LIBS) -o $@. With a static
# libGrid.a the linker is single-pass left-to-right, so -lGrid / -lhdf5 must follow
# the source object or every Grid symbol comes up "undefined reference" (this is
# why the old LIBS-before-SRC order failed to link on dane's clang-14/GNU-ld even
# though it happened to work under tuolumne's more lenient hipcc/lld).
# CXX / *FLAGS are intentionally unquoted so they expand into multiple args.
${CXX} ${CXXFLAGS} ${LDFLAGS} "${SRC}" ${LIBS} -o "${BIN_DIR}/${BIN_NAME}"

echo "Done: ${BIN_DIR}/${BIN_NAME}"

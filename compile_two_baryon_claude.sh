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

# Mirror the Makefile rule:  $(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $< -o $@
# CXX / *FLAGS are intentionally unquoted so they expand into multiple args.
${CXX} ${CXXFLAGS} ${LDFLAGS} ${LIBS} "${SRC}" -o "${BIN_DIR}/${BIN_NAME}"

echo "Done: ${BIN_DIR}/${BIN_NAME}"

#!/bin/bash
# tmp_claude.sh
#
# Chunk 3: compile the contraction and run it on the (point-source) cold q00 file,
# producing the HDF5 output. Reports the single-baryon C_B print and lists the h5
# datasets so Claude can confirm the write.
#
# Run:  ./tmp_claude.sh
# Then Claude reads: chunk3_gevp_contract_claude.log
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${ROOT}/chunk3_gevp_contract_claude.log"
BINDIR="${ROOT}/src/gauge_gen_Nc4"

{
  echo "=================================================================="
  echo " (1) COMPILE two_baryon_gevp_contract_claude.cc"
  date
  echo "=================================================================="
  "${ROOT}/compile_two_baryon_claude.sh" two_baryon_gevp_contract_claude.cc

  echo
  echo "=================================================================="
  echo " (2) RUN contraction on ColdConfig.prop.lime (8^4) -> ColdConfig.gevp.h5"
  echo "=================================================================="
  ( cd "${BINDIR}" && ./bin/two_baryon_gevp_contract_claude --grid 8.8.8.8 --mpi 1.1.1.1 )

  echo
  echo "=================================================================="
  echo " (3) HDF5 dataset listing"
  echo "=================================================================="
  ( cd "${BINDIR}" && h5ls -r ColdConfig.gevp.h5 2>/dev/null | head -60 )

  echo
  echo "=================================================================="
  echo " chunk 3 run finished"
  echo "=================================================================="
} 2>&1 | tee "${LOG}"

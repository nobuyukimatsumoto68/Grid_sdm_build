#!/bin/bash
# tmp_claude.sh
#
# Momentum-projected bvar contraction, Chunk 3 (baryon_variants_mom_impl_plan_claude.md):
#   (1) compile baryon_variants_mom_corr_claude.cc
#   (2) COLD 8^4 (ColdConfig.pdfull.lime, 2 source corners): same args as the production
#       ColdConfig.bvar.h5 run (--pairs block --sink-smear 2.0 20) + all-pairs terms + --momops
#       --mom-raw --phase-check 3 -> ColdConfig.bvarmom.h5
#   (3) COLD with --no-p0 (point sink only)                    -> ColdConfig.bvarmom_nop0.h5
#   (4) REAL GAUGE 16^3x8 (lat16c.pdfull.lime): --phase-check 3, point + smeared sink
#       (3.0 40, --config ckpoint_lat.10000), timing of the 1225-pair pass -> lat16c.bvarmom.h5
#   (5) python checks: legacy vs production files, rest-frame relation, component equality,
#       hermiticity, --no-p0, effective energies
#
# NOTE: uses GPU device 0 -- run only when device 0 is free.
# Run:  ./tmp_claude.sh
# Then Claude reads: chunk_mom_validate_claude.log
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${ROOT}/chunk_mom_validate_claude.log"
BINDIR="${ROOT}/src/gauge_gen_Nc4"
CONF16="/mnt/baracuda_14/grid_claude/16c/ckpoint_lat.10000"
TERMS="baryon_variants_terms_all_claude.txt"
MOMOPS="baryon_mom_ops_claude.txt"
CHECK="baryon_variants_mom_check_claude.py"

# Open MPI rebuilt for the upgraded OS (build_openmpi_cuda129_claude.sh): the old
# openmpi_cuda hangs at launch. Same version 4.1.8, so LD_LIBRARY_PATH (which beats the
# binaries' RUNPATH) is enough to switch the runtime.
export PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/bin:${PATH}
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/lib:${LD_LIBRARY_PATH:-}

{
  echo "=================================================================="
  echo " (0) install the rebuilt Grid library + headers (Grid/ subdir only; the full"
  echo "     'make' stopped in examples/ on disc.cc / chcond3.cc: missing reweight.h)"
  date
  echo "=================================================================="
  echo "MPI runtime: $(which mpirun)"
  export OMPI_CXX=g++-12
  export OMPI_CC=gcc-12
  make -C /mnt/baracuda_14/grid_claude/build_cuda129/Grid install > /dev/null || { echo "GRID INSTALL FAILED"; exit 1; }
  ls -la /mnt/baracuda_14/grid_claude/build_cuda129/lib/libGrid.a || { echo "libGrid.a NOT INSTALLED"; exit 1; }

  echo
  echo "=================================================================="
  echo " (1) COMPILE baryon_variants_mom_corr_claude.cc"
  date
  echo "=================================================================="
  # ORIGINAL: "${ROOT}/compile_two_baryon_claude.sh" baryon_variants_mom_corr_claude.cc || { echo "COMPILE FAILED"; exit 1; }
  # Grid rebuilt with CUDA 12.9 + g++-12 (build_grid_cuda129_claude.sh); the old build/ is CUDA 12.6
  GRID=/mnt/baracuda_14/grid_claude/build_cuda129 "${ROOT}/compile_two_baryon_claude.sh" baryon_variants_mom_corr_claude.cc || { echo "COMPILE FAILED"; exit 1; }

  echo
  echo "=================================================================="
  echo " (2) COLD 8^4: --pairs block --sink-smear 2.0 20 --momops --mom-raw --phase-check 3 -> ColdConfig.bvarmom.h5"
  date
  echo "=================================================================="
  ( cd "${BINDIR}" && ./bin/baryon_variants_mom_corr_claude --grid 8.8.8.8 --mpi 1.1.1.1 \
      --src ColdConfig.pdfull.lime --out ColdConfig.bvarmom.h5 \
      --terms ${TERMS} --ops baryon_variants_ops_claude.txt --op all --pairs block \
      --sink-smear 2.0 20 --momops ${MOMOPS} --mom-raw --phase-check 3 ) || { echo "COLD RUN FAILED"; exit 1; }

  echo
  echo "=================================================================="
  echo " (3) COLD 8^4: --no-p0 (point sink only) -> ColdConfig.bvarmom_nop0.h5"
  date
  echo "=================================================================="
  ( cd "${BINDIR}" && ./bin/baryon_variants_mom_corr_claude --grid 8.8.8.8 --mpi 1.1.1.1 \
      --src ColdConfig.pdfull.lime --out ColdConfig.bvarmom_nop0.h5 \
      --terms ${TERMS} --momops ${MOMOPS} --no-p0 ) || { echo "COLD NO-P0 RUN FAILED"; exit 1; }

  echo
  echo "=================================================================="
  echo " (4) REAL GAUGE 16^3x8: --pairs block --sink-smear 3.0 40 --momops --phase-check 3 -> lat16c.bvarmom.h5"
  date
  echo "=================================================================="
  ( cd "${BINDIR}" && ./bin/baryon_variants_mom_corr_claude --grid 16.16.16.8 --mpi 1.1.1.1 \
      --src lat16c.pdfull.lime --out lat16c.bvarmom.h5 \
      --terms ${TERMS} --ops baryon_variants_ops_claude.txt --op all --pairs block \
      --config "${CONF16}" --sink-smear 3.0 40 --momops ${MOMOPS} --phase-check 3 ) || { echo "16c RUN FAILED"; exit 1; }
  date

  echo
  echo "=================================================================="
  echo " (5) python checks"
  echo "=================================================================="
  cd "${BINDIR}" || exit 1
  export OMP_NUM_THREADS=4
  echo "---- cold: legacy p=0 datasets vs production ColdConfig.bvar.h5"
  python3 ${CHECK} legacy ColdConfig.bvarmom.h5 ColdConfig.bvar.h5
  echo "---- cold: rest-frame relation Aaron ops vs our Cop"
  python3 ${CHECK} rest ColdConfig.bvarmom.h5
  echo "---- cold: component equality within (frame, irrep)"
  python3 ${CHECK} comps ColdConfig.bvarmom.h5
  echo "---- cold: hermiticity"
  python3 ${CHECK} herm ColdConfig.bvarmom.h5
  echo "---- cold: --no-p0 file vs full file"
  python3 ${CHECK} nop0 ColdConfig.bvarmom_nop0.h5 ColdConfig.bvarmom.h5
  echo "---- cold: effective energies"
  python3 ${CHECK} eff ColdConfig.bvarmom.h5
  echo "---- 16c: legacy p=0 datasets vs production lat16c.bvar.h5 (point sink only there)"
  python3 ${CHECK} legacy lat16c.bvarmom.h5 lat16c.bvar.h5
  echo "---- 16c: effective energies (one config, orientation only)"
  python3 ${CHECK} eff lat16c.bvarmom.h5
  echo "---- C++ built-in checks (grep)"
  grep -h -E "phase-check summary|momentum check:|wrote " "${LOG}" 2>/dev/null | tail -n 20

  echo
  echo "=================================================================="
  echo " chunk 3 (momentum) run finished"
  date
  echo "=================================================================="
} 2>&1 | tee "${LOG}"

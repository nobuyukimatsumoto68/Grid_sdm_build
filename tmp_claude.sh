#!/bin/bash
# tmp_claude.sh
#
# Chunk 6a: real-gauge test of the full-PD pipeline on a local SU(4) HMC checkpoint
# (16^3x8, /mnt/baracuda_14/grid_claude/16c/ckpoint_lat.10000, plaquette 0.58; study
# defaults M5 = 1.5, mass = 0.1 -- a code test, not physics):
#   (1) compile both binaries (header changed after the cleanup),
#   (2) dump: 16 pre-rotated solves at the origin, point source, --rottest (16 extra solves)
#       -> src/gauge_gen_Nc4/lat16c.pdfull.lime,
#   (3) contract: 9 operators, --pairs block, --naive-check 256 (sampled sites; the rough
#       gauge field exercises the full colour structure of the minors, unlike the cold config)
#       -> lat16c.bvar.h5,
#   (4) contract: all 35 states, --pairs all  -> lat16c.bvar_all.h5,
#   (5) python: hermiticity / positivity report and effective masses per channel.
#
# NOTE: uses GPU device 0.
# Run:  ./tmp_claude.sh
# Then Claude reads: chunk6a_16c_claude.log
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${ROOT}/chunk6a_16c_claude.log"
BINDIR="${ROOT}/src/gauge_gen_Nc4"
CONF="/mnt/baracuda_14/grid_claude/16c/ckpoint_lat.10000"
GRIDARGS="--grid 16.16.16.8 --mpi 1.1.1.1"

{
  echo "=================================================================="
  echo " (1) COMPILE point2all_prop_dumper_full_claude.cc and baryon_variants_corr_claude.cc"
  date
  echo "=================================================================="
  "${ROOT}/compile_two_baryon_claude.sh" point2all_prop_dumper_full_claude.cc || { echo "COMPILE FAILED (dumper)"; exit 1; }
  "${ROOT}/compile_two_baryon_claude.sh" baryon_variants_corr_claude.cc || { echo "COMPILE FAILED (driver)"; exit 1; }

  echo
  echo "=================================================================="
  echo " (2) DUMP ${CONF}: M5 1.5 mass 0.1, --width 0 --nsrc 1 --rottest -> lat16c.pdfull.lime"
  echo "=================================================================="
  ( cd "${BINDIR}" && ./bin/point2all_prop_dumper_full_claude "${CONF}" 1.5 0.1 lat16c.pdfull.lime ${GRIDARGS} --width 0 --nsrc 1 --rottest ) || { echo "DUMP FAILED"; exit 1; }
  ( cd "${BINDIR}" && ls -la lat16c.pdfull.lime && echo "binary records: $(grep -a -c ildg-binary-data lat16c.pdfull.lime)" )

  echo
  echo "=================================================================="
  echo " (3) CONTRACT 9 ops, --pairs block, --naive-check 256 -> lat16c.bvar.h5"
  echo "=================================================================="
  ( cd "${BINDIR}" && ./bin/baryon_variants_corr_claude ${GRIDARGS} \
      --src lat16c.pdfull.lime --out lat16c.bvar.h5 \
      --op all --pairs block --naive-check 256 ) || { echo "CONTRACT A FAILED"; exit 1; }

  echo
  echo "=================================================================="
  echo " (4) CONTRACT 35 states, --pairs all -> lat16c.bvar_all.h5"
  echo "=================================================================="
  ( cd "${BINDIR}" && ./bin/baryon_variants_corr_claude ${GRIDARGS} \
      --src lat16c.pdfull.lime --out lat16c.bvar_all.h5 \
      --ops baryon_variants_ops_all_claude.txt --op all --pairs all ) || { echo "CONTRACT B FAILED"; exit 1; }

  echo
  echo "=================================================================="
  echo " (5) python report"
  echo "=================================================================="
  ( cd "${BINDIR}" && OMP_NUM_THREADS=4 python3 baryon_variants_report_claude.py lat16c.bvar.h5 lat16c.bvar_all.h5 ) || { echo "REPORT FAILED"; exit 1; }

  echo
  echo "=================================================================="
  echo " chunk 6a run finished"
  echo "=================================================================="
} 2>&1 | tee "${LOG}"

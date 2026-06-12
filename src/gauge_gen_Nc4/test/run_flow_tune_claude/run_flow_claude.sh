#!/bin/bash
#FLUX: -t 240m
#FLUX: --output=flow_tune_{{id}}
#FLUX: -N 8
#FLUX: --exclusive
#
# Eigenvalue-flow scan for M5 tuning.
# For each thermalized config, run eye4_anti over a grid of M5 and collect the
# lowest eigenvalues of the Hermitian Wilson operator H_W(M5)^2 (= MdagM).
# Watch the lowest eigenvalue vs M5: pick M5 where the low spectrum is cleanest.
#
# eye4_anti CLI:  eye4_anti <config> <alpha> <beta> <M5>  [Grid flags]
#   alpha,beta = Chebyshev window for the Lanczos acceleration (TUNE THESE, see README).
#   Output is ALWAYS written to ./evals.dat (truncated each call) -> we cat-append it.

set -u
date

# -------- edit these for your site --------
GRID_SDM_ROOT=${GRID_SDM_ROOT:-/path/to/Grid_sdm_build}
ENV_SH=${ENV_SH:-${GRID_SDM_ROOT}/tioga/env.sh}
APP=${APP:-${GRID_SDM_ROOT}/src/gauge_gen_Nc4/bin/eye4_anti}

GRID_GEOM="32.32.32.64"
MPI_GEOM="2.2.2.4"        # = 32 ranks -> 8 nodes at tasks-per-node=4 (4D op, lighter than HMC)
NODES=8
TPN=4
THREADS=8

BETAS="10.7 10.8 10.9"
M5_LIST=$(seq 1.0 0.1 2.0)   # scan around the 1.5 ballpark
ALPHA=0.5                    # Chebyshev lo edge (just above the low cluster you want to resolve)
BETA_CHEBY=12.0              # Chebyshev hi edge (>= lambda_max of MdagM); see README
CONFIG_GLOB="ckpoint_lat.*"  # which saved configs per stream (default: all)
# ------------------------------------------

source "${ENV_SH}"
OPTIONS="--comms-concurrent --comms-overlap --shm 2048 --shm-mpi 1"

for BETA in ${BETAS}; do
  RUNDIR="b${BETA}"
  if [ ! -d "${RUNDIR}" ]; then
    echo "skip beta=${BETA}: no ${RUNDIR}/ (run HMC first)"
    continue
  fi

  for CFG in $(ls -1 ${RUNDIR}/${CONFIG_GLOB} 2>/dev/null); do
    CFG_ABS="$(pwd)/${CFG}"
    TAG=$(basename "${CFG}")
    OUTDIR="${RUNDIR}/flow_${TAG}"
    mkdir -p "${OUTDIR}"
    COLLECT="${OUTDIR}/flow_b${BETA}_${TAG}_claude.dat"
    : > "${COLLECT}"   # start fresh (truncate via redirect, no rm)

    echo "============================================================"
    echo "flow scan  beta=${BETA}  config=${CFG}"
    echo "============================================================"

    for M5 in ${M5_LIST}; do
      echo "--- M5=${M5} ---"
      PARAMS="--grid ${GRID_GEOM} --mpi ${MPI_GEOM} --threads ${THREADS} --accelerator-threads ${THREADS} ${OPTIONS}"
      # positional args first, Grid flags after (Grid_init strips its own)
      ( cd "${OUTDIR}" && \
        flux run -N ${NODES} --tasks-per-node=${TPN} --verbose --exclusive \
            --setopt=mpibind=verbose:1 \
            "${APP}" "${CFG_ABS}" "${ALPHA}" "${BETA_CHEBY}" "${M5}" ${PARAMS} ) \
        2>&1 | tee "${OUTDIR}/eye4_b${BETA}_${TAG}_M5_${M5}_claude.log"

      # eye4_anti always writes OUTDIR/evals.dat ("M5 <eigvals>"); accumulate it.
      if [ -f "${OUTDIR}/evals.dat" ]; then
        cat "${OUTDIR}/evals.dat" >> "${COLLECT}"
      fi
    done
    echo "collected flow table: ${COLLECT}"
  done
done

date
echo "ALL FLOW SCANS DONE"

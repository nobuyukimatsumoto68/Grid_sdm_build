#!/bin/bash
#FLUX: -t 240m
#FLUX: --output=flow_tune_{{id}}
#FLUX: -q pbatch
#FLUX: -N 8
#FLUX: -n 32
#FLUX: -g 1
#FLUX: --exclusive
#
# Eigenvalue-flow scan for M5 tuning (SU(4) SDM, new 32^3 x 64 tune ensemble).
# For each thermalized config in the HMC tune dirs, run eye4_anti over a grid of
# M5 and collect the lowest eigenvalues of the Hermitian Wilson operator
# H_W(M5)^2 (= MdagM). Pick M5 where the low spectrum is cleanest (largest, most
# stable smallest eigenvalue), which minimizes residual chiral breaking.
#
# eye4_anti CLI (verified against eye4_anti.cc):
#   eye4_anti <config> <alpha> <beta_cheby> <M5>  [Grid flags]
#   - builds WilsonFermionD with mass = -M5, antiperiodic time BC {1,1,1,-1};
#   - alpha,beta_cheby = ChebyParams.alpha/beta (Npoly=101 hardcoded) Lanczos window;
#   - ALWAYS writes ./evals.dat (one line "M5 [eigvals]", truncated each call) -> we
#     cd into a per-call dir and cat-append it into one table per config.
#
# Submit from this dir:  flux batch run_flow_claude.sh
# NOTE: the tune streams are still thermalizing (plaquette not yet plateaued), so
# treat early output as a wiring / qualitative check until the streams equilibrate.

set -u
date

# -------- site settings (su4_32c / tuolumne) --------
ENV_SH=/usr/workspace/lsd/matsumoto5/su4_32c/env.sh
APP=/usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build/src/gauge_gen_Nc4/bin/eye4_anti

GRID_GEOM="32.32.32.64"
MPI_GEOM="2.2.2.4"           # 32 ranks -> 8 nodes at tasks-per-node=4 (4D Wilson op)
NODES=8
TPN=4
THREADS=8

MASSES="0.01 0.05"
BETAS="10.7 10.8 10.9"
M5_LIST=$(seq 1.0 0.1 2.0)   # scan around the 1.5 ballpark
ALPHA=0.05                   # Chebyshev lo edge; calibrated via tune_cheby_claude.sh
BETA_CHEBY=12                # Chebyshev hi edge; 12 (or up to ~30) converges fast on the
                             # current configs (|H_W|min ~ 5.9, no near-zero modes). RE-CHECK
                             # the window at any M5 where the low spectrum drops toward zero.
CONFIG_GLOB="ckpoint_lat.*"  # which saved configs per stream (default: all)
# ----------------------------------------------------

export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off
source "${ENV_SH}"

OPTIONS="--comms-concurrent --comms-overlap --shm 2048 --shm-mpi 1"

for MASS in ${MASSES}; do
  for BETA in ${BETAS}; do
    STREAMDIR="/p/lustre5/matsumoto5/tune_32_64_m${MASS}/beta${BETA}m${MASS}"
    if [ ! -d "${STREAMDIR}" ]; then
      echo "skip m=${MASS} beta=${BETA}: no ${STREAMDIR} (run HMC first)"
      continue
    fi

    for CFG in $(ls -1 ${STREAMDIR}/${CONFIG_GLOB} 2>/dev/null | sort --version-sort); do
      TAG=$(basename "${CFG}")
      OUTDIR="${STREAMDIR}/flow_${TAG}"
      mkdir -p "${OUTDIR}"
      COLLECT="${OUTDIR}/flow_m${MASS}_b${BETA}_${TAG}_claude.dat"
      : > "${COLLECT}"   # start fresh (truncate via redirect, no rm)

      echo "============================================================"
      echo "flow scan  m=${MASS}  beta=${BETA}  config=${TAG}"
      echo "============================================================"

      for M5 in ${M5_LIST}; do
        echo "--- M5=${M5} ---"
        PARAMS="--grid ${GRID_GEOM} --mpi ${MPI_GEOM} --threads ${THREADS} --accelerator-threads ${THREADS} ${OPTIONS}"
        # eye4_anti writes ./evals.dat in cwd -> run inside OUTDIR. Positional args
        # first, Grid flags after (Grid_init strips its own).
        ( cd "${OUTDIR}" && \
          flux run -N ${NODES} --tasks-per-node=${TPN} --verbose --exclusive \
              --setopt=mpibind=verbose:1 \
              "${APP}" "${CFG}" "${ALPHA}" "${BETA_CHEBY}" "${M5}" ${PARAMS} ) \
          2>&1 | tee "${OUTDIR}/eye4_m${MASS}_b${BETA}_${TAG}_M5_${M5}_claude.log"

        # eye4_anti always writes OUTDIR/evals.dat ("M5 [eigvals]"); accumulate it.
        if [ -f "${OUTDIR}/evals.dat" ]; then
          cat "${OUTDIR}/evals.dat" >> "${COLLECT}"
        fi
      done
      echo "collected flow table: ${COLLECT}"
    done
  done
done

date
echo "ALL FLOW SCANS DONE"

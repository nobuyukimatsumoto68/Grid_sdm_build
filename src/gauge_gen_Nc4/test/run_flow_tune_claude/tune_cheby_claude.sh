#!/bin/bash
#FLUX: -t 120m
#FLUX: --output=cheby_tune_{{id}}
#FLUX: -q pbatch
#FLUX: -N 8
#FLUX: -n 32
#FLUX: -g 1
#FLUX: --exclusive
#
# Chebyshev-window calibration for eye4_anti (the M5-flow Lanczos), on the LATEST
# (still-thermalizing) checkpoint of ONE tune stream. Sweeps (alpha, beta_cheby) at
# fixed M5 and reports, per window, the converged-eigenvalue count Nconv and the
# lowest |H_W| eigenvalue. Goal: pick the window that resolves the low cluster with
# Nconv >= Nstop (=10) and a STABLE smallest eigenvalue, then plug alpha/beta_cheby
# into run_flow_claude.sh.
#
# Why a window at all: eye4_anti uses Chebyshev-accelerated Lanczos (Npoly=101 fixed).
# beta_cheby must be >~ lambda_max of MdagM (= H_W^2) so the high spectrum is damped;
# alpha sits just above the low cluster you want to resolve. Too-small beta_cheby
# under-damps and converges few/wrong modes; too-large slows convergence (Npoly fixed).
#
# eye4_anti CLI: eye4_anti <config> <alpha> <beta_cheby> <M5>  [Grid flags]
#   ALWAYS writes ./evals.dat (truncated each call) -> run each window in its own dir.
#   Diagnostics in stdout: Grid IRL "-- Nconv = N" / "#modes converged" / "NOT converged",
#   and eye4_anti's "i sqrt(lambda)= <|H_W eig|> ,l= <signed>" block (i=0 is the smallest).
#
# Submit from this dir:  flux batch tune_cheby_claude.sh

set -u
date

# -------- site settings (su4_32c / tuolumne) --------
ENV_SH=/usr/workspace/lsd/matsumoto5/su4_32c/env.sh
APP=/usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build/src/gauge_gen_Nc4/bin/eye4_anti

GRID_GEOM="32.32.32.64"
MPI_GEOM="2.2.2.4"
NODES=8
TPN=4
THREADS=8
# ----------------------------------------------------

# -------- what to calibrate on --------
MASS=0.01                    # stream to calibrate on (window is ~config-independent)
BETA=10.7
M5=1.5                       # fixed M5 (target ballpark) for the window sweep

# Chebyshev window grid to sweep (alpha = lo edge, beta_cheby = hi edge >~ lambda_max).
ALPHAS="0.05 0.1 0.5"
BETA_CHEBYS="12 30 60"
# --------------------------------------

STREAMDIR="/p/lustre5/matsumoto5/tune_32_64_m${MASS}/beta${BETA}m${MASS}"

export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off
source "${ENV_SH}"
OPTIONS="--comms-concurrent --comms-overlap --shm 2048 --shm-mpi 1"

# latest checkpoint in the stream dir (the "last checkpointer")
CFG=$(ls -1 ${STREAMDIR}/ckpoint_lat.* 2>/dev/null | sort --version-sort | tail -n1)
if [ -z "${CFG}" ]; then
  echo "ERROR: no ckpoint_lat.* in ${STREAMDIR}"
  exit 1
fi
TAG=$(basename "${CFG}")
echo "calibrating Chebyshev window on ${CFG}  (M5=${M5})"

TUNEROOT="${STREAMDIR}/cheby_tune_${TAG}_M5_${M5}"
mkdir -p "${TUNEROOT}"
SUMMARY="${TUNEROOT}/cheby_tune_summary_claude.dat"
: > "${SUMMARY}"
echo "# alpha  beta_cheby  Nconv  NOTconv  lowest_sqrtlambda  log" >> "${SUMMARY}"

PARAMS="--grid ${GRID_GEOM} --mpi ${MPI_GEOM} --threads ${THREADS} --accelerator-threads ${THREADS} ${OPTIONS}"

for A in ${ALPHAS}; do
  for B in ${BETA_CHEBYS}; do
    OUTDIR="${TUNEROOT}/a${A}_b${B}"
    mkdir -p "${OUTDIR}"
    LOG="${OUTDIR}/eye4_a${A}_b${B}_claude.log"

    echo "============================================================"
    echo "window  alpha=${A}  beta_cheby=${B}"
    echo "============================================================"

    ( cd "${OUTDIR}" && \
      flux run -N ${NODES} --tasks-per-node=${TPN} --verbose --exclusive \
          --setopt=mpibind=verbose:1 \
          "${APP}" "${CFG}" "${A}" "${B}" "${M5}" ${PARAMS} ) 2>&1 | tee "${LOG}"

    # definitive Nconv from Grid IRL ("... -- Nconv       = N"); lowest |H_W| eig
    # from eye4_anti's i=0 "sqrt(lambda)=" line.
    nconv=$(grep -- "-- Nconv" "${LOG}" | tail -n1 | sed 's/.*= *//' | tr -d ' ')
    [ -z "${nconv}" ] && nconv="?"
    notconv=$(grep -c "NOT converged" "${LOG}")
    low=$(grep "sqrt(lambda)=" "${LOG}" | head -n1 | sed 's/.*sqrt(lambda)= //; s/ ,l=.*//')
    [ -z "${low}" ] && low="NA"

    echo "${A}  ${B}  ${nconv}  ${notconv}  ${low}  ${LOG}" | tee -a "${SUMMARY}"
  done
done

echo "------------------------------------------------------------"
echo "summary table: ${SUMMARY}"
echo "pick the window with Nconv>=10, NOTconv=0, and a stable lowest_sqrtlambda;"
echo "put its alpha/beta_cheby into run_flow_claude.sh (ALPHA / BETA_CHEBY)."
date
echo "CHEBY TUNE DONE"

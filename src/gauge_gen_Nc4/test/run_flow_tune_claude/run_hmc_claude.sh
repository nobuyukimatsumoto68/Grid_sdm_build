#!/bin/bash
#FLUX: -t 360m
#FLUX: --output=hmc_tune_{{id}}
#FLUX: -N 16
#FLUX: --exclusive
#
# Thermalize a few configs at each gauge_beta for L5/M5 tuning.
# Streams: beta = 10.7, 10.8, 10.9 ; 32^3 x 64, m = 0.01, Ls = 16, M5 = 1.5.
# Configs land as <beta dir>/ckpoint_lat.<traj> (NERSC IEEE64BIG), read later by eye4_anti.
#
# Run sequentially inside ONE flux allocation. Size -N / MPI / -t before submitting.

set -u
date

# -------- edit these for your site --------
GRID_SDM_ROOT=${GRID_SDM_ROOT:-/path/to/Grid_sdm_build}      # repo root on the cluster
ENV_SH=${ENV_SH:-${GRID_SDM_ROOT}/tioga/env.sh}              # tioga/lassen env
APP=${APP:-${GRID_SDM_ROOT}/src/gauge_gen_Nc4/bin/dweofa_mobius_HSDM_v3}
XML_TEMPLATE=${XML_TEMPLATE:-$(pwd)/ip_hmc_mobius_claude.xml}

GRID_GEOM="32.32.32.64"
MPI_GEOM="2.2.2.8"        # = 64 ranks ; with tasks-per-node=4 -> 16 nodes. Keep prod(MPI) = N*4.
NODES=16
TPN=4                     # GPUs (tasks) per node
THREADS=8

BETAS="10.7 10.8 10.9"
# ------------------------------------------

export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off
source "${ENV_SH}"

OPTIONS="--comms-concurrent --comms-overlap --shm 2048 --shm-mpi 1"

for BETA in ${BETAS}; do
  RUNDIR="b${BETA}"
  mkdir -p "${RUNDIR}"
  XML="${RUNDIR}/ip_hmc_mobius_b${BETA}_claude.xml"
  sed "s/@BETA@/${BETA}/" "${XML_TEMPLATE}" > "${XML}"

  echo "============================================================"
  echo "HMC stream beta=${BETA}  dir=${RUNDIR}  xml=${XML}"
  echo "============================================================"

  # XML uses relative ./ckpoint_* prefixes, so launch from RUNDIR.
  PARAMS="--grid ${GRID_GEOM} --mpi ${MPI_GEOM} --threads ${THREADS} --accelerator-threads ${THREADS} ${OPTIONS} --ParameterFile $(basename ${XML})"

  ( cd "${RUNDIR}" && \
    flux run -N ${NODES} --tasks-per-node=${TPN} --verbose --exclusive \
        --setopt=mpibind=verbose:1 \
        "${APP}" ${PARAMS} ) 2>&1 | tee "${RUNDIR}/hmc_b${BETA}_claude.log"

  echo "stream beta=${BETA} done; saved configs:"
  ls -1 "${RUNDIR}"/ckpoint_lat.* 2>/dev/null
done

date
echo "ALL HMC STREAMS DONE"

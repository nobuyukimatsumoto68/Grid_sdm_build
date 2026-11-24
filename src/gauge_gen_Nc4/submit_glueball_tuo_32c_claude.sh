#!/bin/bash
# submit_glueball_tuo_32c_claude.sh
# 32^3 x 64 COPY of submit_glueball_tuo_claude.sh (kept unchanged); plan glueball_32c_impl_plan_claude.md.
# Changes: Lattice 32.32.32.64; GLUE_TRACE passed through (export GLUE_TRACE=1 -> per-path trace lines, Chunk 8.6).
# FLUX batch launcher for glueball_h5_claude (plan glueball_grid_impl_plan_claude.md, Chunk 8.4): NPARTS independent
# single-rank processes on one node (1 GPU each, no MPI comms), sharing the config range START..END step SKIP by
# --part p --nparts NPARTS. One glue.<conf>.h5 per config in OUTDIR (self-skip -> rerun = resume). Modeled on
# Grid/scripts/submit_wflow_q_h5_tuolumne_claude.sh.
# Env: CONF_PATH CONF_PREFIX START END SKIP OUTDIR [BIN] [NPARTS=4] [DEPTH=4] [ALPHA=0.33] [NLEV=4]
#
# Build first (tuolumne login node, after git pull):
#   source /usr/workspace/lsd/matsumoto5/su4_32c/env.sh
#   make -C /usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build/src/gauge_gen_Nc4 -f Makefile_claude glueball_h5_claude
# 32c VALIDATION (4 m0.4 configs, one per GPU, debug queue; Chunk C2):
#   CONF_PATH=/p/lustre5/matsumoto5/32_64_0.4/beta11.045m0.4 CONF_PREFIX=ckpoint_lat START=1940 END=2000 SKIP=20 \
#   OUTDIR=/p/lustre5/matsumoto5/obs_glue_3264_b11p045_m0p4000/bench_claude GLUE_TRACE=1 \
#   flux batch -q pdebug -t 60m --env=CONF_PATH --env=CONF_PREFIX --env=START --env=END --env=SKIP --env=OUTDIR \
#       --env=GLUE_TRACE submit_glueball_tuo_32c_claude.sh
# Claude does NOT submit -- run this yourself.
#FLUX: -t 240m
#FLUX: --output=glue_h5_{{id}}
#FLUX: -q pbatch
#FLUX: -N 1
#FLUX: --exclusive

set -u

# Lattice="24.24.24.48"
Lattice="32.32.32.64"
NPARTS=${NPARTS:-4}
DEPTH=${DEPTH:-4}
ALPHA=${ALPHA:-0.33}
NLEV=${NLEV:-4}
BIN=${BIN:-/usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build/src/gauge_gen_Nc4/glueball_h5_claude}

date
hostname
export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off
export HDF5_USE_FILE_LOCKING=FALSE
export GLUE_TRACE=${GLUE_TRACE:-0}
echo "GLUE_TRACE=${GLUE_TRACE}"
source /usr/workspace/lsd/matsumoto5/su4_32c/env.sh

if [ ! -x "${BIN}" ]; then
    echo "ERROR: binary not found: ${BIN} (make -C .../gauge_gen_Nc4 -f Makefile_claude glueball_h5_claude)" >&2
    exit 1
fi
mkdir -p "${OUTDIR}/logs"

OPTIONS="--decomposition --debug-mem --shm 2048"
PARAMS=" --grid ${Lattice} --mpi 1.1.1.1 --threads 8 --accelerator-threads 8 ${OPTIONS}"
APPPARAMS="--dir ${CONF_PATH} --prefix ${CONF_PREFIX} --start ${START} --end ${END} --step ${SKIP} --outdir ${OUTDIR} --depth ${DEPTH} --alpha ${ALPHA} --nlev ${NLEV}"

echo "--start " $(date) $(date +%s)
for((p=0;p<NPARTS;p++)); do
    LOG=${OUTDIR}/logs/glue_h5_${CONF_PREFIX}_${START}_${END}_p${p}_claude.log
    flux run -N 1 -n 1 -g 1 --verbose ${BIN} ${APPPARAMS} --part ${p} --nparts ${NPARTS} ${PARAMS} > ${LOG} 2>&1 &
done
wait
echo "--end " $(date) $(date +%s)
grep -h "# conf\|level 0\|APE x" ${OUTDIR}/logs/glue_h5_${CONF_PREFIX}_${START}_${END}_p*_claude.log

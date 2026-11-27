#!/bin/bash
# submit_glueball_stress_tuo_claude.sh
# Hang reproduction on tuolumne (plan glueball_grid_impl_plan_claude.md, Chunk 8.7): 4 single-GPU processes on one MI300A
# node, each running glueball_h5_claude --stress with one mode on the SAME 32c config (default m0.4 traj 1680 = the config
# whose production process hung at level 2 path 211 on 2026-11-25). Modes per GPU (default): full kernel slice split.
#   full  : pathTrace + momKernel + momSlice, cycling over the 241 paths (production work, markers K / S between halves)
#   kernel: phase kernel only;  slice: sliceSum of the 28-component object only;  split: 28 single-component sliceSums
# NCALLS is large on purpose; each process stops after MAXTIME (45 min); heartbeat every 1000 calls. A hang = a log whose
# last line stays "# st <i> <K|S|P|X|M>" (no "D") while the others keep going.
# Logs (job id in every name -> a resubmission never overwrites earlier logs):
#   ${OUTDIR}/stress_<jobid>_<mode>_lev<L>_conf<C>_claude.log   per process (markers, heartbeats)
#   ${OUTDIR}/stress_<jobid>_summary_claude.log                 job summary (settings, per-mode last heartbeat
#                                                                and last line), also in the flux output glue_stress_<id>
#
# Build first (tuolumne login node, after the user commits glueball_h5_claude.cc and this script, then git pull):
#   cd /usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build
#   source ../env.sh
#   make -C src/gauge_gen_Nc4 -f Makefile_claude glueball_h5_claude
# Submit (from src/gauge_gen_Nc4):
#   flux batch -q pdebug -t 60m submit_glueball_stress_tuo_claude.sh
# Options via env: CONF_PATH CONF_PREFIX CONF LEVEL MODES NCALLS MAXTIME OUTDIR
# Claude does NOT submit -- run this yourself.
#FLUX: -t 60m
#FLUX: --output=glue_stress_{{id}}
#FLUX: -q pdebug
#FLUX: -N 1
#FLUX: --exclusive

set -u

Lattice="32.32.32.64"
CONF_PATH=${CONF_PATH:-/p/lustre5/matsumoto5/32_64_0.4/beta11.045m0.4}
CONF_PREFIX=${CONF_PREFIX:-ckpoint_lat}
CONF=${CONF:-1680}
LEVEL=${LEVEL:-2}
MODES=${MODES:-"full kernel slice split"}
NCALLS=${NCALLS:-100000000}
# each process stops cleanly after MAXTIME s (so the summary below is written before the 60 min wall); a HUNG process
# never reaches it -> its log ends with "# st <i> <stage>" and the job runs to the wall.
MAXTIME=${MAXTIME:-2700}
OUTDIR=${OUTDIR:-/p/lustre5/matsumoto5/glue32_stress_claude}
BIN=${BIN:-/usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build/src/gauge_gen_Nc4/glueball_h5_claude}

date
hostname
export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off
source /usr/workspace/lsd/matsumoto5/su4_32c/env.sh

if [ ! -x "${BIN}" ]; then
    echo "ERROR: binary not found: ${BIN}" >&2
    exit 1
fi
mkdir -p "${OUTDIR}"
JOB=${FLUX_JOB_ID:-nojob$(date +%s)}
SUM=${OUTDIR}/stress_${JOB}_summary_claude.log
echo "stress job ${JOB} on $(hostname) $(date): conf ${CONF_PATH}/${CONF_PREFIX}.${CONF} level ${LEVEL} modes ${MODES} ncalls ${NCALLS} maxtime ${MAXTIME} s bin ${BIN} ($(ls -la --time-style=+%Y-%m-%d_%H:%M ${BIN} | awk '{print $6}'))" 2>&1 | tee ${SUM}

OPTIONS="--decomposition --debug-mem --shm 2048"
PARAMS=" --grid ${Lattice} --mpi 1.1.1.1 --threads 8 --accelerator-threads 8 ${OPTIONS}"

echo "--start " $(date) $(date +%s) 2>&1 | tee -a ${SUM}
for mode in ${MODES}; do
    LOG=${OUTDIR}/stress_${JOB}_${mode}_lev${LEVEL}_conf${CONF}_claude.log
    echo "mode ${mode} -> ${LOG}" 2>&1 | tee -a ${SUM}
    flux run -N 1 -n 1 -g 1 --verbose ${BIN} --stress ${mode} --dir ${CONF_PATH} --prefix ${CONF_PREFIX} --conf ${CONF} \
        --ncalls ${NCALLS} --maxtime ${MAXTIME} --level ${LEVEL} --depth 4 --alpha 0.33 ${PARAMS} > ${LOG} 2>&1 &
done
wait
echo "--end " $(date) $(date +%s) 2>&1 | tee -a ${SUM}
for mode in ${MODES}; do
    LOG=${OUTDIR}/stress_${JOB}_${mode}_lev${LEVEL}_conf${CONF}_claude.log
    echo "=== ${mode}" 2>&1 | tee -a ${SUM}
    grep "# stress heartbeat" ${LOG} | tail -1 2>&1 | tee -a ${SUM}
    grep "# stress WARNING" ${LOG} | head -3 2>&1 | tee -a ${SUM}
    tail -1 ${LOG} 2>&1 | tee -a ${SUM}
done

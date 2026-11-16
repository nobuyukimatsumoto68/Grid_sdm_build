#!/bin/bash
# run_bvar_dumper_heavy_point_flux_claude.sh
# HEAVY-MASS POINT-SOURCE copy of test/run_bvar_dumper_light/run_bvar_dumper_light_flux_claude.sh (itself a copy of the
# production test/run_bvar_dumper/run_bvar_dumper_flux_claude.sh, which exists only on tuolumne); both untouched.
# Plan: bvar_heavy_point_impl_plan_claude.md (repo root), user 2026-11-16. FLUX batch launcher (tuolumne / MI300A) for STAGE 1 of the baryon-operator-variants
# study: the full-PD point-to-all propagator dump point2all_prop_dumper_full_claude. Per config: the 16 PD
# source columns (4 PD spins x 4 colours, pre-rotated) of the origin corner are solved in 16/nrhs
# split-grid batches; the 16 PD blocks q<s><s'> are written as 16 single-prec SciDAC records per corner
# (~1.36 GB per config, atomic .inprogress -> rename).
#
# Differences from the light script (everything else is identical):
#   - ensembles: the heavy streams m0.1 (b10p865), m0.2 (b10p990), m0.3 (b11p035), m0.4 (b11p045)
#   - config sets = EXACTLY the existing heavy smeared-source dumps, so the two sources pair config by config
#     (user 2026-11-16): m0.1 1000..4640 step 40 (92), m0.2 1000..4680 step 40 (93), m0.3 1000..7960 step 80 (88),
#     m0.4 1000..9960 step 80 (113)
#   - SRC default = point (the smeared heavy dumps already exist); output .../bvar_2448_b<beta>_m<mass>_point/
#   - per-config wall budget PROP_TPT_SECONDS 150 s (as the production heavy dumps)
# NO deflation (user decision 2026-10-03, bvar_point_source_handoff_claude.md).
# Report CG iterations and wall time of the first config of each ensemble.
#
# ONE ENSEMBLE PER JOB (select with ENS). Self-skips configs whose final .pdfull.lime exists (never
# rm/overwrite). Graceful wall blocker (shell + in-binary PROP_DEADLINE_EPOCH). Claude does NOT submit.
#
# Build first (ON tuolumne; already built for the production dumps):
#   cd ${ROOT}/Grid_sdm_build && source ${ROOT}/env.sh
#   GRID=${ROOT}/build ./compile_two_baryon_claude.sh point2all_prop_dumper_full_claude.cc
# Submit yourself, VALIDATION first (one config), then the sets:
#   ENS=m0p1000 CLIST=1000 flux batch run_bvar_dumper_heavy_point_flux_claude.sh
#   ENS=m0p1000 flux batch run_bvar_dumper_heavy_point_flux_claude.sh
#   ENS=m0p2000 flux batch run_bvar_dumper_heavy_point_flux_claude.sh
#   ENS=m0p3000 flux batch run_bvar_dumper_heavy_point_flux_claude.sh
#   ENS=m0p4000 flux batch run_bvar_dumper_heavy_point_flux_claude.sh
# (self-skips finished configs: resubmit the same line to continue after the wall)
#FLUX: -t 480m
#FLUX: --output=bvar_dumper_heavy_point_{{id}}.out
#FLUX: -q pbatch
#FLUX: -N 8
#FLUX: -n 32
#FLUX: -g 1
#FLUX: --exclusive

set -u

date; hostname
export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off

ROOT=/usr/workspace/lsd/matsumoto5/su4_32c
source ${ROOT}/env.sh
APP=${ROOT}/Grid_sdm_build/src/gauge_gen_Nc4/bin/point2all_prop_dumper_full_claude

M5=1.5

# ---------------- ensemble (ONE per job): heavy streams only ----------------
ENS=${ENS:-m0p1000}
case "${ENS}" in
    *m0p1*|*10p865*)  betastr=10p865; massstr=0p1000; mass=0.1000; CMIN_DEF=1000; CMAX_DEF=4640; CSTRIDE_DEF=40 ;;
    *m0p2*|*10p990*)  betastr=10p990; massstr=0p2000; mass=0.2000; CMIN_DEF=1000; CMAX_DEF=4680; CSTRIDE_DEF=40 ;;
    *m0p3*|*11p035*)  betastr=11p035; massstr=0p3000; mass=0.3000; CMIN_DEF=1000; CMAX_DEF=7960; CSTRIDE_DEF=80 ;;
    *m0p4*|*11p045*)  betastr=11p045; massstr=0p4000; mass=0.4000; CMIN_DEF=1000; CMAX_DEF=9960; CSTRIDE_DEF=80 ;;
    *) echo "ERROR: ENS='${ENS}' not one of m0p1000/m0p2000/m0p3000/m0p4000 (light masses: run_bvar_dumper_light_flux_claude.sh)" >&2; exit 1 ;;
esac
cfgfilename=conf_nc4nf1_2448_b${betastr}_m${massstr}
CFGPATH=${CFGPATH:-/p/lustre5/matsumoto5/conf_nc4nf1_2448/${cfgfilename}}

# ---------------- source type: smeared (= q00 production w3/N40) or point ----------------
SRC=${SRC:-point}
case "${SRC}" in
    smeared) W=${W:-3.0}; N=${N:-40}; TAG="" ;;
    point)   W=0.0;       N=0;        TAG="_point" ;;
    *) echo "ERROR: SRC must be 'smeared' or 'point' (got '${SRC}')" >&2; exit 1 ;;
esac
NSRC=${NSRC:-1}               # corners: 1 = origin only (production)
ROTTEST=${ROTTEST:-0}         # 1 = append --rottest (16 extra unit solves/corner): validation only

OUTDIR=${OUTDIR:-/p/lustre5/matsumoto5/bvar_2448_b${betastr}_m${massstr}${TAG}}
mkdir -p "${OUTDIR}"

# ---------------- config set: the existing heavy smeared-source dump configs ----------------
CONF_MIN=${CONF_MIN:-${CMIN_DEF}}
CONF_MAX=${CONF_MAX:-${CMAX_DEF}}
CONF_STRIDE=${CONF_STRIDE:-${CSTRIDE_DEF}}
if [ -z "${CLIST:-}" ]; then
    CLIST=""
    i=${CONF_MIN}
    while [ "${i}" -le "${CONF_MAX}" ]; do CLIST="${CLIST} ${i}"; i=$(( i + CONF_STRIDE )); done
fi

# ---------------- decomposition (q00 production layout: nrhs = 32/4 = 8 -> 2 batches per corner) ----------------
LATT="24.24.24.48"
MPIGRID="2.2.2.4"
SPLIT="1 1 1 4"
OPTIONS="--decomposition --comms-concurrent --comms-overlap --debug-mem --shm 2048 --shm-mpi 1"
PARAMS_GRID="--grid ${LATT} --mpi ${MPIGRID} --split ${SPLIT} --width ${W} --niter ${N} --nsrc ${NSRC} --threads 8 --accelerator-threads 8 ${OPTIONS}"
[ "${ROTTEST}" = "1" ] && PARAMS_GRID="${PARAMS_GRID} --rottest"

if [ ! -x "${APP}" ]; then
    echo "ERROR: binary not found/executable: ${APP} -- build on tuolumne first (see header)" >&2
    exit 1
fi

echo "--start " "$(date)" "$(date +%s)"
echo "ensemble = ${cfgfilename}  (M5=${M5}, mass=${mass})   src=${SRC} W=${W} N=${N}  nsrc=${NSRC}  rottest=${ROTTEST}"
echo "cfgdir   = ${CFGPATH}"
echo "outdir   = ${OUTDIR}"
echo "configs  = ${CLIST}"

# ---------------- graceful wall blocker (shell + in-binary), as in the q00 dumper ----------------
# heavy masses: 150 s per-config budget (as the production heavy dumps; light used 600 s)
export PROP_TPT_SECONDS=${PROP_TPT_SECONDS:-150}
BLOCKER_OVERHEAD=${BLOCKER_OVERHEAD:-300}
TL0=$(flux job timeleft 2>/dev/null)
case "${TL0}" in ''|*[!0-9.]*) TL0=0 ;; esac
TL0=${TL0%.*}
if [ "${TL0}" -gt 0 ]; then
    export PROP_DEADLINE_EPOCH=$(( $(date +%s) + TL0 - BLOCKER_OVERHEAD ))
    echo "blocker: timeleft=${TL0}s -> PROP_DEADLINE_EPOCH=${PROP_DEADLINE_EPOCH} (tpt=${PROP_TPT_SECONDS}s)"
fi

for i in ${CLIST}; do
    cfg="${CFGPATH}/${cfgfilename}_lat.${i}"
    outfile="${OUTDIR}/${cfgfilename}.${i}.pdfull.lime"
    log="${OUTDIR}/${cfgfilename}.${i}.pdfull.log"

    if [ ! -f "${cfg}" ]; then echo "WARNING: config not found, skipping: ${cfg}"; continue; fi
    if [ -s "${outfile}" ]; then echo "SKIP (exists): ${outfile}"; continue; fi

    TIMELEFT=$(flux job timeleft 2>/dev/null)
    case "${TIMELEFT}" in ''|*[!0-9.]*) TIMELEFT=0 ;; esac
    TIMELEFT=${TIMELEFT%.*}
    if [ "${TIMELEFT}" -gt 0 ] && [ "${TIMELEFT}" -lt $(( PROP_TPT_SECONDS + BLOCKER_OVERHEAD )) ]; then
        echo "blocker: timeleft=${TIMELEFT}s -> stopping before config ${i} (resubmit to continue)"
        break
    fi

    echo "==== config ${i} -> ${outfile} ====" | tee "${log}"
    cstart=$(date +%s)
    # positional: <config> <M5> <mass> <outfile>
    flux run -N 8 --tasks-per-node=4 --verbose --exclusive \
        --setopt=mpibind=verbose:1 \
        "${APP}" "${cfg}" "${M5}" "${mass}" "${outfile}" ${PARAMS_GRID} \
        2>&1 | tee -a "${log}"
    cend=$(date +%s)
    echo "==== config ${i} done in $(( cend - cstart ))s ====" | tee -a "${log}"
done

echo "--end " "$(date)" "$(date +%s)"

#!/bin/bash
# run_bvar_dumper_light_flux_claude.sh
# LIGHT-MASS copy of test/run_bvar_dumper/run_bvar_dumper_flux_claude.sh (production, untouched; it exists
# only on tuolumne). FLUX batch launcher (tuolumne / MI300A) for STAGE 1 of the baryon-operator-variants
# study: the full-PD point-to-all propagator dump point2all_prop_dumper_full_claude. Per config: the 16 PD
# source columns (4 PD spins x 4 colours, pre-rotated) of the origin corner are solved in 16/nrhs
# split-grid batches; the 16 PD blocks q<s><s'> are written as 16 single-prec SciDAC records per corner
# (~1.36 GB per config, atomic .inprogress -> rename).
#
# Differences from the production script (everything else is identical):
#   - ensembles: only the light streams m0.05 (b10p840) and m0.01 (b10p800)
#   - config sets = the connected-meson configs after the thermalization cut chosen from the plaquette
#     history (plaq_light_claude/, 2026-10-07): m0.05 500..2032 step 4 (384), m0.01 160..756 step 4 (150)
#   - per-config wall budget PROP_TPT_SECONDS 600 s (light masses need more CG iterations; production 150 s)
# NO deflation (user decision 2026-10-03, bvar_point_source_handoff_claude.md).
# Report CG iterations and wall time of the first config of each ensemble.
#
# ONE ENSEMBLE PER JOB (select with ENS). Self-skips configs whose final .pdfull.lime exists (never
# rm/overwrite). Graceful wall blocker (shell + in-binary PROP_DEADLINE_EPOCH). Claude does NOT submit.
#
# Build first (ON tuolumne; already built for the production dumps):
#   cd ${ROOT}/Grid_sdm_build && source ${ROOT}/env.sh
#   GRID=${ROOT}/build ./compile_two_baryon_claude.sh point2all_prop_dumper_full_claude.cc
# Submit yourself (both source types per ensemble):
#   ENS=m0p0500 SRC=point   flux batch run_bvar_dumper_light_flux_claude.sh
#   ENS=m0p0500 SRC=smeared flux batch run_bvar_dumper_light_flux_claude.sh
#   ENS=m0p0100 SRC=point   flux batch run_bvar_dumper_light_flux_claude.sh
#   ENS=m0p0100 SRC=smeared flux batch run_bvar_dumper_light_flux_claude.sh
#FLUX: -t 480m
#FLUX: --output=bvar_dumper_light_{{id}}.out
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

# ---------------- ensemble (ONE per job): light streams only ----------------
ENS=${ENS:-m0p0500}
case "${ENS}" in
    *m0p05*|*10p840*) betastr=10p840; massstr=0p0500; mass=0.0500; CMIN_DEF=500; CMAX_DEF=2032; CSTRIDE_DEF=4 ;;
    *m0p01*|*10p800*) betastr=10p800; massstr=0p0100; mass=0.0100; CMIN_DEF=160; CMAX_DEF=756;  CSTRIDE_DEF=4 ;;
    *) echo "ERROR: ENS='${ENS}' not one of m0p0500/m0p0100 (heavy masses: use run_bvar_dumper_flux_claude.sh)" >&2; exit 1 ;;
esac
cfgfilename=conf_nc4nf1_2448_b${betastr}_m${massstr}
CFGPATH=${CFGPATH:-/p/lustre5/matsumoto5/conf_nc4nf1_2448/${cfgfilename}}

# ---------------- source type: smeared (= q00 production w3/N40) or point ----------------
SRC=${SRC:-smeared}
case "${SRC}" in
    smeared) W=${W:-3.0}; N=${N:-40}; TAG="" ;;
    point)   W=0.0;       N=0;        TAG="_point" ;;
    *) echo "ERROR: SRC must be 'smeared' or 'point' (got '${SRC}')" >&2; exit 1 ;;
esac
NSRC=${NSRC:-1}               # corners: 1 = origin only (production)
ROTTEST=${ROTTEST:-0}         # 1 = append --rottest (16 extra unit solves/corner): validation only

OUTDIR=${OUTDIR:-/p/lustre5/matsumoto5/bvar_2448_b${betastr}_m${massstr}${TAG}}
mkdir -p "${OUTDIR}"

# ---------------- config set: connected-meson configs after the thermalization cut ----------------
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
# light masses: 600 s per-config budget (production heavy masses: 150 s)
export PROP_TPT_SECONDS=${PROP_TPT_SECONDS:-600}
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

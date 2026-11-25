#!/bin/bash
# run_bvar_dumper_32c_flux_claude.sh
# 32^3 x 64 copy of test/run_bvar_dumper_heavy_point/run_bvar_dumper_heavy_point_flux_claude.sh and
# test/run_bvar_dumper_light/run_bvar_dumper_light_flux_claude.sh (both untouched). Plan: grid_claude/
# bvar_32c_port_impl_plan_claude.md (user go-ahead 2026-11-24). FLUX batch launcher (tuolumne / MI300A) for STAGE A of
# the baryon-operator-variants measurement: the full-PD point-to-all propagator dump point2all_prop_dumper_full_claude.
# Per config: the 16 PD source columns (4 PD spins x 4 colours, pre-rotated) of the origin corner, solved in split-grid
# batches; the 16 PD blocks q<s><s'> are written as 16 single-prec SciDAC records (~4.3 GB per config at 32c,
# 3.16 x the 24c 1.36 GB; atomic .inprogress -> rename inside the binary).
#
# Differences from the 24c scripts (everything else identical):
#   - lattice 32.32.32.64; MPI 2.2.2.4 on 8 nodes x 4 ranks (local 16^4 per rank), --split 1 1 1 4 (as 24c);
#   - all six 32c streams in ONE script: ENS = m0p0100 m0p0500 m0p1000 m0p2000 m0p3000 m0p4000;
#   - configs: the Grid HMC checkpoints /p/lustre5/matsumoto5/32_64_<m>/beta<beta>m<m>/ckpoint_lat.<traj> (NERSC format,
#     LoadNerscCheckpointer; read by NerscIO::readConfiguration); the stream directory is found by globbing
#     32_64_<m>/beta*m<m> (exactly one match required);
#   - DEFAULT CONFIG LIST = EVERY ckpoint_lat.<traj> present (user 2026-11-24: measure every config), CLIST overrides;
#   - output /p/lustre5/matsumoto5/bvar_3264_b<betastr>_m<massstr>[_point]/conf_nc4nf1_3264_b<betastr>_m<massstr>.<traj>.pdfull.lime;
#   - per-config wall budget PROP_TPT_SECONDS: 24c heavy 150 s / light 600 s scaled by ~3.2 -> 600 s heavy, 2400 s light
#     (to be replaced by the validation timings).
# Valence action = 32c sea action: Mobius Ls 16 (hardcoded in the binary), M5 1.5, b 1.5, c 0.5 (Grid/scripts/
# ip_hmc_mobius.xml). Source: SRC=smeared (w 3.0, N 40) or SRC=point. NO deflation (user 2026-10-03).
#
# ONE ENSEMBLE + ONE SOURCE TYPE PER JOB. Self-skips configs whose final .pdfull.lime exists (never rm/overwrite).
# Graceful wall blocker (shell + in-binary PROP_DEADLINE_EPOCH). Claude does NOT submit.
#
# Build (ON tuolumne; same binary as 24c, lattice size from --grid):
#   cd ${ROOT}/Grid_sdm_build && source ${ROOT}/env.sh
#   GRID=${ROOT}/build ./compile_two_baryon_claude.sh point2all_prop_dumper_full_claude.cc
# Submit yourself, VALIDATION first (validate_bvar_32c_flux_claude.sh), then per ensemble and source:
#   ENS=m0p4000 SRC=point   flux batch run_bvar_dumper_32c_flux_claude.sh
#   ENS=m0p4000 SRC=smeared flux batch run_bvar_dumper_32c_flux_claude.sh
# (self-skips finished configs: resubmit the same line to continue after the wall)
#FLUX: -t 480m
#FLUX: --output=bvar_dumper_32c_{{id}}.out
#FLUX: -q pbatch
#FLUX: -N 8
#FLUX: -n 32
#FLUX: -g 1
#FLUX: --exclusive

set -u

date
hostname
export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off

ROOT=/usr/workspace/lsd/matsumoto5/su4_32c
source ${ROOT}/env.sh
APP=${ROOT}/Grid_sdm_build/src/gauge_gen_Nc4/bin/point2all_prop_dumper_full_claude
LUSTRE=/p/lustre5/matsumoto5

M5=1.5

# ---------------- ensemble (ONE per job) ----------------
ENS=${ENS:-m0p4000}
case "${ENS}" in
    *m0p01*)
        m=0.01
        massstr=0p0100
        mass=0.0100
        TPT_DEF=2400
        ;;
    *m0p05*)
        m=0.05
        massstr=0p0500
        mass=0.0500
        TPT_DEF=2400
        ;;
    *m0p1*)
        m=0.1
        massstr=0p1000
        mass=0.1000
        TPT_DEF=600
        ;;
    *m0p2*)
        m=0.2
        massstr=0p2000
        mass=0.2000
        TPT_DEF=600
        ;;
    *m0p3*)
        m=0.3
        massstr=0p3000
        mass=0.3000
        TPT_DEF=600
        ;;
    *m0p4*)
        m=0.4
        massstr=0p4000
        mass=0.4000
        TPT_DEF=600
        ;;
    *)
        echo "ERROR: ENS='${ENS}' not one of m0p0100/m0p0500/m0p1000/m0p2000/m0p3000/m0p4000" >&2
        exit 1
        ;;
esac
# stream directory by glob (exactly one), beta from its name: beta<beta>m<m>
nmatch=$(ls -d ${LUSTRE}/32_64_${m}/beta*m${m} 2>/dev/null | wc -l)
if [ "${nmatch}" -ne 1 ]; then
    echo "ERROR: expected exactly one stream dir ${LUSTRE}/32_64_${m}/beta*m${m}, found ${nmatch}" >&2
    exit 1
fi
CFGDIR=${CFGDIR:-$(ls -d ${LUSTRE}/32_64_${m}/beta*m${m})}
beta=$(basename "${CFGDIR}" | sed -n "s/^beta\(.*\)m${m}$/\1/p")
betastr=$(printf "%.3f" "${beta}" | sed 's/\./p/')
cfgfilename=conf_nc4nf1_3264_b${betastr}_m${massstr}

# ---------------- source type ----------------
SRC=${SRC:-smeared}
case "${SRC}" in
    smeared)
        W=${W:-3.0}
        N=${N:-40}
        TAG=""
        ;;
    point)
        W=0.0
        N=0
        TAG="_point"
        ;;
    *)
        echo "ERROR: SRC must be 'smeared' or 'point' (got '${SRC}')" >&2
        exit 1
        ;;
esac
NSRC=${NSRC:-1}               # corners: 1 = origin only (production)
ROTTEST=${ROTTEST:-0}         # 1 = append --rottest (16 extra unit solves/corner): validation only

OUTDIR=${OUTDIR:-${LUSTRE}/bvar_3264_b${betastr}_m${massstr}${TAG}}
mkdir -p "${OUTDIR}"

# ---------------- config list: every ckpoint_lat.<traj> present (CLIST overrides) ----------------
if [ -z "${CLIST:-}" ]; then
    CLIST=$(ls "${CFGDIR}" | sed -n 's/^ckpoint_lat\.\([0-9][0-9]*\)$/\1/p' | sort -n | tr '\n' ' ')
fi

# ---------------- decomposition: 32.32.32.64 on 8 nodes x 4 ranks, local 16^4 ----------------
LATT="32.32.32.64"
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
echo "ensemble = ${cfgfilename}  (beta=${beta}, M5=${M5}, mass=${mass})   src=${SRC} W=${W} N=${N}  nsrc=${NSRC}  rottest=${ROTTEST}"
echo "cfgdir   = ${CFGDIR}"
echo "outdir   = ${OUTDIR}"
echo "configs  = ${CLIST}"

# ---------------- graceful wall blocker (shell + in-binary) ----------------
export PROP_TPT_SECONDS=${PROP_TPT_SECONDS:-${TPT_DEF}}
BLOCKER_OVERHEAD=${BLOCKER_OVERHEAD:-300}
TL0=$(flux job timeleft 2>/dev/null)
case "${TL0}" in
    ''|*[!0-9.]*)
        TL0=0
        ;;
esac
TL0=${TL0%.*}
if [ "${TL0}" -gt 0 ]; then
    export PROP_DEADLINE_EPOCH=$(( $(date +%s) + TL0 - BLOCKER_OVERHEAD ))
    echo "blocker: timeleft=${TL0}s -> PROP_DEADLINE_EPOCH=${PROP_DEADLINE_EPOCH} (tpt=${PROP_TPT_SECONDS}s)"
fi

for i in ${CLIST}; do
    cfg="${CFGDIR}/ckpoint_lat.${i}"
    outfile="${OUTDIR}/${cfgfilename}.${i}.pdfull.lime"
    log="${OUTDIR}/${cfgfilename}.${i}.pdfull.log"

    if [ ! -f "${cfg}" ]; then
        echo "WARNING: config not found, skipping: ${cfg}"
        continue
    fi
    if [ -s "${outfile}" ]; then
        echo "SKIP (exists): ${outfile}"
        continue
    fi

    TIMELEFT=$(flux job timeleft 2>/dev/null)
    case "${TIMELEFT}" in
        ''|*[!0-9.]*)
            TIMELEFT=0
            ;;
    esac
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

#!/bin/bash
# run_bvar_mom_contract_32c_flux_claude.sh
# 32^3 x 64 copy of test/run_bvar_mom_contract/run_bvar_mom_contract_flux_claude.sh (untouched). Plan: grid_claude/
# bvar_32c_port_impl_plan_claude.md (user go-ahead 2026-11-24). FLUX batch launcher (tuolumne / MI300A) for STAGE B:
# the momentum-projected single-baryon contraction baryon_variants_mom_corr_claude over the full-PD dumps of
# run_bvar_dumper_32c_flux_claude.sh. One config per invocation, SINGLE rank (--mpi 1.1.1.1, one GCD: the FFT host copy
# needs one rank), NSLOTS configs fanned out concurrently. Operators: A. S. Meyer's cubic / little-group irrep 4q ops
# (baryon_mom_ops_claude.txt), total momentum up to (2,1,1); all 1225 elementary pairs (baryon_variants_terms_all_claude.txt);
# smeared-sink copies (--sink-smear 3.0 40, MUST match the source smearing).
#
# Differences from the 24c script (everything else identical):
#   - lattice 32.32.32.64; all six 32c streams (ENS = m0p0100 ... m0p4000); stream dir by glob 32_64_<m>/beta*m<m>;
#   - dumps   /p/lustre5/matsumoto5/bvar_3264_b<betastr>_m<massstr>[_point]/conf_nc4nf1_3264_..<traj>.pdfull.lime
#   - gauge   <stream dir>/ckpoint_lat.<traj> (sink smearing)
#   - output  /p/lustre5/matsumoto5/obs_bvar_3264_b<betastr>_m<massstr>[_point]_ss_mom/bvarmom.all.<traj>.h5
#     (no 32c p = 0 production exists, so the p = 0 Cel_/Cop_ datasets are written too)
#   - MEMCHECK=1: /usr/bin/time -v around each config (max RSS in the per-config log; step-B memory check at 32c)
#   - wall 480 m (32c ~3.2 x the 24c ~170 s/config)
# SKIP if the output exists (never rm/overwrite); written as <out>.inprogress.h5, renamed on success.
# ONE ENSEMBLE + ONE SOURCE TYPE PER JOB. Claude does NOT submit.
#
# Build (ON tuolumne; same binary as 24c, lattice size from --grid):
#   cd ${ROOT}/Grid_sdm_build && source ${ROOT}/env.sh
#   GRID=${ROOT}/build ./compile_two_baryon_claude.sh baryon_variants_mom_corr_claude.cc
# Submit yourself AFTER the dumps of that ensemble/source exist (validation first: validate_bvar_32c_flux_claude.sh):
#   ENS=m0p4000 SRC=point   flux batch run_bvar_mom_contract_32c_flux_claude.sh
#   ENS=m0p4000 SRC=smeared flux batch run_bvar_mom_contract_32c_flux_claude.sh
#FLUX: -t 480m
#FLUX: --output=bvar_mom_contract_32c_{{id}}.out
#FLUX: -q pbatch
#FLUX: -N 2
#FLUX: -n 8
#FLUX: -g 1
#FLUX: --exclusive

set -u

date
hostname
export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off
export HDF5_USE_FILE_LOCKING=FALSE

ROOT=/usr/workspace/lsd/matsumoto5/su4_32c
source ${ROOT}/env.sh
SRCDIR=${ROOT}/Grid_sdm_build/src/gauge_gen_Nc4
APP=${SRCDIR}/bin/baryon_variants_mom_corr_claude
LUSTRE=/p/lustre5/matsumoto5

# runtime text tables (read by the binary)
TERMS=${TERMS:-${SRCDIR}/baryon_variants_terms_all_claude.txt}   # ALL 1225 pairs (momentum path needs them)
OPS=${OPS:-${SRCDIR}/baryon_variants_ops_claude.txt}             # legacy p=0 operators (as production)
MOMOPS=${MOMOPS:-${SRCDIR}/baryon_mom_ops_claude.txt}             # Aaron's irrep operators, 81 momenta

# ---------------- ensemble (ONE per job) ----------------
ENS=${ENS:-m0p4000}
case "${ENS}" in
    *m0p01*)
        m=0.01
        massstr=0p0100
        ;;
    *m0p05*)
        m=0.05
        massstr=0p0500
        ;;
    *m0p1*)
        m=0.1
        massstr=0p1000
        ;;
    *m0p2*)
        m=0.2
        massstr=0p2000
        ;;
    *m0p3*)
        m=0.3
        massstr=0p3000
        ;;
    *m0p4*)
        m=0.4
        massstr=0p4000
        ;;
    *)
        echo "ERROR: ENS='${ENS}' not one of m0p0100/m0p0500/m0p1000/m0p2000/m0p3000/m0p4000" >&2
        exit 1
        ;;
esac
nmatch=$(ls -d ${LUSTRE}/32_64_${m}/beta*m${m} 2>/dev/null | wc -l)
if [ "${nmatch}" -ne 1 ]; then
    echo "ERROR: expected exactly one stream dir ${LUSTRE}/32_64_${m}/beta*m${m}, found ${nmatch}" >&2
    exit 1
fi
CFGDIR=${CFGDIR:-$(ls -d ${LUSTRE}/32_64_${m}/beta*m${m})}
beta=$(basename "${CFGDIR}" | sed -n "s/^beta\(.*\)m${m}$/\1/p")
betastr=$(printf "%.3f" "${beta}" | sed 's/\./p/')
cfgfilename=conf_nc4nf1_3264_b${betastr}_m${massstr}

# ---------------- knobs ----------------
SRC=${SRC:-smeared}                 # which dump: point | smeared  (-> bvar_3264_<ens>[_point])
case "${SRC}" in
    smeared)
        TAG=""
        ;;
    point)
        TAG="_point"
        ;;
    *)
        echo "ERROR: SRC smeared|point" >&2
        exit 1
        ;;
esac
PAIRS=${PAIRS:-sector}              # legacy p=0 pair set, as production
SINKSMEAR=${SINKSMEAR:-1}           # 1 = add smeared-sink copies (--config gauge --sink-smear 3.0 40)
SS_W=${SS_W:-3.0}
SS_N=${SS_N:-40}                    # MUST match the source smearing (as production)
MOMRAW=${MOMRAW:-0}                 # 1 = also write per-component Cmomraw_*
PHASECHECK=${PHASECHECK:-0}         # >0 = --phase-check <n> on the first launched config
NSLOTS=${NSLOTS:-8}                 # concurrent single-GCD configs (= -n above)
MEMCHECK=${MEMCHECK:-0}             # 1 = /usr/bin/time -v (max RSS) around each config

DUMPDIR=${DUMPDIR:-${LUSTRE}/bvar_3264_b${betastr}_m${massstr}${TAG}}
SSTAG=""
if [ "${SINKSMEAR}" = "1" ]; then
    SSTAG="_ss"
fi
P0DIR=${P0DIR:-${LUSTRE}/obs_bvar_3264_b${betastr}_m${massstr}${TAG}${SSTAG}}     # no 32c p=0 production (p0 written)
OUTDIR=${OUTDIR:-${LUSTRE}/obs_bvar_3264_b${betastr}_m${massstr}${TAG}${SSTAG}_mom}
mkdir -p "${OUTDIR}"

LATT="32.32.32.64"
MPIGRID="1.1.1.1"                   # REQUIRED: the FFT host copy needs a single rank
PARAMS_GRID="--grid ${LATT} --mpi ${MPIGRID} --threads 8 --accelerator-threads 8 --shm 2048 --shm-mpi 1"

for f in "${APP}" "${TERMS}" "${OPS}" "${MOMOPS}"; do
    if [ ! -e "${f}" ]; then
        echo "ERROR: missing: ${f}" >&2
        exit 1
    fi
done
if [ ! -x "${APP}" ]; then
    echo "ERROR: binary not executable: ${APP} -- build on tuolumne first" >&2
    exit 1
fi
TIMECMD=""
if [ "${MEMCHECK}" = "1" ]; then
    TIMECMD="/usr/bin/time -v"
fi

echo "--start " "$(date)" "$(date +%s)"
echo "ensemble=${cfgfilename}  beta=${beta}  src=${SRC}  pairs=${PAIRS}  sinksmear=${SINKSMEAR} (w=${SS_W} N=${SS_N})  momraw=${MOMRAW}  phasecheck=${PHASECHECK}  NSLOTS=${NSLOTS}  memcheck=${MEMCHECK}"
echo "dumps  = ${DUMPDIR}"
echo "gauge  = ${CFGDIR}"
echo "outdir = ${OUTDIR}"

# config selection: CLIST override, else every dump present
if [ -n "${CLIST:-}" ]; then
    dumps=""
    for i in ${CLIST}; do
        dumps="${dumps} ${DUMPDIR}/${cfgfilename}.${i}.pdfull.lime"
    done
else
    dumps=$(ls "${DUMPDIR}"/${cfgfilename}.*.pdfull.lime 2>/dev/null)
fi
if [ -z "${dumps}" ]; then
    echo "no dumps found in ${DUMPDIR} (run stage A first)"
    exit 0
fi

running=0
first=1
for dump in ${dumps}; do
    base=$(basename "${dump}")
    stem=${base%.pdfull.lime}
    idx=${stem##*.}
    gauge="${CFGDIR}/ckpoint_lat.${idx}"
    out="${OUTDIR}/bvarmom.all.${idx}.h5"
    tmp="${OUTDIR}/bvarmom.all.${idx}.h5.inprogress.h5"
    log="${OUTDIR}/bvarmom.all.${idx}.log"
    p0="${P0DIR}/bvar.all.${idx}.h5"

    if [ ! -f "${dump}" ]; then
        echo "  WARN no dump, skip idx ${idx}: ${dump}"
        continue
    fi
    if [ -s "${out}" ]; then
        echo "  SKIP (exists): ${out}"
        continue
    fi

    extra=""
    if [ "${SINKSMEAR}" = "1" ]; then
        if [ ! -f "${gauge}" ]; then
            echo "  WARN no gauge for sink smear, skip idx ${idx}: ${gauge}"
            continue
        fi
        extra="--config ${gauge} --sink-smear ${SS_W} ${SS_N}"
    fi
    p0mode="p0 written"
    if [ -s "${p0}" ]; then
        extra="${extra} --no-p0"
        p0mode="--no-p0 (production ${p0} exists)"
    fi
    if [ "${MOMRAW}" = "1" ]; then
        extra="${extra} --mom-raw"
    fi
    if [ "${PHASECHECK}" -gt 0 ] && [ "${first}" = "1" ]; then
        extra="${extra} --phase-check ${PHASECHECK}"
    fi
    first=0

    echo "  -> config ${idx}  [${p0mode}]"
    # one task, one GCD; rename the finished file only if the binary exited cleanly
    ( flux run -n1 -g1 \
        ${TIMECMD} "${APP}" --src "${dump}" --out "${tmp}" \
        --terms "${TERMS}" --ops "${OPS}" --op all --pairs "${PAIRS}" \
        --momops "${MOMOPS}" ${extra} \
        ${PARAMS_GRID} > "${log}" 2>&1 \
      && mv "${tmp}" "${out}" \
      || echo "  FAILED config ${idx} (see ${log})" ) &

    running=$(( running + 1 ))
    if [ "${running}" -ge "${NSLOTS}" ]; then
        wait -n
        running=$(( running - 1 ))
    fi
done
wait
echo "--end " "$(date)" "$(date +%s)"

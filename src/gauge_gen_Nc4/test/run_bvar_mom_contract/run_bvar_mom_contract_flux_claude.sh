#!/bin/bash
# run_bvar_mom_contract_flux_claude.sh
# FLUX batch launcher (tuolumne / MI300A) for the MOMENTUM-PROJECTED single-baryon contraction
# baryon_variants_mom_corr_claude over the full-PD dumps (Stage 1 of the baryon-variants study).
# Mirrors test/run_bvar_contract/run_bvar_contract_flux_claude.sh (production, untouched): one config
# per invocation, single rank (--mpi 1.1.1.1, one GCD), NSLOTS configs fanned out concurrently.
# Operators: A. S. Meyer's cubic / little-group irrep 4q ops (baryon_mom_ops_claude.txt), total
# momentum up to (2,1,1); all 1225 elementary pairs (baryon_variants_terms_all_claude.txt).
# See bvar_mom_contract_handoff_claude.md and baryon_variants_mom_impl_plan_claude.md.
#
# Per config:
#   - output  ${LUSTRE}/obs_bvar_2448_<ens>[_point][_ss]_mom/bvarmom.all.<idx>.h5
#     SKIP if it exists (never rm/overwrite); written as <out>.inprogress.h5, renamed on success
#   - if the production file ${LUSTRE}/obs_bvar_2448_<ens>[_point][_ss]/bvar.all.<idx>.h5 exists,
#     --no-p0 (only the momentum datasets); otherwise the p=0 Cel_/Cop_ datasets are written too
#   - PHASECHECK>0 adds --phase-check <n> to the FIRST launched config only
#
# ONE ENSEMBLE + ONE SOURCE TYPE PER JOB. Claude does NOT submit.
#
# Build first (ON tuolumne):
#   cd ${ROOT}/Grid_sdm_build && source ${ROOT}/env.sh
#   GRID=${ROOT}/build ./compile_two_baryon_claude.sh baryon_variants_mom_corr_claude.cc
# Submit yourself (validation first: one config with the phase check, then the set):
#   ENS=m0p1000 SRC=smeared CLIST=1000 PHASECHECK=3 flux batch run_bvar_mom_contract_flux_claude.sh
#   ENS=m0p1000 SRC=smeared                         flux batch run_bvar_mom_contract_flux_claude.sh
#   ENS=m0p1000 SRC=point ; ENS=m0p0500 SRC=point|smeared ; ENS=m0p0100 SRC=point|smeared
#FLUX: -t 240m
#FLUX: --output=bvar_mom_contract_{{id}}.out
#FLUX: -q pbatch
#FLUX: -N 2
#FLUX: -n 8
#FLUX: -g 1
#FLUX: --exclusive

set -u

date; hostname
export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off
export HDF5_USE_FILE_LOCKING=FALSE

ROOT=/usr/workspace/lsd/matsumoto5/su4_32c
source ${ROOT}/env.sh
SRCDIR=${ROOT}/Grid_sdm_build/src/gauge_gen_Nc4
APP=${SRCDIR}/bin/baryon_variants_mom_corr_claude

# runtime text tables (read by the binary)
TERMS=${TERMS:-${SRCDIR}/baryon_variants_terms_all_claude.txt}   # ALL 1225 pairs (momentum path needs them)
OPS=${OPS:-${SRCDIR}/baryon_variants_ops_claude.txt}             # legacy p=0 operators (as production)
MOMOPS=${MOMOPS:-${SRCDIR}/baryon_mom_ops_claude.txt}             # Aaron's irrep operators, 81 momenta

# ---------------- ensemble (ONE per job): light masses of the point-source update ----------------
ENS=${ENS:-m0p1000}
case "${ENS}" in
    *m0p1*|*10p865*)  betastr=10p865; massstr=0p1000 ;;
    *m0p05*|*10p840*) betastr=10p840; massstr=0p0500 ;;
    *m0p01*|*10p800*) betastr=10p800; massstr=0p0100 ;;
    *) echo "ERROR: ENS='${ENS}' not one of m0p1000/m0p0500/m0p0100 (momentum pass = light masses only)" >&2; exit 1 ;;
esac
cfgfilename=conf_nc4nf1_2448_b${betastr}_m${massstr}
LUSTRE=/p/lustre5/matsumoto5
CONFROOT=${CONFROOT:-${LUSTRE}/conf_nc4nf1_2448}

# ---------------- knobs ----------------
SRC=${SRC:-point}                   # which dump: point | smeared  (-> bvar_2448_<ens>[_point])
case "${SRC}" in smeared) TAG="" ;; point) TAG="_point" ;; *) echo "ERROR: SRC smeared|point" >&2; exit 1 ;; esac
PAIRS=${PAIRS:-sector}              # legacy p=0 pair set, as production (only used without --no-p0)
SINKSMEAR=${SINKSMEAR:-1}           # 1 = add smeared-sink copies (--config gauge --sink-smear 3.0 40)
SS_W=${SS_W:-3.0}; SS_N=${SS_N:-40} # MUST match the source smearing (as production)
MOMRAW=${MOMRAW:-0}                 # 1 = also write per-component Cmomraw_* (~30 MB/config/sink)
PHASECHECK=${PHASECHECK:-0}         # >0 = --phase-check <n> on the first launched config
NSLOTS=${NSLOTS:-8}                 # concurrent single-GCD configs (= -n above)

DUMPDIR=${DUMPDIR:-${LUSTRE}/bvar_2448_b${betastr}_m${massstr}${TAG}}
SSTAG=""; [ "${SINKSMEAR}" = "1" ] && SSTAG="_ss"
P0DIR=${P0DIR:-${LUSTRE}/obs_bvar_2448_b${betastr}_m${massstr}${TAG}${SSTAG}}     # production p=0 output
OUTDIR=${OUTDIR:-${LUSTRE}/obs_bvar_2448_b${betastr}_m${massstr}${TAG}${SSTAG}_mom}
mkdir -p "${OUTDIR}"

LATT="24.24.24.48"
MPIGRID="1.1.1.1"                   # REQUIRED: the FFT host copy needs a single rank
PARAMS_GRID="--grid ${LATT} --mpi ${MPIGRID} --threads 8 --accelerator-threads 8 --shm 2048 --shm-mpi 1"

for f in "${APP}" "${TERMS}" "${OPS}" "${MOMOPS}"; do
    [ -e "${f}" ] || { echo "ERROR: missing: ${f}" >&2; exit 1; }
done
[ -x "${APP}" ] || { echo "ERROR: binary not executable: ${APP} -- build on tuolumne first" >&2; exit 1; }

echo "--start " "$(date)" "$(date +%s)"
echo "ensemble=${cfgfilename}  src=${SRC}  pairs=${PAIRS}  sinksmear=${SINKSMEAR} (w=${SS_W} N=${SS_N})  momraw=${MOMRAW}  phasecheck=${PHASECHECK}  NSLOTS=${NSLOTS}"
echo "dumps  = ${DUMPDIR}"
echo "p0 ref = ${P0DIR}"
echo "outdir = ${OUTDIR}"

# config selection: CLIST override, else every dump present.
if [ -n "${CLIST:-}" ]; then
    dumps=""; for i in ${CLIST}; do dumps="${dumps} ${DUMPDIR}/${cfgfilename}.${i}.pdfull.lime"; done
else
    dumps=$(ls "${DUMPDIR}"/${cfgfilename}.*.pdfull.lime 2>/dev/null)
fi
[ -n "${dumps}" ] || { echo "no dumps found in ${DUMPDIR} (run Stage 1 first)"; exit 0; }

running=0
first=1
for dump in ${dumps}; do
    base=$(basename "${dump}"); stem=${base%.pdfull.lime}; idx=${stem##*.}
    gauge="${CONFROOT}/${cfgfilename}/${cfgfilename}_lat.${idx}"
    out="${OUTDIR}/bvarmom.all.${idx}.h5"
    tmp="${OUTDIR}/bvarmom.all.${idx}.h5.inprogress.h5"
    log="${OUTDIR}/bvarmom.all.${idx}.log"
    p0="${P0DIR}/bvar.all.${idx}.h5"

    if [ ! -f "${dump}" ]; then echo "  WARN no dump, skip idx ${idx}: ${dump}"; continue; fi
    if [ -s "${out}" ]; then echo "  SKIP (exists): ${out}"; continue; fi

    extra=""
    if [ "${SINKSMEAR}" = "1" ]; then
        [ -f "${gauge}" ] || { echo "  WARN no gauge for sink smear, skip idx ${idx}: ${gauge}"; continue; }
        extra="--config ${gauge} --sink-smear ${SS_W} ${SS_N}"
    fi
    p0mode="p0 written"
    if [ -s "${p0}" ]; then
        extra="${extra} --no-p0"
        p0mode="--no-p0 (production ${p0} exists)"
    fi
    [ "${MOMRAW}" = "1" ] && extra="${extra} --mom-raw"
    if [ "${PHASECHECK}" -gt 0 ] && [ "${first}" = "1" ]; then
        extra="${extra} --phase-check ${PHASECHECK}"
    fi
    first=0

    echo "  -> config ${idx}  [${p0mode}]"
    # one task, one GCD (per-task opts only: mixing --tasks-per-node with -g is rejected by flux);
    # rename the finished file only if the binary exited cleanly
    ( flux run -n1 -g1 \
        "${APP}" --src "${dump}" --out "${tmp}" \
        --terms "${TERMS}" --ops "${OPS}" --op all --pairs "${PAIRS}" \
        --momops "${MOMOPS}" ${extra} \
        ${PARAMS_GRID} > "${log}" 2>&1 \
      && mv "${tmp}" "${out}" \
      || echo "  FAILED config ${idx} (see ${log})" ) &

    running=$(( running + 1 ))
    if [ "${running}" -ge "${NSLOTS}" ]; then wait -n; running=$(( running - 1 )); fi
done
wait
echo "--end " "$(date)" "$(date +%s)"

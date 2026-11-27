#!/bin/bash
# validate_bvar_32c_flux_claude.sh
# 32c VALIDATION of the bvar measurement chain (plan grid_claude/bvar_32c_port_impl_plan_claude.md, chunk 3; user
# go-ahead 2026-11-24): ONE config per mass, both source types, in ONE 8-node job.
#   phase 1: stage A (run_bvar_dumper_32c_flux_claude.sh) for every mass x {point, smeared}, sequentially (8 nodes each);
#   phase 2: stage B (run_bvar_mom_contract_32c_flux_claude.sh) for all of them CONCURRENTLY (one GCD each),
#            with PHASECHECK=3 (FFT vs explicit projection) and MEMCHECK=1 (/usr/bin/time -v, max RSS);
#   summary: dump sizes, stage-A wall time + CG lines, stage-B max RSS + status, per mass and source.
# Config per mass: VAL_TRAJ_<mass> if set, else the LAST ckpoint_lat.<traj> present. The outputs are the production
# files of that config (the production jobs later self-skip it), so nothing is wasted.
# Afterwards (local): pull the 12 bvarmom.all.<traj>.h5 + logs and run grid_claude/bvar_32c_validate_check_claude.py.
# Claude does NOT submit.
#
# Submit yourself (from this directory, on tuolumne):
#   flux batch validate_bvar_32c_flux_claude.sh
# Options: ENS_LIST="m0p4000 m0p0100" (default all six), VAL_TRAJ_m0p0100=400 etc.
#FLUX: -t 480m
#FLUX: --output=validate_bvar_32c_{{id}}.out
#FLUX: -q pbatch
#FLUX: -N 8
#FLUX: -n 32
#FLUX: -g 1
#FLUX: --exclusive

set -u

date
hostname
# ORIGINAL: HERE=$(cd "$(dirname "$0")" && pwd)
# flux batch runs a SPOOLED COPY of this script ($0 in a temporary directory), so $0 does not locate the drivers
# (validation f3nDujNjgJFR, 2026-11-26: nothing ran). Use the submit directory (flux keeps the cwd) or BVAR32C_DIR.
HERE=${BVAR32C_DIR:-$(pwd)}
for f in run_bvar_dumper_32c_flux_claude.sh run_bvar_mom_contract_32c_flux_claude.sh; do
    if [ ! -f "${HERE}/${f}" ]; then
        echo "ERROR: ${HERE}/${f} not found -- submit from the run_bvar_32c directory or set BVAR32C_DIR" >&2
        exit 1
    fi
done
echo "driver dir = ${HERE}"
LUSTRE=/p/lustre5/matsumoto5
ENS_LIST=${ENS_LIST:-"m0p0100 m0p0500 m0p1000 m0p2000 m0p3000 m0p4000"}

# mass label -> m (stream directory 32_64_<m>)
mval() {
    case "$1" in
        *m0p01*)
            echo 0.01
            ;;
        *m0p05*)
            echo 0.05
            ;;
        *m0p1*)
            echo 0.1
            ;;
        *m0p2*)
            echo 0.2
            ;;
        *m0p3*)
            echo 0.3
            ;;
        *m0p4*)
            echo 0.4
            ;;
    esac
}

# ---------------- choose one config per mass ----------------
declare -A TRAJ
for e in ${ENS_LIST}; do
    m=$(mval "${e}")
    dir=$(ls -d ${LUSTRE}/32_64_${m}/beta*m${m} 2>/dev/null | head -1)
    var="VAL_TRAJ_${e}"
    t=${!var:-}
    if [ -z "${t}" ]; then
        t=$(ls "${dir}" | sed -n 's/^ckpoint_lat\.\([0-9][0-9]*\)$/\1/p' | sort -n | tail -1)
    fi
    TRAJ[${e}]=${t}
    echo "validation config: ${e} -> ${dir}/ckpoint_lat.${t}"
done

# ---------------- phase 1: stage A (8 nodes each, sequential) ----------------
echo "=== phase 1: stage A dumps  $(date)"
for e in ${ENS_LIST}; do
    for s in point smeared; do
        echo "--- stage A ${e} ${s} traj ${TRAJ[${e}]}"
        ENS=${e} SRC=${s} CLIST=${TRAJ[${e}]} bash "${HERE}/run_bvar_dumper_32c_flux_claude.sh" 2>&1 | tee "${HERE}/val32c_A_${e}_${s}_claude.log"
    done
done

# ---------------- phase 2: stage B (one GCD each, all concurrently) ----------------
echo "=== phase 2: stage B contractions  $(date)"
for e in ${ENS_LIST}; do
    for s in point smeared; do
        ENS=${e} SRC=${s} CLIST=${TRAJ[${e}]} NSLOTS=1 PHASECHECK=3 MEMCHECK=1 \
            bash "${HERE}/run_bvar_mom_contract_32c_flux_claude.sh" > "${HERE}/val32c_B_${e}_${s}_claude.log" 2>&1 &
    done
done
wait

# ---------------- summary ----------------
echo "=== summary  $(date)"
for e in ${ENS_LIST}; do
    m=$(mval "${e}")
    for s in point smeared; do
        tag=""
        if [ "${s}" = "point" ]; then
            tag="_point"
        fi
        t=${TRAJ[${e}]}
        echo "--- ${e} ${s} traj ${t}"
        ls -l ${LUSTRE}/bvar_3264_b*_m${e#m}${tag}/*_m${e#m}.${t}.pdfull.lime 2>/dev/null | awk '{printf "    dump %s  %.2f GB\n", $9, $5 / 1e9}'
        grep -h "done in" "${HERE}/val32c_A_${e}_${s}_claude.log" | sed 's/^/    A: /'
        grep -h -i "converged\|iterations" ${LUSTRE}/bvar_3264_b*_m${e#m}${tag}/*_m${e#m}.${t}.pdfull.log 2>/dev/null | tail -2 | sed 's/^/    A CG: /'
        blog=$(ls ${LUSTRE}/obs_bvar_3264_b*_m${e#m}${tag}_ss_mom/bvarmom.all.${t}.log 2>/dev/null | head -1)
        if [ -n "${blog}" ]; then
            grep -h "Maximum resident set size\|Elapsed (wall clock)\|Exit status" "${blog}" | sed 's/^/    B: /'
            grep -h -i "phase" "${blog}" | tail -3 | sed 's/^/    B phase: /'
        fi
        grep -h "FAILED\|SKIP\|ERROR" "${HERE}/val32c_B_${e}_${s}_claude.log" | sed 's/^/    B status: /'
    done
done
echo "=== done  $(date)"

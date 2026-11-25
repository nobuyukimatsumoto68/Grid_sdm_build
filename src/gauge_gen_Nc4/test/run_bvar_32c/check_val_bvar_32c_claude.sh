#!/bin/bash
# check_val_bvar_32c_claude.sh
# READ-ONLY check of the 32c bvar validation job (validate_bvar_32c_flux_claude.sh), run ON tuolumne from this
# directory after the job finished. Plan: grid_claude/bvar_32c_port_impl_plan_claude.md (Sec. 7). Prints, per mass and
# source: NERSC read of ckpoint_lat (checksum / plaquette vs header, mismatch), step-A wall time per config + last CG
# lines, dump size; step-B max RSS, wall, exit status, phase-check summary, FAILED lines; then every error-like line.
# No writes, no deletions.
#
# Usage (on tuolumne):
#   cd /usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build/src/gauge_gen_Nc4/test/run_bvar_32c
#   bash check_val_bvar_32c_claude.sh 2>&1 | tee check_val_bvar_32c_claude.log

set -u
LUSTRE=/p/lustre5/matsumoto5
HERE=$(cd "$(dirname "$0")" && pwd)
ENS_LIST=${ENS_LIST:-"m0p0100 m0p0500 m0p1000 m0p2000 m0p3000 m0p4000"}

echo "=== validation job output (summary section)"
out=$(ls -t "${HERE}"/validate_bvar_32c_*.out 2>/dev/null | head -1)
echo "file: ${out}"
if [ -n "${out}" ]; then
    sed -n '/=== summary/,$p' "${out}"
fi

for e in ${ENS_LIST}; do
    ms=${e#m}
    for s in point smeared; do
        tag=""
        if [ "${s}" = "point" ]; then
            tag="_point"
        fi
        echo
        echo "=== ${e} ${s}"
        for dlog in ${LUSTRE}/bvar_3264_b*_m${ms}${tag}/*.pdfull.log; do
            if [ ! -f "${dlog}" ]; then
                echo "  no step-A log"
                continue
            fi
            echo "  step-A log: ${dlog}"
            grep -h "NERSC Configuration" "${dlog}" | head -3 | sed 's/^/    /'
            grep -h -i "Plaquette mismatch" "${dlog}" | sed 's/^/    MISMATCH: /'
            grep -h "done in" "${dlog}" | sed 's/^/    /'
            grep -h -i "converged\|true residual\|iterations" "${dlog}" | tail -3 | sed 's/^/    CG: /'
            ls -l "${dlog%.pdfull.log}.pdfull.lime" 2>/dev/null | awk '{printf "    dump %.2f GB\n", $5 / 1e9}'
        done
        for blog in ${LUSTRE}/obs_bvar_3264_b*_m${ms}${tag}_ss_mom/bvarmom.all.*.log; do
            if [ ! -f "${blog}" ]; then
                echo "  no step-B log"
                continue
            fi
            echo "  step-B log: ${blog}"
            grep -h "Maximum resident set size\|Elapsed (wall clock)\|Exit status" "${blog}" | sed 's/^/    /'
            grep -h "phase-check summary\|phase-check: max rel" "${blog}" | tail -2 | sed 's/^/    /'
            h5="${blog%.log}.h5"
            if [ -s "${h5}" ]; then
                ls -l "${h5}" | awk '{printf "    output %s %.1f MB\n", $9, $5 / 1e6}'
            else
                echo "    NO OUTPUT h5 (inprogress or failed)"
            fi
        done
    done
done

echo
echo "=== error-like lines in all validation logs"
grep -l -i "error\|abort\|segmentation\|mismatch\|FAILED\|out of memory" \
    "${HERE}"/val32c_*_claude.log ${LUSTRE}/bvar_3264_*/*.pdfull.log ${LUSTRE}/obs_bvar_3264_*_ss_mom/*.log 2>/dev/null \
    | while read -r f; do
        echo "--- ${f}"
        grep -h -i "error\|abort\|segmentation\|mismatch\|FAILED\|out of memory" "${f}" | head -5
    done
echo "=== done"

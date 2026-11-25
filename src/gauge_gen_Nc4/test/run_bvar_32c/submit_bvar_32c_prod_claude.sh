#!/bin/bash
# submit_bvar_32c_prod_claude.sh
# PRODUCTION submitter for the 32c bvar measurement (plan grid_claude/bvar_32c_port_impl_plan_claude.md, Sec. 7). The USER
# runs it on tuolumne (Claude never submits). For every mass in ENS_LIST and source in SRC_LIST (all in PARALLEL):
#   step A (run_bvar_dumper_32c_flux_claude.sh, 8 nodes, 480 m): a CHAIN of NA jobs, each --dependency=afterany on the
#     previous one (the driver self-skips finished configs, so every link continues where the last stopped);
#   step B (run_bvar_mom_contract_32c_flux_claude.sh, 2 nodes, 480 m): ONE job --dependency=afterany on the LAST step-A
#     job of that chain (contracts every dump present; resubmit later for any config A had not finished).
# Every config: the drivers take every ckpoint_lat.<traj> present when the step-A job STARTS.
# Chain length per job wall: NA = ceil(100 * t_A / (28800 - 300)), t_A = validation seconds per config.
#   Defaults (scaled 24c guesses, REPLACE with the validation timings): NA_HEAVY=3 (t_A ~ 10 min), NA_LIGHT=10.
# PROP_TPT_SECONDS (per-config wall budget of the blocker) is passed through if set (e.g. 1.3 x the validation time).
# DRYRUN=1 prints the commands without submitting.
#
# Usage (on tuolumne, from this directory):
#   DRYRUN=1 bash submit_bvar_32c_prod_claude.sh
#   NA_HEAVY=3 NA_LIGHT=10 bash submit_bvar_32c_prod_claude.sh 2>&1 | tee submit_bvar_32c_prod_claude.log

set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ENS_LIST=${ENS_LIST:-"m0p0100 m0p0500 m0p1000 m0p2000 m0p3000 m0p4000"}
SRC_LIST=${SRC_LIST:-"point smeared"}
NA_HEAVY=${NA_HEAVY:-3}
NA_LIGHT=${NA_LIGHT:-10}
DRYRUN=${DRYRUN:-0}

submit() {
    # $1 = dependency jobid or "", rest = script; prints the new jobid
    dep="$1"
    shift
    if [ "${DRYRUN}" = "1" ]; then
        if [ -n "${dep}" ]; then
            echo "DRYRUN: ENS=${ENS} SRC=${SRC} flux batch --dependency=afterany:${dep} $*" >&2
        else
            echo "DRYRUN: ENS=${ENS} SRC=${SRC} flux batch $*" >&2
        fi
        echo "f_dryrun_${RANDOM}"
        return
    fi
    if [ -n "${dep}" ]; then
        flux batch --dependency=afterany:${dep} "$@"
    else
        flux batch "$@"
    fi
}

for e in ${ENS_LIST}; do
    na=${NA_HEAVY}
    case "${e}" in
        *m0p01*|*m0p05*)
            na=${NA_LIGHT}
            ;;
    esac
    for s in ${SRC_LIST}; do
        export ENS=${e}
        export SRC=${s}
        prev=""
        n=1
        while [ "${n}" -le "${na}" ]; do
            prev=$(submit "${prev}" "${HERE}/run_bvar_dumper_32c_flux_claude.sh")
            echo "${e} ${s} step A link ${n}/${na}: ${prev}"
            n=$(( n + 1 ))
        done
        jb=$(submit "${prev}" "${HERE}/run_bvar_mom_contract_32c_flux_claude.sh")
        echo "${e} ${s} step B (after ${prev}): ${jb}"
    done
done

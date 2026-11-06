#!/bin/bash
# run_glueball_tuo_claude.sh
# Production driver for the Grid glueball code (plan glueball_grid_impl_plan_claude.md, Chunk 8.4), modeled on
# Grid/scripts/run_wflow_q_h5_tuo_claude.sh. For each stream (tag, conf dir, start, step [, end]) the existing configs
# on the start/step grid are cut into ranges of NPER configs; one 1-node job per range (4 single-GPU processes,
# submit_glueball_tuo_claude.sh). A range whose configs all have OUTDIR/glue.<conf>.h5 is skipped; otherwise it is
# (re)submitted -- the binary skips finished configs, so a rerun only does what is missing.
# OUTDIR = /p/lustre5/matsumoto5/obs_glue_2448_<tag>/glue_h5/ (orig + cont share it: config numbers do not overlap).
# Benchmark 2026-11-06: ~49 s per config per GPU -> NPER = 600 is ~2 h per job (4 h wall).
# Usage (tuolumne login, from su4_32c/Grid_sdm_build/src/gauge_gen_Nc4, after building glueball_h5_claude):
#   DRYRUN=1 bash run_glueball_tuo_claude.sh
#   ONLY=b11p045_m0p4000 bash run_glueball_tuo_claude.sh
#   bash run_glueball_tuo_claude.sh                    # all streams
# Claude does NOT submit -- run this yourself.

set -u

NPER=${NPER:-600}
WALL=${WALL:-240m}
QUEUE=${QUEUE:-pbatch}
ONLY=${ONLY:-}
DRYRUN=${DRYRUN:-0}
OBSTOP=/p/lustre5/matsumoto5
ORIG=/p/lustre5/matsumoto5/conf_nc4nf1_2448
CONT=/p/lustre5/matsumoto5/conf_nc4nf1_2448_cont
script=submit_glueball_tuo_claude.sh

# stream table: tag  topdir  start  step  [end]   (end empty = last existing config)
# heavy: traj >= 1000 every 20 (as the meson analysis); m0.05 every 20, m0.01 every 4 (user 2026-11-06).
# START POINTS / _cont STREAMS OF m0.05 AND m0.01: PROPOSED, to be confirmed by the user.
streams=(
  "b11p045_m0p4000 ${ORIG} 1000 20"
  "b11p035_m0p3000 ${ORIG} 1000 20"
  "b10p990_m0p2000 ${ORIG} 1000 20"
  "b10p865_m0p1000 ${ORIG} 1000 20"
  "b10p865_m0p1000 ${CONT} 4660 20"
  "b10p840_m0p0500 ${ORIG} 500 20"
  "b10p840_m0p0500 ${CONT} 2140 20"
  "b10p800_m0p0100 ${ORIG} 4 4"
  "b10p800_m0p0100 ${CONT} 760 4"
)

basedir=$(pwd)
for line in "${streams[@]}"; do
    read -r tag top start step end <<< "${line}"
    end=${end:-}
    if [ -n "${ONLY}" ] && [ "${ONLY}" != "${tag}" ]; then
        continue
    fi
    cfgname=conf_nc4nf1_2448_${tag}
    cdir=${top}/${cfgname}
    prefix=${cfgname}_lat
    outdir=${OBSTOP}/obs_glue_2448_${tag}/glue_h5
    confs=$(ls ${cdir}/${prefix}.* 2>/dev/null | sed 's/.*\.//' | grep -E '^[0-9]+$' | sort -n | \
            awk -v s="${start}" -v k="${step}" -v e="${end}" '$1 >= s && ($1 - s) % k == 0 && (e == "" || $1 <= e + 0)')
    if [ -z "${confs}" ]; then
        echo "no configs: ${cdir} from ${start} step ${step}"
        continue
    fi
    ntot=$(echo "${confs}" | wc -l)
    echo "stream: ${tag} $(basename ${top}) from ${start} step ${step}: ${ntot} configs ($(echo "${confs}" | head -1)..$(echo "${confs}" | tail -1))"
    # ranges of NPER existing configs: "first last" (gaps inside a range are fine, the binary skips them)
    ranges=$(echo "${confs}" | awk -v n="${NPER}" 'NR % n == 1 { a = $1 } { p = $1 } NR % n == 0 { print a, p } END { if (NR % n != 0) print a, p }')
    while read -r a b; do
        todo=$(echo "${confs}" | awk -v a="${a}" -v b="${b}" '$1 >= a && $1 <= b' | while read -r c; do
                   [ -f ${outdir}/glue.${c}.h5 ] || echo ${c}
               done | wc -l)
        if [ "${todo}" -eq 0 ]; then
            echo "done : ${tag} $(basename ${top}) ${a}..${b}"
            continue
        fi
        echo "range: ${tag} $(basename ${top}) ${a}..${b} step ${step} (${todo} to do) -> ${outdir}"
        if [ "${DRYRUN}" = "1" ]; then
            continue
        fi
        mkdir -p ${outdir}
        flux batch --job-name=glue_${tag} -q ${QUEUE} -t ${WALL} \
            --env=CONF_PATH=${cdir} --env=CONF_PREFIX=${prefix} --env=START=${a} --env=END=${b} \
            --env=SKIP=${step} --env=OUTDIR=${outdir} ${basedir}/${script}
    done <<< "${ranges}"
done

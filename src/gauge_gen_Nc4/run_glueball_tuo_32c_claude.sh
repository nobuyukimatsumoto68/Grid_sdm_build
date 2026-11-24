#!/bin/bash
# run_glueball_tuo_32c_claude.sh
# 32^3 x 64 COPY of run_glueball_tuo_claude.sh (kept unchanged); plan glueball_32c_impl_plan_claude.md. EVERY existing
# ckpoint_lat.<traj> of the six 32c streams (user 2026-11-24), /p/lustre5/matsumoto5/32_64_<m>/beta<b>m<m>/;
# output /p/lustre5/matsumoto5/obs_glue_3264_<tag>/glue_h5/. Estimate ~150 s/config/GPU -> NPER 200 = one ~2 h job per
# ensemble. GLUE_TRACE=1 in the environment is passed to the jobs (needs the rebuilt binary, plan (h)).
# Production driver for the Grid glueball code (plan glueball_grid_impl_plan_claude.md, Chunk 8.4), modeled on
# Grid/scripts/run_wflow_q_h5_tuo_claude.sh. For each stream (tag, conf dir, start, step [, end]) the existing configs
# on the start/step grid are cut into ranges of NPER configs; one 1-node job per range (4 single-GPU processes,
# submit_glueball_tuo_32c_claude.sh). A range whose configs all have OUTDIR/glue.<conf>.h5 is skipped; otherwise it is
# (re)submitted -- the binary skips finished configs, so a rerun only does what is missing.
# OUTDIR = /p/lustre5/matsumoto5/obs_glue_3264_<tag>/glue_h5/ (one stream dir per ensemble on 32c).
# 24c benchmark 2026-11-06: ~49 s per config per GPU; 32c estimate x3.2 -> NPER = 200 is ~2 h per job (4 h wall).
# Usage (tuolumne login, from su4_32c/Grid_sdm_build/src/gauge_gen_Nc4, after building glueball_h5_claude; validate with
# the pdebug example in submit_glueball_tuo_32c_claude.sh first):
#   DRYRUN=1 bash run_glueball_tuo_32c_claude.sh
#   ONLY=b11p045_m0p4000 bash run_glueball_tuo_32c_claude.sh
#   GLUE_TRACE=1 bash run_glueball_tuo_32c_claude.sh     # all streams, traced
# Claude does NOT submit -- run this yourself.

set -u

# NPER=${NPER:-600}
NPER=${NPER:-200}
WALL=${WALL:-240m}
QUEUE=${QUEUE:-pbatch}
ONLY=${ONLY:-}
DRYRUN=${DRYRUN:-0}
GLUE_TRACE=${GLUE_TRACE:-0}
OBSTOP=/p/lustre5/matsumoto5
CONFROOT=/p/lustre5/matsumoto5
script=submit_glueball_tuo_32c_claude.sh

# stream table: tag  confdir (under CONFROOT)  start  step  [end]   (end empty = last existing config)
# Every config (user 2026-11-24; same list as the 32c disc, run_disc_mixedprec_32c_tuo_claude.sh).
streams=(
  "b11p045_m0p4000 32_64_0.4/beta11.045m0.4 20 20"
  "b11p035_m0p3000 32_64_0.3/beta11.035m0.3 20 20"
  "b10p990_m0p2000 32_64_0.2/beta10.99m0.2 20 20"
  "b10p865_m0p1000 32_64_0.1/beta10.865m0.1 20 20"
  "b10p840_m0p0500 32_64_0.05/beta10.84m0.05 4 4"
  "b10p800_m0p0100 32_64_0.01/beta10.8m0.01 4 4"
)

basedir=$(pwd)
for line in "${streams[@]}"; do
    read -r tag sub start step end <<< "${line}"
    end=${end:-}
    if [ -n "${ONLY}" ] && [ "${ONLY}" != "${tag}" ]; then
        continue
    fi
    top=${CONFROOT}/${sub}
    cdir=${top}
    prefix=ckpoint_lat
    outdir=${OBSTOP}/obs_glue_3264_${tag}/glue_h5
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
        flux batch --job-name=glue32_${tag} -q ${QUEUE} -t ${WALL} --env=GLUE_TRACE=${GLUE_TRACE} \
            --env=CONF_PATH=${cdir} --env=CONF_PREFIX=${prefix} --env=START=${a} --env=END=${b} \
            --env=SKIP=${step} --env=OUTDIR=${outdir} ${basedir}/${script}
    done <<< "${ranges}"
done

#!/bin/bash
# rsync_pull_glue24c_check_claude.sh
# Pull ONE tuolumne Grid glueball output (m0.4 config 1000, benchmark) and its gauge config (~680 MB) for the local
# 24c NumPy cross-check of the HIP build (glueball_grid_impl_plan_claude.md, Chunk 8.4).
# Run LOCALLY (this machine), same login as rsync_pull_wflow_q_claude.sh:
#   bash /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/rsync_pull_glue24c_check_claude.sh
set -u
REMOTE=matsumoto5@oslic.llnl.gov
DEST=/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/glue24c_check_claude
mkdir -p "${DEST}"
rsync -avhP "${REMOTE}:/p/lustre5/matsumoto5/obs_glue_2448_b11p045_m0p4000/bench_claude/glue.1000.h5" "${DEST}/" || exit 1
rsync -avhP "${REMOTE}:/p/lustre5/matsumoto5/conf_nc4nf1_2448/conf_nc4nf1_2448_b11p045_m0p4000/conf_nc4nf1_2448_b11p045_m0p4000_lat.1000" "${DEST}/" || exit 1
ls -la "${DEST}"

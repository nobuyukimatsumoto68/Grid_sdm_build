#!/bin/bash
#
# Quick local recheck after adding --maxtime to the stress mode (glueball_grid_impl_plan_claude.md, Chunk 8.7):
#   PHASE 1 compile glueball_h5_claude.cc
#   PHASE 2 stress prod, NCALLS huge, --maxtime 5: must stop by itself ("# stress maxtime 5 s reached") with rel diff 0
#   PHASE 3 stress split, 200 calls: must finish with rel diff 0
# ~3 min. Uses ONE GPU (default 0; GPU=1 to switch).
# Run:  bash /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/tmp_glueball_stress2_claude.sh
# Log:  /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/tmp_glueball_stress2_claude.log
#
set -u
set -o pipefail

HERE=/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4
ROOT=/mnt/baracuda_14/grid_claude/Grid_sdm_build
LOG=$HERE/tmp_glueball_stress2_claude.log
EXE=$HERE/bin/glueball_h5_claude
CDIR=/mnt/baracuda_14/grid_claude/16c
GOPT="--grid 16.16.16.8 --mpi 1.1.1.1"

export PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/bin:${PATH}
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/lib:/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
export OMPI_CXX=g++-12
export OMPI_CC=gcc-12
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=${GPU:-0}

echo "================ PHASE 1: compile ($(date)) ================" 2>&1 | tee "$LOG"
cd "$HERE"
GRID=/mnt/baracuda_14/grid_claude/build_cuda129 "$ROOT/compile_two_baryon_claude.sh" glueball_h5_claude.cc > $HERE/tmp_glueball_stress2_build_claude.log 2>&1
RC=$?
grep -nE -A3 " error|error:" $HERE/tmp_glueball_stress2_build_claude.log | head -40 2>&1 | tee -a "$LOG"
echo "compile exit code: $RC" 2>&1 | tee -a "$LOG"
if [ "$RC" -ne 0 ]; then
    tail -40 $HERE/tmp_glueball_stress2_build_claude.log 2>&1 | tee -a "$LOG"
    echo "COMPILE FAILED" 2>&1 | tee -a "$LOG"
    exit 1
fi

echo "================ PHASE 2: prod with --maxtime 5 ($(date)) ================" 2>&1 | tee -a "$LOG"
SLOG=$HERE/tmp_glueball_stress2_prod_claude.log
"$EXE" --stress prod --dir "$CDIR" --prefix ckpoint_lat --conf 1000 --ncalls 100000000 --maxtime 5 --level 2 $GOPT > $SLOG 2>&1
echo "exit code: $?" 2>&1 | tee -a "$LOG"
grep -E "# stress (maxtime|done|WARNING)" $SLOG 2>&1 | tee -a "$LOG"

echo "================ PHASE 3: split, 200 calls ($(date)) ================" 2>&1 | tee -a "$LOG"
SLOG=$HERE/tmp_glueball_stress2_split_claude.log
"$EXE" --stress split --dir "$CDIR" --prefix ckpoint_lat --conf 1000 --ncalls 200 --level 2 $GOPT > $SLOG 2>&1
echo "exit code: $?" 2>&1 | tee -a "$LOG"
grep -E "# stress (done|WARNING)" $SLOG 2>&1 | tee -a "$LOG"
echo "DONE ($(date))" 2>&1 | tee -a "$LOG"

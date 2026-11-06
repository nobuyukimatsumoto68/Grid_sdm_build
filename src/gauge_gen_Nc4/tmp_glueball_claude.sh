#!/bin/bash
#
# Chunk 8.2/8.3 local test of the Grid glueball code (glueball_grid_impl_plan_claude.md):
#   PHASE 1 compile glueball_h5_claude.cc against build_cuda129 (Nc = 4, CUDA 12.9, HDF5)
#   PHASE 2 --selftest: HiRep projection kernel vs an Eigen SVD polar factor on random 8^4 matrices
#   PHASE 3 16c configs 1000..4000 (16^3 x 8, 1 GPU) -> glue16c_test_claude/glue.<conf>.h5; rerun (must skip all)
#   PHASE 4 aggregate + compare with the HiRep reference (hirep_su4 Chunk 8.0), tolerance 1e-10
# Uses ONE GPU (default device 0; override: GPU=1 bash tmp_glueball_claude.sh). Run when that device is not crowded.
# Run:   bash /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/tmp_glueball_claude.sh
# Log:   /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/tmp_glueball_claude.log
#
set -u
set -o pipefail

HERE=/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4
ROOT=/mnt/baracuda_14/grid_claude/Grid_sdm_build
LOG=$HERE/tmp_glueball_claude.log
EXE=$HERE/bin/glueball_h5_claude
CDIR=/mnt/baracuda_14/grid_claude/16c
ODIR=$HERE/glue16c_test_claude
REF=/mnt/baracuda_14/hirep_su4/QuenchedSU4_claude/runs_claude/nersc16c_test/glue_nersc16c_claude.h5
AGG=$ODIR/agg_$(date +%Y%m%d_%H%M%S)_claude.h5

export PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/bin:${PATH}
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/lib:/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
export OMPI_CXX=g++-12
export OMPI_CC=gcc-12
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=${GPU:-0}

echo "================ PHASE 1: compile ($(date)) ================" 2>&1 | tee "$LOG"
cd "$HERE"
GRID=/mnt/baracuda_14/grid_claude/build_cuda129 "$ROOT/compile_two_baryon_claude.sh" glueball_h5_claude.cc > $HERE/tmp_glueball_build_claude.log 2>&1
RC=$?
grep -nE -A3 " error|error:" $HERE/tmp_glueball_build_claude.log | head -40 2>&1 | tee -a "$LOG"
echo "compile exit code: $RC" 2>&1 | tee -a "$LOG"
if [ "$RC" -ne 0 ]; then
    tail -40 $HERE/tmp_glueball_build_claude.log 2>&1 | tee -a "$LOG"
    echo "COMPILE FAILED" 2>&1 | tee -a "$LOG"
    exit 1
fi

echo "================ PHASE 2: selftest, GPU $CUDA_VISIBLE_DEVICES ($(date)) ================" 2>&1 | tee -a "$LOG"
"$EXE" --selftest --grid 8.8.8.8 --mpi 1.1.1.1 > $HERE/tmp_glueball_selftest_claude.log 2>&1
RC=$?
grep "selftest" $HERE/tmp_glueball_selftest_claude.log 2>&1 | tee -a "$LOG"
echo "selftest exit code: $RC" 2>&1 | tee -a "$LOG"
if [ "$RC" -ne 0 ]; then
    tail -30 $HERE/tmp_glueball_selftest_claude.log 2>&1 | tee -a "$LOG"
    echo "SELFTEST FAILED" 2>&1 | tee -a "$LOG"
    exit 1
fi

echo "================ PHASE 3: 16c configs 1000..4000 ($(date)) ================" 2>&1 | tee -a "$LOG"
for run in r1 r2; do
    t0=$SECONDS
    "$EXE" --dir "$CDIR" --prefix ckpoint_lat --start 1000 --end 4000 --step 1000 --outdir "$ODIR" \
        --depth 4 --alpha 0.33 --nlev 4 --grid 16.16.16.8 --mpi 1.1.1.1 > $HERE/tmp_glueball_16c_${run}_claude.log 2>&1
    RC=$?
    echo "---- [$run] exit $RC, $((SECONDS - t0)) s" 2>&1 | tee -a "$LOG"
    grep -E "# glueball_h5|# conf|# skipping|level 0|APE x|done" $HERE/tmp_glueball_16c_${run}_claude.log | head -24 2>&1 | tee -a "$LOG"
    if [ "$RC" -ne 0 ]; then
        tail -30 $HERE/tmp_glueball_16c_${run}_claude.log 2>&1 | tee -a "$LOG"
        echo "16c RUN FAILED" 2>&1 | tee -a "$LOG"
        exit 1
    fi
done

echo "================ PHASE 4: aggregate + compare with HiRep ($(date)) ================" 2>&1 | tee -a "$LOG"
python3 $HERE/glueball_h5_aggregate_claude.py "$ODIR" "$AGG" --L 16 --compare "$REF" --tol 1e-10 2>&1 | tee -a "$LOG"
RC=${PIPESTATUS[0]}
echo "compare exit code: $RC" 2>&1 | tee -a "$LOG"
echo "DONE ($(date))" 2>&1 | tee -a "$LOG"
exit "$RC"

#!/bin/bash
#
# Grid glueball hang investigation (glueball_grid_impl_plan_claude.md, Chunk 8.6), local check of the tracing build:
#   PHASE 1 compile glueball_h5_claude.cc (timing split + GLUE_TRACE) against build_cuda129
#   PHASE 2 16c config 1000 with GLUE_TRACE=1 -> glue16c_trace_claude/glue.1000.h5
#   PHASE 3 output must be IDENTICAL to the validated glue16c_test_claude/glue.1000.h5 (only logging changed);
#           show the new timing lines and the first / last trace lines
# Uses ONE GPU (default device 0; GPU=1 to switch). Run when that device is not crowded.
# Run:  bash /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/tmp_glueball_trace_claude.sh
# Log:  /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/tmp_glueball_trace_claude.log
#
set -u
set -o pipefail

HERE=/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4
ROOT=/mnt/baracuda_14/grid_claude/Grid_sdm_build
LOG=$HERE/tmp_glueball_trace_claude.log
EXE=$HERE/bin/glueball_h5_claude
CDIR=/mnt/baracuda_14/grid_claude/16c
ODIR=$HERE/glue16c_trace_claude

export PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/bin:${PATH}
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/lib:/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
export OMPI_CXX=g++-12
export OMPI_CC=gcc-12
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=${GPU:-0}

echo "================ PHASE 1: compile ($(date)) ================" 2>&1 | tee "$LOG"
cd "$HERE"
GRID=/mnt/baracuda_14/grid_claude/build_cuda129 "$ROOT/compile_two_baryon_claude.sh" glueball_h5_claude.cc > $HERE/tmp_glueball_trace_build_claude.log 2>&1
RC=$?
grep -nE -A3 " error|error:" $HERE/tmp_glueball_trace_build_claude.log | head -40 2>&1 | tee -a "$LOG"
echo "compile exit code: $RC" 2>&1 | tee -a "$LOG"
if [ "$RC" -ne 0 ]; then
    tail -40 $HERE/tmp_glueball_trace_build_claude.log 2>&1 | tee -a "$LOG"
    echo "COMPILE FAILED" 2>&1 | tee -a "$LOG"
    exit 1
fi

echo "================ PHASE 2: 16c config 1000 with GLUE_TRACE=1 ($(date)) ================" 2>&1 | tee -a "$LOG"
if [ -e "$ODIR/glue.1000.h5" ]; then
    echo "$ODIR/glue.1000.h5 exists: move it away yourself, then rerun" 2>&1 | tee -a "$LOG"
    exit 1
fi
GLUE_TRACE=1 "$EXE" --dir "$CDIR" --prefix ckpoint_lat --start 1000 --end 1000 --step 1000 --outdir "$ODIR" \
    --depth 4 --alpha 0.33 --nlev 4 --grid 16.16.16.8 --mpi 1.1.1.1 > $HERE/tmp_glueball_trace_run_claude.log 2>&1
RC=$?
echo "run exit code: $RC" 2>&1 | tee -a "$LOG"
grep -E "level [0-3]:|# conf" $HERE/tmp_glueball_trace_run_claude.log 2>&1 | tee -a "$LOG"
echo "trace lines: $(grep -c '# trace' $HERE/tmp_glueball_trace_run_claude.log)" 2>&1 | tee -a "$LOG"
grep "# trace" $HERE/tmp_glueball_trace_run_claude.log | head -3 2>&1 | tee -a "$LOG"
grep "# trace" $HERE/tmp_glueball_trace_run_claude.log | tail -3 2>&1 | tee -a "$LOG"
if [ "$RC" -ne 0 ]; then
    tail -20 $HERE/tmp_glueball_trace_run_claude.log 2>&1 | tee -a "$LOG"
    echo "RUN FAILED" 2>&1 | tee -a "$LOG"
    exit 1
fi

echo "================ PHASE 3: identical to the validated output? ================" 2>&1 | tee -a "$LOG"
python3 -I - "$ODIR/glue.1000.h5" "$HERE/glue16c_test_claude/glue.1000.h5" <<'EOF_PY' 2>&1 | tee -a "$LOG"
import sys
import h5py
import numpy as np
a = h5py.File(sys.argv[1], "r")["glue"][()]
b = h5py.File(sys.argv[2], "r")["glue"][()]
d = np.abs(a - b).max()
print("max |new - validated| = %.3e over %d numbers" % (d, a.size))
print("IDENTICAL" if d == 0.0 else "DIFFERENT (expected identical: only logging changed)")
EOF_PY
echo "DONE ($(date))" 2>&1 | tee -a "$LOG"

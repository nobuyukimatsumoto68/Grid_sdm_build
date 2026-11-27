#!/bin/bash
#
# Grid glueball hang reproduction, LOCAL (glueball_grid_impl_plan_claude.md, Chunk 8.7):
#   PHASE 1 compile glueball_h5_claude.cc (adds --stress; production path unchanged)
#   PHASE 2 regression: normal 16c config 1000 measurement must be IDENTICAL to the validated glue16c_test_claude/glue.1000.h5
#   PHASE 3 short stress, every mode (prod full kernel slice split), level 2, NSHORT calls each: must finish, max rel diff
#           0 for prod/full/kernel/slice and ~1e-15 for split (summation order)
#   PHASE 4 (only with LONG=1) long stress, mode ${LMODE:-prod}, level ${LLEV:-2}, NLONG calls (~6 ms per call on 16c)
#           -> does a hang reproduce on CUDA? A hang shows as a log that stops growing: then tell Claude, do not wait.
# Uses ONE GPU (default device 0; GPU=1 to switch). Run when that device is not crowded (nvidia-smi printed first).
# Run:  bash /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/tmp_glueball_stress_claude.sh
#       LONG=1 NLONG=300000 bash .../tmp_glueball_stress_claude.sh      (adds the ~30 min long run)
# Log:  /mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4/tmp_glueball_stress_claude.log
#       (stress output per mode: tmp_glueball_stress_<mode>_claude.log next to it)
#
set -u
set -o pipefail

HERE=/mnt/baracuda_14/grid_claude/Grid_sdm_build/src/gauge_gen_Nc4
ROOT=/mnt/baracuda_14/grid_claude/Grid_sdm_build
LOG=$HERE/tmp_glueball_stress_claude.log
EXE=$HERE/bin/glueball_h5_claude
CDIR=/mnt/baracuda_14/grid_claude/16c
ODIR=$HERE/glue16c_stress_regr_claude
NSHORT=${NSHORT:-2000}
NLONG=${NLONG:-300000}
LMODE=${LMODE:-prod}
LLEV=${LLEV:-2}
GOPT="--grid 16.16.16.8 --mpi 1.1.1.1"

export PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/bin:${PATH}
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/openmpi_cuda129/lib:/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
export OMPI_CXX=g++-12
export OMPI_CC=gcc-12
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=${GPU:-0}

echo "================ GPU status ($(date)) ================" 2>&1 | tee "$LOG"
nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv 2>&1 | tee -a "$LOG"

echo "================ PHASE 1: compile ($(date)) ================" 2>&1 | tee -a "$LOG"
cd "$HERE"
GRID=/mnt/baracuda_14/grid_claude/build_cuda129 "$ROOT/compile_two_baryon_claude.sh" glueball_h5_claude.cc > $HERE/tmp_glueball_stress_build_claude.log 2>&1
RC=$?
grep -nE -A3 " error|error:" $HERE/tmp_glueball_stress_build_claude.log | head -40 2>&1 | tee -a "$LOG"
echo "compile exit code: $RC" 2>&1 | tee -a "$LOG"
if [ "$RC" -ne 0 ]; then
    tail -40 $HERE/tmp_glueball_stress_build_claude.log 2>&1 | tee -a "$LOG"
    echo "COMPILE FAILED" 2>&1 | tee -a "$LOG"
    exit 1
fi

echo "================ PHASE 2: production regression, 16c config 1000 ($(date)) ================" 2>&1 | tee -a "$LOG"
if [ -e "$ODIR/glue.1000.h5" ]; then
    echo "$ODIR/glue.1000.h5 exists: move it away yourself, then rerun" 2>&1 | tee -a "$LOG"
    exit 1
fi
"$EXE" --dir "$CDIR" --prefix ckpoint_lat --start 1000 --end 1000 --step 1000 --outdir "$ODIR" \
    --depth 4 --alpha 0.33 --nlev 4 $GOPT > $HERE/tmp_glueball_stress_regr_claude.log 2>&1
RC=$?
echo "run exit code: $RC" 2>&1 | tee -a "$LOG"
grep -E "# conf" $HERE/tmp_glueball_stress_regr_claude.log 2>&1 | tee -a "$LOG"
if [ "$RC" -ne 0 ]; then
    tail -20 $HERE/tmp_glueball_stress_regr_claude.log 2>&1 | tee -a "$LOG"
    echo "RUN FAILED" 2>&1 | tee -a "$LOG"
    exit 1
fi
python3 -I - "$ODIR/glue.1000.h5" "$HERE/glue16c_test_claude/glue.1000.h5" <<'EOF_PY' 2>&1 | tee -a "$LOG"
import sys
import h5py
import numpy as np
a = h5py.File(sys.argv[1], "r")["glue"][()]
b = h5py.File(sys.argv[2], "r")["glue"][()]
d = np.abs(a - b).max()
print("max |new - validated| = %.3e over %d numbers" % (d, a.size))
print("IDENTICAL" if d == 0.0 else "DIFFERENT (expected identical: production path unchanged)")
EOF_PY

echo "================ PHASE 3: short stress, all modes, level 2, $NSHORT calls ($(date)) ================" 2>&1 | tee -a "$LOG"
for mode in prod full kernel slice split; do
    SLOG=$HERE/tmp_glueball_stress_${mode}_claude.log
    "$EXE" --stress $mode --dir "$CDIR" --prefix ckpoint_lat --conf 1000 --ncalls $NSHORT --level 2 $GOPT > $SLOG 2>&1
    RC=$?
    echo "--- mode $mode: exit code $RC" 2>&1 | tee -a "$LOG"
    grep -E "# stress (mode|done|WARNING)" $SLOG | head -5 2>&1 | tee -a "$LOG"
    grep -E "# stress heartbeat" $SLOG | tail -1 2>&1 | tee -a "$LOG"
    if [ "$RC" -ne 0 ]; then
        tail -15 $SLOG 2>&1 | tee -a "$LOG"
        echo "STRESS RUN FAILED (mode $mode)" 2>&1 | tee -a "$LOG"
        exit 1
    fi
done

if [ "${LONG:-0}" = "1" ]; then
    echo "================ PHASE 4: long stress, mode $LMODE, level $LLEV, $NLONG calls ($(date)) ================" 2>&1 | tee -a "$LOG"
    SLOG=$HERE/tmp_glueball_stress_long_${LMODE}_claude.log
    echo "watch: tail -1 $SLOG   (a log that stops growing = hang; tell Claude)" 2>&1 | tee -a "$LOG"
    "$EXE" --stress $LMODE --dir "$CDIR" --prefix ckpoint_lat --conf 1000 --ncalls $NLONG --level $LLEV $GOPT > $SLOG 2>&1
    RC=$?
    echo "--- long mode $LMODE: exit code $RC" 2>&1 | tee -a "$LOG"
    grep -E "# stress (done|WARNING)" $SLOG | head -5 2>&1 | tee -a "$LOG"
    grep -E "# stress heartbeat" $SLOG | tail -1 2>&1 | tee -a "$LOG"
fi
echo "DONE ($(date))" 2>&1 | tee -a "$LOG"

#!/bin/bash
# run_meson_momproj_claude.sh
#
# Production runner for the momentum-projected connected meson
# (baryons_0000_dirac_claude). For each gauge config in CLIST it runs one
# inversion + meson/baryon contraction and writes a per-config HDF5 file
# containing the new "<channel>_t_p<idx>" datasets (17 momenta) that Chunk 2
# combines with the disconnected corr_d.p* files as 2D - C.
#
# This is a CLUSTER script (FLUX launcher, multi-GPU). The local TITAN V cannot
# hold the full 24^3x48 double-precision solve (~18-20 GB); on mi300a a single
# GCD (128 GB HBM) fits it easily, so a small per-node mpi grid is plenty.
#
# >>> EDIT THE SITE-SPECIFIC SETTINGS below (marked TODO) before running. <<<
# No file is ever deleted/overwritten by this script; outputs go to OUTDIR and
# each config writes its own h5. Re-running for a config that already has an h5
# will have the binary overwrite that one h5 (the program's own behaviour).

#FLUX: -t 240m
#FLUX: -q pbatch
#FLUX: -N 1
#FLUX: -n 4
#FLUX: -g 1
#FLUX: --exclusive

set -eu

# ---------------- site-specific settings (TODO: adjust) ----------------
# Cluster environment that sets compiler/runtime/HDF5 paths (as in run.sh).
ENV_SH="/usr/workspace/lsd/matsumoto5/su4_32c/env.sh"   # TODO: confirm path

# Cluster-built binary (must be rebuilt on the cluster for its accelerator;
# the local ./bin build is CUDA/TITAN V only).
BIN="./bin/baryons_0000_dirac_claude"                   # TODO: cluster build path

# Ensemble parameters (this dir: obs_nc4nf1_2448_b11p045_m0p4000).
mass=0.4000
M5=1.5
# Gauge configuration directory and file-name stem. The Nersc config file for
# index i is "${CFGPATH}/${cfgfilename}_lat.${i}" (as in run.sh / run_old.sh).
cfgfilename="conf_nc4nf1_2448_b11p045_m0p4000"          # TODO: confirm stem
CFGPATH="/p/lustre5/park49/SU4_sdm/gauge_conf/${cfgfilename}"  # TODO: confirm

# Config indices to process (matches the disc tracelist: 1500..1580 step 20).
CLIST="1500 1520 1540 1560 1580"                        # TODO: confirm list

# Output directory for per-config connected h5 files.
OUTDIR="./meson_momproj_out_claude"

# Lattice / decomposition. 24^3x48 fits a single mi300a GCD; 4 GPUs/node here.
LATT="24.24.24.48"
MPIGRID="2.2.1.1"
# -----------------------------------------------------------------------

source "${ENV_SH}"
export HDF5_USE_FILE_LOCKING=FALSE

OPTIONS="--decomposition --comms-concurrent --comms-overlap --debug-mem --shm 2048 --shm-mpi 1"
PARAMS_GRID="--grid ${LATT} --mpi ${MPIGRID} --threads 8 --accelerator-threads 8 ${OPTIONS}"

mkdir -p "${OUTDIR}"

if [ ! -x "${BIN}" ]; then
    echo "ERROR: binary not found/executable: ${BIN}" >&2
    echo "       Build baryons_0000_dirac_claude.cc on the cluster first." >&2
    exit 1
fi

echo "--start " "$(date)" "$(date +%s)"

for i in ${CLIST}; do

    cfg="${CFGPATH}/${cfgfilename}_lat.${i}"
    outfile="${OUTDIR}/meson_momproj_${cfgfilename}_lat.${i}.h5"
    log="${OUTDIR}/run_${i}.log"

    if [ ! -f "${cfg}" ]; then
        echo "WARNING: config not found, skipping: ${cfg}" | tee -a "${log}"
        continue
    fi

    echo "==== config ${i} -> ${outfile} ====" | tee "${log}"
    # binary positional args: <config> <M5> <mass> <outfile>
    flux run -N 1 --tasks-per-node=4 --verbose --exclusive \
        --setopt=mpibind=verbose:1 \
        "${BIN}" "${cfg}" "${M5}" "${mass}" "${outfile}" ${PARAMS_GRID} \
        2>&1 | tee -a "${log}"

done

echo "--end " "$(date)" "$(date +%s)"

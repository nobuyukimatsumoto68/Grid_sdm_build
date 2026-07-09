#!/bin/bash
# run_two_baryon_flux_claude.sh
#
# FLUX batch launcher (tuolumne / LLNL, MI300A) for the two-baryon correlator
# two_baryon_corr_claude, on the INTERACTING ensemble b11.045 m0.4 (24^3 x 48).
# Per config: two point-to-all solves (sources O=(0,0,0), M=(12,12,12)) + the
# 8x8-determinant contraction, writing a per-config HDF5 with datasets:
#   bar_0000_t          - single baryon, momentum-projected  (gives M_B)
#   two_baryon_0000_t   - two baryons at fixed points O',M'   (gives E_2B)
#
# Conventions follow the recent measurement scripts run_meson_momproj_claude.sh
# (same <config> <M5> <mass> <outfile> interface, -N 1/4 GPUs, 24^3x48 fits one
# MI300A GCD) and submit_disc_lma_prod_claude.sh (env / paths / flux options).
#
# Submit yourself:  flux batch run_two_baryon_flux_claude.sh   (Claude does not submit.)
# No file is deleted/overwritten by this script; each config writes its own h5.

#FLUX: -t 120m
#FLUX: --output=two_baryon_{{id}}.out
#FLUX: -q pbatch
#FLUX: -N 1
#FLUX: -n 4
#FLUX: -g 1
#FLUX: --exclusive

set -u

date; hostname
export FASTLOAD_VERBOSE=1
export SPINDLE_FLUXOPT=off
export HDF5_USE_FILE_LOCKING=FALSE

# ---------------- environment / binary ----------------
ROOT=/usr/workspace/lsd/matsumoto5/su4_32c
source ${ROOT}/env.sh
# Cluster-built binary (build ON tuolumne against the Grid build at ${ROOT}/build):
#   cd ${ROOT}/Grid_sdm_build && GRID=${ROOT}/build ./compile_two_baryon_claude.sh
APP=${ROOT}/Grid_sdm_build/src/gauge_gen_Nc4/bin/two_baryon_corr_claude

# ---------------- ensemble (interacting: b11.045, m0.4) ----------------
mass=0.4000
M5=1.5
betastr=11p045
massstr=0p4000
cfgfilename=conf_nc4nf1_2448_b${betastr}_m${massstr}
CFGPATH=/p/lustre5/matsumoto5/conf_nc4nf1_2448/${cfgfilename}

# Config indices to process. NERSC file for index i is ${CFGPATH}/${cfgfilename}_lat.${i}.
# Available i ~ 1000..9920 (step 20). A handful of well-separated (decorrelated)
# thermalized configs for the check; edit as desired.
CLIST=${CLIST:-"5000 5500 6000 6500 7000"}

# Output directory for per-config two-baryon h5 files (lustre obs dir, matching
# the disc/meson obs convention: one <obs>.<conf>.h5 per config).
OUTDIR=${OUTDIR:-/p/lustre5/matsumoto5/obs_nc4nf1_2448/obs_nc4nf1_2448_b${betastr}_m${massstr}}

# ---------------- lattice / GPU decomposition ----------------
LATT="24.24.24.48"
MPIGRID="2.2.1.1"        # 4 ranks -> 4 GCDs on one node; 24^3x48 fits one GCD

OPTIONS="--decomposition --comms-concurrent --comms-overlap --debug-mem --shm 2048 --shm-mpi 1"
PARAMS_GRID="--grid ${LATT} --mpi ${MPIGRID} --threads 8 --accelerator-threads 8 ${OPTIONS}"

mkdir -p "${OUTDIR}"

if [ ! -x "${APP}" ]; then
    echo "ERROR: binary not found/executable: ${APP}" >&2
    echo "       Build on tuolumne first:" >&2
    echo "       cd ${ROOT}/Grid_sdm_build && GRID=${ROOT}/build ./compile_two_baryon_claude.sh" >&2
    exit 1
fi

echo "--start " "$(date)" "$(date +%s)"
echo "ensemble = ${cfgfilename}  (M5=${M5}, mass=${mass})"
echo "cfgdir   = ${CFGPATH}"
echo "outdir   = ${OUTDIR}"
echo "configs  = ${CLIST}"

for i in ${CLIST}; do

    cfg="${CFGPATH}/${cfgfilename}_lat.${i}"
    outfile="${OUTDIR}/two_baryon.${i}.h5"
    log="${OUTDIR}/two_baryon.${i}.log"

    if [ ! -f "${cfg}" ]; then
        echo "WARNING: config not found, skipping: ${cfg}" | tee -a "${log}"
        continue
    fi

    echo "==== config ${i} -> ${outfile} ====" | tee "${log}"
    # binary positional args: <config> <M5> <mass> <outfile>
    # Watch the SelfTestTwoBaryon() print at startup: rel diff should be ~1e-12.
    flux run -N 1 --tasks-per-node=4 --verbose --exclusive \
        --setopt=mpibind=verbose:1 \
        "${APP}" "${cfg}" "${M5}" "${mass}" "${outfile}" ${PARAMS_GRID} \
        2>&1 | tee -a "${log}"

done

echo "--end " "$(date)" "$(date +%s)"

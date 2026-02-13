#!/bin/bash

export HDF5_USE_FILE_LOCKING=FALSE
# run if there's any output missing


### Launch parallel executable
# GRID_DIR=/usr/WS2/lsd/sungwoo/SU4_sdm/Grid_sdm_build/mi300a_test_6.2.0

OPTIONS="--decomposition --comms-concurrent --comms-overlap --debug-mem  --shm 2048 --shm-mpi 1"

PARAMS_GRID=" --grid 24.24.24.48 --mpi 2.1.1.2 --threads 8 --accelerator-threads 8 ${OPTIONS}"


#cfg=${CFGPATH}/${cfgfilename}_lat.$i
cfg=conf_nc4nf1_2448_b11p035_m0p4000_lat.10000
# LOG=./log/fixlog.${cfgfilename}_lat.$i
# date > ${LOG}
# tmp=./log/tmp.$i

#--------
# 1. pbp
#APP="$GRID_DIR/build/build_gauge_gen_unified_Nc4/bin/pbp"
# APP=../bin/pbp_test021125
# APP=../../build_gauge_gen_unified_Nc4/bin/pbp
outfile="./outfile.h5"
PARAMS="${cfg} $outfile"
echo $PARAMS

#flux run -N 1 --tasks-per-node=4 --verbose --exclusive --setopt=mpibind=verbose:1 $APP ${PARAMS} $PARAMS_GRID
flux run -N 1 --tasks-per-node=4 --verbose --exclusive --setopt=mpibind=verbose:1 ./glueball ${PARAMS} $PARAMS_GRID

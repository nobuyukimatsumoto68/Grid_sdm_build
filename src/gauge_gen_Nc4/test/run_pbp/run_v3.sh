#!/bin/bash

# export HDF5_USE_FILE_LOCKING=FALSE
# # run if there's any output missing

# # gauge configuration index
# #i=$1
# i=5000


# # here=`pwd`
# # cfgfilename=${here##*/}
# # label=${here##*conf_nc4nf1_}

# cfgfilename=conf_nc4nf1_248_b10p80_m0p1000
# label=${cfgfilename##*conf_nc4nf1_}
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:/opt/cray/pe/lib64/:${LD_LIBRARY_PATH}
echo $LD_LIBRARY_PATH

CFGPATH="/p/lustre5/park49/SU4_sdm/gauge_conf/conf_nc4nf1_248_b10p90c_m0p2000"
ID="conf_nc4nf1_248_b10p90c_m0p2000_lat."
outfilebasename="./conf_nc4nf1_248_b10p90c_m0p2000"
beta="10.9"

mass=0.2
#M5=1.8
M5=1.5

# ### Launch parallel executable
# # GRID_DIR=/usr/WS2/lsd/sungwoo/SU4_sdm/Grid_sdm_build/mi300a_test_6.2.0

OPTIONS="--decomposition --comms-concurrent --comms-overlap --debug-mem  --shm 2048 --shm-mpi 1"

PARAMS_GRID=" --grid 24.24.24.8 --mpi 2.2.1.1 --threads 8 --accelerator-threads 8 ${OPTIONS}"

#export HDF5_USE_FILE_LOCKING=FALSE

cfg=${CFGPATH}/${cfgfilename}_lat.$i
# LOG=./log/fixlog.${cfgfilename}_lat.$i
# date > ${LOG}
# tmp=./log/tmp.$i

#--------
# 1. pbp
#APP="$GRID_DIR/build/build_gauge_gen_unified_Nc4/bin/pbp"
# APP=../bin/pbp_test021125
# APP=../../build_gauge_gen_unified_Nc4/bin/pbp
# outfile="./outfile.h5"
# PARAMS="${cfg} $M5 $mass $outfile"
i=1000
f=1100
interval=4
PARAMS="$CFGPATH $ID $outfilebasename $M5 $mass $i $f $interval"
echo $PARAMS

#flux run -N 1 --tasks-per-node=4 --verbose --exclusive --setopt=mpibind=verbose:1 $APP ${PARAMS} $PARAMS_GRID
flux run -N 1 -q pdebug --tasks-per-node=4 --verbose --exclusive --setopt=mpibind=verbose:1 ./pbp_list ${PARAMS} $PARAMS_GRID

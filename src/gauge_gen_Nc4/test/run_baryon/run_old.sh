#!/bin/bash
OLD_DIR=/g/g20/park49/ws_lsd/sungwoo/SU4_sdm/Grid_sdm_build/mi300a_test_6.2.0/build/test_021125/run_baryon
source /g/g20/park49/ws_lsd/sungwoo/SU4_sdm/Grid_sdm_build/mi300a_test_6.2.0/env.sh
export HDF5_USE_FILE_LOCKING=FALSE
# run if there's any output missing

# gauge configuration index
#i=$1
i=1000


# here=`pwd`
# cfgfilename=${here##*/}
# label=${here##*conf_nc4nf1_}

cfgfilename=conf_nc4nf1_248_b10p90c_m0p2000
label=${cfgfilename##*conf_nc4nf1_}

CFGPATH="/p/lustre5/park49/SU4_sdm/gauge_conf/${cfgfilename}"

#export LD_LIBRARY_PATH=/usr/WS2/lsd/sungwoo/SU4_sdm/Grid_sdm_build/mi300a_rocm6.4.1_hdf5/install/hdf5/lib:${LD_LIBRARY_PATH}

mass=0.${label##*m0p}
#M5=1.8
M5=1.5

### Launch parallel executable
# GRID_DIR=/usr/WS2/lsd/sungwoo/SU4_sdm/Grid_sdm_build/mi300a_test_6.2.0

OPTIONS="--decomposition --comms-concurrent --comms-overlap --debug-mem  --shm 2048 --shm-mpi 1"

#PARAMS_GRID=" --grid 24.24.24.8 --mpi 2.2.1.1 --threads 8 --accelerator-threads 8 ${OPTIONS}"
PARAMS_GRID=" --grid 24.24.24.8 --mpi 2.2.1.1 --threads 8 --accelerator-threads 8 ${OPTIONS}"


cfg=${CFGPATH}/${cfgfilename}_lat.$i
# LOG=./log/fixlog.${cfgfilename}_lat.$i
# date > ${LOG}
# tmp=./log/tmp.$i

#--------
# 1. pbp
#APP="$GRID_DIR/build/build_gauge_gen_unified_Nc4/bin/pbp"
# APP=../bin/pbp_test021125
# APP=../../build_gauge_gen_unified_Nc4/bin/pbp
outfile="./outfile_oldrun.h5"
PARAMS="${cfg} $M5 $mass $outfile"
echo $PARAMS

#flux run -N 1 --tasks-per-node=4 --verbose --exclusive --setopt=mpibind=verbose:1 $APP ${PARAMS} $PARAMS_GRID
#flux run -N 1 --tasks-per-node=2 --verbose --exclusive --setopt=mpibind=verbose:1 ./baryons ${PARAMS} $PARAMS_GRID
flux run -N 1 --tasks-per-node=4 --verbose --exclusive --setopt=mpibind=verbose:1 ${OLD_DIR}/baryons ${PARAMS} $PARAMS_GRID

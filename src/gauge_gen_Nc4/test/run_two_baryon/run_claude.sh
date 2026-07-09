#!/bin/bash
# run_claude.sh - two-baryon correlator C_2B(t) on one configuration.
#
# Mirrors test/run_baryon/run.sh. With no PARAMS the binary uses a cold
# configuration (M5=1.5, mass=0.1, outfile=ColdConfig.h5); set PARAMS to
# "<config> <M5> <mass> <outfile>" to run on a real NERSC gauge field.

#export HDF5_USE_FILE_LOCKING=FALSE

i=1000
cfgfilename=conf_nc4nf1_248_b10p90c_m0p2000
label=${cfgfilename##*conf_nc4nf1_}
CFGPATH="/p/lustre5/park49/SU4_sdm/gauge_conf/${cfgfilename}"

mass=0.${label##*m0p}
M5=1.5

OPTIONS="--decomposition --comms-concurrent --comms-overlap --debug-mem --shm 2048 --shm-mpi 1"
# Production geometry for the 24^3 x 48 ensembles:
#PARAMS_GRID="--grid 24.24.24.48 --mpi 1.1.1.1 --threads 8 --accelerator-threads 8 ${OPTIONS}"
PARAMS_GRID="--grid 8.8.8.8 --mpi 1.1.1.1 --threads 8 --accelerator-threads 8 ${OPTIONS}"

cfg=${CFGPATH}/${cfgfilename}_lat.$i
outfile="./two_baryon.h5"
#PARAMS="${cfg} $M5 $mass $outfile"
PARAMS=""   # empty -> cold configuration smoke test
echo "PARAMS=${PARAMS}"

../../bin/two_baryon_corr_claude ${PARAMS} $PARAMS_GRID

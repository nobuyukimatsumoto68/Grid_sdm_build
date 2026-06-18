#!/bin/bash
#FLUX: -t 120m
#FLUX: --output=hmc_tune_{{id}}
#FLUX: -q pbatch
#FLUX: -N 8
#FLUX: -n 32
#FLUX: -g 1
#FLUX: --exclusive
#
# FLUX batch launcher for one HMC stream of the M5/Ls tuning ensemble.
#   32^3 x 64, m = 0.01, Ls = 16, M5 = 1.5, b/c = 1.5/0.5 (Mobius EOFA, SU(4) SDM).
# One gauge_beta per job; run_tuo_tune_claude.sh seds @BETA@ into the XML and
# submits one of these per stream (beta = 10.7, 10.8, 10.9).
# Modeled on ../../../../submit_hmc_tuolumne.sh.


### four gpus per node
date
here=`pwd`

# To verify fastload2 is loading or not, set env variable FASTLOAD_VERBOSE=1 and rank 0
# (or serial tasks) will print out something if fastload2 is being used (email 10/30/24).
export FASTLOAD_VERBOSE=1

# email notice on 12/11/24
export SPINDLE_FLUXOPT=off


XML=ip_hmc_mobius_claude.xml


### Launch parallel executable
source /usr/workspace/lsd/matsumoto5/su4_32c/env.sh
APP=/usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build/src/gauge_gen_Nc4/bin/dweofa_mobius_HSDM_v4

echo "--start " `date` `date +%s`

OPTIONS="--decomposition --comms-concurrent --comms-overlap --debug-mem  --shm 2048 --shm-mpi 1"
PARAMS=" --grid 32.32.32.64 --mpi 2.2.2.4 --threads 8 --accelerator-threads 8 ${OPTIONS} --ParameterFile ${XML}"

flux run -N 8 --tasks-per-node=4 --verbose --exclusive --setopt=mpibind=verbose:1 $APP $PARAMS

echo "--end " `date` `date +%s`

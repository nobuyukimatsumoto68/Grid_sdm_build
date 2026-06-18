#!/bin/bash
# set -e
#
# Driver: thermalize HMC streams for the M5/Ls tuning ensemble.
#   32^3 x 64, Ls = 16, M5 = 1.5, b/c = 1.5/0.5 (Mobius EOFA, SU(4) SDM).
#   Scan two quark masses (0.01, 0.05) x three gauge_beta (10.7, 10.8, 10.9):
#   6 streams total, each its own flux batch job (scheduled independently / parallel).
#
# Modeled on ../../../../run_tuo.sh + submit_hmc_tuolumne.sh. Resumable thermalization:
# the FIRST submission of a stream HotStarts (no ckpoint present); every later
# submission auto-detects the latest ckpoint_lat.<traj> in the stream dir and continues
# with StartingType=CheckpointStart, StartTrajectory=<latest>, NoMetropolisUntil=0.
# So just re-run this driver to keep thermalizing toward a plaquette plateau; at the
# 120m wall (~11 traj/job at ~646 s/traj) several resubmissions are expected.
#
# All physics parameters except mass, gauge_beta and the trajectory bookkeeping live in
# the XML template; this driver seds those. NTRAJ is the TOTAL trajectory target (upper
# limit of the HMC loop); raise it to extend thermalization.

masses=(0.01 0.05)
betas=(10.7 10.8 10.9)
NTRAJ=200                 # total trajectory target (Trajectories in the XML)

basedir=$(pwd)
Nt=64
xml=ip_hmc_mobius_claude.xml
script=submit_hmc_tuolumne_claude.sh

jmax=${#masses[@]}
for((j=0;j<$jmax;j++))
do
    m=${masses[$j]}
    rundir=/p/lustre5/matsumoto5/tune_32_64_m${m}
    mkdir -p ${rundir}

    imax=${#betas[@]}
    echo $imax
    for((i=0;i<$imax;i++))
    do
        beta=${betas[$i]}
        dir=${rundir}/beta${beta}m${m}
        mkdir -p $dir
        # copy the pristine template (placeholders intact) into each subdir
        cp -f ${basedir}/$xml $dir
        cp -f ${basedir}/$script $dir
        echo $dir
        cd $dir
        sed -i "/<gauge_beta>/{s/@BETA@/${beta}/}" $xml
        sed -i "/<mass>/{s/0.01/${m}/}" $xml
        sed -i "/<Trajectories>/{s/40/${NTRAJ}/}" $xml

        # resume: if a checkpoint already exists, continue from the latest trajectory
        st=$(ls ckpoint_lat.* 2>/dev/null | sed 's/.*ckpoint_lat.//' | sort --version-sort | tail -n1)
        if [ -n "$st" ]; then
            sed -i "/StartTrajectory/{s/0/${st}/}" $xml
            sed -i "/StartingType/{s/HotStart/CheckpointStart/}" $xml
            sed -i "/NoMetropolisUntil/{s/5/0/}" $xml
            echo "  resume from ckpoint_lat.${st} (CheckpointStart)"
        else
            echo "  fresh HotStart"
        fi

        flux batch $script # mybatchscript
        cd ${basedir}
    done
done

#!/bin/bash
#
# Driver: PRODUCTION HMC streams for the 32^3x64 SU(4) SDM ensemble at the chosen
# production beta(s). Modeled on run_tuo_tune_claude.sh and RESUMABLE the same way
# (auto CheckpointStart from the latest ckpoint; first submission HotStarts).
#
# Generates a few configs for the mres check: saveInterval=20, NTRAJ=100 -> saves at
# traj 20,40,60,80,100 = 5 configs per stream. Output is a SEPARATE tree from the
# tuning streams: /p/lustre5/matsumoto5/prod_32_64_m<m>/beta<beta>m<m>/.
# Reuses ip_hmc_mobius_claude.xml (@BETA@, m/Ls/M5/b/c) + submit_hmc_tuolumne_claude.sh.

masses=(0.01 0.05)
betas=(BETA_TODO)        # <<< FILL IN the production beta value(s), e.g. betas=(10.85)
NTRAJ=100                # total trajectory target -> 5 saved configs at saveInterval=20
SAVEINT=20               # saveInterval (save a config every 20 trajectories)

# guard: refuse to run with the placeholder still in place
if [[ " ${betas[*]} " == *" BETA_TODO "* ]]; then
    echo "ERROR: set the production beta(s) in 'betas=(...)' before running." >&2
    exit 1
fi

basedir=$(pwd)
Nt=64
xml=ip_hmc_mobius_claude.xml
script=submit_hmc_tuolumne_claude.sh

jmax=${#masses[@]}
for((j=0;j<$jmax;j++))
do
    m=${masses[$j]}
    rundir=/p/lustre5/matsumoto5/prod_32_64_m${m}
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
        sed -i "/<saveInterval>/{s/5/${SAVEINT}/}" $xml

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

#!/bin/bash
#
# Driver: PRODUCTION HMC streams for the 32^3x64 SU(4) SDM ensemble.
# Two INDEX-ALIGNED streams (beta paired with mass), mimicking the 24c directory/file
# naming under /p/lustre5/matsumoto5 (conf_nc4nf1_3264 = the 32^3x64 analog of the
# 24^3x48 conf_nc4nf1_2448 tree):
#   /p/lustre5/matsumoto5/conf_nc4nf1_3264/conf_nc4nf1_3264_b<betastr>_m<massstr>/
#     conf_nc4nf1_3264_b<betastr>_m<massstr>_lat.<traj>   (NERSC config, IEEE64BIG)
#     conf_nc4nf1_3264_b<betastr>_m<massstr>_rng.<traj>   (RNG state)
#
# RESUMABLE (auto CheckpointStart from the latest saved config; first submission HotStarts).
# saveInterval=4, NTRAJ=100 -> configs at traj 4,8,...,100. Submitted with -t 360m
# (overrides the batch script's 120m directive) so whole saveInterval blocks fit the wall;
# WALL_SECONDS is passed so the blocker's fallback matches the 360m allocation.
# Reuses ip_hmc_mobius_claude.xml + submit_hmc_tuolumne_claude.sh (graceful wall blocker).

# index-aligned arrays: stream j uses masses[j]/betas[j] and the matching name strings.
masses=(0.01     0.05)
betas=(10.8      10.84)
massstrs=(0p0100 0p0500)
betastrs=(10p800 10p840)

NTRAJ=100                # total trajectory target
SAVEINT=4                # saveInterval (save a config every 4 trajectories)
WALL=360m                # flux batch wall override for production
WALL_SECONDS=21600       # same wall in seconds, for the blocker fallback (360m)

basedir=$(pwd)
Nt=64
xml=ip_hmc_mobius_claude.xml
script=submit_hmc_tuolumne_claude.sh
topdir=/p/lustre5/matsumoto5/conf_nc4nf1_3264

mkdir -p ${topdir}

jmax=${#masses[@]}
for((j=0;j<$jmax;j++))
do
    m=${masses[$j]}
    beta=${betas[$j]}
    massstr=${massstrs[$j]}
    betastr=${betastrs[$j]}

    cfgname=conf_nc4nf1_3264_b${betastr}_m${massstr}
    dir=${topdir}/${cfgname}
    mkdir -p $dir

    # copy the pristine template (placeholders intact) into the stream dir
    cp -f ${basedir}/$xml $dir
    cp -f ${basedir}/$script $dir
    echo $dir
    cd $dir

    sed -i "/<gauge_beta>/{s/@BETA@/${beta}/}" $xml
    sed -i "/<mass>/{s/0.01/${m}/}" $xml
    sed -i "/<Trajectories>/{s/40/${NTRAJ}/}" $xml
    sed -i "/<saveInterval>/{s/5/${SAVEINT}/}" $xml
    # config/rng file names mimic the 24c pattern: <cfgname>_lat.<traj>, <cfgname>_rng.<traj>
    sed -i "/<config_prefix>/{s|./ckpoint_lat|./${cfgname}_lat|}" $xml
    sed -i "/<rng_prefix>/{s|./ckpoint_rng|./${cfgname}_rng|}" $xml

    # resume: continue from the latest saved config if present
    st=$(ls ${cfgname}_lat.* 2>/dev/null | sed "s/.*_lat.//" | sort --version-sort | tail -n1)
    if [ -n "$st" ]; then
        sed -i "/StartTrajectory/{s/0/${st}/}" $xml
        sed -i "/StartingType/{s/HotStart/CheckpointStart/}" $xml
        sed -i "/NoMetropolisUntil/{s/5/0/}" $xml
        echo "  resume from ${cfgname}_lat.${st} (CheckpointStart)"
    else
        echo "  fresh HotStart"
    fi

    flux batch -t ${WALL} --env=WALL_SECONDS=${WALL_SECONDS} $script
    cd ${basedir}
done

#!/bin/bash
source env.sh

pushd ${BUILDDIR}

HMC=gauge_gen_unified_Nc4
if [ -d ./build_${HMC} ];
then
    rm -rf ./build_${HMC}
fi
if [ -d ../install/${HMC} ];
then
    rm -rf ../install/${HMC}
fi
if [ ! -d ./build_grid_omp_unified_Nc4 ];
then
    echo "Need build_grid_omp_unified_Nc4"
    exit 1
fi

cp -a ../src/gauge_gen_Nc4 build_${HMC}
cd ./build_${HMC}

sed -i 's/_Nc4/_unified_Nc4/g' Makefile
make
make install
popd

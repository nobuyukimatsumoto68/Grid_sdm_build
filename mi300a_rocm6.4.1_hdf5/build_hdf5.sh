#!/bin/bash

source ./env.sh

pushd ${BUILDDIR}

HDF5=hdf5
if [ -d ./build_${HDF5} ];
then
  rm -rf ./build_${HDF5}
fi
if [ -d ../install/${HDF5} ];
then
  rm -rf ../install/${HDF5}
fi

mkdir  ./build_${HDF5}
cd ./build_${HDF5}

cmake ${SRCDIR}/${HDF5} -G "Unix Makefiles"\
       -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
       -DCMAKE_C_COMPILER=/usr/bin/gcc  \
       -DCMAKE_BUILD_TYPE:STRING=Release \
       -DBUILD_SHARED_LIBS:BOOL=ON \
       -DBUILD_TESTING:BOOL=OFF \
       -DHDF5_BUILD_TOOLS:BOOL=OFF \
       -DHDF5_BUILD_CPP_LIB:BOOL=ON \
       -DHDF5_ENABLE_PARALLEL:BOOL=OFF \
       -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/${HDF5}

cmake --build . -j 32 -v
cmake --install .

popd

       # -DCMAKE_CXX_FLAGS="${MPI_CFLAGS}" \
       # -DCMAKE_C_FLAGS="${MPI_CFLAGS}" \
       # -DCMAKE_SHARED_LINKER_FLAGS="${MPI_LDFLAGS}" \
       # -DCMAKE_EXE_LINKER_FLAGS="${MPI_LDFLAGS}" \

       # -DCMAKE_CXX_STANDARD=11 \
       # -DCMAKE_CXX_FLAGS="${MPI_CFLAGS} -std=c++11" \

# Note that following happens if both CPP_LIB and PARALLEL turned on 
# CMake Error at CMakeLists.txt:1099 (message):
#    **** Parallel and C++ options are mutually exclusive, override with ALLOW_UNSUPPORTED option ****

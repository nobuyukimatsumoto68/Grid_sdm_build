#GRID=../build_grid_omp_Nc4
GRID=../../../build/ # ../../install/Grid_omp_Nc4
CONFIG=$(GRID)/bin/grid-config
CXX=hipcc -std=c++17 # $(shell $(CONFIG) --cxx)
CXXFLAGS=-I/usr/workspace/lsd/matsumoto5/opt/lime-1.3.2//include -O3 -fPIC -I/opt/cray/pe/mpich/9.0.1/ofi/crayclang/20.0/include -I/opt/rocm-6.4.1/include -L/lib64 -fgpu-sanitize --offload-arch=gfx942 -fopenmp -fno-strict-aliasing -I/usr/workspace/lsd/matsumoto5/su4_32c/build/include # $(shell $(CONFIG) --cxxflags)
LDFLAGS=-L/usr/workspace/lsd/matsumoto5/opt/lime-1.3.2//lib -L/lib64 -lamdhip64 -lhipblas -lrocblas -Wl,-rpath=/opt/cray/pe/mpich/9.0.1/ofi/crayclang/20.0/lib -L/opt/cray/pe/mpich/9.0.1/ofi/crayclang/20.0/lib -lmpi -Wl,-rpath=/opt/cray/pe/mpich/9.0.1/gtl/lib/ -L/opt/cray/pe/mpich/9.0.1/gtl/lib/ -L/opt/rocm-6.4.1/llvm/lib -Wl,-rpath=/opt/rocm-6.4.1/llvm/lib -L/opt/rocm-6.4.1/lib -lmpi -L/opt/cray/pe/mpich/9.0.1/gtl/lib -lmpi_gtl_hsa -L/usr/workspace/lsd/matsumoto5/su4_32c/build/lib # $(shell $(CONFIG) --ldflags)
LIBS=-lhdf5_cpp -lz -lcrypto -lunwind-x86_64 -llime -lfftw3f -lfftw3 -lmpfr -lgmp -lhipblas -lrocblas -lhipfft -lz -lhdf5 -lGrid # $(shell $(CONFIG) --libs)

BIN_DIR=./bin
INSTALL_DIR=../../install/gauge_gen_Nc4

SRC=$(wildcard *.cc)
BIN=$(basename $(SRC))

#.PRECIOUS: hmc_SDM.o
.PHONY: clean

all: $(BIN)

%: %.cc
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $< -o $@
	mv $@ $(BIN_DIR)

install:
	mkdir -p $(INSTALL_DIR)
	mv $(BIN_DIR) $(INSTALL_DIR)/.

clean:
	$(RM) -r $(BIN_DIR)

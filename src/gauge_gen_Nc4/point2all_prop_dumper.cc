#include <Grid/Grid.h>
#include <algorithm>
#include <vector>
#include <random>

using namespace std;
using namespace Grid;


int main(int argc, char ** argv)
{
  const int Ls = 16;

  Grid_init(&argc, &argv);

  // Double precision grids
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                       GridDefaultSimd(Nd, vComplex::Nsimd()),
                                       GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  LatticeGaugeField Umu(UGrid);
  std::string config;
  std::string outfile;
  RealD M5, mass;
  if( argc > 1 && argv[1][0] != '-' )
  {
    config  = argv[1];
    M5      = stod(argv[2]);
    mass    = stod(argv[3]);
    outfile = argv[4];
    std::cout << GridLogMessage << "Loading configuration from " << config << std::endl;
    std::cout << GridLogMessage << "M5="   << M5   << std::endl;
    std::cout << GridLogMessage << "mass=" << mass << std::endl;
    std::cout << GridLogMessage << "output: " << outfile << std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, config);
  }
  else
  {
    std::cout << GridLogMessage << "Using cold configuration" << std::endl;
    SU<Nc>::ColdConfiguration(Umu);
    config  = "ColdConfig";
    M5      = 1.5;
    mass    = 0.1;
    outfile = config + ".h5";
  }

  // Mobius domain-wall action as scaled Shamir kernel, anti-periodic in time.
  RealD b = 1.5; // b+c = 2, b-c = 1
  RealD c = 0.5;
  std::cout << GridLogMessage << "==============================================" << std::endl;
  std::cout << GridLogMessage << "MobiusFermion (scaled Shamir kernel), AP BC {1,1,1,-1}" << std::endl;
  std::cout << GridLogMessage << "==============================================" << std::endl;
  std::vector<Complex> boundary = {1,1,1,-1};
  MobiusFermionD::ImplParams Params(boundary);
  MobiusFermionD D(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, c, Params);


  // ------------------------------------------------------------------
  // point source on the source timeslice t_src = 0:
  // ------------------------------------------------------------------
  Coordinate latt = GridDefaultLatt();

  for(int i=0; i<std::pow(2,3); i++){
    Coordinate pt;
    int j=i;
    for(int mu=0; mu<3; mu++){
      const int k=j%2;
      j/=2;
      pt[mu] = k * (latt[mu]/2);
    }
    pt[3] = 0;

    std::cout << GridLogMessage << "Sources " << i << " = "
              << pt[0] << "," << pt[1] << "," << pt[2] << "," << pt[3] << std::endl;

    LatticePropagator src(UGrid);
    PointSource(pt, src);

    LatticePropagator Prop(UGrid);
    std::cout << GridLogMessage << "Solving for source ..." << std::endl;
    Solve(D, src, Prop);
  }

  Grid_finalize();
}

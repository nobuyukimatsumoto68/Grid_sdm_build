#include <Grid/Grid.h>
#include <filesystem>

using namespace std;
using namespace Grid;


void PointSource(Coordinate &coor,LatticePropagator &source)
{
  source=Zero();
  SpinColourMatrix kronecker; kronecker=1.0;
  pokeSite(kronecker, source, coor);
}

void StochasticSource(GridParallelRNG &RNG, LatticePropagator &source)
{
  GridBase *grid = source.Grid();
  LatticeComplex noise(grid);
  bernoulli(RNG, noise); // 0,1 50:50 in cplx

  RealD nrm = 1.0/sqrt(2.0);
  noise = ( 2.0*noise - Complex(1.0,1.0) )*nrm;

  source = 1.0;
  source = source*noise;
}

template<class Action>
void Solve(Action &D,LatticePropagator &source,LatticePropagator &propagator)
{
  GridBase *UGrid = D.GaugeGrid();
  GridBase *FGrid = D.FermionGrid();

  LatticeFermion src4  (UGrid);
  LatticeFermion src5  (FGrid);
  LatticeFermion result5(FGrid);
  LatticeFermion result4(UGrid);

  ConjugateGradient<LatticeFermion> CG(1.0e-8,100000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion> schur(CG);
  ZeroGuesser<LatticeFermion> ZG; // Could be a DeflatedGuesser if have eigenvectors
  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){
      PropToFerm<Action>(src4,source,s,c);

      D.ImportPhysicalFermionSource(src4,src5);

      result5=Zero();
      schur(D,src5,result5,ZG);
      std::cout<<GridLogMessage
               <<"spin "<<s<<" color "<<c
               <<" norm2(src5d) "   <<norm2(src5)
               <<" norm2(result5d) "<<norm2(result5)<<std::endl;

      D.ExportPhysicalFermionSolution(result5,result4);

      FermToProp<Action>(propagator,result4,s,c);
    }
  }
}

void ChCondPtSrc(std::string file, LatticePropagator &q, Coordinate& coord)
{
  LatticeComplex meson_CF( q.Grid() );

  // Gamma G5(Gamma::Algebra::Gamma5);
  meson_CF = trace(q);
  TComplex ChCond;
  peekSite(ChCond, meson_CF, coord);

  Complex res = TensorRemove(ChCond); // Yes this is ugly, not figured a work around
  std::cout << res << std::endl;

  std::cout << "chcond, " << file << ", " << ChCond << std::endl;
  {
    XmlWriter WR(file);
    write(WR, "MesonFile", res);
  }
}

Complex ChCondStochSrc(LatticePropagator &psi, LatticePropagator &eta)
{
  LatticeComplex meson_CF( eta.Grid() );
  meson_CF = trace(psi*adj(eta));
  auto ChCond = sum( meson_CF );

  LatticeComplex identity_CF( eta.Grid() );
  identity_CF = trace(adj(eta)*eta);
  auto norm = sum( identity_CF );

  auto res = ChCond()()/norm()();
  // std::cout << "chcond, " << file << ", " << res << std::endl;
  // {
  //   XmlWriter WR(file);
  //   write(WR, "MesonFile", res);
  // }
  return res;
}


int main (int argc, char ** argv)
{
  const int Ls=16;
  Grid_init(&argc,&argv);

  // Double precision grids
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                   GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // std::vector<int> seeds4({1,2,3,4});
  // GridParallelRNG  RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  std::string config;
  std::string outfile;
  
  std::string infilebasedir;
  std::string infileid;
  std::string outfilebasename;
  RealD M5, mass;
  int conf_min;
  int conf_max;
  int interval;
  if( argc > 1 && argv[1][0] != '-' )
  {
    infilebasedir=argv[1];
    infileid=argv[2];
    outfilebasename=argv[3];
    // config=argv[1];
    M5=stod(argv[4]);
    mass=stod(argv[5]);
    // outfile=argv[4];
    conf_min=atoi(argv[6]);
    conf_max=atoi(argv[7]);
    interval=atoi(argv[8]);
    
    // std::cout << GridLogMessage << "Loading configuration from " << config << std::endl;
    std::cout << GridLogMessage << "mass=" << mass << std::endl;
    std::cout << GridLogMessage << "M5=" << M5 << std::endl;
    std::cout << GridLogMessage << "output: " << outfile << std::endl;
  }
  else
  {
    // std::cout<<GridLogMessage <<"Using cold configuration"<<std::endl;
    // SU<Nc>::ColdConfiguration(Umu);
    // //    SU<Nc>::HotConfiguration(RNG4,Umu);
    // config="ColdConfig";
    std::cout << GridLogMessage << "./pbp <infilebasedir> <infileid> <outfilebasename> <M5> <mass> <conf_min> <conf_max> <interval> --grid xx.xx.xx.x " << std::endl;
    exit(1);
  }

  const std::string infilebasename=infilebasedir + "/" + infileid;
  std::cout << "infilebasename = " << infilebasename << std::endl;
  std::cout << "outfilebasename = " << outfilebasename << std::endl;

  for(int conf=conf_min; conf<conf_max; conf+=interval){
    std::string pathO = outfilebasename+"chcond"+std::to_string(conf)+".xml";

    FieldMetaData header;
    const std::string config = infilebasename+std::to_string(conf);
    if( !std::filesystem::exists( config ) ) continue;
    NerscIO::readConfiguration(Umu, header, config);
  
    // random seed from the configuration name string
    GridParallelRNG          RNG4(UGrid);  RNG4.SeedUniqueString(config);

    RealD b=1.5;// Scale factor b+c=2, b-c=1
    RealD c=0.5;
    std::vector<Complex> boundary = {1,1,1,-1};
    typedef MobiusFermionD FermionAction;
    FermionAction::ImplParams Params(boundary);
    FermionAction FermAct(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, c, Params);

    std::cout<<GridLogMessage <<"======================"<<std::endl;
    std::cout<<GridLogMessage <<"MobiusFermion action as Scaled Shamir kernel"<<std::endl;
    std::cout<<GridLogMessage <<"======================"<<std::endl;

    // noisy estimator with Z4
    LatticePropagator stochastic_source(UGrid);
    StochasticSource( RNG4, stochastic_source );

    LatticePropagator StochProp(UGrid);
    Solve(FermAct, stochastic_source, StochProp);

    Complex chcond = ChCondStochSrc( StochProp, stochastic_source );
  
    {
      // std::unique_ptr<XmlWrite> WR;
      // if (UGrid->IsBoss()){
      XmlWriter WR( pathO );
      WR.scientificFormat( true ); // @@@ double check
      WR.setPrecision( std::numeric_limits<Real>::digits10 + 1 );
      write(WR, "chcond", chcond );
    }
  }
  Grid_finalize();
}




// {
//   SpinColourMatrix tmp;
//   peekSite(tmp, point_source, Origin);
//   std::cout << tmp << std::endl;
// }

// TComplex tmp;
// Coordinate Origin({0,0,0,0});
// peekSite(tmp, check_CF, Origin);
// std::cout << tmp << std::endl;


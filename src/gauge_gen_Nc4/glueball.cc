/* 
 * from Example_plaquette.cc                                                               
 */

#include <Grid/Grid.h>


using namespace std;
using namespace Grid;



void MakePhase(Coordinate mom,LatticeComplex &phase)
{
  GridBase *grid = phase.Grid();
  auto latt_size = grid->GlobalDimensions();
  ComplexD ci(0.0,1.0);
  phase=Zero();

  LatticeComplex coor(phase.Grid());
  for(int mu=0;mu<Nd;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    phase = phase + (TwoPiL * mom[mu]) * coor;
  }
  phase = exp(phase*ci);
}

class MesonFile: Serializable {
public:
  // GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFile, std::vector<std::vector<Complex> >, data);
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFile, std::vector<Complex>,  data);
};


// This copies what already exists in WilsonLoops.h. The point here is to be pedagogical and explain in
// detail what everything does so we can see how GRID works.
template <class Gimpl> class WLoops : public Gimpl {
public:
    // Gimpl seems to be an arbitrary class. Within this class, it is expected that certain types are
    // already defined, things like Scalar and Field. This macro includes a bunch of #typedefs that
    // implement this equivalence at compile time.
    INHERIT_GIMPL_TYPES(Gimpl);

    // Some example Gimpls can be found in GaugeImplementations.h, at the bottom. These are in turn built
    // out of GaugeImplTypes, which can be found in GaugeImplTypes.h. The GaugeImplTypes contain the base
    // field/vector/link/whatever types. These inherit from iScalar, iVector, and iMatrix objects, which
    // are sort of the building blocks for gerenal math objects. The "i" at the beginning of these names
    // indicates that they should be for internal use only. It seems like these base types have the
    // acceleration, e.g. SIMD or GPU or what-have-you, abstracted away. How you accelerate these things
    // appears to be controlled through a template parameter called vtype.

    // The general math/physics objects, such as a color matrix, are built up by nesting these objects.
    // For instance a general color matrix has two color indices, so it's built up like
    //     iScalar<iScalar<iMatrix<vtype ...
    // where the levels going from the inside out are color, spin, then Lorentz indices. Scalars have
    // no indices, so it's what we use when such an index isn't needed. Lattice objects are made by one
    // higher level of indexing using iVector.

    // These types will be used for U and U_mu objects, respectively.
    typedef typename Gimpl::GaugeLinkField GaugeMat;
    typedef typename Gimpl::GaugeField GaugeLorentz;

    // U_mu_nu(x)
    static void dirPlaquette(GaugeMat &plaq, const std::vector<GaugeMat> &U, const int mu, const int nu) {
        // Calls like CovShiftForward and CovShiftBackward have 3 arguments, and they multiply together
        // the first and last argument. (Second arg gives the shift direction.) The CovShiftIdentityBackward
        // has meanwhile only two arguments; it just returns the shifted (adjoint since backward) link. 
        plaq = Gimpl::CovShiftForward(U[mu],mu,
                   // Means Link*Cshift(field,mu,1), arguments are Link, mu, field in that order.
                   Gimpl::CovShiftForward(U[nu],nu,
                       Gimpl::CovShiftBackward(U[mu],mu,
                           // This means Cshift(adj(Link), mu, -1)
                           Gimpl::CovShiftIdentityBackward(U[nu], nu))));
    }

    // tr U_mu_nu(x)
    static void traceDirPlaquette(ComplexField &plaq, const std::vector<GaugeMat> &U, const int mu, const int nu) {
        // This .Grid() syntax seems to get the pointer to the GridBase. Apparently this is needed as argument
        // to instantiate a Lattice object.
        GaugeMat sp(U[0].Grid());
        dirPlaquette(sp, U, mu, nu);
        plaq = trace(sp);
    }

    // sum_mu_nu tr U_mu_nu(x)
    static void sitePlaquette(ComplexField &Plaq, const std::vector<GaugeMat> &U) {
        ComplexField sitePlaq(U[0].Grid());
        Plaq = Zero();
        // Nd=4 and Nc=3 are set as global constants in QCD.h
        for (int mu = 1; mu < Nd; mu++) {
            for (int nu = 0; nu < mu; nu++) {
                traceDirPlaquette(sitePlaq, U, mu, nu);
                Plaq = Plaq + sitePlaq;
            }
        }
    }

    // sum_i_j tr U_i_j(x)
    static void sitespatialPlaquette(ComplexField &Plaq, const std::vector<GaugeMat> &U) {
        ComplexField sitePlaq(U[0].Grid());
        Plaq = Zero();
	// assuming Tdir = Nd-1
        for (int mu = 1; mu < Nd-1; mu++) {
            for (int nu = 0; nu < mu; nu++) {
                traceDirPlaquette(sitePlaq, U, mu, nu);
                Plaq = Plaq + sitePlaq;
            }
        }
    }

    // sum_mu_nu_x Re tr U_mu_nu(x)
    static RealD sumPlaquette(const GaugeLorentz &Umu) {
        std::vector<GaugeMat> U(Nd, Umu.Grid());
        for (int mu = 0; mu < Nd; mu++) {
            // Umu is a GaugeLorentz object, and as such has a non-trivial Lorentz index. We can
            // access the element in the mu Lorentz index with this PeekIndex syntax.
            U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
        }
        ComplexField Plaq(Umu.Grid());
        sitePlaquette(Plaq, U);
        // I guess this should be the line that sums over all space-time sites.
        auto Tp = sum(Plaq);
        // Until now, we have been working with objects inside the tensor nest. This TensorRemove gets
        // rid of the tensor nest to return whatever is inside.
        auto p  = TensorRemove(Tp);
        return p.real();
    }

    // op(t) = sum_mu_nu_{x, spatial} Re tr U_mu_nu(x)
    static void sumspPlaquette(Hdf5Writer &WR,const GaugeLorentz &Umu) {
        std::vector<GaugeMat> U(Nd, Umu.Grid());
        for (int mu = 0; mu < Nd; mu++) {
            // Umu is a GaugeLorentz object, and as such has a non-trivial Lorentz index. We can
            // access the element in the mu Lorentz index with this PeekIndex syntax.
            U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
        }
        ComplexField Plaq(Umu.Grid());
        sitespatialPlaquette(Plaq, U);

	std::vector<TComplex> Plaq_T;
	sliceSum(Plaq,Plaq_T, Tdir);

	int nt=Plaq_T.size();
	double vol = Umu.Grid()->gSites();
        // The number of orientations. xy, xz, yz
        double faces = 3;

	std::vector<Complex> op(nt);
	for(int t=0;t<nt;t++){
	  op[t] = TensorRemove(Plaq_T[t]) * ((double)nt / vol) / faces / Nc; // Yes this is ugly, not figured a work around
	}
	for(int t=0;t<nt;t++){
	  // print real, but save complex later in the file
	  std::cout << " t "<<t<<" " <<op[t].real()<<std::endl;
	}

	// save into hdf file
	MesonFile MF;
	MF.data=op;

	write(WR,"glueballop_S",MF);

    }


    // op(t) = sum_mu_nu_{x, spatial} Re tr U_mu_nu(x)
  static void sumspWilson(Hdf5Writer &WR,const GaugeLorentz &Umu, const int R1, const int R2, Coordinate &mom) {
        std::vector<GaugeMat> U(Nd, Umu.Grid());
        for (int mu = 0; mu < Nd; mu++) {
            // Umu is a GaugeLorentz object, and as such has a non-trivial Lorentz index. We can
            // access the element in the mu Lorentz index with this PeekIndex syntax.
            U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
        }
        ComplexField Wl(Umu.Grid());
        //sitespatialPlaquette(Plaq, U);
	WilsonLoops<Gimpl>::siteSpatialWilsonLoop(Wl, U, R1, R2);

	LatticeComplex phase(Umu.Grid());
	MakePhase(mom,phase);

        ComplexField Wl_phase(Umu.Grid());
	Wl_phase=Wl*phase;
	
	std::vector<TComplex> Wl_T;
	sliceSum(Wl_phase,Wl_T, Tdir);

	int nt=Wl_T.size();
	double vol = Umu.Grid()->gSites();
        // The number of orientations. xy, xz, yz
        double faces = 3;

	std::vector<Complex> op(nt);
	for(int t=0;t<nt;t++){
	  op[t] = TensorRemove(Wl_T[t]) * ((double)nt / vol) / faces / Nc; // Yes this is ugly, not figured a work around
	}
	for(int t=0;t<nt;t++){
	  // print real, but save complex later in the file
	  std::cout << " t "<<t<<" " <<op[t].real()<<std::endl;
	}

	// save into hdf file
	MesonFile MF;
	MF.data=op;

	
	std::string label="glueballop_S_Wl"+std::to_string(R1)+"x"+std::to_string(R2)+"_"
	  +std::to_string(mom[0])+std::to_string(mom[1])+std::to_string(mom[2]);
	write(WR,label,MF);

    }

    // op(t) = sum_mu_nu_{x, spatial} Re tr U_mu_nu(x)
  static void sumspWilson_nsq(Hdf5Writer &WR,const GaugeLorentz &Umu, const int R1, const int R2, const int nsq) {
        std::vector<GaugeMat> U(Nd, Umu.Grid());
        for (int mu = 0; mu < Nd; mu++) {
            // Umu is a GaugeLorentz object, and as such has a non-trivial Lorentz index. We can
            // access the element in the mu Lorentz index with this PeekIndex syntax.
            U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
        }
        ComplexField Wl(Umu.Grid());
        //sitespatialPlaquette(Plaq, U);
	WilsonLoops<Gimpl>::siteSpatialWilsonLoop(Wl, U, R1, R2);

	LatticeComplex phase(Umu.Grid());

	int nt=Umu.Grid()->_fdimensions[3];
	double vol = Umu.Grid()->gSites();
	// The number of orientations. xy, xz, yz
	double faces = 3;

	std::vector<Complex> op(nt);
	for(int t=0;t<nt;t++){
	  op[t]=0.0;		// initialization
	}

	
	int NMAX=2;
	for(int nx=-NMAX;nx<=NMAX;nx++)
	for(int ny=-NMAX;ny<=NMAX;ny++)
	for(int nz=-NMAX;nz<=NMAX;nz++){

	  if (nx*nx+ny*ny+nz*nz!=nsq) continue;
	  
	  Coordinate mom({nx,ny,nz,0});
	  MakePhase(mom,phase);

	  ComplexField Wl_phase(Umu.Grid());
	  Wl_phase=Wl*phase;
	
	  std::vector<TComplex> Wl_T;
	  sliceSum(Wl_phase,Wl_T, Tdir);

	  for(int t=0;t<nt;t++){
	    op[t] += TensorRemove(Wl_T[t]) * ((double)nt / vol) / faces / Nc; // Yes this is ugly, not figured a work around
	  }
	}
	// for(int t=0;t<nt;t++){
	//   op[t]/=6.0;		// normalization over spatial mom
	// }

	for(int t=0;t<nt;t++){
	  // print real, but save complex later in the file
	  std::cout << " t "<<t<<" " <<op[t].real()<<std::endl;
	}

	// save into hdf file
	MesonFile MF;
	MF.data=op;

	
	std::string label="glueballop_S_Wl"+std::to_string(R1)+"x"+std::to_string(R2)+"_nsq"+std::to_string(nsq);
	write(WR,label,MF);

    }

    // < Re tr U_mu_nu(x) >
    static RealD avgPlaquette(const GaugeLorentz &Umu) {
        // Real double type
        RealD sumplaq = sumPlaquette(Umu);
        // gSites() is the number of global sites. there is also lSites() for local sites.
        double vol = Umu.Grid()->gSites();
        // The number of orientations. 4*3/2=6 for Nd=4, as known.
        double faces = (1.0 * Nd * (Nd - 1)) / 2.0;
        return sumplaq / vol / faces / Nc;
    }
};




int main (int argc, char ** argv)
{

  Grid_init(&argc,&argv);

  // Double precision grids
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  // GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  // GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  //////////////////////////////////////////////////////////////////////
  // You can manage seeds however you like.
  // Recommend SeedUniqueString.
  //////////////////////////////////////////////////////////////////////
  // std::vector<int> seeds4({1,2,3,4}); 
  // GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  std::string config;
  std::string outfile;
  if( argc > 1 && argv[1][0] != '-' )
  {
    config=argv[1];
    outfile=argv[2];
    std::cout << GridLogMessage << "Loading configuration from " << config << std::endl;
    std::cout << GridLogMessage << "output: " << outfile << std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, config);
  }
  else
  {
    std::cout<<GridLogMessage <<"Using hot configuration"<<std::endl;
    SU<Nc>::ColdConfiguration(Umu);
    //    SU<Nc>::HotConfiguration(RNG4,Umu);
    config="HotConfig";
    outfile=config+".h5";
  }

  
  // Let's see what we find.
  // RealD plaq = WLoops<PeriodicGimplD>::avgPlaquette(Umu);
  Hdf5Writer WR(outfile);
  WLoops<PeriodicGimplD>::sumspPlaquette(WR,Umu);

  // Coordinate mom({1,0,0,0});

  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,1,1,mom);
  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,1,2,mom);
  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,2,2,mom);
  
  // Coordinate mom2({-1,0,0,0});

  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,1,1,mom2);
  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,1,2,mom2);
  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,2,2,mom2);
  

  // Coordinate mom3({0,1,0,0});

  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,1,1,mom3);
  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,1,2,mom3);
  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,2,2,mom3);
  
  // Coordinate mom4({0,-1,0,0});

  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,1,1,mom4);
  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,1,2,mom4);
  // WLoops<PeriodicGimplD>::sumspWilson(WR,Umu,2,2,mom4);
  

  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,1,1,0);
  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,1,1,1);
  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,1,1,2);
  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,1,2,0);
  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,1,2,1);
  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,1,2,2);
  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,2,2,0);
  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,2,2,1);
  WLoops<PeriodicGimplD>::sumspWilson_nsq(WR,Umu,2,2,2);
  



  
  // // This is how you make log messages.
  // std::cout << GridLogMessage << std::setprecision(std::numeric_limits<Real>::digits10 + 1) << "Plaquette = " << plaq << std::endl;

  

  Grid_finalize();
}




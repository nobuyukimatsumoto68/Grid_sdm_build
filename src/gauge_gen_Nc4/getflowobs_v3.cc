#include "reweight2.h"
#include <sstream>
#include "fftw3.h"


template <class Gimpl, class GaugeField>
void PolyakovField( LatticeComplexD& out, const GaugeField& Umu) {  //assume Nd=4
  typedef typename Gimpl::GaugeLinkField GaugeMat;
  GaugeMat Ut(Umu.Grid()), P(Umu.Grid());
  int T = Umu.Grid()->GlobalDimensions()[3];
  int X = Umu.Grid()->GlobalDimensions()[0];
  int Y = Umu.Grid()->GlobalDimensions()[1];
  int Z = Umu.Grid()->GlobalDimensions()[2];

  Ut = peekLorentz(Umu,3); //Select temporal direction
  P = Ut;
  for (int t=1;t<T;t++){
    P = Gimpl::CovShiftForward(Ut,3,P);
  }
  RealD norm = 1.0/Nc;
  out = trace(P)*norm;
}


template <class Field>
void fft3d( Field& Fv, const Field& v, const int sign=FFT::forward ){
  std::vector<int> dim_mask({1,1,1,0}); // 3d FFT
  FFT theFFT((GridCartesian *) v.Grid());
  theFFT.FFT_dim_mask( Fv, v, dim_mask, sign );
}




template <class Impl>
class FlowObs {
  int nsteps;
  double step_size;
  int meas_interval;
  double t0;

public:
  // here forces the Impl to be of gauge fields
  // if not the compiler will complain
  INHERIT_GIMPL_TYPES(Impl);

  // necessary for HmcObservable compatibility
  typedef typename Impl::Field Field;

  FlowObs(const int nsteps,
          const double step_size,
          const int meas_interval )
    : nsteps(nsteps)
    , step_size(step_size)
    , meas_interval(meas_interval)
    , t0(step_size*nsteps)
  {
  }

  void operator()( std::vector<ComplexD>& all_polyakov,
                   std::vector<RealD>& all_T0,
                   std::vector<RealD>& all_Q,
                   std::vector<RealD>& all_t,
                   Field& Usmear, const Field &U ){
    int def_prec = std::cout.precision();

    all_polyakov.clear();
    all_T0.clear();
    all_Q.clear();

    WilsonFlow<PeriodicGimplR> WF(step_size, nsteps, meas_interval);
    WF.resetActions();		// as in https://github.com/paboyle/Grid/blob/41f344bbd3376fc6ef38a818154c3f80e971aa50/tests/smearing/Test_WilsonFlow_adaptive.cc#L100
    // don't know if it's necessary, but this simplifies stdout prints of "Energy density (plaq)" every iteration..
    // of course it doesn't change any result
    
    WF.addMeasurement(meas_interval,
                      [&all_polyakov,&WF,&all_T0,&all_Q, &all_t](int step, RealD t, const typename Impl::GaugeField &U){
                        all_t.push_back( t );
			
                        const ComplexD p1 = WilsonLoops<Impl>::avgPolyakovLoop(U);
                        all_polyakov.push_back( p1 );
                        // all_polyakov_im.push_back( imag(p1) );
                        std::cout << GridLogMessage
                                  << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                                  << "[WilsonFlow] Polyakov loop           : "  << step << "  "
                                  << p1
                                  << std::endl;
			
                        const RealD p2 = WF.energyDensityPlaquette(t, U);
                        all_T0.push_back( p2 );
                        std::cout << GridLogMessage
                                  << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                                  << "[WilsonFlow] T0           : "  << step << "  "
                                  << p2
                                  << std::endl;

                        const RealD p3 = WilsonLoops<Impl>::TopologicalCharge(U);
                        all_Q.push_back( p3 );
                        std::cout << GridLogMessage
                                  << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                                  << "[WilsonFlow] topo charge        : "  << step << "  "
                                  << p3
                                  << std::endl;
                      });

    WF.smear(Usmear, U);
    std::cout.precision(def_prec);
  }

};



template <typename Field>
void limeWrite(const std::string filestem, Field &vec)
{
  emptyUserRecord   record;
  ScidacWriter binWriter(vec.Grid()->IsBoss());

  binWriter.open(filestem + ".lime.bin");
  binWriter.writeScidacFieldRecord(vec, record);
  binWriter.close();
}

template <typename Field>
void limeRead(Field &vec, const std::string filestem)
{
  emptyUserRecord   record;
  ScidacReader binReader;

  // binReader.open(filestem + ".lime.bin");
  binReader.open(filestem);
  binReader.readScidacFieldRecord(vec, record);
  binReader.close();
}




int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                  GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                  GridDefaultMpi());
  // GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
  //                                                                 GridDefaultSimd(Nd,vComplex::Nsimd()),
  //                                                                 GridDefaultMpi());

  // -------------------------------------
  // for reading data
  // const int conf_min=atoi(argv[1]);

  // const std::string base_dir(argv[2]); // directory of lattice config
  // // -> string path

  // const std::string basedir(argv[3]); // not used

  // const std::string basedir2(argv[4]); // directory for output .bin
  // // -> string pathO

  // const std::string mass(argv[5]);
  // const int interval=atoi(argv[6]);
  // // -> 10, 20

  // const int nbeta=atoi(argv[7]);
  // const int runtype=atoi(argv[8]);
  // //
  // const int Nt=8;
  // const int cinv=2;//atoi(argv[8]);

  // int type = 1;
  // if(mass=="0") type = 0;

  // -----------------------
  // for reading data
  const std::string infilebasedir=argv[1];
  const std::string infileid=argv[2];
  const std::string outfilebasename=argv[3];
  const RealD beta=stod(argv[4]);
  const int conf_min=atoi(argv[5]);
  const int conf_max=atoi(argv[6]);
  const int interval=atoi(argv[7]);
  // -----------------------
  
  const std::string infilebasename=infilebasedir + "/" + infileid;
  std::cout << "infilebasename = " << infilebasename << std::endl;
  std::cout << "outfilebasename = " << outfilebasename << std::endl;
  std::cout << "beta = " << beta << std::endl;

  const int Nt=8;
  const int cinv=2;//atoi(argv[8]);

  // std::vector<std::string> betas;
  // {
  //   for(int i=9; i<9+nbeta; i++) {
  //     std::string str(argv[i]);
  //     betas.push_back(str);
  //   }
  // }

  // for(auto elem : betas) std::cout << elem << std::endl;

  using Impl = PeriodicGimplR;

  const double step_size = 0.01;
  const int nsteps = 200; // (int)(Nt*Nt/cinv/cinv/8/step_size);
  const int meas_interval = 20;
  FlowObs<PeriodicGimplD> flowobs(nsteps, step_size, meas_interval);


  // for(std::string beta : betas){
  std::cout << "beta = " << beta << std::endl;
  WilsonGaugeActionR Waction(beta);
  // WilsonGaugeActionR Waction(std::stod(beta));

  // int conf_max=conf_min;
  // for(int conf=conf_min; conf<conf_min+10000000*interval; conf+=interval){
  //   char f[200];
  //   std::string prefix = get_prefix(beta, mass);
  //   std::sprintf(f,
  // 		 "%s.%d",
  // 		 prefix.data(), conf);

  //   const std::string f_str = f;
  //   const std::string path = base_dir+get_configname(beta, mass, type)+"/"+f_str;
  //   // std::cout << "path = " << path << std::endl;

  //   conf_max = conf-interval;
  //   if( !std::filesystem::exists( path ) ) {
  //     std::cout << "skpping path = " << path << std::endl;
  //     break;
  //   }
  // }

  for(int conf=conf_min; conf<conf_max; conf+=interval){
    Obs obs;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3f", beta);
    std::string betastr(buf);
    obs.beta = betastr;

    obs.base_dir=infilebasedir;
    obs.id=infileid;

    // obs.beta = beta;
    // obs.base_dir=base_dir;
    // obs.id=get_configname(beta, mass, type);

    // std::cout << "pathO" << std::endl;
    // std::cout << "pathO = " << pathO << std::endl;
    std::string pathO = outfilebasename+std::to_string(conf)+".h5";
    std::string pathO2 = outfilebasename+"P"+std::to_string(conf)+".bin";
    std::string pathO3 = outfilebasename+"PF"+std::to_string(conf)+".bin";
    // {
    //   bool exists = false;
    //   if( std::filesystem::exists( pathO )
    // 	  && std::filesystem::exists( pathO2 )
    // 	  && std::filesystem::exists( pathO3 )
    // 	  && (runtype==1) ) {
    // 	exists = true;
    //   }
    //   if(exists) continue;
    // }

    // char f[200];
    // std::string prefix = get_prefix(beta, mass);
    // std::sprintf(f,
    // 		 "%s.%d",
    // 		 prefix.data(), conf);

    // std::cout << "path" << std::endl;
    // const std::string f_str = f;
    // const std::string path = obs.base_dir+obs.id+"/"+f_str;
    // std::cout << "path2 = " << path << std::endl;

    
    //============================================
    // Read conf
    FieldMetaData header;
    LatticeGaugeField U(UGrid);
    const std::string path = infilebasename+std::to_string(conf);
    if( !std::filesystem::exists( path ) ) continue;
    NerscIO::readConfiguration(U, header, path);
    std::cout << "read " << path << std::endl;


    //============================================
    Impl::Field Usmear = U;
    std::vector<ComplexD> all_polyakov;
    std::vector<RealD> all_T0;
    std::vector<RealD> all_Q;
    std::vector<RealD> all_t;

    const ComplexD S    = Waction.S(U);
    obs.energy = real(S)/beta;

    // std::cout << "debug 1." << std::endl;
    flowobs( all_polyakov, all_T0, all_Q, all_t, Usmear, U );
    // std::cout << "debug 2." << std::endl;

#ifdef HAVE_LIME
    // std::cout << "debug 3." << std::endl;
    LatticeComplexD pfield( UGrid );
    PolyakovField<PeriodicGimplD>( pfield, Usmear );
    // std::cout << "debug 4." << std::endl;
    LatticeComplexD pfieldF( UGrid );
    fft3d<LatticeComplexD>( pfieldF, pfield );
    // std::cout << "debug 5." << std::endl;
    // {
    //   // if (UGrid->IsBoss()){
    //   // std::cout << pfieldF << std::endl;
    //   // if (UGrid->IsBoss()){
    //   std::string pathO3 = basedir2+obs.id+"PF"+std::to_string(conf)+".bin";
    //   emptyUserRecord record;
    //   std::cout << "debug 5.1" << std::endl;
    //   GridLimeWriter wr(UGrid->IsBoss());
    //   std::cout << "debug 5.2" << std::endl;
    //   // Grid::ScidacWriter wr;
    //   wr.open( pathO3 );
    //   std::cout << "debug 5.3" << std::endl;
    //   wr.writeLimeLatticeBinaryObject( pfieldF, "pfieldF" );
    //   std::cout << "debug 5.4" << std::endl;
    //   wr.close();
    //   std::cout << "debug 5.5" << std::endl;
    //   // }
    // }
    {
      // LatticeGaugeField Usave(UGrid);
      // LatticeRealD pfieldR( UGrid );
      // pfieldR = toReal( real( pfield ) );
      // if (UGrid->IsBoss()){
      // std::cout << pfieldF << std::endl;
      // if (UGrid->IsBoss()){
      emptyUserRecord record;
      // std::cout << "debug 5.1" << std::endl;
      Grid::ScidacWriter wr(UGrid->IsBoss());
      // std::cout << "debug 5.2" << std::endl;
      // Grid::ScidacWriter wr;
      wr.open( pathO3 );
      // std::cout << "debug 5.3" << std::endl;
      wr.writeScidacFieldRecord( pfieldF, record );
      // wr.writeScidacFieldRecord( pfieldR, record );
      // std::cout << "debug 5.4" << std::endl;
      wr.close();
      // std::cout << "debug 5.5" << std::endl;
      // }
    }
    // std::cout << "debug 5.6" << std::endl;
    {
      // if (UGrid->IsBoss()){
      // std::cout << "debug 5.7" << std::endl;
      emptyUserRecord record;
      // std::cout << "debug 5.8" << std::endl;
      Grid::ScidacWriter wr(UGrid->IsBoss());
      // std::cout << "debug 5.9" << std::endl;
      wr.open( pathO2 );
      // std::cout << "debug 5.11" << std::endl;
      wr.writeScidacFieldRecord( pfield, record );
      // std::cout << "debug 5.12" << std::endl;
      wr.close();
      // }
    }
#endif
    {
      std::unique_ptr<Hdf5Writer> WR;
      if (UGrid->IsBoss()){
	WR = std::make_unique<Hdf5Writer>( pathO );
	write(*WR, "obs", obs );
	write(*WR, "polyakov", all_polyakov );
	write(*WR, "T0", all_T0 );
	write(*WR, "Q", all_Q );
	write(*WR, "t", all_t );
	// Hdf5Writer WR(pathO);
	// write(WR, "obs", obs );
	// write(WR, "polyakov", all_polyakov );
	// write(WR, "T0", all_T0 );
	// write(WR, "Q", all_Q );
      }
    }
    // }

  }

  Grid_finalize();
}




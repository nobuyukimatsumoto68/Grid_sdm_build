/*
 * Warning: This code illustrative only: not well tested, and not meant for production use
 * without regression / tests being applied
 */
// 090924 sungwoo
// calculates mesons 2pts for G5,GiG5 and additional PJ5q
// along both time and spatial (x) direction
// 
// 010424 sungwoo
// add scalar meson

#include <Grid/Grid.h>
// HDF5 C API for the output-completeness check in main (H5Fopen / H5Lexists).
#include <hdf5.h>
#include <filesystem>  // directory_iterator for the --dir loop mode (like disc)
#include <algorithm>   // std::sort / std::all_of for config auto-detect
#include <cassert>     // assert in config auto-detect
#include <ctime>       // std::time for the in-binary wall blocker
#include <cstdlib>     // std::getenv / std::atol for the wall blocker
#include <fstream>     // std::ifstream for the NERSC header pre-check


using namespace std;
using namespace Grid;

// 17-momentum table matching average_trace2.f90 (where entry 3 denotes -1).
// Columns are (px,py,pz) in units of 2\pi/NS; the connected meson is projected
// onto each of these at the sink to match the disconnected loop momenta.
static const int n_mom = 17;
static const int mom_table[n_mom][3] = {
  { 0, 0, 0}, //  0
  { 1, 0, 0}, //  1
  { 0, 1, 0}, //  2
  { 0, 0, 1}, //  3
  { 1, 1, 0}, //  4
  { 1, 0, 1}, //  5
  { 0, 1, 1}, //  6
  { 1,-1, 0}, //  7
  { 1, 0,-1}, //  8
  { 0, 1,-1}, //  9
  { 1, 1, 1}, // 10
  {-1, 1, 1}, // 11
  { 1,-1, 1}, // 12
  { 1, 1,-1}, // 13
  { 2, 0, 0}, // 14
  { 0, 2, 0}, // 15
  { 0, 0, 2}, // 16
};

// Channel names, single source of truth shared by MesonTrace_hdf (which writes
// the datasets) and output_complete (which checks they are all present). The
// order and entries must match the Gammas table in MesonTrace_hdf.
static const int n_channel = 10;
static const char *channel_names[n_channel] = {
  "I_I","G5_G5",
  "GTG5_GTG5",
  "GXG5_GXG5",
  "GYG5_GYG5",
  "GZG5_GZG5",
  "GT_GT",
  "GX_GX",
  "GY_GY",
  "GZ_GZ",
};

// Build the phase exp(i p.x). Copied from Mobius_mesons_xt.cc:50.
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

void PointSource(Coordinate &coor,LatticePropagator &source)
{
  //  Coordinate coor({0,0,0,0});
  source=Zero();
  SpinColourMatrix kronecker; kronecker=1.0;
  pokeSite(kronecker,source,coor);
}

// template<class Action>
// void Solve(Action &D,LatticePropagator &source,LatticePropagator &propagator,LatticePropagator &prop5)
template<class Action>
void Solve(Action &D,LatticePropagator &source,LatticePropagator &propagator)
{
  GridBase *UGrid = D.GaugeGrid();
  GridBase *FGrid = D.FermionGrid();

  LatticeFermion src4  (UGrid); 
  LatticeFermion src5  (FGrid); 
  LatticeFermion result5(FGrid);
  LatticeFermion result4(UGrid);
  // LatticePropagator prop5(FGrid);
  
  ConjugateGradient<LatticeFermion> CG(1.0e-9,100000);
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

      // FermToProp<Action>(prop5,result5,s,c);
      FermToProp<Action>(propagator,result4,s,c);
    }
  }
  // LatticePropagator Axial_mu(UGrid); 

  // LatticeComplex    PA (UGrid); 
  // // LatticeComplex    PJ5q(UGrid);

  // std::vector<TComplex> sumPA;
  // // std::vector<TComplex> sumPJ5q;

  // Gamma g5(Gamma::Algebra::Gamma5);
  // D.ContractConservedCurrent(prop5,prop5,Axial_mu,source,Current::Axial,Tdir);
  // PA       = trace(g5*Axial_mu);      // Pseudoscalar-Axial conserved current
  // sliceSum(PA,sumPA,Tdir);

  // int Nt{static_cast<int>(sumPA.size())};

  // for(int t=0;t<Nt;t++) std::cout<<GridLogMessage <<"PAc["<<t<<"] "<<real(TensorRemove(sumPA[t]))*LCscale<<std::endl;
  
  // // D.ContractJ5q(prop5,PJ5q);
  // // sliceSum(PJ5q,sumPJ5q,Tdir);
  // // for(int t=0;t<Nt;t++) std::cout<<GridLogMessage <<"PJ5q["<<t<<"] "<<real(TensorRemove(sumPJ5q[t]))<<std::endl;


}


class MesonFile: Serializable {
public:
  // GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFile, std::vector<std::vector<Complex> >, data);
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFile, std::vector<Complex>,  data);
};

void MesonTrace_hdf(Hdf5Writer *WR,LatticePropagator &q1,LatticePropagator &q2)
{
  // const int nchannel=4;
  // Gamma::Algebra Gammas[nchannel][2] = {
  //   {Gamma::Algebra::Gamma5      ,Gamma::Algebra::Gamma5},
  //   {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::GammaTGamma5},
  //   {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::Gamma5},
  //   {Gamma::Algebra::Gamma5      ,Gamma::Algebra::GammaTGamma5}
  // };

  // std::vector<std::string> channel_name({"G5_G5",
  // 					 "GTG5_GTG5",
  // 					 "GTG5_G5",
  // 					 "G5_GTG5"});
  
  const int nchannel=n_channel;
  Gamma::Algebra Gammas[nchannel][2] = {
    {Gamma::Algebra::Identity    ,Gamma::Algebra::Identity},
    {Gamma::Algebra::Gamma5      ,Gamma::Algebra::Gamma5},
    {Gamma::Algebra::GammaTGamma5,Gamma::Algebra::GammaTGamma5},
    {Gamma::Algebra::GammaXGamma5,Gamma::Algebra::GammaXGamma5},
    {Gamma::Algebra::GammaYGamma5,Gamma::Algebra::GammaYGamma5},
    // original had GammaYGamma5 at the sink (copy-paste bug): sink must match
    // source GammaZGamma5 so GZG5_GZG5 is consistent with trace_gzg5.
    // {Gamma::Algebra::GammaZGamma5,Gamma::Algebra::GammaYGamma5},
    {Gamma::Algebra::GammaZGamma5,Gamma::Algebra::GammaZGamma5},
    {Gamma::Algebra::GammaT      ,Gamma::Algebra::GammaT},
    {Gamma::Algebra::GammaX      ,Gamma::Algebra::GammaX},
    {Gamma::Algebra::GammaY      ,Gamma::Algebra::GammaY},
    // original had GammaY at the sink (copy-paste bug): sink must match source
    // GammaZ so GZ_GZ is consistent with trace_gz.
    // {Gamma::Algebra::GammaZ      ,Gamma::Algebra::GammaY},
    {Gamma::Algebra::GammaZ      ,Gamma::Algebra::GammaZ},
  };

  std::vector<std::string> channel_name(channel_names, channel_names + n_channel);
  
  Gamma G5(Gamma::Algebra::Gamma5);

  LatticeComplex meson_CF(q1.Grid());
  // Hdf5Writer WR(file);

  // std::vector<MesonFile> MF;
  // MF.resize(nchannel);
  for(int ch=0;ch<nchannel;ch++){

    MesonFile MF;
    
    Gamma Gsrc(Gammas[ch][0]);
    Gamma Gsnk(Gammas[ch][1]);

    meson_CF = trace(G5*adj(q1)*G5*Gsnk*q2*adj(Gsrc));

    std::vector<TComplex> meson_X;
    sliceSum(meson_CF,meson_X, Xdir);

    int nx=meson_X.size();

    std::vector<Complex> corr_X(nx);
    for(int t=0;t<nx;t++){
      corr_X[t] = TensorRemove(meson_X[t]); // Yes this is ugly, not figured a work around
      std::cout << " channel "<<ch<<" x "<<t<<" " <<corr_X[t]<<std::endl;
    }
    // MF[ch].data.push_back(corr);
    MF.data=corr_X;

    // of << channel_name[ch] << " " << MF << std::endl;
    if (WR) write(*WR,channel_name[ch]+"_x",MF);

    

    std::vector<TComplex> meson_T;
    sliceSum(meson_CF,meson_T, Tdir);

    int nt=meson_T.size();

    std::vector<Complex> corr(nt);
    for(int t=0;t<nt;t++){
      corr[t] = TensorRemove(meson_T[t]); // Yes this is ugly, not figured a work around
      std::cout << " channel "<<ch<<" t "<<t<<" " <<corr[t]<<std::endl;
    }
    // MF[ch].data.push_back(corr);
    MF.data=corr;

    // of << channel_name[ch] << " " << MF << std::endl;
    if (WR) write(*WR,channel_name[ch]+"_t",MF);

    // momentum-projected connected correlator C(p_idx,t), one dataset per
    // momentum in mom_table. Index 0 (zero momentum) reproduces the "_t"
    // dataset above as a cross-check. Phase projects the sink onto p.
    LatticeComplex phase(q1.Grid());
    LatticeComplex meson_CFp(q1.Grid());

    for(int ip=0;ip<n_mom;ip++){

      Coordinate mom({mom_table[ip][0],mom_table[ip][1],mom_table[ip][2],0});
      MakePhase(mom,phase);

      meson_CFp = meson_CF*phase;

      std::vector<TComplex> meson_Tp;
      sliceSum(meson_CFp,meson_Tp, Tdir);

      int ntp=meson_Tp.size();

      std::vector<Complex> corr_p(ntp);
      for(int t=0;t<ntp;t++){
        corr_p[t] = TensorRemove(meson_Tp[t]);
      }

      MesonFile MFp;
      MFp.data=corr_p;

      std::stringstream ss;
      ss<<channel_name[ch]<<"_t_p"<<ip;
      if (WR) write(*WR,ss.str(),MFp);

    }

  }
  // XmlWriter WR(file);
  // write(WR,channel_name[ch],MF);

}


// Baryon sector not needed for the meson study: disabled (kept for reference).
#if 0
void BaryonTrace_hdf(Hdf5Writer *WR,LatticePropagator &q)
{
  const int nchannel=2;
  std::vector<std::string> channel_name({"bar_0000","bar_2222"});
  std::vector<int> channel_spinidx({0,2});
  
  LatticeComplex baryon_CF(q.Grid());
  
  for(int ch=0;ch<nchannel;ch++){

    std::cout << GridLogMessage << " channel_name "<<channel_name[ch] <<std::endl;
    MesonFile MF;


    // Unitary matrix that relates chiral (Weyl) basis with Dirac-Pauli basis
    // gamma(W) = U * gamma(DP) * U^+
    // q(DP) = U^+ * q(W)
    // prop(W) = U^+ * prop(DP) * U
    // U= 1/sqrt(2) [ 1 1]
    //              [-1 1]
    //  = 1/sqrt(2) [ 1  0 1 0 ]
    //              [ 0  1 0 1 ]
    //              [ -1 0 1 0 ]
    //              [ 0 -1 0 1 ]
    SpinMatrix U=Zero();
    U()(0,0) = 1.0;
    U()(1,1) = 1.0;
    U()(2,0) = -1.0;
    U()(3,1) = -1.0;
    U()(0,2) = 1.0;
    U()(1,3) = 1.0;
    U()(2,2) = 1.0;
    U()(3,3) = 1.0;

    SpinMatrix Udagger=transpose(U);

    LatticePropagator q_DP = 0.5 * U * q * Udagger;

    int spinidx=channel_spinidx[ch];
    auto q00 = peekSpin(q_DP,spinidx,spinidx); // pick one diagonal component in spin space

    baryon_CF=0.;

    //===========================
    // sink colors {a,b,c,d}
    for (int a = 0; a < 4; a++) {
    for (int b = 0; b < 4; b++) { if ( b == a) continue;
    for (int c = 0; c < 4; c++) { if ( c == a || c == b) continue;
    for (int d = 0; d < 4; d++) { if ( d == a || d == b || d == c) continue;
      // Compute the parity of (a, b, c, d)
      RealD parity = ((a - b) * (a - c) * (a - d) * (b - c) * (b - d) * (c - d) < 0) ? -1.0 : 1.0;
      
      //=============================
      // source colors {ap,bp,cp,dp}
      for (int ap = 0; ap < 4; ap++) {
      for (int bp = 0; bp < 4; bp++) { if ( bp == ap) continue;
      for (int cp = 0; cp < 4; cp++) { if ( cp == ap || cp == bp) continue;
      for (int dp = 0; dp < 4; dp++) { if ( dp == ap || dp == bp || dp == cp) continue;
	    
	std::cout << GridLogMessage << "sink {"<<a<<b<<c<<d<<"}, source {"<<ap<<bp<<cp<<dp<<"}"<<std::endl;
      
	// Compute the parity of (ap, bp, cp, dp)
	RealD parity_p = ((ap - bp) * (ap - cp) * (ap - dp) * (bp - cp) * (bp - dp) * (cp - dp) < 0) ? -1.0 : 1.0;
	
	baryon_CF += +1*parity*parity_p*peekColour(q00,a,ap) * peekColour(q00,b,bp) * peekColour(q00,c,cp) * peekColour(q00,d,dp);
	
      }
      }
      }
      }
    }
    }
    }
    }
    
    std::vector<TComplex> baryon_T;
    sliceSum(baryon_CF,baryon_T, Tdir);

    int nt=baryon_T.size();

    std::vector<Complex> corr(nt);
    for(int t=0;t<nt;t++){
      corr[t] = TensorRemove(baryon_T[t]); // Yes this is ugly, not figured a work around
      std::cout << " channel "<<ch<<" t "<<t<<" " <<corr[t]<<std::endl;
    }
    // MF[ch].data.push_back(corr);
    MF.data=corr;

    // // of << channel_name[ch] << " " << MF << std::endl;
    if (WR) write(*WR,channel_name[ch]+"_t",MF);

  }
  // XmlWriter WR(file);
  // write(WR,channel_name[ch],MF);

}
#endif


// Completeness check for an existing output file.
//
// A finished run writes, for every channel, the datasets <channel>_x,
// <channel>_t and <channel>_t_p0 .. _t_p{n_mom-1} (each stored by Grid as a
// group holding a "data" dataset). output_complete returns true only if the
// file exists, is a valid HDF5 file, and contains ALL of these datasets. A
// missing or partially-written file therefore returns false and gets
// recomputed, so an interrupted job is resumed without trusting a half-filled
// file. Called on the boss rank only (single reader of the shared file).
static bool output_complete(const std::string &outfile)
{
  std::ifstream probe(outfile.c_str());
  if (!probe.good()){
    return false;
  }
  probe.close();

  if (H5Fis_hdf5(outfile.c_str()) <= 0){
    return false;
  }

  hid_t file = H5Fopen(outfile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0){
    return false;
  }

  bool complete = true;
  for(int ch=0; ch<n_channel && complete; ch++){

    std::vector<std::string> keys;
    keys.push_back(std::string(channel_names[ch]) + "_x");
    keys.push_back(std::string(channel_names[ch]) + "_t");
    for(int ip=0; ip<n_mom; ip++){
      std::stringstream ss;
      ss << channel_names[ch] << "_t_p" << ip;
      keys.push_back(ss.str());
    }

    for(size_t k=0; k<keys.size(); k++){
      std::string grp = "/" + keys[k];
      std::string dat = grp + "/data";
      // Check the group and its data dataset separately so an older HDF5 that
      // errors on a multi-part path with a missing intermediate still works.
      if (H5Lexists(file,grp.c_str(),H5P_DEFAULT) <= 0 ||
          H5Lexists(file,dat.c_str(),H5P_DEFAULT) <= 0){
        complete = false;
        break;
      }
    }
  }

  H5Fclose(file);
  return complete;
}


// Internal config-loop mode, mirroring disc_multipleGamma_binary_claude.cc: ONE
// invocation processes the whole ensemble in `dir`, looping configs and skipping
// complete outputs via output_complete (cheap -- no per-config flux run / MPI
// relaunch, unlike the positional single-config mode). The grids and the
// (gauge-independent) point source are built once by the caller / here; only the
// gauge field and fermion action are rebuilt per config. An in-binary graceful
// wall blocker (env MESON_DEADLINE_EPOCH, epoch s; 0 disables) stops on a config
// boundary before the wall, measuring per-config time as it goes;
// MESON_TPT_SECONDS optionally bootstraps the estimate for the first config. The
// stop decision and measured duration are taken on the boss and broadcast
// (GlobalSum) so every rank breaks together (the loop body is collective).
static void run_loop_mode(GridCartesian *UGrid, GridCartesian *FGrid,
                          GridRedBlackCartesian *FrbGrid, GridRedBlackCartesian *UrbGrid,
                          LatticeGaugeField &Umu,
                          const std::string &dir, const std::string &obsdir,
                          RealD mass, RealD M5, int stride_mult)
{
  std::filesystem::create_directories(obsdir);

  // auto-detect conf_min / interval / conf_max and the lat filename prefix
  // (same logic as the disc binary).
  int conf_min, conf_max, interval;
  std::string lat_prefix;
  {
    std::vector<int> confs;
    const std::string suffix = "_lat.";
    for(const auto &entry : std::filesystem::directory_iterator(dir)){
      const std::string fname = entry.path().filename().string();
      const auto pos = fname.rfind(suffix);
      if(pos == std::string::npos) continue;
      const std::string numstr = fname.substr(pos + suffix.size());
      if(numstr.empty() || !std::all_of(numstr.begin(), numstr.end(), ::isdigit)) continue;
      if(lat_prefix.empty()) lat_prefix = fname.substr(0, pos + suffix.size());
      confs.push_back(std::stoi(numstr));
    }
    assert(!confs.empty());
    std::sort(confs.begin(), confs.end());
    conf_min = confs.front();
    interval = (confs.size() >= 2) ? confs[1] - confs[0] : 1;
    conf_max = confs.back() + interval;
  }
  // --stride multiplies the native interval (1 = every config, matches disc).
  const int step = interval * (stride_mult > 0 ? stride_mult : 1);
  std::cout << GridLogMessage << "loop mode: dir=" << dir << " conf_min=" << conf_min
            << " conf_max=" << conf_max << " interval=" << interval
            << " step=" << step << std::endl;

  // in-binary wall blocker inputs (see disc); MESON_DEADLINE_EPOCH 0 disables.
  long deadline = 0;
  if(const char* e = std::getenv("MESON_DEADLINE_EPOCH")) deadline = std::atol(e);
  double max_dur = 0.0;       // largest per-config wall time seen so far (s)
  if(const char* e = std::getenv("MESON_TPT_SECONDS")) max_dur = std::atof(e); // optional bootstrap
  const double margin = 1.2;  // generous safety factor on the estimate

  // fixed action params; the point source (delta at origin) is gauge-independent
  // -> build it ONCE outside the loop. Only Umu + the action change per config.
  const RealD b = 1.5, c = 0.5;
  std::vector<Complex> boundary = {1,1,1,-1};
  MobiusFermionD::ImplParams Params(boundary);

  LatticePropagator point_source(UGrid);
  Coordinate Origin({0,0,0,0});
  PointSource(Origin, point_source);

  LatticePropagator prop(UGrid);

  for(int conf = conf_min; conf < conf_max; conf += step){

    const std::string outfile = obsdir + "/mesons_conn." + std::to_string(conf) + ".h5";

    // cheap skip of a complete output (boss decides, broadcast so ranks agree).
    uint64_t done = 0;
    if(UGrid->IsBoss()) done = output_complete(outfile) ? 1 : 0;
    UGrid->GlobalSum(done);
    if(done){
      std::cout << GridLogMessage << "skipping conf " << conf << " (output complete)" << std::endl;
      continue;
    }

    // graceful wall stop before a config we cannot finish (needs one measurement
    // first unless MESON_TPT_SECONDS bootstrapped max_dur).
    if(deadline > 0 && max_dur > 0.0){
      uint64_t stop = 0;
      if(UGrid->IsBoss()){
        if((double)std::time(nullptr) + margin*max_dur > (double)deadline) stop = 1;
      }
      UGrid->GlobalSum(stop);
      if(stop){
        std::cout << GridLogMessage << "blocker: est " << (long)(margin*max_dur)
                  << "s for next config exceeds deadline; stopping gracefully before conf "
                  << conf << std::endl;
        break;
      }
    }

    const long t_start = (long)std::time(nullptr);

    {
      const std::string cfgpath = dir + "/" + lat_prefix + std::to_string(conf);

      // guard against partially-written configs (e.g. HMC still running): read
      // first line on boss, skip rather than abort if BEGIN_HEADER is missing.
      uint64_t bad_cfg = 0;
      if(UGrid->IsBoss()){
        std::ifstream chk(cfgpath);
        std::string first_line;
        if(!chk || !std::getline(chk, first_line) || first_line != "BEGIN_HEADER") bad_cfg = 1;
      }
      UGrid->GlobalSum(bad_cfg);
      if(bad_cfg){
        std::cout << GridLogMessage << "skipping conf " << conf
                  << " (NERSC header check failed; config may be incomplete)" << std::endl;
        continue;
      }

      std::cout << GridLogMessage << "conf " << conf << " -> " << outfile << std::endl;
      FieldMetaData header;
      NerscIO::readConfiguration(Umu, header, cfgpath);
    }

    MobiusFermionD FermAct(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, c, Params);
    Solve(FermAct, point_source, prop);

    std::unique_ptr<Hdf5Writer> WR;
    if(UGrid->IsBoss()) WR = std::make_unique<Hdf5Writer>(outfile);
    MesonTrace_hdf(WR.get(), prop, prop);

    // update the per-config wall-time estimate (boss clock; broadcast so max_dur
    // is identical on every rank for the next iteration's decision).
    uint64_t dur = 0;
    if(UGrid->IsBoss()) dur = (uint64_t)((long)std::time(nullptr) - t_start);
    UGrid->GlobalSum(dur);
    if((double)dur > max_dur) max_dur = (double)dur;
    std::cout << GridLogMessage << "conf " << conf << " took " << dur
              << "s (max so far " << (long)max_dur << "s)" << std::endl;
  }
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

  //////////////////////////////////////////////////////////////////////
  // You can manage seeds however you like.
  // Recommend SeedUniqueString.
  //////////////////////////////////////////////////////////////////////
  // std::vector<int> seeds4({1,2,3,4}); 
  // GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);

  // ---- directory loop mode (like disc): if --dir is given, process the whole
  // ensemble in one invocation and return. Otherwise fall through to the
  // positional single-config / cold-config path below (preserved unchanged).
  {
    std::string dir, obsdir;
    RealD lm_M5 = 1.5, lm_mass = 0.1;
    int stride_mult = 1;
    bool loop_mode = false;
    for(int i=1; i<argc; i++){
      std::string a = argv[i];
      if      (a=="--dir"    && i+1<argc){ dir    = argv[++i]; loop_mode = true; }
      else if (a=="--obsdir" && i+1<argc){ obsdir = argv[++i]; }
      else if (a=="--mass"   && i+1<argc){ lm_mass= std::stod(argv[++i]); }
      else if (a=="--M5"     && i+1<argc){ lm_M5  = std::stod(argv[++i]); }
      else if (a=="--stride" && i+1<argc){ stride_mult = std::stoi(argv[++i]); }
    }
    if(loop_mode){
      run_loop_mode(UGrid,FGrid,FrbGrid,UrbGrid,Umu,dir,obsdir,lm_mass,lm_M5,stride_mult);
      Grid_finalize();
      return 0;
    }
  }

  std::string config;
  std::string outfile;
  RealD M5, mass;
  bool cold;
  // Parse arguments first WITHOUT touching the gauge field, so that an
  // already-complete output can be skipped below before paying for the (costly)
  // NERSC read of a config we are not going to use.
  if( argc > 1 && argv[1][0] != '-' )
  {
    cold=false;
    config=argv[1];
    M5=stod(argv[2]);
    mass=stod(argv[3]);
    outfile=argv[4];
  }
  else
  {
    cold=true;
    config="ColdConfig";
    M5=1.5;			// give some default numbers
    mass=0.1;			// give some default numbers
    outfile=config+".h5";
  }

  // Skip recomputation only if a previous run already wrote a COMPLETE output
  // file (all per-channel keys present). A missing or partially-written file is
  // recomputed and overwritten below. The boss is the sole reader of the shared
  // file; broadcast its verdict via GlobalSum so every rank agrees before the
  // collective solve, then exit cleanly. This runs BEFORE the gauge read so a
  // complete config is skipped without reading it.
  uint64_t already_done = 0;
  if (UGrid->IsBoss()){
    already_done = output_complete(outfile) ? 1 : 0;
  }
  UGrid->GlobalSum(already_done);
  if (already_done){
    std::cout << GridLogMessage << "output complete, skipping: " << outfile << std::endl;
    Grid_finalize();
    return 0;
  }

  // Now load the gauge field for the configs we actually need to compute.
  if( !cold )
  {
    std::cout << GridLogMessage << "Loading configuration from " << config << std::endl;
    std::cout << GridLogMessage << "M5=" << M5 << std::endl;
    std::cout << GridLogMessage << "mass=" << mass << std::endl;
    std::cout << GridLogMessage << "output: " << outfile << std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, config);
  }
  else
  {
    std::cout<<GridLogMessage <<"Using hot configuration"<<std::endl;
    SU<Nc>::ColdConfiguration(Umu);
    //    SU<Nc>::HotConfiguration(RNG4,Umu);
  }

  std::vector<RealD> masses({ mass} ); // u/d, s, c ??
  // put just a single mass from input

  int nmass = masses.size();

  std::vector<MobiusFermionD *> FermActs;
  
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion action as Scaled Shamir kernel"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;

  for(auto mass: masses) {

    RealD b=1.5;// Scale factor b+c=2, b-c=1
    RealD c=0.5;

    std::cout << GridLogMessage << "==============================================" << std::endl;
    std::cout << GridLogMessage << "anti-periodic boundary condition {1,1,1,-1}" << std::endl;
    std::vector<Complex> boundary = {1,1,1,-1};
    MobiusFermionD::ImplParams Params(boundary);

    FermActs.push_back(new MobiusFermionD(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,Params));
   
  }

  LatticePropagator point_source(UGrid);
  // LatticePropagator wall_source(UGrid);
  // LatticePropagator gaussian_source(UGrid);

  Coordinate Origin({0,0,0,0});
  PointSource   (Origin,point_source);
  // WallSource  (0,wall_source);
  // Z2WallSource  (RNG4,0,wall_source);
  // GaussianSource(Origin,Umu,gaussian_source);
  
  std::vector<LatticePropagator> PointProps(nmass,UGrid);
  // std::vector<LatticePropagator> PointProps5(nmass,FGrid);
  // std::vector<LatticePropagator> GaussProps(nmass,UGrid);
  // std::vector<LatticePropagator> Z2Props   (nmass,UGrid);
  // std::vector<LatticePropagator> zeromomProps   (nmass,UGrid);


  // Hdf5Writer WR(outfile);
  // UGrid->Barrier();
  // if (UGrid->IsBoss()){
  // Hdf5Writer WRtest("test.h5");
  // }
  // UGrid->Barrier();
  // std::cout << GridLogMessage << "TEST SUNGWOO" << std::endl;

  // Hdf5Writer WR(outfile);
  std::unique_ptr<Hdf5Writer> WR;
  if (UGrid->IsBoss()){
    WR = std::make_unique<Hdf5Writer>(outfile);
  }


  
  for(int m=0;m<nmass;m++) {
    
    Solve(*FermActs[m],point_source   ,PointProps[m]); //, PointProps5[m]);
    // Solve(*FermActs[m],point_source   ,PointProps[m], PointProps5[m]);
    // Meson5Trace_hdf(*FermActs[m],WR,point_source  ,PointProps5[m]);
    // Meson5Trace(*FermActs[m],point_source  ,PointProps5[m]);
    // Solve(*FermActs[m],gaussian_source,GaussProps[m]);
    // Solve(*FermActs[m],wall_source    ,Z2Props[m]);
    // Solve(*FermActs[m],wall_source    ,zeromomProps[m]);
  
  }

  // LatticeComplex phase(UGrid);
  // Coordinate mom({0,0,0,0});
  // MakePhase(mom,phase);
  

  // std::ofstream of;
  // of.open( outfile, std::ios::out | std::ios::trunc);
  // if(!of) assert(false);
  // of << std::scientific << std::setprecision(15);

  

  for(int m1=0 ;m1<nmass;m1++) {
  for(int m2=m1;m2<nmass;m2++) {
    std::stringstream ssp,ssg,ssz;

    ssp<<config<< "_m" << m1 << "_m"<< m2 << "_point_meson.xml";
    // ssg<<config<< "_m" << m1 << "_m"<< m2 << "_smeared_meson.xml";
    // ssz<<config<< "_m" << m1 << "_m"<< m2 << "_wall_meson.xml";

    // sungwoo: phase not even being used
    MesonTrace_hdf(WR.get(),PointProps[m1],PointProps[m2]);
    // baryon sector disabled for the meson study
    // BaryonTrace_hdf(WR.get(),PointProps[m1]);
    // MesonTrace(ssp.str(),PointProps[m1],PointProps[m2]);
    // MesonTrace(ssg.str(),GaussProps[m1],GaussProps[m2],phase);
    // MesonTrace(ssz.str(),Z2Props[m1],Z2Props[m2],phase);
    // MesonTrace(ssz.str(),zeromomProps[m1],zeromomProps[m2],phase);
  }}

  Grid_finalize();
}




// hasenbusch (062626): Hasenbusch mass-preconditioned variant, built on v5
//     (inherits MX_inner=10000). The single EOFA factor det[D(mass)/D(pv)] is
//     split into k=2 factors det[D(mass)/D(m1)] * det[D(m1)/D(pv)] with
//     m1 = sqrt(mass*pv), each a one-flavour EOFA ratio, distributed over a
//     3-level integrator (light / heavy / gauge). The single-factor reference
//     is dweofa_mobius_HSDM_v5_claude.cc. RNG unchanged (continues from ckpoint).
//     Refs: Hasenbusch hep-lat/0107019; EOFA Chen-Chiu arXiv:1403.1683.
// v3: 033024 sungwoo
//     made compatible with Grid version 03/26/24 (da59379)
// 041224 sungwoo: time anti-periodic BC

/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_EODWFRatio.cc

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>
#include <cstdlib>   // getenv/atof for the M1_HASEN tuning override

namespace Grid{
  struct FermionParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionParameters,
				    int, Ls,
				    double, mass,
				    double, M5,
				    double, b,
				    double, c,
				    double, StoppingCondition,
				    int, MaxCGIterations,
				    bool, ApplySmearing);

    //template <class ReaderClass >
    //FermionParameters(Reader<ReaderClass>& Reader){
    //  read(Reader, "Mobius", *this);
    //}

  };

  
  struct MobiusHMCParameters: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(MobiusHMCParameters,
				  double, gauge_beta,
				  FermionParameters, Mobius)

  template <class ReaderClass >
  MobiusHMCParameters(Reader<ReaderClass>& Reader){
    read(Reader, "Action", *this);
  }

};

  struct SmearingParam: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(SmearingParam,
				    double, rho,
				    Integer, Nsmear)

    template <class ReaderClass >
    SmearingParam(Reader<ReaderClass>& Reader){
      read(Reader, "StoutSmearing", *this);
    }

  };
  
  
}






NAMESPACE_BEGIN(Grid);


#define MIXED_PRECISION

/*
 * Need a plan for gauge field update for mixed precision in HMC                      (2x speed up)
 *    -- Store the single prec action operator.
 *    -- Clone the gauge field from the operator function argument.
 *    -- Build the mixed precision operator dynamically from the passed operator and single prec clone.
 */

template<class FermionOperatorD, class FermionOperatorF, class SchurOperatorD, class  SchurOperatorF> 
class MixedPrecisionConjugateGradientOperatorFunction : public OperatorFunction<typename FermionOperatorD::FermionField> {
public:
  typedef typename FermionOperatorD::FermionField FieldD;
  typedef typename FermionOperatorF::FermionField FieldF;

  using OperatorFunction<FieldD>::operator();

  RealD   Tolerance;
  RealD   InnerTolerance; //Initial tolerance for inner CG. Defaults to Tolerance but can be changed
  Integer MaxInnerIterations;
  Integer MaxOuterIterations;
  GridBase* SinglePrecGrid4; //Grid for single-precision fields
  GridBase* SinglePrecGrid5; //Grid for single-precision fields
  RealD OuterLoopNormMult; //Stop the outer loop and move to a final double prec solve when the residual is OuterLoopNormMult * Tolerance

  FermionOperatorF &FermOpF;
  FermionOperatorD &FermOpD;;
  SchurOperatorF &LinOpF;
  SchurOperatorD &LinOpD;

  Integer TotalInnerIterations; //Number of inner CG iterations
  Integer TotalOuterIterations; //Number of restarts
  Integer TotalFinalStepIterations; //Number of CG iterations in final patch-up step

  MixedPrecisionConjugateGradientOperatorFunction(RealD tol, 
                                                  Integer maxinnerit, 
                                                  Integer maxouterit, 
                                                  GridBase* _sp_grid4, 
                                                  GridBase* _sp_grid5, 
                                                  FermionOperatorF &_FermOpF,
                                                  FermionOperatorD &_FermOpD,
                                                  SchurOperatorF   &_LinOpF,
                                                  SchurOperatorD   &_LinOpD): 
  LinOpF(_LinOpF),
    LinOpD(_LinOpD),
    FermOpF(_FermOpF),
    FermOpD(_FermOpD),
    Tolerance(tol), 
    InnerTolerance(tol), 
    MaxInnerIterations(maxinnerit), 
    MaxOuterIterations(maxouterit), 
    SinglePrecGrid4(_sp_grid4),
    SinglePrecGrid5(_sp_grid5),
    OuterLoopNormMult(100.) 
  { 
    /* Debugging instances of objects; references are stored
       std::cout << GridLogMessage << " Mixed precision CG wrapper LinOpF " <<std::hex<< &LinOpF<<std::dec <<std::endl;
       std::cout << GridLogMessage << " Mixed precision CG wrapper LinOpD " <<std::hex<< &LinOpD<<std::dec <<std::endl;
       std::cout << GridLogMessage << " Mixed precision CG wrapper FermOpF " <<std::hex<< &FermOpF<<std::dec <<std::endl;
       std::cout << GridLogMessage << " Mixed precision CG wrapper FermOpD " <<std::hex<< &FermOpD<<std::dec <<std::endl;
    */
  };

  void operator()(LinearOperatorBase<FieldD> &LinOpU, const FieldD &src, FieldD &psi) {

    std::cout << GridLogMessage << " Mixed precision CG wrapper operator() "<<std::endl;

    SchurOperatorD * SchurOpU = static_cast<SchurOperatorD *>(&LinOpU);
      
    //      std::cout << GridLogMessage << " Mixed precision CG wrapper operator() FermOpU " <<std::hex<< &(SchurOpU->_Mat)<<std::dec <<std::endl;
    //      std::cout << GridLogMessage << " Mixed precision CG wrapper operator() FermOpD " <<std::hex<< &(LinOpD._Mat) <<std::dec <<std::endl;
    // Assumption made in code to extract gauge field
    // We could avoid storing LinopD reference alltogether ?
    assert(&(SchurOpU->_Mat)==&(LinOpD._Mat));

    ////////////////////////////////////////////////////////////////////////////////////
    // Must snarf a single precision copy of the gauge field in Linop_d argument
    ////////////////////////////////////////////////////////////////////////////////////
    typedef typename FermionOperatorF::GaugeField GaugeFieldF;
    typedef typename FermionOperatorF::GaugeLinkField GaugeLinkFieldF;
    typedef typename FermionOperatorD::GaugeField GaugeFieldD;
    typedef typename FermionOperatorD::GaugeLinkField GaugeLinkFieldD;

    GridBase * GridPtrF = SinglePrecGrid4;
    GridBase * GridPtrD = FermOpD.Umu.Grid();
    GaugeFieldF     U_f  (GridPtrF);
    GaugeLinkFieldF Umu_f(GridPtrF);
    //      std::cout << " Dim gauge field "<<GridPtrF->Nd()<<std::endl; // 4d
    //      std::cout << " Dim gauge field "<<GridPtrD->Nd()<<std::endl; // 4d

    ////////////////////////////////////////////////////////////////////////////////////
    // Moving this to a Clone method of fermion operator would allow to duplicate the 
    // physics parameters and decrease gauge field copies
    ////////////////////////////////////////////////////////////////////////////////////
    GaugeLinkFieldD Umu_d(GridPtrD);
    for(int mu=0;mu<Nd*2;mu++){ 
      Umu_d = PeekIndex<LorentzIndex>(FermOpD.Umu, mu);
      precisionChange(Umu_f,Umu_d);
      PokeIndex<LorentzIndex>(FermOpF.Umu, Umu_f, mu);
    }
    pickCheckerboard(Even,FermOpF.UmuEven,FermOpF.Umu);
    pickCheckerboard(Odd ,FermOpF.UmuOdd ,FermOpF.Umu);

    ////////////////////////////////////////////////////////////////////////////////////
    // Could test to make sure that LinOpF and LinOpD agree to single prec?
    ////////////////////////////////////////////////////////////////////////////////////
    /*
      GridBase *Fgrid = psi._grid;
      FieldD tmp2(Fgrid);
      FieldD tmp1(Fgrid);
      LinOpU.Op(src,tmp1);
      LinOpD.Op(src,tmp2);
      std::cout << " Double gauge field "<< norm2(FermOpD.Umu)<<std::endl;
      std::cout << " Single gauge field "<< norm2(FermOpF.Umu)<<std::endl;
      std::cout << " Test of operators "<<norm2(tmp1)<<std::endl;
      std::cout << " Test of operators "<<norm2(tmp2)<<std::endl;
      tmp1=tmp1-tmp2;
      std::cout << " Test of operators diff "<<norm2(tmp1)<<std::endl;
    */

    ////////////////////////////////////////////////////////////////////////////////////
    // Make a mixed precision conjugate gradient
    ////////////////////////////////////////////////////////////////////////////////////
    MixedPrecisionConjugateGradient<FieldD,FieldF> MPCG(Tolerance,MaxInnerIterations,MaxOuterIterations,SinglePrecGrid5,LinOpF,LinOpD);
    std::cout << GridLogMessage << "Calling mixed precision Conjugate Gradient" <<std::endl;
    MPCG(src,psi);
  }
};

NAMESPACE_END(Grid);




int main(int argc, char **argv) {
  using namespace Grid;
   ;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  // typedef GenericHMCRunner<ForceGradient> HMCWrapper;  // Uses the default minimum norm
  // typedef WilsonImplR FermionImplPolicy;
  // typedef MobiusFermionD FermionAction;
  // typedef typename FermionAction::FermionField FermionField;
  typedef WilsonImplR FermionImplPolicy;
  typedef MobiusFermionD FermionAction;
  typedef MobiusFermionF FermionActionF;
  typedef MobiusEOFAFermionD FermionEOFAAction;
  typedef MobiusEOFAFermionF FermionEOFAActionF;
  typedef typename FermionAction::FermionField FermionField;
  typedef typename FermionActionF::FermionField FermionFieldF;
  // Serialiser
  typedef Grid::XmlReader       Serialiser;
  
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;
  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
 
  // Reader, file should come from command line
  if (TheHMC.ParameterFile.empty()){
    std::cout << "Input file not specified."
              << "Use --ParameterFile option in the command line.\nAborting" 
              << std::endl;
    exit(1);
  }
  Serialiser Reader(TheHMC.ParameterFile);

  MobiusHMCParameters MyParams(Reader);  
  // Apply smearing to the fermionic action
  bool ApplySmearing = MyParams.Mobius.ApplySmearing;
  
  

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader

  // EOFA parameters
  OneFlavourRationalParams OFRp;
  OFRp.lo       = 0.1;
  OFRp.hi       = 16.0;
  OFRp.MaxIter  = 10000;
  OFRp.tolerance= 1.0e-7;
  OFRp.degree   = 8;
  OFRp.precision= 40;
  // OFRp.degree   = 8;
  // OFRp.precision= 30;

  double ActionStoppingCondition     = 1e-10;
  double DerivativeStoppingCondition = 1e-6;
  double MaxCGIterations = 30000;


  // Checkpointer definition (Name: Checkpointer)
  CheckpointerParameters CPparams(Reader);

  TheHMC.Resources.LoadNerscCheckpointer(CPparams);
  //  TheHMC.Resources.LoadBinaryCheckpointer(CPparams);

  // RNG definition (Name: RandomNumberGenerator)
  RNGModuleParameters RNGpar(Reader);
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // Plaquette and Polyakov loop
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
    
  typedef PolyakovMod<HMCWrapper::ImplPolicy> PolyakovObs;
  TheHMC.Resources.AddObservable<PolyakovObs>();

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard

  WilsonGaugeActionR Waction(MyParams.gauge_beta);
//   SymanzikGaugeActionR Syzaction(MyParams.gauge_beta);
    
  const int Ls   = MyParams.Mobius.Ls;
  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  Coordinate latt  = GridDefaultLatt();
  Coordinate mpi   = GridDefaultMpi();
  Coordinate simdF = GridDefaultSimd(Nd,vComplexF::Nsimd());
  Coordinate simdD = GridDefaultSimd(Nd,vComplexD::Nsimd());
  auto GridPtrF   = SpaceTimeGrid::makeFourDimGrid(latt,simdF,mpi);
  auto GridRBPtrF = SpaceTimeGrid::makeFourDimRedBlackGrid(GridPtrF);
  auto FGridF     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtrF);
  auto FrbGridF   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtrF);


  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);
  LatticeGaugeFieldF UF(GridPtrF);

  typedef SchurDiagMooeeOperator<FermionActionF,FermionFieldF> LinearOperatorF;
  typedef SchurDiagMooeeOperator<FermionAction ,FermionField > LinearOperatorD;
  typedef SchurDiagMooeeOperator<FermionEOFAActionF,FermionFieldF> LinearOperatorEOFAF;
  typedef SchurDiagMooeeOperator<FermionEOFAAction ,FermionField > LinearOperatorEOFAD;

  typedef MixedPrecisionConjugateGradientOperatorFunction<MobiusFermionD,MobiusFermionF,LinearOperatorD,LinearOperatorF> MxPCG;
  typedef MixedPrecisionConjugateGradientOperatorFunction<MobiusEOFAFermionD,MobiusEOFAFermionF,LinearOperatorEOFAD,LinearOperatorEOFAF> MxPCG_EOFA;


  Real mass = MyParams.Mobius.mass; // light quark mass, e.g. 0.01
  Real pv   = 1.0;                  // Pauli-Villars
  RealD M5  = MyParams.Mobius.M5; //1.5;
  RealD b   = MyParams.Mobius.b; //  3./2.;
  RealD c   = MyParams.Mobius.c; //  1./2.;

  // Hasenbusch mass preconditioning (k=2): split the one-flavour determinant
  //   det[ D(mass)/D(pv) ] = det[ D(mass)/D(m1) ] * det[ D(m1)/D(pv) ],
  // intermediate mass m1 = sqrt(mass*pv) (= 0.1 for mass=0.01, pv=1.0). The
  // LIGHT factor det[D(mass)/D(m1)] carries the IR; the HEAVY factor
  // det[D(m1)/D(pv)] is well conditioned. Ref: Hasenbusch, hep-lat/0107019.
  Real m1 = std::sqrt(mass*pv);
  // tuning override: env M1_HASEN sets the intermediate mass (fallback: geometric mean
  // sqrt(mass*pv)). Used by the run_hmc_scan_claude.sh m1 scan; the print below echoes it.
  if (const char* e_m1 = std::getenv("M1_HASEN")) { m1 = std::atof(e_m1); }

  // These lines are unecessary if BC are all periodic
  std::cout << GridLogMessage << "boundary condition {1,1,1,-1}" << std::endl;
  std::cout << GridLogMessage << "Hasenbusch k=2: light det[D(" << mass << ")/D(" << m1
            << ")] * heavy det[D(" << m1 << ")/D(" << pv << ")]" << std::endl;
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);

  // EOFA operator (AbstractEOFAFermion.h): D(mq1) + shift*g5*R5*Delta_pm(mq2,mq3)*P_pm.
  // For a factor det[D(m_lo)/D(m_hi)]:
  //   Op_L(mq1=m_lo, mq2=m_lo, mq3=m_hi, shift= 0.0, pm=-1)
  //   Op_R(mq1=m_hi, mq2=m_lo, mq3=m_hi, shift=-1.0, pm=+1)
  // (single-prec twins _F carry the same masses for the mixed-prec inner CG.)

  // ---------------- LIGHT factor: det[D(mass)/D(m1)] ----------------
  MobiusEOFAFermionD Op_L_light (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , mass, mass, m1, 0.0, -1, M5, b, c, Params);
  MobiusEOFAFermionF Op_LF_light(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, mass, mass, m1, 0.0, -1, M5, b, c, Params);
  MobiusEOFAFermionD Op_R_light (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , m1,   mass, m1, -1.0, 1, M5, b, c, Params);
  MobiusEOFAFermionF Op_RF_light(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, m1,   mass, m1, -1.0, 1, M5, b, c, Params);

  // ---------------- HEAVY factor: det[D(m1)/D(pv)] ----------------
  MobiusEOFAFermionD Op_L_heavy (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , m1, m1, pv, 0.0, -1, M5, b, c, Params);
  MobiusEOFAFermionF Op_LF_heavy(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, m1, m1, pv, 0.0, -1, M5, b, c, Params);
  MobiusEOFAFermionD Op_R_heavy (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , pv, m1, pv, -1.0, 1, M5, b, c, Params);
  MobiusEOFAFermionF Op_RF_heavy(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, pv, m1, pv, -1.0, 1, M5, b, c, Params);


  // Shared double-precision heatbath solver (stateless between calls, reused by
  // both Hasenbusch factors).
  ConjugateGradient<FermionField> ActionCG(ActionStoppingCondition,MaxCGIterations);

  // const int MX_inner = 1000;                                                  // v4 original
  const int MX_inner = 10000;  // inherited from v5: inner single-prec CG cap (1000 truncated at m=0.01)

  // ---- LIGHT factor: linear operators + mixed-prec CGs ----
  LinearOperatorEOFAD LinOp_L_light (Op_L_light);
  LinearOperatorEOFAD LinOp_R_light (Op_R_light);
  LinearOperatorEOFAF LinOp_LF_light(Op_LF_light);
  LinearOperatorEOFAF LinOp_RF_light(Op_RF_light);

  MxPCG_EOFA ActionCGL_light(ActionStoppingCondition,
                             MX_inner, MaxCGIterations, GridPtrF, FrbGridF,
                             Op_LF_light, Op_L_light, LinOp_LF_light, LinOp_L_light);
  MxPCG_EOFA ActionCGR_light(ActionStoppingCondition,
                             MX_inner, MaxCGIterations, GridPtrF, FrbGridF,
                             Op_RF_light, Op_R_light, LinOp_RF_light, LinOp_R_light);
  MxPCG_EOFA DerivativeCGL_light(DerivativeStoppingCondition,
                                 MX_inner, MaxCGIterations, GridPtrF, FrbGridF,
                                 Op_LF_light, Op_L_light, LinOp_LF_light, LinOp_L_light);
  MxPCG_EOFA DerivativeCGR_light(DerivativeStoppingCondition,
                                 MX_inner, MaxCGIterations, GridPtrF, FrbGridF,
                                 Op_RF_light, Op_R_light, LinOp_RF_light, LinOp_R_light);

  // ---- HEAVY factor: linear operators + mixed-prec CGs ----
  LinearOperatorEOFAD LinOp_L_heavy (Op_L_heavy);
  LinearOperatorEOFAD LinOp_R_heavy (Op_R_heavy);
  LinearOperatorEOFAF LinOp_LF_heavy(Op_LF_heavy);
  LinearOperatorEOFAF LinOp_RF_heavy(Op_RF_heavy);

  MxPCG_EOFA ActionCGL_heavy(ActionStoppingCondition,
                             MX_inner, MaxCGIterations, GridPtrF, FrbGridF,
                             Op_LF_heavy, Op_L_heavy, LinOp_LF_heavy, LinOp_L_heavy);
  MxPCG_EOFA ActionCGR_heavy(ActionStoppingCondition,
                             MX_inner, MaxCGIterations, GridPtrF, FrbGridF,
                             Op_RF_heavy, Op_R_heavy, LinOp_RF_heavy, LinOp_R_heavy);
  MxPCG_EOFA DerivativeCGL_heavy(DerivativeStoppingCondition,
                                 MX_inner, MaxCGIterations, GridPtrF, FrbGridF,
                                 Op_LF_heavy, Op_L_heavy, LinOp_LF_heavy, LinOp_L_heavy);
  MxPCG_EOFA DerivativeCGR_heavy(DerivativeStoppingCondition,
                                 MX_inner, MaxCGIterations, GridPtrF, FrbGridF,
                                 Op_RF_heavy, Op_R_heavy, LinOp_RF_heavy, LinOp_R_heavy);

  // One EOFA pseudofermion action per Hasenbusch factor (shared ActionCG heatbath).
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy>
    EOFA_light(Op_L_light, Op_R_light,
               ActionCG,
               ActionCGL_light, ActionCGR_light,
               DerivativeCGL_light, DerivativeCGR_light,
               OFRp, true);

  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy>
    EOFA_heavy(Op_L_heavy, Op_R_heavy,
               ActionCG,
               ActionCGL_heavy, ActionCGR_heavy,
               DerivativeCGL_heavy, DerivativeCGR_heavy,
               OFRp, true);

  // Set smearing (true/false), default: false
  EOFA_light.is_smeared = ApplySmearing;
  EOFA_heavy.is_smeared = ApplySmearing;

  // Collect actions on a 3-level integrator. TheAction[0] is the OUTERMOST
  // (coarsest) level; the LAST pushed is innermost (highest force-eval frequency).
  // Multipliers are relative to the parent level and are TUNABLE:
  //   Level1 (mult 1) : LIGHT factor  -- most expensive solve, smallest force
  //   Level2 (mult 2) : HEAVY factor  -- cheap, well conditioned
  //   Level3 (mult 2) : gauge         -- stiff, cheapest force (finest steps)
  // Gauge runs 2 x 2 = 4 steps per light step, matching the v4 fermion:gauge
  // ratio. Ref idiom: HMC/Mobius2p1f_DD_EOFA_96I_3level.cc.
  ActionLevel<HMCWrapper::Field> Level1(1); // light fermion (coarsest)
  Level1.push_back(&EOFA_light);

  ActionLevel<HMCWrapper::Field> Level2(2); // heavy fermion (middle)
  Level2.push_back(&EOFA_heavy);

  ActionLevel<HMCWrapper::Field> Level3(2); // gauge (finest)
  Level3.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.TheAction.push_back(Level3);

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable
  TheHMC.Parameters.initialize(Reader);

  // Reset performance counters 

  if (ApplySmearing){
    SmearingParam SmPar(Reader);
    //double rho = 0.1;  // smearing parameter
    //int Nsmear = 3;    // number of smearing levels
    Smear_Stout<HMCWrapper::ImplPolicy> Stout(SmPar.rho);
    SmearedConfiguration<HMCWrapper::ImplPolicy> SmearingPolicy(GridPtr, SmPar.Nsmear, Stout);
    TheHMC.Run(SmearingPolicy); // for smearing
  } else {
    TheHMC.Run();  // no smearing
  }


  Grid_finalize();
} 

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


  Real mass = MyParams.Mobius.mass; //0.04;
  Real pv   = 1.0;
  RealD M5  = MyParams.Mobius.M5; //1.5;
  RealD b   = MyParams.Mobius.b; //  3./2.;
  RealD c   = MyParams.Mobius.c; //  1./2.;

  // These lines are unecessary if BC are all periodic
  std::cout << GridLogMessage << "boundary condition {1,1,1,-1}" << std::endl;
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);

  // ConjugateGradient<FermionField>  CG(MyParams.Mobius.StoppingCondition,MyParams.Mobius.MaxCGIterations);
  // // DJM: setup for EOFA ratio (Mobius)
  // MobiusEOFAFermionD Strange_Op_L(U, *FGrid, *FrbGrid, *GridPtr, *GridRBPtr, mass, mass, pv,  0.0, -1, M5, b, c, Params);
  // MobiusEOFAFermionD Strange_Op_R(U, *FGrid, *FrbGrid, *GridPtr, *GridRBPtr, pv,   mass, pv, -1.0,  1, M5, b, c, Params);
  // ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> EOFA(Strange_Op_L, Strange_Op_R, CG, OFRp, true);
  // MobiusEOFAFermionD Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , mass, mass, pv, 0.0, -1, M5, b, c);
  // MobiusEOFAFermionF Op_LF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, mass, mass, pv, 0.0, -1, M5, b, c);
  // MobiusEOFAFermionD Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , pv, mass,      pv, -1.0, 1, M5, b, c);
  // MobiusEOFAFermionF Op_RF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, pv, mass,      pv, -1.0, 1, M5, b, c);
  MobiusEOFAFermionD Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , mass, mass, pv, 0.0, -1, M5, b, c, Params);
  MobiusEOFAFermionF Op_LF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, mass, mass, pv, 0.0, -1, M5, b, c, Params);
  MobiusEOFAFermionD Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , pv, mass,      pv, -1.0, 1, M5, b, c, Params);
  MobiusEOFAFermionF Op_RF(UF, *FGridF, *FrbGridF, *GridPtrF, *GridRBPtrF, pv, mass,      pv, -1.0, 1, M5, b, c, Params);


  ConjugateGradient<FermionField>      ActionCG(ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DerivativeCG(DerivativeStoppingCondition,MaxCGIterations);

  const int MX_inner = 1000;
  // Mixed precision EOFA
  LinearOperatorEOFAD LinOp_L (Op_L);
  LinearOperatorEOFAD LinOp_R (Op_R);
  LinearOperatorEOFAF LinOp_LF(Op_LF);
  LinearOperatorEOFAF LinOp_RF(Op_RF);

  MxPCG_EOFA ActionCGL(ActionStoppingCondition,
                       MX_inner,
                       MaxCGIterations,
                       GridPtrF,
                       FrbGridF,
                       Op_LF,Op_L,
                       LinOp_LF,LinOp_L);

  MxPCG_EOFA DerivativeCGL(DerivativeStoppingCondition,
                           MX_inner,
                           MaxCGIterations,
                           GridPtrF,
                           FrbGridF,
                           Op_LF,Op_L,
                           LinOp_LF,LinOp_L);
  
  MxPCG_EOFA ActionCGR(ActionStoppingCondition,
                       MX_inner,
                       MaxCGIterations,
                       GridPtrF,
                       FrbGridF,
                       Op_RF,Op_R,
                       LinOp_RF,LinOp_R);
  
  MxPCG_EOFA DerivativeCGR(DerivativeStoppingCondition,
                           MX_inner,
                           MaxCGIterations,
                           GridPtrF,
                           FrbGridF,
                           Op_RF,Op_R,
                           LinOp_RF,LinOp_R);

  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> 
    EOFA(Op_L, Op_R, 
         ActionCG, 
         ActionCGL, ActionCGR,
         DerivativeCGL, DerivativeCGR,
         OFRp, true);

  //   FermionAction DenOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,mass,M5,b,c, Params);
  //   FermionAction NumOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,pv,  M5,b,c, Params);
  //   TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> Nf2a(NumOp, DenOp,CG,CG);

  // Set smearing (true/false), default: false
  EOFA.is_smeared = ApplySmearing;

  // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1); // fermion
  Level1.push_back(&EOFA);

  ActionLevel<HMCWrapper::Field> Level2(4); // gauge
  Level2.push_back(&Waction);
  //   Level2.push_back(&Syzaction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

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

/*
 * point2all_prop_dumper_full_claude.cc
 *
 * Dump the FULL Pauli-Dirac (PD) point-to-all quark propagator for the SU(Nc=4)
 * one-flavour SDM baryon programs: all 16 PD spin blocks q[s][s'] (LatticeColourMatrix,
 * row = sink colour, column = source colour) per source point, so that any local baryon
 * operator with arbitrary PD components at source and sink can be contracted downstream
 * (baryon_variants_corr_claude.cc: the 9 highest-weight B-form operators and their
 * J^P blocks; later the spin-general two-baryon operators).
 *
 * Relation to point2all_prop_dumper_claude.cc (the q00 dumper, left untouched): same
 * action, grids, split-grid MRHS batching, smearing, record format and atomic write.
 * The pre-rotated source trick is applied to all four PD source spins: source (s',c') =
 * (row s' of U) (x) e_{c'} in the Weyl basis, sink spin s = (1/2) sum_k U(s,k) x_k, so 16
 * solves per source point give the 16 blocks (4x the solves and 16x the storage of the
 * q00 dumper). Block (0,0) equals the q00 dumper's record bit for bit (same code path).
 *
 * Source points: the first --nsrc (default 1 = origin) of the 8 spatial corners
 * {0, L_mu/2}^3 at t = 0 (production: origin only; more corners later if wanted).
 *
 * Output: one .lime file, 16 single-precision LatticeColourMatrixF records per source
 * point, each tagged with a PropRecord {srcIndex = corner, store = "q<s><s'>", ...},
 * written corner-major, then s, then s'. Default name <config>.pdfull.lime.
 *
 * Multi-RHS via split-grid (P. Boyle et al., Grid, arXiv:1512.03487; template
 * Grid/tests/solver/Test_dwf_mrhs_cg.cc): the 16 columns of one source point are solved
 * in 16/nrhs batches; nrhs = prod(mpi_layout / mpi_split) must divide 16.
 *
 * --rottest: per source point, build the reference propagator from 16 unit Weyl solves
 * (identically smeared), rotate G^{PD} = (1/2) U G^{W} U^T, and compare every block with
 * the pre-rotated solve (identity for any gauge field; the cold field is the cheap
 * choice). See baryon_variants_impl_plan_claude.md.
 */

#include <Grid/Grid.h>
#include <algorithm>
#include <vector>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <filesystem>

#include "baryon_pd_common_claude.h"

using namespace std;
using namespace Grid;

int main(int argc, char ** argv)
{
  const int Ls   = 16;
  const int ncol = NPDSpin * Nc;   // 16 source columns (s', c') per source point

  Grid_init(&argc, &argv);

  // ------------------------------------------------------------------
  // Full grids.
  // ------------------------------------------------------------------
  Coordinate latt        = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt, simd_layout, mpi_layout);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  // ------------------------------------------------------------------
  // Configuration + action parameters (positional argv, as in the baryon codes):
  //   argv: <config> <M5> <mass> <outfile>
  // ------------------------------------------------------------------
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
    outfile = config + ".pdfull.lime";
  }

  // ------------------------------------------------------------------
  // Options: --split sx sy sz st, --width w, --niter N, --nsrc n, --rottest, --out f.
  // ------------------------------------------------------------------
  RealD smearWidth = 3.0;
  int   smearNiter = 40;
  int   nsrc       = 1;
  bool  rottest    = false;
  for(int i=0; i<argc; i++){
    if( std::string(argv[i]) == "--width" )   smearWidth = stod(argv[i+1]);
    if( std::string(argv[i]) == "--niter" )   smearNiter = std::stoi(argv[i+1]);
    if( std::string(argv[i]) == "--nsrc" )    nsrc = std::stoi(argv[i+1]);
    if( std::string(argv[i]) == "--rottest" ) rottest = true;
    if( std::string(argv[i]) == "--out" )     outfile = argv[i+1];
  }
  assert( nsrc >= 1 && nsrc <= 8 && "--nsrc must be in 1..8 (spatial corners)" );
  std::cout << GridLogMessage << "Gaussian source: width = " << smearWidth
            << " (non-overlap ceiling L/4 = " << latt[0]/4.0 << "), niter = " << smearNiter
            << " (stability N >~ 3 w^2 = " << 3.0*smearWidth*smearWidth << ")" << std::endl;
  std::cout << GridLogMessage << "source points: first " << nsrc << " corner(s) of {0,L/2}^3, t = 0" << std::endl;
  std::cout << GridLogMessage << "solves per source point: " << ncol << " (all PD source spins x colours)" << std::endl;

  Coordinate mpi_split;
  ParseSplit(argc, argv, mpi_layout, mpi_split);

  // ------------------------------------------------------------------
  // Mobius domain-wall action (scaled Shamir kernel, b+c=2, b-c=1), anti-periodic
  // in time {1,1,1,-1}; full action for the 4D<->5D import/export, split action
  // inside SplitGridSolver for the batched solve.
  // ------------------------------------------------------------------
  RealD b = 1.5;
  RealD c = 0.5;
  std::cout << GridLogMessage << "==============================================" << std::endl;
  std::cout << GridLogMessage << "MobiusFermion (scaled Shamir kernel), AP BC {1,1,1,-1}" << std::endl;
  std::cout << GridLogMessage << "==============================================" << std::endl;
  std::vector<Complex> boundary = {1,1,1,-1};
  MobiusFermionD::ImplParams Params(boundary);

  MobiusFermionD D(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, c, Params);

  SplitGridSolver solver(latt, simd_layout, mpi_layout, mpi_split, UGrid, FGrid, Ls,
                         Umu, mass, M5, b, c, Params, 1.0e-9, 100000);
  int nrhs = solver.nrhs;
  std::cout << GridLogMessage << "mpi_layout = "
            << mpi_layout[0] << "." << mpi_layout[1] << "."
            << mpi_layout[2] << "." << mpi_layout[3] << std::endl;
  std::cout << GridLogMessage << "mpi_split  = "
            << mpi_split[0] << "." << mpi_split[1] << "."
            << mpi_split[2] << "." << mpi_split[3] << std::endl;
  std::cout << GridLogMessage << "nrhs = " << nrhs << " (columns per source point = " << ncol << ")" << std::endl;
  assert( (ncol % nrhs) == 0 && "nrhs must divide the 16 columns of a source point" );
  int nbatch = ncol / nrhs;
  std::cout << GridLogMessage << "nbatch = " << nbatch << " split-solves per source point" << std::endl;

  // ------------------------------------------------------------------
  // Graceful wall blocker (same env knobs as the q00 dumper):
  //   PROP_DEADLINE_EPOCH : job deadline (epoch s); 0 / unset disables the check.
  //   PROP_TPT_SECONDS    : estimated per-config wall (s) for this layout.
  // Refuses to START when the remaining walltime cannot cover the run; the atomic
  // write below then guarantees no partial-looking <outfile> ever appears.
  // ------------------------------------------------------------------
  {
    long deadline = 0;
    if(const char* e = std::getenv("PROP_DEADLINE_EPOCH")) deadline = std::atol(e);
    double tpt = 0.0;
    if(const char* e = std::getenv("PROP_TPT_SECONDS")) tpt = std::atof(e);
    const double margin = 1.2;
    if(deadline > 0 && tpt > 0.0){
      uint64_t stop = 0;
      if(UGrid->IsBoss()){
        if((double)std::time(nullptr) + margin*tpt > (double)deadline) stop = 1;
      }
      UGrid->GlobalSum(stop);
      if(stop){
        std::cout << GridLogMessage << "blocker: est " << (long)(margin*tpt)
                  << "s to finish exceeds deadline; stopping gracefully WITHOUT writing "
                  << outfile << std::endl;
        Grid_finalize();
        return 0;
      }
    }
  }

  // ------------------------------------------------------------------
  // Source points and smearing links.
  // ------------------------------------------------------------------
  std::vector<Coordinate> srcPts;
  CornerSourcePoints(latt, nsrc, srcPts);
  for(int i=0; i<nsrc; i++){
    std::cout << GridLogMessage << "source " << i << " = ("
              << srcPts[i][0] << "," << srcPts[i][1] << ","
              << srcPts[i][2] << "," << srcPts[i][3] << ")" << std::endl;
  }
  std::vector<LatticeColourMatrix> U(Nd, UGrid);
  SpatialLinks(Umu, U);

  // ------------------------------------------------------------------
  // Output: single-precision SciDAC/LIME, atomic (<outfile>.inprogress, renamed on
  // completion). The writer stays open across the source-point loop so each point's
  // 16 blocks are written as soon as they are solved (memory = one point's blocks).
  // ------------------------------------------------------------------
  GridCartesian * UGridF = SpaceTimeGrid::makeFourDimGrid(latt,
                             GridDefaultSimd(Nd, vComplexF::Nsimd()), mpi_layout);
  const std::string tmpfile = outfile + ".inprogress";
  std::cout << GridLogMessage << "Writing " << nsrc * NPDSpin * NPDSpin
            << " PD block records (single precision) to " << tmpfile
            << " (atomic; rename to " << outfile << " on completion)" << std::endl;
  ScidacWriter WR(UGridF->IsBoss());
  WR.open(tmpfile);

  // Per-point block storage q[s][s'] and reusable buffers.
  std::vector<std::vector<LatticeColourMatrix> > q(NPDSpin, std::vector<LatticeColourMatrix>(NPDSpin, UGrid));
  std::vector<LatticeFermion> src4(nrhs, UGrid);
  std::vector<LatticeFermion> res4(nrhs, UGrid);
  LatticeColourVector cv(UGrid);

  // Full-grid solver for the --rottest reference (not split, for simplicity).
  ConjugateGradient<LatticeFermion>            CGref(1.0e-9, 100000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion>  schurRef(CGref);
  ZeroGuesser<LatticeFermion>                  ZGref;
  RealD rotMaxDiff = 0.0;

  for(int i=0; i<nsrc; i++){
    std::cout << GridLogMessage << "=== source point " << i << " ===" << std::endl;
    for(int s=0; s<NPDSpin; s++){
      for(int sp=0; sp<NPDSpin; sp++){
        q[s][sp] = Zero();
      }
    }

    // ----------------------------------------------------------------
    // 16 pre-rotated solves, batched nrhs at a time. Column index col = 4 s' + c'.
    // ----------------------------------------------------------------
    for(int bt=0; bt<nbatch; bt++){
      int col0 = bt * nrhs;
      for(int r=0; r<nrhs; r++){
        int col    = col0 + r;
        int sprime = col / Nc;
        int cprime = col % Nc;
        MakePreRotatedSmearedSource(srcPts[i], sprime, cprime, U, smearWidth, smearNiter, src4[r]);
      }
      solver.Solve(D, src4, res4);
      for(int r=0; r<nrhs; r++){
        int col    = col0 + r;
        int sprime = col / Nc;
        int cprime = col % Nc;
        for(int s=0; s<NPDSpin; s++){
          PDSinkBlock(res4[r], s, cv);
          PokeColourColumn(cv, cprime, q[s][sprime]);
        }
        std::cout << GridLogMessage << "  point " << i << " column (s'=" << sprime
                  << ", c'=" << cprime << ") norm2(res4) = " << norm2(res4[r]) << std::endl;
      }
    }
    for(int s=0; s<NPDSpin; s++){
      for(int sp=0; sp<NPDSpin; sp++){
        std::cout << GridLogMessage << "  q[" << s << "][" << sp << "] norm2 = " << norm2(q[s][sp]) << std::endl;
      }
    }

    // ----------------------------------------------------------------
    // --rottest: reference from 16 unit Weyl solves, rotated (1/2) U G U^T.
    // ----------------------------------------------------------------
    if( rottest ){
      std::cout << GridLogMessage << "  rotation self-test for point " << i << std::endl;
      SpinMatrix Urot = WeylToPauliDiracU();
      SpinMatrix Udag = transpose(Urot);
      LatticeFermion    usrc4(UGrid);
      LatticeFermion    usrc5(FGrid);
      LatticeFermion    ures5(FGrid);
      LatticeFermion    ures4(UGrid);
      LatticePropagator Gref(UGrid);
      for(int sw=0; sw<Ns; sw++){
        for(int cc=0; cc<Nc; cc++){
          MakeSmearedUnitSource(srcPts[i], sw, cc, U, smearWidth, smearNiter, usrc4);
          D.ImportPhysicalFermionSource(usrc4, usrc5);
          ures5 = Zero();
          schurRef(D, usrc5, ures5, ZGref);
          D.ExportPhysicalFermionSolution(ures5, ures4);
          FermToProp<MobiusFermionD>(Gref, ures4, sw, cc);
        }
      }
      LatticePropagator Gref_PD = 0.5 * Urot * Gref * Udag;
      for(int s=0; s<NPDSpin; s++){
        for(int sp=0; sp<NPDSpin; sp++){
          LatticeColourMatrix ref  = peekSpin(Gref_PD, s, sp);
          LatticeColourMatrix diff = ref - q[s][sp];
          RealD d = norm2(diff);
          rotMaxDiff = std::max(rotMaxDiff, d);
          std::cout << GridLogMessage << "    block (" << s << "," << sp << ")"
                    << " norm2(ref) = " << norm2(ref)
                    << " norm2(prerot) = " << norm2(q[s][sp])
                    << " norm2(diff) = " << d << std::endl;
        }
      }
    }

    // ----------------------------------------------------------------
    // Write the 16 blocks of this point (s outer, s' inner).
    // ----------------------------------------------------------------
    for(int s=0; s<NPDSpin; s++){
      for(int sp=0; sp<NPDSpin; sp++){
        LatticeColourMatrixF qF(UGridF);
        precisionChange(qF, q[s][sp]);

        PropRecord rec;
        rec.config     = config;
        rec.srcCoord   = std::vector<int>({srcPts[i][0], srcPts[i][1], srcPts[i][2], srcPts[i][3]});
        rec.srcIndex   = i;
        rec.mass       = mass;
        rec.M5         = M5;
        rec.b          = b;
        rec.c          = c;
        rec.Ls         = Ls;
        rec.smearWidth = smearWidth;
        rec.smearNiter = smearNiter;
        rec.store      = PDBlockTag(s, sp);
        rec.prec       = "single";

        WR.writeScidacFieldRecord(qF, rec);
      }
    }
    std::cout << GridLogMessage << "  point " << i << ": 16 records written." << std::endl;
  }

  if( rottest ){
    std::cout << GridLogMessage << "rotation self-test: max norm2(diff) over all blocks and points = "
              << rotMaxDiff << std::endl;
    assert( rotMaxDiff < 1.0e-8 && "rotation self-test failed: pre-rotated blocks != reference" );
    std::cout << GridLogMessage << "rotation self-test PASSED." << std::endl;
  }

  WR.close();
  UGridF->Barrier();
  if(UGridF->IsBoss()) std::filesystem::rename(tmpfile, outfile);
  UGridF->Barrier();
  std::cout << GridLogMessage << "full PD propagator dumped to " << outfile << std::endl;

  Grid_finalize();
}

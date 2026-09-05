/*
 * point2all_prop_dumper_claude.cc
 *
 * Dump point-to-all quark propagators for the SU(Nc=4) one-flavour SDM two-baryon
 * GEVP program. Sources sit on the 8 spatial corners {0, L_mu/2}^3 at t_src = 0;
 * the propagators are solved once per configuration and serialised (SciDAC/LIME,
 * single precision) for reuse by downstream two-baryon contractions.
 *
 * v1: q00 mode only. Only the Pauli-Dirac source-spin-0 component is needed for the
 * "0000" baryon, so we pre-rotate the source (Weyl spin (e_0 + e_2) per colour) and
 * solve 4 times per source point instead of 16; the stored object is the PD spin-0
 * colour matrix q00 (16x smaller than the full propagator). See
 * point2all_prop_dumper_impl_plan_claude.md for the derivation and the GEVP framing.
 *
 * Multi-RHS via split-grid: the 8 sources are solved concurrently on MPI
 * sub-communicators (Grid_split / Grid_unsplit). This is the Grid idiom for batching
 * many propagator solves on GPUs.
 *   Method / reference: P. Boyle et al., Grid, arXiv:1512.03487; usage template
 *   Grid/tests/solver/Test_dwf_mrhs_cg.cc.
 *
 * Chunk 1 (this commit): grids + parameters + split-grid setup (full and split
 * Mobius actions, gauge field replicated to the sub-grids). Sources, MRHS solve,
 * time-window extraction and SciDAC output follow in later chunks.
 */

#include <Grid/Grid.h>
#include <algorithm>
#include <vector>
#include <sstream>

using namespace std;
using namespace Grid;

// HDF5/SciDAC-serialisable per-record metadata attached to each q00 field record.
class PropRecord: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(PropRecord,
    std::string,      config,
    std::vector<int>, srcCoord,
    int,              srcIndex,
    double,           mass,
    double,           M5,
    double,           b,
    double,           c,
    int,              Ls,
    double,           smearWidth,
    int,              smearNiter,
    std::string,      store,
    std::string,      prec);
};

// --------------------------------------------------------------------------
// Pre-rotated, Gaussian-smeared point source for q00 mode.
//
// The "0000" baryon needs only the Pauli-Dirac source-spin-0 component. Column
// (0,c') of S^{PD} = (1/2) U S^{W} Udagger is obtained by solving with the Weyl
// source Udagger e_{(0,c')} = (row 0 of U) \otimes e_{c'} = (e_0 + e_2) \otimes e_{c'}
// (from WeylToPauliDiracU in two_baryon_corr_prod_claude.cc). So the 4D source is a
// point at \site with spinor components {0,2} set in colour c'; a single covariant
// (Wuppertal) Gaussian smear shapes the spatial profile. The 1/2 and the sink-spin-0
// projection are applied when q00 is assembled (Chunk 3).
//
// width w sets the radius (rms ~ 1.22 w in 3D); non-overlap wants w <~ L/4. Niter N
// is a stability/convergence knob, N >~ 3 w^2. Smearing (spatial, orthog = Tdir) is
// gauge covariant, so it needs no gauge fixing.
static void MakePreRotatedSmearedSource(const Coordinate &site, int cprime,
                                        const std::vector<LatticeColourMatrix> &U,
                                        RealD width, int niter,
                                        LatticeFermion &src)
{
  src = Zero();
  SpinColourVector scv;
  scv = Zero();
  scv()(0)(cprime) = 1.0;
  scv()(2)(cprime) = 1.0;
  pokeSite(scv, src, site);

  if( niter > 0 && width > 0.0 ){
    CovariantSmearing<PeriodicGimplD>::GaussianSmear(U, src, width, niter, Tdir);
  }
}

// Single spin-colour unit source at \site (identical smearing), for the rotation
// self-test reference (full 16-component propagator).
static void MakeSmearedUnitSource(const Coordinate &site, int s, int cc,
                                  const std::vector<LatticeColourMatrix> &U,
                                  RealD width, int niter,
                                  LatticeFermion &src)
{
  src = Zero();
  SpinColourVector scv;
  scv = Zero();
  scv()(s)(cc) = 1.0;
  pokeSite(scv, src, site);

  if( niter > 0 && width > 0.0 ){
    CovariantSmearing<PeriodicGimplD>::GaussianSmear(U, src, width, niter, Tdir);
  }
}

// Weyl (chiral) -> Pauli-Dirac rotation U, S^{PD} = (1/2) U S^{W} Udagger (same as
// two_baryon_corr_prod_claude.cc:293). Row 0 = (e_0 + e_2), the basis of the q00 trick.
static SpinMatrix WeylToPauliDiracU()
{
  SpinMatrix U = Zero();
  U()(0,0) =  1.0;
  U()(1,1) =  1.0;
  U()(2,0) = -1.0;
  U()(3,1) = -1.0;
  U()(0,2) =  1.0;
  U()(1,3) =  1.0;
  U()(2,2) =  1.0;
  U()(3,3) =  1.0;
  return U;
}

int main(int argc, char ** argv)
{
  const int Ls    = 16;
  const int n_src = 8; // 8 spatial corners {0, L_mu/2}^3 at t_src = 0

  Grid_init(&argc, &argv);

  // ------------------------------------------------------------------
  // Full grids (the physical problem lives here).
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
    outfile = config + ".prop.lime";
  }

  // ------------------------------------------------------------------
  // Split-grid layout: mpi_split = sub-communicator geometry (each sub-grid holds
  // the full lattice on prod(mpi_split) ranks). nrhs = prod(mpi_layout/mpi_split)
  // sources are solved concurrently. Parsed from "--split sx sy sz st"; default
  // 1^Nd (no split, nrhs = prod(mpi_layout)).
  // ------------------------------------------------------------------
  Coordinate mpi_split(mpi_layout.size(), 1);
  for(int i=0; i<argc; i++){
    if( std::string(argv[i]) == "--split" ){
      for(int k=0; k<mpi_layout.size(); k++){
        std::stringstream ss;
        ss << argv[i+1+k];
        ss >> mpi_split[k];
      }
      break;
    }
  }

  int nrhs = 1;
  for(int i=0; i<mpi_layout.size(); i++) nrhs *= (mpi_layout[i] / mpi_split[i]);

  std::cout << GridLogMessage << "mpi_layout = "
            << mpi_layout[0] << "." << mpi_layout[1] << "."
            << mpi_layout[2] << "." << mpi_layout[3] << std::endl;
  std::cout << GridLogMessage << "mpi_split  = "
            << mpi_split[0] << "." << mpi_split[1] << "."
            << mpi_split[2] << "." << mpi_split[3] << std::endl;
  std::cout << GridLogMessage << "nrhs = " << nrhs << " (n_src = " << n_src << ")" << std::endl;
  // Sources are solved in n_src/nrhs split-solve batches; nrhs must divide n_src.
  // nrhs = n_src (=8) is the production case (all sources concurrent); nrhs = 1 is a
  // single-rank sequential run (local validation).
  assert( (n_src % nrhs) == 0 && "nrhs must divide the number of sources (8)" );
  int nbatch = n_src / nrhs;
  std::cout << GridLogMessage << "nbatch = " << nbatch << " (sources per split-solve = " << nrhs << ")" << std::endl;

  // Split grids (the "S" family). The GridCartesian split constructor returns this
  // rank's sub-grid index in "me".
  int me;
  GridCartesian         * SGrid   = new GridCartesian(latt, simd_layout, mpi_split, *UGrid, me);
  GridRedBlackCartesian * SUrbGrid= SpaceTimeGrid::makeFourDimRedBlackGrid(SGrid);
  GridCartesian         * SFGrid  = SpaceTimeGrid::makeFiveDimGrid(Ls, SGrid);
  GridRedBlackCartesian * SFrbGrid= SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, SGrid);

  // ------------------------------------------------------------------
  // Mobius domain-wall action (scaled Shamir kernel, b+c=2, b-c=1), anti-periodic
  // in time {1,1,1,-1}. Built on both the full and the split grids: the full action
  // does the 4D<->5D physical import/export, the split action does the batched solve.
  // ------------------------------------------------------------------
  RealD b = 1.5;
  RealD c = 0.5;
  std::cout << GridLogMessage << "==============================================" << std::endl;
  std::cout << GridLogMessage << "MobiusFermion (scaled Shamir kernel), AP BC {1,1,1,-1}" << std::endl;
  std::cout << GridLogMessage << "==============================================" << std::endl;
  std::vector<Complex> boundary = {1,1,1,-1};
  MobiusFermionD::ImplParams Params(boundary);

  MobiusFermionD D(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, b, c, Params);

  // Replicate the gauge field onto every sub-grid, then build the split action.
  LatticeGaugeField s_Umu(SGrid);
  Grid_split(Umu, s_Umu);
  MobiusFermionD sD(s_Umu, *SFGrid, *SFrbGrid, *SGrid, *SUrbGrid, mass, M5, b, c, Params);

  std::cout << GridLogMessage << "Chunk 1 complete: grids + split-grid actions ready." << std::endl;
  std::cout << GridLogMessage << "  norm2(Umu)   = " << norm2(Umu)   << std::endl;
  std::cout << GridLogMessage << "  norm2(s_Umu) = " << norm2(s_Umu) << std::endl;

  // ------------------------------------------------------------------
  // Chunk 2: sources.
  //   (a) 8 spatial-corner source points {0, L_mu/2}^3 at t_src = 0.
  //   (b) spatial gauge links for covariant Gaussian smearing.
  //   (c) pre-rotated (PD source-spin-0) smeared source builder (above); a sanity
  //       pass builds the 8 sources for colour 0 and prints their norms.
  // ------------------------------------------------------------------

  // Smearing knobs: --width <w> (radius, ~L/4 non-overlap ceiling), --niter <N>
  // (stability/convergence, N >~ 3 w^2). Placeholder defaults; set from t_0 scale.
  RealD smearWidth = 3.0;
  int   smearNiter = 40;
  bool  rottest    = false;
  for(int i=0; i<argc; i++){
    if( std::string(argv[i]) == "--width" )   smearWidth = stod(argv[i+1]);
    if( std::string(argv[i]) == "--niter" )   smearNiter = std::stoi(argv[i+1]);
    if( std::string(argv[i]) == "--rottest" ) rottest = true;
  }
  std::cout << GridLogMessage << "Gaussian source: width = " << smearWidth
            << " (non-overlap ceiling L/4 = " << latt[0]/4.0 << "), niter = " << smearNiter
            << " (stability N >~ 3 w^2 = " << 3.0*smearWidth*smearWidth << ")" << std::endl;

  // (a) 8 corners: spatial coords from the low 3 bits of the source index, t = 0.
  std::vector<Coordinate> srcPts(n_src, Coordinate(Nd, 0));
  for(int i=0; i<n_src; i++){
    for(int mu=0; mu<3; mu++){
      int bit = (i >> mu) & 1;
      srcPts[i][mu] = bit * (latt[mu] / 2);
    }
    srcPts[i][Tdir] = 0;
    std::cout << GridLogMessage << "source " << i << " = ("
              << srcPts[i][0] << "," << srcPts[i][1] << ","
              << srcPts[i][2] << "," << srcPts[i][3] << ")" << std::endl;
  }

  // (b) spatial gauge links (mu = Tdir is ignored by GaussianSmear via orthog).
  std::vector<LatticeColourMatrix> U(Nd, UGrid);
  for(int mu=0; mu<Nd; mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
  }

  // (c) sanity: build the 8 colour-0 sources, report norms (verifies smearing runs).
  {
    LatticeFermion src(UGrid);
    for(int i=0; i<n_src; i++){
      MakePreRotatedSmearedSource(srcPts[i], 0, U, smearWidth, smearNiter, src);
      std::cout << GridLogMessage << "  src(corner " << i << ", colour 0) norm2 = "
                << norm2(src) << std::endl;
    }
  }
  std::cout << GridLogMessage << "Chunk 2 complete: sources ready." << std::endl;

  // ------------------------------------------------------------------
  // Chunk 3: MRHS split-grid solve, q00 assembly.
  //   For each of the Nc source colours c', the 8 pre-rotated sources are solved
  //   concurrently on the split grid (Grid_split -> solve on sD -> Grid_unsplit).
  //   q00[i](c,c') = (1/2)(x_{spin0,c} + x_{spin2,c}) is the PD spin-0 colour matrix
  //   (row = sink colour c, col = source colour c') of source point i.
  //   Nc = 4 split-solves total.
  // ------------------------------------------------------------------
  ConjugateGradient<LatticeFermion>            CG(1.0e-9, 100000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion>  schur(CG);
  ZeroGuesser<LatticeFermion>                  ZG;

  std::vector<LatticeColourMatrix> q00(n_src, UGrid);
  for(int i=0; i<n_src; i++) q00[i] = Zero();

  // Reusable buffers.
  LatticeFermion src4(UGrid);
  LatticeFermion res4(UGrid);
  std::vector<LatticeFermion> f5  (nrhs, FGrid);    // imported 5D sources for one batch
  std::vector<LatticeFermion> sol5(nrhs, FGrid);    // unsplit 5D solutions for one batch
  LatticeFermion s_src(SFGrid);                     // split-grid source (per sub-grid)
  LatticeFermion s_res(SFGrid);                     // split-grid solution

  for(int cprime=0; cprime<Nc; cprime++){
    std::cout << GridLogMessage << "=== source colour c' = " << cprime << " ===" << std::endl;

    for(int b=0; b<nbatch; b++){
      int i0 = b * nrhs; // global source index of this batch's first source

      // Build + import the nrhs pre-rotated sources of this batch for colour c'.
      for(int r=0; r<nrhs; r++){
        MakePreRotatedSmearedSource(srcPts[i0+r], cprime, U, smearWidth, smearNiter, src4);
        D.ImportPhysicalFermionSource(src4, f5[r]);
      }

      // Distribute across sub-communicators, solve concurrently, gather back.
      Grid_split(f5, s_src);
      s_res = Zero();
      schur(sD, s_src, s_res, ZG);
      Grid_unsplit(sol5, s_res);
      std::cout << GridLogMessage << "  batch " << b
                << " norm2(s_src) = " << norm2(s_src)
                << " norm2(s_res) = " << norm2(s_res) << std::endl;

      // Export to 4D and assemble the c'-th column of each q00 in the batch.
      for(int r=0; r<nrhs; r++){
        int i = i0 + r;
        D.ExportPhysicalFermionSolution(sol5[r], res4);
        LatticeColourVector cv0 = peekSpin(res4, 0);
        LatticeColourVector cv2 = peekSpin(res4, 2);
        for(int c=0; c<Nc; c++){
          LatticeComplex e = 0.5 * (peekColour(cv0, c) + peekColour(cv2, c));
          pokeColour(q00[i], e, c, cprime);
        }
        std::cout << GridLogMessage << "  corner " << i << " c'=" << cprime
                  << " norm2(res4) = " << norm2(res4) << std::endl;
      }
    }
  }

  for(int i=0; i<n_src; i++){
    std::cout << GridLogMessage << "q00[" << i << "] norm2 = " << norm2(q00[i]) << std::endl;
  }
  std::cout << GridLogMessage << "Chunk 3 complete: q00 colour matrices assembled." << std::endl;

  // ------------------------------------------------------------------
  // Chunk 6 (optional, --rottest): validate the pre-rotation / 4-solve trick.
  //   Reference q00 from a full 16-component propagator with the SAME (smeared)
  //   source, rotated the reference way q00_ref = peekSpin(0.5 U S Udagger, 0,0),
  //   compared to the 4-solve q00. The identity holds for any gauge field; the
  //   free (cold) field is just the cheap, clean choice. Full-grid solves (not
  //   split) for simplicity.
  // ------------------------------------------------------------------
  if( rottest ){
    std::cout << GridLogMessage << "=== rotation self-test (--rottest) ===" << std::endl;
    SpinMatrix Urot  = WeylToPauliDiracU();
    SpinMatrix Udag  = transpose(Urot);

    LatticeFermion    usrc4(UGrid);
    LatticeFermion    usrc5(FGrid);
    LatticeFermion    ures5(FGrid);
    LatticeFermion    ures4(UGrid);
    RealD maxdiff = 0.0;

    for(int i=0; i<n_src; i++){
      // Full 16-component Weyl propagator from the identically-smeared source.
      LatticePropagator Sref(UGrid);
      for(int s=0; s<Nd; s++){
        for(int cc=0; cc<Nc; cc++){
          MakeSmearedUnitSource(srcPts[i], s, cc, U, smearWidth, smearNiter, usrc4);
          D.ImportPhysicalFermionSource(usrc4, usrc5);
          ures5 = Zero();
          schur(D, usrc5, ures5, ZG);
          D.ExportPhysicalFermionSolution(ures5, ures4);
          FermToProp<MobiusFermionD>(Sref, ures4, s, cc);
        }
      }
      // Reference rotation, then PD spin-(0,0) colour matrix.
      LatticePropagator   Sref_PD = 0.5 * Urot * Sref * Udag;
      LatticeColourMatrix q00_ref = peekSpin(Sref_PD, 0, 0);

      LatticeColourMatrix diff = q00_ref - q00[i];
      RealD d = norm2(diff);
      maxdiff = std::max(maxdiff, d);
      std::cout << GridLogMessage << "  corner " << i
                << " norm2(q00_ref) = " << norm2(q00_ref)
                << " norm2(q00_4solve) = " << norm2(q00[i])
                << " norm2(diff) = " << d << std::endl;
    }
    std::cout << GridLogMessage << "rotation self-test: max norm2(diff) = "
              << maxdiff << std::endl;
    assert( maxdiff < 1.0e-8 && "rotation self-test failed: 4-solve q00 != reference" );
    std::cout << GridLogMessage << "rotation self-test PASSED." << std::endl;
  }

  // ------------------------------------------------------------------
  // Chunk 4+5: single-precision demote and SciDAC/LIME output.
  //   v1 stores the full time extent of q00 (no plateau window; q00 files are
  //   small, and localCopyRegion cannot repartition a window across the time
  //   decomposition -- select the plateau downstream at read time).
  //   One .lime file, n_src field records, each a LatticeColourMatrixF (single
  //   precision) tagged with a PropRecord (source coord + action parameters).
  // ------------------------------------------------------------------
  GridCartesian * UGridF = SpaceTimeGrid::makeFourDimGrid(latt,
                             GridDefaultSimd(Nd, vComplexF::Nsimd()), mpi_layout);

  std::cout << GridLogMessage << "Writing " << n_src
            << " q00 records (single precision) to " << outfile << std::endl;

  ScidacWriter WR(UGridF->IsBoss());
  WR.open(outfile);
  for(int i=0; i<n_src; i++){
    LatticeColourMatrixF q00F(UGridF);
    precisionChange(q00F, q00[i]);

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
    rec.store      = "q00";
    rec.prec       = "single";

    WR.writeScidacFieldRecord(q00F, rec);
    std::cout << GridLogMessage << "  record " << i << " (corner "
              << srcPts[i][0] << "," << srcPts[i][1] << ","
              << srcPts[i][2] << "," << srcPts[i][3] << ") written." << std::endl;
  }
  WR.close();

  std::cout << GridLogMessage << "Chunk 4+5 complete: q00 dumped to " << outfile << std::endl;

  Grid_finalize();
}

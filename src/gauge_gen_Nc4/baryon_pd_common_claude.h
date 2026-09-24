/*
 * baryon_pd_common_claude.h
 *
 * Shared pieces of the SU(Nc=4), one-flavour SDM baryon programs that work in the
 * Pauli-Dirac (PD) spin basis:
 *
 *   PropRecord               per-record metadata of the SciDAC/LIME propagator dumps
 *                            (field order = on-disk format; identical to
 *                            point2all_prop_dumper_claude.cc, never reorder)
 *   CorrFile                 one correlator series C(t) for Hdf5Writer
 *   WeylToPauliDiracU        Weyl (chiral, Grid default) -> PD rotation U,
 *                            G^{PD} = (1/2) U G^{W} U^T
 *   MakePreRotatedSmearedSource   point source for PD source spin s' (row s' of U) in
 *                            colour c', optionally covariant-Gaussian smeared
 *   MakeSmearedUnitSource    single Weyl spin-colour unit source (rotation self-test)
 *   PDSinkBlock              PD sink spin s of a Weyl solution: (1/2) sum_k U(s,k) x_k
 *   CornerSourcePoints       the spatial corners {0, L_mu/2}^3 at t = 0
 *   SpatialLinks             gauge links for CovariantSmearing
 *   SplitGridSolver          split-grid multi-RHS Mobius solve (Grid_split / Grid_unsplit)
 *
 * Conventions (baryon_theory_claude.md Sec. 4.0): PD index 0 = upper spin up, 1 = upper
 * spin down (parity +), 2 = lower spin up, 3 = lower spin down (parity -);
 * gamma_4^{PD} = diag(1,1,-1,-1). Grid's chiral basis has gamma_5 = diag(1,1,-1,-1) and
 * gamma_4 swapping 0<->2, 1<->3; U below maps them onto the PD forms.
 *
 * Pre-rotated source trick (point2all_prop_dumper_claude.cc, generalised to all four PD
 * source spins): column (s',c') of G^{PD} = (1/2) U G^{W} U^T is obtained by solving with
 * the Weyl source sum_l U(s',l) e_l (x) e_{c'} (row s' of U); the PD sink spin s is then
 * (1/2) sum_k U(s,k) x_k of the Weyl solution x. 16 solves per source point give all
 * 16 blocks q[s][s'] (LatticeColourMatrix, row = sink colour, column = source colour).
 *
 * Multi-RHS: Grid split-grid idiom, P. Boyle et al., arXiv:1512.03487,
 * template Grid/tests/solver/Test_dwf_mrhs_cg.cc.
 *
 * This header expects to be included after <Grid/Grid.h> by single-translation-unit
 * drivers that already use namespace Grid.
 */

#pragma once

#include <Grid/Grid.h>
#include <vector>
#include <string>
#include <memory>
#include <sstream>

using namespace Grid;

// --------------------------------------------------------------------------
// Serialisable records
// --------------------------------------------------------------------------

// Per-record metadata attached to each propagator field record (SciDAC/LIME).
// store = "q<s><s'>" names the PD (sink spin s, source spin s') colour-matrix block;
// the q00 dumper writes only store = "q00".
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

// One correlator series C(t) (Hdf5Writer writes group <name>/data as compound (re,im)).
class CorrFile: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CorrFile, std::vector<Complex>, data);
};

// --------------------------------------------------------------------------
// Weyl -> Pauli-Dirac rotation
// --------------------------------------------------------------------------

// G^{PD} = (1/2) U G^{W} U^T  (same U as two_baryon_corr_prod_claude.cc:293).
// Rows of U: (1,0,1,0), (0,1,0,1), (-1,0,1,0), (0,-1,0,1); U U^T = 2.
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

// Number of PD spin blocks per source point.
static const int NPDSpin = Ns;

// Record tag of block (s, s'): "q00", "q01", ..., "q33".
static std::string PDBlockTag(int s, int sprime)
{
  std::string tag = "q";
  tag += std::to_string(s);
  tag += std::to_string(sprime);
  return tag;
}

// --------------------------------------------------------------------------
// Sources
// --------------------------------------------------------------------------

// Pre-rotated, optionally Gaussian-smeared point source for PD source spin s' in
// colour c': Weyl spinor components l set to U(s',l) at \site. Smearing (spatial,
// orthog = Tdir) is gauge covariant, so it needs no gauge fixing. width w sets the
// radius (rms ~ 1.22 w in 3D, non-overlap w <~ L/4); niter N >~ 3 w^2 for stability.
static void MakePreRotatedSmearedSource(const Coordinate &site, int sprime, int cprime,
                                        const std::vector<LatticeColourMatrix> &U,
                                        RealD width, int niter,
                                        LatticeFermion &src)
{
  SpinMatrix Urot = WeylToPauliDiracU();
  src = Zero();
  SpinColourVector scv;
  scv = Zero();
  for(int l=0; l<Ns; l++){
    scv()(l)(cprime) = TensorRemove(Urot()(sprime,l));
  }
  pokeSite(scv, src, site);

  if( niter > 0 && width > 0.0 ){
    CovariantSmearing<PeriodicGimplD>::GaussianSmear(U, src, width, niter, Tdir);
  }
}

// Single Weyl spin-colour unit source at \site (identical smearing), for the rotation
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

// PD sink spin s of a Weyl solution x: cv = (1/2) sum_k U(s,k) x_k (a colour vector
// field). For s = 0 this is the q00 dumper's (1/2)(x_0 + x_2).
static void PDSinkBlock(const LatticeFermion &x, int s, LatticeColourVector &cv)
{
  SpinMatrix Urot = WeylToPauliDiracU();
  cv = Zero();
  for(int k=0; k<Ns; k++){
    Complex u = TensorRemove(Urot()(s,k));
    if( u != Complex(0.0, 0.0) ){
      LatticeColourVector xk = peekSpin(x, k);
      cv = cv + u * xk;
    }
  }
  cv = 0.5 * cv;
}

// Insert the sink colour vector cv as column c' of the colour matrix q.
static void PokeColourColumn(const LatticeColourVector &cv, int cprime, LatticeColourMatrix &q)
{
  for(int c=0; c<Nc; c++){
    LatticeComplex e = peekColour(cv, c);
    pokeColour(q, e, c, cprime);
  }
}

// The first n of the 8 spatial corners {0, L_mu/2}^3 at t = 0; corner index bits =
// (x, y, z) halves (0 = origin, 7 = (L/2, L/2, L/2)).
static void CornerSourcePoints(const Coordinate &latt, int n, std::vector<Coordinate> &pts)
{
  pts.assign(n, Coordinate(Nd, 0));
  for(int i=0; i<n; i++){
    for(int mu=0; mu<3; mu++){
      int bit = (i >> mu) & 1;
      pts[i][mu] = bit * (latt[mu] / 2);
    }
    pts[i][Tdir] = 0;
  }
}

// Gauge links per direction (mu = Tdir is ignored by GaussianSmear via orthog).
static void SpatialLinks(const LatticeGaugeField &Umu, std::vector<LatticeColourMatrix> &U)
{
  for(int mu=0; mu<Nd; mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
  }
}

// --------------------------------------------------------------------------
// Split-grid multi-RHS Mobius solver
// --------------------------------------------------------------------------

// Holds the split grids, the replicated gauge field, the split action and the batch
// buffers. Solve() takes nrhs 4D sources and returns nrhs 4D solutions, solved
// concurrently on the MPI sub-communicators (nrhs = prod(mpi_layout / mpi_split)).
// The full action D (on the full grids) does the 4D<->5D import/export.
class SplitGridSolver {
public:
  int nrhs;
  Coordinate mpi_split;
  GridCartesian         *SGrid;
  GridRedBlackCartesian *SUrbGrid;
  GridCartesian         *SFGrid;
  GridRedBlackCartesian *SFrbGrid;
  std::unique_ptr<LatticeGaugeField> s_Umu;
  std::unique_ptr<MobiusFermionD>    sD;
  std::unique_ptr<ConjugateGradient<LatticeFermion> >           CG;
  std::unique_ptr<SchurRedBlackDiagMooeeSolve<LatticeFermion> > schur;
  ZeroGuesser<LatticeFermion> ZG;
  std::vector<LatticeFermion> f5;      // imported 5D sources of one batch (full grid)
  std::vector<LatticeFermion> sol5;    // unsplit 5D solutions of one batch (full grid)
  std::unique_ptr<LatticeFermion> s_src;   // split-grid source (per sub-grid)
  std::unique_ptr<LatticeFermion> s_res;   // split-grid solution

  SplitGridSolver(const Coordinate &latt, const Coordinate &simd_layout,
                  const Coordinate &mpi_layout, const Coordinate &split,
                  GridCartesian *UGrid, GridCartesian *FGrid, int Ls,
                  const LatticeGaugeField &Umu,
                  RealD mass, RealD M5, RealD b, RealD c,
                  const MobiusFermionD::ImplParams &Params,
                  RealD tol, int maxit)
  {
    mpi_split = split;
    nrhs = 1;
    for(int i=0; i<mpi_layout.size(); i++){
      nrhs *= (mpi_layout[i] / mpi_split[i]);
    }

    int me;
    SGrid    = new GridCartesian(latt, simd_layout, mpi_split, *UGrid, me);
    SUrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(SGrid);
    SFGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, SGrid);
    SFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, SGrid);

    s_Umu.reset(new LatticeGaugeField(SGrid));
    Grid_split(const_cast<LatticeGaugeField&>(Umu), *s_Umu);
    sD.reset(new MobiusFermionD(*s_Umu, *SFGrid, *SFrbGrid, *SGrid, *SUrbGrid, mass, M5, b, c, Params));

    CG.reset(new ConjugateGradient<LatticeFermion>(tol, maxit));
    schur.reset(new SchurRedBlackDiagMooeeSolve<LatticeFermion>(*CG));

    f5.resize(nrhs, LatticeFermion(FGrid));
    sol5.resize(nrhs, LatticeFermion(FGrid));
    s_src.reset(new LatticeFermion(SFGrid));
    s_res.reset(new LatticeFermion(SFGrid));
  }

  // src4[r], res4[r] for r < nrhs, all on the full 4D grid.
  void Solve(MobiusFermionD &D, const std::vector<LatticeFermion> &src4, std::vector<LatticeFermion> &res4)
  {
    assert( (int)src4.size() == nrhs && "Solve: need exactly nrhs sources" );
    assert( (int)res4.size() == nrhs && "Solve: need exactly nrhs results" );
    for(int r=0; r<nrhs; r++){
      D.ImportPhysicalFermionSource(src4[r], f5[r]);
    }
    Grid_split(f5, *s_src);
    *s_res = Zero();
    (*schur)(*sD, *s_src, *s_res, ZG);
    Grid_unsplit(sol5, *s_res);
    for(int r=0; r<nrhs; r++){
      D.ExportPhysicalFermionSolution(sol5[r], res4[r]);
    }
  }
};

// Parse "--split sx sy sz st" (default 1^Nd = no split).
static void ParseSplit(int argc, char **argv, const Coordinate &mpi_layout, Coordinate &mpi_split)
{
  mpi_split = Coordinate(mpi_layout.size(), 1);
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
}

// ==========================================================================
// Contraction of B-form baryon components (baryon_variants_impl_plan_claude.md Sec. 2)
//
//   C_{AA'}(t) = sum_x < B_A(x,t) B^\dagger_{A'}(0) >
//              = P(A') sum_{a,b in S_4} \epsilon_a \epsilon_b det_{ij}[ G^{a_i b_j}_{A_i A'_j}(x,0) ]
//
// with G the 16x16 PD spin-colour site propagator (composite index 4*spin + colour).
// The term tables (weights + row/column composite indices of each 4x4 minor) are
// generated by baryon_variants_terms_gen_claude.py; the operator coefficient tables by
// baryon_variants_ops_gen_claude.py. Both are read from text here.
// ==========================================================================

#include <fstream>
#include <map>

// Composite spin-colour index helpers (row/column index of the 16x16 site matrix).
static const int NSC = Ns * Nc;   // 16

// --------------------------------------------------------------------------
// 4x4 determinant of a row-major matrix M[16] of generic element type E
// (Complex on the GPU lane, vComplex on a CPU build). Laplace expansion along
// rows 0-1: sum over column pairs (j<k) of the 2x2 minor of rows 0-1 times the
// complementary 2x2 minor of rows 2-3 with sign (-1)^{1+j+k}. 30 multiplies, no
// branches, no division.
// --------------------------------------------------------------------------
template<class E>
accelerator_inline E Det4(const E *M)
{
  E m01 = M[0]*M[5]  - M[1]*M[4];
  E m02 = M[0]*M[6]  - M[2]*M[4];
  E m03 = M[0]*M[7]  - M[3]*M[4];
  E m12 = M[1]*M[6]  - M[2]*M[5];
  E m13 = M[1]*M[7]  - M[3]*M[5];
  E m23 = M[2]*M[7]  - M[3]*M[6];
  E n01 = M[8]*M[13] - M[9]*M[12];
  E n02 = M[8]*M[14] - M[10]*M[12];
  E n03 = M[8]*M[15] - M[11]*M[12];
  E n12 = M[9]*M[14] - M[10]*M[13];
  E n13 = M[9]*M[15] - M[11]*M[13];
  E n23 = M[10]*M[15] - M[11]*M[14];
  return m01*n23 - m02*n13 + m03*n12 + m12*n03 - m13*n02 + m23*n01;
}

// --------------------------------------------------------------------------
// Term tables of the elementary pairs (A, A'), read from baryon_variants_terms_claude.txt.
// Flattened for the device: idx holds 8 ints per term (4 row then 4 column composite
// indices), w one weight per term; pair p occupies terms [offset[p], offset[p]+nterm[p]).
// --------------------------------------------------------------------------
class PDTermTable {
public:
  std::vector<std::string> A;        // per pair: sink tuple, e.g. "0002"
  std::vector<std::string> Ap;       // per pair: source tuple
  std::vector<int>         parity;   // per pair: P(A')
  std::vector<int>         offset;   // per pair: first term
  std::vector<int>         nterm;    // per pair: number of terms
  std::vector<int>         idx;      // 8 ints per term
  std::vector<double>      w;        // weight per term (includes P(A'))
  std::string              fname;
  deviceVector<int>        d_idx;    // device copies (uploaded once)
  deviceVector<double>     d_w;
  bool                     uploaded = false;

  int NPairs() const { return (int)A.size(); }

  // Parse "s:c" -> 4*s + c.
  static int ParseSC(const std::string &tok)
  {
    size_t colon = tok.find(':');
    assert( colon != std::string::npos && "term token is not s:c" );
    int s = std::stoi(tok.substr(0, colon));
    int c = std::stoi(tok.substr(colon + 1));
    assert( s >= 0 && s < Ns && c >= 0 && c < Nc && "s:c out of range" );
    return Ns * s + c;
  }

  // P(A') = (-1)^{number of lower (2,3) indices}.
  static int ParityOf(const std::string &tuple)
  {
    int nminus = 0;
    for(size_t i=0; i<tuple.size(); i++){
      if( tuple[i] == '2' || tuple[i] == '3' ) nminus++;
    }
    if( nminus % 2 == 0 ) return 1;
    return -1;
  }

  void Read(const std::string &file)
  {
    fname = file;
    std::ifstream in(file);
    assert( in.good() && "cannot open term table file" );
    std::string line;
    int npairs = -1;
    int pairsRead = 0;
    int termsLeft = 0;
    long wsum = 0;
    while( std::getline(in, line) ){
      if( line.empty() ) continue;
      if( line[0] == '#' ) continue;
      std::stringstream ss(line);
      std::string key;
      ss >> key;
      if( key == "npairs" ){
        ss >> npairs;
        continue;
      }
      if( key == "pair" ){
        assert( termsLeft == 0 && "previous pair has missing term lines" );
        if( pairsRead > 0 ){
          assert( wsum == 576 && "sum|w| of previous pair != 576" );
        }
        std::string a, ap;
        int p, n;
        ss >> a >> ap >> p >> n;
        assert( a.size() == 4 && ap.size() == 4 && "tuple must have 4 digits" );
        assert( p == ParityOf(ap) && "P(A') in file inconsistent with A'" );
        A.push_back(a);
        Ap.push_back(ap);
        parity.push_back(p);
        offset.push_back((int)w.size());
        nterm.push_back(n);
        termsLeft = n;
        wsum = 0;
        pairsRead++;
        continue;
      }
      // term line: <w> s:c s:c s:c s:c | s:c s:c s:c s:c
      assert( termsLeft > 0 && "term line outside a pair block" );
      std::stringstream ts(line);
      long wt;
      ts >> wt;
      std::string tok;
      for(int i=0; i<4; i++){
        ts >> tok;
        idx.push_back(ParseSC(tok));
      }
      ts >> tok;
      assert( tok == "|" && "expected | between rows and columns" );
      for(int j=0; j<4; j++){
        ts >> tok;
        idx.push_back(ParseSC(tok));
      }
      w.push_back((double)wt);
      wsum += std::labs(wt);
      termsLeft--;
    }
    assert( termsLeft == 0 && "last pair has missing term lines" );
    if( pairsRead > 0 ) assert( wsum == 576 && "sum|w| of last pair != 576" );
    assert( npairs == pairsRead && "npairs header != number of pair blocks" );
    assert( (int)w.size() * 8 == (int)idx.size() );
    std::cout << GridLogMessage << "PDTermTable: read " << pairsRead << " pairs, "
              << w.size() << " terms from " << file << std::endl;
  }

  int Find(const std::string &a, const std::string &ap) const
  {
    for(int p=0; p<NPairs(); p++){
      if( A[p] == a && Ap[p] == ap ) return p;
    }
    return -1;
  }

  void Upload()
  {
    if( uploaded ) return;
    d_idx.resize(idx.size());
    d_w.resize(w.size());
    acceleratorCopyToDevice((void*)&idx[0], (void*)&d_idx[0], idx.size() * sizeof(int));
    acceleratorCopyToDevice((void*)&w[0],   (void*)&d_w[0],   w.size()   * sizeof(double));
    uploaded = true;
    // round-trip check of the device copies
    std::vector<int>    idx2(idx.size(), -1);
    std::vector<double> w2(w.size(), 0.0);
    acceleratorCopyFromDevice((void*)&d_idx[0], (void*)&idx2[0], idx.size() * sizeof(int));
    acceleratorCopyFromDevice((void*)&d_w[0],   (void*)&w2[0],   w.size()   * sizeof(double));
    long bad = 0;
    for(size_t i=0; i<idx.size(); i++){
      if( idx2[i] != idx[i] ) bad++;
    }
    for(size_t i=0; i<w.size(); i++){
      if( w2[i] != w[i] ) bad++;
    }
    std::cout << GridLogMessage << "PDTermTable: device table round-trip: " << bad << " mismatches of "
              << idx.size() + w.size() << " entries; d_idx at " << (void*)&d_idx[0]
              << " d_w at " << (void*)&d_w[0] << std::endl;
    assert( bad == 0 && "device term table round-trip failed" );
  }

  // Site density dens(x) = sum_terms w * det[ G(x)[rows, cols] ] for pair p:
  // one fused kernel over sites; the 16 entries of each minor are gathered directly
  // from the propagator view (coalescedRead on a sub-element), the determinant is
  // Det4. C_{AA'}(t) = sliceSum(dens, Tdir).
  int kernelMode = 0;   // 0 = gather sub-elements from the view, 1 = read the whole site object

  void ContractPair(const LatticePropagator &G, int p, LatticeComplex &dens)
  {
    assert( uploaded && "call Upload() before ContractPair" );
    GridBase *grid = G.Grid();
    const int    *idx_p = &d_idx[0] + 8 * offset[p];
    const double *w_p   = &d_w[0] + offset[p];
    // NOTE: Grid's accelerator_for macro declares a local 'nt' (thread count) around the
    // lambda; a variable named nt here would be shadowed inside the kernel. Use ntermP.
    const int     ntermP = nterm[p];
    autoView(vG, G,    AcceleratorRead);
    autoView(vd, dens, AcceleratorWrite);
    if( kernelMode == 0 ){
      accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
        typedef decltype(coalescedRead(vd[0])) cVec;
        typedef decltype(coalescedRead(vG[0]()(0,0)(0,0))) Elem;
        // accumulator zeroed as (acc - acc): portable for Complex (GPU lane) and vComplex (CPU)
        Elem acc = coalescedRead(vG[ss]()(0,0)(0,0));
        acc = acc - acc;
        for(int k=0; k<ntermP; k++){
          const int *t = idx_p + 8*k;
          Elem M[16];
          for(int i=0; i<4; i++){
            int sr = t[i] / Nc;
            int cr = t[i] % Nc;
            for(int j=0; j<4; j++){
              int sc = t[4+j] / Nc;
              int cc = t[4+j] % Nc;
              M[4*i+j] = coalescedRead(vG[ss]()(sr,sc)(cr,cc));
            }
          }
          Elem dv = Det4(M);
          acc = acc + Complex(w_p[k], 0.0) * dv;
        }
        cVec res;
        res()()() = acc;
        coalescedWrite(vd[ss], res);
      });
    } else {
      accelerator_for(ss, grid->oSites(), grid->Nsimd(), {
        typedef decltype(coalescedRead(vd[0])) cVec;
        typedef decltype(coalescedRead(vG[0]()(0,0)(0,0))) Elem;
        auto Gsite = vG(ss);
        Elem acc = Gsite()(0,0)(0,0);
        acc = acc - acc;
        for(int k=0; k<ntermP; k++){
          const int *t = idx_p + 8*k;
          Elem M[16];
          for(int i=0; i<4; i++){
            int sr = t[i] / Nc;
            int cr = t[i] % Nc;
            for(int j=0; j<4; j++){
              int sc = t[4+j] / Nc;
              int cc = t[4+j] % Nc;
              M[4*i+j] = Gsite()(sr,sc)(cr,cc);
            }
          }
          Elem dv = Det4(M);
          acc = acc + Complex(w_p[k], 0.0) * dv;
        }
        cVec res;
        res()()() = acc;
        coalescedWrite(vd[ss], res);
      });
    }
  }
};

// --------------------------------------------------------------------------
// Operator coefficient tables, read from baryon_variants_ops_claude.txt (9 highest-weight
// operators) or baryon_variants_ops_all_claude.txt (all 35 states).
// --------------------------------------------------------------------------
struct PDOperator {
  std::string name;               // e.g. "1m_31"
  std::string JP;                 // e.g. "1-"
  int nplus;
  int nminus;
  int J;
  int M;
  std::vector<std::string> comp;  // A of each component
  std::vector<double>      coeff; // c_A, same order
};

class PDOperatorTable {
public:
  std::vector<PDOperator>               ops;
  std::vector<std::string>              blockJP;   // per block
  std::vector<int>                      blockM;    // per block
  std::vector<std::vector<std::string> > blockOps; // per block: operator names
  std::string                           fname;

  void Read(const std::string &file)
  {
    fname = file;
    std::ifstream in(file);
    assert( in.good() && "cannot open operator table file" );
    std::string line;
    int nops = -1;
    int nblocks = -1;
    int compLeft = 0;
    while( std::getline(in, line) ){
      if( line.empty() ) continue;
      if( line[0] == '#' ) continue;
      std::stringstream ss(line);
      std::string key;
      ss >> key;
      if( key == "nops" ){
        ss >> nops;
        continue;
      }
      if( key == "nblocks" ){
        ss >> nblocks;
        continue;
      }
      if( key == "op" ){
        assert( compLeft == 0 && "previous operator has missing component lines" );
        PDOperator op;
        int ncomp;
        ss >> op.name >> op.JP >> op.nplus >> op.nminus >> op.J >> op.M >> ncomp;
        ops.push_back(op);
        compLeft = ncomp;
        continue;
      }
      if( key == "block" ){
        std::string jp;
        int M, nop;
        ss >> jp >> M >> nop;
        std::vector<std::string> names(nop);
        for(int i=0; i<nop; i++) ss >> names[i];
        blockJP.push_back(jp);
        blockM.push_back(M);
        blockOps.push_back(names);
        continue;
      }
      // component line: <A> <c_A>
      assert( compLeft > 0 && "component line outside an operator block" );
      double c;
      ss >> c;
      assert( key.size() == 4 && "component tuple must have 4 digits" );
      ops.back().comp.push_back(key);
      ops.back().coeff.push_back(c);
      compLeft--;
    }
    assert( compLeft == 0 && "last operator has missing component lines" );
    assert( nops == (int)ops.size() && "nops header != number of operators" );
    assert( nblocks == (int)blockJP.size() && "nblocks header != number of blocks" );
    std::cout << GridLogMessage << "PDOperatorTable: read " << ops.size() << " operators, "
              << blockJP.size() << " blocks from " << file << std::endl;
  }

  int Find(const std::string &name) const
  {
    for(int i=0; i<(int)ops.size(); i++){
      if( ops[i].name == name ) return i;
    }
    return -1;
  }

  // Index of the block containing operator name (-1 if none).
  int BlockOf(const std::string &name) const
  {
    for(int b=0; b<(int)blockJP.size(); b++){
      for(int i=0; i<(int)blockOps[b].size(); i++){
        if( blockOps[b][i] == name ) return b;
      }
    }
    return -1;
  }
};

// --------------------------------------------------------------------------
// Naive reference (host, one site): the Wick form with the 24 quark pairings and the
// full \epsilon_a \epsilon_b colour sum, no determinant:
//   dens(x) = P(A') sum_{pi} sum_{a,b} \epsilon_a \epsilon_b prod_i G[(A_i,a_i),(A'_{pi(i)},b_i)]
// Used by --naive-check to validate the term tables + kernel on a real lattice.
// --------------------------------------------------------------------------
class PDNaive {
public:
  int perms[24][4];
  int sign[24];

  PDNaive()
  {
    int p[4] = {0, 1, 2, 3};
    int n = 0;
    do {
      for(int i=0; i<4; i++) perms[n][i] = p[i];
      sign[n] = PermSign(p);
      n++;
    } while( std::next_permutation(p, p + 4) );
    assert( n == 24 );
  }

  static int PermSign(const int *p)
  {
    int inv = 0;
    for(int i=0; i<4; i++){
      for(int j=i+1; j<4; j++){
        if( p[i] > p[j] ) inv++;
      }
    }
    if( inv % 2 == 0 ) return 1;
    return -1;
  }

  static void ParseTuple(const std::string &s, int *A)
  {
    for(int i=0; i<4; i++) A[i] = s[i] - '0';
  }

  Complex SiteDensity(const SpinColourMatrix &Gs, const std::string &As, const std::string &Aps) const
  {
    int A[4];
    int Ap[4];
    ParseTuple(As, A);
    ParseTuple(Aps, Ap);
    Complex val(0.0, 0.0);
    for(int ip=0; ip<24; ip++){
      for(int ia=0; ia<24; ia++){
        for(int ib=0; ib<24; ib++){
          Complex prod(1.0, 0.0);
          for(int i=0; i<4; i++){
            int sr = A[i];
            int cr = perms[ia][i];
            int sc = Ap[perms[ip][i]];
            int cc = perms[ib][i];
            prod = prod * Gs()(sr,sc)(cr,cc);
          }
          val = val + (RealD)(sign[ia] * sign[ib]) * prod;
        }
      }
    }
    return (RealD)PDTermTable::ParityOf(Aps) * val;
  }
};

// --------------------------------------------------------------------------
// Reader for the full PD dump (point2all_prop_dumper_full_claude.cc): records are read
// sequentially; the 16 blocks with srcIndex == corner are promoted to double and poked
// into G as spin block (s, s'). Records of other corners are skipped. Returns the
// PropRecord of block (0,0) of that corner in rec, and the number of corners seen.
// --------------------------------------------------------------------------
static uint64_t CountLimeBinaryRecords(const std::string &infile)
{
  const std::string tag = "ildg-binary-data";
  const size_t chunk = 16 * 1024 * 1024;
  std::ifstream f(infile, std::ios::binary);
  assert( f.good() && "cannot open propagator dump" );
  std::string buf;
  std::string carry;
  buf.resize(chunk);
  uint64_t n = 0;
  while( f ){
    f.read(&buf[0], chunk);
    std::streamsize got = f.gcount();
    if( got <= 0 ) break;
    std::string window = carry + buf.substr(0, (size_t)got);
    size_t pos = 0;
    while( (pos = window.find(tag, pos)) != std::string::npos ){
      n++;
      pos += tag.size();
    }
    // keep the last tag.size()-1 bytes so a tag split across chunks is still found
    if( window.size() >= tag.size() ){
      carry = window.substr(window.size() - (tag.size() - 1));
    } else {
      carry = window;
    }
  }
  return n;
}

static int ReadPDFullCorner(const std::string &infile, GridCartesian *UGrid, GridCartesian *UGridF,
                            int corner, LatticePropagator &G, PropRecord &rec)
{
  G = Zero();
  std::vector<bool> got(NPDSpin * NPDSpin, false);
  std::map<int,int> cornersSeen;

  // Count binary records first (ScidacReader has no record iterator we can query):
  // chunked scan for the LIME type string on the boss rank, then broadcast.
  uint64_t nrec = 0;
  if( UGrid->IsBoss() ){
    nrec = CountLimeBinaryRecords(infile);
  }
  UGrid->GlobalSum(nrec);
  assert( nrec > 0 && (nrec % (NPDSpin * NPDSpin)) == 0 && "record count is not a multiple of 16" );

  ScidacReader RD;
  RD.open(infile);
  LatticeColourMatrixF qF(UGridF);
  LatticeColourMatrix  q(UGrid);
  for(uint64_t r=0; r<nrec; r++){
    PropRecord rr;
    RD.readScidacFieldRecord(qF, rr);
    cornersSeen[rr.srcIndex] = 1;
    if( rr.srcIndex != corner ) continue;
    assert( rr.store.size() == 3 && rr.store[0] == 'q' && "record store tag is not q<s><s'>" );
    int s  = rr.store[1] - '0';
    int sp = rr.store[2] - '0';
    assert( s >= 0 && s < NPDSpin && sp >= 0 && sp < NPDSpin && "store tag out of range" );
    assert( !got[NPDSpin*s + sp] && "duplicate block record" );
    precisionChange(q, qF);
    pokeSpin(G, q, s, sp);
    got[NPDSpin*s + sp] = true;
    if( s == 0 && sp == 0 ) rec = rr;
  }
  RD.close();
  for(int k=0; k<NPDSpin*NPDSpin; k++){
    assert( got[k] && "missing a PD block record for the requested corner" );
  }
  return (int)cornersSeen.size();
}

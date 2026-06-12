/*
 * two_baryon_corr_claude.cc
 *
 * Two-baryon correlator for SU(Nc=4) one-flavour SDM theory.
 *
 * Geometry (see two_baryon_impl_plan_claude.md):
 *   C_2B(t) = < B(O',t) B(M',t) Bbar(O,0) Bbar(M,0) >
 * with O = (0,0,0) and M = (L/2,L/2,L/2) the two most distant spatial points,
 * source timeslice fixed at t_src = 0, sink timeslice t scanned over all slices.
 * Both baryons are the spin-2, all-spin-up, positive-parity "0000" operator
 *   B(x) = \varepsilon_{abcd} q_0^a q_0^b q_0^c q_0^d
 * built from the upper Pauli-Dirac (nonrelativistic) spinor component 0.
 *
 * Chunk 1: sources, two point-to-all solves, Weyl -> Pauli-Dirac rotation,
 *          extraction of the spin-0 colour matrices q00_O, q00_M.
 */

#include <Grid/Grid.h>
#include <algorithm>
#include <vector>
#include <random>

using namespace std;
using namespace Grid;

typedef std::vector<std::vector<ComplexD>> CMat; // small dense complex matrix

// Complex determinant via Gaussian elimination with partial pivoting, O(n^3).
// Host-only std::complex arithmetic (no Eigen / thrust). Takes A by value.
static ComplexD detLU(CMat A)
{
  int n = A.size();
  ComplexD det(1.0, 0.0);
  for(int k=0;k<n;k++){
    int piv = k; double best = abs(A[k][k]);
    for(int i=k+1;i<n;i++){ double v = abs(A[i][k]); if(v>best){ best=v; piv=i; } }
    if(best == 0.0) return ComplexD(0.0,0.0);
    if(piv != k){ std::swap(A[piv], A[k]); det = det * ComplexD(-1.0,0.0); }
    det = det * A[k][k];
    ComplexD inv = ComplexD(1.0,0.0) / A[k][k];
    for(int i=k+1;i<n;i++){
      ComplexD f = A[i][k] * inv;
      for(int j=k+1;j<n;j++) A[i][j] = A[i][j] - f * A[k][j];
    }
  }
  return det;
}

// HDF5-serialisable container for a correlator C(t) (same as the baryon code).
class MesonFile: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFile, std::vector<Complex>, data);
};

// Scalar complex determinant of an n x n matrix via the Leibniz formula
// (sum over permutations with fermion sign). n is small (4 for a single SU(4)
// baryon, 8 for two baryons), so brute-force permutation enumeration is fine.
// This is the kernel reused by the two-baryon contraction in Chunk 3.
static ComplexD detLeibniz(const std::vector<std::vector<ComplexD>> &A)
{
  int n = A.size();
  std::vector<int> perm(n);
  for(int i=0;i<n;i++) perm[i]=i;
  ComplexD det = 0.0;
  do {
    int inv=0;
    for(int i=0;i<n;i++)
      for(int j=i+1;j<n;j++)
        if(perm[i]>perm[j]) inv++;
    double sign = (inv%2==0) ? 1.0 : -1.0;
    ComplexD term = sign;
    for(int i=0;i<n;i++) term *= A[i][perm[i]];
    det += term;
  } while(std::next_permutation(perm.begin(), perm.end()));
  return det;
}

// --------------------------------------------------------------------------
// Two-baryon "0000" contraction.
//
// Each baryon is \varepsilon q_0^4 (SU(4) color singlet, PD spin-0). The full
// 8-quark, fully-antisymmetrised correlator
//   C_2B = \sum_{4 \varepsilon's} \varepsilon\varepsilon\varepsilon\varepsilon
//          det(Q_8x8[colors])
// collapses (det is alternating in rows/cols, each \varepsilon sum -> a factor
// 4! = 24 at fixed colors 0..3) to a single 8x8 determinant:
//   C_2B = 24^4 * det( [[P_OO, P_OM],[P_MO, P_MM]] ),
// where P_{X'Y}[i][j] is the PD spin-0 propagator colour matrix from source
// point Y to sink point X' (i = sink colour, j = source colour):
//   P_OO = q00_O(O'),  P_OM = q00_M(O'),  P_MO = q00_O(M'),  P_MM = q00_M(M').
// --------------------------------------------------------------------------

static const RealD POW24_4 = 24.0*24.0*24.0*24.0; // 331776 = (4!)^4

// Assemble the 8x8 block matrix Q from the four 4x4 colour blocks.
static void assemble8x8(const CMat &OO, const CMat &OM,
                        const CMat &MO, const CMat &MM, CMat &Q8)
{
  Q8.assign(8, std::vector<ComplexD>(8));
  for(int i=0;i<Nc;i++){
    for(int j=0;j<Nc;j++){
      Q8[i]   [j]    = OO[i][j];
      Q8[i]   [j+Nc] = OM[i][j];
      Q8[i+Nc][j]    = MO[i][j];
      Q8[i+Nc][j+Nc] = MM[i][j];
    }
  }
}

// Two-baryon correlator from the four colour blocks via the 8x8 determinant.
static ComplexD TwoBaryonCorr(const CMat &OO, const CMat &OM,
                              const CMat &MO, const CMat &MM)
{
  CMat Q8;
  assemble8x8(OO, OM, MO, MM, Q8);
  return POW24_4 * detLeibniz(Q8);
}

// Extract the 4x4 colour block of a LatticeColourMatrix at global site X.
static CMat colourBlockAt(const LatticeColourMatrix &q, const Coordinate &X)
{
  ColourMatrix cm;
  peekSite(cm, q, X);
  CMat P(Nc, std::vector<ComplexD>(Nc));
  for(int i=0;i<Nc;i++)
    for(int j=0;j<Nc;j++)
      P[i][j] = static_cast<ComplexD>(cm()()(i,j));
  return P;
}

// Brute-force two-baryon correlator: explicit sum over the four \varepsilon
// colour assignments (24^4 terms) with an Eigen 8x8 determinant per term.
// Used only as a validation reference for TwoBaryonCorr.
static ComplexD TwoBaryonBruteForce(const CMat &OO, const CMat &OM,
                                    const CMat &MO, const CMat &MM)
{
  // 24 permutations of {0,1,2,3} with their signs.
  std::vector<std::array<int,4>> perms;
  std::vector<double> signs;
  std::array<int,4> p = {0,1,2,3};
  do {
    int inv=0;
    for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(p[i]>p[j]) inv++;
    perms.push_back(p);
    signs.push_back((inv%2==0)?1.0:-1.0);
  } while(std::next_permutation(p.begin(), p.end()));

  // Pointer table: which 4x4 block for (sink half, source half).
  // sink half: 0 = O', 1 = M';  source half: 0 = O, 1 = M.
  const CMat* blk[2][2] = { { &OO, &OM }, { &MO, &MM } };

  CMat M(8, std::vector<ComplexD>(8));
  ComplexD C = 0.0;
  for(size_t ia=0; ia<perms.size(); ia++){           // sink O' colours (abcd)
  for(size_t ie=0; ie<perms.size(); ie++){           // sink M' colours (efgh)
  for(size_t ip=0; ip<perms.size(); ip++){           // source O colours (a'b'c'd')
  for(size_t iq=0; iq<perms.size(); iq++){           // source M colours (e'f'g'h')
    double sgn = signs[ia]*signs[ie]*signs[ip]*signs[iq];

    for(int i=0;i<8;i++){
      int sinkHalf  = (i<Nc)?0:1;
      int cs        = (i<Nc) ? perms[ia][i] : perms[ie][i-Nc];     // sink colour
      for(int j=0;j<8;j++){
        int srcHalf = (j<Nc)?0:1;
        int csrc    = (j<Nc) ? perms[ip][j] : perms[iq][j-Nc];     // source colour
        M[i][j] = (*blk[sinkHalf][srcHalf])[cs][csrc];
      }
    }
    C += sgn * detLU(M);
  }}}}
  return C;
}

// Self-test of the two-baryon contraction on random colour blocks:
//   (1) TwoBaryonCorr (24^4 det) vs TwoBaryonBruteForce (explicit \varepsilon sum)
//   (2) factorisation limit P_OM=P_MO=0 -> product of single-baryon correlators
//       (single baryon = 576 det P = 24 * weight-1 density).
static void SelfTestTwoBaryon()
{
  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  auto randMat = [&](){
    CMat P(Nc, std::vector<ComplexD>(Nc));
    for(int i=0;i<Nc;i++) for(int j=0;j<Nc;j++) P[i][j] = ComplexD(u(rng), u(rng));
    return P;
  };
  CMat OO = randMat(), OM = randMat(), MO = randMat(), MM = randMat();

  ComplexD c_formula = TwoBaryonCorr(OO, OM, MO, MM);
  ComplexD c_brute   = TwoBaryonBruteForce(OO, OM, MO, MM);
  RealD diff1 = abs(c_formula - c_brute);

  CMat Z(Nc, std::vector<ComplexD>(Nc, ComplexD(0.0,0.0)));
  ComplexD c_decoupled = TwoBaryonCorr(OO, Z, Z, MM);
  ComplexD c_O = 576.0 * detLeibniz(OO);   // single-baryon correlator (eq 1.2.17)
  ComplexD c_M = 576.0 * detLeibniz(MM);
  RealD diff2 = abs(c_decoupled - c_O * c_M);

  std::cout << GridLogMessage << "=== two-baryon contraction self-test ===" << std::endl;
  std::cout << GridLogMessage << "  formula  = " << c_formula << std::endl;
  std::cout << GridLogMessage << "  brute    = " << c_brute   << std::endl;
  std::cout << GridLogMessage << "  |formula - brute|       = " << diff1 << std::endl;
  std::cout << GridLogMessage << "  |decoupled - C_O*C_M|   = " << diff2 << std::endl;
  RealD scale = std::max(abs(c_brute), abs(c_O*c_M));
  std::cout << GridLogMessage << "  max relative diff       = "
            << std::max(diff1, diff2)/scale << std::endl;
}

// Single-baryon "0000" density field (full lattice), weight +1, identical to
// BaryonTrace_hdf bar_0000 in baryons_0000_dirac.cc:
//   cf(x) = \varepsilon_{abcd} \varepsilon_{a'b'c'd'}
//           q00[a,a'] q00[b,b'] q00[c,c'] q00[d,d']
// where q00 is the Pauli-Dirac spin-0 colour matrix of the propagator from the
// source point to x (row = sink colour, col = source colour).
static LatticeComplex BaryonSingleDensity(const LatticeColourMatrix &q00)
{
  LatticeComplex cf(q00.Grid());
  cf = Zero();

  // sink colours {a,b,c,d}, all distinct -> \varepsilon_{abcd}
  for(int a=0;a<Nc;a++){
  for(int b=0;b<Nc;b++){ if(b==a) continue;
  for(int c=0;c<Nc;c++){ if(c==a||c==b) continue;
  for(int d=0;d<Nc;d++){ if(d==a||d==b||d==c) continue;
    RealD parity = ((a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d) < 0) ? -1.0 : 1.0;
    // source colours {ap,bp,cp,dp}, all distinct -> \varepsilon_{a'b'c'd'}
    for(int ap=0;ap<Nc;ap++){
    for(int bp=0;bp<Nc;bp++){ if(bp==ap) continue;
    for(int cp=0;cp<Nc;cp++){ if(cp==ap||cp==bp) continue;
    for(int dp=0;dp<Nc;dp++){ if(dp==ap||dp==bp||dp==cp) continue;
      RealD parity_p = ((ap-bp)*(ap-cp)*(ap-dp)*(bp-cp)*(bp-dp)*(cp-dp) < 0) ? -1.0 : 1.0;
      cf += parity*parity_p
            * peekColour(q00,a,ap) * peekColour(q00,b,bp)
            * peekColour(q00,c,cp) * peekColour(q00,d,dp);
    }}}}
  }}}}
  return cf;
}

// Point source: Kronecker SpinColourMatrix at a single site \coor.
void PointSource(Coordinate &coor, LatticePropagator &source)
{
  source = Zero();
  SpinColourMatrix kronecker; kronecker = 1.0;
  pokeSite(kronecker, source, coor);
}

// Solve D q = source for all 12 spin-colour source components, building the
// point-to-all propagator. Same structure as baryons_0000_dirac.cc.
template<class Action>
void Solve(Action &D, LatticePropagator &source, LatticePropagator &propagator)
{
  GridBase *UGrid = D.GaugeGrid();
  GridBase *FGrid = D.FermionGrid();

  LatticeFermion src4   (UGrid);
  LatticeFermion src5   (FGrid);
  LatticeFermion result5(FGrid);
  LatticeFermion result4(UGrid);

  ConjugateGradient<LatticeFermion> CG(1.0e-9, 100000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion> schur(CG);
  ZeroGuesser<LatticeFermion> ZG;

  for(int s=0; s<Nd; s++){
    for(int c=0; c<Nc; c++){
      PropToFerm<Action>(src4, source, s, c);

      D.ImportPhysicalFermionSource(src4, src5);

      result5 = Zero();
      schur(D, src5, result5, ZG);
      std::cout << GridLogMessage
                << "spin " << s << " color " << c
                << " norm2(src5d) "    << norm2(src5)
                << " norm2(result5d) " << norm2(result5) << std::endl;

      D.ExportPhysicalFermionSolution(result5, result4);

      FermToProp<Action>(propagator, result4, s, c);
    }
  }
}

// Unitary matrix U that rotates a quark propagator from the Weyl (chiral) basis,
// in which Grid works, to the Pauli-Dirac (PD) basis used for spectroscopy with
// definite parity (thermo_gw_details eq. 1.2.3-1.2.6):
//   gamma(W) = U gamma(PD) U^\dagger,  S(W) = U S(PD) U^\dagger
//   => S(PD) = U^\dagger S(W) U,  with the 1/sqrt(2) normalisation absorbed
//      into the 0.5 prefactor below (q_DP = 0.5 * U * q * Udagger).
// In the PD basis the upper two spinor components (0,1) are positive parity,
// the lower (2,3) negative parity; component 0 = positive parity, spin up.
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

  // Validate the two-baryon contraction kernel before touching the lattice.
  SelfTestTwoBaryon();

  // ------------------------------------------------------------------
  // Two point sources on the source timeslice t_src = 0:
  //   O = (0,0,0,0),  M = (L/2,L/2,L/2,0)  (most distant spatial points)
  // ------------------------------------------------------------------
  Coordinate latt = GridDefaultLatt();
  Coordinate Opt({0, 0, 0, 0});
  Coordinate Mpt({latt[0]/2, latt[1]/2, latt[2]/2, 0});

  std::cout << GridLogMessage << "Source O = "
            << Opt[0] << "," << Opt[1] << "," << Opt[2] << "," << Opt[3] << std::endl;
  std::cout << GridLogMessage << "Source M = "
            << Mpt[0] << "," << Mpt[1] << "," << Mpt[2] << "," << Mpt[3] << std::endl;

  LatticePropagator srcO(UGrid), srcM(UGrid);
  PointSource(Opt, srcO);
  PointSource(Mpt, srcM);

  // Point-to-all propagators from O and from M.
  LatticePropagator PropO(UGrid), PropM(UGrid);
  std::cout << GridLogMessage << "Solving for source O ..." << std::endl;
  Solve(D, srcO, PropO);
  std::cout << GridLogMessage << "Solving for source M ..." << std::endl;
  Solve(D, srcM, PropM);

  // ------------------------------------------------------------------
  // Rotate Weyl -> Pauli-Dirac and extract the spin-0 colour matrix.
  // q00_Y(x) is a LatticeColourMatrix: row = sink colour, col = source colour,
  // for the propagator from source point Y with both spinor indices = 0.
  // ------------------------------------------------------------------
  SpinMatrix U       = WeylToPauliDiracU();
  SpinMatrix Udagger = transpose(U);

  LatticePropagator PropO_DP = 0.5 * U * PropO * Udagger;
  LatticePropagator PropM_DP = 0.5 * U * PropM * Udagger;

  LatticeColourMatrix q00_O = peekSpin(PropO_DP, 0, 0);
  LatticeColourMatrix q00_M = peekSpin(PropM_DP, 0, 0);

  std::cout << GridLogMessage << "norm2(q00_O) = " << norm2(q00_O) << std::endl;
  std::cout << GridLogMessage << "norm2(q00_M) = " << norm2(q00_M) << std::endl;

  // ------------------------------------------------------------------
  // Chunk 2: single-baryon validation.
  //   (a) explicit \varepsilon\varepsilon density (matches baryons_0000_dirac),
  //       sliceSum over T -> C^{0000}(t), written to HDF5 as "bar_0000_t".
  //   (b) cross-check the detLeibniz kernel via the per-site identity
  //       density(x) = 24 * det( q00(x) )  [4x4 colour determinant].
  // ------------------------------------------------------------------
  std::unique_ptr<Hdf5Writer> WR;
  if (UGrid->IsBoss()){
    WR = std::make_unique<Hdf5Writer>(outfile);
  }

  LatticeComplex bar0000_O = BaryonSingleDensity(q00_O);

  std::vector<TComplex> bar_T;
  sliceSum(bar0000_O, bar_T, Tdir);
  int nt = bar_T.size();
  std::vector<Complex> corr(nt);
  for(int t=0;t<nt;t++){
    corr[t] = TensorRemove(bar_T[t]);
    std::cout << GridLogMessage << "bar_0000 t " << t << " " << corr[t] << std::endl;
  }
  {
    MesonFile MF; MF.data = corr;
    if (WR) write(*WR, std::string("bar_0000_t"), MF);
  }

  // (b) detLeibniz cross-check at a few sink sites: 24*det(q00(X)) == density(X).
  std::vector<Coordinate> checkX = {
    Coordinate({1,0,0,0}),
    Coordinate({0,1,0,0}),
    Coordinate({1,2,3,4}),
    Coordinate({latt[0]/2, 0, 0, latt[3]/2}),
    Coordinate({latt[0]-1, latt[1]-1, latt[2]-1, latt[3]-1}),
  };
  RealD maxdiff = 0.0;
  std::cout << GridLogMessage << "=== detLeibniz vs explicit single-baryon density ===" << std::endl;
  for(auto &X : checkX){
    // explicit per-site density
    typename LatticeComplex::scalar_object dens_s;
    peekSite(dens_s, bar0000_O, X);
    ComplexD dens = TensorRemove(dens_s);

    // 24 * det( q00(X) )  (full 4x4 colour matrix, rows=sink, cols=source)
    ColourMatrix cm;
    peekSite(cm, q00_O, X);
    std::vector<std::vector<ComplexD>> Q(Nc, std::vector<ComplexD>(Nc));
    for(int i=0;i<Nc;i++)
      for(int j=0;j<Nc;j++)
        Q[i][j] = static_cast<ComplexD>(cm()()(i,j));
    ComplexD dens_det = 24.0 * detLeibniz(Q);

    RealD d = abs(dens - dens_det);
    maxdiff = std::max(maxdiff, d);
    std::cout << GridLogMessage
              << "X=(" << X[0] << "," << X[1] << "," << X[2] << "," << X[3] << ")"
              << " explicit=" << dens
              << " 24*det="   << dens_det
              << " |diff|="   << d << std::endl;
  }
  std::cout << GridLogMessage << "single-baryon detLeibniz check: max|diff| = "
            << maxdiff << std::endl;

  // ------------------------------------------------------------------
  // Chunk 4: two-baryon correlator C_2B(t) = 24^4 det Q_8x8 (physical, fully
  // antisymmetrised normalisation) at each sink timeslice, source slice fixed
  // at t_src = 0. Sink points O' = (0,0,0,t), M' = (L/2,L/2,L/2,t). Written to
  // HDF5 as "two_baryon_0000_t".
  // ------------------------------------------------------------------
  std::cout << GridLogMessage << "=== two-baryon correlator C_2B(t) ===" << std::endl;
  std::vector<Complex> corr2b(latt[3]);
  for(int t=0; t<latt[3]; t++){
    Coordinate Op({0, 0, 0, t});
    Coordinate Mp({latt[0]/2, latt[1]/2, latt[2]/2, t});

    CMat P_OO = colourBlockAt(q00_O, Op); // source O -> sink O'
    CMat P_OM = colourBlockAt(q00_M, Op); // source M -> sink O'
    CMat P_MO = colourBlockAt(q00_O, Mp); // source O -> sink M'
    CMat P_MM = colourBlockAt(q00_M, Mp); // source M -> sink M'

    corr2b[t] = TwoBaryonCorr(P_OO, P_OM, P_MO, P_MM);
    std::cout << GridLogMessage << "C_2B t " << t << " " << corr2b[t] << std::endl;
  }
  {
    MesonFile MF; MF.data = corr2b;
    if (WR) write(*WR, std::string("two_baryon_0000_t"), MF);
  }

  Grid_finalize();
}

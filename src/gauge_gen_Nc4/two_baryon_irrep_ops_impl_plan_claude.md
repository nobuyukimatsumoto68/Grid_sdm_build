# two_baryon_irrep_ops -- implementation plan (irrep-projected momentum operators)

Status: Chunk 0 DONE (2026-09-22). Chunks 1-3 in progress in this order; chunk log at the end.

Owner: this session ("scattering projection", grid-claude-6e). Contraction driver
`two_baryon_gevp_contract_claude.cc` and the q00 data are owned by grid-claude-18/ef; they
integrate the C++ header from Chunk 2 into their driver. On approval this plan is copied to
`Grid_sdm_build/src/gauge_gen_Nc4/two_baryon_irrep_ops_impl_plan_claude.md` (Chunk 0) and
kept up to date there.

## Context

The current two-baryon operators (corner pairs, fixed sink points, plus $P=0$ sink smearing) are
a redundant basis: the 3x3 GEVP and a Prony analysis both collapse onto the lowest level
$\approx 2M_B$. The first excited two-baryon level is back-to-back $p=2\pi/L$; with
$M_B\approx3.0$ and $L=24$ its free energy is $2\sqrt{M_B^2+p^2}-2M_B\approx0.023$, i.e. at
threshold, so no rest-frame position-space operator resolves it. Distinct levels need sink
operators with definite total momentum $P$, definite relative-momentum shell, and definite
little-group irrep; moving frames ($P\ne0$) add kinematic points between the rest-frame levels
(Rummukainen-Gottlieb). Output: irrep-projected correlator matrices for the unmatched GEVP,
then finite-volume levels for a Luescher analysis (downstream).

## Physics facts (verified, not assumed)

1. **Bosonic baryon.** $N_c=4$, four quarks, integer spin. Lattice symmetry group is the
   single-valued $O_h$ (48 elements, 10 irreps $A_{1,2},E,T_{1,2}\times g/u$), not the double
   cover. Verified numerically: 48 signed permutation matrices, class structure and characters
   reproduce the known content of momentum shells (table below).
2. **The `0000` operator** $B=\varepsilon_{abcd}q_0^aq_0^bq_0^cq_0^d$ is spin 2, $s_z=+2$,
   positive parity (thermo notes eq. 1.2.8-1.2.13: PD component 0 = positive parity, spin up).
   Spin 2 subduces to $E_g\oplus T_{2g}$ in $O_h$; the $m=+2$ component alone is
   $\tfrac1{\sqrt2}(|E\rangle+|T_2\rangle)$, not an irrep. $E$ vs $T_2$ separation needs the
   $m=-2$ component (all quarks PD spin-down), which needs the sink-spin-1 propagator block.
3. **Stored data limits the spin projection.** The dumper keeps only $q_{00}$: PD sink spin 0,
   source spin 0 (`point2all_prop_dumper_claude.cc:371-375`, sink spins 1 and 3 are discarded).
   So every baryon in every correlator has $s_z=+2$ at both ends, and the two-baryon spin state is
   $|2,2\rangle|2,2\rangle=|S{=}4,S_z{=}4\rangle$ exactly.
4. **Consequence: orbital-only projection is exact for the S-wave, S=4 channel.** With the spin
   state fixed to $|4,4\rangle$ (spectator), projecting the momentum labels only,
   $$\tilde O_\Lambda(p_1,p_2)=\frac{d_\Lambda}{|G|}\sum_{g\in G}\chi_\Lambda(g)^*\,B_{+2}(gp_1)B_{+2}(gp_2),$$
   gives an orbital wavefunction in irrep $\Lambda$ times $|4,4\rangle$. Bose symmetry (symmetric
   spin) forces even $L$; the rest-frame orbital $A_{1g}$ contains $L=0,4,6,\dots$ and no $L=2$, so
   the created state is $|L{=}0,S{=}4;J{=}4\rangle$ up to G-wave, spread over the lattice irreps
   $J{=}4\to A_1+E+T_1+T_2$ which are degenerate up to $O(a^2)$ cubic breaking. Only the
   "spin-only" numerical splittings are lost; this is the same $L$ content a full projection has.
   In moving frames the orbital $A_1$ of the little group contains $L=2$ ($m_L{=}0$), the standard
   S-D mixing, again identical to the full projection. Full spin projection (Tier B, optional)
   needs $q_{10}$ stored in addition: zero extra solves, 2x storage, but a re-dump.
5. **Little groups** of $P$ (units $2\pi/L$): $P=(0,0,n)$: $C_{4v}$ (order 8); $(n,n,0)$:
   $C_{2v}$ (4); $(n,n,n)$: $C_{3v}$ (6); $P=0$: $O_h$. The module computes the stabiliser of
   $P$ in $O_h$ and asserts these orders (self-test), rather than hardcoding membership.
6. **Rest-frame orbital content of shells of pairs $\{p,-p\}$** (computed this session):

   | $p^2$ | #pairs | irreps |
   |---|---|---|
   | 1 | 3 | $A_{1g}+E_g$ |
   | 2 | 6 | $A_{1g}+E_g+T_{2g}$ |
   | 3 | 4 | $A_{1g}+T_{2g}$ |
   | 4 | 3 | $A_{1g}+E_g$ |
   | 5 | 12 | $A_{1g}+A_{2g}+2E_g+T_{1g}+T_{2g}$ |

   Only $g$ irreps appear (even $L$), consistent with point 4.
7. **Momentum-projected sink contraction factorises.** Rows 0-3 of the $8\times8$ block matrix
   $Q(x_1,x_2)$ depend on sink point $x_1$ only, rows 4-7 on $x_2$ only. Laplace expansion along
   rows 0-3:
   $$\det Q(x_1,x_2)=\sum_{S\subset\{0..7\},|S|=4}\epsilon_S\,D_S(x_1)\,D_{\bar S}(x_2),\qquad
   D_S(x)=\det\big([\,q^a_{00}(x)\,|\,q^b_{00}(x)\,]_{\,\cdot,S}\big),\quad
   \epsilon_S=(-1)^{\sum_{i\in S}i+6},$$
   70 terms; $D_S(x)$ is a "baryon block" (4 sink quarks at $x$ contracted with 4 of the 8 source
   quark slots). Verified numerically (rel. err. 3e-15). Hence
   $$C_{(a,b)}(p_1,p_2;t)=\sum_{x_1,x_2}e^{-ip_1x_1-ip_2x_2}\,24^4\det Q
   =24^4\sum_S\epsilon_S\,\tilde D_S(p_1,t)\,\tilde D_{\bar S}(p_2,t),$$
   with $\tilde D_S$ the spatial FFT of $D_S$: cost $V\log V$ per block instead of $V^2$.
   $S=\{0,1,2,3\}$ is the single-baryon density: $D_{0123}(x)=\det q^a_{00}(x)=$
   `BaryonSingleDensity`$/24$, giving $C_B(p,t)$ for free. The $x_1=x_2$ term vanishes
   automatically (identical rows).
8. **Source side.** Sources stay at the 8 corners (unprojected pairs = unmatched GEVP; sink
   projection alone isolates the irrep since the projector commutes with the transfer matrix).
   The corner set $\{0,L/2\}^3$ is $O_h$-invariant about the origin, so source pairs can also be
   orbit-averaged. When combining source pairs related by a translation $s\in\{0,L/2\}^3$ the
   phase $e^{iP\cdot s}=\pm1$ must be applied (matters only for $P$ with odd components).

## References (to cite in the plan header and in code comments)

- Detmold, Savage, "A method to study complex systems of mesons in lattice QCD" / multi-baryon
  determinant contraction: W. Detmold, K. Orginos, PRD 87 (2013) 114512, arXiv:1207.1452;
  W. Detmold, M. J. Savage, PRD 82 (2010) 014511, arXiv:1001.2768.
- T. Doi, M. G. Endres, "Unified contraction algorithm for multi-baryon correlators",
  CPC 184 (2013) 117, arXiv:1205.0585.
- C. Morningstar et al., "Extended hadron and two-hadron operators of definite momentum for
  spectrum calculations in lattice QCD", PRD 88 (2013) 014511, arXiv:1303.6816.
- T. Luu, M. J. Savage, "Extracting scattering phase shifts in higher partial waves from lattice
  QCD calculations", PRD 83 (2011) 114508, arXiv:1101.3347 (cubic-group projection tables of
  two-particle momentum operators, rest and moving frames).
- M. Goeckeler et al., "Scattering phases for meson and baryon resonances on general
  moving-frame lattices", PRD 86 (2012) 094513, arXiv:1206.4141 (little groups, tables).
- K. Rummukainen, S. Gottlieb, NPB 450 (1995) 397, hep-lat/9503028 (moving-frame Luescher).
- S. Basak et al., PRD 72 (2005) 094506, hep-lat/0506029 and PRD 72 (2005) 074501,
  hep-lat/0508018 (group-theoretical baryon operator construction; fermionic case).
- J. Dudek et al., PRD 82 (2010) 034508, arXiv:1004.4930 (integer-spin subduction tables).
- Grid: Boyle et al., arXiv:1512.03487 (FFT class `Grid/algorithms/FFT.h`, `FFT_dim_mask`).
- GEVP: Blossier et al., JHEP 0904:094, arXiv:0902.1265.

## Files

- **Create** `Grid_sdm_build/src/gauge_gen_Nc4/two_baryon_irrep_ops_impl_plan_claude.md` (this
  plan, Chunk 0).
- **Create** `Grid_sdm_build/src/gauge_gen_Nc4/lattice_irreps_claude.py` -- group theory module
  (Chunk 1). Pure numpy.
- **Create** `Grid_sdm_build/src/gauge_gen_Nc4/two_baryon_momblock_claude.h` -- C++ baryon
  blocks + FFT + HDF5 momentum-block output (Chunk 2). Header-only so the peer's driver includes
  it; no edits to `two_baryon_gevp_contract_claude.cc` by this session.
- **Create** `Grid_sdm_build/src/gauge_gen_Nc4/two_baryon_momblock_test_claude.cc` -- small
  standalone validation driver (brute-force double sum vs blocks on a tiny cold q00 dump).
- **Create** `Grid_sdm_build/src/gauge_gen_Nc4/two_baryon_irrep_corr_claude.py` -- builds the
  irrep-projected correlator matrices from the stored momentum blocks (Chunk 3).
- **Create** `tmp_claude.sh` (repo root) for the peer/user to compile + run the GPU test.

## Implementation chunks

### Chunk 0 -- plan file
Files: `two_baryon_irrep_ops_impl_plan_claude.md`
Copy of this plan into the source dir; send the peer a summary + the two open questions.

### Chunk 1 -- group theory module `lattice_irreps_claude.py`
Files: `lattice_irreps_claude.py`
- `Oh()`: the 48 signed permutation $3\times3$ integer matrices; `det` = parity; class label from
  (trace of the proper part, axis type): $E$, $C_4$, $C_2$, $C_3$, $C_2'$ and their improper
  partners $I$, $S_4$, $\sigma_h$, $S_6$, $\sigma_d$. Character table hardcoded for $O_h$ (10
  irreps) with g/u by parity; orthogonality checked in a self-test.
- `little_group(P)`: stabiliser $\{g: gP=P\}$; identifies $O_h$, $C_{4v}$, $C_{2v}$, $C_{3v}$ by
  order + generators; hardcoded character tables ($C_{4v}$: $A_1,A_2,B_1,B_2,E$; $C_{2v}$:
  $A_1,A_2,B_1,B_2$; $C_{3v}$: $A_1,A_2,E$) with class identification done geometrically (axis
  of rotation parallel to $P$; mirror plane contains $P$ and a coordinate axis => $\sigma_v$,
  else $\sigma_d$). The $B_1/B_2$ and $\sigma_v/\sigma_v'$ labelling convention is fixed to
  Goeckeler et al. 1206.4141 and documented in the module docstring.
- `subduction(J, P)`: multiplicities of little-group irreps in continuum spin $J$ (rest frame:
  $\chi^{(J)}(\theta)=\sin((2J+1)\theta/2)/\sin(\theta/2)$, parity $\eta(-1)^{?}$ handled
  explicitly for proper/improper elements; moving frames: helicity subduction). Reproduces
  $J{=}2\to E_g+T_{2g}$, $J{=}4\to A_{1g}+E_g+T_{1g}+T_{2g}$ as self-tests.
- `pair_basis(P, shell)`: unordered pairs $\{p_1,p_2\}$ with $p_1+p_2=P$ and
  $p_1^2+p_2^2=$ shell (integer, units $(2\pi/L)^2$), as a list; `regular_rep(g)` = permutation
  matrix on that basis.
- `projector(P, irrep, basis)`: $P^\Lambda=\frac{d_\Lambda}{|G|}\sum_g\chi_\Lambda(g)^*R(g)$
  (orbital only, spin spectator per fact 4); returns an orthonormal basis of its column space
  (rank = multiplicity $\times d_\Lambda$), so each projected sink operator is a coefficient vector
  $c_{(p_1,p_2)}$ over the pair basis. Self-tests: $P^\Lambda$ idempotent, $\sum_\Lambda P^\Lambda=1$,
  ranks reproduce the shell table above.
- `source_orbit(P, pair)`: orbit of a corner pair $(a,b)$ under the little group about the
  origin, with the translation phase $e^{iP\cdot s}$ (fact 8), for optional source averaging.
- Tier-B hook (documented, not implemented): `wignerD(2, g)` for rotating the spin index once
  $q_{10}$ exists.
- `python3 lattice_irreps_claude.py` runs all self-tests and prints the operator tables for
  $P^2\in\{0,1,2,3\}$, shells up to $p_1^2+p_2^2\le 6$.

### Chunk 2 -- C++ momentum blocks `two_baryon_momblock_claude.h` + test driver
Files: `two_baryon_momblock_claude.h`, `two_baryon_momblock_test_claude.cc`, `tmp_claude.sh`
- `SubsetTable`: the 70 4-subsets of $\{0..7\}$ in lexicographic order, complement index, sign
  $\epsilon_S$. Written to HDF5 (`subsets` [70][4], `subsetSign` [70]).
- `BaryonBlocks(qa, qb, D)`: for a source pair $(a,b)$, `unvectorizeToLexOrdArray` both
  `LatticeColourMatrix` to host, per site build the $4\times8$ matrix $[q^a|q^b]$, compute the 70
  $4\times4$ minors with a fixed-size LU (`det4`), `vectorizeFromLexOrdArray` into 70
  `LatticeComplex`. Uses `Nc` from Grid (asserted $=4$). Check `D[0123] == BaryonSingleDensity(qa)/24`.
- `MomentumBlocks(D, plist, MB)`: `FFT theFFT(UGrid); FFT_dim_mask(Dt, D, {1,1,1,0}, FFT::forward)`
  (FFTW forward = $e^{-ipx}$), then `sliceSum`-free peek: for each $p$ in `plist` and each $t$,
  `peekSite` at $(p_x,p_y,p_z,t)$ with negative components wrapped mod $L$. Output
  `MB[S][ip][t]` complex.
- Momentum list: all $p$ with $p^2\le p^2_{\max}$ (default 4: 33 momenta), written as
  `momenta` [Np][3].
- HDF5 datasets (same `Hdf5Writer` compound (re,im) convention as the driver):
  `MB_set<s>_pair<a><b>` shape [70][Np][T]; `CBmom_set<s>_corner<Y>` [Np][T] (=
  $24\tilde D_{0123}$ per corner, with the corner phase left to python); tables `subsets`,
  `subsetSign`, `momenta`, `pairs`. Everything else (projection, $P$ selection) is python.
- Cost: per source pair 70 4x4-dets/site + 70 FFTs of a $24^3\times48$ complex field. 3 pairs
  ~210 FFTs; all 28 pairs ~2000 FFTs, still minutes. Memory: process one pair at a time.
- Test driver: reads a q00 `.lime`, runs the blocks, and for $p_1,p_2\in\{0,\pm e_z\}$ compares
  against the brute-force $\sum_{x_1,x_2}e^{-ip_1x_1-ip_2x_2}24^4\det Q$ (`TwoBaryonCorr` per
  site pair) on a $4^4$ cold dump. Also checks $C(p_1{=}0,p_2{=}0)$ with source pair $(a,a)$ is 0
  and that the momentum-0 single-baryon block equals `sliceSum(BaryonSingleDensity)`.
- `tmp_claude.sh`: compile via `compile_two_baryon_claude.sh two_baryon_momblock_test_claude.cc`,
  dump a $4^4$ cold q00 with the dumper (`--grid 4.4.4.4 --width 0`, its cold output name is
  hardcoded `ColdConfig.prop.lime`), run the test, tee `two_baryon_momblock_test_claude.log`.
  Peer/user runs it (GPU binary; device-0 co-tenancy rule).
- Integration (peer): include the header in `two_baryon_gevp_contract_claude.cc`, call
  `BaryonBlocks`+`MomentumBlocks` per (set, pair) after the existing Chunk-2 loop, add
  `--pmax2 <n>` and `--pairs all|anchored|class` CLI. Coordinate via SendMessage.

### Chunk 3 -- projected correlators `two_baryon_irrep_corr_claude.py`
Files: `two_baryon_irrep_corr_claude.py`
- Reads per-config `MB_*` datasets (or the peer's combined per-ensemble h5, schema agreed in
  Chunk 2), builds for each requested $(P,\Lambda,\text{shell},\text{mult})$ sink operator
  $$C_{\Lambda;(a,b)}(t)=24^4\sum_{\{p_1,p_2\}}c^\Lambda_{(p_1,p_2)}\sum_S\epsilon_S\tilde D_S(p_1,t)\tilde D_{\bar S}(p_2,t),$$
  for every source op (set, pair), optionally source-orbit-averaged with phases. Writes
  `C2Birrep_<P>_<irrep>_shell<n>_m<k>_src<...>` into an h5 the peer's `gevp_h5_claude.py` can
  aggregate, plus $C_B(p^2,t)$ shell-averaged single-baryon correlators for $E_B(p^2)$.
- Default operator set (rest frame): $A_{1g}$ at shells $p^2=0,1,2,3$ (4 sink ops) x source
  (3 anchored pairs x sets). Moving frames: $P=(0,0,1)$ $A_1$ at shells $(0,1)$, $(1,2)$;
  $P=(1,1,0)$ $A_1$; $P=(1,1,1)$ $A_1$.
- Symmetry check on real data: correlators of pairs in the same orbit agree within errors.

### Chunk 4 -- (after data) GEVP + levels
Files: python, with the peer (their GEVP code). Effective energies per irrep, $E_B(p^2)$
dispersion (needed because $aM_B\approx3$: use the measured lattice dispersion, not continuum,
in the Luescher $q^{*2}$), $\Delta E$ per level. Luescher zeta functions are a later plan.

## Verification

1. `python3 lattice_irreps_claude.py`: self-tests pass (group orders 48/8/4/6, character
   orthogonality, projector idempotence/completeness, shell table above, subduction of $J\le4$).
2. `tmp_claude.sh` (peer/user runs): momentum blocks vs brute-force double sum on $4^4$ cold
   agree to ~1e-10 relative; block $0123$ matches `BaryonSingleDensity`.
3. On one real config (lat.1000, already local): $C_{A_{1g}}(p^2{=}0)$ sink op vs the existing
   Chunk-8 $P=0$ smeared-sink correlator have consistent plateaus; single-baryon $E_B(p^2)$
   rises with $p^2$; orbit-mates agree within noise.

## Open questions (for the user / peer)

1. **Tier B (full spin projection)?** Storing $q_{10}$ (sink PD spin 1, source spin 0) in the
   dumper is zero extra solves and 2x file size but requires re-dumping. Needed only to separate
   $E/T_2$-type spin splittings, not for S-wave $S=4$ levels. Recommendation: defer; design the
   python module with the hook.
2. **Source-pair table for the momentum blocks:** 3 anchored class reps (as v1), 7 anchored, or
   all 28 pairs (allows full source orbit averaging: 12 edge + 12 face + 4 body)? Cost is FFTs
   only. Recommendation: all 28 in production, 3 for validation.
3. $p^2_{\max}$ for stored momenta: 4 (33 momenta) suffices for rest-frame shells $\le4$ and
   moving frames $|P|^2\le3$ with shells $\le6$. OK?

# Two-baryon correlator: implementation plan

## Physics / goal

SU($N_c=4$), one-flavour SDM theory, $24^3\times 48$ ensembles. Measure the
correlator of a **two-baryon** state in which the two baryons sit at the two
most distant spatial points of the periodic box:

- $O$ = spatial origin $(0,0,0)$,
- $M$ = spatial center $(12,12,12)$ (the furthest point from $O$ under periodic
  boundary conditions).

Each baryon is the spin-2, all-spin-up, positive-parity operator (the `0000`
operator of `thermo_gw_details` eq. 1.2.13), built from four quarks in the
upper Pauli-Dirac (nonrelativistic) component `0`:

$$
B(x) = \varepsilon_{abcd}\, q_0^a(x)\, q_0^b(x)\, q_0^c(x)\, q_0^d(x) .
$$

We loop over a time separation $t$ between a source set and a sink set of the
two baryons and measure the correlator $C_{2B}(t)$, to extract the two-baryon
energy at the maximal separation $r=|O-M|$ and hence (with $2M_B$ from the
single-baryon code) the two-baryon interaction at that separation.

## Conventions reused from the single-baryon code

- Action: `MobiusFermionD`, $L_s=16$, scaled Shamir $b=1.5$, $c=0.5$,
  anti-periodic temporal BC `{1,1,1,-1}`. $M_5$, mass from command line.
- Solver: `SchurRedBlackDiagMooeeSolve` + `ConjugateGradient` (tol $10^{-9}$),
  one point-to-all propagator per point source (12 inversions).
- Weyl$\to$Pauli-Dirac rotation `q_DP = 0.5 * U * q * Udagger` (eq. 1.2.3-1.2.6),
  then the (spin) `0` component is a `LatticeColourMatrix`
  `q00 = peekSpin(q_DP,0,0)`. For the `0000` baryon all four quark lines use
  `q00`.
- Single-baryon validation target (eq. 1.2.17):
  $$
  C^{(0000)}(t) = 24\,\varepsilon_{abcd}\varepsilon_{a'b'c'd'}\,
  S_{00}^{aa'} S_{00}^{bb'} S_{00}^{cc'} S_{00}^{dd'} ,
  $$
  exactly the weight-24 loop in `baryons_0000_dirac.cc:254`.

## Geometry (CONFIRMED)

Standard two-baryon (NN-type) correlator. Notation: unprimed = source-side
spatial point on the source timeslice $t_{\rm src}$ (fixed, $=0$); primed =
sink-side spatial point on the sink timeslice $t_{\rm snk}=t$.

- Source at $t_{\rm src}=0$: $\bar B(O)\,\bar B(M)$, with $O=(0,0,0)$,
  $M=(12,12,12)$.
- Sink at $t_{\rm snk}=t$: $B(O')\,B(M')$, with $O'=(0,0,0)$, $M'=(12,12,12)$.
- Loop $t$ over all sink timeslices $1\ldots N_t-1$.

$$
C_{2B}(t) = \big\langle B(O',t)\,B(M',t)\;\bar B(O,0)\,\bar B(M,0) \big\rangle .
$$

Every quark line runs $0\!\to\!t$. The four propagator types that appear are
$O\!\to\!O'$, $M\!\to\!M'$ (direct) and $O\!\to\!M'$, $M\!\to\!O'$ (exchange);
**no same-timeslice line** (the earlier "(O-M) same timeslice" wording was a
mistake). Propagators needed: point-to-all from $O$ and from $M$ on $t_{\rm src}$,
evaluated at $O'$ and $M'$ on every sink timeslice. **Two solves total**, reused
for the whole $t$ range.

Resolved choices:
- $M$ = spatial center $(12,12,12)$ on the source timeslice (same slice as $O$).
- Fixed source timeslice $t_{\rm src}=0$; scan all sink timeslices.
- Position-space fixed sink points $O',M'$ (no momentum projection).
- Channel: `0000` ($J^P=2^+$, all spin up) for both baryons only.
- Full antisymmetrisation of the 8 identical quarks (all exchange terms).

## Contraction (determinant form)

All eight quark lines use the Pauli-Dirac spin-`0` color matrix. Let $P_{XY}$ be
the $4\times4$ color matrix (row = sink color, col = source color) of the
spin-`0` propagator from source point $Y\in\{O,M\}$ to sink point
$X\in\{O',M'\}$ — i.e. `peekColour(q00_Y, ., .)` evaluated at site $X$. There are
four such matrices: $P_{O'O},P_{O'M},P_{M'O},P_{M'M}$.

The sink baryon at $O'$ carries $\varepsilon_{abcd}$, the one at $M'$ carries
$\varepsilon_{efgh}$; the source baryon at $O$ carries $\varepsilon_{a'b'c'd'}$,
the one at $M$ carries $\varepsilon_{e'f'g'h'}$. With all color indices fixed by
these four $\varepsilon$'s, Wick's theorem for the eight identical quarks is the
determinant of the $8\times8$ propagator matrix $Q$ over the 8 quark lines:

$$
C_{2B}(t) = \sum_{\substack{abcd\\efgh}}\sum_{\substack{a'b'c'd'\\e'f'g'h'}}
\varepsilon_{abcd}\,\varepsilon_{efgh}\,
\varepsilon_{a'b'c'd'}\,\varepsilon_{e'f'g'h'}\;\det Q ,
$$

where $Q$ is indexed by the 8 sink quark lines (rows: 4 at $O'$ with colors
$a,b,c,d$, then 4 at $M'$ with colors $e,f,g,h$) and the 8 source quark lines
(cols: 4 at $O$ with colors $a',b',c',d'$, then 4 at $M$ with colors
$e',f',g',h'$):

$$
Q_{(\text{sink }X,\,c_{\rm sink})\,(\text{source }Y,\,c_{\rm src})}
 = (P_{XY})_{c_{\rm sink}\,c_{\rm src}} .
$$

The determinant sums over all $8!$ quark pairings with the correct fermion
signs, automatically including the direct and all exchange ($1$- to $4$-quark)
topologies. The four $\varepsilon$'s project each vertex onto the SU(4) color
singlet.

### Key simplification (no $24^4$ color loop)

$\det Q$ is alternating in its 8 columns and 8 rows. Summing over each
$\varepsilon$ just permutes a block of 4 columns (or rows), so each
$\varepsilon$ sum collapses to a factor $4! = 24$ times the determinant
evaluated at the fixed color assignment $\{0,1,2,3\}$ for that block:

$$
C_{2B}(t) = 24^4 \,\det
\begin{pmatrix} P_{O'O} & P_{O'M} \\ P_{M'O} & P_{M'M} \end{pmatrix},
$$

where each $P_{XY}$ is the full $4\times4$ color matrix (rows = sink color
$0..3$, cols = source color $0..3$). So the entire contraction is a single
$8\times8$ complex determinant per sink timeslice — cheap. With
$P_{O'M}=P_{M'O}=0$ (no exchange) this factorizes to
$24^4 \det P_{O'O}\det P_{M'M} = C^{(0000)}_O \cdot C^{(0000)}_M$, the product of
the two single-baryon correlators (each $= 576\det P = 24\times$ the weight-1
`baryons_0000_dirac` value), fixing the normalization.

**Validation (Chunk 3).** A brute-force reference sums the four $\varepsilon$'s
explicitly over $24^4$ color assignments (inner $8\times8$ determinant via
Eigen) and is checked against the $24^4\det Q_{8\times8}$ formula on random
$4\times4$ blocks; they must agree to machine precision. The factorization limit
is checked against the validated single-baryon result.

**Self-check (single baryon).** Dropping the $M$ baryon, the same construction
reduces to $\varepsilon_{abcd}\varepsilon_{a'b'c'd'}\det Q_{4\times4}$; the
determinant's $4!=24$ permutations are all identical under $\varepsilon_{a'b'c'd'}$,
giving $24\,\varepsilon_{abcd}\varepsilon_{a'b'c'd'}S_{00}S_{00}S_{00}S_{00}$ =
eq. (1.2.17). So the determinant form reproduces the published single-baryon
weight-24 normalisation, which is the Chunk 2 validation.

**Cost.** The contraction is evaluated only at the two sink points per
timeslice (not over the full volume): $\sim 24^4$ nonzero color combinations
$\times$ an $8\times8$ determinant $\times N_t$ timeslices $\approx 10^{10}$
complex ops per config — acceptable for a naive first version. Chunk 7 may
pre-antisymmetrise the source columns (collapsing the $24^2$ source sum into
baryon blocks) if speed matters.

## Files

All new files, `_claude` suffix, in `Grid_sdm_build/src/gauge_gen_Nc4/`:

- `two_baryon_corr_claude.cc` — main driver (sources, solves, contraction, IO).
- (optional) `two_baryon_contract_claude.h` — the contraction kernel as a
  reusable header if it grows large.
- `test/run_two_baryon/run_claude.sh` — run script (mirrors
  `test/run_baryon/run.sh`).
- `test/run_two_baryon/compare_claude.py` — HDF5 regression compare (reuse
  existing `compare.py`).

## Ordered implementation chunks (Candidate A/C)

### Chunk 1 — Skeleton, sources, propagators
Build the driver: parse `config M5 mass outfile`, set up grids/action exactly as
`baryons_0000_dirac.cc`, create point sources at $O$ and at $M$, run two
`Solve()` calls to get `PropO`, `PropM` (point-to-all). Rotate both to PD basis
and extract the spin-`0` color matrices `q00_O`, `q00_M`.
Files: `two_baryon_corr_claude.cc`.

### Chunk 2 — Single-baryon validation block
From `q00_O` reproduce the single-baryon $C^{(0000)}(t)$ via the
$\varepsilon_{abcd}\varepsilon_{a'b'c'd'}\det Q_{4\times4}$ form and check it
against `baryons_0000_dirac` output on the same config (must match including the
weight-24 factor). Confirms the propagator, PD-rotation, $\varepsilon$ machinery,
and the determinant kernel before the two-baryon contraction.
Files: `two_baryon_corr_claude.cc`.

### Chunk 3 — Two-baryon contraction kernel
Implement $C_{2B}$ via the $8\times8$ determinant form using the four propagator
color matrices $P_{O'O},P_{O'M},P_{M'O},P_{M'M}$ obtained by `peekSite` of
`q00_O`,`q00_M` at sink points $O',M'$. Validate on a cold config (where the
two-baryon answer factorises into a product of single-baryon correlators up to
exchange) and against a small-volume brute-force pairing reference.
Files: `two_baryon_corr_claude.cc` (+ `two_baryon_contract_claude.h` if large).

### Chunk 4 — Time loop and correlator assembly
For fixed source timeslice, evaluate $C_{2B}$ at every sink timeslice to get
$C_{2B}(t)$ vs $t$ (sink points $O',M'$ at each $t$). Assemble the $t$ array.
Files: `two_baryon_corr_claude.cc`.

### Chunk 5 — HDF5 output + run script
Write $C_{2B}(t)$ to HDF5 with the `MesonFile`-style serialisable container;
add `test/run_two_baryon/run_claude.sh`.
Files: `two_baryon_corr_claude.cc`, `test/run_two_baryon/run_claude.sh`.

### Chunk 6 — Test / regression
Run on cold config and one real config; sanity-check $t$-dependence
(two-baryon effective mass $\to 2M_B$ at large $t$, large separation).
Files: `test/run_two_baryon/`.

### Chunk 7 (optional) — Optimise contraction
If the naive contraction is too slow, replace with the block / unified
contraction algorithm. Deferred until correctness is established.

## Open questions

All resolved (see "Geometry (CONFIRMED)"). Remaining minor item, can default
without blocking: exact HDF5 dataset name/layout for $C_{2B}(t)$ — I will reuse
the `MesonFile` serialisable container with key `two_baryon_0000_t`.

# two_baryon_gevp_contract — implementation plan

Reads the stored `q00` point-to-all propagator records (from
`point2all_prop_dumper_claude.cc`) and builds the two-baryon **GEVP correlator matrix**
$C_{ij}(t)$ plus the single-baryon correlator $C_B(t)$, for the SU($N_c{=}4$) SDM
scattering study. New file: `two_baryon_gevp_contract_claude.cc`.

## Goal / physics

The 8 corner sources $\{0,L_\mu/2\}^3$ (indices 0..7, spatial coords = bits of the index;
displacement class of a pair $\{a,b\}$ = `popcount(a XOR b)`: 1 = edge $L/2$, 2 = face
diagonal $L/\sqrt2$, 3 = body diagonal $\sqrt3 L/2$). A two-baryon **operator** is a pair
of corners. The correlator matrix
$$
C_{ij}(t) = \langle\, O_i(t)\, O_j^\dagger(0)\,\rangle,\qquad
O_j^\dagger(0)\sim \bar B(a_j)\,\bar B(b_j),\ \ O_i(t)\sim B(c_i)\,B(d_i),
$$
feeds the GEVP $C(t)v_n=\lambda_n(t,t_0)C(t_0)v_n$, $\lambda_n\to e^{-E_n(t-t_0)}$, to
resolve the low-lying two-baryon spectrum. Source and sink operator sets need not match
(unmatched GEVP, see `point2all_prop_dumper_impl_plan_claude.md`).

Interaction energy needs the single baryon: $\Delta E = E_{2B}-2M_B$, $M_B$ from
$C_B(t)$. Luescher analysis is downstream (python).

## References (mandatory citations)

Cite these in the implementing code comments as well.

- **Grid library**: P. Boyle, A. Yamaguchi, G. Cossu, A. Portelli, "Grid: A next generation
  data parallel C++ QCD library", arXiv:1512.03487.
- **Two-baryon "0000" contraction (v1)**: fully-antisymmetrised correlator collapses to a
  single $8\times8$ determinant of the four propagator colour blocks,
  $C_{2B}=24^4\det[[P_{c a},P_{c b}],[P_{d a},P_{d b}]]$, $P_{X'Y}=q_{00}^Y(X')$. Derivation
  in this project: `two_baryon_impl_plan_claude.md`. Reuse
  `detLeibniz`/`TwoBaryonCorr`/`assemble8x8`/`colourBlockAt`/`BaryonSingleDensity` verbatim
  from `two_baryon_corr_prod_claude.cc`.
- **GEVP**: B. Blossier, M. Della Morte, G. von Hippel, T. Mendes, R. Sommer (Della Morte,
  Garcia i Tormo, Sommer, Weber), "On the generalized eigenvalue method for energies and
  matrix elements in lattice field theory", JHEP 0904:094 (2009), arXiv:0902.1265. Unmatched
  (asymmetric, source set $\ne$ sink set) extension: `point2all_prop_dumper_impl_plan_claude.md`.
- **Multi-baryon block / determinant contraction (v2 momentum projection)**:
  W. Detmold, K. Orginos / M. J. Savage et al., "Nuclear physics from lattice QCD", the
  determinant contraction method, Phys. Rev. D 82 (2010) 014511, arXiv:1001.2768;
  T. Doi, M. G. Endres, "Unified contraction algorithm for multi-baryon correlators on the
  lattice", Comput. Phys. Commun. 184 (2013) 117, arXiv:1205.0585. The double sink sum
  factorises through the block structure of the $8\times8$ det into per-sink-point
  momentum-projected baryon blocks (FFT, $V\log V$), avoiding the $V^2$ brute force.
- **SciDAC/LIME IO**: Grid `ScidacReader` (`Grid/parallelIO/IldgIO.h`).

## Staging (per user)

- **v1 (this plan): fixed-point sink.** Sink baryons at fixed corner points $(c_i,d_i)$ at
  time $t$; no momentum projection. Reuses the $8\times8$ det directly. Validates against the
  old `two_baryon_corr_prod_claude.cc` cold-config number for the $\{0,7\}$ element.
- **v2 (later): momentum-projected sink** via single-baryon momentum blocks (FFT over sink
  $x$), definite total $P{=}0$ + relative momentum. Separate plan.

## Input format (from the dumper, confirmed)

One `.lime` per config, **8 records**, each a `LatticeColourMatrixF` (single precision;
row = sink colour, col = source colour) with a `PropRecord` metadata
{config, srcCoord[4], srcIndex, mass, M5, b, c, Ls, smearWidth, smearNiter, store="q00",
prec="single"}. Records written in order `srcIndex` 0..7. `q00` holds all sink times.

## Operator basis (v1) and asymmetric / multi-smear structure

Independent **source** and **sink** operator tables (unmatched GEVP: source set need not
equal sink set; $C$ may be rectangular `nSnk x nSrc x T`).

A **source operator = (source set s, corner pair {a,b})**, where a "source set" is one input
`.lime` file = one source-smearing variant. This makes source smearing a first-class knob: v1
production has one set (Gaussian $W{=}3$), and a **point source (`--width 0`) set can be added
later just by passing a second file** -- no code change. So `q00[s][corner]`, and the source
basis rank = (#sets) x (#pairs), the natural way to grow rank and resolve more levels.

A **sink operator = corner pair** at fixed sink points (v1 point/unsmeared sink; single set).

Default corner-pair list (both src and snk) = one representative per displacement class,
anchored at corner 0: `{0,1}` (edge), `{0,3}` (face diag), `{0,7}` (body diag) -> 3 pairs.
With one source set that is a `3x3` matrix; with a second (point) source set it becomes
`3 (snk) x 6 (src)`. `{0,7}` is O/M and reproduces the old two-baryon number (validation
anchor). Tables are explicit in the code, easily extended.

Note: v1 real-config data is **smeared-source / point-sink** = inherently asymmetric; energies
still correct via the unmatched GEVP (non-normal $C$; check $t_0$-stability).

## Files

- **Create** `two_baryon_gevp_contract_claude.cc`.
- **Create** this plan (done).
- **Create later** `rsync_pull_prop_dumper_claude.sh` (pull a few real `.lime` for dev) and a
  cluster run script `run_gevp_contract_*_claude.sh` (after v1 validates).
- No edits to the dumper or the existing two-baryon driver.

## Implementation chunks

### Chunk 1 — grids + SciDAC reader
Files: `two_baryon_gevp_contract_claude.cc`
- Grid init, 4D UGrid (single-precision grid `UGridF` to read the `LatticeColourMatrixF`
  records; promote to double `LatticeColourMatrix q00[8]` for the host-side det).
- `PropRecord` struct (copy from dumper) to read metadata.
- `readQ00(infile, q00, recs)`: `ScidacReader`, loop 8 `readScidacFieldRecord`, map by
  `srcIndex`, sanity-print corner coords + norms. Assert 8 records, store=="q00".
- CLI: `<infile.lime> <outfile.h5>` (+ `--grid`/`--mpi`). Cold default = `ColdConfig.prop.lime`.

### Chunk 2 — two-baryon fixed-point C_ij(t)
Files: `two_baryon_gevp_contract_claude.cc`
- Copy `detLeibniz`, `assemble8x8`, `TwoBaryonCorr`, `colourBlockAt` from
  `two_baryon_corr_prod_claude.cc`.
- Operator table `std::vector<std::array<int,2>> ops` (default 3 pairs above).
- For each sink time $t$, each op pair $(i,j)$: sink points $(c_i,t),(d_i,t)$, source corners
  $(a_j,b_j)$; four colour blocks $P_{c_i a_j}=$`colourBlockAt(q00[a_j],(corner c_i,t))` etc.;
  `C2B[i][j][t] = TwoBaryonCorr(...)`.
- **Validation**: on `ColdConfig.prop.lime`, the $\{0,7\}\times\{0,7\}$ element must match
  the old code's cold `two_baryon_0000_t` (both unsmeared). Print max abs diff.

### Chunk 3 — single-baryon C_B(t) + HDF5 output
Files: `two_baryon_gevp_contract_claude.cc`
- `BaryonSingleDensity(q00[Y])` per corner Y -> `sliceSum` over Tdir -> $C_B^Y(t)$ (p=0);
  also store the corner-averaged $\bar C_B(t)$.
- HDF5 (`Hdf5Writer`) datasets: `C2B_i_j_t` per op pair (or a flattened `[nop*nop][T]` +
  the `ops` table), `CB_Y_t` per corner, `CB_avg_t`, and a metadata record (config, mass,
  M5, smear, ops). Matches the two-baryon h5 conventions.

## v1.5 — zero-momentum ($P=0$) projected sink operators (in addition to point sink)

Add total-momentum-zero sink operators alongside the fixed-point ones. Two possible
meanings of "zero momentum":

**(A) COM-projected, fixed relative displacement $d$ (RECOMMENDED, tractable now).**
$$
O_d(t) = \sum_{\mathbf X} B(\mathbf X)\,B(\mathbf X+\mathbf d),\qquad
C_d^{P=0}(t) = \sum_{\mathbf X} 24^4\,\det Q(\mathbf X,\mathbf X+\mathbf d),
$$
where $Q$ is the same $8\times8$ block matrix as the fixed-point case, now with sink baryons
at $\mathbf X$ and $\mathbf X+\mathbf d$ (get $q_{00}(\mathbf X+\mathbf d)$ by `Cshift`). A
**single** COM sum $\Rightarrow$ total momentum $P=0$; the relative separation $\mathbf d$ is
fixed (a position-space relative wavefunction $\delta(\mathbf r-\mathbf d)$). Different
$\mathbf d$ = different operators = a clean $P=0$ variational basis. **No plane-wave phases,
no baryon-block factorisation** -- just a per-site $8\times8$ determinant summed over the
lattice. This is the natural $P=0$ version of the current displacement operators and should
resolve the $P=0$ tower better than the fixed-point (all-total-momenta) operators.

**(B) Both baryons individually at rest, $B(\mathbf p{=}0)\,B(\mathbf p{=}0)$.** This is a
DOUBLE sink sum $\sum_{\mathbf X_1,\mathbf X_2}$, which needs the baryon-block factorisation
(the $8\times8$ det factorises into per-sink-point blocks) -- i.e. the same machinery as the
general momentum projection (v2). Deferred with v2.

We do **(A)**. Sink displacements $\mathbf d$ = the 3 corner classes:
$(L/2,0,0)$ edge, $(L/2,L/2,0)$ face, $(L/2,L/2,L/2)$ body. Keep the 3 fixed-point sink ops
too, so the sink basis grows to 6 (3 point + 3 $P{=}0$); source ops stay the 3 fixed corner
pairs (source can't be projected -- it's where the solves were placed). $C$ becomes
$6_{\rm snk}\times3_{\rm src}$; the GEVP downstream picks a square block (e.g. the 3 $P{=}0$
sink vs 3 source).

### Implementation (Chunk 7)
Files: `two_baryon_gevp_contract_claude.cc`
- Fast `detLU8` (8x8 complex LU, partial pivoting, fixed array) -- Leibniz is too slow per
  site.
- `TwoBaryonP0(qA, qB, dcoord)`: `Cshift` $q_{00}$ blocks by $\mathbf d$; unvectorise the 4
  fields to host lex-order arrays (`unvectorizeToLexOrdArray` / `peekLocalSite`); loop sites,
  assemble $8\times8$, `detLU8`, accumulate $24^4\det$ into $C[t]$ (t = site time). Global-sum
  $C[t]$ over ranks (single-rank in dev).
- For each source op $(a,b)$ and each sink $\mathbf d$: `C2Bp0_snkd<k>_src<jo>` dataset.
  Keep existing `C2B_snk*_src*` (point). Record the $P{=}0$ sink displacement list in `meta`.
- Cost: ~$V\times$(#src ops)$\times$(#d) $8\times8$ LU dets/config -- ~seconds, host-side OK.

## Chunk 8 — smeared sink (DONE)

Add covariant Gaussian (Wuppertal) smearing on the **sink** index of $q_{00}$, giving
smeared-sink two-baryon and single-baryon correlators alongside the point-sink ones.

- Sink smearing = smear each quark propagator at the sink: $q_{00}(x)\to\sum_y\phi(x,y)q_{00}(y)$,
  gauge-covariant on the sink (row) colour index. `CovariantSmearing<PeriodicGimplD>::
  GaussianSmear(U, q00, w, N, Tdir)` does exactly this (same call as the source smearing, now on
  the stored $q_{00}$). The baryon is then $\varepsilon\,q_{00}^{\rm smeared}(x)^4$ (quark-smeared
  sink, standard). Needs the **gauge field**, so the contraction now loads the NERSC config.
- CLI: `--config <nersc>` (gauge field for this config; cold if omitted) and
  `--sink-smear <w> <N>`. When both given, also emit `C2Bss_snk<i>_src<jo>`, `CBss_set<s>_*`.
- Implementation: the Chunk-2 fixed-point loop and Chunk-3 single-baryon loop are refactored into
  `ContractFixedPoint(q,...)` and `SingleBaryonP0(q,...)`, called once for the raw $q_{00}$ (point
  sink) and once for the sink-smeared $q_{00}^s$ -> one validated code path.
- Gives the source x sink smearing matrix (SS/SP/PS/PP across the source sets + sink types) = a
  smearing variational basis for the GEVP.
- Logistics (handoff): the contraction now needs the gauge config (lustre5). Run where both
  $q_{00}$ and the gauge configs are visible (tuolumne/lustre5), or stage the configs; avoids the
  ~525 GB gauge stage to dane if run on lustre5.

## Open questions

1. **Operator basis size** — 3 class-representative pairs (v1 default) or a larger set
   (all pairs anchored at 0 = 7 ops, or all 28)? Bigger = richer GEVP but noisier / more
   redundant. Start with 3?
2. **Single-baryon momenta** — v1 does p=0 only; nonzero-p single baryon waits for the v2
   FFT machinery. OK?
3. **Precision** — read single, contract the det in double (peek -> ComplexD). Fine.

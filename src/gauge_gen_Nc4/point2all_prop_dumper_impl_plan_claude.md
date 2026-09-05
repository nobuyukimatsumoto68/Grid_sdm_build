# point2all_prop_dumper — implementation plan

Working file: `point2all_prop_dumper_claude.cc` (renamed from `point2all_prop_dumper.cc`,
which is kept untouched as the seed).

## Goal / physics summary

Compute and dump **point-to-all propagators** from a set of point sources on the
source timeslice $t_\text{src}=0$, for later reuse in two-baryon (and single-baryon)
contractions of the SU($N_c=4$) one-flavour SDM theory. The dumper does no
contraction itself — it only solves $D\,q = \delta$ and serialises the propagators.

Baryon operator downstream is the spin-2, positive-parity "0000" operator
$B(x)=\varepsilon_{abcd}\,q_0^a q_0^b q_0^c q_0^d$; the two-baryon correlator places
baryons at two distant spatial points (origin $O$ and box centre $M$). This dumper
generalises to the **8 spatial corners** $\{0,L_\mu/2\}^3$ (all with $t=0$), so any
$O/M$ pair (or corner-averaged statistics) can be contracted downstream.

Action: Mobius DWF, scaled Shamir kernel $b=1.5,\ c=0.5$ ($b+c=2,\ b-c=1$),
$L_s=16$, $M5$/mass from argv, anti-periodic in time $\{1,1,1,-1\}$. Schur red-black
preconditioned CG, tol $10^{-9}$.

## Method / algorithm source

**Split-grid multi-RHS (MRHS)**: partition the MPI layout into `nrhs` sub-communicators
via `Grid_split` / `Grid_unsplit`; each sub-grid holds the *full* lattice on fewer
ranks and solves one RHS, all concurrently. This is the Grid idiom for batching many
propagator solves on GPUs.
- Reference: P. Boyle et al., "Grid: data parallel C++ ...", arXiv:1512.03487; canonical
  usage in `Grid/tests/solver/Test_dwf_mrhs_cg.cc` (`Grid_split`/`Grid_unsplit`, `SGrid`).
- Rationale (this project): at the production layouts the full-layout local volume is tiny
  ($\sim12^4$/GPU at 24c, $\sim16^3\times8$/GPU at 32c), so a whole-machine DWF solve is
  past the strong-scaling knee (comms-bound). Splitting to fewer GPUs/solve gives a much
  larger local volume and better MI300A utilisation. Per-RHS CG stays independent →
  bit-reproducible per source (no shared-Krylov coupling; `BlockCGrQ` deliberately NOT used).

## Run configurations (Livermore tuolumne, Flux, 1 GPU/task)

| lattice        | mpi (`MPIGRID`) | `MPISPLIT` | ranks/subgrid | nrhs | local vol / GPU (split) |
|----------------|-----------------|------------|---------------|------|-------------------------|
| $24^3\times48$ | 2.2.2.4 (32)    | 1.1.1.4    | 4             | 8    | $24^3\times12$          |
| $32^3\times64$ | 2.2.2.8 (64)    | 2.2.2.1    | 8             | 8    | $16^3\times64$          |

`nrhs = prod(MPIGRID)/prod(MPISPLIT)` must equal the number of sources (8). Both are env
knobs (defaults set per lattice); code asserts `nrhs == n_src`.

## IO format

- **SciDAC / LIME** via Grid `ScidacWriter` (`Grid/parallelIO/IldgIO.h`).
- **One `.lime` file per config, 8 field records** (one per source), appended by repeated
  `writeScidacFieldRecord`. Each record carries a `userRecord` XML with metadata:
  source coordinate, config name, mass, M5, b, c, Ls, t-window [tmin,tmax], precision.
- **Single precision** payload (`LatticePropagatorF` demoted copy): halves size, ample for
  spectroscopy. Container is standard LIME (readable by `lime_contents`, QIO, Chroma, other
  Grid builds); byte-exact physics reuse is Grid->Grid (index/layout conventions match).
- Per-config storage (single prec): ~10.9 GB (24c) / ~34.3 GB (32c) for the FULL time
  extent; the plateau window (below) cuts this by $T/T_\text{win}$.

## Time-window ("plateau") dump

We only need the sink-time plateau region, so dump slices $t\in[t_\text{min},t_\text{max}]$
instead of all $T$:
- Build a thin output grid `WGrid` with dims $(L_x,L_y,L_z,T_\text{win})$,
  $T_\text{win}=t_\text{max}-t_\text{min}+1$, MPI layout with **time undivided**
  (e.g. spatial = MPISPLIT spatial, time mpi = 1) so $T_\text{win}$ need not divide the
  time decomposition.
- Copy the window with `localCopyRegion(Prop, PropW, {0,0,0,tmin}, {0,0,0,0},
  {Lx,Ly,Lz,Twin})` (`Grid/lattice/Lattice_transfer.h:741`), then demote to single and
  SciDAC-write `PropW`.
- `tmin`/`tmax` are argv/env knobs; default = full range (`0..T-1`) so nothing is lost
  until a plateau is chosen.

## Files

- **Create** `point2all_prop_dumper_claude.cc` (from the seed `point2all_prop_dumper.cc`).
- **Create** this plan (done).
- **Create later** a Flux run script `run_prop_dumper_flux_claude.sh` under
  `test/run_two_baryon/` (or a new `test/run_prop_dumper/`) — separate task, after the code.
- No edits to `two_baryon_corr_claude.cc` / `two_baryon_corr_prod_claude.cc`.

## Source-side collapse: 4 solves instead of 16 (q00 mode)

The "0000" operator uses PD **source-spin 0** only, so we never need the other source spin
components. Pre-rotate the source instead of rotating the full propagator afterward.

Column $(0,c')$ of $S^{(PD)}=\tfrac12\,U\,S^{(W)}\,U^\dagger$ is obtained by solving with the
Weyl-basis source $U^\dagger e_{(0,c')}$. From `WeylToPauliDiracU` (`two_baryon_corr_prod_claude.cc:293`),
$$
U=\begin{pmatrix}1&0&1&0\\0&1&0&1\\-1&0&1&0\\0&-1&0&1\end{pmatrix},\qquad
U^\dagger e_{(0,c')} = (\text{row 0 of }U)\otimes e_{c'} = (e_0+e_2)\otimes e_{c'} .
$$
So:
- **Source** per colour $c'$: point (or Gaussian-smeared) source with spinor $(e_0+e_2)$,
  colour $c'$ → **one solve per colour = 4 solves per source point** (vs 16).
- **Sink**: $q_{00}[c,c'] = \tfrac12\big(x^{(c')}_{\text{spin }0,\,c}+x^{(c')}_{\text{spin }2,\,c}\big)$;
  both sink spins are already in the one solution.
- Floor is 4 (the baryon needs all 4 source colours for $\varepsilon_{abcd}$).
- Exact — pure linearity of $D^{-1}$; identical to `0.5*U*Prop*Udagger` then `peekSpin(0,0)`.
- Gaussian smearing commutes (spatial/colour), so it only shapes the source profile before
  the spin structure is imposed.

This 4x reduction is **tied to `q00` mode** (fixed 0000 channel). `full` mode still needs 16
solves. MRHS batching: outer loop over 4 colours, each an nrhs=8 split-solve (32 RHS total),
instead of 16 colours.

## Decisions (this session)

- **24c first** ($24^3\times48$) — moderate file size for exploration before 32c production.
- **Single Gaussian-source width** (one $w$), so $n_w=1$ and `nrhs` stays 8. Point vs single
  smeared is an operator choice; likely dump the smeared set (one width) as the source basis.
- **Q2 plateau window**: defer; fix $[t_\text{min},t_\text{max}]$ later from existing
  single-baryon data. Code keeps it a knob (default full range).
- **Q3 (batch width) + Q4 (MPISPLIT shape)**: deferred to a remote tuning agent.
- **Q6 stored object**: **v1 = `q00` mode only** (pre-rotated 4-solve, store $q_{00}$ colour
  matrices, single precision). `full` mode documented but not implemented in v1.
- **Source profile**: **single Gaussian width** per corner (point = $w{\to}0$ not stored);
  width value + thin-vs-smeared links deferred as knobs.
- **Split control**: `--split sx sy sz st` argv (Grid idiom, as in `Test_dwf_mrhs_cg`); run
  script maps its `MPISPLIT` env onto it.

## Implementation chunks

### Chunk 1 — grids + params + split-grid setup
Files: `point2all_prop_dumper_claude.cc`
- Keep existing UGrid/UrbGrid/FGrid/FrbGrid + config/M5/mass argv parsing + Mobius `D`.
- Read `MPISPLIT` (Coordinate), `tmin`/`tmax` (default full), `outfile`.
- Build `SGrid, SUrbGrid, SFGrid, SFrbGrid` at `GridDefaultLatt()` with mpi = MPISPLIT;
  derive `nrhs = prod(mpi)/prod(split)`; assert `n_src % nrhs == 0` and set
  `nbatch = n_src/nrhs`. Sources are solved in `nbatch` split-solve batches, so nrhs = 8
  (production, all concurrent) and nrhs = 1 (single-rank local validation, sequential) both
  work.
- Split the gauge field `Umu -> s_Umu` (`Grid_split`) and build the split-grid Mobius `sD`.

### Chunk 2 — sources (pre-rotated PD-spin-0, optional Gaussian smear)
Files: `point2all_prop_dumper_claude.cc`
- Fix the 8-corner source list: size each `Coordinate` to `Nd`, set spatial from bits of `i`,
  `pt[Tdir]=0`. Print all 8.
- `q00` mode source builder: for each source point and colour $c'$, a `LatticeFermion` with a
  unit at the source site in spin components $\{0,2\}$, colour $c'$ (the $(e_0+e_2)\otimes e_{c'}$
  pre-rotated source). Optional Gaussian smear of the spatial profile
  (`CovariantSmearing::GaussianSmear(U,src,width,Niter,Tdir)`) before/around the point.
- `full` mode source builder: the usual 16 spin-colour delta sources (identity
  `SpinColourMatrix`), retained behind the `--store full` path.

### Chunk 3 — MRHS solve (split-grid)
Files: `point2all_prop_dumper_claude.cc`
- `q00` mode: outer loop over **4 source colours** $c'$. For each, gather the pre-rotated
  fermion source of all `n_src` sources -> `std::vector<LatticeFermion> f(n_src,FGrid)`
  (`ImportPhysicalFermionSource`), `Grid_split(f,s_f)`, solve on the split grid
  (`SchurRedBlackDiagMooeeSolve` + `ConjugateGradient` tol 1e-9, zero guess, `sD`),
  `Grid_unsplit`, `ExportPhysicalFermionSolution`. Form
  $q_{00}[c,c'](x)=\tfrac12(x_{\text{spin }0,c}+x_{\text{spin }2,c})$ into
  `std::vector<LatticeColourMatrix> q00(n_src, UGrid)`. **4 split-solves.**
- `full` mode: 16-component loop (s,c) -> `std::vector<LatticePropagator> Prop(n_src,UGrid)`.
  **16 split-solves.**
- Log per-component norms as a sanity trace.

### Chunk 4+5 — single-precision demote + SciDAC write (DONE, full-T v1)
Files: `point2all_prop_dumper_claude.cc`
- **v1 decision: full time extent, no plateau window.** In `q00` mode the files are small
  (~0.68 GB/config at 24c, ~2.1 GB at 32c), and `localCopyRegion` asserts identical MPI
  processor layout on every dim and does a purely local copy (`Lattice_transfer.h:764,786`),
  so it CANNOT repartition a time window across the time decomposition. Correct windowing
  would need a communicating `ExtractSlice`/`InsertSlice` per-slice loop — deferred; select
  the plateau downstream at read time. (tmin/tmax knobs dropped from v1.)
- Build single-precision grid `UGridF` (same latt/mpi, `vComplexF` simd).
- `PropRecord : Serializable` metadata {config, srcCoord[4], srcIndex, mass, M5, b, c, Ls,
  smearWidth, smearNiter, store="q00", prec="single"}.
- `ScidacWriter WR(UGridF->IsBoss()); WR.open(outfile)`; per source: `precisionChange(q00F,
  q00[i])` then `WR.writeScidacFieldRecord(q00F, rec)`; `WR.close()`. Field =
  `LatticeColourMatrixF`.

### Chunk 6 (optional) — readback verification
Files: `point2all_prop_dumper_claude.cc` (guarded by a `--verify` flag) or a tiny separate reader.
- `ScidacReader` re-reads records, checks `norm2` matches and `userRecord` round-trips.

## Variational / GEVP framing (Q1 resolved)

The 8 corner sources are a **variational basis for a two-baryon GEVP** to remove
excited-state contamination and resolve the low-lying finite-volume two-baryon spectrum
(the multiple energy levels feeding a Lüscher analysis).

Build a correlator matrix from two-baryon interpolators $O_i$ that differ in the
**relative displacement** $d_i$ of the two baryon source points:
$$
C_{ij}(t) = \langle\, O_i(t)\, O_j^\dagger(0)\,\rangle ,
\qquad
O_j^\dagger(0) \sim \bar B(y_1^{(j)})\,\bar B(y_2^{(j)}),\quad d_j = y_2^{(j)}-y_1^{(j)} .
$$
Solve the GEVP $C(t)\,v_n = \lambda_n(t,t_0)\,C(t_0)\,v_n$; the eigenvalues give
$\lambda_n(t,t_0)\to e^{-E_n(t-t_0)}$, projecting onto the $n$-th finite-volume eigenstate.

### Unmatched (asymmetric) GEVP — source basis need not equal sink basis

(User note `scattering/unmatched_gevp.pdf`.) Let the **sink** operator set be $\{O_i\}$ and the
**source** set $\{\tilde O_\kappa\}$, possibly different (different displacements, smearings,
momenta, even different cardinality $\ge N$). With $g_{in}=\langle0|O_i|n\rangle$,
$\tilde g_{\kappa n}=\langle0|\tilde O_\kappa|n\rangle$,
$$
C_{i\kappa}(t)=\langle O_i(t)\,\tilde O_\kappa(0)\rangle
=\sum_n g_{in}\,e^{-E_n t}\,\tilde g^*_{\kappa n}
\;\Longrightarrow\; C(t)=g\,\Lambda(t)\,\tilde g^\dagger,\quad \Lambda(t)=\mathrm{diag}(e^{-E_n t}).
$$
The source factor cancels in the GEVP:
$$
C(t)\,C(t_0)^{-1}=g\,\Lambda(t-t_0)\,g^{-1},
\qquad \lambda_n(t,t_0)=e^{-E_n(t-t_0)} .
$$
So $C(t)v=\lambda\,C(t_0)v$ yields the correct energies with **any** invertible source basis;
only requirement is $g,\tilde g$ full-rank over the target $N$-state subspace.

Reference for the standard (symmetric) GEVP and its $O(e^{-\Delta E_{N}(t-t_0)})$ excited-state
bound: Blossier, Della Morte, Garcia i Tormo, Sommer, Weber, JHEP 0904:094 (arXiv:0902.1265).

**Design consequence (important for how much we dump).** The sink side is *free*
post-processing on the stored point-to-all propagators (momentum projection, sink smearing,
sink displacement — as many operators as we like). So make the **sink basis rich** downstream
and keep the **expensive source basis** only as large as the number of two-baryon levels we
want to resolve:
- source-basis rank $N$ = number of GEVP levels targeted;
- 8 corners give $\sim3$–$4$ independent source displacement classes $\Rightarrow \sim3$–$4$
  levels; add Gaussian-source widths (Q5) to grow $N$ cheaply on the source side;
- we are NOT obliged to match sink operators to sources, so the number of solves is set purely
  by the desired source rank, not by the sink operator count.

Caveat: unmatched $C$ is non-symmetric/non-normal — eigenvalues still give $E_n$, but the
rigorous symmetric-GEVP contamination bound no longer strictly applies; keep source/sink bases
well-conditioned and check plateau stability in $t_0$.

**Corner displacement classes.** From the 8 corners $\{0,L_\mu/2\}^3$, the source-pair
separations fall into 3 cubic-symmetry classes (from a fixed reference corner):

| class | coords differing | separation | multiplicity from a corner |
|-------|------------------|-----------|----------------------------|
| edge          | 1 | $L/2$              | 3 |
| face diagonal | 2 | $L/\sqrt2$         | 3 |
| body diagonal | 3 | $\sqrt3\,L/2$ (antipodal $O/M$) | 1 |

So the variational dimension is $\sim3$–$4$ (by separation class); the 8 corners let us
**average over symmetry-equivalent pairs** to cut noise. All 8 point-to-all propagators are
the ingredients for every corner pair (sink separation is free — sink is all-points), so
**keep nrhs = 8**.

Implications for the dumper (mostly downstream, but noted):
- Store the **corner index** in each record's `userRecord` so the contraction code can index
  the GEVP matrix and average within a class.
- **Source/sink asymmetry.** Sources are fixed points; the sink side is summed over position
  (momentum projection). So $C_{ij}$ is not manifestly Hermitian — downstream symmetrise
  $C\to(C+C^\dagger)/2$ (standard) or use the non-Hermitian GEVP. Not a dumper concern.
- **Richer knobs later.** Displacement-only is a modest variational basis; the stronger knobs
  for a two-hadron GEVP are usually **relative momentum** (several $|p|$) and **source
  smearing** (radii). Momentum projection is a sink/contraction operation (no new solves);
  smearing needs new sources (Q5). v1 = 8 point corners; widen later.

### Chunk 6 — rotation self-test (`--rottest`, DONE)
Files: `point2all_prop_dumper_claude.cc`
- Validates the pre-rotation / 4-solve trick: builds the reference `q00` from a full
  16-component Weyl propagator with the SAME (smeared) source, rotates
  `q00_ref = peekSpin(0.5 U S Udagger, 0,0)`, compares to the 4-solve `q00`, asserts
  `max norm2(diff) < 1e-8`. Gauge-independent identity; run on the cold (free) field first,
  ideally with `--width 0` (point) for the cleanest check, then smeared.

## Open questions / ideas to explore

1. **Source count vs averaging.** Variational dimension is ~3 (separation classes), but 8
   corners buy symmetry-averaging. Could drop to 4 (origin + one representative per class) if
   storage bites, losing the averaging. Keep 8 for v1.
2. **Plateau window choice.** What $[t_\text{min},t_\text{max}]$? Depends on the effective-mass
   plateau seen in the earlier 5-config two-baryon run ($M_B\sim3.0$, $E_{2B}\sim6.4$).
   Keep it a run knob; pick per ensemble after a first look.
3. **`nrhs` / batch width tuning.** nrhs=8 (=sources) is the default; is a wider batch
   (fold spin-colour in too) worth it on MI300A, or does memory / util favour 8? Benchmark.
4. **MPISPLIT shape.** For 32c, `2.2.2.1` (cubic-ish local $16^3\times64$) vs `1.1.1.8`
   (local $32^3\times8$) — pick by comms/util. Tunable.
5. **Gaussian (smeared) sources — planned variational knob.** Use Grid's gauge-covariant
   Gaussian/Wuppertal smearing `CovariantSmearing<Gimpl>::GaussianSmear(U, chi, width,
   Iterations, orthog=Tdir)` (`Grid/qcd/utils/CovariantSmearing.h:41`), profile
   $\sim e^{-x^2/w^2}$, spatial-only (orthog = time). A few widths $\{w_k\}$ per corner add
   GEVP operators (point = $w{\to}0$); smear the sink identically for a clean $C_{ij}$.
   Each (corner, width) is a new RHS → `nrhs = 8 n_w` (retune `MPISPLIT`) or loop widths
   outer. v1 = point sources; add a `--widths` list + smeared-source step as v1.1. Open:
   thin vs APE/Stout spatial links for the covariant Laplacian; width values.
6. **Sink projection.** Currently full point-to-all sink stored. If storage bites, could we
   store only PD-rotated $q_{00}$ colour matrices in the window (16x smaller) at the cost of
   locking the operator to "0000"? Deferred; full propagator kept for flexibility.

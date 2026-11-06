# Grid (GPU) glueball 1-pt functions for SU(4), validated against HiRep -- implementation plan

Status 2026-11-06: user go-ahead, Q1-Q7 RESOLVED (all proposals accepted); Chunks 8.0 + 8.1 DONE; Chunk 8.2 DONE + VALIDATED vs HiRep (16c, 4 cfg, max rel 1.6e-12); 24c timing moved to tuolumne (user): `submit_glueball_tuo_claude.sh` written, benchmark = 4 m0.4 configs, 1 per GPU (not yet run). Parent plan (all physics decisions, HiRep
side): `/mnt/baracuda_14/hirep_su4/QuenchedSU4_claude/quenched_su4_glueball_impl_plan_claude.md`, Chunk 8.

## Algorithm sources (cite in the code too)

- Glueball operators: Wilson loops projected onto $O_h$ irreps, B. Berg and A. Billoire, Nucl. Phys. B221 (1983) 109;
  R. C. Johnson, Phys. Lett. B114 (1982) 147; moving-frame little groups $D_4$, $D_2$, $D_3$; operator tables generated
  by A. Rago's `OpProjection.m` (patched copy `hirep_su4/GridBridge_claude/glueball_gen_claude/OpProjection_claude.m`).
- APE smearing: M. Albanese et al. (APE), Phys. Lett. B192 (1987) 163.
- Spatial blocking: M. Teper, Phys. Lett. B183 (1987) 345 (HiRep variant: no decimation, step doubles per level).
- SU(N) projection as HiRep `covariant_project_to_suNg` (`hirep_su4/LibHR/Utils/suN_utils.c:334`): unitary polar factor
  (N. J. Higham, Functions of Matrices, SIAM 2008, ch. 8); eigendecomposition by the cyclic Jacobi method, W. H. Press
  et al., Numerical Recipes, 2nd ed. (CUP 1992), sec. 11.1, as HiRep `jacobi2` (`hirep_su4/LibHR/Inverters/jacobi.c`).
- Analysis (unchanged, Python): GEVP, M. Luscher and U. Wolff, Nucl. Phys. B339 (1990) 222.

## Goal

HiRep's glueball measurement is CPU only (~620 s per $24^3 \times 48$ config on 4 cores). Write a Grid code that computes
the SAME 1292 glueball 1-pt functions (all levels 0-3, APE depth 4, 14 momenta) on a GPU, from Grid NERSC configs, so
that the 6 dynamical 24c ensembles (m0.3 + m0.4 first) can be measured on tuolumne MI300A. The output must agree with
HiRep's measure-only mode (`suN_glueballs_h5_claude`, `configlist` key) on the same configs to round-off, so the
existing analysis (`hirep_su4/QuenchedSU4_claude/analyze_glueball_h5_claude.py`) reads either output unchanged.

## What HiRep computes (the spec the Grid code reproduces)

Per config, for each timeslice $t$ (all spatial operations only; time links unused):

1. Level 0 = raw spatial links $U_k(x)$, step $s_0 = 1$.
2. Level $\ell \to \ell + 1$ (`LibHR/Observables/glueballs.c:27`): $n = 4$ APE steps on the level-$\ell$ geometry, then one
   blocking step; $s_{\ell+1} = 2 s_\ell$, lattice size unchanged.
   - APE step (`LibHR/Utils/spatial_transformations.c:284`), $\alpha = 0.33$, staples over the 2 other spatial
     directions $j \neq k$, all neighbours at distance $s_\ell$:

$$
V_k(x) = (1-\alpha)\, U_k(x) + \frac{\alpha}{4} \sum_{j \neq k} \Big[ U_j(x)\, U_k(x+s\hat j)\, U_j^\dagger(x+s\hat k)
  + U_j^\dagger(x-s\hat j)\, U_k(x-s\hat j)\, U_j(x-s\hat j+s\hat k) \Big]
$$

  Followed by the projection $V \to P(V)$ below. Order: the step reads only the old field (Jacobi-style update).
   - Projection (`suN_utils.c:334`): with the principal branch $z^{-1/N} = |z|^{-1/N} e^{-i \arg(z)/N}$, $\arg \in (-\pi, \pi]$,

$$
W = V \det(V)^{-1/N}, \qquad P(V) = W\, (W^\dagger W)^{-1/2}
$$

  $(W^\dagger W)^{-1/2}$ via the Jacobi eigendecomposition $W^\dagger W = E\, \mathrm{diag}(\lambda)\, E^\dagger$. The polar
  factor is unique, so any exact method agrees to round-off; $\det P(V) = 1$ exactly in exact arithmetic.
   - Blocking (`spatial_transformations.c:111`): $U'_k(x) = U_k(x)\, U_k(x + s\hat k)$ at EVERY site.
3. At each level: 241 closed spatial paths (`glueballs_op.c`, `path0..path240`; lengths 4/6/8/10 links: 3/22/120/96 paths),
   step $s_\ell$ per link, trace $\mathrm{tr}_i(x)$ with the path starting at $x$.
4. Momentum projection for the 14 momenta $\vec p \in \{(-1,1,1),(0,0,0),(0,0,1),\dots,(1,1,1)\}$ (order of
   `eval_all_glueball_ops`), with $L = $ `GLB_X`:

$$
C^+_i(\vec p, t) = \sum_{\vec x} e^{2\pi i\, \vec p \cdot \vec x / L}\, \mathrm{Re\,tr}_i(\vec x, t), \qquad
C^-_i(\vec p, t) = i \sum_{\vec x} e^{2\pi i\, \vec p \cdot \vec x / L}\, \mathrm{Im\,tr}_i(\vec x, t)
$$

5. Operator $a$ (1292 of them, index = call order in `eval_all_glueball_ops`) $= \sum$ (complex coefficient) $\times$
   (shift phase) $\times\, C^\pm_i(\vec p_a, t)$, the shift phase being a product of $c_m = e^{i \pi m\, s_\ell / L}$,
   $m \in \{\pm 2, \pm 4, \pm 6, 8\}$ (`c0..c6` in `glueballs_op.c`).
6. Output per config: complex array [level 4][op 1292][t 48] (HiRep h5 `/glue_ape4` row = (4, 1292, T, 2)).

## Design (Grid)

- Fields: `LatticeGaugeField` read by `NerscIO::readConfiguration` (checksum + plaquette checked by Grid); three spatial
  `LatticeColourMatrix` per level; everything 4D (all timeslices at once), double precision.
- Path traces: for path steps $(d_1, \pm), \dots, (d_n, \pm)$, accumulate right to left with one `Cshift` by $\pm s_\ell$ and one
  matmul per link: forward $W \leftarrow U_d\, \mathrm{Cshift}(W, d, +s)$; backward
  $W \leftarrow \mathrm{Cshift}(U_d^\dagger W, d, -s)$. Same trace as HiRep's left-to-right product up to round-off.
  ~2064 link steps per level, ~8300 per config (+ APE: 4 x 3 levels x 3 dirs x 8 products).
- Momentum projection: one `sliceSum` per path over a 28-component object (14 momenta x {Re, Im}) built from
  precomputed phase fields, instead of 28 separate `sliceSum`s.
- Projection kernel: `accelerator_for` over sites, `coalescedRead` -> per-site 4x4 complex matrix, determinant, principal
  root, 4x4 Hermitian Jacobi (port of `jacobi2`, fixed sweep cap + convergence test), `coalescedWrite`. GPU only
  (CUDA locally, HIP on tuolumne); a CPU SIMD build is not targeted.
- Operator table: generated from the INSTALLED `hirep_su4/LibHR/Observables/glueballs_op.c` (single source of truth, the
  file HiRep compiles) by a Python parser: paths as step lists (the loop closure is checked), operators as sparse lists
  of (C+/C-, path, complex coefficient, shift-phase exponent $m_\mathrm{tot}$). The coefficient expressions are
  evaluated symbolically (numbers, `I`, `c0..c6`) so no hand transcription.
- Config loop + restart as `wflow_q_h5_claude.cc` (`--dir --prefix --start --end --step --part --nparts`, skip existing
  output, write to `.tmp` then rename).

## Files

New (all in `Grid_sdm_build/src/gauge_gen_Nc4/` unless the location question says otherwise):
- `glueball_grid_impl_plan_claude.md` (this file)
- `glueball_ops_table_gen_claude.py` (parser of `glueballs_op.c` -> table)
- `glueball_ops_table_claude.h` (generated table; or a text table, open question Q4)
- `glueball_ref_numpy_claude.py` (independent NumPy reference on 16c, Chunk 8.1)
- `glueball_h5_claude.cc` (the Grid measurement code)
- `glueball_h5_aggregate_claude.py` (per-config files -> one HiRep-layout h5, if Q3 = per-config)
- local handoff `tmp_glueball_claude.sh` + log; tuolumne flux script `run_glueball_flux_claude.sh`
Modified: `hirep_su4/QuenchedSU4_claude/quenched_su4_glueball_impl_plan_claude.md` (Chunk 8 pointer + results),
`Makefile_claude` untouched (its `%: %.cc` rule already builds any `.cc` here on tuolumne).

## Chunks

### Chunk 8.0 -- HiRep reference on 16c (user runs, already written)
`hirep_su4/QuenchedSU4_claude/tmp_nersc_test_claude.sh`: HiRep measure-only on `grid_claude/16c/ckpoint_lat.{1000..4000}`
($16^3 \times 8$, 2 ranks, minutes) -> `runs_claude/nersc16c_test/glue_nersc16c_claude.h5` (`/glue_ape4`, all levels).
Files: none new.

### Chunk 8.1 -- operator table + NumPy reference
Parser -> table; NumPy re-implementation of steps 1-6 (APE + projection via SVD polar factor, blocking, `np.roll` paths,
phases, ops) on 16c config 1000; compare with the HiRep h5 at all 4 levels. This checks the table (and our reading of
HiRep) before any GPU code, and gives a second reference.
Files: `glueball_ops_table_gen_claude.py`, `glueball_ops_table_claude.h`, `glueball_ref_numpy_claude.py`.

DONE 2026-11-06. Chunk 8.0 (`tmp_nersc_test_claude.sh`, run by the user): ALL NERSC TESTS PASSED, reference
`hirep_su4/QuenchedSU4_claude/runs_claude/nersc16c_test/glue_nersc16c_claude.h5` (`/glue_ape4` (4 cfg, 4 lev, 1292, 8, 2),
traj 1000-4000; HiRep ~240 s per 16c config on 2 ranks). Parser output: 241 paths (lengths 4/6/8/10: 3/22/120/96, all
closed), 14 momenta, 1292 ops, 11692 terms, total shift exponents $m \in \{-6,-4,-2,0,2,4,6,8\}$; header 7248 lines.
NumPy reference (`python3 glueball_ref_numpy_claude.py <nersc> <glueballs_op.c> <h5> 1000`, 2 min) vs HiRep, config 1000,
per-op $\max_t|\Delta| / \max_t|O|$: level 0 7.2e-14, level 1 1.6e-13, level 2 3.9e-14, level 3 1.1e-14 (medians ~3e-15).
=> table, APE/projection/blocking spec, path directions and phases all confirmed.

### Chunk 8.2 -- Grid code + local build
`glueball_h5_claude.cc`: read, levels, paths, projection, ops, h5 output, config loop. Local build against
`grid_claude/build_cuda129` (Nc = 4, CUDA, HDF5), run on 16c via `tmp_glueball_claude.sh` (GPU: only when device 0 is free).
Files: `glueball_h5_claude.cc`, `tmp_glueball_claude.sh`.

CODED 2026-11-06 (go-ahead): `glueball_h5_claude.cc` = class `GlueballMeasure` (buildPhases, apeStep, blockLinks, pathTrace,
momSums, levelOps, measure, static projectField) + device helpers gb_det (pivoted LU), gb_jacobi (port of HiRep jacobi2,
incl. its initial abs_sum parenthesisation), gb_project_site; `--selftest` (projection vs Eigen SVD polar factor on random
8^4 matrices: near-unitary, gaussian, mixed). `glueball_h5_aggregate_claude.py` (per-config -> HiRep layout + `--compare`;
tested on a mock file built from the NumPy result: PASSED, same numbers as 8.1). Handoff `tmp_glueball_claude.sh`
(compile vs build_cuda129, selftest, 16c 1000..4000 twice (2nd must skip), aggregate + compare at 1e-10).
Analysis of Grid output: `analyze_glueball_h5_claude.py ... --gevp-only` (no /tor in Grid files).

RESULT 2026-11-06 (`tmp_glueball_claude.log`): compile fix = `momSums` public (CUDA: no extended lambda in a private
member). Selftest: max $|P - P_\mathrm{svd}|$ 3.1e-15 (near-unitary), 4.0e-12 (pure gaussian), 4.8e-13 (mixed); 0 Jacobi
failures. 16c (16^3 x 8, TITAN V shared with other jobs at ~90% util): 8.7 s/config (level: paths + momentum sums ~2.0 s,
APE x 4 + blocking ~0.12 s, operators ~1 ms) vs HiRep ~240 s on 2 CPU ranks; rerun skipped all. Aggregator fix: Grid
stores scalars as 1-element attribute arrays (`get_scalar`). COMPARE PASSED (tol 1e-10): per (level, op)
$\max_t|\Delta|/\max_t|O|$ worst 1.6e-12 (lev 0/1, ops 299/306), lev 2 <= 2.5e-13, lev 3 <= 4e-14, medians ~2-3e-15;
|plaq diff| <= 5e-15. Aggregate: `glue16c_test_claude/agg_20261106_claude.h5`.

### Chunk 8.3 -- validation vs HiRep (16c, 4 configs)
Per (level, op): $\max_t |O^\mathrm{Grid} - O^\mathrm{HiRep}| / \max_t |O^\mathrm{HiRep}|$ below the tolerance (Q6), plus a
projection unit test (random near-unitary and far-from-unitary 4x4 matrices vs the NumPy SVD polar factor). Timing.
Files: `tmp_glueball_claude.sh` (check phase), `glueball_h5_aggregate_claude.py`.

### Chunk 8.4 -- tuolumne: build, 1-config 24c cross-check, production
Build with `Makefile_claude` (hipcc, gfx942), one m0.4 config on 1 GPU vs HiRep CPU on the same config (HiRep tuolumne
CPU build from Chunk 7, or locally on a copied config), then m0.3 + m0.4 (traj >= 1000, every 20: 350 + 450 configs),
then the other 4 ensembles; output `/p/lustre5/matsumoto5/obs_glue_2448_<ens>/` (user OK'd).
Files: `submit_glueball_tuo_claude.sh` (renamed from the planned run_glueball_flux_claude.sh; modeled on
Grid/scripts/submit_wflow_q_h5_tuolumne_claude.sh: NPARTS = 4 single-GPU processes per node, --part/--nparts).

2026-11-06: user moved the 24c timing to tuolumne. Benchmark: m0.4 configs 1000..1060 step 20 (4 configs, one per GPU),
`flux batch -q pdebug -t 30m`, OUTDIR `/p/lustre5/matsumoto5/obs_glue_2448_b11p045_m0p4000/bench_claude`. Build:
`make -C .../su4_32c/Grid_sdm_build/src/gauge_gen_Nc4 -f Makefile_claude glueball_h5_claude` (after `source env.sh`).

### Chunk 8.5 (optional) -- torelons
288 torelon ops (`torellons_op.c`: 1 APE step, no blocking, Polyakov lines in spatial directions via `polyleg`, partial
momentum sums `pf`). Different structure from the glueball paths, separate parser.
Files: TBD.

## Cost estimate

Per level ~2064 (Cshift + SU(4) matmul) on the full lattice; $24^3 \times 48$ = 663552 sites, one SU(4) matmul pass
moves ~0.5 GB -> ~1 ms on MI300A (~5 TB/s), Cshift similar -> ~8300 x 2 ms ~ 15-20 s per config incl. sliceSums,
APE and projection (to be measured in 8.3/8.4). 800 configs (m0.3 + m0.4) ~ 4 GPU-hours. HiRep CPU: ~620 s x 4 cores
per config.

## Open questions -- ALL RESOLVED 2026-11-06 (user accepted every proposal: Q1 location gauge_gen_Nc4, Q2 glueballs only in v1, Q3 per-config h5 + aggregator, Q4 generated header, Q5 Jacobi port, Q6 1e-10, Q7 keep NumPy reference)

- **Q1 (location).** Proposal: `Grid_sdm_build/src/gauge_gen_Nc4/glueball_h5_claude.cc`, next to `wflow_q_h5_claude.cc`
  (same Makefile on tuolumne, git repo Grid_sdm_build).
- **Q2 (torelons).** v1 = glueballs only, torelons as optional Chunk 8.5 (no torelon signal in the quenched data)?
- **Q3 (output).** (a) one h5 per config (`glue.<conf>.h5`, as the wflow code; robust to kills, several processes per
  node) + aggregator to the HiRep single-file layout; or (b) one appended h5 per ensemble written by the Grid code
  (as HiRep). Proposal (a).
- **Q4 (table format).** (a) generated C++ header (compiled in, nothing to ship at run time); or (b) a text table read at
  run time. Proposal (a).
- **Q5 (projection algorithm).** (a) port HiRep's Jacobi (same algorithm, robust for any input); or (b) Newton polar
  iteration (needs a per-site 4x4 inverse kernel). Proposal (a).
- **Q6 (tolerance).** Proposal: relative difference per (level, op) $\le 10^{-10}$ (double precision, summation order and
  12 APE + projection steps differ from HiRep).
- **Q7 (NumPy reference, Chunk 8.1).** Keep it (cheap, catches table/spec errors before GPU debugging)?

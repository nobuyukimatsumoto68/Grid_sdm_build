# Two-baryon correlator — remote handoff (Chunk 6)

Handoff for an agent on the machine where the real $24^3\times48$ SU(4) SDM gauge
configurations live. Chunks 1-4 are done and validated locally (free field /
cold config). Your job is **Chunk 6**: build (if needed), run on a real
confining configuration, and analyze the two-baryon correlator.

This file lives in `src/gauge_gen_Nc4/` of the `Grid_sdm_build` repo. Paths
below are relative to that repo unless absolute. The companion design docs are
`two_baryon_impl_plan_claude.md` (full plan + derivations) and
`two_baryon_overview_claude.md` (the existing single-baryon code).

## What the code computes

Two-baryon correlator for two spin-2, all-spin-up, positive-parity (`0000`)
baryons at the two most distant spatial points $O=(0,0,0)$ and
$M=(L/2,L/2,L/2)$, source timeslice fixed at $t_{src}=0$, scanning all sink
slices $t$:
$$ C_{2B}(t) = \langle B(O',t)B(M',t)\,\bar B(O,0)\bar B(M,0)\rangle . $$
Each baryon $B=\varepsilon_{abcd}q_0^aq_0^bq_0^cq_0^d$ (PD spin component 0).
With full antisymmetrisation the contraction collapses to a single $8\times8$
determinant of the four propagator colour blocks:
$$ C_{2B}(t) = 24^4\,\det\begin{pmatrix}P_{O'O}&P_{O'M}\\P_{M'O}&P_{M'M}\end{pmatrix},
\quad P_{X'Y}=q00_Y(X'). $$
Source: `src/gauge_gen_Nc4/two_baryon_corr_claude.cc`. Key routines: `detLU`,
`TwoBaryonCorr`/`assemble8x8`, `colourBlockAt`, `TwoBaryonBruteForce`,
`SelfTestTwoBaryon`, `BaryonSingleDensity`, `main`.

Action conventions (hard-coded, matching the publication baryon code):
MobiusFermionD, $L_s=16$, scaled Shamir $b=1.5$, $c=0.5$, anti-periodic
temporal BC `{1,1,1,-1}`, CG tol $10^{-9}$. $M_5$, mass, config from CLI.

## Build

Compile script at repo root: `compile_two_baryon_claude.sh`. It mirrors
`src/gauge_gen_Nc4/Makefile` (pulls `--cxx/--cxxflags/--ldflags/--libs` from
`grid-config`) and writes the binary to `src/gauge_gen_Nc4/bin/`.

```bash
cd <repo>/Grid_sdm_build
GRID=/path/to/grid/build ./compile_two_baryon_claude.sh        # builds two_baryon_corr_claude.cc
```

IMPORTANT build notes (these tripped up the local build — verify on the remote):
- Point `GRID` at the Grid **build tree** that has `grid-config`, `include/`,
  `lib/libGrid.a`. The script auto-detects `grid-config` in either `$GRID/` or
  `$GRID/bin/`. Default `GRID` is `<repo>/../build` — override as needed.
- If the build's `include/` is an incomplete header copy, the script adds
  `-I<Grid source root>` (`GRID_SRC`, default `<repo>/../Grid`) and
  `-I$GRID/Grid` (for the bare `#include "Config.h"` in `GridStd.h`). If the
  remote Grid is a full `make install`, these are harmless.
- Compiler is `nvcc -x cu`. Two nvcc gotchas already handled in the code, keep
  in mind if you edit: do NOT use Eigen `determinant()` (thrust/complex
  `common_type` error) — use the host `detLU`; and use unqualified `abs` (not
  `std::abs`) on `ComplexD`.
- Must be an `--enable-Nc=4` Grid build.

## Run (Chunk 6 — real configuration)

```bash
cd <repo>/Grid_sdm_build/src/gauge_gen_Nc4
export HDF5_USE_FILE_LOCKING=FALSE
./bin/two_baryon_corr_claude <config_file> <M5> <mass> <outfile.h5> \
    --grid 24.24.24.48 --mpi <decomp> --threads <N> --accelerator-threads <N>
```
- Positional args `<config M5 mass outfile>` select a real NERSC config; omit
  them for the cold/free field. For the LSD ensembles $M_5=1.5$ and `mass`
  matches the ensemble tag (e.g. `conf_..._m0p2000` -> `mass=0.2`).
- Template run script: `test/run_two_baryon/run_claude.sh` (edit `CFGPATH`,
  `cfgfilename`, `i`, and uncomment the `PARAMS="${cfg} $M5 $mass $outfile"`
  line; set `PARAMS_GRID` to `--grid 24.24.24.48 --mpi ...`).

GPU MEMORY: a single GPU could **not** hold $24^3\times48$, $L_s=16$
(`cudaMalloc 2.7 GB` OOM locally). Use an MPI/GPU decomposition, e.g.
`mpirun -np 4 ./bin/two_baryon_corr_claude ... --mpi 2.2.1.1` across 4 GPUs
(adjust to the remote node), and/or `--shm`. Pick `--mpi` so the local volume
fits the per-GPU memory.

Output HDF5 keys (via the `MesonFile` serialisable container, complex `data`
vector indexed by $t$):
- `bar_0000_t` — single baryon, momentum-projected (spatial sliceSum),
  weight +1 (identical normalization to `baryons_0000_dirac`). Gives $M_B$.
- `two_baryon_0000_t` — two-baryon, physical $24^4\det Q$, fixed sink points.

## Analysis to do

1. Effective masses $m_{\rm eff}(t)=\log[C(t)/C(t+1)]$ for both keys.
2. From `bar_0000_t` (momentum-projected) read the single-baryon mass $M_B$
   from its plateau.
3. From `two_baryon_0000_t` read the two-baryon energy $E_{2B}$ from its plateau.
4. Interaction at maximal separation: $\Delta E = E_{2B} - 2M_B$. With
   confinement this is now well-posed (unlike free field).

PHYSICS CAVEATS (important for interpretation):
- The two-baryon sink is at **fixed spatial points** $O',M'$ (position space, no
  momentum projection). Its asymptotic decay is the lowest two-baryon energy in
  that fixed-geometry channel; it is not a definite-relative-momentum state. The
  single-baryon `bar_0000_t` is momentum-projected, so $M_B$ is the
  zero-momentum mass. Comparing $E_{2B}$ to $2M_B$ is meaningful for the energy,
  but a clean scattering/binding extraction may eventually want a
  momentum-projected (or smeared) two-baryon operator — flag if needed.
- In the FREE field this whole check is ill-posed (no confinement => no bound
  baryon; the correlators are just free-quark products and the fixed-point
  two-baryon underflows). Free field only validated the contraction algebra
  (done). So run Chunk 6 on a real confining config.
- `two_baryon_0000_t` can get very small at large $t$ (it did underflow in free
  field); on a confining config check it stays above double-precision round-off
  in the fit window, else restrict the plateau range.

## Validation already passed (local, cold 8^3 lattice)

- `SelfTestTwoBaryon` runs at startup every run: $24^4\det Q$ vs brute-force
  $24^4$-term $\varepsilon$-sum agree to rel ~4e-12; decoupled limit factorizes
  to the product of single-baryon correlators to ~1e-12. (Watch this print in
  the log — it should pass on the remote too.)
- Single-baryon per-site identity `density == 24*det(q00)` to ~1e-24.
- `bar_0000_t` matches the published `baryons_0000_dirac` binary bit-for-bit
  (`test/run_baryon/compare.py`, rtol 1e-10).

## Suggested deliverables back to the user

- The output `.h5` per config.
- Effective-mass plot/table for `bar_0000_t` and `two_baryon_0000_t`, the
  fitted $M_B$, $E_{2B}$, and $\Delta E = E_{2B}-2M_B$.
- Note any plateau/underflow issues and the chosen fit windows.

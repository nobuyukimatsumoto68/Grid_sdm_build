# point2all_prop_dumper — remote handoff (tuolumne run script + MRHS tuning)

Handoff for an agent on **tuolumne** (LLNL, MI300A) where the real SU(4) SDM gauge
configs live. The dumper code is **written and validated locally** (cold/free field);
your job is to (1) build it on tuolumne, (2) write the Flux run script for the
$24^3\times48$ ensemble, and (3) tune the split-grid MRHS decomposition. **Do not submit
jobs or delete/overwrite anything — the user submits and manages files.**

Paths below are relative to `src/gauge_gen_Nc4/` of the `Grid_sdm_build` repo unless
absolute. Companion design doc (read it): `point2all_prop_dumper_impl_plan_claude.md`
(full plan, GEVP / unmatched-GEVP framing, the 4-solve derivation).

## What the binary does

`point2all_prop_dumper_claude.cc` dumps **point-to-all quark propagators** for the
two-baryon GEVP program. Sources: the 8 spatial corners $\{0,L_\mu/2\}^3$ at $t_{src}=0$.
For the spin-2 positive-parity "0000" baryon only the Pauli-Dirac source-spin-0 component
is needed, so it **pre-rotates the source** (Weyl spinor $(e_0+e_2)$ per colour) and solves
**4 times per source** (one per source colour) instead of 16. The stored object per source
is the PD spin-0 colour matrix $q_{00}$ (`LatticeColourMatrix`, $N_c\times N_c$; row = sink
colour, col = source colour), 16x smaller than the full propagator.

Solves are batched with **split-grid MRHS** (`Grid_split`/`Grid_unsplit`, as in
`tests/solver/Test_dwf_mrhs_cg.cc`): the sources are solved in `nbatch = n_src/nrhs`
concurrent batches, `nrhs = prod(mpi)/prod(mpi_split)`. Output: one **SciDAC/LIME** file
per config, `n_src` field records (single precision), each tagged with a `PropRecord`
(source coord/index, mass, M5, b, c, Ls, smear width/niter, store, prec).

Action (hard-coded): `MobiusFermionD`, $L_s=16$, scaled Shamir $b=1.5$, $c=0.5$,
anti-periodic temporal BC `{1,1,1,-1}`, CG tol $10^{-9}$, Schur red-black.

### Local validation already passed
- Compiles clean against `../build` (`compile_two_baryon_claude.sh point2all_prop_dumper_claude.cc`).
- `--grid 8.8.8.8 --mpi 1.1.1.1 --rottest --width 0` (cold config, nrhs=1, nbatch=8):
  rotation self-test PASSED (4-solve $q_{00}$ == 16-solve reference to $<10^{-8}$), wrote
  `ColdConfig.prop.lime` with 8 records + correct metadata (verified `srcIndex` 0-7,
  `store=q00`).

## CLI reference

```
point2all_prop_dumper_claude <config> <M5> <mass> <outfile> \
    --grid Lx.Ly.Lz.Lt --mpi mx.my.mz.mt \
    --split sx.sy.sz.st           # via 4 argv ints: --split sx sy sz st
    --width <w> --niter <N> \
    [--rottest]                   # cold-config self-test only; omit in production
    [Grid runtime opts: --threads, --accelerator-threads, --shm, ...]
```
- Positional `argv[1..4]`: NERSC `config`, `M5`, `mass`, `outfile` (`.prop.lime`).
- `--split` is read as **4 separate integer args** (`--split 1 1 1 4`), matching
  `Test_dwf_mrhs_cg`. `nrhs = prod(mpi)/prod(split)` must **divide** `n_src = 8`.
- Cold config default (no positional args) uses M5=1.5, mass=0.1.

## Build on tuolumne

Same pattern as the two-baryon code — build **on the cluster** against the cluster Grid
build (the local build/include tree is incomplete; the compile script adds the Grid source
tree include as well):
```bash
cd ${ROOT}/Grid_sdm_build          # ROOT=/usr/workspace/lsd/matsumoto5/su4_32c
GRID=${ROOT}/build ./compile_two_baryon_claude.sh point2all_prop_dumper_claude.cc
# -> src/gauge_gen_Nc4/bin/point2all_prop_dumper_claude
```

## Run configuration — 24c ($24^3\times48$)

From the user's Flux allocation (`-N 8 -n 32 -g 1`, 1 GPU/task = 32 GCDs):

| item | value |
|---|---|
| `#FLUX` | `-N 8  -n 32  -g 1  --exclusive  -q pbatch  -t <...>m` |
| `LATT`  | `24.24.24.48` |
| `MPIGRID` | `2.2.2.4` (32 ranks) |
| `--split` | `1 1 1 4` -> prod=4, `nrhs = 32/4 = 8`, `nbatch = 1` (all 8 sources concurrent) |
| local vol / GPU (split) | $24^3\times12$ (32x the full-layout local vol) |
| `M5`, `b`, `c`, `Ls` | 1.5, 1.5, 0.5, 16 (hard-coded) |
| per-config output | ~0.68 GB (8 single-prec $q_{00}$ records) |

Mirror the existing `test/run_two_baryon/run_two_baryon_flux_claude.sh` for env, paths, and
the `flux run` line. Key differences: positional args are the same
`<config> <M5> <mass> <outfile>`, but add `--split 1 1 1 4 --width <w> --niter <N>` and use
`-N 8 --tasks-per-node=4` (32 tasks over 8 nodes). Output file per config:
`${OUTDIR}/<ens>.<conf>.prop.lime`, with a matching `.log`.

### Task: write `test/run_prop_dumper/run_prop_dumper_flux_claude.sh`
- Same `#FLUX` header scaled to `-N 8 -n 32`.
- `source ${ROOT}/env.sh`; `APP=${ROOT}/Grid_sdm_build/src/gauge_gen_Nc4/bin/point2all_prop_dumper_claude`.
- Ensemble knobs (`betastr`, `massstr`, `CFGPATH`, `mass`) — pick the target $24^3\times48$
  ensemble for the scattering study (template values are in the two-baryon script;
  confirm the intended beta/mass with the user).
- `CLIST` of decorrelated config indices (env-overridable).
- `OUTDIR` on lustre (env-overridable); `mkdir -p`.
- Loop configs: skip-if-missing (do NOT overwrite; if a `.prop.lime` exists, skip or warn —
  never `rm`), `flux run -N 8 --tasks-per-node=4 ... "${APP}" "${cfg}" "${M5}" "${mass}"
  "${outfile}" --grid ${LATT} --mpi ${MPIGRID} --split 1 1 1 4 --width ${W} --niter ${N}
  ${PARAMS_GRID}`, tee to per-config log.
- `PARAMS_GRID`/`OPTIONS` as in the two-baryon script
  (`--threads 8 --accelerator-threads 8 --comms-overlap --shm 2048 --shm-mpi 1` etc.).
- **Do NOT** put `--rottest` in production (it adds the 16-solve reference per config).
- **Do NOT** submit; leave a `# Submit yourself: flux batch run_prop_dumper_flux_claude.sh`
  note. **No `rm`/overwrite/kill** anywhere in the script.

## Smearing width / niter (set from single-baryon data)

`width` $w$ sets the source radius (rms $\approx1.22\,w$ in 3D). Constraints:
- **Non-overlap** ceiling: $w \lesssim L/4$ (= 6 at 24c) so the two baryon sources don't
  overlap at their minimum separation $L/2$.
- **Stability/convergence**: $N \gtrsim 3w^2$ (fewer diverges/under-converges).
Pick $w$ from the single-baryon plateau/overlap optimization (the user will point you at the
existing single-baryon data); then set $N\gtrsim3w^2$. Placeholder defaults in the code are
$w=3.0$, $N=40$. Smearing uses **thin links** (covariant Gaussian, spatial-only, gauge
invariant — no gauge fixing).

## Tuning tasks (deferred here on purpose)

1. **`nrhs` / batch width (plan Q3).** Default `--split 1 1 1 4` gives nrhs=8 (all sources at
   once, nbatch=1). Benchmark MI300A throughput vs memory: is nrhs=8 the sweet spot, or does
   a narrower split (e.g. `1 1 1 2` -> prod=2, nrhs=16 — **but** nrhs must divide 8, so 16 is
   invalid; use nrhs in {1,2,4,8}) fill the GPU better? Note nrhs>8 is not usable with only 8
   sources unless the code is extended to fold source-colour into the batch — do NOT do that
   here; just tune among nrhs $\in\{1,2,4,8\}$.
2. **`mpi_split` shape (plan Q4).** For nrhs=8 the split must satisfy prod(split)=4. Compare
   `1 1 1 4` (local $24^3\times12$) vs `2 2 1 1` (local $12^2\times24\times48$) for
   comms/utilization. Time a single config each way; keep the faster.
Report timings; hand the chosen `--split`/nrhs back to the user.

## Output + verification

Per config, one `.prop.lime` with 8 records. Sanity checks:
- `grep -a -c ildg-binary-data <file>` == 8; `grep -a -o '<srcIndex>[0-7]</srcIndex>'`
  shows 0..7; `<store>q00</store>` x8.
- Optional correctness: one **cold-config** run with `--rottest --width 0` must print
  `rotation self-test PASSED` (this is gauge-independent, so it validates the build).
- `q00[i] norm2` printed per source should be finite and O(1)-ish (cold) / config-dependent.
- rsync the `.prop.lime` (+ logs) back for the downstream contraction, following the
  existing `rsync_pull_*` pattern.

## Downstream (context, not your task)

These $q_{00}$ records feed the **unmatched-GEVP** two-baryon contraction: source operators
= corner-displacement two-baryon interpolators (rank sets the number of resolvable levels);
sink operators (momentum projection, sink displacement/smearing) are built in
post-processing on the stored $q_{00}$ (all commute with the fixed spin structure). Source
and sink bases need not match. See the impl plan's "Unmatched (asymmetric) GEVP" section.

## Repo conventions (follow these)
- Any file you create/edit gets `_claude` before the extension.
- One statement per line (all languages, including shell). No Unicode in comments/text —
  use LaTeX macros (`\pi`, `-`, ...).
- Never `rm`/overwrite/kill in any script; if a rerun is blocked by an existing output,
  flag it and let the user remove it. Never submit the Flux job yourself.

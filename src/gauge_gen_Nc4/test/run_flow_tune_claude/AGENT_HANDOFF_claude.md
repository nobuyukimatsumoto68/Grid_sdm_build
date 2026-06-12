# Agent handoff: $L_5$ / $M_5$ tuning via eigenvalue flow

You are picking up a lattice-QCD tuning task on a GPU cluster (tioga or lassen, flux
scheduler). This document is self-contained: it tells you the goal, where everything is,
what to run, what to decide, and how to know you are done. Read it fully before acting.

## 1. Objective

Tune the domain-wall / Mobius parameters for an SU(4) stealth-dark-matter ensemble:

- **$M_5$** (domain-wall height): tuned by the *eigenvalue flow* of the Hermitian Wilson
  operator $H_W(M_5)$ as a function of $M_5$. Pick $M_5$ where the low spectrum is cleanest
  (largest, most stable smallest eigenvalue), which minimizes residual chiral breaking.
- **$L_5 = L_s$** (fifth-dimension extent): tuned later by driving the residual mass
  $m_\text{res}$ down. That step uses a different binary (`Mobius_mesons_xt`) and is OUT OF
  SCOPE here unless explicitly asked. This task delivers the $M_5$ flow.

Target physics point:

| quantity        | value                          |
|-----------------|--------------------------------|
| volume          | $32^3 \times 64$               |
| quark mass $m$  | $0.01$                         |
| $L_s$           | $16$                           |
| $M_5$ (ballpark)| $1.5$ (scan around it)         |
| Mobius $b, c$   | $1.5,\ 0.5$ (Shamir-like)      |
| gauge $\beta$   | $10.7,\ 10.8,\ 10.9$ (3 streams)|

## 2. What the codes do

- `dweofa_mobius_HSDM_v3.cc` — Mobius EOFA HMC. Reads an XML via `--ParameterFile`,
  generates gauge configs, and the NERSC checkpointer writes `ckpoint_lat.<traj>`
  (format `IEEE64BIG`) plus `ckpoint_rng.<traj>`.
- `eye4_anti.cc` — the eigenvalue-flow code. Builds `WilsonFermionD` with mass $=-M_5$ and
  anti-periodic temporal BC `{1,1,1,-1}`, forms `MdagMLinearOperator` (eigenvalues of
  $H_W^2$), and runs Chebyshev-accelerated Implicitly Restarted Lanczos to get the lowest
  eigenvalues. CLI: `eye4_anti <config> <alpha> <beta> <M5>  [Grid flags]`.
  It **always writes `./evals.dat`** (one line `M5 <eigvals>`) and **truncates it every call**.

## 3. Locations

- Repo root: `GRID_SDM_ROOT` (ask the human / find the checkout of `Grid_sdm_build`).
- Source: `${GRID_SDM_ROOT}/src/gauge_gen_Nc4/` (the `.cc` files, `Makefile`, `bin/`).
- This scratch: `${GRID_SDM_ROOT}/src/gauge_gen_Nc4/test/run_flow_tune_claude/`, containing:
  - `ip_hmc_mobius_claude.xml` — HMC parameter template (`@BETA@` placeholder).
  - `run_hmc_claude.sh` — thermalize one stream per $\beta$.
  - `run_flow_claude.sh` — eigenvalue-flow scan over saved configs.
  - `README_claude.md` — same content in brief; this file is the fuller brief.
- Site env scripts: `${GRID_SDM_ROOT}/tioga/env.sh`, `${GRID_SDM_ROOT}/lassen/env.sh`.

## 4. Prerequisites: build the binaries

Neither binary is built yet. In `src/gauge_gen_Nc4/`:

```bash
source ${GRID_SDM_ROOT}/tioga/env.sh        # or lassen/env.sh, matching the machine
make dweofa_mobius_HSDM_v3                   # -> bin/dweofa_mobius_HSDM_v3
make eye4_anti                              # -> bin/eye4_anti
```

The `Makefile` expects a Grid install at `../../install/Grid_omp_Nc4` (see `grid-config`).
If `grid-config` is missing, Grid must be built/installed first (see repo build scripts,
e.g. `build_grid_ubuntu.sh` / the `tioga`/`lassen` build dirs). Do NOT proceed to runs until
both binaries exist in `bin/`.

## 5. Step 1 - thermalize a few configs (per $\beta$)

Edit the config block at the top of `run_hmc_claude.sh`:

- `GRID_SDM_ROOT`, `ENV_SH`, `APP` (paths for this machine).
- `NODES`, `MPI_GEOM`, `TPN` so that `prod(MPI_GEOM) == NODES * TPN`. Default
  `MPI_GEOM=2.2.2.8` (= 64 ranks) with `TPN=4` -> `NODES=16`; local volume $16^3\times8$.
- flux `-t` wall time.

Submit (runs the three $\beta$ streams sequentially inside one allocation):

```bash
cd ${GRID_SDM_ROOT}/src/gauge_gen_Nc4/test/run_flow_tune_claude
flux batch run_hmc_claude.sh
```

Output: `b10.7/`, `b10.8/`, `b10.9/`, each with `ckpoint_lat.<traj>` and a
`hmc_b<beta>_claude.log`.

**Check thermalization before trusting configs.** Inspect the plaquette history in each
`hmc_*.log`. The scratch XML uses `Trajectories=40`, `saveInterval=5` (a few configs only) -
this is almost certainly too short for real thermalization. Extend `Trajectories` and
continue a stream with `StartingType=CheckpointStart` in the XML if needed. Report the
plaquette trend to the human rather than assuming equilibration.

## 6. Step 2 - eigenvalue flow scan ($M_5$)

Edit the config block at the top of `run_flow_claude.sh`:

- Same path/MPI variables (default `MPI_GEOM=2.2.2.4`, `NODES=8`; the 4D op is lighter).
- `M5_LIST` (default `seq 1.0 0.1 2.0`), `ALPHA`, `BETA_CHEBY`.

Submit:

```bash
flux batch run_flow_claude.sh
```

For every config it loops over `M5_LIST`, runs `eye4_anti` in a per-config subdir
`b<beta>/flow_<cfg>/`, and cat-appends each `evals.dat` into
`b<beta>/flow_<cfg>/flow_b<beta>_<cfg>_claude.dat`. Per-call stdout is in
`eye4_*_M5_*_claude.log`.

### Chebyshev window (the one tuning knob in eye4_anti)

`ALPHA`, `BETA_CHEBY` map to `ChebyParams.alpha/beta` (`Npoly=101`) - the Lanczos
acceleration window. `BETA_CHEBY` should be $\gtrsim \lambda_\max$ of $MdagM$; `ALPHA` just
above the low cluster you want to resolve. In-code reference value was
`Chebyshev(0.0, 10.0, ...)`. If `eye4_anti` reports few converged eigenvalues (`Nconv`
small) or the lowest eigenvalues look wrong, widen `BETA_CHEBY` and/or lower `ALPHA` and
rerun that point. Defaults `ALPHA=0.5`, `BETA_CHEBY=12.0` are a starting guess only.

## 7. Step 3 - interpret and report

Each line of `flow_*.dat` is `M5  <ascending list of smallest eigenvalues of MdagM>`.
Since $MdagM \approx H_W^2$, an eigenvalue near zero means a near-zero mode of $H_W$, which
drives large $m_\text{res}$.

Deliverable to the human:

1. For each $\beta$ and config, the smallest eigenvalue vs $M_5$.
2. A recommended $M_5$ (and a sensible window) where the smallest eigenvalue is largest and
   most stable across configs and across the three $\beta$ values.
3. Any configs/points where Lanczos failed to converge and what window fixed them.

A quick plotting/extraction script (python/numpy) reading the `flow_*.dat` tables is fine to
write; keep it under this scratch dir with a `_claude` suffix.

## 8. Hard rules (must follow)

- **Never put `rm`, `rmdir`, or any destructive delete/overwrite into a script.** If a rerun
  is blocked by a pre-existing output (e.g. an existing `evals.dat`, a checkpoint, a `complete`
  flag), STOP and ask the human to remove it. The flow script intentionally truncates only via
  shell redirection on files it owns; do not add `rm`.
- **Any file you create or edit gets `_claude` before the extension** (e.g. `plot_flow_claude.py`).
- **No Unicode in code/comments/text** - use LaTeX macros (`\beta`, `\lambda`, etc.); ASCII
  hyphen for minus.
- **`prod(MPI_GEOM)` must equal `NODES * tasks-per-node`** or the launch aborts.
- Do not run two `eye4_anti` invocations in the same directory at once (shared `evals.dat`).

## 9. Decisions to confirm with the human before/while running

- Are $\beta = 10.7, 10.8, 10.9$ the intended couplings for $32^3\times64$, $m=0.01$? (This is
  the one physics input that was not derivable.)
- Node count / MPI layout actually available -> set `MPI_GEOM`, `NODES`.
- How many thermalized configs per $\beta$ are wanted, and how long to thermalize.
- Whether to also do the $L_s$ / $m_\text{res}$ side now (`Mobius_mesons_xt.cc`); if so, ask
  for a separate handoff - it is not covered here.

## 10. Definition of done

Both binaries built; at least a few thermalized (plaquette-checked) configs per $\beta$
saved under `b<beta>/`; `flow_*.dat` tables produced across the $M_5$ grid for each; and a
short written recommendation of $M_5$ with the supporting smallest-eigenvalue-vs-$M_5$
evidence handed back to the human.

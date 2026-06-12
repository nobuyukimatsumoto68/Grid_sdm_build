# L5 / M5 tuning scratch

Minimum pipeline to thermalize a few configs and run the eigenvalue-flow code
for tuning $M_5$ (and, with `Mobius_mesons_xt`, $L_5=L_s$ via $m_\text{res}$).

Target point: $32^3\times64$, $m=0.01$, $L_s=16$, $M_5=1.5$ (Mobius $b=1.5,\ c=0.5$),
gauge $\beta \in \{10.7, 10.8, 10.9\}$.

## Files

- `ip_hmc_mobius_claude.xml` — HMC parameter template; `@BETA@` substituted per stream.
- `run_hmc_claude.sh` — thermalizes one stream per $\beta$ into `b<beta>/`, saving NERSC
  configs `b<beta>/ckpoint_lat.<traj>`.
- `run_flow_claude.sh` — for each saved config, scans $M_5$ with `eye4_anti` and collects
  the lowest eigenvalues of $H_W(M_5)^2$ into `b<beta>/flow_<cfg>/flow_*.dat`.

## Build (once, in `src/gauge_gen_Nc4/`)

```bash
source <repo>/tioga/env.sh        # or lassen/env.sh
make dweofa_mobius_HSDM_v3        # HMC -> bin/dweofa_mobius_HSDM_v3
make eye4_anti                    # flow -> bin/eye4_anti
```

## Run

Set `GRID_SDM_ROOT` (and `ENV_SH`/`APP` if non-default) at the top of each script,
then size `NODES` / `MPI_GEOM` / flux `-t` for your allocation. Submit:

```bash
flux batch run_hmc_claude.sh      # 1) thermalize 3 streams
flux batch run_flow_claude.sh     # 2) eigenvalue flow over saved configs
```

## Reading the result

Each line of `flow_*.dat` is `M5  <list of lowest eigenvalues of MdagM>`.
$H_W(M_5)^2$ eigenvalues near zero <=> near-zero modes of $H_W$ <=> large residual
mass. Plot the smallest eigenvalue vs $M_5$ and pick $M_5$ in the region where it is
largest / most stable across configs and across $\beta$. Then fix that $M_5$ and tune
$L_s$ by driving $m_\text{res}$ down with `Mobius_mesons_xt` (separate measurement).

## Gotchas / things to set on-site

- **beta is the only physics input I could not default** — confirm 10.7/10.8/10.9 are
  the intended couplings for this volume/mass.
- **Thermalization length**: `Trajectories=40`, `saveInterval=5` in the XML is a scratch
  value. Check the plaquette history in `b<beta>/hmc_*.log` and extend before trusting
  the configs as "thermalized". Use `StartingType=CheckpointStart` to continue a stream.
- **MPI layout**: `prod(MPI_GEOM)` must equal `NODES * tasks-per-node`. Defaults assume
  4 GPUs/node (tioga/lassen MI300A/V100). Local volume for HMC default is $16^3\times8$.
- **Chebyshev window** (`ALPHA`, `BETA_CHEBY` in `run_flow_claude.sh`): these set the
  Lanczos acceleration in `eye4_anti` (`ChebyParams.alpha/beta`, `Npoly=101`). `BETA_CHEBY`
  should be $\ge \lambda_\max$ of MdagM; `ALPHA` just above the low cluster you want to
  resolve. If the lowest eigenvalues do not converge (`Nconv` small), widen the window.
  The in-code reference value was `Chebyshev(0.0, 10.0, ...)`.
- **`eye4_anti` writes `evals.dat` and truncates it every call** (`eye4_anti.cc:78`), so the
  flow script runs each call in its own `flow_<cfg>/` dir and cat-appends. Do not run two
  `eye4_anti` in the same dir concurrently.
- The current `eye4_anti.cc` takes a single `M5` from `argv[4]`; the internal M5 sweep loop
  is commented out (`eye4_anti.cc:83`). The scan is driven externally by this script. If you
  prefer an in-binary sweep, uncomment that loop instead (note it would append to one
  `evals.dat`).
```

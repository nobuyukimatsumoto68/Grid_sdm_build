# Glueball measurement on 32^3 x 64 -- implementation plan

**Status 2026-11-24: user decisions relayed by the Global organizer [dafddd]: (a) go-ahead for 32c driver copies + one
pdebug validation job (4 configs), user submits; (b) same 4 blocking levels, APE 4 x 0.33 per level as 24c; (c) EVERY
config of each 32c ensemble (same list as the mesons / disc); (f) channels run in parallel; (h) YES (user 2026-11-24 via organizer): rebuild with tracing on
tuolumne, run the 32c jobs with GLUE_TRACE=1.** Parent plan: glueball_grid_impl_plan_claude.md (Chunks 8.0-8.6).

## Goal

Run the validated Grid glueball code (glueball_h5_claude, HiRep-validated to 1e-12; ops: B. Berg and A. Billoire, Nucl.
Phys. B221 (1983) 109; HiRep generator by A. Rago) unchanged on the six 32c streams. The binary reads $L$, $T$ from
the grid (FullDimensions); nothing in it is 24-specific. Only the drivers change.

## Ensembles (tuolumne, /p/lustre5/matsumoto5/32_64_<m>/beta<b>m<m>/ckpoint_lat.<traj>, NERSC format)

| m | beta | dir | configs (every one) |
|---|---|---|---|
| 0.01 | 10.8 | 32_64_0.01/beta10.8m0.01 | 4..400 step 4 (100) |
| 0.05 | 10.84 | 32_64_0.05/beta10.84m0.05 | 4..400 step 4 (100) |
| 0.1 | 10.865 | 32_64_0.1/beta10.865m0.1 | 20..2000 step 20 (100 when the chain is done) |
| 0.2 | 10.99 | 32_64_0.2/beta10.99m0.2 | 20..2000 step 20 (100) |
| 0.3 | 11.035 | 32_64_0.3/beta11.035m0.3 | 20..2000 step 20 (100) |
| 0.4 | 11.045 | 32_64_0.4/beta11.045m0.4 | 20..2000 step 20 (100) |

Betas and dir names as Grid/scripts/run_disc_mixedprec_32c_tuo_claude.sh (which measures the same config list) and
run_prod32c_claude.sh. The driver lists the configs that exist (no thermalization cut at measurement time).
Output: /p/lustre5/matsumoto5/obs_glue_3264_<tag>/glue_h5/glue.<traj>.h5 (tag b10p800_m0p0100 style).

## Cost

24c: 47-49 s/config/GPU. Volume ratio $(32/24)^3 (64/48) = 3.16$ -> ~150 s/config/GPU (estimate), 4 single-GPU
processes per MI300A node. 100 configs per ensemble -> ~1.1 h per 1-node job; NPER = 200 per job (one job per
ensemble, 4 h wall). Total ~6-7 node-hours. Memory per process ~4x 24c, well inside one MI300A.

## Files

- NEW submit_glueball_tuo_32c_claude.sh (copy of submit_glueball_tuo_claude.sh: Lattice 32.32.32.64; GLUE_TRACE passed
  through when set).
- NEW run_glueball_tuo_32c_claude.sh (copy of run_glueball_tuo_claude.sh: 32c stream table, ckpoint_lat prefix,
  obs_glue_3264 output, NPER 200).
- Later: aggregator run (glueball_h5_aggregate_claude.py, unchanged), rsync pull script for obs_glue_3264.

## Chunks

- **C1. DONE 2026-11-24.** Driver copies. Files: the two scripts above (bash -n checked).
- **C2.** pdebug validation (user submits): 4 configs of m0.4 (one per GPU), check 4 glue.<traj>.h5 written, timing per
  config ("# conf" lines in the part logs), no NaN; pull them and aggregate with glueball_h5_aggregate_claude.py
  (expect /glue_ape4 shape (4 cfg, 4, 1292, 64, 2)).
- **C3.** Production: `bash run_glueball_tuo_32c_claude.sh` (6 jobs), user submits; reruns resume (per-config h5 skip).

## (h) GLUE_TRACE rebuild -- steps for the user's decision

The tracing source is already committed (Grid_sdm_build commit 7e740e9, 2026-11-17). On a tuolumne login node:

```
cd /usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build
git pull
source /usr/workspace/lsd/matsumoto5/su4_32c/env.sh
make -C src/gauge_gen_Nc4 -f Makefile_claude glueball_h5_claude 2>&1 | tee src/gauge_gen_Nc4/make_glueball_h5_32c_claude.log
```

Then submit with `GLUE_TRACE=1` in the environment of run_glueball_tuo_32c_claude.sh (passed to the jobs): ~1446 extra
flushed trace lines per config in the per-part logs, measurement unchanged (verified locally 2026-11-17, identical
output). Without the rebuild the old binary ignores GLUE_TRACE.

## Open questions

1. RESOLVED 2026-11-24: (h) rebuild + GLUE_TRACE=1 for pdebug and production.
2. m0.1: start when its HMC chain reaches 2000 (the driver measures what exists; a rerun later picks up the rest).

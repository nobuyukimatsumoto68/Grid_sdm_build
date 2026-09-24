# Baryon operator variants: production handoff (2026-09-24)

End-to-end handoff for the single-baryon correlators of the nine B-form highest-weight operators
($2^+$ (4,0),(2,2),(0,4); $2^-$ (3,1),(1,3); $1^-$ (3,1),(1,3); $1^+$ (2,2); $0^+$ (2,2)) on the four
heavy $24^3\times48$ SU(4) SDM ensembles. Pipeline: **full-PD propagator dump (tuolumne GPU) ->
per-operator contraction (tuolumne GPU, same filesystem) -> pull to local -> python analysis**.
Design and validation: `baryon_variants_impl_plan_claude.md` (same directory). The remote agent
writes and tunes the run scripts; **the user submits all jobs; never rm / overwrite / kill.**

## Ensembles (as in `two_baryon_gevp_production_handoff_claude.md`)

| mass | beta (str) | stride | #cfg | index range   | $2^+$ (4,0) $m_\text{eff}$ (obs) |
|------|-----------|--------|------|---------------|----------------------------------|
| 0.1  | 10p865    | 40     | 92   | 1000 .. 4640  | ~1.56 |
| 0.2  | 10p990    | 40     | 93   | 1000 .. 4680  | |
| 0.3  | 11p035    | 80     | 88   | 1000 .. 7960  | |
| 0.4  | 11p045    | 80     | 113  | 1000 .. 9960  | ~3.0 (validation anchor) |

Gauge configs (lustre5): `/p/lustre5/matsumoto5/conf_nc4nf1_2448/conf_nc4nf1_2448_b<beta>_m<mass>/conf_nc4nf1_2448_b<beta>_m<mass>_lat.<idx>`.
Action: MobiusFermionD, Ls=16, b=1.5, c=0.5, M5=1.5, AP time BC, CG tol 1e-9 (hard-coded, same as the q00 dumper).

**Start with m0.4 / b11p045, lat.1000 alone** (validation: $2^+$ (4,0) plateau ~3.0 must reproduce the
q00-based single baryon; then the $2^-$ against the publication value), then the full m0.4 set, then the others.

## Files to copy to the cluster (`/usr/workspace/lsd/matsumoto5/su4_32c/`, i.e. `su4_32c`)

From `Grid_sdm_build/src/gauge_gen_Nc4/`:
- `baryon_pd_common_claude.h` (shared header; included by both sources)
- `point2all_prop_dumper_full_claude.cc` (stage 1)
- `baryon_variants_corr_claude.cc` (stage 2)
- `baryon_variants_terms_claude.txt`, `baryon_variants_ops_claude.txt`, `baryon_variants_ops_all_claude.txt`
  (text tables read at run time; pass their paths with `--terms` / `--ops` or run in the directory that holds them)
- `baryon_variants_report_claude.py` (optional per-config sanity report; needs h5py)

Build (tuolumne, hipcc, same as the q00 dumper):
```
cd /usr/workspace/lsd/matsumoto5/su4_32c
GRID=${PWD}/build ./compile_two_baryon_claude.sh point2all_prop_dumper_full_claude.cc
GRID=${PWD}/build ./compile_two_baryon_claude.sh baryon_variants_corr_claude.cc
```
(The compile script resolves the source under `src/gauge_gen_Nc4/` and puts binaries in `src/gauge_gen_Nc4/bin/`.)
Portability note: the contraction kernel is plain `accelerator_for` code (validated on CUDA; HIP uses the
same macros). Do not name any kernel variable `nt`: Grid's `accelerator_for` macro declares a local `int nt`
around the lambda (`Grid/threads/Accelerator.h:136,:440`).

## Stage 1 -- full-PD propagator dump (tuolumne, MI300A)

Binary `point2all_prop_dumper_full_claude`. Per config: the 16 PD source columns (4 PD source spins x
4 colours, pre-rotated sources) of the origin corner are solved in 16/nrhs split-grid batches; the four
PD sink blocks of each solution give the 16 blocks `q<s><s'>` (LatticeColourMatrix, row = sink colour,
column = source colour), written as 16 single-precision SciDAC records per corner, tagged
`PropRecord{srcIndex = corner, store = "q<s><s'>", smearWidth, smearNiter, ...}`.

CLI (positional as the q00 dumper, plus knobs):
```
point2all_prop_dumper_full_claude <config> <M5> <mass> <outfile> --grid 24.24.24.48 --mpi 2.2.2.4 \
    --split 1 1 1 4 --width 3.0 --niter 40 --nsrc 1 [--rottest]
```
- `--nsrc n` : first n of the 8 corners {0,L/2}^3 at t=0 (**production: 1 = origin only**; more corners later).
- `--width/--niter` : source smearing, same convention as the q00 dumper (`SRC=smeared` w=3.0 N=40;
  `SRC=point` -> `--width 0`). **Use the same source smearing as the q00 production (smeared w=3 N=40)
  so the $2^+$ (4,0) channel can be checked against the q00 single baryon.**
- `--split sx sy sz st` : split-grid MRHS; nrhs = prod(mpi)/prod(split) **must divide 16**
  (`--mpi 2.2.2.4 --split 1 1 1 4` -> nrhs 8, 2 batches per corner; the q00 production layout).
- `--rottest` : self-test (16 extra unit solves per corner, full-grid); run ONCE on lat.1000, not in production.
- Env `PROP_DEADLINE_EPOCH` / `PROP_TPT_SECONDS`: graceful wall blocker as in the q00 dumper; atomic
  `.inprogress` -> rename write.

Cost per config (estimate from the q00 dumper: 32 solves in 61-78 s at `-N 8 -n 32 -g 1`):
**16 solves, ~35-45 s**, plus the write of **1.36 GB** (16 x 85 MB single-precision records).
Storage: 1.36 GB x 113 configs = **154 GB for m0.4**; all four ensembles = **~525 GB** on lustre5.
Output: `/p/lustre5/matsumoto5/bvar_2448_b<beta>_m<mass>/conf_nc4nf1_2448_b<beta>_m<mass>.<idx>.pdfull.lime`.

Script to write: `test/run_bvar_dumper/run_bvar_dumper_flux_claude.sh` mirroring
`run_prop_dumper_chain_claude.sh` / `submit_prop_dumper_b<beta>_m<mass>_claude.sh`
(`#FLUX -N 8 -n 32 -g 1 -t 480m`, `flux run -N 8 --tasks-per-node=4`, PARAMS_GRID
`--threads 8 --accelerator-threads 8 --comms-overlap --shm 2048 --shm-mpi 1`), knobs ENS / CLIST / SRC,
skip-if-exists on the final `.pdfull.lime`.

## Stage 2 -- contraction (tuolumne, same node type; reads the dump, GPU kernel)

Binary `baryon_variants_corr_claude`. One config per invocation; **one operator block per run
(parallel runs) or all nine in one run** (the contraction itself takes seconds on a GPU; the run is
dominated by reading the 1.36 GB dump, so `--op all` is the economical default and per-op runs are
for parallelism only).
```
baryon_variants_corr_claude --grid 24.24.24.48 --mpi 1.1.1.1 \
    --src <file.pdfull.lime> --out <bvar.<op>.<idx>.h5> \
    --terms baryon_variants_terms_claude.txt --ops baryon_variants_ops_claude.txt \
    --op all|2p_40|2m_31|1m_31|1p_22|0p_22 --pairs block \
    [--config <nersc gauge> --sink-smear 3.0 40]
```
- `--op <name>` computes the whole $J^P$ block that contains `<name>` (all cross-correlators of that
  block); the five block representatives are `2p_40`, `2m_31`, `1m_31`, `1p_22`, `0p_22`.
- `--pairs block|sector|all` : which elementary $C_{AA'}$ to compute and store (`block` = the ones the
  operators need; `sector` = every same-(parity, $J_3$ mod 4) pair, needed later for lower-$M$ states
  and the $E$/$T_2$ split; `all` = all 157). Use `sector` in production if disk is no concern (tiny files).
- `--sink-smear 3.0 40 --config <gauge>` adds the smeared-sink copies (`Celss_*`, `Copss_*`); the smear
  of a full propagator costs 6 shifts of 2.7 GB per iteration, i.e. comparable to the solves: include it
  (the smeared-sink single baryon was what resolved the smeared/point disagreement in the GEVP study).
- Memory: one rank holds the 2.7 GB double propagator (+2.7 GB smeared copy): `--mpi 1.1.1.1` on one
  MI300A is fine.
- `--naive-check <n>` (validation only, lat.1000): compares the kernel against the naive Wick loop at n
  sampled sites; expect `naive-check PASSED`.
Output: `/p/lustre5/matsumoto5/obs_bvar_2448_b<beta>_m<mass>[_ss]/bvar.<op>.<idx>.h5`
(datasets `meta` [HDF5 attributes], `Cel_<A>_<A'>`, `Cop_<O>_<O'>`, `Celss_*`, `Copss_*`; ~50-200 KB).
Script to write: `test/run_bvar_contract/run_bvar_contract_flux_claude.sh` (fan-out over configs,
`HDF5_USE_FILE_LOCKING=FALSE`, skip-if-exists, knobs ENS / CLIST / OP / PAIRS / SINKSMEAR).

## Stage 3 -- pull to local

Repo root: `./rsync_pull_bvar_claude.sh` (oslic, lustre5 by default; `DRYRUN=1` to list) ->
`./obs_bvar_2448/obs_bvar_2448_b<beta>_m<mass>[_ss]/bvar.<op>.<idx>.h5`.

## Stage 4 -- local aggregation + analysis

- `python3 bvar_h5_claude.py` -> `bvar_h5_out_claude/bvar_b<beta>_m<mass>[_ss].h5`
  (attrs: run meta + beta/ens; `confs`; `Cop/<O>_<O'>`, `Cel/<A>_<A'>`, `Copss/`, `Celss/`, shape (Nconf,T)).
- Analysis (to write, notebook, one plot per cell, jackknife as in `analyze_gevp_claude.ipynb`):
  1. $m_\text{eff}(t)$ of every diagonal $C_{OO}$; $2^+$ (4,0) must reproduce the q00 single baryon
     (`gevp_h5_out_claude` `CB` data) exactly up to the factor 24 on the same configs.
  2. $2^\pm$ ground states vs the publication (LSD "thermo gw details" Sec. 1.2 / Fig. 4.6 values in
     flow units; the lattice-unit reference is the output of `baryons_S2_dirac_parity.cc` on these ensembles).
  3. Level coincidence: same plateau in (4,0)/(2,2)/(0,4) for $2^+$, in (3,1)/(1,3) for $2^-$ and $1^-$;
     small GEVP per block (3x3, 2x2) with $t_0$ stability.
  4. Backward propagation: parity partners appear in the backward half of every channel.

## Sanity checks per config (remote agent)

- dump: `grep -a -c ildg-binary-data <file>` = 16 x nsrc; tags `<store>q00</store>` .. `<store>q33</store>` each once per corner.
- contraction log: `PDTermTable: read 157 pairs`, `elementary pairs to compute: 42 (needed by the blocks: 42)`
  for `--pairs block --op all`; `wrote ... 42 Cel + 19 Cop datasets`.
- lat.1000 only: `--naive-check 256` -> `naive-check PASSED`; `Cop_2p_40_2p_40` = 24 x `CB_set0_corner0`
  of the existing `obs_gevp_2448_b11p045_m0p4000/gevp.1000.h5` (same source smearing, point sink).

## Conventions

New/edited files carry `_claude` before the extension; one statement per line; no Unicode in comments
(LaTeX macros); never rm/overwrite/kill in scripts; the user submits all jobs.

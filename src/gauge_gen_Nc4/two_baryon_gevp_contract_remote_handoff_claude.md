# two_baryon_gevp_contract — remote handoff (run over the 386 configs on tuolumne)

Handoff for an agent on **tuolumne** (LLNL, MI300A). The contraction code is **written and
validated locally**: cold point-source matches the old driver to ~7-8 sig figs, and a real
b11.045/m0.4 config (lat.1000) gives a physical single-baryon $C_B(t)$ with $m_\text{eff}$
plateau ~3.0, consistent with the existing $M_B\approx3.0$. Your job: build it on tuolumne
and run it over all 4 ensembles' `q00` files, producing one small HDF5 per config.
**Do not submit jobs or delete/overwrite; the user submits.**

Companion design doc: `two_baryon_gevp_contract_impl_plan_claude.md` (physics, GEVP,
references). Paths are relative to `src/gauge_gen_Nc4/` of `Grid_sdm_build` unless absolute.

## What the binary does (v1)

Reads the stored `q00` point-to-all propagator `.lime` (8 corner records/config, single
prec, from `point2all_prop_dumper_claude.cc`) and builds:
- Two-baryon GEVP matrix `C2B_snk<i>_src<jo>` (fixed-point sink; default 3x3 over corner
  pairs `{0,1},{0,3},{0,7}` = edge/face/body classes; `{0,7}` = O/M),
- Single-baryon `CB_set<s>_corner<Y>` and `CB_set<s>_avg` (p=0), and a `meta` record.
Contraction = `24^4 det(8x8)` of q00 colour blocks (reused verbatim from
`two_baryon_corr_prod_claude.cc`). **Light: I/O-bound** (read 0.68 GB, trivial dets + 8
single-baryon densities). No momentum projection yet (that is v2).

## CLI

```
two_baryon_gevp_contract_claude \
    --src <file.prop.lime>    # repeatable; one per source-smearing SET (v1: one set)
    --out <file.h5> \
    --grid 24.24.24.48 --mpi 1.1.1.1     # single rank holds the whole lattice
```
Single rank / 1 GPU per config (q00[8] double ~1.36 GB fits one MI300A GCD). It IS a Grid
GPU binary but the work is light. (A CPU Grid build would also work and avoid the GPU.)

## Build on tuolumne

```bash
cd ${ROOT}/Grid_sdm_build          # ROOT=/usr/workspace/lsd/matsumoto5/su4_32c
GRID=${ROOT}/build ./compile_two_baryon_claude.sh two_baryon_gevp_contract_claude.cc
# -> src/gauge_gen_Nc4/bin/two_baryon_gevp_contract_claude
```

## Input / output paths + ensembles

Input `q00` (from the dumper) on lustre:
`/p/lustre5/matsumoto5/two_baryon_2448_b<BETA>_m<MASS>/conf_nc4nf1_2448_b<BETA>_m<MASS>.<idx>.prop.lime`

Output (one small ~70 KB HDF5 per config), suggested:
`/p/lustre5/matsumoto5/obs_gevp_2448_b<BETA>_m<MASS>/gevp.<idx>.h5`

Config lists match the dumper's production sampling (`point2all_prop_dumper_benchmarks_claude.md`):

| mass | BETA (str) | MASS str | stride | #cfg | index range   |
|------|-----------|----------|--------|------|---------------|
| 0.1  | 10p865    | 0p1000   | 40     | 92   | 1000 .. 4640  |
| 0.2  | 10p990    | 0p2000   | 40     | 93   | 1000 .. 4680  |
| 0.3  | 11p035    | 0p3000   | 80     | 88   | 1000 .. 7960  |
| 0.4  | 11p045    | 0p4000   | 80     | 113  | 1000 .. 9960  |

(Only process configs whose `.prop.lime` exists; skip missing. ~386 total.)

## Task: write `test/run_gevp_contract/run_gevp_contract_flux_claude.sh`

- Loop the config list per ensemble; per config:
  `./bin/two_baryon_gevp_contract_claude --grid 24.24.24.48 --mpi 1.1.1.1
   --src <indir>/conf_...<idx>.prop.lime --out <outdir>/gevp.<idx>.h5`.
- **Skip if the `.prop.lime` is missing** (warn) and **skip if the output `.h5` already
  exists** (do NOT `rm`/overwrite; if a rerun is needed the user removes it).
- Light per config -> either a serial loop on 1 GPU, or fan out across the allocation's GCDs
  (e.g. `-N 8 -n 32`, one `flux run --tasks-per-node=... -g 1` per config, up to 32
  concurrent). Serial on 1 GPU over 386 configs is ~1-2 h (I/O-bound) and fits one `-t 480m`
  job; fan-out is a straightforward speedup.
- Mirror env/paths from `run_two_baryon_flux_claude.sh` (`source ${ROOT}/env.sh`,
  `HDF5_USE_FILE_LOCKING=FALSE`). `mkdir -p` the outdir. No `--rottest`. No submit / rm / kill.
- Leave a `# Submit yourself: flux batch ...` note.

## Sanity checks per ensemble

- Output HDF5 has 19 datasets: `meta`, 9 `C2B_snk*_src*`, 8 `CB_set0_corner*`, `CB_set0_avg`
  (python h5py; h5ls may be absent).
- Single-baryon `m_eff(t)=log CB_set0_avg(t)/CB_set0_avg(t+1)` should plateau near the known
  $M_B$ for the ensemble (m0.4 ~ 3.0; the others per the existing single-baryon analysis).
- `C2B_snk2_src2` (the `{0,7}` element) is the two-baryon correlator; on a real config it
  decays ~2x faster than $C_B$ (near-non-interacting).

## rsync back

The `.h5` are tiny (~70 KB); pull the whole `obs_gevp_2448_*` tree back for the downstream
python GEVP/Luescher analysis (mirror `rsync_pull_meson_conn_claude.sh` auto-discovery, or a
per-ensemble pull like `rsync_pull_prop_dumper_claude.sh`).

## Conventions (follow)
- Files you create/edit get `_claude` before the extension.
- One statement per line (shell included). No Unicode; LaTeX macros in comments.
- Never `rm`/overwrite/kill in any script; never submit the Flux job yourself.

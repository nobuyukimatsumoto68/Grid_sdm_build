# Two-baryon GEVP production + contraction — thread handoff (2026-09-17)

End-to-end handoff for the SU(4) SDM two-baryon GEVP scattering study on the four heavy $24^3\times48$
ensembles. Pipeline: **q00 dumper (tuolumne GPU) -> stage (oslic) -> GEVP contraction (dane CPU) ->
pull to local -> python GEVP/Luescher analysis**. All scripts live in
`/usr/workspace/lsd/matsumoto5/su4_32c/` (= `su4_32c`); companion design docs sit beside the code in
`Grid_sdm_build/src/gauge_gen_Nc4/`.

## Ensembles + measurement sampling

On-disk HMC checkpoints are stride 20 from lat.1000. Measurement strides chosen per ensemble:

| mass | beta (str) | stride | #cfg | index range   | single-baryon $m_\text{eff}$ (obs) |
|------|-----------|--------|------|---------------|-----------------------------------|
| 0.1  | 10p865    | 40     | 92   | 1000 .. 4640  | ~1.56 |
| 0.2  | 10p990    | 40     | 93   | 1000 .. 4680  | (per single-baryon analysis) |
| 0.3  | 11p035    | 80     | 88   | 1000 .. 7960  | (per single-baryon analysis) |
| 0.4  | 11p045    | 80     | 113  | 1000 .. 9960  | ~3.0 (validation anchor) |

Total = **386 configs**. (m0.05 @ b10p840 / m0.01 @ b10p800 DEFERRED -- need deflation.)

## Stage 1 -- q00 point-to-all dumper (tuolumne, MI300A GPU)

- Binary `point2all_prop_dumper_claude` (built on tuolumne against `su4_32c/build`). Dumps the PD
  spin-0 colour matrix $q_{00}$ (8 corner sources at $t_{src}=0$, 4 solves/source, split-grid MRHS
  `--split 1 1 1 4`, nrhs=8). Output = one SciDAC/LIME per config, 8 single-prec records, ~0.68 GB.
- **Two source operators** via one runtime knob (no code change): `SRC=smeared` (Gaussian width 3.0,
  niter 40) and `SRC=point` (width 0 -> the binary skips smearing). `PropRecord.smearWidth` tags them.
- Scripts (submit yourself; Claude does not submit):
  - per-mass batch: `submit_prop_dumper_b<beta>_m<mass>_claude.sh` (`#FLUX -N 8 -n 32 -g 1 -t 480m`)
  - driver/chain:   `run_prop_dumper_chain_claude.sh`  (all 4 ensembles; `SRC=point` for point pass)
- Hardening: in-binary graceful wall blocker (`PROP_DEADLINE_EPOCH`/`PROP_TPT_SECONDS`, ~150 s/config
  budget) + **atomic write** (`.inprogress` -> rename), so a wall kill never leaves a partial file.
- Per-config wall (measured): 61 s (m0.4) .. 78 s (m0.1); point ~59 s.
- Output (lustre5): `/p/lustre5/matsumoto5/two_baryon_2448_b<beta>_m<mass>[_point]/*.prop.lime`
- STATUS: **smeared 386/386 DONE; point 386/386 DONE.**

## Stage 2 -- stage q00 to a dane-visible filesystem (run on OSLIC)

Dane mounts only lustre1-3; the data is on lustre5; **oslic mounts both**. Copy there.
- Script: `copy_q00_to_dane_claude.sh` (run on oslic; `rsync -a`, .prop.lime only, idempotent, no rm).
  `SRC=smeared|point`, `DSTROOT` default `/p/lustre2/matsumoto5`. ~242 GB per source set.
- Output (lustre2): `/p/lustre2/matsumoto5/two_baryon_2448_b<beta>_m<mass>[_point]/`
- STATUS: **smeared staged 386/386; point staged 386/386.**

## Stage 3 -- GEVP contraction (dane, CPU) -- runs here, NOT tuolumne

I/O-bound + light (read 0.68 GB, trivial $24^4\det(8\times8)$ + single-baryon densities) -> CPU, to
spare MI300A allocation. Binary `two_baryon_gevp_contract_claude`.
- Build ON dane (grid-config native flags; the LINK ORDER must be object-before-libs for static
  libGrid.a -- already fixed in the shared compile script):
  ```
  module load clang/14.0.6-magic mvapich2/2.3.7 hdf5-serial
  cd su4_32c/Grid_sdm_build
  GRID=/usr/workspace/lsd/matsumoto5/su4_16_8_dane/build \
  GRID_SRC=/usr/workspace/lsd/matsumoto5/su4_16_8_dane/Grid \
    ./compile_two_baryon_claude.sh two_baryon_gevp_contract_claude.cc
  ```
  (Use the `su4_16_8_dane/build` CPU Grid: clang++ + mvapich2-clang-14, AVX512, Nc=4, HDF5.)
- CLI: `two_baryon_gevp_contract_claude --src <file.prop.lime> [--src ...] --out <file.h5>
  --grid 24.24.24.48 --mpi 1.1.1.1`. `--src` is repeatable (one per source-smearing set).
- SLURM scripts (account `-A latticgc`, matches the other dane scripts):
  - batch: `submit_gevp_contract_dane_claude.sh` (`#SBATCH -N1 --cpus-per-task=112 -t 12:00:00
    --requeue`; ENS/SRC/DATAROOT env knobs; serial config loop; skip-if-exists).
  - driver: `run_gevp_contract_dane_claude.sh` (one job PER ensemble; knobs `SRC` smeared|point,
    `QUEUE` pbatch|pdebug, `TLIMIT`, `KEYS="m0p2000 m0p3000"` to submit a subset / resume).
- Per-config wall: ~13 s solo, ~80 s under 4-way lustre2 I/O contention.
- Output (lustre2): `/p/lustre2/matsumoto5/obs_gevp_2448_b<beta>_m<mass>[_point]/gevp.<idx>.h5` (~84 KB).
- STATUS (2026-09-18): **point-sink contraction DONE** -- smeared 386/386 AND point 386/386
  (all 4 ensembles complete, both source sets); pulled + analysed locally (dE preliminary).

## Stage 3b -- smeared SINK rerun (Chunk 8, NEW code)

`two_baryon_gevp_contract_claude.cc` now also does a **covariant Gaussian sink smearing** of q00
(in addition to the point sink), emitting `C2Bss_snk<i>_src<jo>` + `CBss_set<s>_*` datasets. This
needs the **gauge field**, so the contraction loads the NERSC config.
- New CLI: `--config <nersc_gauge_file>` and `--sink-smear <w> <N>` (match the source: `3.0 40`).
  Without them it behaves exactly as before (point sink only; point-sink numbers unchanged --
  validated on cold: C2B_07_07 identical after the refactor).
- **WHERE to run:** the gauge configs live on **lustre5**
  (`/p/lustre5/matsumoto5/conf_nc4nf1_2448/conf_nc4nf1_2448_b<beta>_m<mass>/..._lat.<i>`), and the
  q00 also still sit on lustre5 (they were copied, not moved). So run the smeared-sink pass where
  BOTH are visible -- i.e. **tuolumne (reads lustre5)** or oslic-visible compute -- NOT dane
  (dane sees only lustre1-3; staging ~525 GB of gauge to lustre2 is NOT worth it). The contraction
  is still light (CPU-bound work; on a GPU node it just uses one device briefly).
- Run script to write (mirror `run_gevp_contract_dane_claude.sh`): per config, map the gevp index
  to its NERSC gauge file and pass `--config <that> --sink-smear 3.0 40`, `--src` the q00
  set(s), `--out` a NEW obs dir e.g. `obs_gevp_2448_b<beta>_m<mass>_ss/gevp.<idx>.h5` (keep the
  smeared-sink outputs separate from the point-sink ones). Skip-if-exists; no rm/overwrite.
- Output HDF5 then has 37 datasets (19 point-sink + 9 `C2Bss_*` + 9 `CBss_*`); `meta` carries
  `sinkSmearWidth`/`sinkSmearNiter`. Pull with the existing `rsync_pull_gevp_claude.sh` (its
  filter already matches `obs_gevp_2448_*`).
- STATUS: code DONE + validated locally; **rerun to run** (choose smeared and/or point source set).

### Output format -- gevp.<idx>.h5, 19 datasets (Grid Hdf5Writer)
- `meta`
- `C2B_snk<i>_src<jo>` for i,jo in 0..2 (3x3 two-baryon GEVP matrix; corner pairs
  {0,1}=edge, {0,3}=face, {0,7}=body(=O/M)); each a length-48 complex vector `.data` (re,im).
- `CB_set0_corner0..7` (single-baryon per corner) and `CB_set0_avg` (corner-averaged, p=0).
- Sanity: single-baryon $m_\text{eff}$ plateaus at the ensemble $M_B$ (m0.1 ~1.56, m0.4 ~3.0);
  the two-baryon `C2B_snk2_src2` ({0,7}x{0,7}) decays ~2x faster (near-non-interacting). VALIDATED.
- smeared vs point kept in SEPARATE dirs (`_point`); combine as a wider GEVP operator basis downstream.

## Operational notes

- **Preemption:** dane pbatch jobs are frequently `CANCELLED ... DUE TO PREEMPTION` (~every 10-30 min).
  This is NOT a code or memory bug (memory verified: 1 GiB Grid shm via shmget + ~0.76 GB q00 doubles,
  no OOM). Frequency suggests low/standby-priority scheduling (check `sacct -o QOS,Priority`, `sshare -M
  dane`, `scontrol show partition pbatch | grep -i preempt`; a bank with a live dane allocation would
  fix it). WORKAROUND: `--requeue` is set but this QOS cancels rather than requeues, so RESUME by
  re-running the driver -- skip-if-exists continues where it stopped:
  `KEYS="m0p2000 m0p3000" bash run_gevp_contract_dane_claude.sh`. A dependency-chain (CHAINLEN) auto-
  resume can be added if wanted.
- **Never rm/overwrite/kill** in any script; the user removes files + submits jobs. The contraction
  `.h5` write is NON-atomic, so never run two jobs on the SAME ensemble+set concurrently (collision).
- Memory allocation verified correct at all levels (SLURM whole-node, Grid default 1 GiB shm, ~2-3 GB
  peak vs ~256 GB/node).

## Pull to local + analysis

Outputs are tiny (~84 KB each; ~30 MB per full ensemble set). Pull the `obs_gevp_2448_*` trees from
lustre2 (dane-visible; via oslic or directly to your local box):
```
rsync -avP <lc_user>@<host>:/p/lustre2/matsumoto5/'obs_gevp_2448_b*_m*' <local_dest>/
# add the _point dirs once the point contraction is done (they are separate dirs)
```
Analysis (local, python/h5py): read the 19 datasets per config; build the GEVP matrix C(t) from
`C2B_snk*_src*` (+ optionally fold in the point set as extra operators -> asymmetric/unmatched GEVP,
source and sink bases need not match, see `two_baryon_gevp_contract_impl_plan_claude.md`); solve
$C(t)v=\lambda(t,t_0)C(t_0)v$ for the energy levels; single-baryon $M_B$ from `CB_set0_avg`; then
Luescher to get the scattering phase shift.

## Key scripts (all under su4_32c)
- Dumper:  `submit_prop_dumper_b<beta>_m<mass>_claude.sh`, `run_prop_dumper_chain_claude.sh`,
           `submit_prop_dumper_b11p045_m0p4_debug_claude.sh` (pdebug validation)
- Stage:   `copy_q00_to_dane_claude.sh` (oslic)
- Contract:`submit_gevp_contract_dane_claude.sh` (batch), `run_gevp_contract_dane_claude.sh` (driver)
- Design:  `Grid_sdm_build/src/gauge_gen_Nc4/two_baryon_gevp_contract_{impl_plan,remote_handoff}_claude.md`,
           `point2all_prop_dumper_{impl_plan,remote_handoff,benchmarks}_claude.md`

## Immediate next steps
1. Finish smeared contraction: resume m0.2 (-2) / m0.3 (-26) until 93/88 (preemption -> repeat).
2. Run point contraction: `SRC=point bash run_gevp_contract_dane_claude.sh` -> obs_gevp_2448_*_point.
3. Pull obs_gevp_2448_* (+ _point) to local; run the GEVP/Luescher analysis. First look: smeared m0.4
   (113 configs complete) for an early two-baryon signal.

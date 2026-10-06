# Momentum-projected single-baryon contraction: remote handoff (tuolumne)

Adds moving-frame (and finer rest-frame $E/T_2$) single-baryon correlators to the light-mass full-PD
dumps, using A. S. Meyer's (LLNL) cubic / little-group irrep 4q operators (private communication
2026-09-28). **Contraction only; the dumps are the existing / planned Stage 1 outputs.** Plan and
validation record: `baryon_variants_mom_impl_plan_claude.md`. The production contraction
(`baryon_variants_corr_claude.cc`, `test/run_bvar_contract/`) is untouched. **The user submits all jobs;
never rm / overwrite / kill.**

## What it computes

Per config and sink type (point, smeared): for all $35^2=1225$ elementary pairs the site density
$d_{AA'}(x)$ (same kernel as production), FFT in space to 81 momenta (frames 000, 100, 110, 111, 200,
210, 211, all directions), source phase $e^{+i\mathbf p\cdot\mathbf y}$, then Aaron's operator correlators

$$
C_{ij}(\mathbf p,t)=\sum_{A,A'} c^i_A(\mathbf p)\,\overline{c^j_{A'}(\mathbf p)}\,C_{AA'}(\mathbf p,t),
$$

averaged over the components (rows x momentum directions) of each (frame, irrep) class: 28 classes,
2111 correlators per sink type. In flight parity mixes, hence all 1225 pairs (~22x the det terms of
the production block mode).

## Files (git sync of `Grid_sdm_build`, all under `src/gauge_gen_Nc4/`)

- `baryon_variants_mom_corr_claude.cc` (binary source; includes the unchanged `baryon_pd_common_claude.h`)
- `baryon_mom_ops_claude.txt` (operator groups, 81 momenta; from `aaron_ops_convert_claude.py`)
- `baryon_variants_terms_all_claude.txt` (all 1225 pair term tables, 2.4 MB)
- `baryon_variants_ops_claude.txt` (already there; legacy p = 0 operators)
- `test/run_bvar_mom_contract/run_bvar_mom_contract_flux_claude.sh` (run script)
- `baryon_variants_mom_check_claude.py` (optional, local checks)

## Build (on tuolumne)

```
cd /usr/workspace/lsd/matsumoto5/su4_32c/Grid_sdm_build && source ../env.sh
GRID=/usr/workspace/lsd/matsumoto5/su4_32c/build ./compile_two_baryon_claude.sh baryon_variants_mom_corr_claude.cc
```
The FFT uses Grid's `FFT` class (hipFFT backend; Grid's configure adds `-lhipfft` for HIP builds).
The momentum path asserts a single MPI rank (`--mpi 1.1.1.1`, as the production contraction).

## Run

One ensemble and one source type per job; knobs `ENS` (`m0p1000` / `m0p0500` / `m0p0100`), `SRC`
(`point` default / `smeared`), `CLIST` (default: every dump present), `PHASECHECK`, `MOMRAW`, `NSLOTS`.

Per config the script
- skips it if `obs_bvar_2448_<ens>[_point]_ss_mom/bvarmom.all.<idx>.h5` exists;
- passes `--no-p0` if the production `obs_bvar_2448_<ens>[_point]_ss/bvar.all.<idx>.h5` exists (only the
  momentum datasets are written); otherwise it also writes the p = 0 `Cel_`/`Cop_` datasets (same
  `--pairs sector` selection as production, 109 pairs), so no separate production contraction is
  needed for the new dumps;
- writes `<out>.inprogress.h5` and renames it only after a clean exit.

Order:
1. **Validation, m0.1 smeared (dumps and production p = 0 results already exist):**
   `ENS=m0p1000 SRC=smeared CLIST=1000 PHASECHECK=3 flux batch run_bvar_mom_contract_flux_claude.sh`
2. m0.1 smeared, full set (92 configs, all `--no-p0`): `ENS=m0p1000 SRC=smeared flux batch ...`
3. m0.1 point, once the point dumps `bvar_2448_b10p865_m0p1000_point/` exist (per
   `bvar_point_source_handoff_claude.md`): `ENS=m0p1000 SRC=point flux batch ...`
4. m0.05 and m0.01, both source types, once their dumps exist (the production dumper script only knows
   m0p1000..m0p4000 so far; it needs the `b10p840_m0p0500` / `b10p800_m0p0100` cases added, see
   `bvar_point_source_handoff_claude.md`).

## Per-config sanity (log `bvarmom.all.<idx>.log`)

- `combined compute list: 1225 pairs (... legacy p = 0, 1225 momentum)`
- per sink: `max rel |FFT(p=0) - sliceSum| = O(1e-15)` (built-in check on every pair)
- validation config: `phase-check summary: max rel = O(1e-15) PASSED`
- `wrote mom_meta + 2111 Cmom datasets (+ smeared-sink copies)`
- record the `contraction+sliceSum` and `FFT+projection` times of the first config (local TITAN V,
  16^3x8: 1.4 s + 7.1 s per sink pass; the FFT part is per-call overhead, ~6 ms per pair).

## Output

`/p/lustre5/matsumoto5/obs_bvar_2448_b<beta>_m<mass>[_point]_ss_mom/bvarmom.all.<idx>.h5`, about 1.6 MB per
sink type per config. Datasets (CorrFile, compound (re, im), length T = 48):
- `mom_meta` (attributes): `momOpsFile`, `nmom`, `momList` (`k:px,py,pz` joined by `|`), `classList`
  (`<frame>_<irrep>:<ncomp>:<uid>,<uid>,...` joined by `|`), `raw`, `p0Written`
- `Cmom_<frame>_<irrep>_<uid_i>_<uid_j>` (point sink), `Cmomss_...` (smeared sink); component average
- with `MOMRAW=1`: `Cmomraw_<frame>_<irrep>_c<comp>_<uid_i>_<uid_j>`, `Cmomrawss_...`
- `meta` (as production) and, when p = 0 is written, `Cel_*`, `Cop_*`, `Celss_*`, `Copss_*`

Frames: `000` (irreps `a1p eem eep t1m t1p t2m t2p`), `100`/`200` (`a1x a2x b1x b2x eex`),
`110` (`a1x a2x b1x b2x`), `111` (`a1x a2x eex`), `210`/`211` (`a1x a2x`). Operator IDs (`uid`) and their
spin content: `baryon_mom_ops_claude.txt`; rest frame <-> our 9 operators: `0p_22`=a1p.00040,
`1p_22`=t1p.00035, `1m_31/1m_13`=t1m.00029/00033, `2p_40/2p_22/2p_04`=eep+t2p 00036+00039 / 00034+00038 /
00031+00042, `2m_31/2m_13`=eem+t2m 00030+00041 / 00032+00037.

## Local side

Pull with `./rsync_pull_bvar_claude.sh` (its filter `obs_bvar_2448_*/` also matches the `_mom` dirs);
aggregation in `bvar_h5_claude.py` (grid-claude-54 extends its directory regex for `_ss_mom`).

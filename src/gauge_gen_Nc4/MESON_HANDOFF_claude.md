# Meson measurement + 2D-C combination — handoff

Self-contained pointer doc for the momentum-projected meson pipeline. For the
full design rationale and open questions see `meson_momproj_impl_plan_claude.md`.

## Goal

Measure the momentum-projected CONNECTED meson correlator and combine it with the
existing DISCONNECTED loop-loop correlator to form the physical flavor-singlet
correlator, per Rebbi's convention
$$
C_\text{phys}(\vec p, t) = 2\,D(\vec p, t) - C_\text{conn}(\vec p, t)
$$
(the `2*corr - ccorr` of `makecorr_96_d.f90:353`). Ensemble:
`obs_nc4nf1_2448_b11p045_m0p4000` ($24^3\times48$, $NS=24$, $NT=48$, Mobius DWF
$L_s=16$, $b=1.5$, $c=0.5$, $M5=1.5$, mass $=0.4$, antiperiodic in time).

## Files (all in `src/gauge_gen_Nc4/`)

- `baryons_0000_dirac_claude.cc` — the meson driver. Point source at the origin,
  CG/Schur solve, connected meson contraction
  `tr[G5 q1^dag G5 Gsnk q2 Gsrc^dag]` (`:147`) projected onto 17 momenta
  (`MakePhase` + `mom_table` at `:18-61`; momentum loop in `MesonTrace_hdf`
  writing `<channel>_t_p<idx>`). Baryon sector is disabled (`#if 0` at `:269`,
  call commented at `:507`). 10 channels: I_I, G5_G5, GTG5_GTG5, GXG5_GXG5,
  GYG5_GYG5, GZG5_GZG5, GT_GT, GX_GX, GY_GY, GZ_GZ.
- `tmp_claude.sh` — LOCAL build + cold-config smoke test (`--grid 8.8.8.16`,
  fits one GPU), runs `verify_momproj_p0_claude.py`, tees to
  `baryons_momproj_cold_claude.log`.
- `verify_momproj_p0_claude.py` — checks `<channel>_t_p0` == `<channel>_t`
  (since `MakePhase({0,0,0})=1`). PASSED for all 10 channels.
- `run_meson_momproj_claude.sh` — CLUSTER production runner (FLUX, multi-GPU),
  loops the config list, one inversion per config, writes
  `meson_momproj_<cfg>_lat.<i>.h5` into `meson_momproj_out_claude/`.
- `meson_combine_claude.py` — analysis: read connected h5, fold, class-average,
  read disc `corr_d.p*`, form `2D - C`.
- `meson_momproj_impl_plan_claude.md` — design + open questions.

## Momentum table (must match `average_trace2.f90:74`; `3` = -1)

17 momenta, indices 0-16, in units $2\pi/NS$. $p^2$ classes (stored indices):
- p0: {0}; p1: {1,2,3}; p2: {4,5,6,7,8,9}; p3: {10,11,12,13}; p4: {14,15,16}.

Rebbi's `pcomp(m)/2` equals the stored-index count per class, so the connected
class average is a plain mean over those indices (real part).

## Build

LOCAL (CUDA/TITAN V), mirrors `src/gauge_gen_Nc4/Makefile`:
```
bash ../../compile_two_baryon_claude.sh baryons_0000_dirac_claude.cc
```
CLUSTER (mi300a/ROCm): rebuild there with the cluster's `grid-config`; the local
`./bin` binary is CUDA-only.

## Run

LOCAL smoke test (wiring only, ~seconds, <0.5 GB):
```
./tmp_claude.sh        # -> baryons_momproj_cold_claude.log ; expect "ALL MATCH"
```
Single-GPU full $24^3\times48$ needs ~18-20 GB and OOMs an 11.8 GB TITAN V, so
production must be multi-GPU.

CLUSTER production (fill the TODO-marked paths first: `ENV_SH`, `BIN`,
`CFGPATH`, `cfgfilename`, `CLIST`):
```
./run_meson_momproj_claude.sh
```
Positional binary args are `<config> <M5> <mass> <outfile>`; cold config if the
first arg starts with `-`.

## Analyze (2D-C)

`meson_combine_claude.py` API:
- `read_connected_h5(file, channel)` -> `(17, NT)` complex `C(p_idx,t)`.
- `connected_pclass(file, channel)` -> `(5, NT/2+1)` real, folded + class-averaged.
- `read_corr_d(file)` -> `(nconf, NT/2+1)` real disc `D` for one p-class file.
- `combine_2D_minus_C(D, C, rel_norm=1.0)` -> `2*rel_norm*D - C`.
- `CHANNEL_TO_TRACEDIR` maps each meson channel to its `trace_<G>/` dir.

Connected-side smoke test (no production data needed):
```
python3 meson_combine_claude.py test/run_baryon_momproj_claude/ColdConfig.h5
```
Disc files live at `obs_nc4nf1_2448/<ens>/trace_<G>/corr_d.p{0..4}` (per-config
blocks of `t  D(t)`, t=0..24).

## Status

- DONE: momentum projection in Grid (verified); GZ channel sink-Gamma bug fixed
  (`:155-164`, originals commented); baryon sector disabled; production run
  script; combination module with connected-side test passing.
- PENDING: cluster production run on real configs; end-to-end `2D-C` +
  jackknife + plots; integrate into `analyze_corr_claude.ipynb`.

## Open questions (see plan for detail)

- OQ1 (UNRESOLVED): relative D-vs-C normalization (`rel_norm`; the factor-3 at
  `makecorr_96_d.f90:231`). Calibrate before trusting absolute `2D-C`.
- OQ2: confirm disc `ns^3/NT` and `/(pcomp/2)` are matched on the connected side.
- OQ3: confirm channel <-> `trace_<G>` mapping (`CHANNEL_TO_TRACEDIR`).
- OQ4 (RESOLVED): GZ/GZG5 sink Gamma fixed.

# Baryon operator variants: POINT-SOURCE pass (addendum to `baryon_variants_production_handoff_claude.md`)

Goal: add the **point-source** column to the single-baryon variant data, so every operator has the
full source x sink smearing matrix {point, smeared} x {point, smeared}. Locally this enables a
$2n\times2n$ GEVP per $J^P$ block (operators x smearings), e.g. $2^+$: 3 ops x 2 smearings = 6x6.
Everything else (ensembles, action, build, binaries, scripts, sanity checks, conventions) is exactly
as in `baryon_variants_production_handoff_claude.md`. **No code change. The user submits all jobs;
never rm / overwrite / kill.**

## What already exists (smeared source, w=3 N=40)
- Dumps: `/p/lustre5/matsumoto5/bvar_2448_b<beta>_m<mass>/conf_nc4nf1_2448_b<beta>_m<mass>.<idx>.pdfull.lime`
- Contractions (point + smeared sink): `/p/lustre5/matsumoto5/obs_bvar_2448_b<beta>_m<mass>_ss/bvar.all.<idx>.h5`
- All four ensembles complete (92 / 93 / 88 / 113 configs), pulled and analysed locally.

## Stage 1b -- point-source full-PD dump
Same binary, same layout, **`--width 0`** (the dumper then skips the source smearing):
```
point2all_prop_dumper_full_claude <config> 1.5 <mass> <outfile> --grid 24.24.24.48 --mpi 2.2.2.4 \
    --split 1 1 1 4 --width 0 --nsrc 1
```
- Use the existing `SRC=point` knob of `run_bvar_dumper_flux_claude.sh` (add it if missing: `SRC=point`
  -> `--width 0`).
- **OUTPUT MUST GO TO A SEPARATE DIR** (filenames are identical to the smeared pass, so writing into the
  same dir would make skip-if-exists skip every config):
  `/p/lustre5/matsumoto5/bvar_2448_b<beta>_m<mass>_point/conf_nc4nf1_2448_b<beta>_m<mass>.<idx>.pdfull.lime`
- Cost: same as the smeared pass, 16 solves/config (~35-45 s), 1.36 GB/config.
  **Storage: another ~525 GB on lustre5 for all four ensembles.** The dumps are only needed until the
  contraction of a config is done; the user may delete them afterwards (do not delete them yourself).

## Stage 2b -- contraction of the point-source dumps (point + smeared sink)
Exactly the production command, on the `_point` dumps, with the sink smearing on:
```
baryon_variants_corr_claude --grid 24.24.24.48 --mpi 1.1.1.1 \
    --src <..._point/...<idx>.pdfull.lime> --out bvar.all.<idx>.h5 \
    --terms baryon_variants_terms_claude.txt --ops baryon_variants_ops_claude.txt \
    --op all --pairs block --config <nersc gauge of <idx>> --sink-smear 3.0 40
```
- **Separate obs dir:** `/p/lustre5/matsumoto5/obs_bvar_2448_b<beta>_m<mass>_point_ss/bvar.all.<idx>.h5`
- Keep `--pairs block --op all` and `--sink-smear 3.0 40` identical to the smeared pass, so the two
  source columns carry the same datasets (`Cel_*`, `Cop_*`, `Celss_*`, `Copss_*`).

## Validation anchor (run first: **m0.1 / b10p865, lat.1000 only** -- user priority)
The q00 production already has a point source (`obs_gevp_2448_b10p865_m0p1000_point`, and the `_ss`
rerun holds it as set1). The PD source-spin-0 column of the point-source full dump uses the same
pre-rotated source as the q00 point dump, so on the same config:
- `Cop_2p_40_2p_40`   = 24 x `CB_set1_corner0`   (point source, point sink)
- `Copss_2p_40_2p_40` = 24 x `CBss_set1_corner0` (point source, smeared sink)
of `obs_gevp_2448_b10p865_m0p1000_ss/gevp.1000.h5`, to machine precision (the smeared pass agreed
with set0 to 7e-16). **Order: full m0.1 set first (92 configs), then m0.4, m0.3, m0.2.**

## Light ensembles m0.05 (b10p840) and m0.01 (b10p800) -- BOTH source smearings, new
The user wants these too. Nothing exists yet for them in this pipeline (no q00, no full-PD, no bvar).
Run BOTH passes per ensemble: smeared source (`--width 3.0 --niter 40`, dirs `bvar_2448_b<beta>_m<mass>/`,
`obs_bvar_2448_b<beta>_m<mass>_ss/`) and point source (`--width 0`, `_point` / `_point_ss` dirs as above).
- **DECISION (user, 2026-10-03): NO deflation.** Plain Schur-RB CG, tol 1e-9, zero guess, as in the
  heavy runs. Reason: the disc study tuned shift-invert deflation at exactly m0.01 (b10.8, lat.758;
  `Grid/examples/disc_mrhs_defl_bench_results_claude.md`): the 100-mode eigensolve costs ~643 s/config
  and saves only ~62 s per 16 RHS (168 -> 106 s). The disc amortises that over ~1500 solves/config; the
  full-PD dumper has only 16 solves per config per source smearing, so deflation is a net loss.
  Undeflated reference at m0.01: ~3000 CG iterations, ~13 s per solve (plain double CG), i.e. a few
  minutes per config per pass. The connected-meson study ran these ensembles undeflated too (166 / 394
  configs). Still report CG iterations and wall for the first config of each light ensemble, as a check.
- Find the gauge stream on lustre5 (`conf_nc4nf1_2448/conf_nc4nf1_2448_b10p840_m0p0500/`,
  `..._b10p800_m0p0100/`; the meson obs dirs `obs_nc4nf1_2448_b10p840_m0p0500` / `_b10p800_m0p0100`
  show which indices were used), report available indices and propose a stride as in the benchmarks doc.
- Watch memory: the 24^3x48 solve is unchanged in size; light masses only cost iterations.

## Per-config sanity (as for the smeared pass)
- dump: `grep -a -c ildg-binary-data` = 16; `<store>q00</store>` .. `q33` once each;
  `<smearWidth>0</smearWidth>` (confirms the point source).
- contraction log: `wrote ... 42 Cel + 19 Cop datasets` (+ the same counts for the `ss` copies).

## Local side (already prepared)
- Pull: `./rsync_pull_bvar_claude.sh` (its filter `obs_bvar_2448_*/` picks up the `_point_ss` dirs).
- Aggregate: `python3 bvar_h5_claude.py` -> `bvar_h5_out_claude/bvar_b<beta>_m<mass>_point_ss.h5`
  next to the existing `..._ss.h5` (directory regex extended for `_point`).
- Analysis: smearing x operator GEVP per $J^P$ block in `bvar_spectrum_claude.py` (to be added when
  the data lands).

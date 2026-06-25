# Conn binary: internal config-loop mode (like disc) — impl plan

## Goal
Make `baryons_0000_dirac_claude.cc` (CONNECTED meson) process a whole ensemble in
ONE invocation by looping configs INTERNALLY, exactly like
`disc_multipleGamma_binary_claude.cc`. This eliminates the per-config `flux run`
relaunch tax: today each config is its own MPI job, so even a *skip* costs ~6 s
(MPI init + Grid setup + teardown) although the skip decision itself is ~0.65 s.
With an internal loop the Grid/MPI init happens ONCE and skips become ~instant
(an `output_complete` check before any gauge read), matching disc.

## Why disc is cheap and conn is not (measured)
- disc: one `flux run` -> one resident process -> internal `for(conf)` with
  `std::filesystem::exists()` skip BEFORE the NERSC read (`disc:191-201`).
- conn: bash loops, one `flux run` per config; ~6 s launch/teardown per config
  regardless of skip (b11p045 log: conf 1000 rank-up 12:48:47, conf 1020 12:48:53).

## Design (mirror disc)
Per-config loop inside `main`, after building the grids ONCE:
```
detect conf_min/conf_max/interval + lat_prefix from --dir (disc logic)
build grids, point source (config-independent) ONCE
read deadline (env) for the in-binary wall blocker
for(conf = conf_min; conf < conf_max; conf += interval):
    outfile = obsdir + "/mesons_conn." + conf + ".h5"
    if output_complete(outfile) on boss, GlobalSum, continue        # cheap skip
    wall-blocker: boss checks now + 1.2*max_dur > deadline, GlobalSum, break
    t0 = time()
    NerscIO::readConfiguration(Umu, ...)                            # per config
    build MobiusFermionD action (mass, M5, b=1.5, c=0.5, antiperiodic)
    Solve(point_source -> PointProps)
    Hdf5Writer WR(outfile) on boss; MesonTrace_hdf(WR, PointProps, PointProps)
    measure dur, GlobalSum, update max_dur
```
Reuse existing `output_complete`, `MesonTrace_hdf`, `Solve`, `PointSource`,
`MakePhase`. The point source (delta at origin) is gauge-independent -> build once.
The fermion action depends on the gauge -> rebuild per config (same as disc's
FermAct-per-config). Writer scoped per config (open/write/close), like disc's
ScidacWriter.

## CLI
Add disc-style options for loop mode (conn needs no `--beta`; valence mass only):
`--dir <latdir> --obsdir <obsdir> --mass <m> --M5 <M5>` (Ls hardcoded 16 as now).
KEEP the existing positional path `<config> <M5> <mass> <outfile>` (single config /
cold config) as a fallback so the smoke test still works — preserve, don't delete.
Selection: if `--dir` is present -> loop mode; else positional/cold (current code).

## Wall blocker (in-binary, like disc)
Env `MESON_DEADLINE_EPOCH` (epoch s; 0 disables) + optional `MESON_TPT_SECONDS`
bootstrap (~20 s conn compute). Boss decides, `GlobalSum` so all ranks break
together. Mirrors the disc blocker added in `disc_..._claude.cc`.

## Files to modify
1. `baryons_0000_dirac_claude.cc` — add `<getopt.h>`/`<ctime>`/`<cstdlib>` (ctime/cstdlib
   already added), a `ParseArgsDir` for `--dir/--obsdir/--mass/--M5`, the config
   auto-detect block, and the internal loop in `main` (positional path preserved).
   Files: baryons_0000_dirac_claude.cc
2. `workdir/submit_meson_momproj_tuolumne_claude.sh` — becomes a SINGLE `flux run`
   with `--dir <latdir> --obsdir <obsdir> --mass <m> --M5 <M5>` (like submit_disc),
   no bash config loop, no per-config blocker (now in binary); export
   MESON_DEADLINE_EPOCH from `flux job timeleft`. Files: submit_meson_momproj_tuolumne_claude.sh
3. Drivers `run_meson_momproj_tuo_claude.sh` + `run_meson_momproj_newens_tuo_claude.sh`
   — pass `--dir/--obsdir/--mass/--M5` env to the new submit (mostly unchanged;
   they already pass latdir/obsdir/mass/M5). Files: the two run_*.sh
4. OBSOLETED by this change (conn now self-skips cheaply like disc):
   `submit_meson_momproj_finalize_claude.sh` + the conn half of
   `finalize_old_ensembles_claude.sh` (CONFLIST no longer needed). Decide: retire or
   keep. Files: finalize_*.

## Open questions
- OQ1: CLI — add `--dir` loop mode and KEEP the positional single-config/cold mode
  as fallback (recommended), or replace positional entirely?
- OQ2: After loop mode, the conn CONFLIST finalizer is redundant. Retire
  `submit_meson_momproj_finalize_claude.sh` and drop the conn calls from the
  finalizer (disc-style resubmit suffices), or keep them?
- OQ3: STRIDE/thinning — the per-config submit supported a STRIDE env. Loop mode
  uses the native interval (like disc). Keep an optional `--stride`? (disc has none.)
```

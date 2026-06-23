# Momentum-projected connected meson + 2D-C combination — implementation plan

## Goal / physics

Add spatial momentum projection to the connected meson correlator computed in
`baryons_0000_dirac.cc`, so it can be combined with the existing disconnected
loop-loop correlator to form the physical flavor-singlet correlator.

Per Rebbi's convention the physical singlet correlator is, per momentum class,
$$
C_\text{phys}(\vec p, t) = 2 D(\vec p, t) - C_\text{conn}(\vec p, t)
$$
(the `2*corr - ccorr` line in `makecorr_96_d.f90:353`, currently commented out).
Here $D$ is the disconnected loop-loop correlator and $C_\text{conn}$ the connected
meson correlator. The combination will be done in Python (Path B, user choice),
not by re-enabling the Fortran.

This is the SU(4), $N_f=1$ ensemble `obs_nc4nf1_2448_b11p045_m0p4000`, lattice
$24^3\times48$ ($NS=24$, $NT=48$).

## How the existing pieces fit

Disconnected (already working):
- `disc_multipleGamma_binary_claude` (Grid) writes the loop trace
  $T_\Gamma(\vec x,t)=\text{tr}[\Gamma S(x,x)]$ over the whole lattice as Scidac
  binary, one $\Gamma$ per `trace_*/` directory.
- `average_trace2.f90` spatial-projects $\tilde T(\vec p,t)=\frac{1}{NS^3}\sum_{\vec x}
  e^{i\vec p\cdot\vec x}T(\vec x,t)$ over the fixed 17-momentum table
  (`average_trace2.f90:74` ptable, `:130` projection, `:143` the $1/NS^3$).
  Only the real part of $T$ is kept (`:114`).
- `makecorr_96_d.f90` time-correlates into the disc correlator
  $D(\vec p,t)=\frac{ns^3}{NT}\sum_{t'}\tilde T(t')\tilde T(t'+t)$
  (`:240-243`, $ns=NT/2=24$), averages over $p^2$ classes p0..p4 (`:247-255`),
  writes `corr_d.p0 .. corr_d.p4` (real part of disc correlator vs $t=0..NT/2$).

Connected (to add):
- `baryons_0000_dirac.cc:147` forms
  $C(t)=\text{tr}[G_5\,q_1^\dagger G_5\,\Gamma_\text{snk}\,q_2\,\Gamma_\text{src}^\dagger]$
  but `sliceSum(meson_CF, meson_T, Tdir)` at `:168` only projects $\vec p=0$.

## The 17-momentum table (must match Rebbi exactly)

From `average_trace2.f90:74` (`3` denotes $-1$):

```
index  px py pz      index  px py pz
  0     0  0  0         9     0  1 -1
  1     1  0  0        10     1  1  1
  2     0  1  0        11    -1  1  1
  3     0  0  1        12     1 -1  1
  4     1  1  0        13     1  1 -1
  5     1  0  1        14     2  0  0
  6     0  1  1        15     0  2  0
  7     1 -1  0        16     0  0  2
  8     1  0 -1
```

$p^2$ classes and component counts `pcomp=(/1,6,12,8,6/)` (`makecorr_96_d.f90:137`):
- p0: index 0            (1 component)
- p1: indices 1-3        (stored 3, count 6 with c.c. partners)
- p2: indices 4-9        (stored 6, count 12)
- p3: indices 10-13      (stored 4, count 8)
- p4: indices 14-16      (stored 3, count 6)

Momentum unit is $2\pi n/NS$ with $NS=24$. In Grid this is `MakePhase` with
`mom = {nx,ny,nz}` (the integers above); `MakePhase` already exists in
`Mobius_mesons_xt.cc:50-63` and builds $e^{i\vec p\cdot\vec x}$ with
`TwoPiL = 2*pi/L`.

## Files to modify

- `baryons_0000_dirac_claude.cc` (new, copy of `baryons_0000_dirac.cc`):
  add `MakePhase`, momentum-project the connected meson over the 17 momenta,
  write `C(channel, p_index, t)` to HDF5.
- `analyze_corr_claude.ipynb` (Python, Path B): read `corr_d.p*` (disc) + the new
  connected HDF5, fold + class-average the connected, form $2D-C$.

## Implementation chunks

### Chunk 1 — momentum-projected connected meson in Grid  [DONE, verified 2026-06-11]
Files: `baryons_0000_dirac_claude.cc`
Status: implemented (`MakePhase` + `mom_table` at lines 18-61; momentum loop in
`MesonTrace_hdf` writing `<channel>_t_p<idx>` for the 17 momenta). Cold-config
smoke test on `--grid 8.8.8.16` PASSED: `verify_momproj_p0_claude.py` confirms
`<channel>_t_p0` == `<channel>_t` for all 10 channels (ALL MATCH). Build/run via
`tmp_claude.sh` -> `baryons_momproj_cold_claude.log`. Single-GPU full 24^3x48 is
~18-20 GB and OOMs an 11.8 GB TITAN V; production must run multi-GPU on cluster.
- Add `MakePhase(Coordinate mom, LatticeComplex &phase)` (copy from
  `Mobius_mesons_xt.cc:50-63`).
- Add the 17-entry integer momentum table as a `static const int mom_table[17][3]`
  matching `average_trace2.f90:74`.
- In `MesonTrace_hdf`, after forming `meson_CF` (`:147`), loop the 17 momenta:
  build `phase`, `sliceSum(meson_CF * phase, meson_T, Tdir)` -> `C(p_idx, t)`
  for $t=0..NT-1$ (full length, complex). Keep the existing $p=0$ path too.
- Write one dataset per (channel, momentum), e.g. `<channel>_t_p<idx>`, holding
  the complex `C(p_idx, t)`. Preserve the existing `_t` / `_x` datasets
  (comment-out style: leave originals, add the new momentum loop alongside).
- Build + smoke-test on the cold config (writes to an `.h5`), check $p=0$
  dataset reproduces the current `_t` output bit-for-bit.

### Production run (prerequisite for real Chunk 2 data)  [script written]
Files: `run_meson_momproj_claude.sh`
FLUX/multi-GPU runner looping the config list, one inversion per config, writing
`meson_momproj_<cfg>_lat.<i>.h5` (with `<channel>_t_p<idx>`) into
`meson_momproj_out_claude/`. Site-specific paths marked TODO (ENV_SH, BIN,
CFGPATH, cfgfilename, CLIST). Must be (re)built on the cluster's accelerator;
the local `./bin` build is CUDA/TITAN V only.

### Chunk 2 — Python 2D-C combination
Files: `meson_combine_claude.py` [DONE, connected-side smoke-tested],
       `analyze_corr_claude.ipynb` (already has `read_corr_d`, `read_corr_d_p`,
       `Jackknife`; extend, do not clobber) [PENDING: notebook cells].
`meson_combine_claude.py` provides: `read_connected_h5`, `fold_connected`,
`connected_pclass` (fold + class-average to (5, NT/2+1)), `read_corr_d`,
`combine_CminusD(D, C, rel_norm, nf)` and `CHANNEL_TO_TRACEDIR`. Connected-side
smoke test on the cold h5 (8^3x16) PASSED for all 10 channels (p0 class == fold
of `_t`); disc reader validated on a real `corr_d.p0` (39 configs x 25).
Still PENDING (needs production connected data): end-to-end 2D-C, jackknife,
plots; and `rel_norm` calibration (OQ1).
Prerequisite: a production run of `baryons_0000_dirac_claude` on real configs
(multi-GPU, cluster) to produce per-config HDF5 with the `<channel>_t_p<idx>`
datasets. The cold smoke-test h5 (8^3x16) only exercises the code path.

Disc file format (`corr_d.p{0..4}`, e.g.
`obs_nc4nf1_2448/<ens>/trace_<G>/corr_d.p0`): per config a `#File:` + `#Average
trace:` header then 25 lines `t  D(t)`, t=0..NT/2=24, real, already class-averaged
and folded. ~39 configs per file.

Class membership maps directly onto our 17 stored indices, and Rebbi's
`pcomp(m)/2` equals the stored-index count per class, so connected class-average
is a plain mean over the stored indices (real part):
- p0: {0}            (1)
- p1: {1,2,3}        (3 = pcomp1/2)
- p2: {4,5,6,7,8,9}  (6 = pcomp2/2)
- p3: {10,11,12,13}  (4 = pcomp3/2)
- p4: {14,15,16}     (3 = pcomp4/2)
Connected fold: C_fold(t)=(C(t)+C(NT-t))/2, t=1..NT/2-1; endpoints C(0), C(NT/2).
- Read `corr_d.p0..p4` (disc, already in Rebbi's class-averaged, folded,
  $t=0..NT/2$ form).
- Read the connected HDF5; fold $C(\vec p,t)$ as
  $C_\text{fold}(t)=(C(t)+C(NT-t))/2$ and average over each $p^2$ class
  (sum the stored components, divide by the same effective count Rebbi uses,
  `pcomp(m)/2`), keeping only the real part (connected data are real).
- Form $2D-C$ per class; plot / tabulate.

## Open questions

### OQ1 clarification (2026-06-18): the factor-3 is NOT a color factor
Rebbi's `makecorr_96_d.f90:110-115` says the disc trace data for his 96^3x192
runs were normalized a factor of 3 LOWER than his 64^3x128 runs, while the
connected `conn4d` used the 64^3 convention; the `pdata=3*pdata` rescales the
trace data so D and C are coherent. Both ensembles are N_c=3 QCD, so the 3 is a
ratio between two trace-file CODE VERSIONS, not a color trace:
- not 1/N_c (ratio of two N_c=3 ensembles would be 1, not 3);
- does NOT become 4 for SU(4): the relative D-vs-C weight is fixed by flavor/Wick
  combinatorics (the "2" in 2D-C, N_f-related); the color trace tr[Gamma S] sums
  N_c terms identically in D and C, so N_c cancels in the ratio.
For this SU(4) N_f=1 project both D (disc_multipleGamma_binary_claude) and C
(baryons_0000_dirac_claude) come from our own Grid codes, so there is no legacy
mismatch to import: rel_norm = 1 is the expected value. Residual risk is only a
projection/volume convention difference between the two Grid codes -> catch it by
the p=0 calibration below.

### hadron0 vs mesons_conn ~1.5 amplitude (RESOLVED 2026-06-18)
Cross-check of the connected p=0 correlator against the previously measured
`hadron/hadron0_2448_b11p045_m0p4000.npz` showed `ref ~= 1.5 x mesons_conn` at
small t. Investigation conclusion: it is NOT the driver and NOT thermalization.
- hadron0 is jackknife-NOT: variance test shows full per-config spread (would be
  ~100x smaller if resampled) -> bare per-config correlators (100 configs x 48 t).
- thermalization ruled out: mesons_conn C(0) is flat over the whole MC chain
  (corr with config index ~ -0.07; first/last quarter agree <1%).
- driver ruled out via git: hadron0 (Jan 17) predates commit b55ba91 (Feb 12);
  the only meson driver then was Mobius_mesons_xt.cc. Its meson sector
  (contraction trace(G5 adj(q1) G5 Gsnk q2 adj(Gsrc)), b=1.5/c=0.5, M5/mass args,
  PointSource kronecker=1.0, Import/ExportPhysical) is byte-identical to
  baryons_S2_dirac_parity.cc and baryons_0000_dirac.cc. The meson-region diff
  between the baryons_* files is empty. So the committed driver normalization
  never changed.
- signature points to a SOURCE overlap/normalization difference in the hadron0
  production: ground-state mass agrees to ~1% (source-independent), amplitude off
  by ~const ~1.5 (different source Z), and t-dependent ratio + slightly higher ref
  m_eff at mid-t (different excited-state contamination = different source shape).
  The hadron0 npz-builder + production source config live on the cluster, not in
  this tree (Untitled.ipynb only READS hadron0).
Consequences: previously measured MASSES are unaffected (normalization cancels in
m_eff). For 2D-C use mesons_conn (current driver) consistently with disc D; the
hadron0 offset is a separate-production artifact, not an OQ1 input.

### OQ1 — factor-3 trace normalization (UNRESOLVED, revisit later)
`makecorr_96_d.f90:231` does `pdata = 3*pdata`, with the in-code comment stating
this factor is to match the $96^3$ trace normalization to the $64^3$ one. It is
NOT obviously correct for this $24^3\times48$ ensemble. This factor sets the
RELATIVE normalization between $D$ and $C$ and therefore directly affects
$2D-C$. Decision deferred: for now treat the relative $C$-vs-$D$ constant as
to-be-calibrated, e.g. by matching the $p=0$ connected to a trusted reference
correlator, or by checking that $2D-C$ shows clean single-state behavior at
large $t$. Do not hard-code a factor until this is settled.

### OQ1/OQ2 normalization MATCHING derivation (2026-06-18)
Status: rel_norm has NOT been matched. All 2D-C numbers so far used rel_norm=1.0
(placeholder). Derivation of what rel_norm should be:

Disc as written to corr_d (average_trace2 + makecorr, lines cited):
- projection (average_trace2:143):  Stilde(p,t) = (1/NS^3) S(p,t),
  with S(p,t) = sum_x e^{ipx} T(x,t),  T = tr[Gamma S(x,x)].
- factor-3 (makecorr:231, ACTIVE):  pdata = 3 * Stilde.
- correlate (makecorr:240-243):  D_file(p,t) = (ns^3/NT) sum_t' pdata(t')pdata(t'+t),
  ns = NT/2 = 24 = NS.
Substituting (pdata = 3 Stilde = 3 S/NS^3):
  D_file(p,t) = (NS^3/NT) * 9 * (1/NS^6) sum_t' S(p,t') S(p,t'+t)
              = 9 * [ (1/(NS^3 NT)) sum_t' S(p,t') S(p,t'+t) ].
By translation invariance the bracket equals
  D_match(p,t) = sum_x e^{ipx} <T(x,t) T(0,0)>,
i.e. the point-source-at-origin / sum-at-sink disconnected correlator -- the SAME
convention as the connected C from a point source. Hence
  D_file = 9 * D_match,
and the physical singlet 2*D_match - C is
  2*D_match - C = 2*(D_file/9) - C   =>   rel_norm = 1/9
PROVIDED the two Grid codes normalize the single-quark object identically (see
caveat). The good news: Rebbi's 1/NS^3 and ns^3/NT factors are exactly what
reconstruct D_match, so structurally D and C are the same convention; only the
spurious factor-3 (-> 9) and the cross-code caveat remain.

Caveat (cross-code per-quark normalization): D comes from
disc_multipleGamma_binary_claude, C from baryons_0000_dirac_claude.

STEP (a) DONE 2026-06-18: read both codes. The cross-code constant is EXACTLY 1.
Both use:
- MobiusFermionD, same mass arg, M5=1.5, b=1.5, c=0.5;
- boundary {1,1,1,-1} (antiperiodic in time);
- ImportPhysicalFermionSource / ExportPhysicalFermionSolution (same physical
  single-quark normalization).
Disc loop estimator (disc_multipleGamma_binary_claude.cc):
StochasticDilutedSource uses Z4 noise nrm=1/sqrt(2) (E[xi xi*]=1) with FULL
spin/color dilution and full (t,eo) dilution; res[ig] += trace(Gamma psi adj(eta))
accumulated once per (t,eo) so every site is covered exactly once -> unbiased
T(x)=tr[Gamma S_phys(x,x)] with NO 1/Nhit or dilution factor (constant 1).
Connected uses an exact point source kronecker=1.0. Both C and D are second order
in the same S_phys. => cross-code constant = 1.

STEP (c) DONE 2026-06-18: validated at p0 on the LIGHT ensemble b10p865_m0p1000
with rel_norm=1/9. G5_G5 disc contribution (2D/9)/C grows 0.8% -> ~48% across
t=0..6 (eta'-like), m_eff(2D-C) departs upward from m_eff(C) -> disc has real
signal and rel_norm=1/9 is the right scale. I_I (scalar) needs VACUUM
SUBTRACTION (disc dominated by const condensate VEV ~240); separate from rel_norm.

STEP (b) DONE (tooling) 2026-06-18: centralized the Fortran into
obs_nc4nf1_2448/fortran_claude/ (average_trace2_claude.f90,
makecorr_96_d_claude.f90 with the factor-3 line commented out,
build_fortran_claude.sh, regenerate_corr_d_claude.sh, README_claude.md).
Verified: factor-3-free makecorr reproduces existing corr_d / 9 to ~1e-11 (all 5
classes). REGEN DONE + CONFIRMED 2026-06-18: user ran CONFIRM=1
./regenerate_corr_d_claude.sh. Verified 2*D_new(p0,t0)=3.908e-3 == earlier
2*(1/9)*D_old -> corr_d are now D_match. rel_norm = 1 everywhere
(meson_CminusD_claude.py default); the 1/9 / factor-3 is fully retired, no magic
factor remains.

Vacuum subtraction for I_I (scalar): LEFT OUT per user (2026-06-18). I_I p0 stays
condensate-VEV-dominated and is not a usable singlet without it; intentionally
not implemented.

CONCLUSION: rel_norm = 1/9 EXACTLY (only the spurious factor-3^2 from makecorr).
Apply rel_norm=1/9 in the Python combine (no makecorr rerun needed). Step (b)
(rerun makecorr without the factor-3 so corr_d files are correct at the source)
is optional cleanup. Step (c): validate at p0 on a LIGHT ensemble where D has
signal, confirming 2D-C single-state behavior with rel_norm=1/9.

At the heavy ensemble (b11p045 m=0.4) D is noise-dominated so rel_norm barely
affects 2D-C; it is decisive on the light ensembles where D has real signal.

### OQ2 — disc time-correlation normalization to replicate in Python
$D$ as written by `makecorr_96_d.f90` already includes $ns^3/NT$ (`:243`) and the
per-class $/(pcomp/2)$ (`:253`). The connected $C$ from a point source is the
standard $\sum_{\vec x}e^{i\vec p\cdot\vec x}\langle O(\vec x,t)O^\dagger(0)\rangle$
with no spatial-volume or source-time average. Confirm in Python that the two are
put on a common footing before subtracting (tied to OQ1).

### OQ4 — GZ channel definitions copy-paste bugs  [RESOLVED 2026-06-11]
Fixed in `baryons_0000_dirac_claude.cc:155-164`: sink set to GammaZGamma5 /
GammaZ (originals left commented). The p0==_t projection invariant is unaffected,
so no re-run needed for wiring; only the GZ/GZG5 correlator values change.
Original description:
In the ORIGINAL `baryons_0000_dirac.cc` the sink Gamma for two channels does not
match the source:
- `:115` `{Gamma::Algebra::GammaZGamma5, Gamma::Algebra::GammaYGamma5}` -- GZG5
  channel uses GammaYGamma5 at the sink (should be GammaZGamma5).
- `:119` `{Gamma::Algebra::GammaZ, Gamma::Algebra::GammaY}` -- GZ channel uses
  GammaY at the sink (should be GammaZ).
These propagate into `baryons_0000_dirac_claude.cc`. For the connected/disc
matching, `trace_gz` is tr[GammaZ S(x,x)] and `trace_gzg5` is tr[GammaZGamma5
S(x,x)], so the connected GZ_GZ / GZG5_GZG5 must use GammaZ / GammaZGamma5 at
BOTH ends or 2D-C is inconsistent for those classes. Ask user before changing
(do not silently "fix" the reference code).

### OQ3 — channel/directory mapping
Disc runs one $\Gamma$ per `trace_*/` dir (g5, gt, gtg5, gx..gz, gxg5..gzg5, id).
The connected channel list in `baryons_0000_dirac.cc:108-131` must use the SAME
$\Gamma$ at source and sink for each of these to be combined. Confirm the channel
names line up before forming $2D-C$.

## Algorithm source
Disconnected pipeline + momentum table + 2D-C convention: C. Rebbi, internal LSD
codes (`average_trace2.f90`, `makecorr_96_d.f90`), see this ensemble's `README`.

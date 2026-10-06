# Analysis of corr_data.tar.gz (LLNL copy) -- impl plan

## Goal
Analyze the meson correlators in `corr_data.tar.gz` pulled from LLNL
(`rsync_pull_corr_data_claude.sh` -> `/mnt/baracuda_14/grid_claude/corr_data_llnl_claude/extract/corr_h5_out_claude/`).
Six SU(4) $N_f=1$ ensembles, $24^3\times 48$, 10 channels, $p^2$ classes p0..p4,
contributions `C` (connected), `Dsub` (disc, per-config DC subtracted),
`CminusDsub` $= C - N_f D_\text{sub}$ (flavor singlet).

Methods (standard):
- Jackknife errors (M. H. Quenouille, Biometrika 43 (1956) 353; J. W. Tukey 1958),
  reusing `jackknife()` in `meson_combine_claude.py`.
- cosh effective mass on the folded correlator,
  $$m_\text{eff}(t) = \text{arccosh}\frac{C(t-1)+C(t+1)}{2C(t)},$$
  exact for a single $\cosh(m(t-T/2))$, independent of $T$.
- Plateau mass = uncorrelated weighted constant fit of $m_\text{eff}(t)$ over
  $[t_\text{min}, t_\text{max}]$, inside the jackknife.

## Step 1 (DONE 2026-10-06)
LLNL tarball is byte-identical to the June 23 local push (md5 93c8077e...),
all 6 h5 identical to local `corr_h5_out_claude/`. Layout = our own
`<contr>/<ch>/<mom>` (Nconf, 25) + `confs/<ch>`.
Configs per ensemble: b10.800/m0.01: 11, b10.840/m0.05: 24, b10.865/m0.1: 178,
b10.990/m0.2: 185, b11.035/m0.3: 339, b11.045/m0.4: 385.

## Step 2 -- analysis script
Files: `analyze_corr_data_llnl_claude.py`
- loop ensemble x channel x p-class x contr; jackknife mean/err of C(t) and $m_\text{eff}(t)$
- per-(contr, channel) default plateau windows (by eye, editable dict at top)
- outputs: `corr_llnl_analysis_claude/meff_<ens>.npz` + plateau summary table
  `corr_llnl_analysis_claude/plateau_summary_claude.txt`

## Step 3 -- plots
Files: `analyze_corr_data_llnl_claude.ipynb` (new; one plot per cell)
- C, Dsub, CminusDsub correlators (log |C|) for a chosen (ens, ch, psq)
- effective mass overlay C vs CminusDsub (marker + color distinct)
- plateau mass vs quark mass m across ensembles

## Step 2 status (2026-10-06): script DONE, run, outputs in
`corr_llnl_analysis_claude/` (npz per ensemble + `plateau_summary_claude.txt`).
Connected plateaus [12,20] clean, e.g. G5 p0: 0.5657(8), 0.8576(7), 1.1193(5),
1.3629(4) for m=0.1..0.4. Singlet `CminusDsub` [4,7] is weight-dominated by
t=4 -> NOT a plateau (excited states); signal dies at t~6-8.

## Open questions
- Fit windows: start by eye, user to refine.

- **Q1 (BLOCKING for z channels): GZ / GZG5 disc anomaly.** Cubic symmetry
  requires $D_{\gamma_x} = D_{\gamma_y} = D_{\gamma_z}$ statistically. Observed
  $D_\text{sub}(t=0)$, p0:

  | ens | GX | GY | GZ | GT | GXG5 | GZG5 |
  |---|---|---|---|---|---|---|
  | b10.800/m0.01 | 1.18e-2 | 1.20e-2 | 1.52e-3 | 1.28e-2 | 1.28e-2 | 1.77e-3 |
  | b10.865/m0.1  | 1.15e-2 | 1.12e-2 | 7.39e-4 | 1.22e-2 | 1.15e-2 | 9.59e-4 |
  | b11.045/m0.4  | 9.51e-3 | 9.32e-3 | 3.48e-4 | 9.57e-3 | 9.32e-3 | 3.92e-4 |

  GZ/GZG5 are 8-27x smaller AND mass-dependent; the others are a
  mass-independent ~1e-2 floor (looks like estimator noise). Already present in
  raw `trace_gz/averages/averages2.NNNN` (per-timeslice loop ~4x smaller than
  gx), so it is upstream of h5/makecorr: Grid trace output
  (`Grid/examples/disc_multipleGamma_binary_claude.cc:179-191`, Gamma list
  looks correct) or average_trace2's read of the complex trace.
  Hypothesis to test: $\text{tr}[\gamma_\mu S(x,x)]$ is purely imaginary for a
  $\gamma_5$-hermitian $S$
  ($\text{tr}[\gamma_\mu S]^* = \text{tr}[\gamma_5 S \gamma_5 \gamma_\mu] = -\text{tr}[\gamma_\mu S]$),
  so if average_trace2 keeps only the REAL part, the vector loops carry only
  noise, with a Gamma-matrix-dependent (Grid basis: $\gamma_x,\gamma_z$
  imaginary entries vs $\gamma_y,\gamma_t$ real) noise level. Which of
  GX/GY/GT vs GZ is "right" is then open.

  **Q1 RESOLVED (2026-10-06).** Read the raw SciDAC traces directly
  (`disc_loop_cache_claude.py`, data block at byte 1992, validated: Re-only
  projection == `averages2` to 1e-14, D == `corr_d.p0` to 5e-12).
  1. BUG: `average_trace2_claude.f90` keeps only Re tr (`tr(x,y,z,0) = a`).
     $\text{tr}[\gamma_\mu S]$ is imaginary, so ALL vector disc channels
     GX/GY/GZ/GT in corr_d / h5 are pure noise (zero signal). id, g5 (exactly
     real: Im slice sum ~1e-19) and $\gamma_\mu\gamma_5$ (real) use the right part.
  2. NOT a bug: the GZ/GZG5 "anomaly" is estimator variance. The spin-diagonal
     $Z_2\otimes Z_2$ noise in Grid's chiral basis makes the off-diagonal noise
     land differently in Re/Im depending on which spins $\Gamma$ connects:
     $\gamma_z,\gamma_t$ (spins 0-2, 1-3) have one low-noise part (~1e-5),
     $\gamma_x,\gamma_y$ (0-3, 1-2) have both parts noisy (~1e-4). The signal
     part is low-noise for gzg5 and gt, high-noise for gxg5, gyg5, gtg5, gx,
     gy, gz.
  3. Signal test b10.865 p0 (|mean/err| at t=1,2,3): g5 26/23/19, gzg5 13/6/0.6,
     all others unresolved. Correct D: real loops $+\langle r r\rangle$,
     vector loops $-\langle l l\rangle$ (`disc_corr_from_loops_claude.py`).
  => reliable disc: id, g5, gzg5 (stands for gx/gy g5 by cubic symmetry).
     Vector disc unresolvable (gt ~0 as expected at p=0).

- **Q3 (NEW): topology.** b10.865: total g5 loop per config $\propto Q/m$ has
  ensemble mean $-3.3(5)\times10^{-3}$ (7 $\sigma$ from 0); raw $D_{g5}(t)$
  has a minimum at $t\approx13$ and rises to $1.6\times10^{-4}$ at $T/2$.
  A finite-volume/poorly-sampled-topology constant in the eta' disc
  (cf. S. Aoki, H. Fukaya, S. Hashimoto, T. Onogi, PRD 76 (2007) 054508,
  arXiv:0707.0396). Subtraction choice for the final singlet numbers is OPEN.

- **Q4 (item 3, 2026-10-06): vector disc recomputed from the Im part**
  ($D_V = -\langle l\,l\rangle$, all 4 directions, p0..p4, all 6 ensembles,
  caches refreshed after the 2026-10-06 pull: disc 12/24/180/185/343/385).
  No direction / p-class shows a stable signal (all within ~2 $\sigma$);
  $D_{\gamma_t}$ = 0 to ~3e-6 (current conservation). Noise vs connected,
  spatial average, p0, matched configs:
  $\sigma_D/C$ at $t=2,3,4,5$ = 0.004/0.013/0.044/0.095 (b10.865),
  0.003/0.019/0.10/0.51 (b11.045); ~100% by $t\approx 6$, while the
  connected plateau starts at $t\approx 12$. So including $D_V$ kills the
  vector-singlet plateau; $D_V$ only bounds OZI violation:
  $|D_V/C| \lesssim 1$-$2\%$ at $t \le 3$.
  Options: (a) vector singlet = connected, quote the OZI bound;
  (b) improve the estimator in a new measurement (more hits / dilution
  adapted to the $\gamma_x,\gamma_y$ noise structure). Decision OPEN.

- **Q2: per-config DC subtraction biases the singlet.** Stored
  $$D_\text{sub}(t) = D(t) - \frac{1}{N_T}\sum_{t'} D(t')$$
  (all channels, all p-classes, `meson_CminusD_h5_claude.py:141`). This forces
  $\sum_t D_\text{sub}(t)=0$ per config, i.e. it also removes the ensemble
  signal $\frac{1}{N_T}\sum_{t'}\langle D(t')\rangle_\text{signal}$, giving a
  negative offset in $D_\text{sub}$ at intermediate $t$ (seen: $D_\text{sub}<0$
  from $t\approx 2$). Consequence: G5 singlet $m_\text{eff}$ drops BELOW the
  connected one (eta' lighter than the non-singlet, opposite of the
  anomaly expectation), I_I singlet goes above. For p>0 there is no VEV, so the
  subtraction is unmotivated there. Alternatives: (a) subtract ensemble-mean
  $V\,\langle\bar\psi\psi\rangle^2$ only for I_I p0 and nothing elsewhere;
  (b) fit with an explicit constant. Needs raw D (corr_d.p* locally), not in h5.

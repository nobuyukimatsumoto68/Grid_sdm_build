# EOFA HMC speedup: Hasenbusch + mixed-prec tuning + (optional) forecasting

## STATUS (2026-06-26)
- Chunk 0 (MX_inner 1000->10000): DONE -> `dweofa_mobius_HSDM_v5_claude.cc`. Local
  8^4 CUDA debug run validated (mixed-prec EOFA converges, dH O(1)). Cluster perf
  test PENDING.
- Chunks 1-2 (Hasenbusch k=2 + 3-level): DONE -> `dweofa_mobius_HSDM_hasenbusch_claude.cc`
  (built on v5, inherits MX_inner=10000). m1=sqrt(mass*pv)=0.1; LIGHT det[D(0.01)/D(0.1)]
  on Level1(mult1, coarsest), HEAVY det[D(0.1)/D(1.0)] on Level2(mult2), gauge on
  Level3(mult2, finest). Local debug: `run_local_8c_hasenbusch_claude.sh`. NOT yet
  built/run. Multipliers (1,2,2) are starting guesses, TUNABLE.
- Chunk 3 (forecasting): DEFERRED.

Target: `dweofa_mobius_HSDM_v4.cc` (SU(4) SDM, Mobius DWF EOFA HMC). Pain point:
m=0.01, 32c trajectories are slow. Scope agreed with user:
- ADD Hasenbusch mass preconditioning (new code).
- REVIEW/tune the EXISTING mixed-precision numbers (user flagged tol/iter caps as
  crucial).
- INCLUDE chronological forecasting -- but see Open Question Q3, it is NOT a
  drop-in for EOFA.
- Multi-timescale already exists (gauge Level2 / fermion Level1); reuse it.
- OUT OF SCOPE (user): Ls/ZMobius tuning, mixed-prec heatbath variant, standalone
  deflation (valence deflation eigvecs are tied to fixed configs, do not transfer
  to HMC).

## Physics / goal summary

The single-flavour determinant currently handled by ONE EOFA pseudofermion is
$$
\det\frac{D(m)}{D(\text{pv})}, \qquad m=0.01,\ \text{pv}=1.0 .
$$
At small $m$ the even/odd Schur operator $\hat H=$ `SchurDiagMooee` is
ill-conditioned. Valence bench (same op, $m=0.01$, $24^3\times48$, $L_s=16$):
$\lambda_\text{min}\approx4.8\times10^{-4}$, $\lambda_\text{max}\approx81.8$,
$\kappa\approx1.7\times10^5$, $\sqrt\kappa\approx410$ -> ~3000-iter CG baseline.
CG count scales like $1/m$, so the wins come from REDUCING the condition number
of each solve, not from faster mat-vec.

Hasenbusch: factor the determinant through intermediate masses
$m=m_0 < m_1 < \dots < m_{k}=\text{pv}$,
$$
\det\frac{D(m_0)}{D(m_k)} = \prod_{i=0}^{k-1}\det\frac{D(m_i)}{D(m_{i+1})} .
$$
Each factor $\det[D(m_i)/D(m_{i+1})]$ is a well-conditioned one-flavour object
handled by its own EOFA pseudofermion. The light factor ($m_0\to m_1$) still
carries the IR but with a much smaller ratio; the heavier factors are cheap.
Combined with the existing multi-timescale integrator, heavy/cheap factors ride
coarser MD timescales (smaller force, larger step) and the light factor rides the
finest fermion timescale.

## Algorithm sources (mandatory citations -- put in code comments too)

- Hasenbusch mass preconditioning: M. Hasenbusch, "Speeding up the Hybrid-Monte-
  Carlo algorithm for dynamical fermions," Phys. Lett. B519 (2001) 177,
  hep-lat/0107019.
- EOFA: Y.-C. Chen and T.-W. Chiu (TWQCD), "Exact Pseudofermion Action for Monte
  Carlo Simulation of Domain-Wall Fermion," Phys. Lett. B738 (2014) 55,
  arXiv:1403.1683. RBC/UKQCD / Grid implementation: D. J. Murphy.
- Mixed-precision reliable-update CG: M. A. Clark et al., Comput. Phys. Commun.
  181 (2010) 1517, arXiv:0911.3191.
- Chronological forecasting (initial-guess extrapolation): R. C. Brower, T.
  Ivanenko, A. R. Levi, K. N. Orginos, Nucl. Phys. B484 (1997) 353,
  hep-lat/9509012. (Grid: `Grid/algorithms/approx/Forecast.h`.)

## EOFA operator mass map (verified against Grid)

EOFA fermion op (`AbstractEOFAFermion.h:48`):
`D(mq1) + shift*\gamma_5 R_5 \Delta_\pm(mq2,mq3) P_\pm`.
Current code for $\det[D(m)/D(\text{pv})]$:
- `Op_L (mq1=mass, mq2=mass, mq3=pv, shift= 0.0, pm=-1)`
- `Op_R (mq1=pv,   mq2=mass, mq3=pv, shift=-1.0, pm=+1)`

Generalise factor $i$, $\det[D(m_i)/D(m_{i+1})]$, by replacing (mass,pv) ->
($m_i$, $m_{i+1}$):
- `Op_L_i (m_i,     m_i, m_{i+1}, 0.0, -1)`
- `Op_R_i (m_{i+1}, m_i, m_{i+1},-1.0, +1)`
(plus single-prec twins `Op_LF_i/Op_RF_i` for the mixed-prec CG, exactly as the
current `Op_LF/Op_RF`).

## Mixed-precision tuning review (user-flagged "crucial numbers")

Grid `MixedPrecisionConjugateGradient` (`ConjugateGradientMixedPrec.h`):
- Inner single-prec tol is SELF-RELAXING (line 126:
  `while(norm*inner_tol^2 < stop*1.01) inner_tol*=2`), so the inner tol knob is
  secondary.
- The governing knobs are `MaxInnerIterations` and `OuterLoopNormMult`.

Current HMC settings vs valence:
- `MX_inner = 1000` (`dweofa_mobius_HSDM_v4.cc:375`). Valence used `MAXINNER=10000`.
  At $m=0.01$ the first inner solve needs ~3000 iters, so `1000` TRUNCATES it and
  forces extra expensive double-prec outer restarts. PRIME suspect for slow
  mixed-prec. Action outer tol `1e-10`, force outer tol `1e-6`.
- `OuterLoopNormMult = 100.` hardcoded in the wrapper (`:144`); standard, leave.
- Action heatbath solve uses pure-double `ActionCG` (`:372`,`:416`); leave (out of
  scope: mixed-prec heatbath).

Valence reference record (the "other agent" conversation residue):
`Grid/examples/disc_mixedprec_impl_plan_claude.md`,
`Grid/examples/disc_mrhs_defl_bench_results_claude.md`.

## Files

- `dweofa_mobius_HSDM_v4.cc` -> NEW variant `dweofa_mobius_HSDM_hasenbusch_claude.cc`
  (keep v4 as the working reference for A/B). All changes are in this driver --
  Hasenbusch is just multiple EOFA actions; no Grid-core change needed UNLESS we
  do forecasting via Q3 option (a).
- `Makefile` / build target: add the new binary (user builds; Claude never
  compiles).
- XML parameter file: add a Hasenbusch mass list + per-factor timescale (Q1/Q2).

## Ordered implementation chunks

### Chunk 0 -- mixed-prec number tuning (cheap, standalone, do FIRST)
Files: `dweofa_mobius_HSDM_v4.cc` (in place is fine, single-line) OR fold into the
new variant.
- Bump `MX_inner` 1000 -> 10000 (match valence). This alone may give a large
  speedup with zero algorithmic change.
- (Optional) expose `MX_inner`, action/deriv outer tols via XML instead of
  hardcoded.
- Validate: one short run, compare iteration counts / wall vs current.

### Chunk 1 -- Hasenbusch ladder (core new code)
Files: `dweofa_mobius_HSDM_hasenbusch_claude.cc`
- Read a mass ladder `m_1 < ... < m_{k-1}` from XML (endpoints m and pv fixed).
- Loop $i=0..k-1$: build `Op_L_i/Op_R_i` (double) + `Op_LF_i/Op_RF_i` (single)
  with the mass map above; build the 4 mixed-prec CGs per factor (ActionCGL/R,
  DerivativeCGL/R) exactly mirroring the current single-factor wiring
  (`:382-412`); construct one `ExactOneFlavourRatioPseudoFermionAction` per factor.
- Store actions in a container (must outlive the run; reserve to avoid realloc /
  dangling refs -- EOFA holds references).
- `use_fc=true` per factor (heatbath forecasting over poles, as now).

### Chunk 2 -- integrator level assignment (reuse existing multi-timescale)
Files: same
- Put the light factor ($m_0\to m_1$) on the finest fermion level; heavier factors
  on a coarser fermion level (or share). Decide level counts (Q2).
- Gauge stays on its own level (currently Level2, 4 substeps).

### Chunk 3 -- chronological forecasting for MD solves (SEE Q3, may be deferred)
Files: TBD by Q3 (Grid-core EOFA edit, or driver-side `_claude` subclass).
- Add across-trajectory initial-guess forecasting to the action/force solves
  (`deriv` currently zeros the guess, `ExactOneFlavourRatio.h:434,446`).

### Chunk 4 -- build + validation
Files: build/run handoff via `tmp_claude.sh` teeing to `*_claude.log` (user runs).
- Reversibility / dH check, acceptance, and wall-time per trajectory vs v4
  baseline on one short stream.

## Open questions (RESOLVED 2026-06-26)

Q1. Mass ladder -> k=2, ONE intermediate $m_1=\sqrt{m\cdot\text{pv}}=\sqrt{0.01}=0.1$.
  Two factors:
  - Factor LIGHT: $\det[D(0.01)/D(0.1)]$  -- IR-sensitive, finest level.
  - Factor HEAVY: $\det[D(0.1)/D(1.0)]$   -- cheap, coarser level.
  (Hardcode for now; can promote to XML later.)

Q2. Integrator -> 3 levels:
  - Level1 (finest, most substeps): LIGHT factor.
  - Level2 (middle):                HEAVY factor.
  - Level3 (coarsest):              gauge (Wilson action).
  Ref `HMC/Mobius2p1f_DD_EOFA_96I_3level.cc`. Step ratios tunable; start from the
  current MD steps and split (e.g. gauge ~4x relative to its parent, as now).

Q3. Forecasting -> DEFERRED. Do Chunk 0 + Hasenbusch (Chunks 1-2), measure, revisit
  only if the LIGHT factor still dominates wall time. (Chunk 3 stays in this plan as
  future work; option (b) driver-side subclass preferred if revisited.)

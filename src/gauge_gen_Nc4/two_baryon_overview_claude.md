# Two-baryon correlator project: overview of the scratch / publication code

Goal of the new project: compute the correlator between **two-baryon states**
in the SU($N_c=4$) "stealth dark matter" (SDM) theory, in a naive way (point
source, single propagator, all contractions done explicitly).

This document lists the files in
`Grid_sdm_build/src/gauge_gen_Nc4/` that are relevant as a starting point and
describes what each does. The `.cc` files were copied from
`https://github.com/vmos1/su4_dm_grid_lsd/tree/main/run_gauge_gen/run_tioga/code`
(see `README`).

## Physics context (SU(4) baryon)

For gauge group SU($N_c$) the (single) baryon is a colour singlet built from
$N_c$ quarks antisymmetrised with the rank-$N_c$ Levi-Civita tensor
$\varepsilon_{a_1 \dots a_{N_c}}$. For $N_c=4$ this is a **4-quark** object:

$$
B \sim \varepsilon_{abcd}\, q_a q_b q_c q_d .
$$

The two-point function of $B$ is therefore a contraction of four quark
propagators with an $\varepsilon$ at the source and an $\varepsilon$ at the
sink. In the existing code these contractions are done by explicit nested loops
over the four colour indices $\{a,b,c,d\}$ at the sink and $\{a',b',c',d'\}$ at
the source, each carrying the permutation sign (parity) of the index tuple,
which reproduces $\varepsilon\,\varepsilon$. Spin structure is selected by
rotating the propagator into the Dirac-Pauli basis and picking specific
diagonal spin components.

The **two-baryon** correlator extends this: two such baryon interpolators at the
source and two at the sink, i.e. an 8-quark contraction. The single-baryon code
below is the natural template to generalise.

### Operator construction from the papers

`thermo_gw_details-5.pdf` Section 1.2 ("Baryon mass") and
`Lattice_Strong_Dynamics-19.pdf` End Matter (ii) fix the single-baryon
conventions that the code implements:

- For the SU($N_c=4$) **one-flavour** theory the flavour wavefunction is totally
  symmetric, so (local operator, no displacement) the spin must be fully
  symmetric, uniquely selecting a **total spin $S=2$** baryon.
- Spectroscopy at rest with definite parity is done in the **Pauli-Dirac (PD)**
  basis, where $\gamma_4 = \text{diag}(1,-1)$ is diagonal. Grid works in the
  **Weyl** basis, so the propagator is rotated
  $S^{(W)} = U S^{(PD)} U^\dagger$ with
  $U = (1/\sqrt2)(1 + \gamma_4^{(DP)}\gamma_5^{(DP)})$ — this is the
  `q_DP = 0.5 * U * q * Udagger` step in the code (`U` matches eq. 1.2.3, up to
  the $1/\sqrt2$ absorbed into the $0.5$).
- Dirac components in PD basis: `0` = positive parity / spin up, `1` = pos /
  down, `2` = negative parity / spin up, `3` = neg / down (eq. 1.2.7).
- **Parity-even $J^P=2^+$ operator** $B^{(0000)} = \varepsilon_{abcd}\,q_0^a q_0^b q_0^c q_0^d$
  (eq. 1.2.13). All $4!=24$ permutations give the identical term, so
  $$
  C^{(0000)} = 24\,\varepsilon_{abcd}\varepsilon_{a'b'c'd'}\,
  S_{00}^{aa'} S_{00}^{bb'} S_{00}^{cc'} S_{00}^{dd'}
  $$
  (eq. 1.2.17) — exactly the weight-24 contraction in `baryons_S2_dirac_parity.cc:261`.
- **Parity-odd $J^P=2^-$ operator** $B^{(0002)} = \sqrt4\,\varepsilon_{abcd}\,q_0^a q_0^b q_0^c q_2^d$
  (eq. 1.2.19), giving
  $$
  C^{(0002)} = 4\,\varepsilon_{abcd}\varepsilon_{a'b'c'd'}
  [\,6\,S_{00}^{aa'}S_{00}^{bb'}S_{00}^{cc'}S_{22}^{dd'}
  + 18\,S_{00}^{aa'}S_{00}^{bb'}S_{02}^{cc'}S_{20}^{dd'}\,]
  $$
  (eq. 1.2.22): 6 direct ($q_2\!\to\!q_2$) + 18 cross ($q_2\!\to\!q_0$ via
  $S_{02},S_{20}$) contractions — exactly the weight-6 and weight-18 terms in
  `baryons_S2_dirac_parity.cc:319,321`.

So the existing publication code computes the **single-baryon** $2^+$ and $2^-$
masses (the $M_B$ used to set the line of constant physics in the LSD paper).
The two-baryon correlator $\langle (BB)(t)\,\overline{(BB)}(t_0)\rangle$ is a
genuinely new, 8-quark contraction that does **not** yet exist in this code base
(see search note below).

## Primary reference files (single baryon + meson spectroscopy)

These are the publication-style spectrum codes. Both share an identical
`main()`, `Solve()`, `PointSource()`, and `MesonTrace_hdf()`; they differ only
in `BaryonTrace_hdf()`.

### `baryons_S2_dirac_parity.cc` (most complete single-baryon reference)
- `main()` (line 353): loads a NERSC gauge config (`config M5 mass outfile`
  command-line args), builds a `MobiusFermionD` domain-wall action ($L_s=16$,
  scaled Shamir kernel $b=1.5$, $c=0.5$), anti-periodic temporal BC
  `{1,1,1,-1}`. Default with no args = cold config, $M_5=1.5$, $m=0.1$.
- `PointSource()` (line 18): point source = Kronecker `SpinColourMatrix` at a
  site, poked into a `LatticePropagator`.
- `Solve()` (line 28): inverts the Dirac operator for all 12 spin-colour source
  components via `SchurRedBlackDiagMooeeSolve` + `ConjugateGradient`
  (tol 1e-9), building the full point-to-all `LatticePropagator`.
- `MesonTrace_hdf()` (line 93): 10 meson channels (I, $\gamma_5$,
  $\gamma_i\gamma_5$, $\gamma_\mu$), `meson_CF = trace(G5 adj(q1) G5 Gsnk q2 adj(Gsrc))`,
  `sliceSum` along both $X$ and $T$, written to HDF5.
- `BaryonTrace_hdf()` (line 190): the **key routine**. Rotates the propagator to
  Dirac-Pauli basis `q_DP = 0.5 * U * q * Udagger` (U defined at line 206),
  peeks spin components `q00,q02,q20,q22` (line 220), then does the
  $\varepsilon\varepsilon$ colour contraction over 4 sink + 4 source colours
  with parity signs. Two channels:
  - `ch=0` ("bar_0000"): all four quarks use spin component `q00`
    (weight 24, line 261).
  - `ch=1` ("bar_0002"): mixed spin components `q00,q22` (weight 6) and
    `q00,q02,q20` (weight 18) (lines 319, 321) — different spin coupling of the
    4-quark operator.
  Result `sliceSum` along $T$, written to HDF5.

### `baryons_0000_dirac.cc` (simpler single-baryon variant)
- Same structure, but `BaryonTrace_hdf()` (line 190) loops over a channel array
  `channel_spinidx={0,2}` ("bar_0000","bar_2222"), each channel using a single
  diagonal spin component `q_{ii}` for all four quarks (weight +1, line 254).
  This is the cleanest, most readable version of the $\varepsilon\varepsilon$
  4-quark contraction and is the best starting template.

## Test / run scaffolding

### `test/run_baryon/`
- `run.sh`: runs the baryon binary on one config. Sets `--grid/--mpi/--threads`,
  config path, `M5`, `mass`, output `outfile.h5`. Currently points at
  `../bin/baryons_S2_dirac_parity` with no config args (cold-config smoke test).
- `run_old.sh`: same but launches a reference ("old") binary via `flux run` on a
  24.24.24.8 lattice for regression comparison.
- `compare.py`: reads two HDF5 files and `np.allclose`-compares every
  `<channel>/data` dataset (re+im) — the regression check against the old code.

### `test/Makefile`
Builds every `*.cc` in `test/` against the installed Grid
(`../../install/Grid_omp_Nc4`) via `grid-config`, output to `test/bin/`.

## Build

### `Makefile` (top level of `gauge_gen_Nc4/`)
Builds each `*.cc` into `./bin/` using `grid-config` from
`../../install/Grid_omp_Nc4`; `make install` moves binaries to
`../../install/gauge_gen_Nc4`.

## Other files in the directory (context, not directly baryon-related)

- `Mobius_mesons_xt.cc`: fuller meson spectroscopy (includes a covariant
  Laplacian / smearing helper `CovariantLaplacianCshift`).
- `eye4_anti.cc`: Wilson-fermion utility (RNG sources); name suggests an
  eye/disconnected or identity test for the 4-quark setup.
- `pbp.cc`, `pbp_list.cc`: $\bar\psi\psi$ (chiral condensate) measurements.
- `glueball.cc`, `glueball_imp.cc`: glueball operators.
- `WilsonFlow.cc`, `getflowobs_v3.cc`: Wilson-flow scale setting / flowed
  observables.
- `dweofa_mobius_HSDM_v*.cc`, `checkBC_MobiusEOFA.cc`, `hmc_WilsonGauge.cc`,
  `reweight2.h`: HMC gauge generation (EOFA Mobius DWF) and reweighting — the
  configuration-generation side, not measurement.
- `freeze/`: older frozen versions of the HMC / DWF drivers.

## Suggested starting point for the two-baryon code

Generalise `BaryonTrace_hdf()` from `baryons_0000_dirac.cc`:
- Build the single point-to-all propagator `q` exactly as now.
- For two baryons at rest, the simplest naive construction contracts two
  $\varepsilon_{abcd}$ blocks at the source and two at the sink (8 quark lines),
  with the appropriate spin projection per baryon. Open question to settle with
  the user: which interpolating operators / spin-isospin channel, momentum
  projection (zero vs back-to-back), and whether both baryons share the same
  point source or use displaced/wall sources.

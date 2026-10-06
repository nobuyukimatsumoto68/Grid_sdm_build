# Momentum-projected single-baryon correlators from A. Meyer's cubic-irrep operators

Status: PLAN, awaiting go-ahead (2026-10-03). Session grid-claude-02; peer grid-claude-54 owns the
local bvar analysis side (`rsync_pull_bvar_claude.sh`, `bvar_h5_claude.py`, `bvar_spectrum_claude.py`)
and the production files `baryon_variants_corr_claude.cc` / `baryon_pd_common_claude.h` (NOT edited).

## Sources (algorithms / operators)

- Operators: A. S. Meyer (LLNL), 4-quark local operators projected onto $O_h$ irreps at rest and
  onto little-group irreps in flight, total momentum up to $(2,1,1)$, private communication
  2026-09-28 (LC give/take, `baryon_group_theory/aaron_4q_ops/`). Construction in the spirit of
  Basak et al., hep-lat/0506029, and Thomas, Edwards, Dudek, arXiv:1107.1930 (in-flight operators);
  little-group conventions Goeckeler et al., arXiv:1206.4141.
- Contraction: unchanged determinant form of `baryon_variants_impl_plan_claude.md` Sec. 2
  (Detmold-Savage-style composite-index determinants, arXiv:1001.2768).
- Momentum projection by spatial FFT of each site density (Grid FFT, cuFFT/hipFFT backend).

## 1. Physics / goal

Add the single-baryon spectrum in moving frames, i.e. the lattice dispersion $E_B(\mathbf P)$ per
little-group irrep, and the finer rest-frame split $E^\pm$ / $T_2^\pm$ of our $J=2$ channels, by
changing only the contraction (no new propagator dumps). Scope: the light ensembles of the
point-source update (m0.1, m0.05, m0.01).

Verified so far (`baryon_group_theory/compare_aaron_ops_claude.py`): Aaron's single-flavour (`i2c2`)
operators are local ($B_A(x)$ with one momentum per operator, all quarks at $y0$), his projectors give our
PD basis ($(+,\uparrow)=0$, $(+,\downarrow)=1$, $(-,\uparrow)=2$, $(-,\downarrow)=3$), and every rest-frame
operator lies inside exactly one of our 9 $(J^P,n_+,n_-)$ operators.

Correlator. With the point source at $y$ (the dump's `srcCoord`, origin in production) and sink momentum
$\mathbf p$,

$$
C_{AA'}(\mathbf p,t) = \sum_{\mathbf x} e^{-i\mathbf p\cdot(\mathbf x-\mathbf y)}\,
  \langle B_A(\mathbf x,t)\, B^\dagger_{A'}(y)\rangle ,
$$

i.e. the existing site density $d_{AA'}(x)$ (`PDTermTable::ContractPair`) Fourier transformed in space
instead of `sliceSum`. For Aaron's operators $O_i(\mathbf p)=\sum_A c^i_A(\mathbf p) B_A(\mathbf p)$ of one
(frame, irrep, component) group, all at the same $\mathbf p$,

$$
C_{ij}(\mathbf p,t) = \sum_{A,A'} c^i_A(\mathbf p)\, \overline{c^j_{A'}(\mathbf p)}\, C_{AA'}(\mathbf p,t).
$$

The source operator is $O_j^\dagger$ built from the conjugated `.op` coefficients. Aaron's `.cc` files
(charge-conjugate creation operators, `DOWNBART`) are NOT used; they belong to his two-flavour contraction
code. Phase convention: Aaron's `MOM` is taken as the sink momentum with $e^{-i\mathbf p\cdot\mathbf x}$
(Grid `FFT::forward`). The opposite convention would give his operator at $-\mathbf p$, which is still a
valid operator of the same little-group irrep, and all directions of the star are computed anyway.

Cost drivers:
- In flight, every one of the $35^2=1225$ elementary pairs is needed, including parity-mixing pairs
  (parity is not a symmetry in flight). That is 65536 det terms per site in total, vs 3024 for the
  production 42-pair block mode: about 22x the kernel work.
- 81 momenta (frames 000, 100, 110, 111, 200, 210, 211 with all directions). One 3D FFT per pair per
  sink type, then a host copy of the transformed density, picking the 81 x T values.

Output sizes per config per sink type: component-averaged operator correlators 2111 x T (about 1.6 MB);
raw per-component 39193 x T (about 30 MB); all elementary $C_{AA'}(\mathbf p)$ 1225 x 81 x T (about 76 MB).

## 2. Files

New (none of the production files is edited):
- `Grid_sdm_build/src/gauge_gen_Nc4/aaron_ops_convert_claude.py`: Aaron `i2c2/*.op` ->
  `baryon_mom_ops_claude.txt`; also writes `baryon_variants_terms_all_claude.txt` (all 1225 pairs) by
  importing `pair_terms` / `write_terms` from `baryon_variants_terms_gen_claude.py` (read-only use).
- `Grid_sdm_build/src/gauge_gen_Nc4/baryon_mom_ops_claude.txt`, `baryon_variants_terms_all_claude.txt`
  (generated tables).
- `Grid_sdm_build/src/gauge_gen_Nc4/baryon_variants_mom_corr_claude.cc`: copy of
  `baryon_variants_corr_claude.cc` plus the momentum path (class `PDMomOperatorTable` in the .cc).
- `Grid_sdm_build/src/gauge_gen_Nc4/run_bvar_mom_contract_flux_claude.sh` (remote, skip-if-exists) and
  `bvar_mom_contract_handoff_claude.md`.
- Local test script `tmp_claude.sh` + log (user runs; GPU binary).

## 3. Chunks

### Chunk 1: operator + term tables (Python)
Files: `aaron_ops_convert_claude.py`, `baryon_mom_ops_claude.txt`, `baryon_variants_terms_all_claude.txt`
- Parse every `i2c2/pcom*/qqqq*/*.op`; key = (frame, irrep, component); within a group the ops are the
  unique IDs; each group has one momentum (checked: true for all 265 groups).
- Table format (text, C++-friendly):
  `nmom` + `mom <k> <px> <py> <pz>`; `ngroup` + `group <frame> <irrep> <comp> <k> <nops>`, then per op
  `op <uid> <ncomp>` + `<A> <re> <im>` lines (A = sorted PD tuple, coefficients merged per sorted A).
- Self-checks: rest-frame ops reproduce `compare_aaron_ops_claude.py`; all pairs needed by the groups
  exist in the terms file.

### Chunk 2: new contraction binary
Files: `baryon_variants_mom_corr_claude.cc`
- Same CLI as production plus `--momops <file>` (enables the momentum path), `--mom-raw` (also store
  per-component correlators) and `--no-p0` (omit the legacy p=0 `Cel_`/`Cop_` datasets; set by the run
  script when the production `bvar.all.<idx>.h5` already exists).
- With `--momops`: compute all pairs the groups need (all 1225, taken from `--terms` = the all-pairs file).
  For each pair and sink type: density -> `sliceSum` (legacy p=0 `Cel`/`Cop` for the `--ops`/`--pairs block`
  set, unchanged dataset names, so the file is a superset of the production output) -> spatial FFT ->
  host copy -> accumulate $C_{AA'}(\mathbf p_k,t)$ for the 81 momenta (source phase $e^{+i\mathbf p\cdot\mathbf y}$).
  The host copy requires a single MPI rank (the production contraction already runs `--mpi 1.1.1.1`).
- Operator correlators $C_{ij}$ per group; averaged over the components of each (frame, irrep).
- HDF5: `Cmom_<frame>_<irrep>_<uid_i>_<uid_j>` (+ `Cmomss_` smeared sink), meta extended with
  frames, irreps, op IDs and the source table name; with `--mom-raw` also `Cmomraw_<frame>_<irrep>_c<comp>_...`.

### Chunk 3: local validation (user runs `tmp_claude.sh`)
Files: `tmp_claude.sh`, `chunk_mom_validate_claude.log`, small check script
- Cold `ColdConfig.pdfull.lime`: legacy `Cel_*`/`Cop_*` equal the production `ColdConfig.bvar.h5` to
  machine precision; FFT at $\mathbf p=0$ equals `sliceSum`; rest-frame `a1p`, `t1p`, `eep`/`t2p` equal the
  matching combinations of our `Cop_*`; all components of one (frame, irrep) agree exactly (the cold
  field with the source at the origin is cubic-symmetric), which also tests Aaron's relative phases
  before we trust the component average; hermiticity $C_{ij}=\overline{C_{ji}}$.
- Real gauge `lat16c.pdfull.lime` ($16^3\times8$): FFT path vs explicit phase x `sliceSum` at a few momenta.
- Timing of the 1225-pair pass on the TITAN V (scaled estimate for tuolumne).

### Chunk 4: remote run script + handoff
Files: `run_bvar_mom_contract_flux_claude.sh`, `bvar_mom_contract_handoff_claude.md`
- Fan-out over configs, `HDF5_USE_FILE_LOCKING=FALSE`, skip a config when its output file exists.
- Output dirs `obs_bvar_2448_b<beta>_m<mass>[_point]_ss_mom/bvarmom.all.<idx>.h5` (distinct suffix, agreed
  with grid-claude-54, who extends the aggregator regex). Dataset names sent to grid-claude-54.

## 4. Decisions (user, 2026-10-03)

1. Skip test = OLD p=0 output: if the production `bvar.all.<idx>.h5` exists in the matching
   `obs_bvar_2448_..._ss/` (or `_point_ss/`) dir, the run script passes `--no-p0` and the binary writes
   only the momentum datasets; otherwise it also writes the p=0 `Cel_`/`Cop_` datasets. In addition, a
   config whose own `bvarmom.all.<idx>.h5` exists is skipped (restartability).
2. Scope: every new light dump, i.e. m0.1 point, m0.05 and m0.01 with both point and smeared source.
3. Storage: component-averaged `Cmom_*` by default; per-component `Cmomraw_*` only with `--mom-raw`.

## 5. Chunk log

- Chunk 1 DONE (2026-10-03): `aaron_ops_convert_claude.py` -> `baryon_mom_ops_claude.txt` (2835 .op files,
  265 groups, 81 momenta; frames 000/100/110/111/200/210/211 = 17/36/48/32/36/48/48 groups) and
  `baryon_variants_terms_all_claude.txt` (1225 pairs, 65536 terms). Self-checks passed: one momentum per
  group, every ID present in every component of its (frame, irrep), every op at fixed $(n_+,n_-)$.
- Chunk 2 WRITTEN, NOT YET COMPILED (2026-10-03): `baryon_variants_mom_corr_claude.cc` (copy of production +
  `MomMeta`, `PDMomOperatorTable`, `MomProjector` (spatial FFT, host lex copy, single rank),
  `ContractAllMom` (one density per pair: legacy sliceSum + FFT projection, built-in FFT(p=0) vs sliceSum
  check), `PhaseCheck` (explicit phase x sliceSum), `GroupCorrelators`; CLI `--momops --mom-raw --no-p0
  --phase-check`; changed production calls kept commented as `// ORIGINAL:`). Datasets: `mom_meta`,
  `Cmom_<frame>_<irrep>_<uid_i>_<uid_j>`, `Cmomss_*`, optional `Cmomraw_<frame>_<irrep>_c<comp>_*`.
- Chunk 3 SCRIPTED (2026-10-03): `Grid_sdm_build/tmp_claude.sh` -> `chunk_mom_validate_claude.log` (compile; cold 8^4
  full + `--no-p0`; 16c real gauge with sink smearing and timing), checks in `baryon_variants_mom_check_claude.py`.
  Rest-frame relation tested on the cold field (cubic symmetry, positive overlaps from
  `compare_aaron_ops_claude.py`, $n$ = ket norm$^2$ $\sum_A |c_A|^2/k'_A$ of Aaron's operator):

$$
C^{\rm ours}_{XY} = \frac{C^{\Gamma}_{XY}}{\sqrt{n_X n_Y}} \quad (J=0,1),
\qquad
C^{\rm ours}_{XY} = \frac12\left[\frac{C^{E}_{XY}}{\sqrt{n^E_X n^E_Y}} + \frac{C^{T_2}_{XY}}{\sqrt{n^{T_2}_X n^{T_2}_Y}}\right] \quad (J=2),
$$

  since our $M=J$ state is a $T_1$ row ($J=1$) or $|2,2\rangle=(E_1+T_{2,z})/\sqrt2$ ($J=2$; the cross terms vanish by
  cubic symmetry). Also: legacy datasets vs production files, every component equal to its class average
  (cold), hermiticity, `--no-p0` consistency, FFT vs explicit phase (`--phase-check`), FFT$(p=0)$ vs `sliceSum`.
- Chunk 3 DONE + VALIDATED (2026-10-06, `Grid_sdm_build/chunk_mom_validate_claude.log`; local toolchain rebuilt first:
  Grid `build_cuda129` + Open MPI `openmpi_cuda129`, see memory local_toolchain_cuda129_gotcha). Cold 8^4 (2 corners):
  legacy 122 datasets identical to production (0.0); FFT vs explicit phase 5e-16; FFT(p=0) vs sliceSum 5e-16; rest-frame
  relation 9e-14 (point) / 5e-14 (smeared); all 78386 per-component correlators equal their class average 7e-13 (Aaron's
  conventions = our PD basis, also in flight); hermiticity 9e-16; `--no-p0` Cmom identical (0.0) with no legacy datasets.
  16c real gauge: legacy 61 datasets identical (0.0); phase check 9e-16; FFT(p=0) 1.5e-15. Cold E_eff(t=2) rises with |P|
  (3.31 rest -> 3.36/3.41/3.45/3.44/3.48/3.52). Fixes during the chunk: check metric now scaled by the pair's max over
  (p,t) (per-entry scale was meaningless for parity-vanishing p=0 entries); empty legacy meta strings -> "none" (HDF5
  rejects zero-length string attributes). Timing (TITAN V, per sink pass, 1225 pairs): 8^4 contraction 0.33 s + FFT 7.5 s;
  16c contraction 1.37 s + FFT 7.1 s (FFT cost is per-call overhead, ~6 ms/pair).
- Chunk 4 DONE (2026-10-06): `test/run_bvar_mom_contract/run_bvar_mom_contract_flux_claude.sh` (mirrors the remote production
  `test/run_bvar_contract/run_bvar_contract_flux_claude.sh`, pasted by the user: NSLOTS fan-out, `flux run -n1 -g1`, gauge
  `${CONFROOT}/${cfgfilename}/${cfgfilename}_lat.<idx>`, dumps `bvar_2448_<ens>[_point]/<cfgfilename>.<idx>.pdfull.lime`;
  ENS adds b10p840_m0p0500 / b10p800_m0p0100; skip if own output exists, `--no-p0` if production `bvar.all.<idx>.h5` exists,
  `.inprogress.h5` -> rename; PHASECHECK on first config) + `bvar_mom_contract_handoff_claude.md`. Production uses
  `--pairs sector`; on the all-pairs term file it selects the same 109 pairs (checked). m0.1 SMEARED dumps still exist on
  lustre5 -> first validation/production target (all `--no-p0`). Names sent to grid-claude-54.

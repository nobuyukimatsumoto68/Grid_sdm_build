# Single-baryon correlators of the 9 B-form highest-weight operators: implementation plan

## Sources / citations

- Operator definition, PD basis, $J^P$ decomposition and the highest-weight table:
  `baryon_group_theory/baryon_theory_claude.md` Sec. 4.0 (lines 385-417) and Sec. 4.1 (lines 421-511),
  the user's own note.
- PD rotation of Grid's chiral-basis propagator, $G^{PD}=\tfrac12 U G^W U^T$, and the pre-rotated
  source trick: `point2all_prop_dumper_claude.cc:106-118` (`WeylToPauliDiracU`), `:58-84`
  (pre-rotated smeared source), `:346-381` (sink projection), `:396-436` (`--rottest`).
- Determinant form of identical-quark contractions (all quark pairings of one flavour collapse
  into a determinant over composite spin-colour indices): W. Detmold and M. J. Savage,
  arXiv:1001.2768; T. Doi and M. G. Endres, arXiv:1205.0585 (unified contraction algorithm).
  The single-baryon case used here is the $4\times4$ special case, written out in Sec. 2 below.
- $576\det q_{00}$ normalisation of the $0000$ correlator: `two_baryon_corr_prod_claude.cc:175-197`.
- Gaussian smearing: `Grid/qcd/utils/CovariantSmearing.h:41` (Chroma convention).
- Split-grid MRHS batching: `Grid/tests/solver/Test_dwf_mrhs_cg.cc`, as used by the dumper.

## 1. Goal

SU(4), $N_f=1$ SDM ensembles ($24^3\times48$, $m\in\{0.1,0.2,0.3,0.4\}$). Compute the zero-momentum
two-point functions of the nine highest-weight local baryon operators (one per $(J^P,n_+,n_-)$
sector), including the cross-correlators between operators of the same $J^P$, so that

- the published $2^\pm$ ground-state masses are reproduced,
- the same levels appear in every channel with the same quantum numbers
  ($2^+$: (4,0), (2,2), (0,4); $2^-$: (3,1), (1,3); $1^-$: (3,1), (1,3)).

Only the highest-weight states are in scope. No correlators between different quantum numbers on
real configurations.

## 2. Physics

### 2.1 Operators

PD Dirac index: 0 = upper, $s_z=+\tfrac12$, $P=+$; 1 = upper, $s_z=-\tfrac12$, $P=+$;
2 = lower, $s_z=+\tfrac12$, $P=-$; 3 = lower, $s_z=-\tfrac12$, $P=-$. With $A:=(\alpha,\beta,\gamma,\delta)$ a sorted 4-tuple of PD
indices (written as digits, e.g. `0002`; $A_i$ denotes its $i$-th entry),

$$
B_A(x)=\epsilon_{abcd}\,q^a_{\alpha}(x)\,q^b_{\beta}(x)\,q^c_{\gamma}(x)\,q^d_{\delta}(x),
\qquad P(A)=(-1)^{n_-(A)},
$$

totally symmetric in the four indices (35 distinct sorted tuples $A$). Highest-weight operators (Sec. 4.1 of the
note, verbatim):

| name | $J^P$ | $(n_+,n_-)$ | $B$ form |
|---|---|---|---|
| `2p_40` | $2^+$ | (4,0) | $B_{0000}$ |
| `2m_31` | $2^-$ | (3,1) | $B_{0002}$ |
| `1m_31` | $1^-$ | (3,1) | $\tfrac{\sqrt3}{2}(B_{0003}-B_{0012})$ |
| `2p_22` | $2^+$ | (2,2) | $B_{0022}$ |
| `1p_22` | $1^+$ | (2,2) | $B_{0023}-B_{0122}$ |
| `0p_22` | $0^+$ | (2,2) | $\tfrac1{\sqrt3}(B_{0033}-2B_{0123}+B_{1122})$ |
| `2m_13` | $2^-$ | (1,3) | $B_{0222}$ |
| `1m_13` | $1^-$ | (1,3) | $\tfrac{\sqrt3}{2}(B_{0223}-B_{1222})$ |
| `2p_04` | $2^+$ | (0,4) | $B_{2222}$ |

Same-$J^P$ blocks: $2^+$ {`2p_40`,`2p_22`,`2p_04`}, $2^-$ {`2m_31`,`2m_13`}, $1^-$ {`1m_31`,`1m_13`},
$1^+$ {`1p_22`}, $0^+$ {`0p_22`}.

### 2.2 Conjugate operator

Osterwalder-Schrader reflection: $\Theta q=\bar q\gamma_4$, antilinear, order reversing;
$\gamma_4^{PD}=\mathrm{diag}(1,1,-1,-1)$, so $(\gamma_4)_{\alpha\alpha}=\eta_\alpha$ with
$\eta=(+,+,-,-)$ and $\prod_i\eta_{A_i}=P(A)$:

$$
B^\dagger_{A'}(0)\equiv\Theta B_{A'}(0)
=P(A')\,\epsilon_{abcd}\,\bar q^d_{\delta'}\,\bar q^c_{\gamma'}\,\bar q^b_{\beta'}\,\bar q^a_{\alpha'}(0).
$$

For an operator $O=\sum_A c_A B_A$ the reflection-positive correlator is
$C_{OO'}(t)=\sum_{AA'}c^O_A\,\bar c^{O'}_{A'}\,C_{AA'}(t)$, with
$C_{AA'}(t)=\sum_{\vec x}\langle B_A(\vec x,t)\,B^\dagger_{A'}(0)\rangle$.

### 2.3 Wick contraction (one flavour)

Write $q_i=q^{a_i}_{A_i}$, $\bar q_{j}=\bar q^{b_j}_{A'_j}$ and
$G^{aa'}_{\alpha \alpha'}(x,0)=\langle q^a_\alpha(x)\bar q^{a'}_{\alpha'}(0)\rangle$ (the
$16\times16$ spin-colour site propagator; colour up, spin down). The ordering
$q_1q_2q_3q_4\,\bar q_4\bar q_3\bar q_2\bar q_1$ gives Grassmann sign $+1$ for the identity pairing
and $\mathrm{sgn}(\pi)$ for the pairing $q_i\leftrightarrow\bar q_{\pi(i)}$. Relabelling the source
colours by $\pi$ turns $\epsilon_{a'}$ into $\mathrm{sgn}(\pi)\epsilon_b$, so the two signs cancel:
$$
\langle B_A\,B^\dagger_{A'}\rangle
=P(A')\sum_{\pi\in S_4}\ \sum_{a,b\in S_4}\epsilon_a\,\epsilon_b\ \prod_{i=1}^4
G^{a_i\,b_i}_{A_i\,A'_{\pi(i)}} .
\tag{naive form}
$$

Substituting $b\to b\circ\pi$ inside the $b$ sum turns the pairing sum into a determinant over
composite indices:
$$
\langle B_A\,B^\dagger_{A'}\rangle
=P(A')\sum_{a,b\in S_4}\epsilon_a\,\epsilon_b\ \det_{ij}\Big[G^{a_i\,b_j}_{A_i\,A'_j}\Big].
\tag{det form}
$$

The naive form is used by the Python checker and by the driver's `--naive-check`; the det form is
what the term tables encode.

### 2.4 Term tables

Viewing $G$ as a $16\times16$ matrix with composite index $4\times\text{spin}+\text{colour}$, the rows
of the $4\times4$ minor are $r_i=4A_i+a_i$ and the columns $c_j=4A'_j+b_j$.
Sorting the rows introduces the sign of the sorting permutation, so each colour ordering $a$
contributes $\epsilon_a\,\mathrm{sgn}(\mathrm{sort}_a)$ to the sorted key. Orderings that only
permute slots with equal spin index give the same sorted key and the same sign, so the
$4!$ orderings collapse to $4!/\prod_j k_j!$ distinct keys with integer weight $\prod_j k_j!$
($k_j$ = multiplicity of index $j$ in $A$). The table for a pair is the list of

$$
\big(w,\ (r_1,r_2,r_3,r_4),\ (c_1,c_2,c_3,c_4)\big),\qquad
w=P(A')\,w_R\,w_C,\qquad \sum|w|=576 .
$$

Checks: `0000 0000` has one term with $w=576$ (so $C=576\det q_{00}$); `0123 0123` has $576$ terms
with $w=\pm1$.

### 2.5 Selection rules used only for validation

After averaging, $C_{AA'}=0$ unless $P(A)=P(A')$ and $J_3(A)\equiv J_3(A')\pmod 4$. On the cold
(free) configuration with a point source these hold exactly, which validates the coefficient signs
of the operator table (Chunk 5).

## 3. Pipeline

1. **Dump** (`point2all_prop_dumper_full_claude.cc`): per config, for each of `--nsrc` corners
   (default 1 = origin) solve the 16 pre-rotated PD sources (spin $s'$ = row $s'$ of $U$, colour
   $c'$), project the sink to the four PD spins, store 16 `LatticeColourMatrixF` records per corner
   tagged `store="q<s><s'>"`. Output `<config>.pdfull.lime`.
2. **Contract** (`baryon_variants_corr_claude.cc`): read the 16 blocks into one `LatticePropagator`
   per corner, evaluate the term tables of the selected operator's $J^P$ block with one fused site
   kernel per elementary pair, `sliceSum` over `Tdir`, form $C_{OO'}$, write HDF5. `--op <name>`
   selects the block (parallel runs); `--op all` does all nine (local validation).
3. **Analyse** locally from the rsynced `.h5` files (effective masses per channel; comparison with
   the published $2^\pm$ masses; level coincidence across same-$J^P$ channels).

## 4. Files

| file | role |
|---|---|
| `baryon_variants_impl_plan_claude.md` | this file |
| `baryon_variants_common_claude.py` | index helpers shared by the two generators (tuple naming, parity, $2J_3$, ops-file reader) |
| `baryon_variants_ops_gen_claude.py` | group theory only: writes `baryon_variants_ops_claude.txt` (9 highest-weight operators) and `baryon_variants_ops_all_claude.txt` (all 35 $\lvert J^P,(n_+,n_-),J,M\rangle$ states by repeated $J_-$, ket-normalised, orthonormality asserted) |
| `baryon_variants_terms_gen_claude.py` | contraction combinatorics only: writes `baryon_variants_terms_claude.txt`; `--pairs block` reads the pairs it needs from an ops file |
| `baryon_variants_terms_check_claude.py` | brute-force check of the term file on a random $16\times16$ matrix |
| `baryon_pd_common_claude.h` | shared header (records, PD rotation, sources, split-grid solver, table readers, kernel, naive reference) |
| `point2all_prop_dumper_full_claude.cc` | dump stage |
| `baryon_variants_corr_claude.cc` | contraction driver |

Not edited: `point2all_prop_dumper_claude.cc`, `two_baryon_gevp_contract_claude.cc` (owned by
another session), `baryon_theory_claude.md`.

### File formats

`baryon_variants_terms_claude.txt`:
```
# comment lines
npairs <N>
pair <A> <A'> <P(A')> <nterms>
<w> <s:c> <s:c> <s:c> <s:c> | <s:c> <s:c> <s:c> <s:c>
...
```
Each `s:c` is (PD spin, colour) of one row (before `|`) or column (after `|`) of the
$4\times4$ minor of $G$; composite index $=4s+c$. `w` is an integer (includes $P(A')$).

`baryon_variants_ops_claude.txt` and `baryon_variants_ops_all_claude.txt` (same format):
```
nops <N>
op <name> <JP> <n+> <n-> <J> <M> <ncomp>
<A> <c_A>
...
nblocks <Nb>
block <JP> <M> <nop> <name> ...
```
A block lists the operators that cross-correlate (same $J^P$ and same $M$). Coefficients are
normalised in the note's ket convention, $\sum_A c_A^2/k'_A=1$ with
$k'_A=\binom{n_+}{k_0}\binom{n_-}{k_2}$ (the number of terms of the product ket that $B_A$
stands for); the nine table entries have unit norm in this convention as written. Names in the
all-states file are `<base>_Mp2`, `_Mp1`, `_M0`, `_Mm1`, `_Mm2`.

## 5. Chunks

1. DONE 2026-09-24. Impl plan (this file) + generators + checker. Files:
   `baryon_variants_impl_plan_claude.md`, `baryon_variants_common_claude.py`,
   `baryon_variants_ops_gen_claude.py`, `baryon_variants_terms_gen_claude.py`,
   `baryon_variants_terms_check_claude.py`. Checker passes on the 157 sector pairs (8e-14) and
   the 42 block pairs; all-states table (35) orthonormal per $(n_+,n_-,M)$ group.
2. DONE 2026-09-24. Shared header + full dumper, cold $8^4$ `--rottest` (2 corners, point
   source): max norm2(diff) 1.2e-32 over all 32 blocks; 32 records, each `q<s><s'>` tag twice;
   $q_{00}=q_{33}$ on the free field. Files: `baryon_pd_common_claude.h`,
   `point2all_prop_dumper_full_claude.cc`, `Grid_sdm_build/tmp_claude.sh`
   (log `chunk2_pdfull_dumper_claude.log`).
3. DONE 2026-09-24. Contraction kernel + naive reference in the header. Files:
   `baryon_pd_common_claude.h`.
4. DONE 2026-09-24. Driver (reading, blocks, HDF5). Files: `baryon_variants_corr_claude.cc`.
5. DONE 2026-09-24. Cold validation (log `chunk45_bvar_corr_claude.log`,
   `baryon_variants_cold_check_claude.py`): `Cop_2p_40_2p_40` = $24\times$`CB_set0_corner<i>`
   of `ColdConfig.gevp.h5` to 3e-16 (point and smeared sink); `--naive-check 0`: kernel =
   host det form = naive Wick loop to 4e-14 at all 4096 sites of corner 0, both kernel modes;
   all 35 diagonals positive; 50 different-$J$ pairs vanish to 4e-17; hermiticity 3e-17.
   GOTCHA found on the way: Grid's `accelerator_for` macro declares a local `int nt`
   (`Accelerator.h:136`, `:440`) around the kernel lambda, so any kernel variable named `nt`
   is silently shadowed; symptom was wrong values at ~25% of sites and an illegal device
   memory access on the last pair. Also: `Hdf5Writer` stores the `meta` fields as HDF5
   ATTRIBUTES of the `meta` group (read with `f["meta"].attrs[...]`), not as datasets.
   Cleanup after validation (bisection diagnostics removed; backups
   `baryon_pd_common_claude_backup.h`, `baryon_variants_corr_claude_backup.cc`); the cleaned
   build re-validated identically (naive check 4e-14, python checks 0 failures).
6. 6a DONE 2026-09-24: real-gauge code test on the local SU(4) HMC checkpoint
   `/mnt/baracuda_14/grid_claude/16c/ckpoint_lat.10000` ($16^3\times8$, plaquette 0.58, M5 1.5,
   mass 0.1 study defaults): `--rottest` 1e-16 on all 16 blocks; `--naive-check 256` on the rough
   gauge field (non-diagonal colour minors) kernel = host det to 3e-14 = naive Wick to 3e-13; all 9
   ops + 35 states written (`lat16c.bvar.h5`, `lat16c.bvar_all.h5`, log `chunk6a_16c_claude.log`,
   report script `baryon_variants_report_claude.py`). No $24^3\times48$ gauge config exists
   locally (confirmed by the Grid_sdm session; only the q00 dumps are staged), so 6b (lat.1000
   physics: $2^+$ (4,0) ~3.0, $2^\pm$ vs publication, level coincidence) runs on the cluster.
7. Handoff WRITTEN 2026-09-24: `baryon_variants_production_handoff_claude.md` (this dir),
   `rsync_pull_bvar_claude.sh` + `bvar_h5_claude.py` (repo root; aggregation tested on renamed
   cold outputs). PENDING: user/remote agent runs lat.1000, then the m0.4 set, then the others;
   local analysis notebook (m_eff per channel, $2^\pm$ vs publication, level coincidence, small
   GEVP per block, backward parity partners).

## 6. Decisions log

- 2026-09-23: two-stage pipeline (dump full PD propagator once, per-operator contraction runs);
  pre-rotated PD sources; origin corner only in production, `--nsrc` knob; Python generator +
  checker with text tables read by C++; `--naive-check` on a single config; no cross-quantum-number
  contractions on real configs; only highest-weight states.
- 2026-09-24: index tuple named $A:=(\alpha,\beta,\gamma,\delta)$, propagator named $G$; the
  elementary $C_{AA'}$ of every computed pair are stored in the HDF5 output as well; the full
  35-state $(J,M)$ table is generated (group theory script separate from the determinant script)
  so that any state can be formed offline from the stored $C_{AA'}$.

## 7. Open questions

- Handoff assumes the SAME smeared source as the q00 production (w=3, N=40) so the $2^+$ (4,0)
  channel cross-checks against the existing single baryon; switch to point sources if preferred.

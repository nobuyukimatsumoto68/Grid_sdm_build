# two_baryon_gevp_contract v2 -- add SINK smearing (full symmetric smeared/point GEVP)

Spec for extending `two_baryon_gevp_contract_claude.cc` (v1) so the GEVP basis has {point, smeared}
operators at BOTH source and sink -- i.e. the full 2x2 smearing matrix incl. the cross blocks
(smeared-source/point-sink, point-source/smeared-sink) and smeared-smeared. Written for the LOCAL agent
to implement + validate on a cold config (as v1 was), since it is physics/normalization-sensitive and
must reproduce v1 exactly on the point-point block.

## Why v1 is incomplete
v1 computes `C2B[sink pair][ (source set s, source pair) ]` by peeking the point-to-all `q00[s]` at the
sink corner points -> **the sink is always a POINT**; only the SOURCE smearing varies (set by which q00
file is read). So of the 2x2 smearing matrix:

|            | point sink | smeared sink |
|------------|------------|--------------|
| point src  | v1 (SRC=point)  | MISSING |
| smeared src| v1 (SRC=smeared) | MISSING |

The smeared-SINK column is absent. A full variational GEVP needs it.

## Key point: NO new solves
`q00` is point-to-all (sink index = whole volume), so a smeared sink is obtained by covariantly
Gaussian-smearing the SINK index of the stored q00 in post-processing. All data is already on disk
(lustre5 `two_baryon_2448_*[_point]`). Only a contraction code change + the gauge config are needed.

## Design
Keep the det(8x8) kernel (`TwoBaryonCorr`) untouched -- it is agnostic to what q00 blocks you feed it.
The only new ingredient is a sink-smeared copy of q00.

1. **New input: gauge config.** Add a CLI arg, e.g. `--conf <NERSC path>` (or reuse the config index to
   locate it). Load `Umu`, extract spatial links `U[mu]=PeekIndex<LorentzIndex>(Umu,mu)` (mu=0,1,2).
2. **Sink-smear q00.** For each source set s and source corner a, build `q00_smsink[s][a]` by applying the
   SAME covariant Gaussian as the source (mirror the dumper `MakePreRotatedSmearedSource`):
   `CovariantSmearing<PeriodicGimplD>::GaussianSmear(U, field, width=3.0, niter=40, orthog=Tdir)`.
   The smear acts on the SINK colour+space index. q00 is a `LatticeColourMatrix` (row = sink colour,
   col = source colour); smear each SOURCE-colour column as a `LatticeColourVector`
   (`peekColour(q00, .,cprime)` -> smear -> `pokeColour`), or a colour-matrix-aware smear. Width/niter
   MUST equal the source's (else smeared-smeared is mis-normalised).
3. **Operator tables gain a smearing index at BOTH ends:**
   - source op = (source smearing s in {point,smeared}, source pair (a,b))
   - sink   op = (sink   smearing sm in {point,smeared}, sink   pair (c,d))
   - `C[sinkop][srcop] = TwoBaryonCorr( peek(Q(s,sm)[a],Xc), peek(Q(s,sm)[b],Xc),
                                        peek(Q(s,sm)[a],Xd), peek(Q(s,sm)[b],Xd) )`
     where `Q(s,sm=point)=q00[s]` (raw) and `Q(s,sm=smeared)=q00_smsink[s]`.
   With 2 smearings x 3 pairs at each end this is a 6x6 C matrix (= 2x2 smearing blocks x 3x3 pairs).
4. **Single input file is enough per source smearing** but you still pass BOTH source sets
   (`--src smeared.lime --src point.lime`) to populate the source-smearing dimension; the sink-smearing
   dimension is generated internally by (2).

## Output schema (suggestion)
`C2B_snk<S|P><ip>_src<S|P><jp>` for smearing in {S=smeared,P=point}, pair index ip,jp in 0..2 (36
datasets), plus single-baryon `CB_<S|P>_corner<Y>` / `CB_<S|P>_avg`, and `meta`. Keep names explicit so
the downstream python can assemble the smearing-blocked GEVP.

## Validation (local, cold config -- REQUIRED before production)
- **Regression:** the (point src, point sink) sub-block must equal v1's output bit-for-bit (same
  `C2B_snkP*_srcP*` == v1 `C2B_snk*_src*`). This anchors the refactor.
- **Symmetry:** for matched operators, `C_ij(t) ~ C_ji(t)` (point-to-all source/sink asymmetry aside);
  smeared-smeared and point-point blocks real/positive at small t.
- **Cross-block sanity:** smeared-point and point-smeared should agree with each other up to the
  source/sink exchange, and the smeared blocks should have better overlap (flatter early m_eff) than the
  point block.

## Operational -- run v2 on TUOLUMNE (recommended)
v2 now needs BOTH the q00 AND the gauge config, and the covariant smear is heavier than the v1 dets. This
flips the earlier dane-vs-tuolumne call:
- **Tuolumne**: sees lustre5 -> reads q00 AND the NERSC configs (`conf_nc4nf1_2448_*`) directly, NO
  staging; GPU makes the smear fast; runs under an allocated flux bank (no preemption). Build the v2
  binary on tuolumne (HIP, `su4_32c/build`), submit via flux mirroring the dumper scripts.
- Dane would additionally require staging the gauge configs to lustre1-3 and still suffer standby
  preemption -> not recommended for v2.

## Status of v1 outputs (keep -- they are the point-sink column)
`obs_gevp_2448_*` (smeared src / point sink) and `obs_gevp_2448_*_point` (point src / point sink) remain
valid; v2 supersedes them for the full matrix but v1 point-point is the regression anchor.

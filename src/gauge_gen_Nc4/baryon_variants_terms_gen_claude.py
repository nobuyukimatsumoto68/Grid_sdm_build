#!/usr/bin/env python3
r"""
baryon_variants_terms_gen_claude.py   (contraction combinatorics only; no group theory)

Generates the determinant term tables for the single-baryon two-point functions
C_{AA'} of the SU(4), N_f = 1 B-form components B_A.
See baryon_variants_impl_plan_claude.md Sec. 2.3-2.4 for the derivation and Sec. 4 for the
file format. The operator coefficients live in a separate script
(baryon_variants_ops_gen_claude.py); this script only needs the ops FILE for --pairs block.

Physics summary (PD basis: 0 = upper up, 1 = upper down, 2 = lower up, 3 = lower down):

  B_A(x) = \epsilon_{abcd} q^a_\alpha q^b_\beta q^c_\gamma q^d_\delta(x),
  A := (\alpha, \beta, \gamma, \delta) = sorted 4-tuple of PD indices (35 distinct), P(A) = (-1)^{n_-(A)}.

  B^\dagger_{A'} = P(A') \epsilon_{abcd} \bar q^d_{\delta'} \bar q^c_{\gamma'} \bar q^b_{\beta'} \bar q^a_{\alpha'}

  C_{AA'}(t) = \sum_x < B_A(x,t) B^\dagger_{A'}(0) >
             = P(A') \sum_{a,b \in S_4} \epsilon_a \epsilon_b det_{ij}[ G^{a_i b_j}_{A_i A'_j} ]

  where G^{a a'}_{\alpha \alpha'} is the 16 x 16 spin-colour site propagator (colour up, spin
  down), composite index 4 \alpha + a. Sorting the determinant rows/columns and merging colour
  orderings that only permute equal-spin slots gives the term table

     C_{AA'} = \sum_terms w * det[ G[rows, cols] ],   \sum |w| = 576.

Output (same directory):
  baryon_variants_terms_claude.txt   term tables for the selected pair set

Usage:
  python3 baryon_variants_terms_gen_claude.py [--pairs sector|parity|block] [--ops <ops file>]
    sector (default): all (A, A') with equal parity and J_3 equal mod 4   (157 pairs)
    parity          : all (A, A') with equal parity                       (617 pairs)
    block           : only the pairs needed by the operator blocks of the ops file
                      (default ops file baryon_variants_ops_claude.txt: 42 pairs)
"""

import itertools
import math
import os
import sys

from baryon_variants_common_claude import NSPIN, NC, perm_sign, all_tuples, aname, aparse, parity, j3_twice, read_ops_file


# --------------------------------------------------------------------------------------
# term tables
# --------------------------------------------------------------------------------------

def multiplicity_product(A):
  """prod_j k_j! with k_j the multiplicity of index j in A."""
  prod = 1
  for j in range(NSPIN):
    k = A.count(j)
    prod *= math.factorial(k)
  return prod


def side_table(A):
  r"""
  One side (rows or columns) of the determinant.

  For each colour ordering a in S_4 the slot-ordered composite indices are
  r_i = 4 A_i + a_i. Sorting the rows of the determinant multiplies it by the sign of the
  sorting permutation, so ordering a contributes \epsilon_a * sgn(sort) to the sorted key.
  Returns a dict: sorted key (4-tuple of composite indices) -> integer weight.
  """
  table = {}
  for a in itertools.permutations(range(NC)):
    eps = perm_sign(a)
    r = []
    for i in range(4):
      r.append(NSPIN * A[i] + a[i])
    order = sorted(range(4), key=r.__getitem__)
    sgn_sort = perm_sign(order)
    key = tuple(sorted(r))
    if key not in table:
      table[key] = 0
    table[key] += eps * sgn_sort
  # consistency: every surviving key carries weight +-prod k_j!, and the number of keys
  # is 24 / prod k_j!
  mult = multiplicity_product(A)
  keys = list(table.keys())
  for key in keys:
    w = table[key]
    if w == 0:
      del table[key]
      continue
    if abs(w) != mult:
      raise RuntimeError("side_table(%s): weight %d != multiplicity %d" % (aname(A), w, mult))
  if len(table) != 24 // mult:
    raise RuntimeError("side_table(%s): %d keys, expected %d" % (aname(A), len(table), 24 // mult))
  return table


def pair_terms(A, Ap):
  """
  Term list for C_{AA'}: list of (w, rows, cols) with
    w = P(A') * w_rows * w_cols,  rows/cols = sorted composite-index 4-tuples.
  """
  rows = side_table(A)
  cols = side_table(Ap)
  pm = parity(Ap)
  terms = []
  for rkey in sorted(rows.keys()):
    for ckey in sorted(cols.keys()):
      w = pm * rows[rkey] * cols[ckey]
      terms.append((w, rkey, ckey))
  total = 0
  for term in terms:
    total += abs(term[0])
  if total != 576:
    raise RuntimeError("pair_terms(%s,%s): sum|w| = %d != 576" % (aname(A), aname(Ap), total))
  return terms


def same_sector(A, Ap):
  """Selection rule: equal parity and J_3 equal mod 4 (i.e. 2 J_3 equal mod 8)."""
  if parity(A) != parity(Ap):
    return False
  if (j3_twice(A) - j3_twice(Ap)) % 8 != 0:
    return False
  return True


def same_parity(A, Ap):
  return parity(A) == parity(Ap)


def block_pairs(ops, blocks):
  """Ordered pairs (A, A') needed by all operator pairs inside each block of an ops file."""
  byname = {}
  for op in ops:
    byname[op["name"]] = op
  pairs = set()
  for block in blocks:
    for na in block[2]:
      for nb in block[2]:
        for Aa in byname[na]["coeffs"]:
          for Ab in byname[nb]["coeffs"]:
            pairs.add((aparse(Aa), aparse(Ab)))
  return sorted(pairs)


# --------------------------------------------------------------------------------------
# writer
# --------------------------------------------------------------------------------------

def sc(idx):
  """Composite index -> 's:c'."""
  return "%d:%d" % (idx // NSPIN, idx % NSPIN)


def write_terms(fname, pairs, mode):
  with open(fname, "w") as f:
    f.write("# baryon_variants_terms_claude.txt : generated by baryon_variants_terms_gen_claude.py\n")
    f.write("#\n")
    f.write("# C_{AA'} = sum_terms w * det[ G[rows, cols] ],\n")
    f.write("# G = 16x16 PD spin-colour site propagator\n")
    f.write("#\n")
    f.write("# indices = four [ spin (\\alpha) + color (a) ],\n")
    f.write("# written as \\alpha : a\n")
    f.write("#\n")
    f.write("# w includes P(A')\n")
    f.write("#\n")
    f.write("# header per pair:  pair <A> <A'> <P(A')> <#terms>\n")
    f.write("# term line:        <w> <rows of the 4x4 minor of G> | <columns of the minor>\n")
    f.write("#\n")
    f.write("# generated with --pairs %s (sector: equal parity and J_3 mod 4; block: only the pairs the ops file needs; parity: all same-parity pairs)\n" % mode)
    f.write("npairs %d\n" % len(pairs))
    nterm_total = 0
    for (A, Ap) in pairs:
      terms = pair_terms(A, Ap)
      nterm_total += len(terms)
      f.write("#\n")
      f.write("pair %s %s %d %d\n" % (aname(A), aname(Ap), parity(Ap), len(terms)))
      for (w, rows, cols) in terms:
        rs = " ".join(sc(r) for r in rows)
        cs = " ".join(sc(c) for c in cols)
        f.write("%d %s | %s\n" % (w, rs, cs))
  return nterm_total


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------

def main(argv):
  mode = "sector"
  opsfile = "baryon_variants_ops_claude.txt"
  i = 1
  while i < len(argv):
    if argv[i] == "--pairs":
      mode = argv[i + 1]
      i += 2
      continue
    if argv[i] == "--ops":
      opsfile = argv[i + 1]
      i += 2
      continue
    raise RuntimeError("unknown argument %s" % argv[i])

  As = all_tuples()
  if len(As) != 35:
    raise RuntimeError("expected 35 tuples A, got %d" % len(As))

  if mode == "sector":
    pairs = []
    for A in As:
      for Ap in As:
        if same_sector(A, Ap):
          pairs.append((A, Ap))
  elif mode == "parity":
    pairs = []
    for A in As:
      for Ap in As:
        if same_parity(A, Ap):
          pairs.append((A, Ap))
  elif mode == "block":
    ops, blocks = read_ops_file(opsfile)
    pairs = block_pairs(ops, blocks)
  else:
    raise RuntimeError("unknown --pairs mode %s" % mode)

  nterm = write_terms("baryon_variants_terms_claude.txt", pairs, mode)
  print("pairs: %d (mode %s), terms: %d" % (len(pairs), mode, nterm))

  # coverage report against the ops file(s), if present
  for fn in [opsfile, "baryon_variants_ops_all_claude.txt"]:
    if not os.path.exists(fn):
      continue
    ops, blocks = read_ops_file(fn)
    needed = block_pairs(ops, blocks)
    pairset = set(pairs)
    missing = 0
    for p in needed:
      if p not in pairset:
        missing += 1
    print("%s: needs %d pairs, %d missing from the generated set" % (fn, len(needed), missing))

  # spot checks quoted in the plan
  t0 = pair_terms(aparse("0000"), aparse("0000"))
  if len(t0) != 1 or t0[0][0] != 576:
    raise RuntimeError("0000 0000 should be a single term with w = 576")
  t1 = pair_terms(aparse("0123"), aparse("0123"))
  if len(t1) != 576:
    raise RuntimeError("0123 0123 should have 576 terms")
  print("spot checks passed: 0000x0000 -> 1 term w=576; 0123x0123 -> 576 terms")
  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv))

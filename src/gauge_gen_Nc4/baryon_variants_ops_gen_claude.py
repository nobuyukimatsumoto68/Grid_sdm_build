#!/usr/bin/env python3
r"""
baryon_variants_ops_gen_claude.py   (group theory only; no contraction combinatorics)

Writes the operator coefficient tables of the SU(4), N_f = 1 B-form baryon operators:

  baryon_variants_ops_claude.txt      the 9 highest-weight (M = J) operators, one per
                                      (J^P, n_+, n_-) sector, coefficients verbatim from
                                      baryon_group_theory/baryon_theory_claude.md Sec. 4.1
  baryon_variants_ops_all_claude.txt  all 35 states |J^P, (n_+, n_-), J, M>, obtained by
                                      repeated J_- on the highest-weight operators

Conventions (PD basis 0 = upper up, 1 = upper down, 2 = lower up, 3 = lower down):
  B_A(x) = \epsilon_{abcd} q^a_\alpha q^b_\beta q^c_\gamma q^d_\delta(x),
  A := (\alpha, \beta, \gamma, \delta) sorted (B_A is totally symmetric), 35 components.
  An operator is O = \sum_A c_A B_A.
  J_- acts slot-wise: 0 -> 1 and 2 -> 3 (upper and lower spinors lower independently).
  Ket normalisation (the note's convention "a normalised symmetric ket with k terms
  corresponds to sqrt(k) B"): B_A <-> product ket (upper quarks) x (lower quarks) with
  k'_A = binom(n_+, k_0) binom(n_-, k_2) terms, so < B_A | B_A > = 1 / k'_A and
  || O ||^2 = \sum_A c_A^2 / k'_A.  The 9 table entries have unit norm in this convention
  (asserted), and the lowered states are normalised the same way.

Blocks: operators that cross-correlate at zero momentum, i.e. same J^P and same M.

Usage: python3 baryon_variants_ops_gen_claude.py
"""

import math
import sys

from baryon_variants_common_claude import aname, aparse

SQRT3 = math.sqrt(3.0)

# Operator record: (name, J^P, n_+, n_-, J, M, {A (digit string): coefficient}).
# Highest-weight states (M = J) from baryon_theory_claude.md Sec. 4.1,
# "Highest weight construction of the irrep", coefficients verbatim.
OPERATORS = [
  ("2p_40", "2+", 4, 0, 2, 2, {"0000": 1.0}),
  ("2m_31", "2-", 3, 1, 2, 2, {"0002": 1.0}),
  ("1m_31", "1-", 3, 1, 1, 1, {"0003": 0.5 * SQRT3, "0012": -0.5 * SQRT3}),
  ("2p_22", "2+", 2, 2, 2, 2, {"0022": 1.0}),
  ("1p_22", "1+", 2, 2, 1, 1, {"0023": 1.0, "0122": -1.0}),
  ("0p_22", "0+", 2, 2, 0, 0, {"0033": 1.0 / SQRT3, "0123": -2.0 / SQRT3, "1122": 1.0 / SQRT3}),
  ("2m_13", "2-", 1, 3, 2, 2, {"0222": 1.0}),
  ("1m_13", "1-", 1, 3, 1, 1, {"0223": 0.5 * SQRT3, "1222": -0.5 * SQRT3}),
  ("2p_04", "2+", 0, 4, 2, 2, {"2222": 1.0}),
]

OP_NAME = 0
OP_JP = 1
OP_NPLUS = 2
OP_NMINUS = 3
OP_J = 4
OP_M = 5
OP_COEFFS = 6

BLOCK_ORDER = ["2+", "2-", "1-", "1+", "0+"]


# --------------------------------------------------------------------------------------
# SU(2) x SU(2) (upper x lower spin) algebra on the symmetric index space
# --------------------------------------------------------------------------------------

def lowered(coeffs):
  """
  Apply the total lowering operator J_- to sum_A c_A B_A.
  On each slot: 0 -> 1, 2 -> 3; the result is re-sorted because B_A is symmetric.
  """
  out = {}
  for As in coeffs:
    c = coeffs[As]
    A = aparse(As)
    for i in range(4):
      if A[i] == 0 or A[i] == 2:
        low = list(A)
        low[i] = A[i] + 1
        key = aname(tuple(sorted(low)))
        if key not in out:
          out[key] = 0.0
        out[key] += c
  # drop exact cancellations (e.g. J_- on the J = 0 state gives nothing)
  keys = list(out.keys())
  for key in keys:
    if abs(out[key]) < 1e-14:
      del out[key]
  return out


def kprime(A):
  """k'_A = binom(n_+, k_0) binom(n_-, k_2): number of terms of the product ket of B_A."""
  k0 = A.count(0)
  k1 = A.count(1)
  k2 = A.count(2)
  k3 = A.count(3)
  return math.comb(k0 + k1, k0) * math.comb(k2 + k3, k2)


def ket_inner(ca, cb):
  """< sum_A ca_A B_A | sum_A cb_A B_A > in the ket convention (real coefficients)."""
  val = 0.0
  for As in ca:
    if As in cb:
      val += ca[As] * cb[As] / kprime(aparse(As))
  return val


def normalised(coeffs):
  norm = math.sqrt(ket_inner(coeffs, coeffs))
  out = {}
  for As in coeffs:
    out[As] = coeffs[As] / norm
  return out


def mname_suffix(M):
  if M > 0:
    return "_Mp%d" % M
  if M < 0:
    return "_Mm%d" % (-M)
  return "_M0"


def all_states():
  """
  All 35 states |J^P, (n_+, n_-), J, M>: repeated J_- on each highest-weight operator,
  normalised in the ket convention. Names: <base>_Mp2, _Mp1, _M0, _Mm1, _Mm2.
  Asserts: the 9 table entries have unit norm; states in the same (n_+, n_-, M) group
  (different J) are orthonormal.
  """
  out = []
  for op in OPERATORS:
    J = op[OP_J]
    coeffs = dict(op[OP_COEFFS])
    M = J
    while True:
      out.append((op[OP_NAME] + mname_suffix(M), op[OP_JP], op[OP_NPLUS], op[OP_NMINUS], J, M, normalised(coeffs)))
      if M == -J:
        break
      coeffs = lowered(coeffs)
      if len(coeffs) == 0:
        raise RuntimeError("lowering %s below M = %d gave zero" % (op[OP_NAME], M))
      M -= 1
  if len(out) != 35:
    raise RuntimeError("expected 35 states, got %d" % len(out))
  for op in OPERATORS:
    n2 = ket_inner(op[OP_COEFFS], op[OP_COEFFS])
    if abs(n2 - 1.0) > 1e-12:
      raise RuntimeError("highest-weight %s has ket norm^2 %.15g, expected 1" % (op[OP_NAME], n2))
  for sa in out:
    for sb in out:
      same_group = (sa[OP_NPLUS] == sb[OP_NPLUS]) and (sa[OP_M] == sb[OP_M])
      if not same_group:
        continue
      g = ket_inner(sa[OP_COEFFS], sb[OP_COEFFS])
      expect = 0.0
      if sa[OP_NAME] == sb[OP_NAME]:
        expect = 1.0
      if abs(g - expect) > 1e-12:
        raise RuntimeError("Gram(%s, %s) = %.15g, expected %g" % (sa[OP_NAME], sb[OP_NAME], g, expect))
  return out


def operator_blocks(ops):
  """Blocks of operators with the same J^P and the same M (the ones that cross-correlate)."""
  blocks = []
  for jp in BLOCK_ORDER:
    Ms = []
    for op in ops:
      if op[OP_JP] == jp and op[OP_M] not in Ms:
        Ms.append(op[OP_M])
    Ms.sort(reverse=True)
    for M in Ms:
      names = []
      for op in ops:
        if op[OP_JP] == jp and op[OP_M] == M:
          names.append(op[OP_NAME])
      blocks.append((jp, M, names))
  return blocks


# --------------------------------------------------------------------------------------
# writer
# --------------------------------------------------------------------------------------

def write_ops(fname, ops, title):
  with open(fname, "w") as f:
    f.write("# %s : generated by baryon_variants_ops_gen_claude.py\n" % fname)
    f.write("# %s\n" % title)
    f.write("#\n")
    f.write("# O = sum_A c_A B_A, normalised in the ket convention sum_A c_A^2 / k'_A = 1,\n")
    f.write("# k'_A = binom(n+, k_0) binom(n-, k_2)\n")
    f.write("#\n")
    f.write("# header per operator:  op <name> <J^P> <n+> <n-> <J> <M> <#components>\n")
    f.write("# component line:       <A> <c_A>\n")
    f.write("# block line:           block <J^P> <M> <#ops> <names ...>   (operators that cross-correlate)\n")
    f.write("nops %d\n" % len(ops))
    for op in ops:
      coeffs = op[OP_COEFFS]
      f.write("#\n")
      f.write("op %s %s %d %d %d %d %d\n" % (op[OP_NAME], op[OP_JP], op[OP_NPLUS], op[OP_NMINUS], op[OP_J], op[OP_M], len(coeffs)))
      for As in sorted(coeffs.keys()):
        f.write("%s %.17g\n" % (As, coeffs[As]))
    blocks = operator_blocks(ops)
    f.write("#\n")
    f.write("nblocks %d\n" % len(blocks))
    for block in blocks:
      f.write("block %s %d %d %s\n" % (block[0], block[1], len(block[2]), " ".join(block[2])))


def main(argv):
  if len(argv) > 1:
    raise RuntimeError("no arguments expected")
  write_ops("baryon_variants_ops_claude.txt", OPERATORS, "the 9 highest-weight (M = J) operators")
  states = all_states()
  write_ops("baryon_variants_ops_all_claude.txt", states, "all 35 states |J^P, (n+,n-), J, M> by repeated J_- on the highest-weight operators")
  print("wrote baryon_variants_ops_claude.txt (9 operators, %d blocks)" % len(operator_blocks(OPERATORS)))
  print("wrote baryon_variants_ops_all_claude.txt (35 states, %d blocks); orthonormality asserted" % len(operator_blocks(states)))
  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv))

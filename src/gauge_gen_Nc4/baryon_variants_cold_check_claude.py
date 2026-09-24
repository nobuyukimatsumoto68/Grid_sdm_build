#!/usr/bin/env python3
r"""
baryon_variants_cold_check_claude.py

Cold-configuration (free field, 8^4, point source) checks of the baryon_variants_corr
output (Chunk 5 of baryon_variants_impl_plan_claude.md):

  (a) ColdConfig.bvar.h5 (9 operators, --pairs block, --sink-smear 2.0 20):
      Cop_2p_40_2p_40 == 24 * mean_corners CB_set0_corner<i> of ColdConfig.gevp.h5
      (the gevp contraction stores 24 det q00 = C_{0000,0000} / 24), point and smeared sink.
  (b) ColdConfig.bvar_all.h5 (35 states, --pairs all, --cross):
      operator pairs with different J (same parity and M mod 4) vanish exactly on the free
      field; diagonal correlators are positive at small t; hermiticity C_ij = conj(C_ji).

Usage: python3 baryon_variants_cold_check_claude.py [bvar.h5] [bvar_all.h5] [gevp.h5]
Exit status 0 on success.
"""

import sys

import h5py
import numpy as np

from baryon_variants_common_claude import read_ops_file


def read_corr(f, name):
  d = f[name + "/data"][:]
  return np.array([complex(x[0], x[1]) for x in d])


def check_reference(bvar, gevp, ncorner):
  nfail = 0
  fb = h5py.File(bvar, "r")
  fg = h5py.File(gevp, "r")
  for (ours, theirs) in [("Cop_2p_40_2p_40", "CB_set0_corner"), ("Copss_2p_40_2p_40", "CBss_set0_corner")]:
    if ours not in fb:
      print("  %s not in %s (skipped)" % (ours, bvar))
      continue
    c = read_corr(fb, ours)
    ref = np.zeros_like(c)
    for i in range(ncorner):
      ref += read_corr(fg, theirs + str(i))
    ref = 24.0 * ref / ncorner
    rel = np.max(np.abs(c - ref)) / np.max(np.abs(ref))
    print("  %s vs 24*%s<i>: max rel diff %.3e" % (ours, theirs, rel))
    print("    ours: ", np.real(c[:4]))
    print("    ref : ", np.real(ref[:4]))
    if rel > 1e-5:
      print("  FAIL: reference mismatch")
      nfail += 1
  return nfail


def op_quantum_numbers(opsfile):
  ops, blocks = read_ops_file(opsfile)
  qn = {}
  for op in ops:
    qn[op["name"]] = (op["J"], op["JP"][-1], op["M"])
  return qn


def split_pair_name(name, opnames):
  """'Cop_<O>_<O2>' -> (O, O2) using the known operator names (they contain underscores)."""
  body = name[len("Cop_"):]
  for a in opnames:
    if body.startswith(a + "_"):
      b = body[len(a) + 1:]
      if b in opnames:
        return a, b
  return None, None


def check_all_states(bvar_all, opsfile):
  nfail = 0
  f = h5py.File(bvar_all, "r")
  qn = op_quantum_numbers(opsfile)
  opnames = sorted(qn.keys(), key=len, reverse=True)
  corr = {}
  for k in f.keys():
    if not k.startswith("Cop_"):
      continue
    a, b = split_pair_name(k, opnames)
    if a is None:
      raise RuntimeError("cannot split dataset name %s" % k)
    corr[(a, b)] = read_corr(f, k)
  print("  %d operator-pair correlators" % len(corr))
  scale = 0.0
  for key in corr:
    scale = max(scale, np.max(np.abs(corr[key])))
  ncross = 0
  maxcross = 0.0
  ndiag = 0
  for (a, b) in corr:
    Ja, Pa, Ma = qn[a]
    Jb, Pb, Mb = qn[b]
    c = corr[(a, b)]
    if a == b:
      ndiag += 1
      if not (c[1].real > 0.0 and c[2].real > 0.0):
        print("  FAIL: diagonal %s not positive at t=1,2: %s %s" % (a, c[1], c[2]))
        nfail += 1
      continue
    if Ja != Jb:
      ncross += 1
      v = np.max(np.abs(c)) / scale
      maxcross = max(maxcross, v)
      if v > 1e-10:
        print("  FAIL: different-J pair %s x %s not zero: rel %.3e" % (a, b, v))
        nfail += 1
  print("  diagonals: %d (positive at t=1,2)" % ndiag)
  print("  different-J pairs: %d, max |C|/scale = %.3e (must vanish on the free field)" % (ncross, maxcross))
  # hermiticity report
  maxherm = 0.0
  for (a, b) in corr:
    if (b, a) in corr:
      maxherm = max(maxherm, np.max(np.abs(corr[(a, b)] - np.conj(corr[(b, a)]))) / scale)
  print("  hermiticity max |C_ab - conj(C_ba)|/scale = %.3e (report)" % maxherm)
  # E/T2 mixing report for J=2, M=+2 vs M=-2
  for (a, b) in sorted(corr):
    Ja, Pa, Ma = qn[a]
    Jb, Pb, Mb = qn[b]
    if Ja == 2 and Jb == 2 and Ma == 2 and Mb == -2:
      v = np.max(np.abs(corr[(a, b)])) / scale
      print("  M=+2 x M=-2 (E/T2 mixing) %s x %s: max |C|/scale = %.3e (report)" % (a, b, v))
  return nfail


def main(argv):
  bvar = "ColdConfig.bvar.h5"
  bvar_all = "ColdConfig.bvar_all.h5"
  gevp = "ColdConfig.gevp.h5"
  if len(argv) > 1:
    bvar = argv[1]
  if len(argv) > 2:
    bvar_all = argv[2]
  if len(argv) > 3:
    gevp = argv[3]
  nfail = 0
  print("(a) reference check against %s" % gevp)
  fb = h5py.File(bvar, "r")
  ncorner = int(fb["meta"].attrs["nCorner"][0])
  fb.close()
  nfail += check_reference(bvar, gevp, ncorner)
  print("(b) all-states checks on %s" % bvar_all)
  nfail += check_all_states(bvar_all, "baryon_variants_ops_all_claude.txt")
  print("failures: %d" % nfail)
  if nfail > 0:
    return 1
  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv))

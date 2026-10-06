#!/usr/bin/env python3
r"""
aaron_ops_convert_claude.py

Converts A. S. Meyer's (LLNL) single-flavour 4q operators (cubic-group irreps at rest,
little-group irreps in flight, total momentum up to (2,1,1); private communication 2026-09-28,
directory output_quad00_op/i2c2 of his tarball) into the text tables read by
baryon_variants_mom_corr_claude.cc:

  baryon_mom_ops_claude.txt            momenta + operator groups with complex coefficients
  baryon_variants_terms_all_claude.txt determinant term tables of all 35^2 = 1225 pairs
                                       (in flight parity is not a symmetry, so every pair is
                                       needed); tables from pair_terms() of
                                       baryon_variants_terms_gen_claude.py (imported, unchanged)

Aaron's projectors map onto our PD basis (0 = upper up, 1 = upper down, 2 = lower up,
3 = lower down):
  PARITYPLUS SPINUP -> 0, PARITYPLUS SPINDN -> 1, PARITYMINUS SPINUP -> 2, PARITYMINUS SPINDN -> 3.
Each term is a product of four UP quark fields at one site (identifier y0) with momentum MOM;
B_A is totally symmetric in A, so a term is filed under sorted(A) and coefficients are summed.

Operator O_i(p) = \sum_A c^i_A B_A(p). A group = (frame, irrep, component): all operators
(unique IDs) of one irrep row at one momentum direction p; the contraction forms
C_{ij}(p,t) = \sum_{A,A'} c^i_A conj(c^j_{A'}) C_{AA'}(p,t) inside each group.

ops table format:
  nmom <n>
  mom <k> <px> <py> <pz>                      (integer momenta, units 2\pi/L)
  ngroup <n>
  group <frame> <irrep> <comp> <k> <nops>     (k = momentum index)
  op <uid> <ncomp>
  <A> <re> <im>

Self-checks: one momentum per group; every operator ID appears in every component of its
(frame, irrep); every operator has fixed (n_+, n_-) (rotations preserve upper/lower spinors).

Usage: python3 aaron_ops_convert_claude.py [--aaron <i2c2 dir>]
"""

import glob
import os
import re
import sys

from baryon_variants_common_claude import all_tuples, aname, parity
from baryon_variants_terms_gen_claude import pair_terms, sc

PROJ = {
  ("PARITYPLUS", "SPINUP"): 0,
  ("PARITYPLUS", "SPINDN"): 1,
  ("PARITYMINUS", "SPINUP"): 2,
  ("PARITYMINUS", "SPINDN"): 3,
}

FRAME_ORDER = ["000", "100", "110", "111", "200", "210", "211"]

TINY = 1e-14

DEFAULT_AARON = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "../../../baryon_group_theory/aaron_4q_ops/output_quad00_op/i2c2")


def read_op(fname):
  """Returns (momentum string, {A-string: complex}) of one Aaron .op file."""
  coeffs = {}
  mom = None
  fac = None
  idx = []
  with open(fname) as f:
    lines = f.read().splitlines()
  lines.append("")
  for line in lines:
    w = line.split()
    if len(w) == 0:
      if fac is not None:
        if len(idx) != 4:
          raise RuntimeError("%s: term with %d quark lines" % (fname, len(idx)))
        key = aname(tuple(sorted(idx)))
        if key not in coeffs:
          coeffs[key] = 0.0
        coeffs[key] += fac
      fac = None
      idx = []
      continue
    if w[0] == "FACTOR":
      fac = complex(float(w[1]), float(w[2]))
      continue
    if w[0] == "MOM":
      if mom is not None and mom != w[1]:
        raise RuntimeError("%s: more than one momentum (%s, %s)" % (fname, mom, w[1]))
      mom = w[1]
      if w[2] != "y0":
        raise RuntimeError("%s: identifier %s (expected y0)" % (fname, w[2]))
      continue
    if w[0] == "UP":
      if w[3] != "y0":
        raise RuntimeError("%s: identifier %s (expected y0)" % (fname, w[3]))
      idx.append(PROJ[(w[1], w[2])])
      continue
    raise RuntimeError("%s: unexpected line '%s'" % (fname, line))
  clean = {}
  for key in coeffs:
    c = coeffs[key]
    re_ = c.real
    im_ = c.imag
    if abs(re_) < TINY:
      re_ = 0.0
    if abs(im_) < TINY:
      im_ = 0.0
    if re_ != 0.0 or im_ != 0.0:
      clean[key] = complex(re_, im_)
  return mom, clean


def n_minus(A):
  n = 0
  for ch in A:
    if ch in "23":
      n += 1
  return n


def main(argv):
  aaron = DEFAULT_AARON
  i = 1
  while i < len(argv):
    if argv[i] == "--aaron":
      aaron = argv[i + 1]
      i += 2
      continue
    raise RuntimeError("unknown argument %s" % argv[i])

  # groups[(frame, irrep, comp)] = {uid: (mom, coeffs)}
  groups = {}
  pat = re.compile(r"^(\w+?)_(\d{3})\.(\d+)\.(\d+)\.op$")
  nfiles = 0
  for frame in FRAME_ORDER:
    files = sorted(glob.glob(os.path.join(aaron, "pcom" + frame, "*", "*.op")))
    if len(files) == 0:
      raise RuntimeError("no .op files for frame %s under %s" % (frame, aaron))
    for fname in files:
      m = pat.match(os.path.basename(fname))
      if m is None:
        raise RuntimeError("cannot parse file name %s" % fname)
      irrep = m.group(1)
      uid = m.group(3)
      comp = m.group(4)
      mom, coeffs = read_op(fname)
      key = (frame, irrep, comp)
      if key not in groups:
        groups[key] = {}
      groups[key][uid] = (mom, coeffs)
      nfiles += 1

  # momenta, in order of first appearance (frame order, then file order)
  moms = []
  for key in sorted(groups, key=lambda k: (FRAME_ORDER.index(k[0]), k[1], int(k[2]))):
    for uid in sorted(groups[key]):
      mom = groups[key][uid][0]
      if mom not in moms:
        moms.append(mom)

  # self-checks
  uids_of = {}
  for key in groups:
    fi = (key[0], key[1])
    if fi not in uids_of:
      uids_of[fi] = set()
    uids_of[fi] |= set(groups[key].keys())
  for key in groups:
    gmoms = set(groups[key][uid][0] for uid in groups[key])
    if len(gmoms) != 1:
      raise RuntimeError("group %s has %d momenta" % (str(key), len(gmoms)))
    if set(groups[key].keys()) != uids_of[(key[0], key[1])]:
      raise RuntimeError("group %s lacks some operator IDs of its irrep" % str(key))
    for uid in groups[key]:
      coeffs = groups[key][uid][1]
      nms = set(n_minus(A) for A in coeffs)
      if len(nms) != 1:
        raise RuntimeError("op %s of group %s mixes n_- = %s" % (uid, str(key), str(sorted(nms))))

  keys = sorted(groups, key=lambda k: (FRAME_ORDER.index(k[0]), k[1], int(k[2])))
  with open("baryon_mom_ops_claude.txt", "w") as f:
    f.write("# baryon_mom_ops_claude.txt : generated by aaron_ops_convert_claude.py\n")
    f.write("# A. S. Meyer (LLNL) single-flavour (i2c2) local 4q operators, from %s\n" % os.path.normpath(aaron))
    f.write("# O_i(p) = sum_A c^i_A B_A(p), A = sorted PD tuple (0 = upper up, 1 = upper down, 2 = lower up, 3 = lower down)\n")
    f.write("# group = (frame, irrep, component): one irrep row at one momentum direction; ops = unique IDs\n")
    f.write("#\n")
    f.write("# mom <k> <px> <py> <pz>\n")
    f.write("# group <frame> <irrep> <comp> <k> <nops>\n")
    f.write("# op <uid> <ncomp>\n")
    f.write("# <A> <re> <im>\n")
    f.write("nmom %d\n" % len(moms))
    for k in range(len(moms)):
      p = moms[k].strip("[]").split(",")
      f.write("mom %d %s %s %s\n" % (k, p[0], p[1], p[2]))
    f.write("ngroup %d\n" % len(keys))
    for key in keys:
      uids = sorted(groups[key])
      k = moms.index(groups[key][uids[0]][0])
      f.write("#\n")
      f.write("group %s %s %s %d %d\n" % (key[0], key[1], key[2], k, len(uids)))
      for uid in uids:
        coeffs = groups[key][uid][1]
        f.write("op %s %d\n" % (uid, len(coeffs)))
        for A in sorted(coeffs):
          c = coeffs[A]
          f.write("%s %.17g %.17g\n" % (A, c.real, c.imag))

  # all 1225 pairs, same file format as baryon_variants_terms_gen_claude.py
  As = all_tuples()
  nterm_total = 0
  with open("baryon_variants_terms_all_claude.txt", "w") as f:
    f.write("# baryon_variants_terms_all_claude.txt : generated by aaron_ops_convert_claude.py\n")
    f.write("#   (pair_terms() of baryon_variants_terms_gen_claude.py, all 35^2 pairs)\n")
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
    f.write("npairs %d\n" % (len(As) * len(As)))
    for A in As:
      for Ap in As:
        terms = pair_terms(A, Ap)
        nterm_total += len(terms)
        f.write("#\n")
        f.write("pair %s %s %d %d\n" % (aname(A), aname(Ap), parity(Ap), len(terms)))
        for (w, rows, cols) in terms:
          rs = " ".join(sc(r) for r in rows)
          cs = " ".join(sc(c) for c in cols)
          f.write("%d %s | %s\n" % (w, rs, cs))

  nfi = {}
  for key in keys:
    fi = key[0]
    if fi not in nfi:
      nfi[fi] = 0
    nfi[fi] += 1
  print("read %d .op files -> %d groups, %d momenta" % (nfiles, len(keys), len(moms)))
  for fr in FRAME_ORDER:
    print("  frame %s: %d groups" % (fr, nfi.get(fr, 0)))
  print("wrote baryon_mom_ops_claude.txt and baryon_variants_terms_all_claude.txt (%d pairs, %d terms)" % (len(As) * len(As), nterm_total))


if __name__ == "__main__":
  main(sys.argv)

#!/usr/bin/env python3
r"""
baryon_variants_report_claude.py

Single-configuration report for baryon_variants_corr output files:
  (a) <bvar.h5>      9 operators: effective masses m_eff(t) = log C(t)/C(t+1) of every
                     diagonal C_{OO}, and the same for the smeared sink if present;
                     off-diagonal/diagonal ratios per J^P block.
  (b) <bvar_all.h5>  35 states: hermiticity |C_ab - conj(C_ba)| / scale, positivity of the
                     diagonals at small t, per-block summary.

No ensemble average here (one configuration): numbers are for orientation and code checks.
Usage: python3 baryon_variants_report_claude.py [bvar.h5] [bvar_all.h5]
"""

import sys

import h5py
import numpy as np


def read_corr(f, name):
  d = f[name + "/data"][:]
  return np.array([complex(x[0], x[1]) for x in d])


def meta_str(f, key):
  v = f["meta"].attrs[key][0]
  if isinstance(v, bytes):
    return v.decode()
  return str(v)


def effmass(c):
  T = len(c)
  out = []
  for t in range(T - 1):
    a = c[t].real
    b = c[t + 1].real
    if a > 0.0 and b > 0.0:
      out.append(np.log(a / b))
    else:
      out.append(float("nan"))
  return out


def report_ops(fname):
  f = h5py.File(fname, "r")
  names = meta_str(f, "opNames").split(",")
  jps = meta_str(f, "opJP").split(",")
  T = int(f["meta"].attrs["T"][0])
  print("(a) %s : config %s, T = %d, corners %d" % (fname, meta_str(f, "config"), T, int(f["meta"].attrs["nCorner"][0])))
  for prefix in ["Cop_", "Copss_"]:
    have = False
    for i in range(len(names)):
      key = prefix + names[i] + "_" + names[i]
      if key not in f:
        continue
      have = True
      c = read_corr(f, key)
      me = effmass(c)
      print("  %-6s %-6s C(1..3) = %s" % (prefix[:-1], names[i], " ".join("%.4e" % c[t].real for t in range(1, min(4, T)))))
      print("         m_eff(t) = %s" % " ".join("%.3f" % m for m in me))
    if not have:
      continue
    # off-diagonal / sqrt(diag diag) at t = 1, 2 per block
    seen = set()
    for i in range(len(names)):
      for j in range(len(names)):
        if i == j or jps[i] != jps[j]:
          continue
        key = prefix + names[i] + "_" + names[j]
        if key not in f or (i, j) in seen:
          continue
        seen.add((i, j))
        cij = read_corr(f, key)
        cii = read_corr(f, prefix + names[i] + "_" + names[i])
        cjj = read_corr(f, prefix + names[j] + "_" + names[j])
        r = []
        for t in range(1, min(4, T)):
          den = np.sqrt(abs(cii[t].real * cjj[t].real))
          r.append(cij[t].real / den if den > 0 else float("nan"))
        print("  %-6s %s x %s (%s): normalised overlap t=1..3 = %s" % (prefix[:-1], names[i], names[j], jps[i], " ".join("%.3f" % x for x in r)))


def report_all(fname):
  f = h5py.File(fname, "r")
  names = meta_str(f, "opNames").split(",")
  print("(b) %s : %d states" % (fname, len(names)))
  corr = {}
  for i in range(len(names)):
    for j in range(len(names)):
      key = "Cop_" + names[i] + "_" + names[j]
      if key in f:
        corr[(names[i], names[j])] = read_corr(f, key)
  scale = 0.0
  for key in corr:
    scale = max(scale, np.max(np.abs(corr[key])))
  maxherm = 0.0
  for (a, b) in corr:
    if (b, a) in corr:
      maxherm = max(maxherm, np.max(np.abs(corr[(a, b)] - np.conj(corr[(b, a)]))) / scale)
  nneg = 0
  for a in names:
    c = corr[(a, a)]
    if not (c[1].real > 0.0 and c[2].real > 0.0):
      nneg += 1
      print("  diagonal %s not positive at t=1,2: %s %s" % (a, c[1], c[2]))
  print("  %d pair correlators; hermiticity max |C_ab - conj(C_ba)|/scale = %.3e; non-positive diagonals: %d of %d" % (len(corr), maxherm, nneg, len(names)))


def main(argv):
  bvar = "lat16c.bvar.h5"
  bvar_all = "lat16c.bvar_all.h5"
  if len(argv) > 1:
    bvar = argv[1]
  if len(argv) > 2:
    bvar_all = argv[2]
  report_ops(bvar)
  report_all(bvar_all)
  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv))

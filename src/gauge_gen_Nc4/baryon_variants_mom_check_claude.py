#!/usr/bin/env python3
r"""
baryon_variants_mom_check_claude.py

Checks of baryon_variants_mom_corr_claude output (see baryon_variants_mom_impl_plan_claude.md,
Chunk 3). Subcommands:

  legacy NEW REF   every Cel_/Cop_/Celss_/Copss_ dataset of the production file REF is present
                   in NEW and equal (relative to max_t |REF|)
  rest NEW         rest-frame relation between Aaron's irrep operators and our highest-weight
                   operators (exact on the cold field, cubic symmetry):
                     J = 0, 1:  Cop_{XY} = C^{\Gamma}_{XY} / sqrt(n_X n_Y)
                     J = 2:     Cop_{XY} = (1/2) [ C^{E}_{XY} / sqrt(n^E_X n^E_Y) + C^{T2}_{XY} / sqrt(n^{T2}_X n^{T2}_Y) ]
                   n = ket norm^2 \sum_A |c_A|^2 / k'_A of Aaron's operator (positive overlaps, from
                   baryon_group_theory/compare_aaron_ops_claude.py)
  comps NEW        every per-component Cmomraw_ equals its (frame, irrep) average Cmom_ (exact on
                   the cold field with a cubic-symmetric source); needs --mom-raw
  herm NEW         C_{ij} = conj(C_{ji}) for every (frame, irrep) class
  nop0 NOP0 FULL   NOP0 (--no-p0) has no Cel_/Cop_ datasets and the same Cmom_ datasets as FULL
  eff NEW          effective energies log(C(t)/C(t+1)) of the first diagonal operator per class

Usage: python3 baryon_variants_mom_check_claude.py <subcommand> <files...> [--ops baryon_mom_ops_claude.txt]
"""

import math
import sys

import h5py
import numpy as np

TOL = 1e-10

# our highest-weight operator -> Aaron rest-frame (irrep, uid) list (compare_aaron_ops_claude.py)
REST_MAP = {
  "0p_22": [("a1p", "00040")],
  "1p_22": [("t1p", "00035")],
  "1m_31": [("t1m", "00029")],
  "1m_13": [("t1m", "00033")],
  "2p_40": [("eep", "00036"), ("t2p", "00039")],
  "2p_22": [("eep", "00034"), ("t2p", "00038")],
  "2p_04": [("eep", "00031"), ("t2p", "00042")],
  "2m_31": [("eem", "00030"), ("t2m", "00041")],
  "2m_13": [("eem", "00032"), ("t2m", "00037")],
}


def read_corr(f, name):
  d = f[name + "/data"][:]
  return d["re"] + 1j * d["im"]


def kprime(A):
  k0 = A.count("0")
  k1 = A.count("1")
  k2 = A.count("2")
  k3 = A.count("3")
  return math.comb(k0 + k1, k0) * math.comb(k2 + k3, k2)


def read_ketnorms(opsfile):
  """{(frame, irrep, uid): ket norm^2 of the first component found}"""
  out = {}
  key = None
  cur = None
  with open(opsfile) as f:
    for line in f:
      w = line.split()
      if len(w) == 0 or w[0].startswith("#"):
        continue
      if w[0] in ("nmom", "mom", "ngroup"):
        continue
      if w[0] == "group":
        key = (w[1], w[2])
        continue
      if w[0] == "op":
        cur = (key[0], key[1], w[1])
        if cur in out:
          cur = None
        else:
          out[cur] = 0.0
        continue
      if cur is not None:
        c = complex(float(w[1]), float(w[2]))
        out[cur] += abs(c) ** 2 / kprime(w[0])
  return out


def rel(a, ref):
  scale = max(1e-300, np.max(np.abs(ref)))
  return np.max(np.abs(a - ref)) / scale


def verdict(name, worst, tol=TOL):
  status = "PASSED" if worst < tol else "FAILED"
  print("%s: max rel = %.3e  %s" % (name, worst, status))
  return worst < tol


def check_legacy(new, ref):
  fn = h5py.File(new, "r")
  fr = h5py.File(ref, "r")
  worst = 0.0
  n = 0
  for k in fr.keys():
    if not (k.startswith("Cel") or k.startswith("Cop")):
      continue
    if k not in fn:
      print("  MISSING in new file: %s" % k)
      worst = 1.0
      continue
    r = rel(read_corr(fn, k), read_corr(fr, k))
    worst = max(worst, r)
    n += 1
  print("legacy: compared %d datasets" % n)
  return verdict("legacy p=0 datasets vs production", worst)


def check_rest(new, opsfile):
  f = h5py.File(new, "r")
  n = read_ketnorms(opsfile)
  ok = True
  for sink in ["", "ss"]:
    if "Cop%s_0p_22_0p_22" % sink not in f:
      continue
    worst = 0.0
    for X in REST_MAP:
      for Y in REST_MAP:
        name = "Cop%s_%s_%s" % (sink, X, Y)
        if name not in f:
          continue
        if len(REST_MAP[X]) != len(REST_MAP[Y]):
          continue
        pred = 0.0
        for (irX, uX), (irY, uY) in zip(REST_MAP[X], REST_MAP[Y]):
          if irX != irY:
            raise RuntimeError("irrep mismatch %s %s" % (irX, irY))
          c = read_corr(f, "Cmom%s_000_%s_%s_%s" % (sink, irX, uX, uY))
          pred = pred + c / math.sqrt(n[("000", irX, uX)] * n[("000", irY, uY)])
        pred = pred / len(REST_MAP[X])
        r = rel(pred, read_corr(f, name))
        worst = max(worst, r)
        if r > TOL:
          print("  %s: rel %.3e" % (name, r))
    ok = verdict("rest-frame Aaron ops vs our Cop%s" % sink, worst) and ok
  return ok


def classes_of(f):
  """{(frame_irrep): [uids]} from mom_meta.classList"""
  cl = f["mom_meta"].attrs["classList"][0]
  if isinstance(cl, bytes):
    cl = cl.decode()
  out = {}
  for item in cl.split("|"):
    name, ncomp, uids = item.split(":")
    out[name] = (int(ncomp), uids.split(","))
  return out


def check_comps(new):
  f = h5py.File(new, "r")
  cls = classes_of(f)
  keys = list(f.keys())
  worst = 0.0
  nraw = 0
  for sink in ["", "ss"]:
    for k in keys:
      pre = "Cmomraw%s_" % sink
      if not k.startswith(pre):
        continue
      if sink == "" and k.startswith("Cmomrawss_"):
        continue
      rest = k[len(pre):]
      w = rest.split("_")
      # <frame>_<irrep>_c<comp>_<uid_i>_<uid_j>
      avg = "Cmom%s_%s_%s_%s_%s" % (sink, w[0], w[1], w[3], w[4])
      ca = read_corr(f, avg)
      diag = "Cmom%s_%s_%s_%s_%s" % (sink, w[0], w[1], w[3], w[3])
      scale = max(1e-300, np.max(np.abs(read_corr(f, diag))))
      r = np.max(np.abs(read_corr(f, k) - ca)) / scale
      if r > TOL and r > worst:
        print("  %s: rel %.3e (to the diagonal scale)" % (k, r))
      worst = max(worst, r)
      nraw += 1
  print("comps: compared %d per-component datasets with their class average" % nraw)
  return verdict("components equal within each (frame, irrep)", worst)


def check_herm(new):
  f = h5py.File(new, "r")
  cls = classes_of(f)
  worst = 0.0
  for sink in ["", "ss"]:
    for name in cls:
      if "Cmom%s_%s_%s_%s" % (sink, name, cls[name][1][0], cls[name][1][0]) not in f:
        continue
      uids = cls[name][1]
      for i in uids:
        for j in uids:
          cij = read_corr(f, "Cmom%s_%s_%s_%s" % (sink, name, i, j))
          cji = read_corr(f, "Cmom%s_%s_%s_%s" % (sink, name, j, i))
          sii = np.max(np.abs(read_corr(f, "Cmom%s_%s_%s_%s" % (sink, name, i, i))))
          sjj = np.max(np.abs(read_corr(f, "Cmom%s_%s_%s_%s" % (sink, name, j, j))))
          scale = max(1e-300, math.sqrt(sii * sjj))
          r = np.max(np.abs(cij - np.conj(cji))) / scale
          worst = max(worst, r)
  return verdict("hermiticity C_ij = conj(C_ji)", worst, 1e-8)


def check_nop0(nop0, full):
  fa = h5py.File(nop0, "r")
  fb = h5py.File(full, "r")
  legacy = [k for k in fa.keys() if k.startswith("Cel") or k.startswith("Cop")]
  print("nop0: %d legacy datasets in the --no-p0 file (expect 0)" % len(legacy))
  worst = 0.0
  n = 0
  for k in fb.keys():
    if not k.startswith("Cmom_"):
      continue
    if k not in fa:
      print("  MISSING in --no-p0 file: %s" % k)
      worst = 1.0
      continue
    worst = max(worst, rel(read_corr(fa, k), read_corr(fb, k)))
    n += 1
  print("nop0: compared %d Cmom_ datasets" % n)
  return verdict("--no-p0 Cmom_ vs full run", worst) and len(legacy) == 0


def show_eff(new):
  f = h5py.File(new, "r")
  cls = classes_of(f)
  print("%-14s %5s %22s %22s %22s" % ("class", "nops", "E_eff(0)", "E_eff(1)", "E_eff(2)"))
  for name in cls:
    u = cls[name][1][0]
    c = read_corr(f, "Cmom_%s_%s_%s" % (name, u, u))
    es = []
    for t in range(3):
      if c[t + 1] != 0 and (c[t] / c[t + 1]).real > 0:
        es.append("%.6f" % math.log((c[t] / c[t + 1]).real))
      else:
        es.append("nan")
    print("%-14s %5d %22s %22s %22s" % (name, len(cls[name][1]), es[0], es[1], es[2]))
  return True


def main(argv):
  opsfile = "baryon_mom_ops_claude.txt"
  args = []
  i = 1
  while i < len(argv):
    if argv[i] == "--ops":
      opsfile = argv[i + 1]
      i += 2
      continue
    args.append(argv[i])
    i += 1
  cmd = args[0]
  if cmd == "legacy":
    ok = check_legacy(args[1], args[2])
  elif cmd == "rest":
    ok = check_rest(args[1], opsfile)
  elif cmd == "comps":
    ok = check_comps(args[1])
  elif cmd == "herm":
    ok = check_herm(args[1])
  elif cmd == "nop0":
    ok = check_nop0(args[1], args[2])
  elif cmd == "eff":
    ok = show_eff(args[1])
  else:
    raise RuntimeError("unknown subcommand %s" % cmd)
  return 0 if ok else 1


if __name__ == "__main__":
  sys.exit(main(sys.argv))

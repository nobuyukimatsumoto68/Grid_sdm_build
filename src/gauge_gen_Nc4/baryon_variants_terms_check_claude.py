#!/usr/bin/env python3
r"""
baryon_variants_terms_check_claude.py

Brute-force check of baryon_variants_terms_claude.txt on a random complex 16 x 16 matrix G
(composite index 4*spin + colour), independent of the generator's determinant reorganisation.

For every pair (A, A') in the file it compares

  table:  sum_terms w * det[ G[rows, cols] ]

against the naive Wick form (baryon_variants_impl_plan_claude.md Sec. 2.3)

  naive:  P(A') \sum_{\pi \in S_4} \sum_{a,b \in S_4} \epsilon_a \epsilon_b
                 \prod_i G[ 4 A_i + a_i , 4 A'_{\pi(i)} + b_i ]

which is the sum over the 24 quark pairings with the colour epsilons contracted explicitly.
Also checks sum |w| = 576 per pair and that P(A') in the file matches (-1)^{n_-(A')}.

Usage: python3 baryon_variants_terms_check_claude.py [terms_file] [--seed N]
Exit status 0 on success, 1 on any mismatch.
"""

import itertools
import sys

import numpy as np

NSPIN = 4
NC = 4


def perm_sign(p):
  n = len(p)
  inv = 0
  for i in range(n):
    for j in range(i + 1, n):
      if p[i] > p[j]:
        inv += 1
  if inv % 2 == 0:
    return 1
  return -1


def parity(A):
  nminus = 0
  for i in A:
    if i >= 2:
      nminus += 1
  if nminus % 2 == 0:
    return 1
  return -1


def parse_sc(tok):
  s, c = tok.split(":")
  return NSPIN * int(s) + int(c)


def read_terms(fname):
  """Returns list of (A, Ap, pfile, terms) with terms = list of (w, rows, cols)."""
  pairs = []
  with open(fname) as f:
    lines = []
    for line in f:
      line = line.strip()
      if len(line) == 0:
        continue
      if line[0] == "#":
        continue
      lines.append(line)
  pos = 0
  head = lines[pos].split()
  if head[0] != "npairs":
    raise RuntimeError("expected 'npairs', got %s" % lines[pos])
  npairs = int(head[1])
  pos += 1
  for ip in range(npairs):
    tok = lines[pos].split()
    if tok[0] != "pair":
      raise RuntimeError("expected 'pair' at line %d" % pos)
    A = tuple(int(ch) for ch in tok[1])
    Ap = tuple(int(ch) for ch in tok[2])
    pfile = int(tok[3])
    nterm = int(tok[4])
    pos += 1
    terms = []
    for it in range(nterm):
      tok = lines[pos].split()
      w = int(tok[0])
      bar = tok.index("|")
      rows = tuple(parse_sc(t) for t in tok[1:bar])
      cols = tuple(parse_sc(t) for t in tok[bar + 1:])
      if len(rows) != 4 or len(cols) != 4:
        raise RuntimeError("bad term line %d" % pos)
      terms.append((w, rows, cols))
      pos += 1
    pairs.append((A, Ap, pfile, terms))
  return pairs


def table_value(G, terms):
  val = 0.0 + 0.0j
  for (w, rows, cols) in terms:
    sub = G[np.ix_(rows, cols)]
    val += w * np.linalg.det(sub)
  return val


def naive_value(G, A, Ap):
  perms = list(itertools.permutations(range(NC)))
  signs = {}
  for p in perms:
    signs[p] = perm_sign(p)
  val = 0.0 + 0.0j
  for pi in perms:
    for a in perms:
      ea = signs[a]
      for b in perms:
        eb = signs[b]
        prod = 1.0 + 0.0j
        for i in range(4):
          prod *= G[NSPIN * A[i] + a[i], NSPIN * Ap[pi[i]] + b[i]]
        val += ea * eb * prod
  return parity(Ap) * val


def main(argv):
  fname = "baryon_variants_terms_claude.txt"
  seed = 12345
  i = 1
  while i < len(argv):
    if argv[i] == "--seed":
      seed = int(argv[i + 1])
      i += 2
      continue
    fname = argv[i]
    i += 1

  rng = np.random.default_rng(seed)
  G = rng.standard_normal((16, 16)) + 1j * rng.standard_normal((16, 16))

  pairs = read_terms(fname)
  nfail = 0
  maxrel = 0.0
  for (A, Ap, pfile, terms) in pairs:
    tag = "%s %s" % ("".join(str(x) for x in A), "".join(str(x) for x in Ap))
    if pfile != parity(Ap):
      print("FAIL %s: file parity %d != %d" % (tag, pfile, parity(Ap)))
      nfail += 1
    wsum = 0
    for term in terms:
      wsum += abs(term[0])
    if wsum != 576:
      print("FAIL %s: sum|w| = %d" % (tag, wsum))
      nfail += 1
    vt = table_value(G, terms)
    vn = naive_value(G, A, Ap)
    scale = max(1.0, abs(vn))
    rel = abs(vt - vn) / scale
    if rel > maxrel:
      maxrel = rel
    if rel > 1e-10:
      print("FAIL %s: table %s naive %s rel %.3e" % (tag, vt, vn, rel))
      nfail += 1
  print("checked %d pairs from %s, max relative deviation %.3e, failures %d" % (len(pairs), fname, maxrel, nfail))
  if nfail > 0:
    return 1
  return 0


if __name__ == "__main__":
  sys.exit(main(sys.argv))

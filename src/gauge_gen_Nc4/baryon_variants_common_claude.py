#!/usr/bin/env python3
r"""
baryon_variants_common_claude.py

Index helpers shared by the group-theory script (baryon_variants_ops_gen_claude.py) and the
determinant-table script (baryon_variants_terms_gen_claude.py).

PD Dirac index: 0 = upper up, 1 = upper down, 2 = lower up, 3 = lower down.
A := (\alpha, \beta, \gamma, \delta) = sorted 4-tuple of PD indices, written as digits (e.g. 0002).
"""

import itertools

NSPIN = 4
NC = 4


def perm_sign(p):
  """Sign of the permutation p (a tuple/list of 0..n-1), by inversion count."""
  n = len(p)
  inv = 0
  for i in range(n):
    for j in range(i + 1, n):
      if p[i] > p[j]:
        inv += 1
  if inv % 2 == 0:
    return 1
  return -1


def all_tuples():
  """The 35 sorted 4-tuples A over {0,1,2,3}."""
  out = []
  for A in itertools.combinations_with_replacement(range(NSPIN), 4):
    out.append(tuple(A))
  return out


def aname(A):
  """Tuple A as digit string, e.g. (0,0,0,2) -> '0002'."""
  return "".join(str(i) for i in A)


def aparse(s):
  """Inverse of aname."""
  return tuple(int(ch) for ch in s)


def parity(A):
  """P(A) = (-1)^{n_-}, n_- = number of lower (2 or 3) indices."""
  nminus = 0
  for i in A:
    if i >= 2:
      nminus += 1
  if nminus % 2 == 0:
    return 1
  return -1


def j3_twice(A):
  """2 J_3 = (n_0 + n_2) - (n_1 + n_3)."""
  up = 0
  down = 0
  for i in A:
    if i % 2 == 0:
      up += 1
    else:
      down += 1
  return up - down


def read_ops_file(fname):
  """
  Reader for the ops files written by baryon_variants_ops_gen_claude.py.
  Returns (ops, blocks) with
    ops    = list of dicts {name, JP, nplus, nminus, J, M, coeffs (dict A-string -> float)}
    blocks = list of (JP, M, [names])
  """
  lines = []
  with open(fname) as f:
    for line in f:
      line = line.strip()
      if len(line) == 0:
        continue
      if line[0] == "#":
        continue
      lines.append(line)
  pos = 0
  tok = lines[pos].split()
  if tok[0] != "nops":
    raise RuntimeError("%s: expected 'nops'" % fname)
  nops = int(tok[1])
  pos += 1
  ops = []
  for io in range(nops):
    tok = lines[pos].split()
    if tok[0] != "op":
      raise RuntimeError("%s: expected 'op' at line %d" % (fname, pos))
    op = {}
    op["name"] = tok[1]
    op["JP"] = tok[2]
    op["nplus"] = int(tok[3])
    op["nminus"] = int(tok[4])
    op["J"] = int(tok[5])
    op["M"] = int(tok[6])
    ncomp = int(tok[7])
    pos += 1
    coeffs = {}
    for ic in range(ncomp):
      tok = lines[pos].split()
      coeffs[tok[0]] = float(tok[1])
      pos += 1
    op["coeffs"] = coeffs
    ops.append(op)
  tok = lines[pos].split()
  if tok[0] != "nblocks":
    raise RuntimeError("%s: expected 'nblocks'" % fname)
  nblocks = int(tok[1])
  pos += 1
  blocks = []
  for ib in range(nblocks):
    tok = lines[pos].split()
    if tok[0] != "block":
      raise RuntimeError("%s: expected 'block' at line %d" % (fname, pos))
    nop = int(tok[3])
    names = tok[4:4 + nop]
    blocks.append((tok[1], int(tok[2]), names))
    pos += 1
  return ops, blocks

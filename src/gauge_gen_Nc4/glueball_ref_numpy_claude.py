#!/usr/bin/env python3
r"""
glueball_ref_numpy_claude.py

Independent NumPy reference of HiRep's glueball 1-pt functions (measure_1pt_glueballs, APE depth n, levels 0..3) on
one Grid NERSC configuration, compared with the HiRep measure-only h5. Plan: glueball_grid_impl_plan_claude.md,
Chunk 8.1 (checks the operator table and our reading of HiRep before any GPU code).

Methods (as HiRep LibHR/Observables/glueballs.c, LibHR/Utils/spatial_transformations.c, LibHR/Utils/suN_utils.c):
  APE smearing, M. Albanese et al., Phys. Lett. B192 (1987) 163: V = (1-\alpha) U + (\alpha/4) \sum spatial staples;
  projection: W = V det(V)^{-1/N} (principal root), then the unitary polar factor of W (here via SVD, N. J. Higham,
  Functions of Matrices, SIAM 2008, ch. 8; HiRep uses Jacobi, Numerical Recipes sec. 11.1 -- same unique factor);
  blocking (M. Teper, Phys. Lett. B183 (1987) 345, HiRep variant): U'_k(x) = U_k(x) U_k(x + s k), s -> 2 s;
  operators: glueball_ops_table_gen_claude.py (parsed from HiRep's glueballs_op.c).

Usage:
  python3 glueball_ref_numpy_claude.py <nersc file> <glueballs_op.c> [<hirep h5> <traj>] [--depth 4] [--alpha 0.33]
          [--nlev 4] [--save out.npy]
"""
import argparse
import re
import sys

import numpy as np

import glueball_ops_table_gen_claude as gen


def read_nersc(fname):
    """Grid NERSC: header up to END_HEADER, data [t][z][y][x][mu][Nc][Nc] complex IEEE64BIG.
    Returns U[d] for d = 1,2,3 (x,y,z), each (X,Y,Z,T,N,N), plus the header dict."""
    raw = open(fname, "rb").read()
    end = raw.index(b"END_HEADER\n") + len(b"END_HEADER\n")
    hdr = {}
    for line in raw[:end].decode().split("\n"):
        if "=" in line:
            k, v = line.split("=", 1)
            hdr[k.strip()] = v.strip()
    nd = [int(hdr["DIMENSION_%d" % i]) for i in (1, 2, 3, 4)]
    m = re.match(r"4D_SU(\d+)_GAUGE", hdr["DATATYPE"])
    nc = int(m.group(1))
    a = np.frombuffer(raw[end:], dtype=">f8").astype(np.float64)
    a = a.reshape(nd[3], nd[2], nd[1], nd[0], 4, nc, nc, 2)
    c = a[..., 0] + 1j * a[..., 1]
    c = np.transpose(c, (3, 2, 1, 0, 4, 5, 6))
    U = {}
    for d in (1, 2, 3):
        U[d] = np.ascontiguousarray(c[:, :, :, :, d - 1])
    return U, hdr, nd, nc


def shift(A, d, n):
    """A(x + n e_d)"""
    return np.roll(A, -n, axis=d - 1)


def dag(A):
    return np.conj(np.swapaxes(A, -1, -2))


def project(V, nc):
    det = np.linalg.det(V)
    W = V * np.power(det, -1.0 / nc)[..., None, None]
    A, sv, Bh = np.linalg.svd(W)
    return A @ Bh


def ape_step(U, s, alpha, nc):
    out = {}
    for k in (1, 2, 3):
        st = np.zeros_like(U[k])
        for j in (1, 2, 3):
            if j == k:
                continue
            up = U[j] @ shift(U[k], j, s) @ dag(shift(U[j], k, s))
            Ujm = shift(U[j], j, -s)
            dn = dag(Ujm) @ shift(U[k], j, -s) @ shift(shift(U[j], j, -s), k, s)
            st = st + up + dn
        V = (1.0 - alpha) * U[k] + (alpha / 4.0) * st
        out[k] = project(V, nc)
    return out


def block(U, s):
    return {k: U[k] @ shift(U[k], k, s) for k in (1, 2, 3)}


def path_trace(U, steps, s):
    """left-to-right product along the signed steps, starting at every site"""
    off = [0, 0, 0]
    W = None
    for st in steps:
        d = abs(st)
        if st > 0:
            L = U[d]
            for e in (1, 2, 3):
                if off[e - 1] != 0:
                    L = shift(L, e, off[e - 1] * s)
            off[d - 1] += 1
        else:
            off[d - 1] -= 1
            L = U[d]
            for e in (1, 2, 3):
                if off[e - 1] != 0:
                    L = shift(L, e, off[e - 1] * s)
            L = dag(L)
        if W is None:
            W = L
        else:
            W = W @ L
    return np.trace(W, axis1=-2, axis2=-1)


def eval_level(U, tab, s, L, nd):
    X, Y, Z, T = nd
    x = np.arange(X)[:, None, None]
    y = np.arange(Y)[None, :, None]
    z = np.arange(Z)[None, None, :]
    phase = []
    for p in tab["moms"]:
        phase.append(np.exp(2j * np.pi * (x * p[0] + y * p[1] + z * p[2]) / L))
    npth = len(tab["paths"])
    Cp = np.zeros((len(tab["moms"]), npth, T), dtype=complex)
    Cm = np.zeros((len(tab["moms"]), npth, T), dtype=complex)
    for ip, steps in enumerate(tab["paths"]):
        tr = path_trace(U, steps, s)
        for im, ph in enumerate(phase):
            Cp[im, ip] = np.einsum("xyz,xyzt->t", ph, tr.real)
            Cm[im, ip] = 1j * np.einsum("xyz,xyzt->t", ph, tr.imag)
    ops = np.zeros((len(tab["ops"]), T), dtype=complex)
    for io, o in enumerate(tab["ops"]):
        acc = np.zeros(T, dtype=complex)
        for pm, ip, m, v in o["terms"]:
            c = v * np.exp(1j * np.pi * m * s / L)
            if pm == "p":
                acc = acc + c * Cp[o["mom"], ip]
            else:
                acc = acc + c * Cm[o["mom"], ip]
        ops[io] = acc
    return ops


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("nersc")
    ap.add_argument("opsrc")
    ap.add_argument("h5", nargs="?")
    ap.add_argument("traj", nargs="?", type=int)
    ap.add_argument("--depth", type=int, default=4)
    ap.add_argument("--alpha", type=float, default=0.33)
    ap.add_argument("--nlev", type=int, default=4)
    ap.add_argument("--save", default=None)
    a = ap.parse_args()

    tab = gen.parse_glueballs_op(a.opsrc)
    U, hdr, nd, nc = read_nersc(a.nersc)
    L = nd[0]
    print("config", a.nersc, "dims", nd, "Nc", nc, "header plaquette", hdr.get("PLAQUETTE"))
    res = np.zeros((a.nlev, len(tab["ops"]), nd[3]), dtype=complex)
    s = 1
    for lev in range(a.nlev):
        res[lev] = eval_level(U, tab, s, L, nd)
        print("level", lev, "step", s, "done", flush=True)
        if lev + 1 < a.nlev:
            for it in range(a.depth):
                U = ape_step(U, s, a.alpha, nc)
            U = block(U, s)
            s = 2 * s
    if a.save is not None:
        np.save(a.save, res)
        print("saved", a.save)
    if a.h5 is not None:
        import h5py
        f = h5py.File(a.h5, "r")
        traj = list(f["traj"][:])
        n = traj.index(a.traj)
        g = f["glue_ape%d" % a.depth][n]
        ref = g[..., 0] + 1j * g[..., 1]
        print("HiRep h5", a.h5, "traj", a.traj, "row", n, "shape", ref.shape)
        for lev in range(a.nlev):
            num = np.abs(res[lev] - ref[lev]).max(axis=1)
            den = np.abs(ref[lev]).max(axis=1)
            rel = num / np.maximum(den, 1e-300)
            worst = int(np.argmax(rel))
            print("level %d: max_op rel diff %.3e (op %d %s, |ref| %.3e); median rel %.3e; min |ref| over ops %.3e"
                  % (lev, rel.max(), worst, tab["ops"][worst]["name"], den[worst], np.median(rel), den.min()))


if __name__ == "__main__":
    main()

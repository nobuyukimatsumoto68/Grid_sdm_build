#!/usr/bin/env python3
r"""
glueball_h5_aggregate_claude.py

Collect the per-config files <dir>/glue.<conf>.h5 written by glueball_h5_claude.cc into ONE h5 in the HiRep layout of
suN_glueballs_h5_claude (hirep_su4/PureGauge), so hirep_su4/QuenchedSU4_claude/analyze_glueball_h5_claude.py reads it
(use --gevp-only or skip torelons: the Grid code has no /tor). Optionally compare with a HiRep h5 on the common configs.
Plan: glueball_grid_impl_plan_claude.md (Chunks 8.2/8.3).

Output datasets: traj (N,), plaq (N,), glue_ape<d> (N, nlev, nops, T, 2), topo (0, 11, 4), topo_traj (0,);
attrs GLB_T, GLB_X, GLB_Y, GLB_Z, NG, nblkstart = 0, nblkend = nlev - 1, APEsmear, APE_depth_0, n_APE_depths = 1,
total_n_glue_op, source = "grid glueball_h5_claude".

Usage:
  python3 glueball_h5_aggregate_claude.py <dir> <out.h5> [--L 16] [--compare <hirep.h5>] [--tol 1e-10]
  (out.h5 must not exist; the script never overwrites)
"""
import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np


def get(f, key):
    if key in f:
        return f[key][()]
    return f.attrs[key]


def get_scalar(f, key):
    # Grid's Hdf5Writer stores scalars as 1-element attribute arrays
    v = np.asarray(get(f, key)).reshape(-1)
    if v.size != 1:
        raise ValueError("%s is not a scalar" % key)
    return v[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("out")
    ap.add_argument("--L", type=int, required=True, help="spatial extent (Grid files do not store it)")
    ap.add_argument("--nc", type=int, default=4)
    ap.add_argument("--compare", default=None)
    ap.add_argument("--tol", type=float, default=1e-10)
    a = ap.parse_args()

    if os.path.exists(a.out):
        print("output exists, not overwriting:", a.out)
        sys.exit(1)
    files = glob.glob(os.path.join(a.dir, "glue.*.h5"))
    confs = []
    for fn in files:
        m = re.search(r"glue\.(\d+)\.h5$", fn)
        if m:
            confs.append(int(m.group(1)))
    confs = sorted(confs)
    if len(confs) == 0:
        print("no glue.<conf>.h5 in", a.dir)
        sys.exit(1)
    rows = []
    plaq = []
    depth = None
    alpha = None
    shape = None
    for c in confs:
        with h5py.File(os.path.join(a.dir, "glue.%d.h5" % c), "r") as f:
            sh = tuple(int(x) for x in np.atleast_1d(get(f, "glue_shape")))
            g = np.asarray(get(f, "glue"), dtype=np.float64).reshape(sh)
            d = int(get_scalar(f, "ape_depth"))
            al = float(get_scalar(f, "ape_alpha"))
            if int(get_scalar(f, "conf")) != c:
                print("conf mismatch in file", c)
                sys.exit(1)
            if shape is None:
                shape = sh
                depth = d
                alpha = al
            if sh != shape or d != depth or al != alpha:
                print("inconsistent shape/depth/alpha at conf", c)
                sys.exit(1)
            rows.append(g)
            plaq.append(float(get_scalar(f, "plaq")))
    glue = np.stack(rows)
    nlev, nops, T, two = shape
    dname = "glue_ape%d" % depth
    with h5py.File(a.out, "w") as f:
        f.create_dataset("traj", data=np.array(confs, dtype=np.int32))
        f.create_dataset("plaq", data=np.array(plaq))
        f.create_dataset(dname, data=glue)
        f.create_dataset("topo", data=np.zeros((0, 11, 4)))
        f.create_dataset("topo_traj", data=np.zeros((0,), dtype=np.int32))
        f.attrs["GLB_T"] = np.int32(T)
        f.attrs["GLB_X"] = np.int32(a.L)
        f.attrs["GLB_Y"] = np.int32(a.L)
        f.attrs["GLB_Z"] = np.int32(a.L)
        f.attrs["NG"] = np.int32(a.nc)
        f.attrs["nblkstart"] = np.int32(0)
        f.attrs["nblkend"] = np.int32(nlev - 1)
        f.attrs["APEsmear"] = np.float64(alpha)
        f.attrs["APE_depth_0"] = np.int32(depth)
        f.attrs["n_APE_depths"] = np.int32(1)
        f.attrs["total_n_glue_op"] = np.int32(nops)
        f.attrs["source"] = "grid glueball_h5_claude"
    print("wrote", a.out, ":", len(confs), "configs", confs[0], "..", confs[-1], dname, glue.shape)

    if a.compare is None:
        return
    ok = True
    with h5py.File(a.compare, "r") as f:
        htraj = list(f["traj"][:])
        hg = f[dname][:]
        hplaq = f["plaq"][:]
    common = [c for c in confs if c in htraj]
    print("compare with", a.compare, ": common configs", common)
    if len(common) == 0:
        print("COMPARE FAILED: no common configs")
        sys.exit(1)
    for c in common:
        i = confs.index(c)
        j = htraj.index(c)
        G = glue[i][..., 0] + 1j * glue[i][..., 1]
        H = hg[j][..., 0] + 1j * hg[j][..., 1]
        dp = abs(plaq[i] - hplaq[j])
        line = "conf %d: |plaq diff| %.2e;" % (c, dp)
        for lev in range(nlev):
            num = np.abs(G[lev] - H[lev]).max(axis=1)
            den = np.abs(H[lev]).max(axis=1)
            rel = num / np.maximum(den, 1e-300)
            line += " lev%d max %.2e (op %d) med %.2e;" % (lev, rel.max(), int(np.argmax(rel)), np.median(rel))
            if not rel.max() <= a.tol:
                ok = False
        print(line)
    if ok:
        print("COMPARE PASSED (tol %.1e)" % a.tol)
    else:
        print("COMPARE FAILED (tol %.1e)" % a.tol)
        sys.exit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# meson_combine_claude.py
#
# Chunk 2: combine the momentum-projected CONNECTED meson correlator (produced by
# baryons_0000_dirac_claude into per-config HDF5, datasets "<channel>_t_p<idx>")
# with the DISCONNECTED loop-loop correlator (Rebbi's corr_d.p{0..4} text files)
# to form the physical flavor-singlet correlator
#
#     C_phys(p_class, t) = 2 * D(p_class, t) - C_conn(p_class, t)
#
# (the "2*corr - ccorr" of makecorr_96_d.f90:353). See
# meson_momproj_impl_plan_claude.md for the full derivation and open questions.
#
# IMPORTANT (OQ1, unresolved): the RELATIVE normalization between D and C is not
# yet pinned (the factor-3 at makecorr_96_d.f90:231 etc.). It is exposed here as
# `rel_norm` and defaults to 1.0 -- calibrate before trusting absolute 2D-C.

import re
import numpy as np
import h5py

# ----------------------------------------------------------------------
# Momentum table (matches average_trace2.f90 and mom_table in the .cc).
# p^2 classes -> the stored momentum indices belonging to each class.
# Rebbi's pcomp(m)/2 equals the number of stored indices per class, so the
# connected class average is a plain mean over these indices (real part).
# ----------------------------------------------------------------------
PCLASS_INDICES = {
    0: [0],
    1: [1, 2, 3],
    2: [4, 5, 6, 7, 8, 9],
    3: [10, 11, 12, 13],
    4: [14, 15, 16],
}

N_MOM = 17

# Connected channel name (HDF5 dataset prefix) -> disconnected trace_<G> dir.
CHANNEL_TO_TRACEDIR = {
    "I_I": "trace_id",
    "G5_G5": "trace_g5",
    "GTG5_GTG5": "trace_gtg5",
    "GXG5_GXG5": "trace_gxg5",
    "GYG5_GYG5": "trace_gyg5",
    "GZG5_GZG5": "trace_gzg5",
    "GT_GT": "trace_gt",
    "GX_GX": "trace_gx",
    "GY_GY": "trace_gy",
    "GZ_GZ": "trace_gz",
}

# Grid HDF5 complex storage: compound (re, im).
_DT = np.dtype([('re', np.double), ('im', np.double)])


def _read_complex(f, name):
    tmp = np.array(f[name]['data'][:], dtype=_DT)
    return tmp['re'] + 1j * tmp['im']


def read_connected_h5(filename, channel):
    """Read the 17 momentum datasets for one channel from one config's h5.

    Returns a complex array of shape (N_MOM, NT): C(p_idx, t), t = 0 .. NT-1.
    """
    out = None
    with h5py.File(filename, 'r') as f:
        for ip in range(N_MOM):
            c = _read_complex(f, f"{channel}_t_p{ip}")
            if out is None:
                out = np.zeros((N_MOM, c.shape[0]), dtype=complex)
            out[ip] = c
    return out


def fold_connected(c_t):
    """Fold C(t) about t=0 and take the real part.

    c_t: complex array (..., NT). Returns real array (..., NT/2+1):
    Cf[0]=C[0]; Cf[NT/2]=C[NT/2]; Cf[t]=(C[t]+C[NT-t])/2 for t=1..NT/2-1.
    """
    nt = c_t.shape[-1]
    nh = nt // 2
    cf = np.zeros(c_t.shape[:-1] + (nh + 1,), dtype=complex)
    cf[..., 0] = c_t[..., 0]
    cf[..., nh] = c_t[..., nh]
    for t in range(1, nh):
        cf[..., t] = 0.5 * (c_t[..., t] + c_t[..., nt - t])
    return np.real(cf)


def connected_pclass(filename, channel):
    """Per-config connected correlator averaged into p^2 classes.

    Returns real array (5, NT/2+1): C_conn(p_class=0..4, t).
    """
    c_pidx = read_connected_h5(filename, channel)   # (N_MOM, NT)
    c_fold = fold_connected(c_pidx)                 # (N_MOM, NT/2+1)
    out = np.zeros((5, c_fold.shape[-1]), dtype=float)
    for m, idxs in PCLASS_INDICES.items():
        out[m] = np.mean(c_fold[idxs, :], axis=0)
    return out


def read_corr_d(filename):
    """Read a Rebbi corr_d.p{m} file: per-config blocks of 't  D(t)'.

    Returns real array (nconf, NT/2+1): D(p_class fixed by file, t).
    """
    blocks = []
    cur = []
    with open(filename) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                if cur:
                    blocks.append(cur)
                    cur = []
                continue
            parts = s.split()
            # 't  D(t)'
            cur.append(float(parts[1]))
    if cur:
        blocks.append(cur)
    return np.array(blocks, dtype=float)


def combine_CminusD(D_pclass, C_pclass, rel_norm=1.0, nf=1.0):
    """Form the Nf=1 physical singlet C - Nf*rel_norm*D for one p-class, one config.

    rel_norm scales D relative to C (= 1.0 here). The overall sign is irrelevant
    for masses; C-D is used. Both inputs are real arrays of shape (NT/2+1,).
    """
    return C_pclass - nf * rel_norm * D_pclass


def read_corr_d_indexed(filename):
    """Like read_corr_d but also return the per-config indices.

    The index is parsed from each '#File: ... data_files.<idx>' header.
    Returns (indices (nconf,) int, data (nconf, NT/2+1) real).
    """
    idxs = []
    blocks = []
    cur = []
    cur_idx = None
    with open(filename) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                m = re.search(r'data_files\.(\d+)', s)
                if m:
                    if cur:
                        blocks.append(cur)
                        idxs.append(cur_idx)
                        cur = []
                    cur_idx = int(m.group(1))
                continue
            cur.append(float(s.split()[1]))
    if cur:
        blocks.append(cur)
        idxs.append(cur_idx)
    return np.array(idxs, dtype=int), np.array(blocks, dtype=float)


def collect_connected(ens_dir, channel):
    """Read all mesons_conn.<idx>.h5 in ens_dir for one channel.

    Returns (indices (nconf,) int sorted, C (nconf, 5, NT/2+1) real),
    where C[:, m, :] is the folded class-averaged connected correlator for
    p^2 class m.
    """
    import glob
    import os
    files = glob.glob(os.path.join(ens_dir, "mesons_conn.*.h5"))
    pairs = []
    for f in files:
        m = re.search(r'mesons_conn\.(\d+)\.h5$', f)
        if m:
            pairs.append((int(m.group(1)), f))
    pairs.sort()
    idxs = []
    arrs = []
    n_bad = 0
    for idx, f in pairs:
        try:
            a = connected_pclass(f, channel)         # (5, NT/2+1)
        except (OSError, KeyError) as e:
            n_bad += 1
            print(f"WARNING: skipping unreadable/incomplete {os.path.basename(f)}: {e}")
            continue
        idxs.append(idx)
        arrs.append(a)
    if n_bad:
        print(f"WARNING: skipped {n_bad} bad connected file(s)")
    return np.array(idxs, dtype=int), np.array(arrs, dtype=float)


def jackknife(samples, binsize=1):
    """Leave-one-out jackknife over axis 0 (configs).

    samples: (nconf, ...). Bins of `binsize` consecutive configs are averaged
    first (remainder dropped). Returns (mean, err) each with shape samples[0].
    """
    x = np.asarray(samples, dtype=float)
    n = x.shape[0]
    nb = n // binsize
    x = x[:nb * binsize].reshape((nb, binsize) + x.shape[1:]).mean(axis=1)
    total = x.sum(axis=0)
    jk = (total - x) / (nb - 1)            # leave-one-out means, (nb, ...)
    mean = jk.mean(axis=0)
    err = np.sqrt((nb - 1) / nb * np.sum((jk - mean) ** 2, axis=0))
    return mean, err


# ----------------------------------------------------------------------
# Smoke test of the CONNECTED side only (runnable now, no production data).
# Verifies: (a) all 17 momenta read; (b) p0 class fold == fold of the original
# "_t" dataset; (c) class averaging produces finite numbers.
# Usage: python3 meson_combine_claude.py <connected.h5>
# ----------------------------------------------------------------------
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python3 meson_combine_claude.py <connected.h5>")
        sys.exit(1)
    fname = sys.argv[1]
    ok = True
    for ch in CHANNEL_TO_TRACEDIR:
        c_pidx = read_connected_h5(fname, ch)
        nt = c_pidx.shape[1]
        cls = connected_pclass(fname, ch)
        # p0 class is just index 0; compare its fold to fold of "_t".
        with h5py.File(fname, 'r') as f:
            c_t = _read_complex(f, ch + "_t")
        p0_from_t = fold_connected(c_t)
        maxdiff = np.abs(cls[0] - p0_from_t).max()
        finite = np.all(np.isfinite(cls))
        status = "OK" if (maxdiff < 1e-12 and finite) else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"[{ch:11s}] NT={nt} classes shape={cls.shape} "
              f"p0-vs-_t maxdiff={maxdiff:.2e} finite={finite} {status}")
    print("CONNECTED-SIDE SMOKE TEST:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

#!/usr/bin/env python3
# Cross-check: the momentum-projected connected meson at p=0 (dataset
# "<channel>_t_p0") must equal the original zero-momentum dataset
# "<channel>_t" bit-for-bit, since MakePhase({0,0,0}) = 1.
#
# Usage: python3 verify_momproj_p0_claude.py file.h5

import sys
import h5py
import numpy as np

dt = np.dtype([('re', np.double), ('im', np.double)])

channels = ["I_I", "G5_G5", "GTG5_GTG5", "GXG5_GXG5", "GYG5_GYG5",
            "GZG5_GZG5", "GT_GT", "GX_GX", "GY_GY", "GZ_GZ"]

def read_complex(f, name):
    tmp = np.array(f[name]['data'][:], dtype=dt)
    return tmp['re'] + 1j * tmp['im']

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 verify_momproj_p0_claude.py file.h5")
        sys.exit(1)

    fname = sys.argv[1]
    ok = True
    with h5py.File(fname, 'r') as f:
        for ch in channels:
            n_t = ch + "_t"
            n_p0 = ch + "_t_p0"
            if n_t not in f or n_p0 not in f:
                print(f"[{ch}] missing dataset ({n_t} or {n_p0})")
                ok = False
                continue
            a = read_complex(f, n_t)
            b = read_complex(f, n_p0)
            if a.shape != b.shape:
                print(f"[{ch}] shape mismatch {a.shape} vs {b.shape}")
                ok = False
                continue
            maxdiff = np.abs(a - b).max()
            if np.allclose(a, b, rtol=1e-12, atol=1e-14):
                print(f"[{ch}] p0 matches _t   (max diff = {maxdiff:.3e})")
            else:
                print(f"[{ch}] p0 DIFFERS      (max diff = {maxdiff:.3e})")
                ok = False

    print("ALL MATCH" if ok else "MISMATCH FOUND")
    sys.exit(0 if ok else 1)

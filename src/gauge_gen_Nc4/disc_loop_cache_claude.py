#!/usr/bin/env python3
# disc_loop_cache_claude.py
#
# Re-do the average_trace2 momentum projection of the disconnected loops
# tr[Gamma S(x,x)] KEEPING THE FULL COMPLEX VALUE (average_trace2_claude.f90
# keeps only the real part: "tr(x,y,z,0) = a  ! because the traces should be
# real"). gamma5-hermiticity S^\dagger = \gamma_5 S \gamma_5 gives
#   tr[Gamma S]^* = tr[S \gamma_5 Gamma^\dagger \gamma_5]
# so the loop is REAL for Gamma in {1, \gamma_5, \gamma_\mu\gamma_5} and
# IMAGINARY for Gamma = \gamma_\mu.
#
# Reads the raw SciDAC files obsdir/traces.<g>.<conf> (big-endian complex128,
# x fastest, data block at byte offset DATA_OFFSET), and writes
#   loops_<ens>_claude.npz : confs (Nconf,), gammas (10,),
#                            R, I (Nconf, 10, NT, 17) complex128
# R = projection of the real-part field r(x) = Re tr(x), I = projection of the
# imaginary-part field l(x) = Im tr(x), so tr = r + i l and
#   R[..., t, 0]   = (1/NS^3) sum_x r(x)  (= avtr of average_trace2)
#   R[..., t, j>0] = (1/NS^3) sum_x e^{2 \pi i p_j.x/NS} r(x), same 16-momentum
#                    ptable as average_trace2 (index 3 means p=-1); same for I.
# r and l are real fields, so r(-p) = conj(r(p)) and l(-p) = conj(l(p)).
# (The first version stored L = projection of the full complex tr; for p0,
#  Re L = R and Im L = I.)
#
# Usage:
#   OMP_NUM_THREADS=4 python3 disc_loop_cache_claude.py <ens_dir> [out_dir]

import sys
import os
import glob
import numpy as np

NS = 24
NT = 48
DATA_OFFSET = 1992
NBYTES = NS * NS * NS * NT * 16
GAMMAS = ["id", "g5", "gx", "gy", "gz", "gt", "gxg5", "gyg5", "gzg5", "gtg5"]

# (px, py, pz) for momentum index 1..16, copied from average_trace2_claude.f90
PTABLE = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
          (1, 1, 0), (1, 0, 1), (0, 1, 1),
          (1, 3, 0), (1, 0, 3), (0, 1, 3),
          (1, 1, 1), (3, 1, 1), (1, 3, 1), (1, 1, 3),
          (2, 0, 0), (0, 2, 0), (0, 0, 2)]


def make_expk():
    # expk[k, j] = exp(2 \pi i j k / NS) for j=0,1,2 and j=3 -> conj(j=1)
    k = np.arange(NS)
    expk = np.zeros((NS, 4), dtype=complex)
    expk[:, 0] = 1.0
    expk[:, 1] = np.exp(2j * np.pi * k / NS)
    expk[:, 2] = np.exp(2j * np.pi * 2 * k / NS)
    expk[:, 3] = np.conj(expk[:, 1])
    return expk


def read_trace(path):
    with open(path, "rb") as f:
        f.seek(DATA_OFFSET)
        buf = f.read(NBYTES)
    if len(buf) != NBYTES:
        return None
    # memory order t, z, y, x (x fastest)
    return np.frombuffer(buf, dtype=">c16").astype(complex).reshape(NT, NS, NS, NS)


def project(d, expk):
    out = np.zeros((NT, 1 + len(PTABLE)), dtype=complex)
    out[:, 0] = d.mean(axis=(1, 2, 3))
    for j, (px, py, pz) in enumerate(PTABLE):
        out[:, j + 1] = np.einsum("tzyx,x,y,z->t", d, expk[:, px], expk[:, py], expk[:, pz], optimize=True) / NS**3
    return out


def main():
    ens_dir = sys.argv[1].rstrip("/")
    out_dir = "disc_loops_claude"
    if len(sys.argv) > 2:
        out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    ens = os.path.basename(ens_dir)

    # configs that have ALL gammas
    conf_sets = []
    for g in GAMMAS:
        fs = glob.glob(os.path.join(ens_dir, f"traces.{g}.*"))
        conf_sets.append(set(int(f.rsplit(".", 1)[1]) for f in fs))
    confs = sorted(set.intersection(*conf_sets))
    print(f"{ens}: {len(confs)} configs with all {len(GAMMAS)} gammas")

    expk = make_expk()
    # L = np.zeros((len(confs), len(GAMMAS), NT, 1 + len(PTABLE)), dtype=complex)
    R = np.zeros((len(confs), len(GAMMAS), NT, 1 + len(PTABLE)), dtype=complex)
    I = np.zeros((len(confs), len(GAMMAS), NT, 1 + len(PTABLE)), dtype=complex)
    good = np.ones(len(confs), dtype=bool)
    for ic, conf in enumerate(confs):
        for ig, g in enumerate(GAMMAS):
            d = read_trace(os.path.join(ens_dir, f"traces.{g}.{conf}"))
            if d is None:
                print(f"  truncated: traces.{g}.{conf}")
                good[ic] = False
                continue
            # L[ic, ig] = project(d, expk)
            R[ic, ig] = project(d.real.astype(complex), expk)
            I[ic, ig] = project(d.imag.astype(complex), expk)
        if ic % 20 == 0:
            print(f"  {ic}/{len(confs)} conf {conf}", flush=True)

    out = os.path.join(out_dir, f"loops_{ens}_claude.npz")
    # np.savez(out, confs=np.array(confs)[good], gammas=np.array(GAMMAS), L=L[good])
    np.savez(out, confs=np.array(confs)[good], gammas=np.array(GAMMAS), R=R[good], I=I[good])
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

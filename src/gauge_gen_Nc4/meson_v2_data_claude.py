#!/usr/bin/env python3
# meson_v2_data_claude.py
#
# Per-ensemble data file for the corrected meson analysis (v2):
#   - connected C from ALL mesons_conn.<conf>.h5 (folded, p^2-class averaged;
#     meson_combine_claude.connected_pclass)
#   - disconnected D from the complex loop cache (disc_loop_cache_claude.py),
#     built from the SIGNAL part with the correct sign
#     (disc_corr_from_loops_claude.disc_corr_gamma: real loops +<r r>,
#     vector loops -<l l>); NO subtraction of any kind
#   - Lid[conf] = (1/NT) sum_t (slice-averaged id loop), for the vacuum
#     subtraction N_S^3 <Lid>^2 of the scalar disc (done in the analysis)
#
# Output: meson_v2_out_claude/meson_v2_<ens>_claude.npz with
#   channels (10,), gammas (10,),
#   confs_C (NC,), C (10, NC, 5, NT/2+1)    [channel order = CHANNELS]
#   confs_D (ND,), D (10, ND, 5, NT/2+1)    [gamma order = GAMMAS]
#   Lid (ND,), beta, mass
#
# Usage:
#   OMP_NUM_THREADS=4 python3 meson_v2_data_claude.py [ens_name ...]

import sys
import os
import glob
import re
import numpy as np
import meson_combine_claude as mc
import disc_corr_from_loops_claude as dc

OBS_ROOT = "/mnt/baracuda_14/grid_claude/obs_nc4nf1_2448"
LOOP_DIR = "disc_loops_claude"
OUT_DIR = "meson_v2_out_claude"

CHANNELS = ["I_I", "G5_G5", "GTG5_GTG5", "GXG5_GXG5", "GYG5_GYG5", "GZG5_GZG5",
            "GT_GT", "GX_GX", "GY_GY", "GZ_GZ"]
GAMMAS = ["id", "g5", "gtg5", "gxg5", "gyg5", "gzg5", "gt", "gx", "gy", "gz"]


def read_connected_all(ens_dir):
    files = glob.glob(os.path.join(ens_dir, "mesons_conn.*.h5"))
    pairs = []
    for f in files:
        m = re.search(r"mesons_conn\.(\d+)\.h5$", f)
        if m:
            pairs.append((int(m.group(1)), f))
    pairs.sort()
    confs = []
    rows = []
    for conf, f in pairs:
        try:
            per_ch = [mc.connected_pclass(f, ch) for ch in CHANNELS]
        except (OSError, KeyError) as e:
            print(f"  WARNING: skipping {os.path.basename(f)}: {e}")
            continue
        confs.append(conf)
        rows.append(per_ch)
    C = np.array(rows)                       # (NC, 10, 5, nt2+1)
    return np.array(confs, dtype=int), np.transpose(C, (1, 0, 2, 3))


def build(ens, lma=False):
    # lma=True: disc from the LMA-estimator loop cache loops_lma_nc4nf1_2448_<b_m>_claude.npz
    # (light ensembles), output tagged <ens>_lma; connected always from OBS_ROOT/<ens>.
    m = re.search(r"b(\d+p\d+)_m(\d+p\d+)$", ens)
    beta = float(m.group(1).replace("p", "."))
    mass = float(m.group(2).replace("p", "."))

    confs_C, C = read_connected_all(os.path.join(OBS_ROOT, ens))
    print(f"{ens}: connected {len(confs_C)} configs")

    loop_name = f"loops_{ens}_claude.npz"
    out_tag = ens
    if lma:
        loop_name = "loops_lma_nc4nf1_2448_" + ens.replace("obs_nc4nf1_2448_", "") + "_claude.npz"
        out_tag = ens + "_lma"
    # confs_D, gammas, R, I = dc.load_loops(os.path.join(LOOP_DIR, f"loops_{ens}_claude.npz"))
    confs_D, gammas, R, I = dc.load_loops(os.path.join(LOOP_DIR, loop_name))
    D = np.array([dc.disc_corr_gamma(R, I, gammas, g) for g in GAMMAS])
    Lid = R[:, gammas.index("id"), :, 0].real.mean(axis=1)
    print(f"{ens}: disc {len(confs_D)} configs")

    os.makedirs(OUT_DIR, exist_ok=True)
    # out = os.path.join(OUT_DIR, f"meson_v2_{ens}_claude.npz")
    out = os.path.join(OUT_DIR, f"meson_v2_{out_tag}_claude.npz")
    np.savez(out, channels=np.array(CHANNELS), gammas=np.array(GAMMAS),
             confs_C=confs_C, C=C, confs_D=np.array(confs_D, dtype=int), D=D,
             Lid=Lid, beta=beta, mass=mass)
    print(f"wrote {out}")


def main():
    # usage: meson_v2_data_claude.py [--lma] [ens_name ...]
    args = sys.argv[1:]
    lma = False
    if len(args) > 0 and args[0] == "--lma":
        lma = True
        args = args[1:]
    ens_list = args
    if len(ens_list) == 0:
        ens_list = sorted(os.path.basename(f)[len("loops_"):-len("_claude.npz")]
                          for f in glob.glob(os.path.join(LOOP_DIR, "loops_obs_*_claude.npz")))
    for ens in ens_list:
        build(ens, lma)


if __name__ == "__main__":
    main()

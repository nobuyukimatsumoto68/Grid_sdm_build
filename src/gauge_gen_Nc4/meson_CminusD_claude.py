#!/usr/bin/env python3
# meson_CminusD_claude.py
#
# End-to-end physical flavor-singlet meson correlator with jackknife errors:
#
#     C_phys(p_class, t) = C_conn(p_class, t) - Nf * rel_norm * D(p_class, t)
#
# (Nf=1 here; the overall sign is irrelevant for masses, so C-D is used to match
# the file name. See meson_singlet_sign_derivation_claude.md.)
#
# D : disconnected loop-loop correlator, Rebbi's corr_d.p{0..4} (per config,
#     already folded + p^2-class averaged).
# C : connected meson correlator from mesons_conn.<idx>.h5 (folded + class
#     averaged here by meson_combine_claude.connected_pclass).
#
# Configs are matched by index between the disc files and the connected h5
# files; only the common set is used. Errors are leave-one-out jackknife over
# configs (optional binning).
#
# rel_norm defaults to 1.0. Per OQ1 (see meson_momproj_impl_plan_claude.md), the
# factor-3 in makecorr_96_d.f90 is a legacy QCD inter-ensemble artifact, NOT a
# color factor, so for this SU(4) N_f=1 project rel_norm = 1 is expected; the
# p=0 behaviour is the check.
#
# Usage:
#   python3 meson_CminusD_claude.py <ensemble_dir> [--channel G5_G5]
#       [--rel-norm 1.0] [--binsize 1] [--out PREFIX]
# Example:
#   python3 meson_CminusD_claude.py \
#     ../../../obs_nc4nf1_2448/obs_nc4nf1_2448_b11p045_m0p4000 --channel G5_G5

import os
import sys
import argparse
import numpy as np

import meson_combine_claude as mc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ensemble_dir", help="path to obs_nc4nf1_2448_b..._m... dir")
    ap.add_argument("--channel", default="G5_G5",
                    choices=list(mc.CHANNEL_TO_TRACEDIR.keys()))
    ap.add_argument("--rel-norm", type=float, default=1.0,
                    help="relative D-vs-C normalization (OQ1, default 1.0)")
    ap.add_argument("--binsize", type=int, default=1,
                    help="jackknife bin size (consecutive configs)")
    ap.add_argument("--out", default=None,
                    help="output prefix; default derived from ensemble+channel")
    args = ap.parse_args()

    ens_dir = args.ensemble_dir.rstrip("/")
    channel = args.channel
    tracedir = mc.CHANNEL_TO_TRACEDIR[channel]

    # --- connected: (cidx, C[nconf,5,NT/2+1]) ---
    cidx, C_all = mc.collect_connected(ens_dir, channel)
    if cidx.size == 0:
        print(f"ERROR: no mesons_conn.*.h5 in {ens_dir}", file=sys.stderr)
        sys.exit(1)
    nh = C_all.shape[-1]                      # NT/2 + 1

    # --- disconnected: one file per p-class, build D[5, nconf_d, NT/2+1] ---
    # indices may differ per class file in principle; intersect per class.
    D_by_class = {}
    didx_by_class = {}
    for m in range(5):
        dfile = os.path.join(ens_dir, tracedir, f"corr_d.p{m}")
        didx, D = mc.read_corr_d_indexed(dfile)
        D_by_class[m] = D
        didx_by_class[m] = didx

    # Common configs across connected and all disc class files.
    common = set(cidx.tolist())
    for m in range(5):
        common &= set(didx_by_class[m].tolist())
    common = np.array(sorted(common), dtype=int)
    if common.size == 0:
        print("ERROR: no common configs between connected and disc",
              file=sys.stderr)
        sys.exit(1)

    print(f"channel={channel}  tracedir={tracedir}")
    print(f"connected configs={cidx.size}  common configs={common.size}  "
          f"binsize={args.binsize}  rel_norm={args.rel_norm}")

    # index lookups
    cpos = {idx: i for i, idx in enumerate(cidx.tolist())}

    # aligned per-config arrays, shape (5, ncommon, NT/2+1)
    nc = common.size
    D_al = np.zeros((5, nc, nh))
    C_al = np.zeros((5, nc, nh))
    for m in range(5):
        dpos = {idx: i for i, idx in enumerate(didx_by_class[m].tolist())}
        for j, idx in enumerate(common):
            D_al[m, j] = D_by_class[m][dpos[idx]]
            C_al[m, j] = C_all[cpos[idx], m]

    # physical C - D per config, then jackknife (Nf=1; sign irrelevant).
    phys = C_al - args.rel_norm * D_al              # (5, nc, NT/2+1)

    out_prefix = args.out
    if out_prefix is None:
        ens_name = os.path.basename(ens_dir)
        out_prefix = f"CminusD_{ens_name}_{channel}"

    ts = np.arange(nh)
    phys_mean = np.zeros((5, nh))
    phys_err = np.zeros((5, nh))
    D_mean = np.zeros((5, nh))
    D_err = np.zeros((5, nh))
    C_mean = np.zeros((5, nh))
    C_err = np.zeros((5, nh))
    for m in range(5):
        phys_mean[m], phys_err[m] = mc.jackknife(phys[m], args.binsize)
        D_mean[m], D_err[m] = mc.jackknife(D_al[m], args.binsize)
        C_mean[m], C_err[m] = mc.jackknife(C_al[m], args.binsize)

    # save arrays for downstream analysis / plotting
    npz = out_prefix + ".npz"
    np.savez(npz, t=ts, common=common, rel_norm=args.rel_norm,
             binsize=args.binsize,
             phys_mean=phys_mean, phys_err=phys_err,
             D_mean=D_mean, D_err=D_err, C_mean=C_mean, C_err=C_err)
    print(f"saved {npz}")

    # readable table
    for m in range(5):
        print(f"\n# p^2 class p{m}:  t   C-D         err          "
              f"D           err          C           err")
        for t in range(nh):
            print(f"  {t:3d}  {phys_mean[m,t]: .6e} {phys_err[m,t]: .3e}  "
                  f"{D_mean[m,t]: .6e} {D_err[m,t]: .3e}  "
                  f"{C_mean[m,t]: .6e} {C_err[m,t]: .3e}")


if __name__ == "__main__":
    main()

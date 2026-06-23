#!/usr/bin/env python3
# plot_CminusD_h5_claude.py
#
# Read a CminusD_<ens>.h5 (datasets "<conf>/<ch>/<mom>" = C - Nf*D_sub, the
# singlet correlator per config; root keys beta, m; attrs Nf, Nc) and plot the
# jackknifed C-D vs t for one channel / p^2 class.
#
# Usage:
#   python3 plot_CminusD_h5_claude.py CminusD_h5_out_claude/CminusD_<ens>.h5 \
#       [--ch G5_G5] [--mom p0] [--symlog] [--save out.pdf]

import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt


def jackknife(a):
    # a : (nconf, nt). Leave-one-out mean and error per t.
    n = a.shape[0]
    mean = a.mean(0)
    loo = (a.sum(0) - a) / (n - 1)
    err = np.sqrt((n - 1) / n * ((loo - loo.mean(0)) ** 2).sum(0))
    return mean, err


def load_CminusD(h5file, ch, mom):
    rows = []
    with h5py.File(h5file, "r") as h:
        beta = float(h["beta"][()])
        mass = float(h["m"][()])
        nc = int(h.attrs.get("Nc", 0))
        confs = sorted((k for k in h.keys() if k.isdigit()), key=int)
        for c in confs:
            key = f"{c}/{ch}/{mom}"
            if key in h:
                rows.append(np.array(h[key]))
    data = np.array(rows)            # (nconf, nt)
    return data, beta, mass, nc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5file")
    ap.add_argument("--ch", default="G5_G5")
    ap.add_argument("--mom", default="p0")
    ap.add_argument("--symlog", action="store_true",
                    help="symlog y-axis (shows the zero crossing)")
    ap.add_argument("--save", default=None, help="save figure to this path")
    args = ap.parse_args()

    data, beta, mass, nc = load_CminusD(args.h5file, args.ch, args.mom)
    if data.size == 0:
        raise SystemExit(f"no datasets for ch={args.ch} mom={args.mom} in {args.h5file}")

    nconf = data.shape[0]
    mean, err = jackknife(data)
    t = np.arange(data.shape[1])

    plt.errorbar(t, mean, err, marker='o', markersize=5)
    plt.hlines(0.0, t[0], t[-1])
    if args.symlog:
        plt.yscale("symlog", linthresh=1e-5)
    plt.xlabel("$t$")
    plt.ylabel("$C-D$  (singlet)")
    plt.title(f"{args.ch}  {args.mom}   "
              f"$\\beta$={beta}  $m$={mass}  $N_c$={nc}  $N_{{\\rm conf}}$={nconf}")
    if args.save:
        plt.savefig(args.save, bbox_inches="tight")
        print(f"saved {args.save}")
    plt.show()


if __name__ == "__main__":
    main()

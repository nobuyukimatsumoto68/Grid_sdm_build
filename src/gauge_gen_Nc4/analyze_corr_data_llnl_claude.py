#!/usr/bin/env python3
# analyze_corr_data_llnl_claude.py
#
# Jackknife analysis of the per-ensemble correlator h5 files in corr_data.tar.gz
# (pulled from LLNL by rsync_pull_corr_data_claude.sh). Layout of each
# corr_<ens>.h5 (written by meson_CminusD_h5_claude.py):
#   <contr>/<ch>/<mom>  shape (Nconf, NT/2+1), folded, one row per config
#   confs/<ch>          config numbers aligned with the rows
#   contr in {C, Dsub, CminusDsub}, mom in {p0..p4} (p^2 class)
#
# For every (ensemble, contr, channel, mom):
#   - jackknife mean/err of the correlator
#   - cosh effective mass  m_eff(t) = arccosh[(C(t-1)+C(t+1)) / (2 C(t))],
#     exact for a single cosh(m(t-T/2)); nan where the ratio < 1
#   - plateau mass = weighted constant fit of m_eff over [tmin, tmax]
#     (weights 1/err^2 from the full jackknife, fixed across samples)
# Jackknife: M. H. Quenouille, Biometrika 43 (1956) 353; J. W. Tukey (1958).
#
# Outputs (OUT_DIR):
#   meff_<ens>_claude.npz   arrays [contr, ch, mom, t] + plateau [contr, ch, mom]
#   plateau_summary_claude.txt
#
# Usage:
#   OMP_NUM_THREADS=4 python3 analyze_corr_data_llnl_claude.py [data_dir] [out_dir]

import sys
import os
import glob
import numpy as np
import h5py

DATA_DIR = "/mnt/baracuda_14/grid_claude/corr_data_llnl_claude/extract/corr_h5_out_claude"
OUT_DIR = "corr_llnl_analysis_claude"

CONTRS = ["C", "Dsub", "CminusDsub"]
CHANNELS = ["G5_G5", "GTG5_GTG5", "GXG5_GXG5", "GYG5_GYG5", "GZG5_GZG5",
            "I_I", "GT_GT", "GX_GX", "GY_GY", "GZ_GZ"]
MOMS = ["p0", "p1", "p2", "p3", "p4"]

# Bin consecutive configs before the jackknife (1 = no binning).
BINSIZE = 1

# Plateau windows [tmin, tmax] (inclusive, in m_eff time t). By eye from the
# p0 survey: connected plateaus from t~12; the singlet loses signal at t~6-8,
# so its window is SHORT and is not a true plateau. Dsub alone is not a
# single-exponential -> no fit (None).
WINDOWS = {
    "C": (12, 20),
    "Dsub": None,
    "CminusDsub": (4, 7),
}
# Per-(ensemble, contr, channel) overrides, e.g.
#   WINDOW_OVERRIDE[("obs_nc4nf1_2448_b11p045_m0p4000", "CminusDsub", "G5_G5")] = (3, 5)
WINDOW_OVERRIDE = {}


def jackknife_samples(x, binsize):
    # Leave-one-out means over axis 0 after binning.
    n = x.shape[0]
    nb = n // binsize
    xb = x[:nb * binsize].reshape((nb, binsize) + x.shape[1:]).mean(axis=1)
    total = xb.sum(axis=0)
    jk = (total - xb) / (nb - 1)
    return jk


def jk_mean_err(jk):
    # Jackknife mean and error of samples on axis 0, ignoring nan samples.
    nb = jk.shape[0]
    mean = np.nanmean(jk, axis=0)
    err = np.sqrt((nb - 1) / nb * np.nansum((jk - mean) ** 2, axis=0))
    return mean, err


def cosh_meff(c):
    # c: (..., nt2+1) folded correlator. Returns (..., nt2+1) with m_eff at
    # t=1..nt2-1 and nan at the two ends.
    meff = np.full(c.shape, np.nan)
    ratio = (c[..., :-2] + c[..., 2:]) / (2.0 * c[..., 1:-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        meff[..., 1:-1] = np.arccosh(ratio)
    return meff


def plateau_fit(meff_jk, meff_err, tmin, tmax):
    # Weighted constant fit of m_eff over [tmin, tmax] in each jackknife sample.
    w = 1.0 / meff_err[tmin:tmax + 1] ** 2
    w[~np.isfinite(w)] = 0.0
    seg = meff_jk[:, tmin:tmax + 1]
    good = np.isfinite(seg)
    wsum = np.where(good, w, 0.0).sum(axis=1)
    num = np.where(good, w * np.nan_to_num(seg), 0.0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        fit = num / wsum
    nan_frac = 1.0 - good.mean()
    return fit, nan_frac


def analyze_file(fn):
    h5 = h5py.File(fn, "r")
    ens = os.path.basename(fn)[len("corr_"):-len(".h5")]
    beta = float(h5["beta"][()])
    mass = float(h5["m"][()])
    nt2p1 = h5["C/G5_G5/p0"].shape[1]

    shape4 = (len(CONTRS), len(CHANNELS), len(MOMS), nt2p1)
    shape3 = (len(CONTRS), len(CHANNELS), len(MOMS))
    corr_mean = np.full(shape4, np.nan)
    corr_err = np.full(shape4, np.nan)
    meff_mean = np.full(shape4, np.nan)
    meff_err = np.full(shape4, np.nan)
    plat_mean = np.full(shape3, np.nan)
    plat_err = np.full(shape3, np.nan)
    plat_nanfrac = np.full(shape3, np.nan)
    plat_win = np.full(shape3 + (2,), -1, dtype=int)
    nconf = np.zeros(len(CHANNELS), dtype=int)

    for ic, contr in enumerate(CONTRS):
        for ich, ch in enumerate(CHANNELS):
            nconf[ich] = h5[f"confs/{ch}"].shape[0]
            win = WINDOW_OVERRIDE.get((ens, contr, ch), WINDOWS[contr])
            for im, mom in enumerate(MOMS):
                x = h5[f"{contr}/{ch}/{mom}"][()]
                jk = jackknife_samples(x, BINSIZE)
                cm, ce = jk_mean_err(jk)
                corr_mean[ic, ich, im] = cm
                corr_err[ic, ich, im] = ce

                mjk = cosh_meff(jk)
                mm, me = jk_mean_err(mjk)
                meff_mean[ic, ich, im] = mm
                meff_err[ic, ich, im] = me

                if win is None:
                    continue
                tmin, tmax = win
                fit, nan_frac = plateau_fit(mjk, me, tmin, tmax)
                fm, fe = jk_mean_err(fit)
                plat_mean[ic, ich, im] = fm
                plat_err[ic, ich, im] = fe
                plat_nanfrac[ic, ich, im] = nan_frac
                plat_win[ic, ich, im] = (tmin, tmax)
    h5.close()

    res = {
        "ensemble": ens,
        "beta": beta,
        "mass": mass,
        "contrs": np.array(CONTRS),
        "channels": np.array(CHANNELS),
        "moms": np.array(MOMS),
        "nconf": nconf,
        "binsize": BINSIZE,
        "corr_mean": corr_mean,
        "corr_err": corr_err,
        "meff_mean": meff_mean,
        "meff_err": meff_err,
        "plat_mean": plat_mean,
        "plat_err": plat_err,
        "plat_nanfrac": plat_nanfrac,
        "plat_win": plat_win,
    }
    return res


def summary_lines(res):
    lines = []
    lines.append(f"# {res['ensemble']}  beta={res['beta']}  m={res['mass']}  binsize={res['binsize']}")
    lines.append(f"#   {'channel':10s} {'Ncf':>4s} | " + " | ".join(
        f"{c + ' p0 [win]':>28s}" for c in CONTRS if WINDOWS[c] is not None))
    for ich, ch in enumerate(CHANNELS):
        cols = []
        for ic, contr in enumerate(CONTRS):
            if WINDOWS[contr] is None:
                continue
            m = res["plat_mean"][ic, ich, 0]
            e = res["plat_err"][ic, ich, 0]
            tmin, tmax = res["plat_win"][ic, ich, 0]
            nf = res["plat_nanfrac"][ic, ich, 0]
            flag = "*" if nf > 0.0 else " "
            cols.append(f"{m:9.4f} +- {e:7.4f} [{tmin:2d},{tmax:2d}]{flag}")
        lines.append(f"    {ch:10s} {res['nconf'][ich]:4d} | " + " | ".join(f"{c:>28s}" for c in cols))
    lines.append("")
    return lines


def main():
    data_dir = DATA_DIR
    out_dir = OUT_DIR
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    if len(sys.argv) > 2:
        out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(data_dir, "corr_*.h5")))
    if len(files) == 0:
        print(f"no corr_*.h5 in {data_dir}")
        sys.exit(1)

    all_lines = []
    all_lines.append("# plateau masses from cosh m_eff, jackknife errors; p0 class")
    all_lines.append("# '*' = some jackknife samples had nan m_eff inside the window (signal lost)")
    all_lines.append("# CminusDsub window is short: NOT a true plateau, an upper-bound-like estimate")
    all_lines.append("")
    for fn in files:
        res = analyze_file(fn)
        np.savez(os.path.join(out_dir, f"meff_{res['ensemble']}_claude.npz"), **res)
        lines = summary_lines(res)
        all_lines.extend(lines)
        print("\n".join(lines))

    with open(os.path.join(out_dir, "plateau_summary_claude.txt"), "w") as f:
        f.write("\n".join(all_lines) + "\n")
    print(f"wrote {out_dir}/meff_*_claude.npz and plateau_summary_claude.txt")


if __name__ == "__main__":
    main()

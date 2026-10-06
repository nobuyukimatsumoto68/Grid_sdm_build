#!/usr/bin/env python3
# meson_v2_png_claude.py
#
# PNG version of meson_v2_meff_claude.ipynb: for each ensemble, the 7 plots
#   1 G5 (eta')      C / D raw / C - D  m_eff
#   2 GTG5 (eta')    C / D / C - D      m_eff
#   3 I_I (sigma)    C / D vac-sub / C - D
#   4 GiG5 (f1)      C xyz-avg / D gzg5 / C - D
#   5 Gi (omega)     C xyz-avg / D xyz-avg (Im part) / C - D
#   6 disc significance |D|/sigma_D for all 10 loops
#   7 raw g5 disc D(t) to T/2
# Same data and conventions as the notebook (see its header cell).
# Jackknife: M. H. Quenouille, Biometrika 43 (1956) 353; J. W. Tukey (1958).
#
# Usage:
#   OMP_NUM_THREADS=4 python3 meson_v2_png_claude.py [psq] [out_dir]

import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NS = 24
BIN = 1
# TMAX = 16
TMAX = 24
YMAX = 2.6

# Large-t plateau windows [tlo, thi] of the RAW folded disc D(t), p0, picked BY EYE
# from the linear plots (_8_g5_disc_tsumsub, _9_id_disc_tsumsub) + printed tables,
# 2026-10-06. None = no plateau identifiable (light ensembles, plain-disc N=12/24).
#   m0.1 g5: dip to 6e-5 at t=13-14 then flat ~1.7e-4 for t=19-24 (ambiguous)
#   m0.2 g5: flat ~3e-5 from t=6;  m0.3/m0.4 g5: ~0 from t=3-4
#   id: vac-sub D flat (~0) from t=10 (m0.1), t=5 (m0.2-0.4)
# Plateau windows [tlo, thi] for the constant fit of the CONNECTED cosh m_eff
# (plots 14/15), m0.2-0.4 only, picked BY EYE from the plots + m_eff tables 2026-10-06.
CONN_FIT_WIN = {
    ("b10p990_m0p2000", "G5_G5"): (19, 23),
    ("b10p990_m0p2000", "I_I"): (15, 20),
    ("b11p035_m0p3000", "G5_G5"): (18, 23),
    ("b11p035_m0p3000", "I_I"): (14, 21),
    ("b11p045_m0p4000", "G5_G5"): (20, 23),
    ("b11p045_m0p4000", "I_I"): (17, 22),
}

# (tag, channel) -> (quark mass, plateau m_eff, error), filled by the CONN_FIT_WIN fits
FIT_RESULTS = {}

# Single-baryon 2^+ / 2^- GEVP ground states (lattice units), from session grid-claude-54
# (2026-10-06, read-only recompute): bvar_h5_out_claude/bvar_<b_m>_ss.h5 group Copss,
# bvar_spectrum_claude.py; smeared source + smeared sink (Gaussian w=3, N=40);
# 2+: 3x3 GEVP (2p_40, 2p_22, 2p_04), 2-: 2x2 GEVP (2m_31, 2m_13), t0=3, log GEVP m_eff,
# uncorrelated const fit, windows 2+_0 [7,14], 2-_0 [5,14] (chosen for m0.4; m0.1-0.3 PRELIMINARY).
# rows: (m, m_2+, err, m_2-, err)
BARYON_2PM = [
    (0.1, 1.4868, 0.0122, 1.7249, 0.0170),
    (0.2, 1.9920, 0.0076, 2.1739, 0.0132),
    (0.3, 2.4678, 0.0094, 2.6336, 0.0155),
    (0.4, 2.9298, 0.0075, 3.0847, 0.0110),
]

PLAT_WIN = {
    ("b10p800_m0p0100", "g5"): None,
    ("b10p800_m0p0100", "id"): None,
    ("b10p840_m0p0500", "g5"): None,
    ("b10p840_m0p0500", "id"): None,
    ("b10p865_m0p1000", "g5"): (19, 24),
    ("b10p865_m0p1000", "id"): (10, 24),
    ("b10p990_m0p2000", "g5"): (6, 24),
    ("b10p990_m0p2000", "id"): (5, 24),
    ("b11p035_m0p3000", "g5"): (4, 24),
    ("b11p035_m0p3000", "id"): (5, 24),
    ("b11p045_m0p4000", "g5"): (4, 24),
    ("b11p045_m0p4000", "id"): (5, 24),
    # m0.01 LMA disc (160 cfgs): no plateau seen; user: use the last two points
    ("b10p800_m0p0100_lma", "g5"): (23, 24),
    ("b10p800_m0p0100_lma", "id"): (23, 24),
    # m0.05 LMA disc (208 cfgs): same as m0.01, last two points
    ("b10p840_m0p0500_lma", "g5"): (23, 24),
    ("b10p840_m0p0500_lma", "id"): (23, 24),
}


def jk_samples(x, binsize):
    nb = x.shape[0] // binsize
    xb = x[:nb * binsize].reshape((nb, binsize) + x.shape[1:]).mean(axis=1)
    return (xb.sum(axis=0) - xb) / (nb - 1)


def jk_mean_err(s):
    nb = s.shape[0]
    m = s.mean(axis=0)
    e = np.sqrt((nb - 1) / nb * ((s - m) ** 2).sum(axis=0))
    return m, e


def meff_jk(s):
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (s[:, :-2] + s[:, 2:]) / (2.0 * s[:, 1:-1])
        me = np.arccosh(r)
    bad = ~np.all(np.isfinite(me), axis=0)
    m, e = jk_mean_err(np.where(np.isfinite(me), me, 0.0))
    m[bad] = np.nan
    e[bad] = np.nan
    t = np.arange(1, s.shape[1] - 1)
    return t, m, e


def const_fit_meff_jk(s, tlo, thi):
    # uncorrelated weighted constant fit of the cosh m_eff over [tlo, thi] in each
    # jackknife sample; weights 1/err^2 from the full jackknife (fixed across samples)
    with np.errstate(invalid="ignore", divide="ignore"):
        me = np.arccosh((s[:, :-2] + s[:, 2:]) / (2.0 * s[:, 1:-1]))
    seg = me[:, tlo - 1:thi]
    m, e = jk_mean_err(seg)
    w = 1.0 / e**2
    fit = (seg * w).sum(axis=1) / w.sum()
    return jk_mean_err(fit)


def meff_log_jk(s):
    # log m_eff = ln[G(t)/G(t+1)] on each jackknife sample; nan at t where any
    # sample is undefined (ratio <= 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        me = np.log(s[:, :-1] / s[:, 1:])
    bad = ~np.all(np.isfinite(me), axis=0)
    m, e = jk_mean_err(np.where(np.isfinite(me), me, 0.0))
    m[bad] = np.nan
    e[bad] = np.nan
    t = np.arange(0, s.shape[1] - 1)
    return t, m, e


def plot_three(fn, title, r, labC, labD, nC, nD):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    t, m, e = r["C"]
    ax.errorbar(t, m, e, fmt="s", color="tab:blue", ms=5, capsize=2, label=f"{labC} (N={nC})")
    t, m, e = r["D"]
    ax.errorbar(t + 0.1, m, e, fmt="^", color="tab:orange", ms=5, capsize=2, label=f"{labD} (N={nD})")
    t, m, e = r["CmD"]
    ax.errorbar(t + 0.2, m, e, fmt="o", color="tab:red", ms=5, capsize=2, label="C - D")
    ax.set_xlim(0, TMAX + 0.5)
    ax.set_ylim(0, YMAX)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$m_\mathrm{eff}$")
    ax.set_title(title)
    ax.legend()
    fig.savefig(fn, dpi=120, bbox_inches="tight")
    plt.close(fig)


def do_ensemble(npz, psq, out_dir):
    z = np.load(npz)
    ens = os.path.basename(npz).replace("meson_v2_", "").replace("_claude.npz", "")
    tag = ens.replace("obs_nc4nf1_2448_", "")
    beta = float(z["beta"])
    mass = float(z["mass"])
    channels = list(z["channels"])
    gammas = list(z["gammas"])
    confs_C = z["confs_C"]
    confs_D = z["confs_D"]
    common = np.intersect1d(confs_C, confs_D)
    iC = np.searchsorted(confs_C, common)
    iD = np.searchsorted(confs_D, common)
    Lid = z["Lid"]

    def Cof(ch):
        return z["C"][channels.index(ch)][:, psq]

    def Dof(g):
        return z["D"][gammas.index(g)][:, psq]

    groups = {
        "1_G5_etaprime": ("G5 (eta')", Cof("G5_G5"), Dof("g5"), "C", "D raw", False),
        "2_GTG5_etaprime": ("GTG5 (eta')", Cof("GTG5_GTG5"), Dof("gtg5"), "C", "D", False),
        "3_II_sigma": ("I_I (sigma)", Cof("I_I"), Dof("id"), "C", "D vac-sub", True),
        "4_GiG5_f1": ("GiG5 (f1), D from gzg5",
                      (Cof("GXG5_GXG5") + Cof("GYG5_GYG5") + Cof("GZG5_GZG5")) / 3.0,
                      Dof("gzg5"), "C xyz-avg", "D gzg5", False),
        "5_Gi_omega": ("Gi (omega)",
                       (Cof("GX_GX") + Cof("GY_GY") + Cof("GZ_GZ")) / 3.0,
                       (Dof("gx") + Dof("gy") + Dof("gz")) / 3.0, "C xyz-avg", "D xyz-avg", False),
    }
    tail = rf"$\beta$={beta}, $m$={mass}, $p^2$ class {psq}"
    edir = os.path.join(out_dir, tag)
    os.makedirs(edir, exist_ok=True)

    sD_g5 = None
    for key, (name, C, D, labC, labD, vacsub) in groups.items():
        sCm = jk_samples(C[iC], BIN)
        sDm = jk_samples(D[iD], BIN)
        if vacsub and psq == 0:
            sDm = sDm - NS**3 * (jk_samples(Lid[iD], BIN) ** 2)[:, None]
        if key.startswith("1_"):
            sD_g5 = sDm
        r = {"C": meff_jk(jk_samples(C, BIN)), "D": meff_jk(sDm), "CmD": meff_jk(sCm - sDm)}
        plot_three(os.path.join(edir, f"{tag}_p{psq}_{key}_claude.png"), f"{name}: " + tail, r,
                   labC, labD, len(confs_C), len(common))

    # disc significance, all loops
    markers = ["o", "s", "^", "v", "D", "P", "X", "<", ">", "*"]
    colors = plt.cm.tab10(np.arange(10))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    tt = np.arange(1, TMAX + 1)
    for k, g in enumerate(gammas):
        s = jk_samples(Dof(g)[iD], BIN)
        if g == "id" and psq == 0:
            s = s - NS**3 * (jk_samples(Lid[iD], BIN) ** 2)[:, None]
        m, e = jk_mean_err(s)
        ax.semilogy(tt, np.abs(m[tt]) / e[tt], marker=markers[k], color=colors[k], ls="-", lw=0.8, ms=5, label=g)
    ax.axhline(2.0, color="gray", ls="--", lw=1)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$|\bar D| / \sigma_D$")
    ax.set_title(f"disc significance (N={len(common)}): " + tail)
    ax.legend(ncol=2, fontsize=8)
    fig.savefig(os.path.join(edir, f"{tag}_p{psq}_6_disc_significance_claude.png"), dpi=120, bbox_inches="tight")
    plt.close(fig)

    # raw g5 disc correlator to T/2
    m, e = jk_mean_err(sD_g5)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.errorbar(np.arange(len(m)), m, e, fmt="^", color="tab:orange", ms=5, capsize=2, label=f"D raw g5 (N={len(common)})")
    ax.axhline(0.0, color="gray", lw=0.8)
    ax.set_xlabel("t")
    ax.set_ylabel("D(t)")
    ax.set_title(r"raw $\gamma_5$ disc: " + tail)
    ax.legend()
    fig.savefig(os.path.join(edir, f"{tag}_p{psq}_7_g5_disc_raw_claude.png"), dpi=120, bbox_inches="tight")
    plt.close(fig)

    # t-sum-subtracted disc (old h5 "Dsub"), linear scale, vs raw / vac-sub
    #   D_sub(t) = D(t) - (1/NT) sum_t' D(t'),  folded: avg = (D0 + D_{NT/2} + 2 sum_{1}^{NT/2-1} D)/NT
    nt2 = sD_g5.shape[1] - 1
    for key, g, lab0, vacsub in [("8_g5_disc_tsumsub", "g5", "D raw", False),
                                 ("9_id_disc_tsumsub", "id", "D vac-sub", True)]:
        Dm = Dof(g)[iD]
        avg = (Dm[:, 0] + Dm[:, nt2] + 2.0 * Dm[:, 1:nt2].sum(axis=1)) / (2 * nt2)
        s0 = jk_samples(Dm, BIN)
        if vacsub and psq == 0:
            s0 = s0 - NS**3 * (jk_samples(Lid[iD], BIN) ** 2)[:, None]
        s1 = jk_samples(Dm - avg[:, None], BIN)
        m0, e0 = jk_mean_err(s0)
        m1, e1 = jk_mean_err(s1)
        tt = np.arange(nt2 + 1)
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.errorbar(tt, m0, e0, fmt="^", color="tab:orange", ms=5, capsize=2, label=f"{lab0} (N={len(common)})")
        ax.errorbar(tt + 0.15, m1, e1, fmt="D", color="tab:purple", ms=4, capsize=2, label="D t-sum-subtracted")
        ax.axhline(0.0, color="gray", lw=0.8)
        ax.set_xlabel("t")
        ax.set_ylabel("D(t)")
        ax.set_title(rf"{g} disc, linear: " + tail)
        ax.legend()
        fig.savefig(os.path.join(edir, f"{tag}_p{psq}_{key}_claude.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)

    # PLATEAU-subtracted disc, linear scale: per jackknife sample subtract the
    # large-t plateau c = mean_{t in PLAT_WIN} D(t) of the RAW folded D (p0 only).
    #   10: g5   raw / t-sum-sub / plateau-sub
    #   11: id   vac-sub (N_S^3 <Lbar>^2) / plateau-sub of raw D
    for key, g in [("10_g5_disc_platsub", "g5"), ("11_id_disc_platsub", "id")]:
        win = PLAT_WIN.get((tag, g))
        if win is None or psq != 0:
            print(f"  {tag} {g}: no plateau window -> skip {key}")
            continue
        tlo, thi = win
        Dm = Dof(g)[iD]
        sraw = jk_samples(Dm, BIN)
        cplat = sraw[:, tlo:thi + 1].mean(axis=1)
        splat = sraw - cplat[:, None]
        cm, ce = jk_mean_err(cplat)
        tt = np.arange(nt2 + 1)
        fig, ax = plt.subplots(figsize=(7, 4.5))
        if g == "g5":
            avg = (Dm[:, 0] + Dm[:, nt2] + 2.0 * Dm[:, 1:nt2].sum(axis=1)) / (2 * nt2)
            m0, e0 = jk_mean_err(sraw)
            m1, e1 = jk_mean_err(jk_samples(Dm - avg[:, None], BIN))
            ax.errorbar(tt, m0, e0, fmt="^", color="tab:orange", ms=5, capsize=2, label=f"D raw (N={len(common)})")
            ax.errorbar(tt + 0.15, m1, e1, fmt="D", color="tab:purple", ms=4, capsize=2, label="D t-sum-subtracted")
        else:
            svac = sraw - NS**3 * (jk_samples(Lid[iD], BIN) ** 2)[:, None]
            m0, e0 = jk_mean_err(svac)
            ax.errorbar(tt, m0, e0, fmt="^", color="tab:orange", ms=5, capsize=2, label=f"D vac-sub (N={len(common)})")
        m2, e2 = jk_mean_err(splat)
        ax.errorbar(tt + 0.3, m2, e2, fmt="o", color="tab:green", ms=4, capsize=2,
                    label=f"D plateau-sub, c[{tlo},{thi}]={cm:.2e}({ce:.0e})")
        ax.axvspan(tlo - 0.3, thi + 0.3, color="gray", alpha=0.12)
        ax.axhline(0.0, color="gray", lw=0.8)
        ax.set_xlabel("t")
        ax.set_ylabel("D(t)")
        ax.set_title(rf"{g} disc, linear: " + tail)
        ax.legend(fontsize=8)
        fig.savefig(os.path.join(edir, f"{tag}_p{psq}_{key}_claude.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  {tag} {g}: plateau c[{tlo},{thi}] = {cm:.3e} +- {ce:.1e}")

    # C and C - D, log scale, D = the plateau-subtracted curve of plots 10/11,
    # subtracted from C per jackknife sample (matched configs). Overall sign
    # sign(C(0)) applied to both (I_I has C < 0).
    for key, ch, g, lab in [("12_G5_C_CmD_log", "G5_G5", "g5", r"$\eta'$"), ("13_II_C_CmD_log", "I_I", "id", r"$\sigma$")]:
        win = PLAT_WIN.get((tag, g))
        if win is None or psq != 0:
            continue
        tlo, thi = win
        sC = jk_samples(Cof(ch)[iC], BIN)
        sD = jk_samples(Dof(g)[iD], BIN)
        sD = sD - sD[:, tlo:thi + 1].mean(axis=1)[:, None]
        mC, eC = jk_mean_err(sC)
        mX, eX = jk_mean_err(sC - sD)
        sgn = np.sign(mC[0])
        tt = np.arange(nt2 + 1)
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.set_yscale("log", nonpositive="clip")
        ax.errorbar(tt, sgn * mC, eC, fmt="s", color="tab:blue", ms=5, capsize=2, label=f"C (N={len(common)})")
        ax.errorbar(tt + 0.2, sgn * mX, eX, fmt="o", color="tab:red", ms=5, capsize=2, label="C - D")
        # negative values (after the sign(C(0)) factor): |.| with open markers
        negC = sgn * mC < 0
        negX = sgn * mX < 0
        if np.any(negC):
            ax.errorbar(tt[negC], -sgn * mC[negC], eC[negC], fmt="s", mfc="none", color="tab:blue", ms=5,
                        capsize=2, label="|C| (C < 0)")
        if np.any(negX):
            ax.errorbar(tt[negX] + 0.2, -sgn * mX[negX], eX[negX], fmt="o", mfc="none", color="tab:red", ms=5,
                        capsize=2, label="|C - D| (C - D < 0)")
        ax.set_xlabel("t")
        ax.set_ylabel("correlator")
        ax.set_title(f"{lab}: " + tail)
        ax.legend()
        fig.savefig(os.path.join(edir, f"{tag}_p{psq}_{key}_claude.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)

        # log m_eff = ln[G(t)/G(t+1)] of the same C and C - D (plots 14 / 15)
        mkey = key.replace("12_", "14_").replace("13_", "15_").replace("_log", "_meff")
        fig, ax = plt.subplots(figsize=(7, 4.5))
        t, m, e = meff_jk(sC)
        # t, m, e = meff_log_jk(sC)
        ax.errorbar(t, m, e, fmt="s", color="tab:blue", ms=5, capsize=2, label=f"C cosh (N={len(common)})")
        fwin = CONN_FIT_WIN.get((tag, ch))
        if fwin is not None:
            flo, fhi = fwin
            fm, fe = const_fit_meff_jk(sC, flo, fhi)
            ax.axvspan(flo - 0.3, fhi + 0.3, color="tab:blue", alpha=0.10)
            ax.fill_between([flo - 0.3, fhi + 0.3], fm - fe, fm + fe, color="tab:blue", alpha=0.35, lw=0)
            ax.hlines(fm, flo - 0.3, fhi + 0.3, color="tab:blue", lw=1.5)
            ax.text(0.5 * (flo + fhi), fm + 0.08, f"{fm:.4f}({fe * 1e4:.0f}) [{flo},{fhi}]",
                    ha="center", va="bottom", color="tab:blue", fontsize=9)
            print(f"  {tag} {ch}: conn plateau [{flo},{fhi}] = {fm:.4f} +- {fe:.4f}")
            FIT_RESULTS[(tag, ch)] = (mass, fm, fe)
        # t, m, e = meff_jk(sC - sD)
        t, m, e = meff_log_jk(sC - sD)
        ax.errorbar(t + 0.2, m, e, fmt="o", color="tab:red", ms=5, capsize=2, label="C - D")
        ax.set_xlim(0, TMAX + 0.5)
        ax.set_ylim(0, YMAX)
        ax.set_xlabel("t")
        ax.set_ylabel(r"$m_\mathrm{eff}$")
        ax.set_title(f"{lab}: " + tail)
        ax.legend()
        fig.savefig(os.path.join(edir, f"{tag}_p{psq}_{mkey}_claude.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)
    print(f"{tag}: png -> {edir}")


def main():
    psq = 0
    out_dir = "meson_v2_png_claude"
    if len(sys.argv) > 1:
        psq = int(sys.argv[1])
    if len(sys.argv) > 2:
        out_dir = sys.argv[2]
    for npz in sorted(glob.glob("meson_v2_out_claude/meson_v2_*_claude.npz")):
        do_ensemble(npz, psq, out_dir)

    # summary: connected plateau masses m_eta' (G5) and m_sigma (I_I) vs quark mass m
    if len(FIT_RESULTS) > 0:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for ch, lab, fmt, col in [("G5_G5", r"$m_{\eta'}$ (G5)", "o", "tab:red"),
                                  ("I_I", r"$m_\sigma$ (I_I)", "s", "tab:blue")]:
            pts = sorted(v for k, v in FIT_RESULTS.items() if k[1] == ch)
            x = np.array([p[0] for p in pts])
            y = np.array([p[1] for p in pts])
            e = np.array([p[2] for p in pts])
            ax.errorbar(x, y, e, fmt=fmt, color=col, ms=6, capsize=3, label=lab)
            for xi, yi, ei in zip(x, y, e):
                print(f"  summary {ch} m={xi}: {yi:.4f} +- {ei:.4f}")
        bx = np.array([r[0] for r in BARYON_2PM])
        ax.errorbar(bx, [r[1] for r in BARYON_2PM], [r[2] for r in BARYON_2PM], fmt="^", color="tab:green",
                    ms=6, capsize=3, label=r"$2^+$ baryon (GEVP, ss)")
        ax.errorbar(bx, [r[3] for r in BARYON_2PM], [r[4] for r in BARYON_2PM], fmt="v", color="tab:purple",
                    ms=6, capsize=3, label=r"$2^-$ baryon (GEVP, ss)")
        ax.set_xlabel("quark mass $m$")
        # ax.set_ylabel("meson mass (connected plateau)")
        # ax.set_title(rf"connected $\eta'$ and $\sigma$ vs $m$, $p^2$ class {psq}")
        ax.set_ylabel("mass (lattice units)")
        ax.set_title(rf"connected $\eta'$, $\sigma$ ($p^2$ class {psq}) and $2^\pm$ baryons vs $m$")
        ax.legend()
        fig.savefig(os.path.join(out_dir, f"summary_p{psq}_meson_mass_vs_m_claude.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)

        # ratios m_eta'/m_2+ and m_2-/m_2+ vs m; errors in quadrature (treated uncorrelated)
        bar = {r[0]: r for r in BARYON_2PM}
        fig, ax = plt.subplots(figsize=(7, 4.5))
        pts = sorted(v for k, v in FIT_RESULTS.items() if k[1] == "G5_G5" and v[0] in bar)
        x = np.array([p[0] for p in pts])
        r1 = np.array([p[1] / bar[p[0]][1] for p in pts])
        e1 = r1 * np.sqrt(np.array([(p[2] / p[1]) ** 2 + (bar[p[0]][2] / bar[p[0]][1]) ** 2 for p in pts]))
        ax.errorbar(x, r1, e1, fmt="o", color="tab:red", ms=6, capsize=3, label=r"$m_{\eta'}/m_{2^+}$")
        bx = np.array([r[0] for r in BARYON_2PM])
        r2 = np.array([r[3] / r[1] for r in BARYON_2PM])
        e2 = r2 * np.sqrt(np.array([(r[4] / r[3]) ** 2 + (r[2] / r[1]) ** 2 for r in BARYON_2PM]))
        ax.errorbar(bx, r2, e2, fmt="v", color="tab:purple", ms=6, capsize=3, label=r"$m_{2^-}/m_{2^+}$")
        sig = {v[0]: v for k, v in FIT_RESULTS.items() if k[1] == "I_I"}
        eta = {v[0]: v for k, v in FIT_RESULTS.items() if k[1] == "G5_G5"}
        x3 = np.array(sorted(m for m in sig if m in eta))
        r3 = np.array([sig[m][1] / eta[m][1] for m in x3])
        e3 = r3 * np.sqrt(np.array([(sig[m][2] / sig[m][1]) ** 2 + (eta[m][2] / eta[m][1]) ** 2 for m in x3]))
        ax.errorbar(x3, r3, e3, fmt="s", color="tab:blue", ms=6, capsize=3, label=r"$m_\sigma/m_{\eta'}$")
        for xi, yi, ei in zip(x3, r3, e3):
            print(f"  ratio sigma/eta' m={xi}: {yi:.4f} +- {ei:.4f}")
        for xi, yi, ei in zip(x, r1, e1):
            print(f"  ratio eta'/2+ m={xi}: {yi:.4f} +- {ei:.4f}")
        for xi, yi, ei in zip(bx, r2, e2):
            print(f"  ratio 2-/2+  m={xi}: {yi:.4f} +- {ei:.4f}")
        ax.set_xlabel("quark mass $m$")
        ax.set_ylabel("mass ratio")
        # ax.set_title(rf"$m_{{\eta'}}/m_{{2^+}}$ (connected $\eta'$, $p^2$ class {psq}) and $m_{{2^-}}/m_{{2^+}}$")
        ax.set_title(rf"$m_\sigma/m_{{\eta'}}$, $m_{{\eta'}}/m_{{2^+}}$ (connected mesons, $p^2$ class {psq}), $m_{{2^-}}/m_{{2^+}}$")
        ax.legend()
        fig.savefig(os.path.join(out_dir, f"summary_p{psq}_ratios_vs_m_claude.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)

        # relative 2-/2+ splitting (m_2- - m_2+)/m_2+ = r2 - 1, same error as m_2-/m_2+
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.errorbar(bx, r2 - 1.0, e2, fmt="v", color="tab:purple", ms=6, capsize=3,
                    label=r"$(m_{2^-} - m_{2^+})/m_{2^+}$")
        for xi, yi, ei in zip(bx, r2 - 1.0, e2):
            print(f"  splitting (2- - 2+)/2+ m={xi}: {yi:.4f} +- {ei:.4f}")
        ax.axhline(0.0, color="gray", lw=0.8)
        ax.set_xlabel("quark mass $m$")
        ax.set_ylabel(r"$(m_{2^-} - m_{2^+})/m_{2^+}$")
        ax.set_title(r"relative $2^-$ - $2^+$ baryon splitting (GEVP, ss)")
        ax.legend()
        fig.savefig(os.path.join(out_dir, f"summary_p{psq}_2pm_splitting_vs_m_claude.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()

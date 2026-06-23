#!/usr/bin/env python3
# meson_CminusD_batch_claude.py
#
# Run the physical singlet C-D across ALL channels and ALL ensembles with
# rel_norm=1 -- disc corr_d regenerated to D_match. Saves one npz per ensemble
# holding every channel and p^2 class, plus jackknife errors. (Nf=1 singlet is
# C - Nf*D; the overall sign is irrelevant for masses, so C-D is used to match
# the file name. See meson_singlet_sign_derivation_claude.md.)
#
# Output: CminusD_out_claude/CminusD_<ensemble>.npz with arrays indexed
# [channel, p2class, t]:  phys_mean/phys_err, D_mean/D_err, C_mean/C_err,
# channels (names), t, common (config indices), rel_norm.

import os
import glob
import re
import numpy as np
import meson_combine_claude as mc

OBS = "/mnt/baracuda_14/grid_claude/obs_nc4nf1_2448"
REL_NORM = 1.0
OUTDIR = "CminusD_out_claude"

chans = list(mc.CHANNEL_TO_TRACEDIR.keys())


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    ensembles = sorted(d for d in glob.glob(OBS + "/obs_nc4nf1_2448_*")
                       if os.path.isdir(d) and glob.glob(d + "/mesons_conn.*.h5"))

    for ens in ensembles:
        name = os.path.basename(ens)

        # connected: one pass, all channels -> Cc[nconf,10,5,nh]
        files = sorted((int(re.search(r'conn\.(\d+)', f).group(1)), f)
                       for f in glob.glob(ens + "/mesons_conn.*.h5"))
        cidx = []
        Cc = []
        nbad = 0
        for idx, f in files:
            try:
                perch = [mc.connected_pclass(f, ch) for ch in chans]
            except (OSError, KeyError):
                nbad += 1
                continue
            Cc.append(perch)
            cidx.append(idx)
        Cc = np.array(Cc)            # (nconf,10,5,nh)
        cidx = np.array(cidx)
        cpos = {d: i for i, d in enumerate(cidx)}
        nh = Cc.shape[-1]

        nchan = len(chans)
        phys_mean = np.zeros((nchan, 5, nh))
        phys_err = np.zeros((nchan, 5, nh))
        D_mean = np.zeros((nchan, 5, nh))
        D_err = np.zeros((nchan, 5, nh))
        C_mean = np.zeros((nchan, 5, nh))
        C_err = np.zeros((nchan, 5, nh))
        ncommon_ch = []

        for ci, ch in enumerate(chans):
            td = mc.CHANNEL_TO_TRACEDIR[ch]
            # disc per p-class, matched indices
            Didx = {}
            Dval = {}
            for m in range(5):
                di, D = mc.read_corr_d_indexed(os.path.join(ens, td, f"corr_d.p{m}"))
                Didx[m] = {d: i for i, d in enumerate(di)}
                Dval[m] = D
            common = [k for k in cidx if all(k in Didx[m] for m in range(5))]
            ncommon_ch.append(len(common))
            for m in range(5):
                Dal = np.array([Dval[m][Didx[m][k]] for k in common])
                Cal = np.array([Cc[cpos[k], ci, m] for k in common])
                phys = Cal - 1.0 * REL_NORM * Dal  # Nf=1 singlet: C - D
                phys_mean[ci, m], phys_err[ci, m] = mc.jackknife(phys)
                D_mean[ci, m], D_err[ci, m] = mc.jackknife(Dal)
                C_mean[ci, m], C_err[ci, m] = mc.jackknife(Cal)

        out = os.path.join(OUTDIR, f"CminusD_{name}.npz")
        np.savez(out, channels=np.array(chans), t=np.arange(nh),
                 rel_norm=REL_NORM, ncommon=np.array(ncommon_ch),
                 phys_mean=phys_mean, phys_err=phys_err,
                 D_mean=D_mean, D_err=D_err, C_mean=C_mean, C_err=C_err)

        skip = f" ({nbad} bad h5 skipped)" if nbad else ""
        print(f"{name}: connected={cidx.size}  matched~{min(ncommon_ch)}-"
              f"{max(ncommon_ch)}{skip}  -> {out}")


if __name__ == "__main__":
    main()

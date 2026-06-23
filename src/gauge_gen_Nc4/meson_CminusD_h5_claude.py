#!/usr/bin/env python3
# meson_CminusD_h5_claude.py
#
# Store the per-config connected, disconnected, and singlet meson correlators
# in HDF5 (Nf=1; overall sign irrelevant, C-D chosen to match the name).
#
# One file per ensemble: corr_h5_out_claude/corr_<ens>.h5
#   datasets : "<contr>/<ch>/<mom>"  ->  real array shape (Nconf_ch, NT/2+1)
#              one row per config; rows aligned with "confs/<ch>".
#              contr = contribution:
#                        C          = connected meson
#                        Dsub       = sum-subtracted disconnected loop-loop
#                        CminusDsub = C - Nf*Dsub (physical singlet)
#              ch    = channel name (G5_G5, GTG5_GTG5, ..., 10 channels)
#              mom   = p^2 class label p0..p4 (disc exists only per p^2 class)
#   "confs/<ch>" : int array of config numbers for that channel (row order)
#   root keys: beta, m  (scalar datasets = the ensemble parameters)
#   attrs    : ensemble, Nf, Nc, note
#
# This (Nconf, NT/2+1)-array layout (vs one tiny dataset per config) keeps the
# file small: ~150 datasets instead of ~Nconf*150, avoiding HDF5 per-object
# overhead.
#
# Only configs present in BOTH the connected h5 and the disc corr_d are written.
# C : connected meson, mesons_conn.<conf>.h5 (folded + class-averaged here).
# D : disconnected loop-loop, Rebbi's corr_d.p{0..4} (per config, folded), with
#     the per-config DC/zero-mode constant SUM-SUBTRACTED:
#       D_sub(t) = D(t) - (1/NT) sum_t' D(t'),  the average reconstructed from
#       the folded 0..NT/2 with weights (1,2,..,2,1). Removes the non-decaying
#       constant (the g5 loop has no true VEV). Same subtraction as the notebook.

import os
import glob
import re
import numpy as np
import h5py
import meson_combine_claude as mc

OBS = "/mnt/baracuda_14/grid_claude/obs_nc4nf1_2448"
OUTDIR = "corr_h5_out_claude"
NF = 1
NC = 4
NPCLASS = 5

chans = list(mc.CHANNEL_TO_TRACEDIR.keys())


def parse_beta_mass(ens_name):
    mobj = re.match(r'obs_nc4nf1_2448_b(\d+p\d+)_m(\d+p\d+)$', ens_name)
    beta = float(mobj.group(1).replace('p', '.'))
    mass = float(mobj.group(2).replace('p', '.'))
    return beta, mass


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    ensembles = []
    for d in sorted(glob.glob(OBS + "/obs_nc4nf1_2448_b*_m*")):
        if re.match(r'obs_nc4nf1_2448_b\d+p\d+_m\d+p\d+$', os.path.basename(d)):
            ensembles.append(d)

    for ens in ensembles:
        name = os.path.basename(ens)
        beta, mass = parse_beta_mass(name)

        # config number -> connected h5 path
        conn_files = {}
        for f in glob.glob(ens + "/mesons_conn.*.h5"):
            conf = int(re.search(r'conn\.(\d+)', f).group(1))
            conn_files[conf] = f

        # disc[ch][mclass][conf] = D(t) array; skip channels with no corr_d
        disc = {}
        for ch in chans:
            td = mc.CHANNEL_TO_TRACEDIR[ch]
            dd = {}
            ok = True
            for mclass in range(NPCLASS):
                dfile = os.path.join(ens, td, f"corr_d.p{mclass}")
                if not os.path.exists(dfile):
                    ok = False
                    break
                didx, D = mc.read_corr_d_indexed(dfile)
                per_conf = {}
                for i, dval in enumerate(didx):
                    per_conf[int(dval)] = D[i]
                dd[mclass] = per_conf
            if ok:
                disc[ch] = dd

        # configs we actually need connected for = union of disc configs
        needed = set()
        for ch in disc:
            for mclass in range(NPCLASS):
                needed |= set(disc[ch][mclass].keys())
        needed &= set(conn_files.keys())

        # read connected only for the needed configs (all channels per file)
        conn = {}
        for conf in sorted(needed):
            f = conn_files[conf]
            per_ch = {}
            for ch in chans:
                try:
                    per_ch[ch] = mc.connected_pclass(f, ch)   # (NPCLASS, NT/2+1)
                except (OSError, KeyError):
                    continue
            conn[conf] = per_ch

        out = os.path.join(OUTDIR, f"corr_{name}.h5")
        n_written = 0
        with h5py.File(out, "w") as h5:
            h5["beta"] = beta
            h5["m"] = mass
            h5.attrs["ensemble"] = name
            h5.attrs["Nf"] = NF
            h5.attrs["Nc"] = NC
            h5.attrs["note"] = "datasets <contr>/<ch>/<mom> shape (Nconf_ch, NT/2+1) (one row per config, aligned with confs/<ch>); contr in {C, Dsub, CminusDsub}; CminusDsub = C - Nf*Dsub; Dsub = disc with per-config DC/zero-mode subtracted; mom=p^2 class"

            nconf_ch = []
            for ch in chans:
                if ch not in disc:
                    continue
                # matched configs for this channel (present in conn and all p-classes)
                ch_confs = sorted(c for c in conn
                                  if ch in conn[c]
                                  and all(c in disc[ch][m] for m in range(NPCLASS)))
                if not ch_confs:
                    continue
                h5[f"confs/{ch}"] = np.array(ch_confs, dtype=int)
                nconf_ch.append(len(ch_confs))

                for mclass in range(NPCLASS):
                    C_rows = []
                    Dsub_rows = []
                    CmD_rows = []
                    for conf in ch_confs:
                        Dv = disc[ch][mclass][conf]
                        nt2 = Dv.shape[0] - 1
                        avg = (Dv[0] + Dv[nt2] + 2.0 * Dv[1:nt2].sum()) / (2 * nt2)
                        Dsub = Dv - avg
                        C = conn[conf][ch][mclass]
                        C_rows.append(C)
                        Dsub_rows.append(Dsub)
                        CmD_rows.append(C - NF * Dsub)
                    h5[f"C/{ch}/p{mclass}"] = np.array(C_rows)
                    h5[f"Dsub/{ch}/p{mclass}"] = np.array(Dsub_rows)
                    h5[f"CminusDsub/{ch}/p{mclass}"] = np.array(CmD_rows)
                    n_written += 1

        lo = min(nconf_ch) if nconf_ch else 0
        hi = max(nconf_ch) if nconf_ch else 0
        print(f"{name}: beta={beta} m={mass}  Nconf_ch={lo}-{hi}  "
              f"wrote {n_written} (contr,ch,mom) arrays -> {out}")


if __name__ == "__main__":
    main()

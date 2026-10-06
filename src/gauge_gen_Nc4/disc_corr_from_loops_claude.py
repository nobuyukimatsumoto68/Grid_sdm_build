#!/usr/bin/env python3
# disc_corr_from_loops_claude.py
#
# Build the disconnected loop-loop correlator D(t, p^2 class) from the complex
# loop cache written by disc_loop_cache_claude.py, with the SAME normalization
# as Rebbi's makecorr_96_d_claude.f90 (factor-3 free):
#   pcorr_j(t) = (NS^3/NT) sum_i Re[ X_j(i) conj(X_j(i+t)) ]
#   D(t, m)    = pcorr_0(t)                         (m = 0)
#              = mean_{j in class m} pcorr_j(t)     (m > 0; = sum/(pcomp/2))
# X = projection of the real-part field r (R) or of the imaginary-part field
# l (I) of tr[Gamma S] = r + i l.
#
# gamma5-hermiticity: the loop is real for Gamma in {1, g5, g_mu g5} and
# imaginary for Gamma = g_mu. The physical D is the product of the two loops:
#   real loops:      D = +<r r>
#   imaginary loops: D = <(i l)(i l)> = -<l l>
# The other part has zero expectation (stochastic noise only).
# The old Fortran pipeline (average_trace2) always used r -> pure noise for
# the vector channels gx, gy, gz, gt.

import numpy as np

NS = 24
NT = 48
GAMMAS = ["id", "g5", "gx", "gy", "gz", "gt", "gxg5", "gyg5", "gzg5", "gtg5"]

# which part of tr[Gamma S] = r + i l carries the signal
SIGNAL_PART = {"id": "R", "g5": "R", "gx": "I", "gy": "I", "gz": "I", "gt": "I",
               "gxg5": "R", "gyg5": "R", "gzg5": "R", "gtg5": "R"}

# momentum index (0..16, average_trace2 ptable) -> p^2 class
PCLASS = [[0], [1, 2, 3], [4, 5, 6, 7, 8, 9], [10, 11, 12, 13], [14, 15, 16]]

# h5 channel name -> loop gamma
CHANNEL_TO_GAMMA = {"I_I": "id", "G5_G5": "g5", "GX_GX": "gx", "GY_GY": "gy",
                    "GZ_GZ": "gz", "GT_GT": "gt", "GXG5_GXG5": "gxg5",
                    "GYG5_GYG5": "gyg5", "GZG5_GZG5": "gzg5", "GTG5_GTG5": "gtg5"}


def load_loops(fn):
    # Returns confs, gammas, R, I. Old-format cache (key L, full complex
    # projection) only supports p0: R[...,0] = Re L, I[...,0] = Im L; the
    # p>0 entries are set to nan.
    z = np.load(fn)
    confs = z["confs"]
    gammas = list(z["gammas"])
    if "R" in z.files:
        return confs, gammas, z["R"], z["I"]
    L = z["L"]
    R = np.full(L.shape, np.nan, dtype=complex)
    I = np.full(L.shape, np.nan, dtype=complex)
    R[..., 0] = L[..., 0].real
    I[..., 0] = L[..., 0].imag
    return confs, gammas, R, I


def pcorr(X):
    # X: (..., NT) complex projection for one momentum. Returns (..., NT/2+1).
    out = np.zeros(X.shape[:-1] + (NT // 2 + 1,))
    for t in range(NT // 2 + 1):
        out[..., t] = np.sum((X * np.conj(np.roll(X, -t, axis=-1))).real, axis=-1)
    return out * NS**3 / NT


def disc_corr(X):
    # X: (Nconf, NT, 17) projections. Returns D (Nconf, 5, NT/2+1).
    D = np.zeros((X.shape[0], len(PCLASS), NT // 2 + 1))
    for m, idx in enumerate(PCLASS):
        acc = np.zeros((X.shape[0], NT // 2 + 1))
        for j in idx:
            acc = acc + pcorr(X[:, :, j])
        D[:, m] = acc / len(idx)
    return D


def disc_corr_gamma(R, I, gammas, g, part=None):
    # Physical D for loop gamma g, from its signal part (or a forced part).
    ig = gammas.index(g)
    if part is None:
        part = SIGNAL_PART[g]
    if part == "R":
        return disc_corr(R[:, ig])
    # loop = i l -> product (i l)(i l) = -l l
    return -disc_corr(I[:, ig])

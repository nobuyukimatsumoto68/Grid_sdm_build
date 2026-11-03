#!/bin/bash
# meson_v2_overview_merge_claude.sh
# Merge the p = 0 overview (summary_p0_2expfit_overview_claude.pdf, 12 pages: ens x {sigma, eta'})
# and the p^2 = 1-4 overview (summary_pmom_2expfit_overview_claude.pdf, 48 pages: ens x ch x class 1-4)
# into ONE PDF with, per (ensemble, channel), class 0 then classes 1-4 (60 pages):
#   meson_v2_png_claude/summary_pall_2expfit_overview_claude.pdf
# No fits are rerun; rebuild the two inputs first if ranges changed. Page splits go to a scratch dir
# (argument 1, default ./merge_tmp_claude), never deleted by this script.
set -u

cd "$(dirname "$0")" || exit 1
P0=meson_v2_png_claude/summary_p0_2expfit_overview_claude.pdf
PM=meson_v2_png_claude/summary_pmom_2expfit_overview_claude.pdf
OUT=meson_v2_png_claude/summary_pall_2expfit_overview_claude.pdf
TMP=${1:-./merge_tmp_claude}

mkdir -p "${TMP}"
pdfseparate "${P0}" "${TMP}/p0_%02d.pdf" || exit 1
pdfseparate "${PM}" "${TMP}/pm_%02d.pdf" || exit 1

pages=""
for ens in 0 1 2 3 4 5; do
    for ch in 0 1; do
        i0=$(( ens * 2 + ch + 1 ))
        pages="${pages} ${TMP}/p0_$(printf %02d ${i0}).pdf"
        for k in 1 2 3 4; do
            im=$(( ens * 8 + ch * 4 + k ))
            pages="${pages} ${TMP}/pm_$(printf %02d ${im}).pdf"
        done
    done
done
pdfunite ${pages} "${OUT}" || exit 1
echo "pdf -> ${OUT}"

#!/usr/bin/env python3
"""
23_figure_calibration.py
========================
Figure 7: behaviour of the mixed model under phenotypes that carry no
association by construction.

Panel (a) is the empirical distribution of the smallest p-value returned per
permutation, against the two thresholds that matter: the Bonferroni threshold
actually applied in the analysis, and the smallest value the exact test can
return for the design. Reading off the curve at the Bonferroni line gives the
family-wise error rate the analysis would have had.

Panel (b) contrasts that observed error rate with the nominal 5%, and the
floor-crossing rate with the union bound that any calibrated test must satisfy.

The script also prints the numbers needed for the text.

Usage:
  python3 23_figure_calibration.py \\
      --calib   results/calibration_v2 \\
      --n       29 \\
      --m       636 \\
      --outdir  figures
"""

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "figure.dpi": 150, "savefig.dpi": 900, "savefig.bbox": "tight",
})

C = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
     "red": "#D55E00", "grey": "#999999", "dark": "#333333"}
MM = 1 / 25.4
DOUBLE = 180 * MM

# smallest positive value plotted; the model returns exact zeros, which cannot
# be shown on a log axis and are in any case not p-values
FLOOR_PLOT = 1e-300


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calib", required=True)
    ap.add_argument("--n", type=int, default=29)
    ap.add_argument("--m", type=int, default=636,
                    help="number of distinct presence/absence patterns tested")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--outdir", default="figures")
    a = ap.parse_args()

    bonf = a.alpha / a.m
    rows, data = [], {}

    for r in (4, 6, 8):
        f = Path(a.calib) / f"perm_r{r}.tsv"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t")
        floor = 1.0 / math.comb(a.n, r)
        mp = d["min_p_lmm"].to_numpy(dtype=float)

        # Ties with the floor are not crossings. The permutation tables store
        # p-values with 7 significant digits, so a Fisher p exactly equal to
        # the floor can be written marginally below it; the tolerance absorbs
        # that rounding.
        tol = 1e-6
        cross = float(np.mean(mp < floor * (1 - tol)))
        fisher = d["min_p_fisher"].to_numpy(dtype=float)
        fisher_cross = (float(np.nanmean(fisher < floor * (1 - tol)))
                        if np.isfinite(fisher).any() else float("nan"))

        data[r] = (mp, floor)
        rows.append({
            "r": r,
            "permutations": len(d),
            "empirical_FWER_pct": round(100 * float(np.mean(mp < bonf)), 2),
            "nominal_FWER_pct": 100 * a.alpha,
            "floor": floor,
            "cross_floor_pct": round(100 * cross, 2),
            "union_bound_pct": round(100 * min(1.0, a.m * floor), 3),
            "excess_fold": round(cross / min(1.0, a.m * floor), 1),
            "fisher_cross_pct": round(100 * fisher_cross, 3),
            "median_lambda": round(d["lambda"].median(), 3),
            "p90_lambda": round(d["lambda"].quantile(0.90), 3),
            "median_min_p": f"{np.median(mp):.2e}",
            "smallest_min_p": f"{np.min(mp):.2e}",
            "n_exact_zero": int((mp == 0).sum()),
        })

    s = pd.DataFrame(rows)
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    s.to_csv(outdir / "calibration_for_text.tsv", sep="\t", index=False)
    print(s.to_string(index=False))
    print(f"\n  Bonferroni threshold applied in the analysis: {bonf:.2e}")
    print(f"  wrote {outdir/'calibration_for_text.tsv'}\n")

    # ── figure ──────────────────────────────────────────────────────────────
    # Fig. 7 of the manuscript. The x axis starts at 1e-15: values below it
    # (including exact zeros returned by the model) cannot be shown on a log
    # axis and are reported as a proportion in an annotation instead.
    XMIN = 1e-15
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(DOUBLE, 68 * MM),
                                   gridspec_kw={"width_ratios": [1.35, 1]})

    col = {4: C["red"], 6: C["green"], 8: C["blue"]}
    for r, (mp, floor) in sorted(data.items()):
        x = np.sort(np.clip(mp, XMIN, None))
        y = np.arange(1, len(x) + 1) / len(x)
        axA.step(x, y, where="post", lw=1.3, color=col[r],
                 label=f"$r$ = {r} resistant of {a.n}")
        axA.plot([floor], [np.mean(mp < floor)], "o", ms=4,
                 mfc="white", mec=col[r], mew=1.1, zorder=5)

    axA.axvline(bonf, color=C["dark"], lw=1.0, ls="--")
    axA.axhline(a.alpha, color=C["grey"], lw=0.8, ls=":")
    # legend-only proxies for the reference marks
    axA.plot([], [], "o", ms=4, mfc="white", mec=C["dark"], mew=1.1,
             label=f"exact-test floor 1/C({a.n}, $r$)")
    axA.plot([], [], color=C["dark"], lw=1.0, ls="--",
             label=f"Bonferroni threshold ({a.alpha:g}/{a.m})")
    axA.plot([], [], color=C["grey"], lw=0.8, ls=":", label="nominal 5%")

    if 4 in data:
        mp4 = data[4][0]
        frac = 100 * float(np.mean(mp4 < XMIN))
        axA.text(2e-15, 0.30,
                 f"{frac:.1f}% of $r$ = 4 permutations returned\n"
                 "values below 10$^{-15}$, including exact\n"
                 "zeros by numerical underflow",
                 fontsize=6, color="#555555", va="center")

    axA.set_xscale("log")
    axA.set_xlim(XMIN, 1)
    axA.set_ylim(-0.01, 1.02)
    axA.set_xlabel("Smallest $p$-value returned in a permutation")
    axA.set_ylabel("Proportion of permutations at or below")
    axA.set_title("a", loc="left", fontweight="bold")
    axA.legend(frameon=True, framealpha=1, edgecolor="#cccccc",
               loc="upper left", fontsize=6.3)

    # (b) observed against nominal
    lab = [f"$r$ = {r}" for r in s["r"]]
    x = np.arange(len(lab)); w = 0.36
    axB.bar(x - w / 2, s["empirical_FWER_pct"], w, color=C["red"],
            label="observed")
    axB.bar(x + w / 2, s["nominal_FWER_pct"], w, color="#bbbbbb",
            label="nominal (5%)")
    for xi, v in zip(x - w / 2, s["empirical_FWER_pct"]):
        axB.text(xi, v + 1.5, f"{v:.1f}", ha="center", fontsize=7,
                 color=C["red"], fontweight="bold")
    axB.set_xticks(x); axB.set_xticklabels(lab)
    axB.set_ylabel("Permutations yielding a\nsignificant hit (%)")
    axB.set_ylim(0, max(60, s["empirical_FWER_pct"].max() * 1.15))
    axB.set_title("b", loc="left", fontweight="bold")
    axB.legend(frameon=True, framealpha=1, edgecolor="#cccccc",
               loc="upper right")

    fig.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(outdir / f"Fig7_calibration.{ext}")
    plt.close(fig)
    print(f"  wrote {outdir/'Fig7_calibration.pdf'} / .png / .svg")

    print("""
  FOR THE TEXT
    empirical_FWER_pct is the proportion of permutations, none of which
    contains any association, in which the analysis as performed would have
    reported at least one significant result. It should equal the nominal 5%.

    cross_floor_pct and union_bound_pct are the corresponding quantities for
    the exact-test floor. excess_fold is their ratio.

    n_exact_zero counts permutations in which the model returned a p-value of
    exactly zero. With 29 observations no test statistic can support such a
    value; these are not small p-values but a failure of the approximation.
""")


if __name__ == "__main__":
    main()

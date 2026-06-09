#!/usr/bin/env python3
"""
10_figure_phylo_correction.py
------------------------------
Figure showing the effect of each methodological step on SNP-MIC associations.
Three bars per antibiotic:
  1. Fisher only, MIN_STRAINS=2 (original)
  2. Fisher only, MIN_STRAINS=5 (stricter)
  3. After phylogenetic correction (CP > DP)

Usage:
    python3 10_figure_phylo_correction.py

Inputs:
    ../data/processed/MICs_converted.tsv
    ../data/raw/allnsSNP.tsv
    ../results/gene_level_phylo_corrected.tsv

Outputs:
    ../figures/Fig5_phylo_correction.pdf/png

Author: Francesc Montoro (2025)
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

BASE    = Path(__file__).parent.parent
RAW     = BASE / "data" / "raw"
PROC    = BASE / "data" / "processed"
RES     = BASE / "results"
FIGURES = BASE / "figures"
FIGURES.mkdir(exist_ok=True)

REFERENCE = "CECT_5769"
MIN_DELTA = 0.5

SNP_TO_MIC = {
    "11028": "BT3_11028", "12": "BT3_12", "162": "BT3_162", "32": "BT3_32",
    "94385": "94-3-85",   "960426-1_4C": "960426_1/4C",
    "CIP8190": "C1P8190", "V4": "V_4",
    "CECT4604": "CECT_4604", "CECT4606": "CECT_4606", "CECT4608": "CECT_4608",
    "CECT4865": "CECT_4865", "CECT4866": "CECT_4866", "CECT4999": "CECT_4999",
    "CECT5198": "CECT_5198", "CECT5763": "CECT_5763", "CECT5769": "CECT_5769",
    "CECT7030": "CECT_7030", "CECT898":  "CECT_898",  "CG100":    "CG_100",
    "94-9-112": "94-9-112", "94-9-119": "94-9-119",   "95-8-161": "95-8-161",
    "95-8-6":   "95-8-6",   "95-8-7":   "95-8-7",     "A14":      "A14",
    "PD-2-51":  "PD-2-51",  "Ra_3":     "Ra_3",        "RIU1":     "RIU1",
    "YJ016":    "YJ016",
}

AB_NAMES = {"GEN2": "Gentamicin (low range)", "STZ": "Sulphathiazole"}


def bh_correction(pvalues):
    n = len(pvalues)
    if n == 0:
        return pvalues.copy()
    idx = pvalues.argsort()
    adj = np.minimum(1.0, pvalues.iloc[idx].values * n / np.arange(1, n+1))
    for i in range(n-2, -1, -1):
        adj[i] = min(adj[i], adj[i+1])
    result = pd.Series(index=pvalues.index, dtype=float)
    result.iloc[idx] = adj
    return result


def count_significant_genes(snp_raw, mics, ref, genome_strains,
                             antibiotics, min_strains):
    """Run Fisher test and count significant genes per antibiotic."""
    snp_raw = snp_raw.copy()
    snp_raw["snp_id"] = (
        snp_raw["Gene"].str.strip() + "__" +
        snp_raw["Position"].astype(str) + "__" +
        snp_raw["Reference"].str.strip() + "__" +
        snp_raw["SNP"].str.strip()
    )

    mat = (snp_raw.groupby(["snp_id", "strain_mic"])
           .size().unstack(fill_value=0).clip(upper=1))
    for s in genome_strains:
        if s not in mat.columns:
            mat[s] = 0
    mat = mat[genome_strains]

    gene_counts = {}
    for ab in antibiotics:
        if ab not in mics.columns:
            continue
        resistant = set(s for s in genome_strains
                        if s in mics.index
                        and mics.loc[s, ab] >= ref[ab] + MIN_DELTA)
        sensitive = set(genome_strains) - resistant
        n_res = len(resistant)
        if n_res < min_strains:
            continue

        records = []
        for snp_id, row in mat.iterrows():
            with_snp    = set(row[row == 1].index)
            without_snp = set(row[row == 0].index)
            if len(with_snp) < min_strains:
                continue
            a = len(with_snp & resistant)
            b = len(with_snp & sensitive)
            c = len(without_snp & resistant)
            d = len(without_snp & sensitive)
            if a + b == 0 or c + d == 0:
                continue
            _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            gene = snp_id.split("__")[0]
            records.append({"gene": gene, "p": p})

        if not records:
            continue
        df = pd.DataFrame(records)
        df["fdr"] = bh_correction(df["p"])
        gene_counts[ab] = df[df["fdr"] < 0.05]["gene"].nunique()

    return gene_counts


if __name__ == "__main__":
    print("=" * 60)
    print("  Figure 5 — Methodological steps comparison")
    print("=" * 60)

    # Load data
    mics = pd.read_csv(PROC / "MICs_converted.tsv", sep="\t", index_col=0)
    mics.index = mics.index.astype(str).str.strip()
    ref = mics.loc[REFERENCE]

    snp_raw = pd.read_csv(RAW / "allnsSNP.tsv", sep="\t")
    snp_raw.columns = snp_raw.columns.str.strip()
    snp_raw["Strain"] = snp_raw["Strain"].astype(str).str.strip()
    snp_raw["strain_mic"] = snp_raw["Strain"].map(SNP_TO_MIC)
    snp_raw = snp_raw[snp_raw["strain_mic"].notna()].copy()
    genome_strains = sorted(snp_raw["strain_mic"].unique())

    antibiotics = ["GEN2", "STZ"]

    # Step 1: Fisher MIN_STRAINS=2
    print("Step 1: Fisher MIN_STRAINS=2...")
    n_fisher_2 = count_significant_genes(
        snp_raw, mics, ref, genome_strains, antibiotics, min_strains=2)
    print(f"  {n_fisher_2}")

    # Step 2: Fisher MIN_STRAINS=5
    print("Step 2: Fisher MIN_STRAINS=5...")
    n_fisher_5 = count_significant_genes(
        snp_raw, mics, ref, genome_strains, antibiotics, min_strains=5)
    print(f"  {n_fisher_5}")

    # Step 3: After phylogenetic correction
    after = pd.read_csv(RES / "gene_level_phylo_corrected.tsv", sep="\t")
    n_phylo = {ab: after[after["antibiotic"] == ab]["gene"].nunique()
               for ab in antibiotics}
    print(f"Step 3: After phylo correction: {n_phylo}")

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5), facecolor="white")

    colors = ["#cccccc", "#aaaaaa", "#4393c3"]  # GEN2 colors
    colors_stz = ["#cccccc", "#aaaaaa", "#9970ab"]  # STZ colors

    labels = ["Fisher\nMIN=2\n(original)",
              "Fisher\nMIN=5\n(stricter)",
              "Pairwise\ncomparisons\n(CP > DP)"]

    step_labels = ["Step 1:\nFisher\n(MIN=2)",
                   "Step 2:\nFisher\n(MIN=5)",
                   "Step 3:\nPhylo\ncorrection"]

    for ax, ab in zip(axes, antibiotics):
        vals = [
            n_fisher_2.get(ab, 0),
            n_fisher_5.get(ab, 0),
            n_phylo.get(ab, 0),
        ]
        col = ["#cccccc", "#aaaaaa",
               "#4393c3" if ab == "GEN2" else "#9970ab"]

        bars = ax.bar(step_labels, vals, color=col,
                      width=0.55, edgecolor="white", linewidth=1.5)

        # Value labels
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + max(vals) * 0.02,
                    f"{val:,}", ha="center", va="bottom",
                    fontsize=13, fontweight="bold")

        # Reduction annotations
        for i in range(len(vals)-1):
            if vals[i] > 0:
                pct = (vals[i] - vals[i+1]) / vals[i] * 100
                sign = "−" if pct >= 0 else "+"
                ax.annotate("",
                            xy=(i+1, vals[i+1] + max(vals)*0.04),
                            xytext=(i, vals[i] + max(vals)*0.04),
                            arrowprops=dict(arrowstyle="->",
                                            color="#555555", lw=1.5))
                ax.text(i + 0.5, max(vals) * 0.75,
                        f"{sign}{abs(pct):.0f}%",
                        ha="center", fontsize=10, color="#333333",
                        bbox=dict(boxstyle="round,pad=0.25",
                                  fc="white", ec="#cccccc", alpha=0.9))

        ax.set_title(AB_NAMES.get(ab, ab), fontsize=12,
                     fontweight="bold", pad=12)
        ax.set_ylabel("Significant genes", fontsize=10)
        ax.set_ylim(0, max(vals) * 1.3)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=9)
        ax.grid(axis="y", color="#eeeeee", lw=0.8, zorder=0)

    fig.suptitle(
        r"Effect of methodological filters on SNP–MIC associations"
        "\n" r"($\it{Vibrio\ vulnificus}$, n=29 strains with WGS)",
        fontsize=12, fontweight="bold", y=1.01
    )

    legend_handles = [
        mpatches.Patch(color="#cccccc", label="Step 1: Fisher (MIN_STRAINS=2, uncorrected)"),
        mpatches.Patch(color="#aaaaaa", label="Step 2: Fisher (MIN_STRAINS=5, stricter)"),
        mpatches.Patch(color="#4393c3", label="Step 3: Pairwise comparisons — GEN2 (CP > DP)"),
        mpatches.Patch(color="#9970ab", label="Step 3: Pairwise comparisons — STZ (CP > DP)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=2, fontsize=8, bbox_to_anchor=(0.5, -0.10),
               framealpha=0.9)

    fig.tight_layout()

    for ext in ["pdf", "png"]:
        out = FIGURES / f"Fig5_phylo_correction.{ext}"
        fig.savefig(out, dpi=900, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")

    plt.close(fig)
    print("Done.")

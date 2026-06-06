#!/usr/bin/env python3
"""
05_figure_volcano.py
--------------------
Volcano plot of SNP-antibiotic associations, collapsed to gene level.
One point per (gene, antibiotic/family) pair.
X-axis: log2(best odds ratio)
Y-axis: -log10(best FDR)
Points coloured by antibiotic class.
Top significant genes labelled with adjustText.

Usage:
    python3 05_figure_volcano.py

Inputs:
    ../results/gene_level_per_ab.tsv
    ../results/gene_level_per_family.tsv

Outputs:
    ../figures/Fig2a_volcano_per_ab.pdf/png
    ../figures/Fig2b_volcano_per_family.pdf/png

Author: Francesc Montoro (2025)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False
    print("  [INFO] adjustText not found — labels may overlap. "
          "Install with: pip install adjustText --break-system-packages")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = Path(__file__).parent.parent
RES     = BASE / "results"
FIGURES = BASE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
FDR_THRESHOLD = 0.05
MAX_OR        = 999     # placeholder for inf in gene_level files
TOP_LABELS    = 12      # max gene labels

# ── Colours ───────────────────────────────────────────────────────────────────
AB_CLASS = {
    "AMI":"Aminoglycoside",  "GEN":"Aminoglycoside",  "GEN2":"Aminoglycoside",
    "NEO":"Aminoglycoside",  "TOB":"Aminoglycoside",  "STR":"Aminoglycoside",
    "SPE":"Aminoglycoside",
    "AMP":"Beta-lactam",     "FAZ":"Beta-lactam",     "FEP":"Beta-lactam",
    "FOX":"Beta-lactam",     "AXO":"Beta-lactam",     "TAZ":"Beta-lactam",
    "POD":"Beta-lactam",     "PIP":"Beta-lactam",     "AMOX":"Beta-lactam",
    "PEN":"Beta-lactam",     "TIO":"Beta-lactam",     "AZT":"Beta-lactam",
    "A/S2":"Beta-lactam",    "P/T4":"Beta-lactam",    "TIM2":"Beta-lactam",
    "IMI":"Carbapenem",      "MERO":"Carbapenem",
    "CIP":"Fluoroquinolone", "ENRO":"Fluoroquinolone","GAT":"Fluoroquinolone",
    "NOV":"Fluoroquinolone",
    "ERY":"Macrolide",       "TYLT":"Macrolide",
    "CLI":"Lincosamide",
    "TET":"Tetracycline",    "OXY":"Tetracycline",
    "SDM":"Sulfonamide",     "STZ":"Sulfonamide",     "SXT":"Sulfonamide",
    "SXT2":"Sulfonamide",    "TANS":"Sulfonamide",
    "FFN":"Phenicol",
    "FUR":"Other",           "NIT":"Other",
}

CLASS_COLORS = {
    "Aminoglycoside":  "#4393c3",
    "Beta-lactam":     "#4dac26",
    "Beta-lactam combo":"#1a7a1a",
    "Carbapenem":      "#006d2c",
    "Fluoroquinolone": "#e08214",
    "Macrolide":       "#8c510a",
    "Lincosamide":     "#543005",
    "Tetracycline":    "#d6604d",
    "Sulfonamide":     "#9970ab",
    "Phenicol":        "#c994c7",
    "Other":           "#888888",
}

FAMILY_COLORS = {
    "Aminoglycoside":  "#4393c3",
    "Beta-lactam":     "#4dac26",
    "Carbapenem":      "#006d2c",
    "Fluoroquinolone": "#e08214",
    "Macrolide":       "#8c510a",
    "Lincosamide":     "#543005",
    "Tetracycline":    "#d6604d",
    "Sulfonamide":     "#9970ab",
    "Phenicol":        "#c994c7",
    "Other":           "#888888",
}


def prep_gene_level(df: pd.DataFrame,
                    group_col: str,
                    color_map: dict) -> pd.DataFrame:
    """
    Prepare gene-level data for volcano plot.
    best_OR column uses 999 as placeholder for inf — restore and cap for display.
    """
    df = df.copy()

    # best_OR: 999 = inf in our pipeline; cap at 100 for display
    df["or_plot"] = df["best_OR"].clip(upper=100)
    df["log2_or"] = np.log2(df["or_plot"].clip(lower=1e-6))

    # FDR: -log10, floor at 1e-10
    df["neg_log_fdr"] = -np.log10(df["best_fdr"].clip(lower=1e-10))

    # Significance flag
    df["significant"] = df["best_fdr"] < FDR_THRESHOLD

    # Add jitter to points with identical coordinates (linkage groups)
    rng = np.random.default_rng(42)
    df["log2_or_j"]    = df["log2_or"]    + rng.uniform(-0.08, 0.08, len(df))
    df["neg_log_fdr_j"] = df["neg_log_fdr"] + rng.uniform(-0.02, 0.02, len(df))

    # Colour
    df["color"] = df[group_col].map(color_map).fillna("#888888")

    return df


def make_volcano(df: pd.DataFrame,
                 group_col: str,
                 color_map: dict,
                 title: str,
                 out_prefix: Path,
                 figsize=(11, 7.5)):

    data = prep_gene_level(df, group_col, color_map)
    df_ns  = data[~data["significant"]]
    df_sig = data[data["significant"]]

    fig, ax = plt.subplots(figsize=figsize, facecolor="white")

    # ── Non-significant background ────────────────────────────────────────────
    ax.scatter(df_ns["log2_or_j"], df_ns["neg_log_fdr_j"],
               c="#e0e0e0", s=12, alpha=0.5, linewidths=0,
               rasterized=True, zorder=1)

    # ── Significant — coloured by class/family ────────────────────────────────
    for group in sorted(df_sig[group_col].unique()):
        sub   = df_sig[df_sig[group_col] == group]
        color = color_map.get(group, "#888888")
        ax.scatter(sub["log2_or_j"], sub["neg_log_fdr_j"],
                   c=color, s=40, alpha=0.85,
                   linewidths=0.5, edgecolors="white",
                   label=group, zorder=3)

    # ── Threshold line ────────────────────────────────────────────────────────
    fdr_line = -np.log10(FDR_THRESHOLD)
    ax.axhline(fdr_line, color="#d73027", lw=1.2, ls="--",
               alpha=0.8, zorder=2)
    ax.axvline(0, color="#777777", lw=0.8, ls="-", alpha=0.4, zorder=2)

    # ── Gene labels for top hits ──────────────────────────────────────────────
    if len(df_sig) > 0:
        # Pick top genes: highest -log10(FDR), then highest OR
        top = (df_sig
               .sort_values(["neg_log_fdr", "log2_or"], ascending=[False, False])
               .drop_duplicates("gene")
               .head(TOP_LABELS))

        texts = []
        for _, row in top.iterrows():
            gene_label = row["gene"][:28]
            t = ax.text(row["log2_or_j"], row["neg_log_fdr_j"],
                        gene_label, fontsize=6.5, color="#1a1a1a",
                        fontfamily="monospace", zorder=5)
            texts.append(t)

        if HAS_ADJUSTTEXT and texts:
            adjust_text(
                texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.6),
                expand_points=(1.8, 1.8),
                expand_text=(1.5, 1.5),
                force_text=(0.5, 0.8),
            )

    # ── Axes ──────────────────────────────────────────────────────────────────
    ax.set_xlabel("log₂(Odds Ratio)", fontsize=11)
    ax.set_ylabel("−log₁₀(FDR)", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=14)

    # Sensible axis limits
    xmax = max(data["log2_or"].max() * 1.12, 1)
    ymax = max(data["neg_log_fdr"].max() * 1.15, fdr_line * 1.5)
    ax.set_xlim(-0.3, xmax)
    ax.set_ylim(-0.1, ymax)

    # FDR line label
    ax.text(xmax * 0.01, fdr_line + ymax * 0.015,
            f"FDR = {FDR_THRESHOLD}", fontsize=7.5,
            color="#d73027", va="bottom")

    # Stats box
    n_sig  = len(df_sig)
    n_tot  = len(data)
    n_genes_sig = df_sig["gene"].nunique()
    ax.text(0.02, 0.97,
            f"Gene-level tests: {n_tot:,}\n"
            f"Significant genes: {n_genes_sig:,} (FDR<{FDR_THRESHOLD})",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec="#cccccc", alpha=0.9))

    # ── Legend ────────────────────────────────────────────────────────────────
    shown = sorted(df_sig[group_col].unique())
    handles = [mpatches.Patch(color=color_map.get(g, "#888"), label=g)
               for g in shown]
    ax.legend(handles=handles, title=group_col.capitalize(),
              fontsize=7.5, title_fontsize=8.5,
              loc="lower right", framealpha=0.9)

    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=9)
    fig.tight_layout()

    for ext in ["pdf", "png"]:
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")

    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Figure 2 — Volcano Plots (gene-level)")
    print("=" * 60)

    # Fig 2a — per antibiotic
    ab_file = RES / "gene_level_per_ab.tsv"
    if ab_file.exists():
        print("\nFig 2a — per-antibiotic...")
        df_ab = pd.read_csv(ab_file, sep="\t")
        df_ab["class"] = df_ab["antibiotic"].map(AB_CLASS).fillna("Other")
        print(f"  {len(df_ab):,} (gene, antibiotic) pairs loaded")

        make_volcano(
            df         = df_ab,
            group_col  = "class",
            color_map  = CLASS_COLORS,
            title      = (r"SNP–MIC gene-level associations in $\it{Vibrio\ vulnificus}$"
                          "\n(per antibiotic, n=29 strains with WGS)"),
            out_prefix = FIGURES / "Fig2a_volcano_per_ab",
        )
    else:
        print(f"  [SKIP] {ab_file} not found — run 04 first.")

    # Fig 2b — per family
    fam_file = RES / "gene_level_per_family.tsv"
    if fam_file.exists():
        print("\nFig 2b — per-family...")
        df_fam = pd.read_csv(fam_file, sep="\t")
        print(f"  {len(df_fam):,} (gene, family) pairs loaded")

        make_volcano(
            df         = df_fam,
            group_col  = "family",
            color_map  = FAMILY_COLORS,
            title      = (r"SNP–MIC gene-level associations in $\it{Vibrio\ vulnificus}$"
                          "\n(per antibiotic family, n=29 strains with WGS)"),
            out_prefix = FIGURES / "Fig2b_volcano_per_family",
        )
    else:
        print(f"  [SKIP] {fam_file} not found — run 04 first.")

    print("\nDone.")

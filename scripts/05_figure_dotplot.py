#!/usr/bin/env python3
"""
05_figure_dotplot.py
--------------------
Dot plot (lollipop chart) of top significant gene-antibiotic associations.
One row per gene, dot size = number of SNPs, colour = antibiotic class,
X-axis = log2(Odds Ratio), sorted by OR descending.

Also produces a family-level version.

Usage:
    python3 05_figure_dotplot.py

Inputs:
    ../results/gene_level_per_ab.tsv
    ../results/gene_level_per_family.tsv

Outputs:
    ../figures/Fig2a_dotplot_per_ab.pdf/png
    ../figures/Fig2b_dotplot_per_family.pdf/png

Author: Francesc Montoro (2025)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = Path(__file__).parent.parent
RES     = BASE / "results"
FIGURES = BASE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
FDR_THRESHOLD = 0.05
TOP_N         = 30      # top genes to show per plot
MAX_OR_PLOT   = 7.0     # log2 scale cap for display (2^7 = 128)

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
    "Aminoglycoside":   "#4393c3",
    "Beta-lactam":      "#4dac26",
    "Beta-lactam combo":"#1a7a1a",
    "Carbapenem":       "#006d2c",
    "Fluoroquinolone":  "#e08214",
    "Macrolide":        "#8c510a",
    "Lincosamide":      "#543005",
    "Tetracycline":     "#d6604d",
    "Sulfonamide":      "#9970ab",
    "Phenicol":         "#c994c7",
    "Other":            "#888888",
}

FAMILY_COLORS = {
    "Aminoglycoside":   "#4393c3",
    "Beta-lactam":      "#4dac26",
    "Carbapenem":       "#006d2c",
    "Fluoroquinolone":  "#e08214",
    "Macrolide":        "#8c510a",
    "Lincosamide":      "#543005",
    "Tetracycline":     "#d6604d",
    "Sulfonamide":      "#9970ab",
    "Phenicol":         "#c994c7",
    "Other":            "#888888",
}


def clean_gene_name(name: str, max_len: int = 35) -> str:
    """Shorten long gene names for display."""
    name = name.replace("_", " ")
    if len(name) > max_len:
        name = name[:max_len - 1] + "…"
    return name


def make_dotplot(df: pd.DataFrame,
                 group_col: str,
                 color_map: dict,
                 title: str,
                 out_prefix: Path,
                 top_n: int = TOP_N):
    """
    Lollipop / dot plot.
    Rows = top genes (by OR, then FDR).
    X = log2(OR), dot size = n_snps, colour = class/family.
    """
    # Filter significant only
    sig = df[df["best_fdr"] < FDR_THRESHOLD].copy()
    if len(sig) == 0:
        print(f"  [SKIP] No significant genes for {out_prefix.name}")
        return

    # log2 OR — cap inf (stored as 999) at MAX_OR_PLOT
    sig["or_plot"]   = sig["best_OR"].clip(upper=2**MAX_OR_PLOT)
    sig["log2_or"]   = np.log2(sig["or_plot"].clip(lower=1e-6)).clip(upper=MAX_OR_PLOT)
    sig["is_inf"]    = sig["best_OR"] >= 999    # flag capped values

    # Sort: highest OR first, then lowest FDR
    sig = (sig.sort_values(["log2_or", "best_fdr"],
                            ascending=[False, True])
              .drop_duplicates("gene")
              .head(top_n))

    # Reverse for bottom-to-top display
    sig = sig.iloc[::-1].reset_index(drop=True)

    n = len(sig)
    fig_h = max(5, n * 0.38)
    fig, ax = plt.subplots(figsize=(11, fig_h), facecolor="white")

    # Y positions
    y_pos = np.arange(n)

    # ── Lollipop stems ────────────────────────────────────────────────────────
    for i, row in sig.iterrows():
        color = color_map.get(row[group_col], "#888888")
        ax.plot([0, row["log2_or"]], [i, i],
                color=color, lw=1.5, alpha=0.6, zorder=1)

    # ── Dots ──────────────────────────────────────────────────────────────────
    # Size proportional to n_snps, capped
    sizes = np.clip(sig["n_snps"] * 8, 40, 300)
    colors = [color_map.get(row[group_col], "#888888")
              for _, row in sig.iterrows()]

    sc = ax.scatter(sig["log2_or"], y_pos,
                    s=sizes, c=colors, alpha=0.9,
                    edgecolors="white", linewidths=0.8,
                    zorder=3)

    # Mark capped (inf) values with a triangle marker
    inf_mask = sig["is_inf"]
    if inf_mask.any():
        ax.scatter(sig.loc[inf_mask, "log2_or"],
                   y_pos[inf_mask.values],
                   s=sizes[inf_mask.values] * 0.4,
                   marker="^", c="white",
                   zorder=4, linewidths=0)

    # ── FDR annotations ───────────────────────────────────────────────────────
    for i, row in sig.iterrows():
        fdr_str = f"FDR={row['best_fdr']:.3f}"
        ax.text(row["log2_or"] + 0.08, i, fdr_str,
                va="center", ha="left", fontsize=6.5,
                color="#555555", fontfamily="monospace")

    # ── Y labels ─────────────────────────────────────────────────────────────
    gene_labels = [clean_gene_name(g) for g in sig["gene"]]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(gene_labels, fontsize=8, fontfamily="monospace")

    # Colour Y labels by class
    for tick, (_, row) in zip(ax.get_yticklabels(), sig.iterrows()):
        tick.set_color(color_map.get(row[group_col], "#333333"))

    # ── Axes ──────────────────────────────────────────────────────────────────
    ax.set_xlabel("log₂(Odds Ratio)", fontsize=11)
    ax.set_xlim(-0.2, MAX_OR_PLOT + 1.5)
    ax.set_ylim(-0.8, n - 0.2)
    ax.axvline(0, color="#888888", lw=0.8, ls="-", alpha=0.4)

    # Inf boundary line
    ax.axvline(MAX_OR_PLOT, color="#aaaaaa", lw=0.8, ls=":",
               alpha=0.6, zorder=0)
    ax.text(MAX_OR_PLOT + 0.05, n - 1.2, "OR = ∞",
            fontsize=7, color="#888888", rotation=90, va="top")

    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#eeeeee", lw=0.8, zorder=0)

    # ── Legends ───────────────────────────────────────────────────────────────
    # Colour legend (class/family)
    shown = sorted(sig[group_col].unique())
    color_handles = [mpatches.Patch(color=color_map.get(g, "#888"), label=g)
                     for g in shown]

    # Size legend (n_snps)
    size_vals  = [1, 5, 20]
    size_handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="#888888",
               markersize=np.sqrt(np.clip(v * 8, 40, 300) / np.pi),
               label=f"{v} SNPs")
        for v in size_vals
    ]
    size_handles.append(
        Line2D([0], [0], marker="^", color="w",
               markerfacecolor="white",
               markeredgecolor="#888888",
               markersize=7,
               label="OR = ∞ (capped)")
    )

    leg1 = ax.legend(handles=color_handles,
                     title=group_col.capitalize(),
                     fontsize=7.5, title_fontsize=8.5,
                     loc="lower right",
                     bbox_to_anchor=(1.0, 0.0),
                     framealpha=0.9)
    ax.add_artist(leg1)
    ax.legend(handles=size_handles,
              title="Dot size",
              fontsize=7.5, title_fontsize=8.5,
              loc="lower right",
              bbox_to_anchor=(1.0, 0.32),
              framealpha=0.9)

    fig.tight_layout()

    for ext in ["pdf", "png"]:
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")

    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Figure 2 — Dot plots (top significant genes)")
    print("=" * 60)

    # Fig 2a — per antibiotic
    ab_file = RES / "gene_level_per_ab.tsv"
    if ab_file.exists():
        print("\nFig 2a — per-antibiotic...")
        df_ab = pd.read_csv(ab_file, sep="\t")
        df_ab["class"] = df_ab["antibiotic"].map(AB_CLASS).fillna("Other")
        make_dotplot(
            df         = df_ab,
            group_col  = "class",
            color_map  = CLASS_COLORS,
            title      = (r"Top SNP-associated genes in $\it{Vibrio\ vulnificus}$"
                          "\n(per antibiotic, n=29 strains with WGS)"),
            out_prefix = FIGURES / "Fig2a_dotplot_per_ab",
        )

    # Fig 2b — per family
    fam_file = RES / "gene_level_per_family.tsv"
    if fam_file.exists():
        print("\nFig 2b — per-family...")
        df_fam = pd.read_csv(fam_file, sep="\t")
        make_dotplot(
            df         = df_fam,
            group_col  = "family",
            color_map  = FAMILY_COLORS,
            title      = (r"Top SNP-associated genes in $\it{Vibrio\ vulnificus}$"
                          "\n(per antibiotic family, n=29 strains with WGS)"),
            out_prefix = FIGURES / "Fig2b_dotplot_per_family",
        )

    print("\nDone.")

#!/usr/bin/env python3
"""
06_figure_phylo.py
------------------
Publication-quality phylogenetic tree of the 30 V. vulnificus strains
with WGS data, annotated with:
  - Lineage colour (tip labels and tip dots)
  - Resistance score heatmap strip (per antibiotic) alongside the tree

Usage:
    python3 06_figure_phylo.py

Inputs:
    ../data/raw/phylo_tree.nwk
    ../data/processed/MICs_converted.tsv
    ../data/raw/lineage_manual_mapping.tsv

Outputs:
    ../figures/Fig3_phylo_resistance.pdf
    ../figures/Fig3_phylo_resistance.png

Author: Francesc Montoro (2025)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

try:
    from Bio import Phylo
    from io import StringIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("  [ERROR] Biopython not found. Install with:")
    print("  pip install biopython --break-system-packages")
    exit(1)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = Path(__file__).parent.parent
RAW     = BASE / "data" / "raw"
PROC    = BASE / "data" / "processed"
FIGURES = BASE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

# ── Name mapping: tree names → MIC names ─────────────────────────────────────
TREE_TO_MIC = {
    "11028":      "BT3_11028",  "12":        "BT3_12",
    "162":        "BT3_162",    "32":        "BT3_32",
    "94385":      "94-3-85",    "960426-1_4C": "960426_1/4C",
    "CIP8190":    "C1P8190",    "V4":        "V_4",
    "CECT4604":   "CECT_4604",  "CECT4606":  "CECT_4606",
    "CECT4608":   "CECT_4608",  "CECT4865":  "CECT_4865",
    "CECT4866":   "CECT_4866",  "CECT4999":  "CECT_4999",
    "CECT5198":   "CECT_5198",  "CECT5763":  "CECT_5763",
    "CECT5769":   "CECT_5769",  "CECT7030":  "CECT_7030",
    "CECT898":    "CECT_898",   "CG100":     "CG_100",
    "94-9-112":   "94-9-112",   "94-9-119":  "94-9-119",
    "95-8-161":   "95-8-161",   "95-8-6":    "95-8-6",
    "95-8-7":     "95-8-7",     "A14":       "A14",
    "PD-2-51":    "PD-2-51",    "Ra_3":      "Ra_3",
    "RIU1":       "RIU1",       "YJ016":     "YJ016",
}

# Antibiotics to show in heatmap strip (select most informative)
AB_SHOW = ["GEN2", "STZ", "SDM", "NOV", "SPE", "STR",
           "CLI",  "ERY", "PEN", "AMOX", "OXY", "TET",
           "FAZ",  "FOX", "TYLT"]

AB_NAMES = {
    "GEN2":"Gentamicin (low)", "STZ":"Sulphathiazole",
    "SDM":"Sulphadimethoxine", "NOV":"Novobiocin",
    "SPE":"Spectinomycin",     "STR":"Streptomycin",
    "CLI":"Clindamycin",       "ERY":"Erythromycin",
    "PEN":"Penicillin",        "AMOX":"Amoxicillin",
    "OXY":"Oxytetracycline",   "TET":"Tetracycline",
    "FAZ":"Cefazolin",         "FOX":"Cefoxitin",
    "TYLT":"Tylosin tartrate",
}

LINEAGE_COLORS = {
    "L1":"#2166ac", "1":"#2166ac",
    "L2":"#74add1", "2":"#74add1",
    "L3":"#d73027", "3":"#d73027",
    "L4":"#f4a582", "4":"#f4a582",
    "unknown":"#aaaaaa",
}

MIC_CMAP = LinearSegmentedColormap.from_list(
    "mic", ["#f7f7f7", "#fee090", "#fc8d59", "#d73027", "#67001f"], N=256
)


def load_lineages(raw_dir, strains):
    manual = raw_dir / "lineage_manual_mapping.tsv"
    if manual.exists():
        df = pd.read_csv(manual, sep="\t")
        df["Strain"] = df["Strain"].astype(str).str.strip()
        lut = df.set_index("Strain")["Lineage"].astype(str)
        return {s: lut.get(s, "unknown") for s in strains}
    return {s: "unknown" for s in strains}


def get_tip_order(tree):
    """Return tip names in display order (top to bottom) from Biopython tree."""
    tips = []
    def traverse(clade):
        if clade.is_terminal():
            tips.append(clade.name)
        else:
            for child in clade.clades:
                traverse(child)
    traverse(tree.root)
    return tips


def draw_tree(ax, tree, tip_y: dict, lineage_colors: dict):
    """Draw the phylogenetic tree on ax. tip_y maps tip name → y position."""

    def get_y(clade):
        if clade.is_terminal():
            return tip_y[clade.name]
        ys = [get_y(c) for c in clade.clades]
        return np.mean([min(ys), max(ys)])

    x_pos = {}
    def calc_x(clade, x=0):
        x_pos[clade] = x
        if not clade.is_terminal():
            for child in clade.clades:
                bl = child.branch_length or 0.01
                calc_x(child, x + bl)
    calc_x(tree.root)

    def draw_clade(clade):
        x = x_pos[clade]
        y = get_y(clade)
        if not clade.is_terminal():
            child_ys = [get_y(c) for c in clade.clades]
            # Vertical line connecting children
            ax.plot([x, x], [min(child_ys), max(child_ys)],
                    color="#444444", lw=1.0, solid_capstyle="round")
            for child in clade.clades:
                cx = x_pos[child]
                cy = get_y(child)
                # Horizontal branch
                ax.plot([x, cx], [cy, cy],
                        color="#444444", lw=1.0, solid_capstyle="round")
                draw_clade(child)
        else:
            # Tip dot coloured by lineage
            lin = lineage_colors.get(clade.name, "#aaaaaa")
            ax.scatter(x, y, s=30, color=lin, zorder=5,
                       edgecolors="white", linewidths=0.5)

    draw_clade(tree.root)

    # Scale bar
    x_vals = [x_pos[c] for c in x_pos]
    x_max = max(x_vals)
    scale = round(x_max * 0.2, 2)
    ax.plot([0, scale], [-1.5, -1.5], color="#444444", lw=1.5)
    ax.text(scale / 2, -2.0, f"{scale}", ha="center", fontsize=7, color="#444444")
    ax.text(0, -2.8, "Branch length", ha="left", fontsize=6.5, color="#888888")

    return x_max


def make_phylo_figure(tree_file, mics, lineage_map, out_prefix):

    # ── Load tree ─────────────────────────────────────────────────────────────
    with open(tree_file) as f:
        nwk = f.read().strip()
    tree = Phylo.read(StringIO(nwk), "newick")

    # Rename tips to MIC names
    for tip in tree.get_terminals():
        tip.name = TREE_TO_MIC.get(tip.name, tip.name)

    tips = get_tip_order(tree)
    n    = len(tips)
    tip_y = {tip: i for i, tip in enumerate(tips)}

    # ── Lineage colours per tip ───────────────────────────────────────────────
    tip_lineage_color = {
        tip: LINEAGE_COLORS.get(lineage_map.get(tip, "unknown"), "#aaaaaa")
        for tip in tips
    }

    # ── MIC data for these strains ────────────────────────────────────────────
    ab_cols = [ab for ab in AB_SHOW if ab in mics.columns]
    mic_sub = mics.reindex(tips)[ab_cols].fillna(0)

    # ── Figure layout ─────────────────────────────────────────────────────────
    n_ab   = len(ab_cols)
    fig_w  = 6 + n_ab * 0.35
    fig_h  = max(8, n * 0.32)

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    # Tree axes (left 40%)
    ax_tree = fig.add_axes([0.01, 0.08, 0.28, 0.88])
    # Lineage strip
    ax_lin  = fig.add_axes([0.38, 0.08, 0.012, 0.88])
    # Heatmap axes (right)
    ax_heat = fig.add_axes([0.41, 0.08, 0.50, 0.88])
    # Colorbar
    ax_cbar = fig.add_axes([0.93, 0.08, 0.012, 0.88])

    # ── Draw tree ─────────────────────────────────────────────────────────────
    x_max = draw_tree(ax_tree, tree, tip_y, tip_lineage_color)
    ax_tree.set_xlim(-x_max * 0.02, x_max * 1.05)
    ax_tree.set_ylim(-3.5, n - 0.5)
    ax_tree.axis("off")

    # Tip labels coloured by lineage
    for tip in tips:
        y = tip_y[tip]
        color = tip_lineage_color[tip]
        ax_tree.text(x_max * 1.08, y, tip.replace("_", " "),
                     va="center", ha="left", fontsize=6.5,
                     color=color, fontfamily="monospace")

    # ── Lineage strip ─────────────────────────────────────────────────────────
    lin_rgb = np.array(
        [[matplotlib.colors.to_rgb(tip_lineage_color[tip])] for tip in tips]
    ).reshape(n, 1, 3)
    ax_lin.imshow(lin_rgb, aspect="auto", interpolation="nearest",
                  origin="lower")
    ax_lin.set_xticks([])
    ax_lin.set_yticks([])

    # ── Heatmap ───────────────────────────────────────────────────────────────
    mat = mic_sub.values   # shape: n_strains × n_ab
    im  = ax_heat.imshow(mat, aspect="auto", cmap=MIC_CMAP,
                         vmin=0, vmax=5, interpolation="nearest",
                         origin="lower")

    ax_heat.set_xticks(range(n_ab))
    ax_heat.set_xticklabels([AB_NAMES.get(ab, ab) for ab in ab_cols],
                             rotation=45, ha="right", fontsize=7.5)
    ax_heat.set_yticks([])
    ax_heat.tick_params(axis="x", length=0)

    # Minor grid
    ax_heat.set_xticks(np.arange(-0.5, n_ab, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n, 1),    minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.6)
    ax_heat.tick_params(which="minor", length=0)

    # ── Colorbar ──────────────────────────────────────────────────────────────
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label("Resistance score [0–5]", fontsize=8,
                   rotation=270, labelpad=14)
    cbar.set_ticks([0, 1, 2, 3, 4, 5])
    cbar.ax.tick_params(labelsize=7)

    # ── Lineage legend ────────────────────────────────────────────────────────
    lin_labels = {"1":"Lineage 1","2":"Lineage 2",
                  "3":"Lineage 3","4":"Lineage 4"}
    shown_lins = sorted({lineage_map.get(t, "unknown") for t in tips}
                        - {"unknown"})
    patches = [mpatches.Patch(color=LINEAGE_COLORS.get(l, "#aaa"),
                               label=lin_labels.get(l, f"L{l}"))
               for l in shown_lins]
    ax_tree.legend(handles=patches, loc="lower left",
                   fontsize=7, title="Lineage", title_fontsize=8,
                   framealpha=0.9)

    # ── Title ─────────────────────────────────────────────────────────────────
    fig.suptitle(
        r"Phylogeny and antibiotic resistance of $\it{Vibrio\ vulnificus}$"
        " (n=30 strains with WGS)",
        fontsize=11, y=0.99, fontweight="bold"
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext in ["pdf", "png"]:
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=900, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")

    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Figure 3 — Phylogenetic tree + resistance heatmap")
    print("=" * 60)

    # Load data
    mics = pd.read_csv(PROC / "MICs_converted.tsv", sep="\t", index_col=0)
    mics.index = mics.index.astype(str).str.strip()

    lineage_df = pd.read_csv(RAW / "lineage_manual_mapping.tsv", sep="\t")
    lineage_df["Strain"] = lineage_df["Strain"].astype(str).str.strip()
    lineage_map = dict(zip(lineage_df["Strain"],
                           lineage_df["Lineage"].astype(str)))

    tree_file = RAW / "phylo_tree.nwk"
    print(f"Tree: {tree_file}")
    print(f"MICs: {mics.shape}")

    make_phylo_figure(tree_file, mics, lineage_map,
                      FIGURES / "Fig3_phylo_resistance")
    print("\nDone.")

#!/usr/bin/env python3
"""
02_figure_heatmap.py
--------------------
Publication-quality heatmap of converted MIC data.
Layout: strains on X-axis, antibiotics on Y-axis.
Both axes hierarchically clustered.
AB class shown via Y-label colour. Lineage shown via top colour strip.

Usage:
    python3 02_figure_heatmap.py

Inputs:
    ../data/processed/MICs_converted.tsv
    ../data/raw/*lineage*.xlsx   (optional)

Outputs:
    ../figures/Fig1_MIC_heatmap.pdf
    ../figures/Fig1_MIC_heatmap.png

Author: Francesc Montoro (2025)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = Path(__file__).parent.parent
PROC    = BASE / "data" / "processed"
RAW     = BASE / "data" / "raw"
FIGURES = BASE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

# ── Antibiotic full names ──────────────────────────────────────────────────────
AB_NAMES = {
    "AMI":  "Amikacin",           "AMP":  "Ampicillin",
    "A/S2": "Amox/Sulbactam",     "AZT":  "Aztreonam",
    "FAZ":  "Cefazolin",          "FEP":  "Cefepime",
    "TANS": "Trimethoprim",       "AXO":  "Ceftriaxone",
    "TAZ":  "Ceftazidime",        "FUR":  "Nitrofurantoin",
    "CIP":  "Ciprofloxacin",      "GAT":  "Gatifloxacin",
    "MERO": "Meropenem",          "GEN":  "Gentamicin",
    "IMI":  "Imipenem",           "NIT":  "Nitrofurantoin B",
    "FOX":  "Cefoxitin",          "PIP":  "Piperacillin",
    "P/T4": "Pip/Tazobactam",     "TIM2": "Ticar/Clavulanate",
    "SXT":  "Trim/Sulfa",         "TOB":  "Tobramycin",
    "POD":  "Cefpodoxime",        "ENRO": "Enrofloxacin",
    "GEN2": "Gentamicin (low)",   "TIO":  "Ceftiofur",
    "NEO":  "Neomycin",           "ERY":  "Erythromycin",
    "OXY":  "Oxytetracycline",    "TET":  "Tetracycline",
    "AMOX": "Amoxicillin",        "SPE":  "Spectinomycin",
    "SDM":  "Sulphadimethoxine",  "SXT2": "Trim/Sulfa B",
    "FFN":  "Florfenicol",        "STZ":  "Sulphathiazole",
    "PEN":  "Penicillin",         "STR":  "Streptomycin",
    "NOV":  "Novobiocin",         "TYLT": "Tylosin tartrate",
    "CLI":  "Clindamycin",
}

AB_CLASS = {
    "AMI":"Aminoglycoside",  "GEN":"Aminoglycoside",  "GEN2":"Aminoglycoside",
    "NEO":"Aminoglycoside",  "TOB":"Aminoglycoside",  "STR":"Aminoglycoside",
    "SPE":"Aminoglycoside",
    "AMP":"Beta-lactam",     "FAZ":"Beta-lactam",     "FEP":"Beta-lactam",
    "FOX":"Beta-lactam",     "AXO":"Beta-lactam",     "TAZ":"Beta-lactam",
    "POD":"Beta-lactam",     "PIP":"Beta-lactam",     "AMOX":"Beta-lactam",
    "PEN":"Beta-lactam",     "TIO":"Beta-lactam",     "AZT":"Beta-lactam",
    "A/S2":"Beta-lactam combo","P/T4":"Beta-lactam combo","TIM2":"Beta-lactam combo",
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
    "Aminoglycoside":    "#4393c3",
    "Beta-lactam":       "#4dac26",
    "Beta-lactam combo": "#1a7a1a",
    "Carbapenem":        "#006d2c",
    "Fluoroquinolone":   "#e08214",
    "Macrolide":         "#8c510a",
    "Lincosamide":       "#543005",
    "Tetracycline":      "#d6604d",
    "Sulfonamide":       "#9970ab",
    "Phenicol":          "#c994c7",
    "Other":             "#888888",
}

LINEAGE_COLORS = {
    "L1": "#2166ac", "1": "#2166ac",
    "L2": "#74add1", "2": "#74add1",
    "L3": "#d73027", "3": "#d73027",
    "L4": "#f4a582", "4": "#f4a582",
    "unknown": "#dddddd",
}


def hierarchical_order(mat, axis=0):
    data = mat if axis == 0 else mat.T
    dist = pdist(data, metric="euclidean")
    link = linkage(dist, method="ward")
    return leaves_list(link)


def load_lineages(raw_dir, strains):
    # Check for manual mapping first
    manual = raw_dir / "lineage_manual_mapping.tsv"
    if manual.exists():
        df = pd.read_csv(manual, sep="\t")
        df["Strain"] = df["Strain"].astype(str).str.strip()
        lut = df.set_index("Strain")["Lineage"].astype(str)
        return pd.Series([lut.get(s, "unknown") for s in strains], index=strains)    
    candidates = list(raw_dir.glob("*lineage*")) + list(raw_dir.glob("*Lineage*"))
    if not candidates:
        return pd.Series("unknown", index=strains)
    try:
        df = pd.read_excel(candidates[0])
        strain_col  = next((c for c in df.columns if "strain" in c.lower()), df.columns[0])
        lineage_col = next((c for c in df.columns
                            if "lineage" in c.lower() or "linaje" in c.lower()), None)
        if lineage_col is None:
            return pd.Series("unknown", index=strains)
        df[strain_col] = df[strain_col].astype(str).str.strip()
        lut = df.set_index(strain_col)[lineage_col].astype(str)
        # Normalise: remove underscores, spaces, slashes for fuzzy matching
        def norm(x): return str(x).upper().replace("_","").replace("-","").replace("/","").replace("*","").replace(" ","")
        lut_norm = {norm(k): v for k, v in lut.items()}
        def add_prefix(v): return f"L{v}" if v not in ("unknown", "5", "5?") else "unknown"
        return pd.Series([add_prefix(lut.get(s, "unknown")) for s in strains], index=strains)
    except Exception as e:
        print(f"  [INFO] Lineage file not loaded: {e}")
        return pd.Series("unknown", index=strains)


def make_heatmap(df: pd.DataFrame, lineages: pd.Series, out_prefix: Path):

    # ── Cluster ───────────────────────────────────────────────────────────────
    col_order = hierarchical_order(df.values,   axis=0)   # strains  → X
    row_order = hierarchical_order(df.values.T, axis=0)   # antibiotics → Y

    # mat: rows=antibiotics, cols=strains
    mat         = df.iloc[col_order, :].T.iloc[row_order, :]
    lin_ordered = lineages.iloc[col_order]
    ab_ordered  = df.columns[row_order]

    n_strains = mat.shape[1]
    n_abs     = mat.shape[0]

    # ── Figure ────────────────────────────────────────────────────────────────
    fig_w = max(22, n_strains * 0.19)
    fig_h = max(9,  n_abs    * 0.32)

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    # Single heatmap axes — generous left margin for Y labels
    ax_heat = fig.add_axes([0.28, 0.22, 0.60, 0.62])
    ax_lin  = fig.add_axes([0.28, 0.845, 0.60, 0.013])   # lineage strip top
    ax_cbar = fig.add_axes([0.90, 0.22, 0.012, 0.62])    # colorbar right

    # ── Colour map ────────────────────────────────────────────────────────────
    cmap = LinearSegmentedColormap.from_list(
        "mic", ["#f7f7f7", "#fee090", "#fc8d59", "#d73027", "#67001f"], N=256
    )

    # ── Heatmap ───────────────────────────────────────────────────────────────
    im = ax_heat.imshow(mat.values, aspect="auto", cmap=cmap,
                        vmin=0, vmax=5, interpolation="nearest")

    # X: strain labels
    ax_heat.set_xticks(range(n_strains))
    ax_heat.set_xticklabels(mat.columns, rotation=50, ha="right",
                             fontsize=6.5, fontfamily="monospace")
    ax_heat.tick_params(axis="x", length=0)

    # Y: antibiotic labels — coloured by class
    ab_labels = [AB_NAMES.get(ab, ab) for ab in ab_ordered]
    ax_heat.set_yticks(range(n_abs))
    ax_heat.set_yticklabels(ab_labels, fontsize=8)
    ax_heat.tick_params(axis="y", length=0)

    # Colour each Y tick label by antibiotic class
    for tick, ab in zip(ax_heat.get_yticklabels(), ab_ordered):
        cls   = AB_CLASS.get(ab, "Other")
        color = CLASS_COLORS.get(cls, "#888888")
        tick.set_color(color)
        tick.set_fontweight("bold")

    # Minor grid
    ax_heat.set_xticks(np.arange(-0.5, n_strains, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n_abs,     1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.5)
    ax_heat.tick_params(which="minor", length=0)

    # ── Lineage strip (top) ───────────────────────────────────────────────────
    lin_rgb = np.array(
        [[matplotlib.colors.to_rgb(LINEAGE_COLORS.get(str(l), "#dddddd"))]
         for l in lin_ordered]
    ).reshape(1, n_strains, 3)
    #print("Lineage sample:", lin_ordered[:10].values)
    #print("Colors sample:", lin_colors[0, :5, :])
    ax_lin.imshow(lin_rgb, aspect="auto", interpolation="nearest")
    ax_lin.set_xticks([])
    ax_lin.set_yticks([])
    ax_lin.set_ylabel("Lineage", fontsize=7, rotation=0, labelpad=38, va="center")

    # ── Colorbar ──────────────────────────────────────────────────────────────
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label("Resistance score [0–5]", fontsize=8, rotation=270, labelpad=14)
    cbar.set_ticks([0, 1, 2, 3, 4, 5])
    cbar.ax.tick_params(labelsize=7)

    # ── Legends ───────────────────────────────────────────────────────────────
    # Lineage
    unique_lins = sorted(set(lin_ordered) - {"unknown"})
    lin_patches = [mpatches.Patch(color=LINEAGE_COLORS.get(l, "#ddd"), label=l)
                   for l in unique_lins]
    if "unknown" in set(lin_ordered):
        lin_patches.append(mpatches.Patch(color="#dddddd", label="Unknown"))

    # AB class
    shown_classes = sorted({AB_CLASS.get(ab, "Other") for ab in ab_ordered})
    class_patches = [mpatches.Patch(color=CLASS_COLORS.get(c, "#aaa"), label=c)
                     for c in shown_classes]

    leg1 = ax_heat.legend(handles=lin_patches, title="Lineage",
                          loc="upper left", bbox_to_anchor=(1.07, 1.10),
                          fontsize=7, title_fontsize=8, framealpha=0.9)
    ax_heat.add_artist(leg1)
    ax_heat.legend(handles=class_patches, title="Antibiotic class",
                   loc="lower left", bbox_to_anchor=(1.07, 0.0),
                   fontsize=7, title_fontsize=8, framealpha=0.9)

    # ── Title ─────────────────────────────────────────────────────────────────
    fig.suptitle(
        r"Antibiotic resistance profiles of $\it{Vibrio\ vulnificus}$"
        " (n=116 strains, 41 antibiotics)",
        fontsize=12, y=0.99, fontweight="bold"
    )

    # ── Save ──────────────────────────────────────────────────────────────────
    for ext in ["pdf", "png"]:
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")

    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Figure 1 — MIC Heatmap (strains on X, AB on Y)")
    print("=" * 60)

    df = pd.read_csv(PROC / "MICs_converted.tsv", sep="\t", index_col=0)
    df = df.fillna(0)
    print(f"Loaded: {df.shape[0]} strains × {df.shape[1]} antibiotics")

    lineages = load_lineages(RAW, df.index.tolist())
    n_known  = (lineages != "unknown").sum()
    print(f"Lineage annotation: {n_known}/{len(lineages)} strains annotated")

    make_heatmap(df, lineages, FIGURES / "Fig1_MIC_heatmap")
    print("\nDone.")

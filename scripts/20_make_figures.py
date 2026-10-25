#!/usr/bin/env python3
"""
20_make_figures.py
==================
Generates five of the seven main figures of the manuscript (numbering as
published; Fig. 1 is a schematic and Fig. 7 is drawn by
23_figure_calibration.py).

  Fig 2  MIC distributions for gentamicin and sulphathiazole across the full
         116-isolate collection, with both candidate cut-offs marked.
  Fig 3  Acquired determinants by strain, with the corresponding MIC.
  Fig 4  Collapse of nominally significant variants under phylogenetic control.
  Fig 5  Attainable significance: combinatorial and phylogenetic bounds.
  Fig 6  The contig_2 variant: distribution over the tree against phenotype.

Output is vector PDF at 300 dpi equivalent, using a colour-blind-safe palette
and a single typeface throughout. Each figure is also written as PNG for
drafting.

Usage:
  python3 20_make_figures.py \
      --mic-xlsx  data/mic/MICs_all_Vv.xlsx \
      --genopheno results/genotype_phenotype/genotype_phenotype_check.tsv \
      --tree      results/pyseer/core.treefile \
      --matrix    results/summary/variants_matrix.tsv \
      --outdir    figures
"""

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ── house style ─────────────────────────────────────────────────────────────
# pdf.fonttype/ps.fonttype 42 forces TrueType rather than matplotlib's default
# Type 3 fonts, which several publishers reject at the production stage because
# they do not embed reliably.
plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "figure.dpi": 150,
    "savefig.dpi": 900,          # line art, as specified by the co-authors
    "savefig.bbox": "tight",
})

# Smallest type used anywhere in these figures. Journals typically require
# 6-8 pt as an absolute floor at final size; nothing here goes below it.
MINPT = 7.0

# Okabe-Ito palette (colour-blind safe)
C = {
    "blue":   "#0072B2",
    "orange": "#E69F00",
    "green":  "#009E73",
    "red":    "#D55E00",
    "purple": "#CC79A7",
    "sky":    "#56B4E9",
    "yellow": "#F0E442",
    "grey":   "#999999",
    "dark":   "#333333",
}

MM = 1 / 25.4          # mm to inches
SINGLE = 85 * MM       # single-column width
DOUBLE = 180 * MM      # double-column width


def save(fig, outdir, name):
    """PDF for submission (vector), PNG at 900 dpi for drafting and review."""
    fig.savefig(Path(outdir) / f"{name}.pdf")
    fig.savefig(Path(outdir) / f"{name}.png")
    fig.savefig(Path(outdir) / f"{name}.svg")
    plt.close(fig)
    print(f"  wrote {name}.pdf (vector) / {name}.png (900 dpi)")


# Column names of agents tested on both plates. In the workbook the duplicated
# headers are read by pandas as GEN/GEN.1 and SXT/SXT.1; the Zenodo release may
# instead carry the plate as a suffix. Both forms are accepted.
MIC_COLUMNS = {"GEN_GN2F": "GEN", "GEN_AVIAN1F": "GEN.1",
               "SXT_GN2F": "SXT", "SXT_AVIAN1F": "SXT.1"}


def read_mics(mic_xlsx):
    return pd.read_excel(mic_xlsx, sheet_name="MICs").rename(columns=MIC_COLUMNS)


# ════════════════════════════════════════════════════════════════════════════
# Figure 1 — MIC distributions
# ════════════════════════════════════════════════════════════════════════════
def figure1(mic_xlsx, outdir):
    d = read_mics(mic_xlsx)

    gen_levels = ["<0.5", "0.5", "1", "2", "4"]
    stz_levels = ["<32", "64", "128", "256", ">256"]

    gen = d["GEN.1"].astype(str).str.strip().value_counts()
    stz = d["STZ"].astype(str).str.strip().value_counts()

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 62 * MM))

    # ── gentamicin ──
    ax = axes[0]
    counts = [gen.get(l, 0) for l in gen_levels]
    bars = ax.bar(range(len(gen_levels)), counts, width=0.72,
                  color=C["grey"], edgecolor="none")
    # strains called resistant under the reference-anchored definition
    for i, l in enumerate(gen_levels):
        if l in ("1", "2", "4"):
            bars[i].set_color(C["orange"])
    ax.set_xticks(range(len(gen_levels)))
    ax.set_xticklabels(gen_levels)
    ax.set_xlabel("Gentamicin MIC (µg ml$^{-1}$)")
    ax.set_ylabel("Number of isolates")
    ax.set_title("a", loc="left", fontweight="bold")
    # bracket over the bars called resistant by a reference-anchored cut-off
    ymax = max(counts)
    ax.annotate("", xy=(1.62, ymax * 1.03), xytext=(4.38, ymax * 1.03),
                arrowprops=dict(arrowstyle="-", color=C["orange"], lw=1.0))
    ax.text(3.0, ymax * 1.06, "called resistant by a\nreference-anchored cut-off",
            color=C["orange"], fontsize=7.0, ha="center", va="bottom")
    for i, v in enumerate(counts):
        if v:
            ax.text(i, v + 1, str(v), ha="center", fontsize=7.0)
    ax.set_ylim(0, ymax * 1.30)

    # ── sulphathiazole ──
    ax = axes[1]
    counts = [stz.get(l, 0) for l in stz_levels]
    bars = ax.bar(range(len(stz_levels)), counts, width=0.72,
                  color=C["grey"], edgecolor="none")
    for i, l in enumerate(stz_levels):
        if l in ("128", "256", ">256"):
            bars[i].set_color(C["blue"])
        elif l == "64":
            bars[i].set_color(C["sky"])
    ax.set_xticks(range(len(stz_levels)))
    ax.set_xticklabels(stz_levels)
    ax.set_xlabel("Sulphathiazole MIC (µg ml$^{-1}$)")
    ax.set_ylabel("Number of isolates")
    ax.set_title("b", loc="left", fontweight="bold")
    ax.axvline(1.5, color=C["dark"], lw=1.0, ls="-")
    ax.text(1.60, max(counts) * 0.93, "primary cut-off\n(\u2265128)",
            fontsize=7.0, va="top", ha="left")
    ax.axvline(0.5, color=C["dark"], lw=0.8, ls=":")
    ax.text(0.60, max(counts) * 0.62, "sensitivity\n(\u226564)",
            fontsize=7.0, va="top", ha="left", color=C["dark"])
    for i, v in enumerate(counts):
        if v:
            ax.text(i, v + 1, str(v), ha="center", fontsize=7.0)
    ax.set_ylim(0, max(counts) * 1.18)

    fig.tight_layout()
    save(fig, outdir, "Fig2_MIC_distributions")


# ════════════════════════════════════════════════════════════════════════════
# Figure 3 — attainable significance
# ════════════════════════════════════════════════════════════════════════════
def figure3(outdir, n_strains=29, n_tests_core=2614, n_tests_pan=636):
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 65 * MM))

    # ── (a) combinatorial bound ──
    ax = axes[0]
    rs = np.arange(1, 15)
    floor = np.array([1.0 / math.comb(n_strains, int(r)) for r in rs])
    ax.semilogy(rs, floor, "o-", color=C["blue"], ms=3.5, lw=1.2,
                label=f"minimum attainable $p$ (n={n_strains})")

    thr_core = 0.05 / n_tests_core
    thr_pan = 0.05 / n_tests_pan
    ax.axhline(thr_core, color=C["red"], lw=1.0, ls="--",
               label=f"Bonferroni, core ({n_tests_core} patterns)")
    ax.axhline(thr_pan, color=C["orange"], lw=1.0, ls="--",
               label=f"Bonferroni, pan-genome ({n_tests_pan} patterns)")
    ax.axhline(0.05, color=C["grey"], lw=0.8, ls=":", label=r"$\alpha=0.05$")

    # Annotation anchors placed clear of the legend and threshold lines.
    for r, lab, col, xy, ha, va in (
            (4, "this study\n(\u2265128 µg ml$^{-1}$)", C["dark"], (-2, -50), "right", "top"),
            (7, "sensitivity\n(\u226564 µg ml$^{-1}$)", C["grey"], (8, 10), "left", "bottom")):
        f = 1.0 / math.comb(n_strains, r)
        ax.plot([r], [f], "o", ms=7, mfc="none", mec=col, mew=1.4)
        ax.annotate(lab, (r, f), textcoords="offset points", xytext=xy,
                    fontsize=7.0, color=col, ha=ha, va=va,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=col,
                                    shrinkA=0, shrinkB=6))

    ax.set_xlabel("Number of resistant strains")
    ax.set_ylabel("Minimum attainable $p$-value")
    ax.set_title("a", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper right", fontsize=7.0)
    ax.set_xticks(range(0, 15, 2))

    # ── (b) phylogenetic bound ──
    ax = axes[1]
    ks = np.arange(1, 11)
    # Scoary reports the two-sided binomial p-value; with k non-overlapping
    # pairs, all supporting the association, the smallest value is 2 x 0.5^k.
    pmin = np.minimum(1.0, 2 * 0.5 ** ks)
    ax.semilogy(ks, pmin, "s-", color=C["green"], ms=3.5, lw=1.2,
                label=r"minimum attainable $p$ = $2 \times 0.5^{k}$")
    ax.axhline(0.05, color=C["grey"], lw=0.8, ls=":", label=r"$\alpha=0.05$")

    for k, lab, col, xy in (
            (4, "sulphathiazole, $k$ = 4", C["blue"], (10, 10)),
            (6, "gentamicin, $k$ = 6", C["orange"], (12, 24))):
        ax.plot([k], [2 * 0.5 ** k], "o", ms=7, mfc="none", mec=col, mew=1.4)
        ax.annotate(lab, (k, 2 * 0.5 ** k), textcoords="offset points", xytext=xy,
                    fontsize=7.0, color=col, ha="left", va="bottom",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=col,
                                    shrinkA=0, shrinkB=6))

    ax.set_xlabel("Phylogenetically independent pairs available, $k$")
    ax.set_ylabel("Minimum attainable $p$-value")
    ax.set_title("b", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="lower left", fontsize=7.0)
    ax.set_xticks(range(1, 11))

    fig.tight_layout()
    save(fig, outdir, "Fig5_attainable_significance")


# ════════════════════════════════════════════════════════════════════════════
# Figure 4 — collapse under phylogenetic control
# ════════════════════════════════════════════════════════════════════════════
def figure4(outdir, gen_naive=2257, stz_naive=1335):
    fig, ax = plt.subplots(figsize=(SINGLE, 66 * MM))

    labels = ["Gentamicin", "Sulphathiazole"]
    naive = [gen_naive, stz_naive]
    corrected = [0, 0]

    x = np.arange(len(labels))
    w = 0.34

    ax.bar(x - w / 2, naive, w, color=C["orange"],
           label="Fisher + BH, no phylogenetic control")
    ax.bar(x + w / 2, corrected, w, color=C["blue"],
           label="pairwise comparisons over the tree")

    ymax = max(naive)
    for xi, v in zip(x - w / 2, naive):
        ax.text(xi, v + ymax * 0.02, f"{v:,}", ha="center", fontsize=7)
    for xi in x + w / 2:
        ax.text(xi, ymax * 0.02, "0", ha="center", fontsize=7,
                fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Significant variants")
    ax.set_ylim(0, ymax * 1.32)
    ax.legend(frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, 1.20), fontsize=7.0)
    fig.tight_layout()
    save(fig, outdir, "Fig4_naive_vs_corrected")


# ════════════════════════════════════════════════════════════════════════════
# Figure 2 — determinants by strain (requires the genotype-phenotype table)
# ════════════════════════════════════════════════════════════════════════════
def _mic_numeric(v):
    """Numeric value of a censored MIC string; combinations use the first part."""
    s = str(v).strip().lstrip("<>").split("/")[0]
    try:
        return float(s)
    except ValueError:
        return None


# Dilution series present on each plate, ascending, read from the Sensititre
# layouts (GN2F and AVIAN1F). These are properties of the plate, not of the
# isolates tested: normalising against the concentrations that happen to be
# observed in a collection distorts the scale for agents where most isolates
# fall below the lowest tested concentration.
PLATE_SERIES = {
    # GN2F
    "AMP":   [4, 8, 16, 32],
    "PIP":   [16, 32, 64, 128],
    "SXT":   [0.5, 1, 2, 4],          # trimethoprim component
    "GEN":   [2, 4, 8, 16],
    # AVIAN1F
    "AMOX":  [0.25, 0.5, 1, 2, 4, 8, 16],
    "SPE":   [8, 16, 32, 64],
    "SDM":   [32, 64, 128, 256],
    "STZ":   [32, 64, 128, 256],
    "STR":   [8, 16, 32, 64, 128, 256, 512, 1024],
    "TET":   [0.25, 0.5, 1, 2, 4, 8],
    "OXY":   [0.25, 0.5, 1, 2, 4, 8],
    "SXT.1": [0.5, 1, 2],             # trimethoprim component
    "GEN.1": [0.5, 1, 2, 4, 8],
    "ERY":   [0.12, 0.25, 0.5, 1, 2, 4],
    "CLI":   [0.5, 1, 2, 4],
    "TYLT":  [2.5, 5, 10, 20],
    "FFN":   [1, 2, 4, 8],
    "NEO":   [2, 4, 8, 16, 32],
    "AMI":   [8, 16, 32, 64],
    "TOB":   [4, 8],
}


def figure2(genopheno, outdir, mic_xlsx=None):
    """
    Determinant x strain, coloured by where the MIC sits within the dilution
    series tested for that agent on the plate.

    Earlier versions of this figure used the internal 0-5 index, which is not
    defined anywhere in the manuscript, and then three coarse categories, which
    grouped an amoxicillin MIC of 0.5 ug/ml with a tetracycline MIC of 4. The
    position within the plate's dilution series is self-explanatory, comparable
    across agents with different ranges, and separates determinants that have a
    matching phenotype from those that do not.
    """
    if not Path(genopheno).exists():
        print(f"  [SKIP] Figure 2: {genopheno} not found")
        return
    d = pd.read_csv(genopheno, sep="\t")
    d = d[(~d["intrinsic"]) & (d["testable"] == True)]
    d = d[~d["gene"].astype(str).str.lower().str.startswith(("gyra", "parc"))]
    if not len(d):
        print("  [SKIP] Figure 2: no testable acquired determinants")
        return

    RENAME = {
        "BT3_11028": "11028", "BT3_12": "12", "BT3_162": "162",
        "BT3_32": "32", "94-3-85": "94385", "960426_1/4C": "960426-1_4C",
        "C1P8190": "CIP8190", "V_4": "V4", "94-9-112": "94-8-112",
        "RIU1": "Riu1", "Ra_3": "Ra3",
    }
    d = d.copy()
    d["strain_pub"] = d["mic_strain"].replace(RENAME)
    d["strain_pub"] = d["strain_pub"].str.replace("CECT_", "CECT", regex=False)

    if not (mic_xlsx and Path(mic_xlsx).exists()):
        print("  [SKIP] Figure 2: raw MIC file required")
        return
    raw = read_mics(mic_xlsx)
    raw["STRAIN"] = raw["STRAIN"].astype(str).str.strip()
    raw = raw.set_index("STRAIN")

    def position(mic_name, agent):
        """0 = lowest concentration on the plate or below, 1 = highest or above."""
        if mic_name not in raw.index or agent not in raw.columns:
            return None
        s = PLATE_SERIES.get(agent)
        if not s or len(s) < 2:
            return None
        v = str(raw.loc[mic_name, agent]).strip()
        if v in ("nan", ""):
            return None
        if v.startswith("<"):
            return 0.0
        if v.startswith(">"):
            return 1.0
        n = _mic_numeric(v)
        if n is None:
            return None
        # index within the plate series, mapped onto [0, 1]
        nearest = min(range(len(s)), key=lambda i: abs(s[i] - n))
        return nearest / (len(s) - 1)

    genes = sorted(d["gene"].unique())
    strains = sorted(d["strain_pub"].unique())

    val = {}
    for (strain, gene), sub in d.groupby(["strain_pub", "gene"]):
        ps = [position(r["mic_strain"], r["antibiotic"]) for _, r in sub.iterrows()]
        ps = [p for p in ps if p is not None]
        val[(strain, gene)] = max(ps) if ps else None

    cmap = plt.get_cmap("YlOrBr")
    fig, ax = plt.subplots(figsize=(DOUBLE, 0.32 * len(strains) + 1.7))

    for i, s in enumerate(strains):
        for j, g in enumerate(genes):
            if (s, g) not in val:
                continue
            p = val[(s, g)]
            face = "#DDDDDD" if p is None else cmap(0.12 + 0.85 * p)
            ax.add_patch(mpatches.Rectangle(
                (j - 0.42, i - 0.40), 0.84, 0.80,
                facecolor=face, edgecolor="#999999", lw=0.4))

    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(strains) - 0.5, -0.5)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right", style="italic")
    ax.set_yticks(range(len(strains)))
    ax.set_yticklabels(strains)
    ax.tick_params(length=0)
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.12, 0.97))
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.015)
    cb.set_ticks([0.12, 0.97])
    cb.set_ticklabels(["lowest concentration\ntested, or below",
                       "highest concentration\ntested, or above"])
    cb.ax.tick_params(length=0)
    cb.set_label("Position of the MIC within the tested dilution series")
    cb.outline.set_visible(False)

    fig.tight_layout()
    save(fig, outdir, "Fig3_determinants")


# ════════════════════════════════════════════════════════════════════════════
# Figure 5 — the contig_2 variant against the phylogeny
# ════════════════════════════════════════════════════════════════════════════
def figure5(tree_path, matrix_path, outdir):
    if not (Path(tree_path).exists() and Path(matrix_path).exists()):
        print("  [SKIP] Figure 5: tree or matrix not found")
        return
    try:
        from Bio import Phylo
    except ImportError:
        print("  [SKIP] Figure 5: Biopython not available")
        return

    VARIANT = "contig_2:393252:C>T"
    STZ = {  # raw MIC, censored series
        "CECT898": ">256", "CECT4606": ">256", "CECT4608": ">256",
        "YJ016": ">256", "11028": "256", "CG100": "64", "32": "64",
        "162": "64",
    }
    LINEAGE = {
        # L1
        "CECT4608": "L1", "CG100": "L1", "V4": "L1", "YJ016": "L1",
        # L2
        "94-9-112": "L2", "94-8-112": "L2", "94-9-119": "L2",
        "95-8-161": "L2", "95-8-6": "L2", "95-8-7": "L2",
        "960426-1_4C": "L2", "A14": "L2", "CECT4604": "L2",
        "CECT4865": "L2", "CECT4866": "L2", "CECT4999": "L2",
        "CECT5763": "L2", "CECT5769": "L2", "CECT7030": "L2",
        "CECT898": "L2", "CIP8190": "L2", "PD-2-51": "L2", "Ra_3": "L2",
        # L3
        "11028": "L3", "12": "L3", "32": "L3", "162": "L3",
        # L4
        "94385": "L4", "RIU1": "L4",
        # not in Roig et al. (2018); assigned by F. J. Roig Molina (Table S1, note c)
        "CECT4606": "L1*", "CECT5198": "L2*",
    }
    # Tip names in the tree follow the analysis files; display the published
    # strain designations (Table S1).
    DISPLAY = {"94-9-112": "94-8-112", "Ra_3": "Ra3", "RIU1": "Riu1"}

    m = pd.read_csv(matrix_path, sep="\t", index_col=0)
    if VARIANT not in m.index:
        print(f"  [SKIP] Figure 5: {VARIANT} not in matrix")
        return
    carrier = {c for c in m.columns if m.loc[VARIANT, c] == 1}

    t = Phylo.read(tree_path, "newick")
    t.ladderize()
    tips = [x.name for x in t.get_terminals()]

    # Support labels are omitted: all nodes are at 100 except seven shallow
    # nodes (81-96), whose labels overlap on the very short branches. The
    # values are reported in the legend instead.
    for cl in t.get_nonterminals():
        cl.confidence = None

    fig, (axT, axD) = plt.subplots(
        1, 2, figsize=(DOUBLE, 0.22 * len(tips) + 1.2),
        gridspec_kw={"width_ratios": [2.2, 1]})

    Phylo.draw(t, axes=axT, do_show=False,
               label_func=lambda c: DISPLAY.get(c.name, c.name) if c.is_terminal() else "")
    axT.set_xlabel("Substitutions per variable site")
    axT.set_ylabel("")
    axT.set_title("a", loc="left", fontweight="bold")
    for sp in ("top", "right", "left"):
        axT.spines[sp].set_visible(False)
    axT.set_yticks([])

    order = tips[::-1]
    y = np.arange(len(order))
    axD.set_ylim(-0.5, len(order) - 0.5)
    axD.set_yticks(y)
    axD.set_yticklabels([f"{DISPLAY.get(s, s)}  ({LINEAGE.get(s,'--')})" for s in order],
                        fontsize=7.0)
    axD.set_xticks([0, 1])
    axD.set_xticklabels(["variant", "STZ MIC"], fontsize=7)
    axD.set_xlim(-0.5, 1.5)
    axD.set_title("b", loc="left", fontweight="bold")
    for sp in ("top", "right", "left", "bottom"):
        axD.spines[sp].set_visible(False)
    axD.tick_params(length=0)

    for i, s in enumerate(order):
        if s in carrier:
            axD.add_patch(mpatches.Rectangle((-0.28, i - 0.34), 0.56, 0.68,
                                             color=C["blue"]))
        elif s in m.columns:
            axD.add_patch(mpatches.Rectangle((-0.28, i - 0.34), 0.56, 0.68,
                                             facecolor="white",
                                             edgecolor=C["grey"], lw=0.5))
        mic = STZ.get(s)
        if mic:
            col = C["red"] if mic in (">256", "256") else C["yellow"]
            axD.add_patch(mpatches.Rectangle((0.72, i - 0.34), 0.56, 0.68,
                                             color=col))
            axD.text(1.0, i, mic, ha="center", va="center", fontsize=7.0)

    handles = [
        mpatches.Patch(color=C["blue"], label="variant present"),
        mpatches.Patch(facecolor="white", edgecolor=C["grey"],
                       label="variant absent"),
        mpatches.Patch(color=C["red"], label="MIC \u2265128 µg ml$^{-1}$"),
        mpatches.Patch(color=C["yellow"], label="MIC 64 µg ml$^{-1}$"),
    ]
    axD.legend(handles=handles, frameon=False, fontsize=7.0,
               loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2)

    fig.tight_layout()
    save(fig, outdir, "Fig6_lineage_marker")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mic-xlsx", default="data/mic/MICs_all_Vv.xlsx")
    ap.add_argument("--genopheno",
                    default="results/genotype_phenotype/genotype_phenotype_check.tsv")
    ap.add_argument("--tree", default="results/pyseer/core.treefile")
    ap.add_argument("--matrix", default="results/summary/variants_matrix.tsv")
    ap.add_argument("--outdir", default="figures")
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    print("Generating manuscript figures")

    figure1(a.mic_xlsx, outdir)
    figure2(a.genopheno, outdir, mic_xlsx=a.mic_xlsx)
    figure3(outdir)
    figure4(outdir)
    figure5(a.tree, a.matrix, outdir)

    print(f"\nFigures written to {outdir}/")
    print("PDF for submission, PNG for drafting.")


if __name__ == "__main__":
    main()

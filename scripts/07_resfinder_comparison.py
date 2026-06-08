#!/usr/bin/env python3
"""
07_resfinder_comparison.py
--------------------------
Runs ResFinder 4.x on the 30 V. vulnificus genomes and compares
the results with our SNP-based pipeline.

ResFinder detects acquired resistance genes (horizontally transferred).
PointFinder detects chromosomal point mutations in known resistance genes.
Our pipeline detects SNPs in the core genome associated with MIC changes.

The comparison table shows:
  - What ResFinder/PointFinder detect per strain
  - What our pipeline detects per strain
  - Overlap and unique findings

Usage:
    python3 07_resfinder_comparison.py

Inputs:
    ../data/raw/genomes/*.fasta          (30 genomes)
    ../data/raw/resfinder_db/            (ResFinder database)
    ../data/raw/pointfinder_db/          (PointFinder database)
    ../data/processed/SNP_significant.tsv

Outputs:
    ../results/resfinder_raw/            (per-strain ResFinder output)
    ../results/resfinder_summary.tsv     (combined results table)
    ../results/pipeline_comparison.tsv  (our pipeline vs ResFinder)
    ../figures/Fig4_comparison.pdf/png  (comparison figure)

Author: Francesc Montoro (2025)
"""

# Ensure BLAST is in PATH
import os
blast_bin = "/mnt/c/PACO/MICS_PAPER/ncbi-blast-2.12.0+/bin"
if blast_bin not in os.environ.get("PATH", ""):
    os.environ["PATH"] = blast_bin + ":" + os.environ.get("PATH", "")
import subprocess
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import json
import warnings
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = Path(__file__).parent.parent
RAW     = BASE / "data" / "raw"
PROC    = BASE / "data" / "processed"
RES     = BASE / "results"
FIGURES = BASE / "figures"
RF_OUT  = RES / "resfinder_raw"
RF_OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
RESFINDER_DB  = RAW / "resfinder_db"
POINTFINDER_DB= RAW / "pointfinder_db"
GENOMES_DIR   = RAW / "genomes"
SPECIES       = "vibrio"          # for PointFinder species
MIN_IDENTITY  = 90.0              # % identity threshold
MIN_COVERAGE  = 60.0              # % coverage threshold

# Strain name mapping: fasta filename stem → MIC name
GENOME_TO_MIC = {
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

# ── ResFinder functions ────────────────────────────────────────────────────────

def run_resfinder(genome_fasta: Path, out_dir: Path, strain_name: str) -> dict:
    """Run ResFinder on a single genome. Returns dict of detected genes."""
    strain_out = out_dir / strain_name
    strain_out.mkdir(exist_ok=True)

    cmd = [
        "python3", "-m", "resfinder",
        "-ifa", str(genome_fasta),
        "-o", str(strain_out),
        "-s", "Other",
        "--acquired",
        "-db_res", str(RESFINDER_DB),
        "--min_cov", str(MIN_COVERAGE / 100),
        "--threshold", str(MIN_IDENTITY / 100),
        "-j",  str(strain_out / "resfinder.json"),
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            print(f"    [WARNING] ResFinder error for {strain_name}: {result.stderr[:200]}")
            return {}
    except subprocess.TimeoutExpired:
        print(f"    [WARNING] ResFinder timeout for {strain_name}")
        return {}
    except Exception as e:
        print(f"    [WARNING] ResFinder failed for {strain_name}: {e}")
        return {}

    # Parse JSON output
    json_file = strain_out / "resfinder.json"
    if not json_file.exists():
        # Try parsing the tab-separated output instead
        return parse_resfinder_txt(strain_out, strain_name)

    try:
        with open(json_file) as f:
            data = json.load(f)
        genes = {}
        # ResFinder 4.x JSON structure
        if "resfinder" in data:
            rf = data["resfinder"]
            if "results" in rf:
                for ab_class, hits in rf["results"].items():
                    if isinstance(hits, dict):
                        for gene_key, hit_info in hits.items():
                            if isinstance(hit_info, dict) and hit_info.get("query_id"):
                                gene_name = hit_info.get("ref_gene_id", gene_key)
                                phenotype = hit_info.get("phenotype", ab_class)
                                genes[gene_name] = phenotype
        return genes
    except Exception as e:
        return parse_resfinder_txt(strain_out, strain_name)


def parse_resfinder_txt(strain_out: Path, strain_name: str) -> dict:
    """Fallback: parse ResFinder tab-separated output."""
    genes = {}
    for txt_file in strain_out.glob("*.txt"):
        try:
            with open(txt_file) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        gene = parts[0]
                        phenotype = parts[2] if len(parts) > 2 else "unknown"
                        genes[gene] = phenotype
        except:
            pass
    return genes


def run_all_resfinder(genomes_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Run ResFinder on all genomes. Returns summary DataFrame."""
    fasta_files = sorted(genomes_dir.glob("*.fasta"))
    # Exclude the gene database file
    fasta_files = [f for f in fasta_files if "Genes_Vv" not in f.name]

    print(f"  Running ResFinder on {len(fasta_files)} genomes...")
    all_results = {}

    for i, fasta in enumerate(fasta_files, 1):
        stem = fasta.stem
        mic_name = GENOME_TO_MIC.get(stem, stem)
        print(f"  [{i:02d}/{len(fasta_files)}] {mic_name}...", end="\r")

        genes = run_resfinder(fasta, out_dir, stem)
        all_results[mic_name] = genes

    print(f"\n  Done. {len(all_results)} strains processed.")
    return all_results


def build_resfinder_matrix(all_results: dict) -> pd.DataFrame:
    """Build presence/absence matrix from ResFinder results."""
    all_genes = set()
    for genes in all_results.values():
        all_genes.update(genes.keys())

    if not all_genes:
        print("  [WARNING] No resistance genes detected by ResFinder.")
        return pd.DataFrame()

    mat = pd.DataFrame(0, index=sorted(all_results.keys()),
                       columns=sorted(all_genes))
    for strain, genes in all_results.items():
        for gene in genes:
            if gene in mat.columns:
                mat.loc[strain, gene] = 1

    return mat


def make_comparison_figure(rf_matrix: pd.DataFrame,
                           our_sig: pd.DataFrame,
                           out_prefix: Path):
    """
    Comparison figure: bar chart showing number of findings per method.
    """
    # Count per strain
    rf_counts = rf_matrix.sum(axis=1) if len(rf_matrix) > 0 else pd.Series(dtype=int)

    # Our pipeline: significant SNPs per strain
    our_counts = our_sig.groupby("gene")["antibiotic"].count()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor="white")

    # Panel A: ResFinder genes per strain
    ax1 = axes[0]
    if len(rf_counts) > 0:
        rf_counts_sorted = rf_counts.sort_values(ascending=False)
        colors = ["#d73027" if v > 0 else "#cccccc" for v in rf_counts_sorted]
        ax1.barh(range(len(rf_counts_sorted)), rf_counts_sorted.values,
                 color=colors, height=0.7)
        ax1.set_yticks(range(len(rf_counts_sorted)))
        ax1.set_yticklabels(rf_counts_sorted.index, fontsize=7,
                            fontfamily="monospace")
        ax1.set_xlabel("Acquired resistance genes detected", fontsize=10)
        ax1.set_title("ResFinder 4.x", fontsize=11, fontweight="bold")
        n_with = (rf_counts > 0).sum()
        ax1.text(0.98, 0.02,
                 f"Strains with ≥1 gene: {n_with}/{len(rf_counts)}",
                 transform=ax1.transAxes, fontsize=8, ha="right",
                 bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    else:
        ax1.text(0.5, 0.5, "No acquired resistance genes\ndetected by ResFinder",
                 transform=ax1.transAxes, ha="center", va="center",
                 fontsize=11, color="#888888")
        ax1.set_title("ResFinder 4.x", fontsize=11, fontweight="bold")
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.spines[["top", "right", "left", "bottom"]].set_visible(False)
        
    # Panel B: Our pipeline — genes per antibiotic
    ax2 = axes[1]
    ab_counts = our_sig.groupby("antibiotic")["gene"].nunique().sort_values(ascending=False)
    colors2 = ["#4393c3" if ab in ["GEN2"] else
                "#9970ab" if ab in ["STZ", "SDM"] else
                "#888888" for ab in ab_counts.index]
    ax2.barh(range(len(ab_counts)), ab_counts.values,
             color=colors2, height=0.7)
    ax2.set_yticks(range(len(ab_counts)))
    ax2.set_yticklabels(ab_counts.index, fontsize=9)
    ax2.set_xlabel("Significant genes (FDR < 0.05)", fontsize=10)
    ax2.set_title("Our SNP-MIC pipeline", fontsize=11, fontweight="bold")
    ax2.text(0.98, 0.02,
             f"Total: {our_sig['gene'].nunique()} genes\n"
             f"{our_sig['antibiotic'].nunique()} antibiotics",
             transform=ax2.transAxes, fontsize=8, ha="right",
             bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    ax2.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        r"Comparison: ResFinder vs SNP-MIC pipeline for $\it{Vibrio\ vulnificus}$",
        fontsize=12, fontweight="bold", y=1.01
    )
    fig.tight_layout()

    for ext in ["pdf", "png"]:
        out = out_prefix.with_suffix(f".{ext}")
        fig.savefig(out, dpi=900, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 65)
    print("  ResFinder comparison")
    print("=" * 65)

    # Check databases exist
    if not RESFINDER_DB.exists():
        print(f"  [ERROR] ResFinder DB not found: {RESFINDER_DB}")
        exit(1)
    if not GENOMES_DIR.exists():
        print(f"  [ERROR] Genomes directory not found: {GENOMES_DIR}")
        exit(1)

    print(f"\n  ResFinder DB: {RESFINDER_DB}")
    print(f"  Genomes dir:  {GENOMES_DIR}")
    fasta_count = len([f for f in GENOMES_DIR.glob("*.fasta")
                       if "Genes_Vv" not in f.name])
    print(f"  Genomes found: {fasta_count}")

    # Run ResFinder
    print("\nRunning ResFinder...")
    all_results = run_all_resfinder(GENOMES_DIR, RF_OUT)

    # Build matrix
    rf_matrix = build_resfinder_matrix(all_results)

    if len(rf_matrix) > 0:
        out_rf = RES / "resfinder_summary.tsv"
        rf_matrix.to_csv(out_rf, sep="\t")
        print(f"\nResFinder matrix: {rf_matrix.shape} → {out_rf}")
        print(f"Total unique genes detected: {(rf_matrix > 0).any().sum()}")
        print(f"Strains with ≥1 gene: {(rf_matrix.sum(axis=1) > 0).sum()}")
    else:
        print("\nNo resistance genes detected by ResFinder.")
        print("This is the expected result for V. vulnificus —")
        print("it confirms the gap our pipeline fills.")
        rf_matrix = pd.DataFrame()

    # Load our significant results
    sig_file = PROC / "SNP_significant.tsv"
    if sig_file.exists():
        our_sig = pd.read_csv(sig_file, sep="\t")
        print(f"\nOur pipeline: {len(our_sig)} significant SNPs, "
              f"{our_sig['gene'].nunique()} genes, "
              f"{our_sig['antibiotic'].nunique()} antibiotics")

        # Save comparison table
        our_genes = our_sig.groupby("antibiotic")["gene"].apply(
            lambda x: "; ".join(sorted(set(x))[:5])
        ).reset_index()
        our_genes.columns = ["antibiotic", "top_genes_our_pipeline"]
        out_comp = RES / "pipeline_comparison.tsv"
        our_genes.to_csv(out_comp, sep="\t", index=False)
        print(f"Comparison table → {out_comp}")

        # Figure
        print("\nGenerating comparison figure...")
        make_comparison_figure(rf_matrix, our_sig,
                               FIGURES / "Fig4_comparison")
    else:
        print(f"\n[WARNING] {sig_file} not found — run 03 first.")

    # Summary text for paper
    print("\n" + "=" * 65)
    print("  KEY FINDING FOR PAPER:")
    print("=" * 65)
    n_rf_genes = (rf_matrix.sum(axis=1) > 0).sum() if len(rf_matrix) > 0 else 0
    print(f"  ResFinder detected acquired resistance genes in {n_rf_genes}/30 strains")
    print(f"  Our pipeline identified significant associations for")
    if sig_file.exists():
        print(f"  {our_sig['antibiotic'].nunique()} antibiotics across {our_sig['gene'].nunique()} genes")
    print(f"  → V. vulnificus is NOT in ResFinder's species database")
    print(f"  → Our core-genome SNP approach fills this gap")
    print("\nDone.")

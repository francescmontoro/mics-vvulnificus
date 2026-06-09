#!/usr/bin/env python3
"""
run_cleancoreSNP_v2.py
-----------------------
Batch processor: runs cleancoreSNP_v2.py on all .aln files in a directory.
Generates a global summary TSV with statistics per gene.

Usage:
    python3 run_cleancoreSNP_v2.py

Outputs:
    ../data/processed/aln_cleaning_summary.tsv   (per-gene statistics)
    ../data/raw/LOGS/cleancoreSNP_v2.log         (detailed per-gene log)
    ../data/raw/SNP_v2/*.aln                     (cleaned alignments)

Author: Francesc Montoro (2025)
"""

import sys
import os
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# Add scripts dir to path to import cleancoreSNP_v2
sys.path.insert(0, str(Path(__file__).parent))
from cleancoreSNP_v2 import parse_fasta, clean_alignment, write_cleaned_fasta, write_log

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = Path(__file__).parent.parent
ALN_DIR = Path("/mnt/c/PACO/MICS_PAPER/MICs/PIPELINE/TSV/GENES/CORE/ALN")
OUT_DIR = BASE / "data" / "raw" / "SNP_v2_clean_aln"
LOG_DIR = BASE / "data" / "raw" / "LOGS"
PROC    = BASE / "data" / "processed"

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)
PROC.mkdir(parents=True, exist_ok=True)

REFERENCE = "CECT5769"   # reference strain name in alignment files


def main():
    print("=" * 70)
    print("  Batch alignment cleaning — cleancoreSNP_v2")
    print(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    aln_files = sorted(ALN_DIR.glob("*.aln"))
    n_total   = len(aln_files)
    print(f"\n  Alignment files found: {n_total}")
    print(f"  Output dir:  {OUT_DIR}")
    print(f"  Log dir:     {LOG_DIR}\n")

    # Global stats
    summary_rows = []
    n_clean      = 0
    n_with_iupac = 0
    n_with_gaps  = 0
    n_errors     = 0

    log_path = str(LOG_DIR / "cleancoreSNP_v2.log")

    # Write log header
    with open(log_path, "w") as f:
        f.write(f"cleancoreSNP_v2 batch run — {datetime.now()}\n")
        f.write(f"Reference: {REFERENCE}\n")
        f.write(f"Input dir: {ALN_DIR}\n")
        f.write(f"Output dir: {OUT_DIR}\n")
        f.write("=" * 70 + "\n")

    for i, aln_file in enumerate(aln_files, 1):
        gene_name = aln_file.name.split(".fasta")[0]

        # Progress
        if i % 100 == 0 or i == n_total:
            print(f"  [{i:4d}/{n_total}] Processing...", end="\r")

        # Parse
        try:
            seqs = parse_fasta(str(aln_file))
        except Exception as e:
            print(f"\n  [ERROR] {gene_name}: {e}")
            n_errors += 1
            continue

        if not seqs:
            n_errors += 1
            continue

        # Clean
        out_path = str(OUT_DIR / (gene_name + "_clean.aln"))
        cleaned, stats = clean_alignment(seqs, REFERENCE, gene_name)

        # Write cleaned alignment
        write_cleaned_fasta(cleaned, out_path)

        # Write to log (only if has substitutions)
        total_subst = stats["iupac_total"] + stats["gap_total"] + stats["n_total"]
        if total_subst > 0:
            write_log(stats, log_path)

        # Update global counters
        if stats["iupac_total"] > 0:
            n_with_iupac += 1
        if stats["gap_total"] > 0:
            n_with_gaps += 1
        if total_subst == 0:
            n_clean += 1

        # Summary row
        summary_rows.append({
            "gene":          gene_name,
            "n_strains":     stats["n_strains"],
            "seq_length":    stats["seq_length"],
            "iupac_subst":   stats["iupac_total"],
            "gap_subst":     stats["gap_total"],
            "n_preexisting": stats["n_total"],
            "total_subst":   total_subst,
            "pct_subst":     round(total_subst / (stats["seq_length"] * stats["n_strains"]) * 100, 4)
                             if stats["seq_length"] > 0 else 0,
            "has_error":     "warning" in stats or "error" in stats,
        })

    print(f"\n  Done. {n_total} files processed.\n")

    # Write global summary
    import pandas as pd
    summary_df = pd.DataFrame(summary_rows)
    out_summary = PROC / "aln_cleaning_summary.tsv"
    summary_df.to_csv(out_summary, sep="\t", index=False)

    # Print summary
    print("=" * 70)
    print("  GLOBAL SUMMARY")
    print("=" * 70)
    print(f"  Total alignments processed: {n_total}")
    print(f"  Clean (no substitutions):   {n_clean} ({n_clean/n_total*100:.1f}%)")
    print(f"  With IUPAC ambiguities:     {n_with_iupac} ({n_with_iupac/n_total*100:.1f}%)")
    print(f"  With gaps:                  {n_with_gaps} ({n_with_gaps/n_total*100:.1f}%)")
    print(f"  Errors:                     {n_errors}")
    print(f"\n  Total IUPAC substitutions:  {summary_df['iupac_subst'].sum():,}")
    print(f"  Total gap substitutions:    {summary_df['gap_subst'].sum():,}")
    print(f"  Total → N:                  {summary_df['total_subst'].sum():,}")
    print(f"\n  Summary TSV → {out_summary}")
    print(f"  Detailed log → {log_path}")

    # Top genes with most substitutions
    top = summary_df[summary_df["total_subst"] > 0].nlargest(10, "total_subst")
    if len(top) > 0:
        print(f"\n  Top 10 genes with most substitutions:")
        for _, row in top.iterrows():
            print(f"    {row['gene'][:50]:<50} "
                  f"total={row['total_subst']:5d}  "
                  f"gaps={row['gap_subst']:5d}  "
                  f"IUPAC={row['iupac_subst']:4d}")

    print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()

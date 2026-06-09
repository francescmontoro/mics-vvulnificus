#!/usr/bin/env python3
"""
cleancoreSNP_v2.py
------------------
Improved version of cleancoreSNP.py for cleaning MSA alignment files.

Changes from v1:
  - IUPAC ambiguous bases → N (instead of most frequent base)
  - Gaps (-) → N (instead of reference nucleotide)
  - N positions are ignored by checkcoreSNPs.py (not in codon table)
  - Full logging of all substitutions made per gene/strain
  - Summary statistics written to log file

Rationale:
  Substituting ambiguous bases or gaps with the most frequent nucleotide
  or the reference nucleotide can introduce artificial SNPs. Using N
  is more conservative: checkcoreSNPs.py skips any codon containing N,
  so these positions do not contribute to the SNP database.

Usage:
    python3 cleancoreSNP_v2.py <aln_file> <reference_strain>

    Or in batch (from the ALN directory):
    for f in *.aln; do
        python3 cleancoreSNP_v2.py "$f" CECT5769
    done

Inputs:
    <aln_file>         : FASTA-format MSA alignment file
    <reference_strain> : name of the reference strain (e.g. CECT5769)

Outputs:
    SNP/<gene>_clean.aln      : cleaned alignment (IUPAC/gaps → N)
    LOGS/cleancoreSNP_v2.log  : per-run log appended each call

Author: Francesc Montoro (refactored 2025)
"""

import sys
import os
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# ── IUPAC ambiguous bases (excluding N which we use as our sentinel) ──────────
IUPAC_AMBIG = set("WSMKRYBDHVZ")  # N is handled separately, - is gap

SENTINEL = "N"   # replacement for both ambiguous bases and gaps


def parse_fasta(filepath: str) -> dict:
    """Read a FASTA file and return {strain: sequence} dict."""
    seqs = {}
    current = None
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].strip()
                seqs[current] = []
            elif current is not None:
                seqs[current].extend(line.upper())
    return {k: v for k, v in seqs.items()}


def clean_alignment(seqs: dict, reference: str, gene_name: str) -> tuple:
    """
    Clean alignment sequences:
      - IUPAC ambiguous bases → N
      - Gaps (-) → N
      - Existing N → kept as N

    Returns:
        cleaned: {strain: cleaned_sequence_string}
        stats:   summary dict with counts of substitutions
    """
    stats = {
        "gene":          gene_name,
        "n_strains":     len(seqs),
        "seq_length":    0,
        "iupac_total":   0,
        "gap_total":     0,
        "n_total":       0,  # pre-existing N
        "per_strain":    defaultdict(lambda: {"iupac": 0, "gap": 0, "n": 0}),
    }

    if reference not in seqs:
        stats["error"] = f"Reference '{reference}' not found"
        # Fall back: use first strain as reference
        reference = list(seqs.keys())[0]
        stats["warning"] = f"Using '{reference}' as fallback reference"

    ref_seq = seqs[reference]
    stats["seq_length"] = len(ref_seq)

    cleaned = {}
    for strain, seq in seqs.items():
        clean_seq = []
        for i, base in enumerate(seq):
            if base in IUPAC_AMBIG:
                clean_seq.append(SENTINEL)
                stats["iupac_total"] += 1
                stats["per_strain"][strain]["iupac"] += 1
            elif base == "-":
                clean_seq.append(SENTINEL)
                stats["gap_total"] += 1
                stats["per_strain"][strain]["gap"] += 1
            elif base == "N":
                clean_seq.append(SENTINEL)
                stats["n_total"] += 1
                stats["per_strain"][strain]["n"] += 1
            else:
                clean_seq.append(base)
        cleaned[strain] = "".join(clean_seq)

    return cleaned, stats


def write_cleaned_fasta(cleaned: dict, output_path: str):
    """Write cleaned sequences to FASTA file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        for strain, seq in cleaned.items():
            f.write(f">{strain}\n{seq}\n")


def write_log(stats: dict, log_path: str):
    """Append run statistics to log file."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "a") as f:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"\n{'='*70}\n")
        f.write(f"[{ts}] Gene: {stats['gene']}\n")
        f.write(f"  Strains: {stats['n_strains']}  |  "
                f"Length: {stats['seq_length']} bp\n")

        if "error" in stats:
            f.write(f"  ERROR: {stats['error']}\n")
        if "warning" in stats:
            f.write(f"  WARNING: {stats['warning']}\n")

        f.write(f"  Total substitutions → N:\n")
        f.write(f"    IUPAC ambiguous: {stats['iupac_total']}\n")
        f.write(f"    Gaps (-):        {stats['gap_total']}\n")
        f.write(f"    Pre-existing N:  {stats['n_total']}\n")
        f.write(f"    TOTAL:           "
                f"{stats['iupac_total'] + stats['gap_total'] + stats['n_total']}\n")

        # Per-strain detail (only if any substitutions)
        any_subst = any(
            v["iupac"] + v["gap"] + v["n"] > 0
            for v in stats["per_strain"].values()
        )
        if any_subst:
            f.write(f"  Per-strain detail:\n")
            for strain, counts in sorted(stats["per_strain"].items()):
                total = counts["iupac"] + counts["gap"] + counts["n"]
                if total > 0:
                    f.write(f"    {strain:<25} "
                            f"IUPAC={counts['iupac']:4d}  "
                            f"gaps={counts['gap']:4d}  "
                            f"N={counts['n']:4d}\n")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    aln_file  = sys.argv[1]
    reference = sys.argv[2]

    if not os.path.exists(aln_file):
        print(f"ERROR: File not found: {aln_file}")
        sys.exit(1)

    # Derive gene name from filename
    gene_name = Path(aln_file).name.split(".fasta")[0]

    # Parse alignment
    seqs = parse_fasta(aln_file)
    if not seqs:
        print(f"ERROR: No sequences found in {aln_file}")
        sys.exit(1)

    # Clean
    cleaned, stats = clean_alignment(seqs, reference, gene_name)

    # Output paths
    out_dir  = os.path.join(os.path.dirname(aln_file), "..", "..", "SNP_v2")
    log_dir  = os.path.join(os.path.dirname(aln_file), "..", "..", "LOGS")
    out_path = os.path.join(out_dir, gene_name + "_clean.aln")
    log_path = os.path.join(log_dir, "cleancoreSNP_v2.log")

    # Write outputs
    write_cleaned_fasta(cleaned, out_path)
    write_log(stats, log_path)

    # Console summary
    total_subst = stats["iupac_total"] + stats["gap_total"] + stats["n_total"]
    if total_subst > 0:
        print(f"{gene_name}: {total_subst} positions → N "
              f"(IUPAC={stats['iupac_total']}, "
              f"gaps={stats['gap_total']}, "
              f"N={stats['n_total']})")
    else:
        print(f"{gene_name}: clean")


if __name__ == "__main__":
    main()

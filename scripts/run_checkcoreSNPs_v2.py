#!/usr/bin/env python3
"""
run_checkcoreSNPs_v2.py
-----------------------
Batch runner for checkcoreSNPs.py on the new cleaned alignments
(produced by cleancoreSNP_v2.py with IUPAC→N and gaps→N treatment).

Collects all non-synonymous SNPs into a single allnsSNP_v2.tsv file.

Usage:
    python3 run_checkcoreSNPs_v2.py

Inputs:
    ../data/raw/SNP_v2_clean_aln/*_clean.aln   (cleaned alignments)

Outputs:
    ../data/raw/allnsSNP_v2.tsv               (all non-synonymous SNPs)
    ../data/raw/allsSNP_v2.tsv                (all synonymous SNPs)
    ../data/processed/snp_summary_v2.tsv      (per-gene SNP counts)

Author: Francesc Montoro (2025)
"""

import sys
import os
from pathlib import Path
from datetime import datetime
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE      = Path(__file__).parent.parent
CLEAN_DIR = BASE / "data" / "raw" / "SNP_v2_clean_aln"
PROC      = BASE / "data" / "processed"
RAW       = BASE / "data" / "raw"
PROC.mkdir(exist_ok=True)

REFERENCE = "CECT5769"

# ── Codon table ───────────────────────────────────────────────────────────────
CODON_TABLE = {
    'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile', 'ATG': 'Met',
    'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',
    'AAC': 'Asn', 'AAT': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGC': 'Ser', 'AGT': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',
    'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',
    'CAC': 'His', 'CAT': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',
    'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',
    'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',
    'GAC': 'Asp', 'GAT': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly',
    'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',
    'TTC': 'Phe', 'TTT': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
    'TAC': 'Tyr', 'TAT': 'Tyr', 'TAA': 'STOP', 'TAG': 'STOP',
    'TGC': 'Cys', 'TGT': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp',
}


def parse_fasta(filepath: Path) -> dict:
    seqs = {}
    current = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].strip()
                seqs[current] = ""
            elif current:
                seqs[current] += line.upper()
    return seqs


def find_snps(seqs: dict, reference: str, gene_name: str) -> tuple:
    """
    Find synonymous and non-synonymous SNPs by codon comparison.
    Codons containing N are skipped (from gaps/IUPAC treatment).
    Returns (nonsynonymous_list, synonymous_list).
    """
    if reference not in seqs:
        return [], []

    ref_seq = seqs[reference]
    nonsyn = []
    syn    = []

    for strain, seq in seqs.items():
        if strain == reference:
            continue
        if len(seq) != len(ref_seq):
            continue

        i = 0
        while i < len(seq) - 2:
            codon_ref   = ref_seq[i:i+3]
            codon_strain = seq[i:i+3]

            # Skip if either codon contains N (gap/ambiguity sentinel)
            if "N" in codon_ref or "N" in codon_strain:
                i += 3
                continue

            if codon_ref in CODON_TABLE and codon_strain in CODON_TABLE:
                if codon_ref != codon_strain:
                    aa_ref    = CODON_TABLE[codon_ref]
                    aa_strain = CODON_TABLE[codon_strain]
                    pos = i + 1  # 1-based

                    record = {
                        "Gene":      gene_name,
                        "Strain":    strain,
                        "Position":  pos,
                        "Reference": codon_ref,
                        "Aa_ref":    aa_ref,
                        "SNP":       codon_strain,
                        "Aa_mod":    aa_strain,
                    }

                    if aa_ref != aa_strain:
                        record["Type_SNP"] = "Nonsynonymous"
                        nonsyn.append(record)
                    else:
                        record["Type_SNP"] = "Synonymous"
                        syn.append(record)
            i += 3

    return nonsyn, syn


def main():
    print("=" * 70)
    print("  Batch SNP detection — checkcoreSNPs v2")
    print(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    clean_files = sorted(CLEAN_DIR.glob("*_clean.aln"))
    n_total = len(clean_files)
    print(f"\n  Clean alignment files: {n_total}")
    print(f"  Reference: {REFERENCE}\n")

    all_nonsyn = []
    all_syn    = []
    summary    = []

    for i, aln_file in enumerate(clean_files, 1):
        gene_name = aln_file.name.replace("_clean.aln", "")

        if i % 200 == 0 or i == n_total:
            print(f"  [{i:4d}/{n_total}] Processing...", end="\r")

        seqs = parse_fasta(aln_file)
        if not seqs or REFERENCE not in seqs:
            summary.append({
                "gene": gene_name,
                "n_strains": len(seqs),
                "n_nonsyn": 0,
                "n_syn": 0,
                "error": "reference_not_found"
            })
            continue

        nonsyn, syn = find_snps(seqs, REFERENCE, gene_name)
        all_nonsyn.extend(nonsyn)
        all_syn.extend(syn)

        summary.append({
            "gene":      gene_name,
            "n_strains": len(seqs),
            "n_nonsyn":  len(nonsyn),
            "n_syn":     len(syn),
            "error":     "",
        })

    print(f"\n  Done. {n_total} genes processed.\n")

    # Save outputs
    cols_snp = ["Gene", "Strain", "Position", "Reference",
                "Aa_ref", "SNP", "Aa_mod", "Type_SNP"]

    # Non-synonymous
    df_ns = pd.DataFrame(all_nonsyn, columns=cols_snp) \
            if all_nonsyn else pd.DataFrame(columns=cols_snp)
    out_ns = RAW / "allnsSNP_v2.tsv"
    df_ns.to_csv(out_ns, sep="\t", index=False)

    # Synonymous
    df_s = pd.DataFrame(all_syn, columns=cols_snp) \
           if all_syn else pd.DataFrame(columns=cols_snp)
    out_s = RAW / "allsSNP_v2.tsv"
    df_s.to_csv(out_s, sep="\t", index=False)

    # Summary
    df_sum = pd.DataFrame(summary)
    out_sum = PROC / "snp_summary_v2.tsv"
    df_sum.to_csv(out_sum, sep="\t", index=False)

    # Print global stats
    print("=" * 70)
    print("  GLOBAL SUMMARY")
    print("=" * 70)
    print(f"  Genes processed:              {n_total:,}")
    print(f"  Genes with ≥1 nsSNP:          {(df_sum['n_nonsyn']>0).sum():,}")
    print(f"  Total non-synonymous SNPs:    {len(df_ns):,}")
    print(f"  Total synonymous SNPs:        {len(df_s):,}")
    print(f"  Unique nsSNP positions:       {df_ns[['Gene','Position']].drop_duplicates().shape[0]:,}")
    print(f"\n  Non-synonymous SNPs → {out_ns}")
    print(f"  Synonymous SNPs    → {out_s}")
    print(f"  Summary            → {out_sum}")

    # Compare with original
    orig = RAW / "allnsSNP.tsv"
    if orig.exists():
        df_orig = pd.read_csv(orig, sep="\t")
        print(f"\n  Comparison with original allnsSNP.tsv:")
        print(f"    Original: {len(df_orig):,} non-synonymous SNPs")
        print(f"    New v2:   {len(df_ns):,} non-synonymous SNPs")
        diff = len(df_orig) - len(df_ns)
        print(f"    Difference: {diff:,} SNPs removed "
              f"({diff/len(df_orig)*100:.1f}% reduction)")

    print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
13_contamination_check.py
=========================
Proper contamination screening for the assemblies, replacing the crude
"GC content by contig count" heuristic.

Two complementary analyses:

  (A) LENGTH-WEIGHTED GC PROFILE
      GC is computed per contig but weighted by contig LENGTH, and the
      assembly is split into a low-GC fraction (V. vulnificus, ~46.6%) and a
      high-GC fraction. Counting contigs is misleading when N50 is small,
      because thousands of short contigs dominate the count while
      contributing little sequence. What matters is how many BASE PAIRS are
      not Vibrio-like.

  (B) 16S rRNA TAXONOMIC ASSIGNMENT  (the definitive test)
      Bakta already predicted rRNA genes for every genome. This script pulls
      the 16S sequences out of the Bakta annotation and (optionally) BLASTs
      them against NCBI's 16S_ribosomal_RNA database, which gives a
      genus/species-level identification. Multiple distinct 16S species in a
      single assembly = mixed sample, definitively.

      To get the database (small, ~250 MB):
          mamba run -n vv_annot update_blastdb.pl --decompress 16S_ribosomal_RNA

Outputs:
  <outdir>/gc_profile.tsv        per-strain length-weighted GC breakdown
  <outdir>/16S_sequences.fasta   all extracted 16S genes, named <strain>|<locus>
  <outdir>/16S_blast.tsv         BLAST hits (if a database was given)
  <outdir>/16S_assignment.tsv    best hit per strain, with n distinct taxa

Usage:
  python3 13_contamination_check.py \
      --genomes-dir data/genomes \
      --bakta-dir   results/bakta \
      --outdir      results/contamination \
      --blast-db    /path/to/16S_ribosomal_RNA        # optional
      --env-annot   vv_annot
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

VV_GC = 46.6            # V. vulnificus reference GC%
GC_WINDOW = 4.0         # +/- around VV_GC counted as "Vibrio-like"
MIN_CONTIG = 500        # ignore very short contigs for the GC profile


def gc_profile(fasta):
    """Length-weighted GC breakdown of an assembly."""
    total_bp = 0
    vib_bp = 0          # bp in contigs with GC within the Vibrio window
    high_bp = 0         # bp in contigs with GC clearly above (>+8%)
    low_bp = 0          # bp in contigs with GC clearly below (<-8%)
    gc_sum = 0.0
    contigs = 0
    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper()
        L = len(seq)
        if L < MIN_CONTIG:
            continue
        gc = 100.0 * (seq.count("G") + seq.count("C")) / L
        contigs += 1
        total_bp += L
        gc_sum += gc * L
        if abs(gc - VV_GC) <= GC_WINDOW:
            vib_bp += L
        elif gc > VV_GC + 8:
            high_bp += L
        elif gc < VV_GC - 8:
            low_bp += L
    if total_bp == 0:
        return None
    return {
        "contigs_ge_500bp": contigs,
        "total_bp": total_bp,
        "mean_GC_weighted": round(gc_sum / total_bp, 2),
        "bp_vibrio_like": vib_bp,
        "pct_vibrio_like": round(100.0 * vib_bp / total_bp, 1),
        "bp_high_GC": high_bp,
        "pct_high_GC": round(100.0 * high_bp / total_bp, 1),
        "bp_low_GC": low_bp,
        "pct_low_GC": round(100.0 * low_bp / total_bp, 1),
    }


def extract_16S(bakta_dir, strain):
    """
    Pull 16S rRNA nucleotide sequences from Bakta output.
    Bakta writes feature nucleotide sequences to <strain>.ffn and a feature
    table to <strain>.tsv (columns include Type and Product).
    """
    d = Path(bakta_dir) / strain
    tsv = d / f"{strain}.tsv"
    ffn = d / f"{strain}.ffn"
    if not tsv.exists() or not ffn.exists():
        return []

    wanted = set()
    with open(tsv) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            ftype = parts[1].strip().lower()
            product = parts[-1].strip().lower() if parts else ""
            locus = parts[5].strip()
            if "rrna" in ftype and "16s" in product:
                wanted.add(locus)
            elif "16s ribosomal rna" in product:
                wanted.add(locus)

    out = []
    for rec in SeqIO.parse(str(ffn), "fasta"):
        rid = rec.id.split()[0]
        desc = rec.description.lower()
        if rid in wanted or "16s ribosomal rna" in desc:
            out.append((rid, str(rec.seq)))
    return out


def run_blast(fasta, db, out_tsv, env=None, threads=4):
    cmd = ["blastn", "-query", str(fasta), "-db", str(db),
           "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
           "-max_target_seqs", "5", "-num_threads", str(threads),
           "-out", str(out_tsv)]
    if env:
        cmd = ["mamba", "run", "-n", env] + cmd
    print("  $ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--bakta-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--blast-db", default="")
    ap.add_argument("--env-annot", default="vv_annot")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 74)
    print("  Contamination check: length-weighted GC + 16S taxonomy")
    print("=" * 74)

    # ── (A) GC profile ───────────────────────────────────────────────────────
    print("\n[A] Length-weighted GC profile")
    print(f"    Vibrio-like window: {VV_GC} +/- {GC_WINDOW} % GC\n")
    rows = []
    genomes = sorted(Path(a.genomes_dir).glob("*.fna"))
    for g in genomes:
        s = g.stem
        prof = gc_profile(str(g))
        if prof is None:
            continue
        rows.append({"strain": s, **prof})
        flag = ""
        if prof["pct_high_GC"] > 10:
            flag = f"  <-- {prof['pct_high_GC']}% of bp are high-GC"
        print(f"    {s:<14} {prof['total_bp']/1e6:>6.2f} Mb  "
              f"GC={prof['mean_GC_weighted']:>5.1f}%  "
              f"Vibrio-like={prof['pct_vibrio_like']:>5.1f}%{flag}")
    gc_df = pd.DataFrame(rows)
    gc_df.to_csv(outdir / "gc_profile.tsv", sep="\t", index=False)

    suspect = gc_df[gc_df["pct_vibrio_like"] < 85]["strain"].tolist()
    if suspect:
        print(f"\n    Strains with <85% Vibrio-like sequence: {', '.join(suspect)}")

    # ── (B) 16S extraction ───────────────────────────────────────────────────
    print("\n[B] 16S rRNA extraction from Bakta annotation")
    fasta_out = outdir / "16S_sequences.fasta"
    n_seqs = 0
    per_strain = {}
    with open(fasta_out, "w") as fh:
        for g in genomes:
            s = g.stem
            seqs = extract_16S(a.bakta_dir, s)
            per_strain[s] = len(seqs)
            for locus, seq in seqs:
                fh.write(f">{s}|{locus}\n{seq}\n")
                n_seqs += 1
            print(f"    {s:<14} {len(seqs)} 16S copies")
    print(f"\n    wrote {n_seqs} sequences to {fasta_out}")

    if n_seqs == 0:
        print("\n    No 16S found in the Bakta output. Check the .tsv/.ffn "
              "column layout, or run with --blast-db on the raw contigs.")
        return

    # ── (C) BLAST against 16S database ───────────────────────────────────────
    if not a.blast_db:
        print("""
[C] BLAST skipped (no --blast-db given).

    To get a definitive taxonomic assignment, download NCBI's 16S database
    (small, ~250 MB) and re-run with --blast-db:

        cd <somewhere>
        mamba run -n vv_annot update_blastdb.pl --decompress 16S_ribosomal_RNA
        python3 13_contamination_check.py ... --blast-db <somewhere>/16S_ribosomal_RNA
""")
        return

    print("\n[C] BLAST 16S vs NCBI 16S_ribosomal_RNA")
    blast_out = outdir / "16S_blast.tsv"
    run_blast(fasta_out, a.blast_db, blast_out, env=a.env_annot,
              threads=a.threads)

    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
            "stitle"]
    b = pd.read_csv(blast_out, sep="\t", names=cols)
    b["strain"] = b["qseqid"].str.split("|").str[0]
    # best hit per query
    best = b.sort_values("bitscore", ascending=False).groupby("qseqid").first()
    best = best.reset_index()
    best["strain"] = best["qseqid"].str.split("|").str[0]

    def taxon(title):
        parts = str(title).split()
        # stitle usually starts with accession then organism name
        return " ".join(parts[1:3]) if len(parts) > 2 else str(title)

    best["taxon"] = best["stitle"].apply(taxon)
    best.to_csv(outdir / "16S_assignment.tsv", sep="\t", index=False)

    print("\n    -- 16S assignment per strain --")
    for s, sub in best.groupby("strain"):
        taxa = sub["taxon"].value_counts()
        n_distinct = len(taxa)
        summary = "; ".join(f"{t} ({c})" for t, c in taxa.items())
        flag = "   <-- MIXED SAMPLE" if n_distinct > 1 else ""
        print(f"    {s:<14} {summary}{flag}")

    mixed = [s for s, sub in best.groupby("strain")
             if sub["taxon"].nunique() > 1]
    if mixed:
        print(f"\n    Assemblies with more than one 16S taxon: {', '.join(mixed)}")
        print("    These are mixed sequencing samples, not just fragmented "
              "assemblies.")

    print(f"\n  wrote {outdir/'16S_assignment.tsv'}")


if __name__ == "__main__":
    main()

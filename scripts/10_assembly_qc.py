#!/usr/bin/env python3
"""
10_assembly_qc.py
=================
Assembly quality-control table for the manuscript (Fran's comment #44).

For every genome (.fna) computes standard assembly metrics:
    - number of contigs
    - total assembly length (bp)
    - largest contig (bp)
    - N50, N90
    - GC content (%)
    - number of contigs >= 1 kb
    - number of Ns (ambiguous bases)

Optionally merges Bakta annotation counts (CDS, tRNA, rRNA) if the Bakta
output directories are available, and adds the number of AMRFinderPlus /
acquired-resistance hits if the AMR summary is present.

Produces:
    results/qc/assembly_qc.tsv       (machine readable)
    results/qc/assembly_qc_table.md  (paper-ready markdown table)

Usage:
    python3 10_assembly_qc.py \
        --genomes-dir data/genomes \
        --reference   data/reference/YJ016.fna \
        --bakta-dir   results/bakta \
        --amr-summary results/amrfinder/amrfinder_summary.tsv \
        --mics        data/mic/MICs_converted.tsv \
        --outdir      results/qc

Only --genomes-dir and --outdir are strictly required.

Requirements: biopython, pandas.
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# tree/SNP name -> MIC name (to add lineage/MIC context if wanted later)
SNP_TO_MIC = {
    "11028": "BT3_11028", "12": "BT3_12", "162": "BT3_162", "32": "BT3_32",
    "94385": "94-3-85", "960426-1_4C": "960426_1/4C", "CIP8190": "C1P8190",
    "V4": "V_4", "CECT4604": "CECT_4604", "CECT4606": "CECT_4606",
    "CECT4608": "CECT_4608", "CECT4865": "CECT_4865", "CECT4866": "CECT_4866",
    "CECT4999": "CECT_4999", "CECT5198": "CECT_5198", "CECT5763": "CECT_5763",
    "CECT5769": "CECT_5769", "CECT7030": "CECT_7030", "CECT898": "CECT_898",
    "CG100": "CG_100", "94-9-112": "94-9-112", "94-9-119": "94-9-119",
    "95-8-161": "95-8-161", "95-8-6": "95-8-6", "95-8-7": "95-8-7",
    "A14": "A14", "PD-2-51": "PD-2-51", "Ra_3": "Ra_3", "RIU1": "RIU1",
    "YJ016": "YJ016",
}


def assembly_metrics(fasta):
    """Compute assembly stats from a FASTA of contigs."""
    lengths = []
    gc = 0
    ns = 0
    total = 0
    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper()
        L = len(seq)
        lengths.append(L)
        total += L
        gc += seq.count("G") + seq.count("C")
        ns += seq.count("N")
    if not lengths:
        return None
    lengths.sort(reverse=True)

    def n_x(frac):
        target = total * frac
        run = 0
        for L in lengths:
            run += L
            if run >= target:
                return L
        return lengths[-1]

    return {
        "n_contigs": len(lengths),
        "total_bp": total,
        "largest_contig": lengths[0],
        "N50": n_x(0.5),
        "N90": n_x(0.9),
        "GC_pct": round(100.0 * gc / total, 2) if total else 0,
        "n_contigs_ge_1kb": sum(1 for L in lengths if L >= 1000),
        "n_Ns": ns,
    }


def bakta_counts(bakta_dir, gid):
    """Pull CDS/tRNA/rRNA counts from Bakta's TSV summary if present."""
    txt = Path(bakta_dir) / gid / f"{gid}.txt"
    if not txt.exists():
        return {}
    out = {}
    for line in txt.read_text().splitlines():
        low = line.lower()
        for key, tag in (("cds", "CDS"), ("trna", "tRNA"),
                         ("rrna", "rRNA"), ("contigs", "bakta_contigs")):
            if low.strip().startswith(key + ":") or f"{key}:" in low[:20]:
                try:
                    out[tag] = int(line.split(":")[-1].strip().split()[0])
                except Exception:
                    pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--reference", default="")
    ap.add_argument("--bakta-dir", default="")
    ap.add_argument("--amr-summary", default="")
    ap.add_argument("--mics", default="")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(os.path.join(a.genomes_dir, "*.fna")))
    if a.reference and os.path.exists(a.reference):
        files.append(a.reference)

    amr = {}
    if a.amr_summary and os.path.exists(a.amr_summary):
        adf = pd.read_csv(a.amr_summary, sep="\t")
        amr = dict(zip(adf["strain"].astype(str), adf["n_amr_genes"]))

    rows = []
    for f in files:
        gid = os.path.basename(f).rsplit(".", 1)[0]
        m = assembly_metrics(f)
        if m is None:
            print(f"  WARNING: no sequences in {f}", file=sys.stderr)
            continue
        row = {"strain": gid, "mic_name": SNP_TO_MIC.get(gid, gid)}
        row.update(m)
        if a.bakta_dir:
            row.update(bakta_counts(a.bakta_dir, gid))
        if amr:
            row["acquired_AMR_genes"] = amr.get(gid, "NA")
        rows.append(row)
        print(f"  {gid:<28} contigs={m['n_contigs']:>5}  "
              f"N50={m['N50']:>8}  total={m['total_bp']:>9}", file=sys.stderr)

    df = pd.DataFrame(rows)
    # order: reference last, rest alphabetical
    df = df.sort_values("strain").reset_index(drop=True)

    tsv = outdir / "assembly_qc.tsv"
    df.to_csv(tsv, sep="\t", index=False)

    # paper-ready markdown
    show_cols = ["mic_name", "n_contigs", "total_bp", "largest_contig",
                 "N50", "GC_pct", "n_contigs_ge_1kb"]
    show_cols = [c for c in show_cols if c in df.columns]
    md = df[show_cols].copy()
    md.columns = ["Strain", "Contigs", "Total (bp)", "Largest (bp)",
                  "N50 (bp)", "GC (%)", "Contigs ≥1 kb"][:len(show_cols)]
    (outdir / "assembly_qc_table.md").write_text(md.to_markdown(index=False))

    # summary line for the paper
    print("\n  ── Assembly QC summary ──", file=sys.stderr)
    print(f"    genomes           : {len(df)}", file=sys.stderr)
    print(f"    contigs  median   : {int(df['n_contigs'].median())} "
          f"(range {int(df['n_contigs'].min())}–{int(df['n_contigs'].max())})",
          file=sys.stderr)
    print(f"    N50      median   : {int(df['N50'].median()):,} bp", file=sys.stderr)
    print(f"    total_bp median   : {int(df['total_bp'].median()):,} bp", file=sys.stderr)
    print(f"    GC%      median   : {df['GC_pct'].median():.1f}%", file=sys.stderr)
    print(f"\n  wrote {tsv}", file=sys.stderr)
    print(f"  wrote {outdir/'assembly_qc_table.md'}", file=sys.stderr)


if __name__ == "__main__":
    main()

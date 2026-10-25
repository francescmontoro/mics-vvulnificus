#!/usr/bin/env python3
"""
14_decontaminate.py
===================
Removes contaminant contigs from mixed assemblies, conservatively.

CONTEXT
-------
16S screening (script 13) showed 5 of 29 assemblies are mixed samples:

    95-8-6     Pseudomonas aeruginosa      ~46% of bp
    12         Renibacterium salmoninarum  ~36% of bp
    CIP8190    Renibacterium salmoninarum  ~28% of bp
    94-9-112   Renibacterium salmoninarum  ~20% of bp
    RIU1       Pseudomonas aeruginosa      ~11% of bp

Coverage in the SPAdes headers separates the two organisms cleanly in some
strains (95-8-6: Vibrio ~1.60x vs contaminant ~0.03x) but not in others
(strain 12: 1.29x vs 0.81x), so coverage alone is not a safe criterion.

CLASSIFICATION RULE (deliberately conservative)
-----------------------------------------------
A contig is discarded ONLY IF BOTH conditions hold:

    (a) it has NO significant nucleotide homology to the V. vulnificus
        reference (YJ016), AND
    (b) its GC content is outside the Vibrio range

Requiring both prevents deleting genuine V. vulnificus accessory sequence:
integrons carrying sul1/sul2, plasmids and recently acquired islands often
have atypical GC, but they are kept as long as they hit the reference — and
conversely, a divergent Vibrio contig with no reference hit is kept if its GC
is Vibrio-like. Everything ambiguous is KEPT: under-cleaning is safer than
deleting real biology.

Coverage is reported for transparency but is NOT used as a filter.

Outputs:
    <outdir>/<strain>.fna              cleaned assembly
    <outdir>/<strain>.removed.fna      discarded contigs (for inspection)
    <outdir>/decontamination_report.tsv per-strain before/after stats
    <outdir>/contig_calls_<strain>.tsv  per-contig decision and evidence

Usage:
    python3 14_decontaminate.py \
        --genomes-dir data/genomes \
        --reference   data/reference/YJ016.fna \
        --strains     95-8-6 12 CIP8190 94-9-112 RIU1 \
        --outdir      data/genomes_clean \
        --env-annot   vv_annot
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# V. vulnificus GC range. Wide on purpose: only contigs OUTSIDE this range
# AND without a reference hit are discarded.
GC_MIN, GC_MAX = 41.0, 51.5
MIN_CONTIG = 200          # contigs shorter than this are dropped as assembly junk
BLAST_MIN_COV = 20.0      # % of contig length covered by hits to the reference
BLAST_MIN_PID = 80.0      # % identity


def contig_gc(seq):
    s = seq.upper()
    L = len(s)
    return 100.0 * (s.count("G") + s.count("C")) / L if L else 0.0


def contig_cov(header):
    m = re.search(r"cov_([0-9.]+)", header)
    return float(m.group(1)) if m else float("nan")


def blast_against_reference(fasta, ref, env=None, threads=4):
    """
    Returns {contig_id: (pct_of_contig_covered_by_hits, best_pident)}.
    Uses blastn (megablast) of the assembly against the reference genome.
    """
    tmp = tempfile.mkdtemp(prefix="decontam_")
    db = os.path.join(tmp, "ref")
    out = os.path.join(tmp, "hits.tsv")

    def run(cmd):
        if env:
            cmd = ["mamba", "run", "-n", env] + cmd
        subprocess.run(cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    run(["makeblastdb", "-in", str(ref), "-dbtype", "nucl", "-out", db])
    cmd = ["blastn", "-query", str(fasta), "-db", db,
           "-outfmt", "6 qseqid qlen pident length qstart qend",
           "-evalue", "1e-10", "-num_threads", str(threads),
           "-max_target_seqs", "5", "-out", out]
    if env:
        cmd = ["mamba", "run", "-n", env] + cmd
    subprocess.run(cmd, check=True)

    # merge hit intervals per contig to get covered fraction
    spans = {}
    qlen = {}
    best_pid = {}
    with open(out) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            q, ql, pid, ln, qs, qe = p[0], int(p[1]), float(p[2]), int(p[3]), int(p[4]), int(p[5])
            if pid < BLAST_MIN_PID:
                continue
            qlen[q] = ql
            a, b = min(qs, qe), max(qs, qe)
            spans.setdefault(q, []).append((a, b))
            best_pid[q] = max(best_pid.get(q, 0.0), pid)

    result = {}
    for q, ivs in spans.items():
        ivs.sort()
        merged = []
        for a, b in ivs:
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        covered = sum(b - a + 1 for a, b in merged)
        result[q] = (100.0 * covered / qlen[q], best_pid[q])
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--strains", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--env-annot", default="vv_annot")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("  Decontamination of mixed assemblies")
    print("=" * 78)
    print(f"  rule: discard a contig only if NO hit to {Path(a.reference).name}"
          f" AND GC outside [{GC_MIN}, {GC_MAX}]%")
    print(f"  reference hit = >={BLAST_MIN_COV}% of contig covered at "
          f">={BLAST_MIN_PID}% identity\n")

    report = []
    for s in a.strains:
        fasta = Path(a.genomes_dir) / f"{s}.fna"
        if not fasta.exists():
            print(f"  [SKIP] {s}: {fasta} not found")
            continue

        print(f"  --- {s} ---")
        hits = blast_against_reference(fasta, a.reference,
                                       env=a.env_annot, threads=a.threads)

        kept, removed, calls = [], [], []
        for rec in SeqIO.parse(str(fasta), "fasta"):
            seq = str(rec.seq)
            L = len(seq)
            gc = contig_gc(seq)
            cov = contig_cov(rec.id)
            hcov, hpid = hits.get(rec.id, (0.0, 0.0))
            has_hit = hcov >= BLAST_MIN_COV
            gc_ok = GC_MIN <= gc <= GC_MAX

            if L < MIN_CONTIG:
                decision, reason = "removed", "too short"
            elif has_hit or gc_ok:
                decision, reason = "kept", ("reference hit" if has_hit
                                            else "Vibrio-like GC")
            else:
                decision, reason = "removed", "no reference hit + atypical GC"

            calls.append({"contig": rec.id, "length": L, "GC": round(gc, 1),
                          "spades_cov": cov, "ref_cov_pct": round(hcov, 1),
                          "ref_pident": round(hpid, 1),
                          "decision": decision, "reason": reason})
            (kept if decision == "kept" else removed).append(rec)

        SeqIO.write(kept, str(outdir / f"{s}.fna"), "fasta")
        SeqIO.write(removed, str(outdir / f"{s}.removed.fna"), "fasta")
        pd.DataFrame(calls).to_csv(outdir / f"contig_calls_{s}.tsv",
                                   sep="\t", index=False)

        bp_in = sum(len(r.seq) for r in kept) + sum(len(r.seq) for r in removed)
        bp_out = sum(len(r.seq) for r in kept)
        gc_out = (contig_gc("".join(str(r.seq) for r in kept))
                  if kept else float("nan"))
        print(f"    contigs  {len(kept)+len(removed):>6} -> {len(kept):>6}"
              f"   ({len(removed)} removed)")
        print(f"    bases    {bp_in/1e6:>6.2f} Mb -> {bp_out/1e6:.2f} Mb"
              f"   ({100.0*(bp_in-bp_out)/bp_in:.1f}% removed)")
        print(f"    GC       -> {gc_out:.1f}%\n")

        report.append({"strain": s,
                       "contigs_before": len(kept) + len(removed),
                       "contigs_after": len(kept),
                       "bp_before": bp_in, "bp_after": bp_out,
                       "pct_removed": round(100.0*(bp_in-bp_out)/bp_in, 1),
                       "GC_after": round(gc_out, 2)})

    if report:
        df = pd.DataFrame(report)
        df.to_csv(outdir / "decontamination_report.tsv", sep="\t", index=False)
        print(f"  wrote {outdir/'decontamination_report.tsv'}")
        print("""
  NEXT STEPS
    1. Sanity-check the cleaned assemblies: sizes should be ~5.0-5.4 Mb and
       GC ~46-47%. Inspect <strain>.removed.fna if anything looks off.
    2. Copy the cleaned files over the originals (keep a backup), then
       re-run Bakta for these strains only, and re-run Panaroo + Snippy so
       the core genome and the variant matrix are built from clean data.
    3. Only then re-run pyseer / Scoary.
""")


if __name__ == "__main__":
    main()

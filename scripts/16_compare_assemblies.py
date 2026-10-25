#!/usr/bin/env python3
"""
16_compare_assemblies.py
========================
Compares the NCBI-deposited assemblies against our locally decontaminated
versions, for the strains present in both, and decides whether the sequence
missing from the NCBI versions is genuine V. vulnificus.

WHY
---
Several deposited assemblies are much smaller than expected for the species
(Ra3 3.73 Mb, Riu1 4.31 Mb, 94-8-112 4.33 Mb, 95-8-6 4.40 Mb vs ~5.2 Mb for
YJ016). Two explanations are possible and they lead to opposite decisions:

  (a) NCBI's contamination screen removed foreign sequence
      -> the deposited versions are correct and should be used.
  (b) The deposited versions are simply less complete
      -> using them would shrink the Panaroo core genome and cost us sites
         for both the tree and the association analysis.

The test: take the sequence that is present in OUR version but absent from the
NCBI one, and ask whether it looks like V. vulnificus. If that extra sequence
has homology to the YJ016 reference and Vibrio-like GC, explanation (b) holds.
If it is high-GC with no reference homology, explanation (a) holds.

METHOD
------
For each shared strain:
  1. blastn our cleaned assembly against the NCBI assembly
  2. classify each of our contigs as "shared" (well covered by NCBI) or
     "extra" (poorly covered / absent)
  3. for the extra sequence, compute GC and homology to YJ016

Usage:
  python3 16_compare_assemblies.py \
      --ours     data/genomes \
      --ncbi-dir data/genomes_ncbi/ncbi_dataset/data \
      --reference data/reference/YJ016.fna \
      --outdir   results/assembly_comparison \
      --env-annot vv_annot
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# NCBI strain label -> our local file stem
# note 94-8-112 (NCBI, and the 2018 paper) == 94-9-112 (our local naming);
# confirmed by Fran as a typo on our side.
NCBI_TO_LOCAL = {
    "Ra3": "Ra_3",
    "PD-2-51": "PD-2-51",
    "CG100": "CG100",
    "CECT7030": "CECT7030",
    "CECT4865": "CECT4865",
    "A14": "A14",
    "CECT4604": "CECT4604",
    "960426-1_4C": "960426-1_4C",
    "94-9-119": "94-9-119",
    "95-8-6": "95-8-6",
    "94385": "94385",
    "CECT898": "CECT898",
    "CECT4608": "CECT4608",
    "95-8-7": "95-8-7",
    "95-8-161": "95-8-161",
    "94-8-112": "94-9-112",
    "162": "162",
    "32": "32",
    "Riu1": "RIU1",
}

MIN_CONTIG = 500
SHARED_COV = 80.0      # % of our contig covered by NCBI hits -> "shared"
BLAST_PID = 95.0       # same strain, so identity should be very high


def run(cmd, env=None, **kw):
    if env:
        cmd = ["mamba", "run", "-n", env] + cmd
    return subprocess.run(cmd, check=True, **kw)


def coverage_map(query_fa, subject_fa, env, threads, pid_min):
    """% of each query contig covered by hits to subject."""
    tmp = tempfile.mkdtemp(prefix="cmp_")
    db = os.path.join(tmp, "subj")
    out = os.path.join(tmp, "hits.tsv")
    run(["makeblastdb", "-in", str(subject_fa), "-dbtype", "nucl", "-out", db],
        env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    run(["blastn", "-query", str(query_fa), "-db", db,
         "-outfmt", "6 qseqid qlen pident qstart qend",
         "-evalue", "1e-20", "-num_threads", str(threads),
         "-max_target_seqs", "20", "-out", out], env=env)

    spans, qlen = {}, {}
    with open(out) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            q, ql, pid, qs, qe = p[0], int(p[1]), float(p[2]), int(p[3]), int(p[4])
            if pid < pid_min:
                continue
            qlen[q] = ql
            a, b = min(qs, qe), max(qs, qe)
            spans.setdefault(q, []).append((a, b))

    cov = {}
    for q, ivs in spans.items():
        ivs.sort()
        merged = []
        for a, b in ivs:
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        cov[q] = 100.0 * sum(b - a + 1 for a, b in merged) / qlen[q]
    return cov


def gc(seq):
    s = seq.upper()
    return 100.0 * (s.count("G") + s.count("C")) / len(s) if s else 0.0


def find_ncbi_fasta(ncbi_dir, accession_prefix):
    """Locate the .fna for a given accession folder (prefer GCF)."""
    for acc in sorted(os.listdir(ncbi_dir), reverse=True):   # GCF before GCA
        if not acc.startswith(accession_prefix):
            continue
        d = os.path.join(ncbi_dir, acc)
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if f.endswith(".fna"):
                return os.path.join(d, f)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--ncbi-dir", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--env-annot", default="vv_annot")
    ap.add_argument("--threads", type=int, default=6)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # map strain -> ncbi fasta, from the assembly report
    report = Path(a.ncbi_dir) / "assembly_data_report.jsonl"
    strain_to_acc = {}
    if report.exists():
        for line in open(report):
            d = json.loads(line)
            acc = d.get("accession", "")
            st = (d.get("organism", {}).get("infraspecificNames", {}) or {}).get("strain", "")
            if st and (st not in strain_to_acc or acc.startswith("GCF")):
                strain_to_acc[st] = acc

    print("=" * 78)
    print("  NCBI deposited assemblies vs our decontaminated versions")
    print("=" * 78)
    print(f"  a contig of ours is 'shared' if >={SHARED_COV}% of it is covered "
          f"by NCBI hits at >={BLAST_PID}% identity\n")

    rows = []
    for ncbi_name, local in sorted(NCBI_TO_LOCAL.items()):
        acc = strain_to_acc.get(ncbi_name)
        if not acc:
            print(f"  [SKIP] {ncbi_name}: no accession in report")
            continue
        ncbi_fa = find_ncbi_fasta(a.ncbi_dir, acc)
        ours_fa = Path(a.ours) / f"{local}.fna"
        if not ncbi_fa or not ours_fa.exists():
            print(f"  [SKIP] {ncbi_name}: missing file")
            continue

        cov = coverage_map(ours_fa, ncbi_fa, a.env_annot, a.threads, BLAST_PID)

        extra_recs, shared_bp, extra_bp = [], 0, 0
        for rec in SeqIO.parse(str(ours_fa), "fasta"):
            L = len(rec.seq)
            if L < MIN_CONTIG:
                continue
            if cov.get(rec.id, 0.0) >= SHARED_COV:
                shared_bp += L
            else:
                extra_bp += L
                extra_recs.append(rec)

        # is the extra sequence Vibrio?
        extra_fa = outdir / f"{local}.extra.fna"
        SeqIO.write(extra_recs, str(extra_fa), "fasta")
        extra_gc = gc("".join(str(r.seq) for r in extra_recs)) if extra_recs else float("nan")

        ref_cov = {}
        if extra_recs:
            ref_cov = coverage_map(extra_fa, a.reference, a.env_annot,
                                   a.threads, 80.0)
        extra_with_ref_hit = sum(len(r.seq) for r in extra_recs
                                 if ref_cov.get(r.id, 0.0) >= 20.0)

        ncbi_bp = sum(len(r.seq) for r in SeqIO.parse(ncbi_fa, "fasta"))
        ours_bp = shared_bp + extra_bp

        pct_ref = 100.0 * extra_with_ref_hit / extra_bp if extra_bp else 0.0
        verdict = ""
        if extra_bp > 200_000:
            verdict = ("NCBI INCOMPLETO" if pct_ref > 50 else "contaminacion")

        print(f"  {local:<14} NCBI {ncbi_bp/1e6:>5.2f} Mb | nuestra {ours_bp/1e6:>5.2f} Mb"
              f" | extra {extra_bp/1e6:>5.2f} Mb  GC={extra_gc:>5.1f}%"
              f"  con homologia YJ016: {pct_ref:>5.1f}%   {verdict}")

        rows.append({
            "strain": local, "ncbi_accession": acc,
            "ncbi_mb": round(ncbi_bp/1e6, 2), "ours_mb": round(ours_bp/1e6, 2),
            "extra_mb": round(extra_bp/1e6, 2),
            "extra_GC": round(extra_gc, 1) if extra_gc == extra_gc else "",
            "extra_pct_with_YJ016_homology": round(pct_ref, 1),
            "verdict": verdict,
        })

    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(outdir / "assembly_comparison.tsv", sep="\t", index=False)
        print(f"\n  wrote {outdir/'assembly_comparison.tsv'}")
        print("""
  HOW TO READ THIS
    extra = sequence present in our assembly but not in the NCBI one.
    If that extra sequence has high homology to YJ016 and Vibrio-like GC
    (~46-47%), the deposited assembly is simply less complete and ours is
    preferable. If it is high-GC with little YJ016 homology, NCBI removed
    contamination and the deposited version should be used.
""")


if __name__ == "__main__":
    main()

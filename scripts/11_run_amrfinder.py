#!/usr/bin/env python3
"""
11_run_amrfinder.py
===================
Runs NCBI AMRFinderPlus on all genomes, properly this time.

WHY THIS SCRIPT EXISTS
----------------------
The first AMRFinder run inside run_vv_variants.sh silently failed on all 30
genomes ("No valid AMRFinder database is found") and wrote empty files, which
were then counted as "0 acquired resistance genes". That would have been a
false result in the manuscript. This script:

  * points explicitly at the AMRFinderPlus database that Bakta downloaded
  * uses --organism Vibrio_vulnificus, which IS supported (checked with
    `amrfinder -l`), enabling detection of KNOWN POINT MUTATIONS in addition
    to acquired genes
  * runs in COMBINED mode (protein + nucleotide + GFF) when possible, which is
    more sensitive on fragmented assemblies where a gene split across two
    contigs is not called as a complete CDS (our worst assembly has 5194
    contigs) -- falls back to protein-only, then nucleotide-only
  * FAILS LOUDLY: any genome whose run errors out is reported and excluded,
    never written as a silent zero

Outputs:
    <outdir>/<strain>.amr.tsv       per-genome AMRFinder output
    <outdir>/amrfinder_summary.tsv  strain, status, n_genes, n_point_mutations, genes
    <outdir>/amrfinder_failures.txt list of genomes that failed (should be empty)

Usage:
    python3 11_run_amrfinder.py \
        --bakta-dir  results/bakta \
        --db         data/bakta_db_light/db-light/amrfinderplus-db/latest \
        --organism   Vibrio_vulnificus \
        --outdir     results/amrfinder \
        --env-annot  vv_annot \
        --threads    12
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


def run_capture(cmd, env=None):
    """Run a command, return (returncode, stdout, stderr). Never raises."""
    if env:
        cmd = ["mamba", "run", "-n", env] + cmd
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def strip_fasta_from_gff(gff_in, gff_out):
    """
    Bakta's GFF3 has the sequences appended after '##FASTA'.
    AMRFinder wants annotation only, so we write a copy without that block.
    """
    lines = []
    with open(gff_in) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            lines.append(line)
    gff_out.write_text("".join(lines))
    return gff_out


def parse_amr_output(path):
    """
    Parse AMRFinder TSV. Returns (n_acquired, n_point, gene_list).
    AMRFinder marks point mutations in the 'Element subtype' column as
    'POINT'; acquired genes are 'AMR' / other subtypes.
    """
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return 0, 0, []
    if df.empty:
        return 0, 0, []

    # Column names differ slightly across versions; find them robustly
    subtype_col = next((c for c in df.columns
                        if c.lower().replace(" ", "_") in
                        ("element_subtype", "subtype")), None)
    symbol_col = next((c for c in df.columns
                       if c.lower().replace(" ", "_") in
                       ("element_symbol", "gene_symbol", "symbol")), None)

    if subtype_col is not None:
        is_point = df[subtype_col].astype(str).str.upper().str.contains("POINT")
        n_point = int(is_point.sum())
        n_acq = int((~is_point).sum())
    else:
        n_point, n_acq = 0, len(df)

    genes = []
    if symbol_col is not None:
        genes = [str(g) for g in df[symbol_col].dropna().tolist()]
    return n_acq, n_point, genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bakta-dir", required=True)
    ap.add_argument("--db", required=True, help="AMRFinderPlus database dir")
    ap.add_argument("--organism", default="Vibrio_vulnificus")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--env-annot", default="vv_annot")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--no-organism", action="store_true",
                    help="skip --organism (no point-mutation screening)")
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    tmpdir = outdir / "_tmp"; tmpdir.mkdir(exist_ok=True)

    bakta = Path(a.bakta_dir)
    genomes = sorted([d.name for d in bakta.iterdir() if d.is_dir()])
    print("=" * 68)
    print(f"  AMRFinderPlus on {len(genomes)} genomes")
    print(f"  database : {a.db}")
    print(f"  organism : {'(none)' if a.no_organism else a.organism}")
    print("=" * 68)

    rows, failures = [], []

    for gid in genomes:
        faa = bakta / gid / f"{gid}.faa"
        fna = bakta / gid / f"{gid}.fna"
        gff = bakta / gid / f"{gid}.gff3"
        out = outdir / f"{gid}.amr.tsv"

        if not faa.exists():
            print(f"  [SKIP] {gid}: no .faa")
            failures.append(f"{gid}: missing .faa")
            continue

        base = ["amrfinder", "-d", a.db, "--plus",
                "--threads", str(a.threads), "-o", str(out)]
        if not a.no_organism:
            base += ["--organism", a.organism]

        attempts = []
        # 1) combined: protein + nucleotide + gff (most sensitive)
        if fna.exists() and gff.exists():
            gff_clean = strip_fasta_from_gff(gff, tmpdir / f"{gid}.noseq.gff3")
            attempts.append(base + ["-p", str(faa), "-n", str(fna),
                                    "-g", str(gff_clean),
                                    "--annotation_format", "bakta"])
            attempts.append(base + ["-p", str(faa), "-n", str(fna),
                                    "-g", str(gff_clean)])
        # 2) protein only
        attempts.append(base + ["-p", str(faa)])
        # 3) nucleotide only
        if fna.exists():
            attempts.append(base + ["-n", str(fna)])

        ok, mode_used, last_err = False, "", ""
        for i, cmd in enumerate(attempts, 1):
            rc, so, se = run_capture(cmd, env=a.env_annot)
            if rc == 0 and out.exists() and out.stat().st_size > 0:
                ok = True
                mode_used = ("combined" if "-n" in cmd and "-p" in cmd
                             else ("protein" if "-p" in cmd else "nucleotide"))
                break
            last_err = (se or so).strip().splitlines()[-1] if (se or so) else "unknown error"

        if not ok:
            print(f"  [FAIL] {gid}: {last_err[:110]}")
            failures.append(f"{gid}: {last_err}")
            continue

        n_acq, n_point, genes = parse_amr_output(out)
        rows.append({
            "strain": gid,
            "status": "OK",
            "mode": mode_used,
            "n_acquired_genes": n_acq,
            "n_point_mutations": n_point,
            "genes": ";".join(genes),
        })
        flag = "" if (n_acq or n_point) else "  (none)"
        print(f"  [ OK ] {gid:<12} mode={mode_used:<10} "
              f"acquired={n_acq}  point={n_point}{flag}")

    # ── outputs ───────────────────────────────────────────────────────────
    if rows:
        df = pd.DataFrame(rows).sort_values("strain")
        df.to_csv(outdir / "amrfinder_summary.tsv", sep="\t", index=False)
        print("\n  ── Summary ──")
        print(f"    genomes analysed        : {len(df)}")
        print(f"    with acquired AMR genes : {(df['n_acquired_genes']>0).sum()}")
        print(f"    with point mutations    : {(df['n_point_mutations']>0).sum()}")
        allg = sorted({g for s in df['genes'] for g in s.split(';') if g})
        if allg:
            print(f"    distinct elements found : {len(allg)}")
            print(f"      {', '.join(allg[:20])}")
        print(f"\n  wrote {outdir/'amrfinder_summary.tsv'}")

    (outdir / "amrfinder_failures.txt").write_text("\n".join(failures))
    if failures:
        print(f"\n  *** {len(failures)} GENOMES FAILED — see "
              f"{outdir/'amrfinder_failures.txt'} ***")
        print("  Do NOT report these as zero. Fix and re-run.")
        sys.exit(1)
    else:
        print("\n  All genomes completed successfully.")


if __name__ == "__main__":
    main()

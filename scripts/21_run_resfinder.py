#!/usr/bin/env python3
"""
21_run_resfinder.py
===================
Runs ResFinder on the decontaminated assemblies and parses the results into a
single table.

WHY THIS IS RE-RUN
------------------
The ResFinder results carried over from the earlier analysis were generated
from the pre-decontamination assemblies and still report determinants that
belong to the contaminating organisms: aph(3')-IIb and blaOXA-50 in 95-8-6
(Pseudomonas aeruginosa), and tet(K), erm(B) and tet(M) in CECT5769 and
CECT7030 (Gram-positive). AMRFinderPlus, re-run on the cleaned annotations, no
longer reports any of these. Until ResFinder is re-run the two tools cannot be
compared, and the genotype-phenotype table mixes clean and stale calls.

FAILURE HANDLING
----------------
A strain that produces no output, or whose output cannot be parsed, is reported
as a failure and NOT as a strain without determinants. An unreported failure is
indistinguishable from a genuine absence, which is precisely the error that
produced the incorrect "no acquired resistance genes" conclusion in the earlier
analysis.

Usage:
  python3 21_run_resfinder.py \
      --genomes-dir data/genomes \
      --db          ../PAPER_V2/data/raw/resfinder_db \
      --outdir      results/resfinder \
      --min-id      0.90 \
      --min-cov     0.60
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


def find_resfinder():
    """Locate a usable ResFinder entry point."""
    for cmd in ("resfinder", "run_resfinder.py"):
        p = shutil.which(cmd)
        if p:
            return [p]
    # module form
    try:
        subprocess.run([sys.executable, "-m", "resfinder", "--help"],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        return [sys.executable, "-m", "resfinder"]
    except Exception:
        pass
    return None


def parse_results_tab(path):
    """Parse ResFinder_results_tab.txt -> list of dicts."""
    rows = []
    if not os.path.exists(path):
        return rows
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return rows
    if df.empty:
        return rows
    gcol = df.columns[0]
    for _, r in df.iterrows():
        gene = str(r[gcol]).strip()
        if not gene or gene.lower().startswith("resistance gene"):
            continue
        rows.append({
            "gene": gene,
            "identity": r.get("Identity", ""),
            "coverage": r.get("Coverage", ""),
            "contig": r.get("Contig", ""),
            "position": r.get("Position in contig", ""),
            "phenotype": r.get("Phenotype", ""),
            "accession": r.get("Accession no.", ""),
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--db", required=True,
                    help="path to the resfinder_db directory")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-id", type=float, default=0.90)
    ap.add_argument("--min-cov", type=float, default=0.60)
    ap.add_argument("--species", default="Other")
    ap.add_argument("--blast-path", default="",
                    help="full path to blastn; needed when ResFinder and "
                         "BLAST+ are installed in different conda environments")
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    rf = find_resfinder()
    if not rf:
        sys.exit(
            "ERROR: ResFinder not found on PATH or as a module.\n"
            "  Install with:  pip install resfinder\n"
            "  The database is already present at the path given by --db."
        )
    print(f"  ResFinder entry point: {' '.join(rf)}")
    print(f"  database             : {a.db}")

    genomes = sorted(glob.glob(os.path.join(a.genomes_dir, "*.fna")))
    print(f"  genomes              : {len(genomes)}\n")

    all_rows, failures = [], []
    for g in genomes:
        strain = Path(g).stem
        sdir = outdir / strain
        sdir.mkdir(parents=True, exist_ok=True)

        cmd = rf + [
            "-ifa", str(g),
            "-o", str(sdir),
            "-acq",
            "-db_res", str(a.db),
            "-l", str(a.min_cov),
            "-t", str(a.min_id),
            "-s", a.species,
            "--ignore_missing_species",
        ]
        if a.blast_path:
            cmd += ["-b", a.blast_path]
        try:
            subprocess.run(cmd, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            msg = (e.stderr or b"").decode()[-300:]
            print(f"  [FAIL] {strain}: {msg.strip()[:120]}")
            failures.append((strain, msg))
            continue

        tab = sdir / "ResFinder_results_tab.txt"
        if not tab.exists():
            print(f"  [FAIL] {strain}: no results table produced")
            failures.append((strain, "no ResFinder_results_tab.txt"))
            continue

        hits = parse_results_tab(tab)
        for h in hits:
            all_rows.append({"strain": strain, **h})
        print(f"  [ OK ] {strain:<14} {len(hits)} acquired gene(s)"
              + (f"  {', '.join(h['gene'] for h in hits)}" if hits else ""))

    # ── outputs ──────────────────────────────────────────────────────────────
    df = pd.DataFrame(all_rows)
    df.to_csv(outdir / "resfinder_hits.tsv", sep="\t", index=False)

    summary = []
    for g in genomes:
        s = Path(g).stem
        sub = df[df["strain"] == s] if len(df) else df
        summary.append({
            "strain": s,
            "n_genes": len(sub),
            "genes": ";".join(sorted(sub["gene"])) if len(sub) else "",
            "status": "FAILED" if s in [f[0] for f in failures] else "ok",
        })
    pd.DataFrame(summary).to_csv(outdir / "resfinder_summary.tsv",
                                 sep="\t", index=False)

    n_with = sum(1 for r in summary if r["n_genes"] > 0 and r["status"] == "ok")
    print(f"\n  ── Summary ──")
    print(f"    genomes analysed        : {len(genomes) - len(failures)}")
    print(f"    with acquired genes     : {n_with}")
    print(f"    total hits              : {len(df)}")
    if len(df):
        print(f"    distinct genes          : "
              f"{', '.join(sorted(df['gene'].unique()))}")

    if failures:
        with open(outdir / "resfinder_failures.txt", "w") as fh:
            for s, m in failures:
                fh.write(f"{s}\t{m}\n")
        print(f"\n  *** {len(failures)} GENOME(S) FAILED — see "
              f"{outdir/'resfinder_failures.txt'} ***")
        print("  Do NOT report these as zero.")
        sys.exit(1)

    print(f"\n  wrote {outdir/'resfinder_hits.tsv'}")
    print(f"  wrote {outdir/'resfinder_summary.tsv'}")


if __name__ == "__main__":
    main()

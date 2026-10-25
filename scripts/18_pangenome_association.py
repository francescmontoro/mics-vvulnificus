#!/usr/bin/env python3
"""
18_pangenome_association.py
===========================
Runs the association analysis on GENE PRESENCE/ABSENCE from the Panaroo
pan-genome, alongside the core-genome SNP analysis.

WHY THIS MATTERS (Fran's point)
-------------------------------
Until now only core-genome SNPs were tested (Snippy against YJ016). Accessory
gene content was never tested, even though gene gain/loss is the dominant
mechanism of acquired resistance and Panaroo already produces the matrix.

BUILT-IN POSITIVE CONTROL
-------------------------
sul1, sul2, dfrA1 and tet(59) are accessory genes carried by a known subset of
strains whose phenotype is confirmed (all sulphonamide determinants sit at
MIC > 256). If the pan-genome analysis recovers them as the top hits for
sulphathiazole, the whole association pipeline is validated end to end. If it
does NOT recover them, then the method demonstrably lacks the power to detect
a determinant we already know is present and phenotypically expressed — which
is itself a quantitative demonstration of the power ceiling, and arguably a
stronger result than the negative finding on SNPs alone.

Either outcome is informative. That is why this analysis is worth running.

Usage:
  python3 18_pangenome_association.py \
      --rtab      results/panaroo/gene_presence_absence.Rtab \
      --pangenome results/panaroo/gene_presence_absence.csv \
      --kinship   results/pyseer/phylogeny_K.tsv \
      --phenos    results/pyseer/pheno_STZ_bin128.tsv \
                  results/pyseer/pheno_STZ_ordinal.tsv \
                  results/pyseer/pheno_GEN_bin.tsv \
      --outdir    results/pangenome_assoc \
      --env       vv_assoc
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# genes whose presence should track sulphonamide resistance if the method works
POSITIVE_CONTROLS = ["sul1", "sul2", "dfrA", "tet(59)", "tet59", "aadA"]


def safe(name):
    s = str(name).strip().replace("/", "_").replace("-", "_").replace(" ", "_")
    return "S_" + s


def run(cmd, env=None, **kw):
    if env:
        cmd = ["mamba", "run", "-n", env] + cmd
    print("  $ " + " ".join(str(c) for c in cmd), file=sys.stderr)
    return subprocess.run(cmd, check=True, **kw)


def prepare_rtab(rtab_in, outdir):
    """Rename columns to safe sample names; keep pyseer's expected layout."""
    df = pd.read_csv(rtab_in, sep="\t", index_col=0)
    df.columns = [safe(c) for c in df.columns]
    # drop genes present in all or absent in all: no variance, no test
    keep = df[(df.sum(axis=1) > 0) & (df.sum(axis=1) < df.shape[1])]
    out = outdir / "pangenome.Rtab"
    keep.to_csv(out, sep="\t")
    print(f"  pan-genome: {df.shape[0]} genes -> {keep.shape[0]} variable "
          f"across {df.shape[1]} strains", file=sys.stderr)
    print(f"    (removed {df.shape[0]-keep.shape[0]} core/absent genes with no "
          f"variance)", file=sys.stderr)
    return out, keep


def gene_annotations(pangenome_csv):
    """Map Panaroo gene id -> annotation, if the csv is available."""
    ann = {}
    if not pangenome_csv or not os.path.exists(pangenome_csv):
        return ann
    try:
        df = pd.read_csv(pangenome_csv, low_memory=False)
        gcol = df.columns[0]
        acol = next((c for c in df.columns[:4]
                     if "annot" in c.lower() or "description" in c.lower()), None)
        for _, r in df.iterrows():
            ann[str(r[gcol])] = str(r[acol]) if acol else ""
    except Exception as e:
        print(f"    could not read pan-genome csv: {e}", file=sys.stderr)
    return ann


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rtab", required=True)
    ap.add_argument("--pangenome", default="")
    ap.add_argument("--kinship", required=True)
    ap.add_argument("--phenos", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--env", default="vv_assoc")
    ap.add_argument("--alpha", type=float, default=0.05)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # WSL: pyseer's multiprocessing pool deadlocks with threaded OpenBLAS
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
              "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[v] = "1"

    print("=" * 76)
    print("  Pan-genome (gene presence/absence) association")
    print("=" * 76)

    rtab, mat = prepare_rtab(a.rtab, outdir)
    ann = gene_annotations(a.pangenome)

    summary = []
    for ph in a.phenos:
        ph = Path(ph)
        if not ph.exists():
            print(f"  [SKIP] {ph} not found")
            continue
        name = ph.stem.replace("pheno_", "")
        pdf = pd.read_csv(ph, sep="\t")
        col = pdf.columns[1]
        is_binary = set(pdf[col].dropna().unique()) <= {0, 1}

        print(f"\n  ── {name} ({'binary' if is_binary else 'continuous'}) ──")
        if is_binary:
            print(f"    resistant: {int(pdf[col].sum())} / {len(pdf)}")

        out_tsv = outdir / f"pyseer_pan_{name}.tsv"
        pat = outdir / f"patterns_pan_{name}.txt"
        cmd = ["pyseer", "--lmm",
               "--phenotypes", str(ph), "--phenotype-column", col,
               "--pres", str(rtab), "--similarity", str(a.kinship),
               "--output-patterns", str(pat), "--cpu", "1"]
        if not is_binary:
            cmd.append("--continuous")
        try:
            with open(out_tsv, "w") as fh:
                run(cmd, env=a.env, stdout=fh)
        except subprocess.CalledProcessError as e:
            print(f"    pyseer failed: {e}")
            continue

        res = pd.read_csv(out_tsv, sep="\t")
        pcol = "lrt-pvalue" if "lrt-pvalue" in res.columns else res.columns[3]
        res = res.dropna(subset=[pcol]).copy()
        if not len(res):
            print("    no results")
            continue

        # block-level correction on unique patterns
        n_pat = len(res)
        if pat.exists():
            pats = [l.strip() for l in open(pat) if l.strip()]
            if len(pats) == len(res):
                res["pattern"] = pats
                res = res.sort_values(pcol).drop_duplicates("pattern")
                n_pat = len(res)
        rej, q, _, _ = multipletests(res[pcol], alpha=a.alpha, method="fdr_bh")
        res["fdr_bh"] = q
        res["significant"] = rej
        res["annotation"] = res["variant"].map(ann).fillna("")
        res = res.sort_values(pcol)
        res.to_csv(outdir / f"blocks_pan_{name}.tsv", sep="\t", index=False)

        n_sig = int(rej.sum())
        print(f"    unique patterns (tests): {n_pat}")
        print(f"    Bonferroni threshold   : {a.alpha/n_pat:.3e}")
        print(f"    significant after BH   : {n_sig}")

        print("    top 8 hits:")
        for _, r in res.head(8).iterrows():
            mark = " *" if r["significant"] else "  "
            an = str(r["annotation"])[:44]
            print(f"     {mark} {str(r['variant'])[:28]:<30} p={r[pcol]:.2e} "
                  f"q={r['fdr_bh']:.3f}  {an}")

        # ── positive controls ───────────────────────────────────────────────
        hits = res[res["variant"].astype(str).str.contains(
            "|".join(POSITIVE_CONTROLS), case=False, regex=True, na=False)]
        print("    POSITIVE CONTROLS (sul/dfrA/tet59/aadA):")
        if len(hits):
            for _, r in hits.head(10).iterrows():
                rank = int((res[pcol] < r[pcol]).sum()) + 1
                mark = "SIGNIFICANT" if r["significant"] else "not significant"
                print(f"      {str(r['variant'])[:30]:<32} p={r[pcol]:.2e} "
                      f"rank {rank}/{n_pat}  {mark}")
        else:
            print("      none found in the pan-genome by name; check the "
                  "Panaroo annotation column")

        summary.append({"phenotype": name, "tests": n_pat,
                        "significant": n_sig,
                        "min_p": float(res[pcol].min()),
                        "bonferroni": a.alpha / n_pat})

    if summary:
        s = pd.DataFrame(summary)
        s.to_csv(outdir / "pangenome_summary.tsv", sep="\t", index=False)
        print("\n" + "=" * 76)
        print(s.to_string(index=False))
        print(f"\n  wrote {outdir/'pangenome_summary.tsv'}")

    print("""
  HOW TO READ THE POSITIVE CONTROLS
    sul1/sul2 are carried by strains whose sulphathiazole MIC is > 256, so a
    working method should rank them at or near the top for that phenotype.
      - recovered and significant -> the pipeline detects real determinants,
        and the negative SNP result is a statement about chromosomal variation
        specifically, not a failure of the method.
      - recovered but not significant -> quantifies the power ceiling: the
        method cannot reach significance even for a determinant that is present,
        phenotypically expressed and mechanistically established.
      - not recovered at all -> check the Panaroo gene naming before concluding
        anything; accessory genes are often labelled group_XXXX rather than by
        gene symbol.
""")


if __name__ == "__main__":
    main()

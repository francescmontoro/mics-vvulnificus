#!/usr/bin/env python3
"""
15_block_level_fdr.py
=====================
Re-applies multiple-testing correction at the level of LINKAGE BLOCKS rather
than individual variants.

WHY THIS IS THE CORRECT DENOMINATOR
-----------------------------------
With n=29 genomes, variants frequently share an IDENTICAL presence/absence
pattern across strains. Such variants are not independent hypotheses: they
produce byte-for-byte identical test statistics and p-values. Correcting for
19,525 variants when only 2,614 distinct patterns exist inflates the
correction ~7.5-fold and is over-conservative.

pyseer itself documents this: the recommended significance threshold is based
on the number of unique patterns (that is why --output-patterns exists), not
on the number of variants. Collapsing to patterns is also exactly what Fran's
review asked for -- report BLOCKS, not genes.

This script:
  1. reads pyseer's per-variant output and its matching patterns file
  2. groups variants by pattern hash -> linkage blocks
  3. keeps one representative test per block (they are identical anyway)
  4. applies Benjamini-Hochberg FDR across BLOCKS, separately for coding and
     non-coding regions (pooling the two raises the threshold and buries the
     coding signal)
  5. also reports the Bonferroni threshold on patterns, for reference

Note on the coding/non-coding split: a block can contain variants from both
regions. Such blocks are labelled "mixed" and tested in their own stratum, so
no variant is counted twice.

Usage:
  python3 15_block_level_fdr.py \
      --pyseer-dir results/pyseer \
      --var-annot  results/summary/variants_variants.tsv \
      --outdir     results/pyseer_blocks \
      --alpha      0.05
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def load_patterns(patterns_file, n_expected):
    """
    pyseer --output-patterns writes one pattern hash per TESTED variant, in
    the same order as the results file. Returns a list of hashes.
    """
    with open(patterns_file) as fh:
        pats = [ln.strip() for ln in fh if ln.strip()]
    if len(pats) != n_expected:
        print(f"    WARNING: {len(pats)} patterns vs {n_expected} variants; "
              f"cannot match by order, block collapse skipped", file=sys.stderr)
        return None
    return pats


def region_of_block(regions):
    """coding / noncoding / mixed for a set of region labels."""
    s = set(regions) - {""}
    if not s:
        return "unknown"
    if len(s) == 1:
        return s.pop()
    return "mixed"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pyseer-dir", required=True)
    ap.add_argument("--var-annot", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pdir = Path(a.pyseer_dir)

    annot = pd.read_csv(a.var_annot, sep="\t")
    region = dict(zip(annot["variant_id"], annot["region"]))
    gene = dict(zip(annot["variant_id"], annot.get("gene", pd.Series(dtype=str))))
    locus = dict(zip(annot["variant_id"],
                     annot.get("yj016_locus", pd.Series(dtype=str))))
    product = dict(zip(annot["variant_id"],
                       annot.get("product", pd.Series(dtype=str))))
    effect = dict(zip(annot["variant_id"],
                      annot.get("effect", pd.Series(dtype=str))))

    print("=" * 76)
    print("  Block-level multiple-testing correction")
    print("=" * 76)

    summary_rows = []

    for res_file in sorted(pdir.glob("pyseer_*_*.tsv")):
        name = res_file.stem                      # pyseer_GEN2_binary
        if name.endswith("_fdr"):
            continue
        parts = name.split("_")
        if len(parts) < 3:
            continue
        ab, mode = parts[1], parts[2]
        pat_file = pdir / f"patterns_{ab}_{mode}.txt"

        res = pd.read_csv(res_file, sep="\t")
        pcol = "lrt-pvalue" if "lrt-pvalue" in res.columns else res.columns[3]
        res = res.dropna(subset=[pcol]).reset_index(drop=True)

        print(f"\n  ── {ab} ({mode}) ──")
        print(f"    variants tested : {len(res)}")

        if not pat_file.exists():
            print(f"    patterns file missing: {pat_file}")
            continue
        pats = load_patterns(pat_file, len(res))
        if pats is None:
            continue
        res["pattern"] = pats
        res["region"] = res["variant"].map(region).fillna("unknown")

        n_pat = res["pattern"].nunique()
        print(f"    unique patterns : {n_pat}"
              f"   ({100.0*n_pat/len(res):.1f}% of variants)")
        print(f"    Bonferroni threshold on patterns : "
              f"{a.alpha/n_pat:.3e}")

        # ── collapse to blocks ──────────────────────────────────────────────
        blocks = []
        for pat, sub in res.groupby("pattern"):
            # all members share the same p-value; take the first
            r0 = sub.iloc[0]
            blocks.append({
                "block_id": pat[:12],
                "n_variants": len(sub),
                "pvalue": float(r0[pcol]),
                "beta": float(r0.get("beta", np.nan)),
                "variant_h2": float(r0.get("variant_h2", np.nan)),
                "af": float(r0.get("af", np.nan)),
                "region": region_of_block(sub["region"]),
                "representative": r0["variant"],
                "variants": ";".join(sub["variant"].tolist()[:20]),
                "genes": ";".join(sorted({str(gene.get(v, "")) for v in sub["variant"]
                                          if str(gene.get(v, ""))})),
                "yj016_loci": ";".join(sorted({str(locus.get(v, "")) for v in sub["variant"]
                                               if str(locus.get(v, ""))})),
                "effects": ";".join(sorted({str(effect.get(v, "")) for v in sub["variant"]
                                            if str(effect.get(v, ""))})),
                "products": " | ".join(sorted({str(product.get(v, "")) for v in sub["variant"]
                                               if str(product.get(v, ""))}))[:200],
            })
        bdf = pd.DataFrame(blocks)

        # ── BH-FDR per region stratum ───────────────────────────────────────
        out_parts = []
        for reg, sub in bdf.groupby("region"):
            sub = sub.sort_values("pvalue").copy()
            rej, q, _, _ = multipletests(sub["pvalue"], alpha=a.alpha,
                                         method="fdr_bh")
            sub["fdr_bh"] = q
            sub["significant"] = rej
            sub["n_blocks_in_stratum"] = len(sub)
            out_parts.append(sub)
            n_sig = int(rej.sum())
            print(f"      {reg:<10} blocks={len(sub):>5}  "
                  f"min p={sub['pvalue'].min():.3e}  significant={n_sig}")

        final = pd.concat(out_parts).sort_values("pvalue")
        outp = outdir / f"blocks_{ab}_{mode}.tsv"
        final.to_csv(outp, sep="\t", index=False)

        n_sig_total = int(final["significant"].sum())
        print(f"    TOTAL significant blocks : {n_sig_total}")

        if n_sig_total:
            print(f"\n    Significant blocks:")
            for _, r in final[final["significant"]].iterrows():
                print(f"      {r['representative']}  p={r['pvalue']:.2e}  "
                      f"q={r['fdr_bh']:.3f}  n_var={r['n_variants']}  "
                      f"{r['region']}")
                if r["genes"]:
                    print(f"        genes: {r['genes'][:100]}")
        else:
            # report how far the best block was from significance
            best = final.iloc[0]
            stratum_n = int(best["n_blocks_in_stratum"])
            thr = a.alpha / stratum_n
            print(f"    best block: {best['representative']}  "
                  f"p={best['pvalue']:.2e}  (needs <= {thr:.2e} for k=1)")

        summary_rows.append({
            "antibiotic": ab, "mode": mode,
            "variants": len(res), "patterns": n_pat,
            "significant_blocks": n_sig_total,
            "min_pvalue": float(final["pvalue"].min()),
            "bonferroni_threshold": a.alpha / n_pat,
        })

    if summary_rows:
        s = pd.DataFrame(summary_rows)
        s.to_csv(outdir / "block_fdr_summary.tsv", sep="\t", index=False)
        print("\n" + "=" * 76)
        print(s.to_string(index=False))
        print(f"\n  wrote {outdir/'block_fdr_summary.tsv'}")

    print("""
  INTERPRETATION
    Blocks, not variants, are the unit of inference: variants sharing a
    pattern give identical statistics and must not be counted as independent
    findings. Report the number of significant BLOCKS in the manuscript, and
    state the number of unique patterns as the effective number of tests.

    If nothing survives, that is a legitimate result for n=29 with 13 and 7
    resistant strains: the honest conclusion is that no variant reaches
    genome-wide significance after proper correction for population
    structure, and the apparent signal of earlier analyses was attributable
    to lineage effects.
""")


if __name__ == "__main__":
    main()

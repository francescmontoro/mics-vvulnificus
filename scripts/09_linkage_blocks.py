#!/usr/bin/env python3
"""
09_linkage_blocks.py
====================
Groups variants that share an IDENTICAL presence/absence pattern across
strains into linkage blocks, then collapses association results to the block
level.

This directly answers Fran's core criticism (comment #55): the old pipeline
reported "68 GEN2 genes, all with CP-DP=15 and OR=inf" — but those are not 68
independent findings. They are ONE co-segregating block counted 68 times. In
a collection of n=29 with strong clonal structure, many variants ride on the
exact same set of strains and carry exactly the same statistical signal; they
are statistically indistinguishable.

What this script does:
  1. Load the variant presence/absence matrix.
  2. Hash each variant's 0/1 pattern across strains; variants with the same
     hash form a linkage block.
  3. Report, per block: size (n variants), the strains carrying it, the genes
     / loci involved, region composition, and (if given) the best association
     statistic from pyseer and/or Scoary among its members.
  4. Produce a block-level table for the manuscript (report N BLOCKS, not N
     genes, in Abstract / Results / Table 1), plus the variant->block map.

Usage:
    python3 09_linkage_blocks.py \
        --var-matrix results/summary/variants_matrix.tsv \
        --var-annot  results/summary/variants_variants.tsv \
        --pyseer-dir results/pyseer \
        --scoary-dir results/scoary \
        --outdir     results/blocks

Only --var-matrix and --var-annot are required; association dirs are optional
(the block structure itself is the key deliverable).

Requirements: pandas, numpy.
"""

import argparse
import glob
import hashlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def pattern_hash(row_values):
    """Stable hash of a 0/1 pattern (as bytes)."""
    b = bytes(int(v) & 1 for v in row_values)
    return hashlib.md5(b).hexdigest()[:12]


def load_association(pyseer_dir, scoary_dir):
    """
    Optionally load best association p-values per variant from pyseer and/or
    Scoary outputs. Returns dict variant_id -> {source: pvalue}. Best-effort:
    if files are missing or columns differ, we skip gracefully.
    """
    assoc = {}

    # pyseer: files pyseer_<AB>_<mode>_fdr.tsv with columns variant, lrt-pvalue, fdr_bh
    if pyseer_dir and os.path.isdir(pyseer_dir):
        for f in glob.glob(os.path.join(pyseer_dir, "pyseer_*_fdr.tsv")):
            try:
                df = pd.read_csv(f, sep="\t")
            except Exception:
                continue
            vcol = "variant" if "variant" in df.columns else df.columns[0]
            pcol = "fdr_bh" if "fdr_bh" in df.columns else (
                "lrt-pvalue" if "lrt-pvalue" in df.columns else None)
            if pcol is None:
                continue
            tag = "pyseer:" + os.path.basename(f).replace("pyseer_", "").replace("_fdr.tsv", "")
            for _, r in df.iterrows():
                assoc.setdefault(r[vcol], {})[tag] = r[pcol]

    # Scoary: files <TRAIT>.results.csv with a 'Gene' column = our variant_id
    if scoary_dir and os.path.isdir(scoary_dir):
        for f in glob.glob(os.path.join(scoary_dir, "*.results.csv")):
            try:
                df = pd.read_csv(f)
            except Exception:
                continue
            if "Gene" not in df.columns:
                continue
            pcol = None
            for cand in ("Benjamini_H_p", "Bonferroni_p", "Empirical_p",
                         "Best_pairwise_comp_p", "p_value"):
                if cand in df.columns:
                    pcol = cand
                    break
            if pcol is None:
                continue
            trait = os.path.basename(f).split(".")[0]
            for _, r in df.iterrows():
                assoc.setdefault(r["Gene"], {})[f"scoary:{trait}:{pcol}"] = r[pcol]

    return assoc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--var-matrix", required=True)
    ap.add_argument("--var-annot", required=True)
    ap.add_argument("--pyseer-dir", default="")
    ap.add_argument("--scoary-dir", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-block-size", type=int, default=1,
                    help="only report blocks with >= this many variants")
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 68)
    print("  Linkage block analysis (collapse co-segregating variants)")
    print("=" * 68)

    mat = pd.read_csv(a.var_matrix, sep="\t", index_col=0)
    annot = pd.read_csv(a.var_annot, sep="\t").set_index("variant_id")
    strains = list(mat.columns)
    print(f"\n  variants: {mat.shape[0]}  strains: {mat.shape[1]}")

    # 1. Assign each variant to a block by identical pattern
    block_of = {}
    members = {}
    for vid, row in mat.iterrows():
        h = pattern_hash(row.values)
        block_of[vid] = h
        members.setdefault(h, []).append(vid)

    n_blocks = len(members)
    print(f"  linkage blocks (unique patterns): {n_blocks}")
    print(f"  compression: {mat.shape[0]} variants -> {n_blocks} blocks")

    # 2. Optional association
    assoc = load_association(a.pyseer_dir, a.scoary_dir)
    if assoc:
        print(f"  loaded association stats for {len(assoc)} variants")

    # 3. Build block-level table
    assoc_tags = sorted({t for d in assoc.values() for t in d}) if assoc else []
    rows = []
    for h, vids in members.items():
        if len(vids) < a.min_block_size:
            continue
        carriers = [s for s in strains if int(mat.loc[vids[0], s]) == 1]
        # collect gene / locus / region info across members
        genes, loci, regions, effects, products = set(), set(), set(), set(), set()
        for v in vids:
            if v in annot.index:
                an = annot.loc[v]
                if an.get("gene"): genes.add(str(an["gene"]))
                if an.get("yj016_locus"): loci.add(str(an["yj016_locus"]))
                if an.get("region"): regions.add(str(an["region"]))
                if an.get("effect"): effects.add(str(an["effect"]))
                if an.get("product"): products.add(str(an["product"]))
        # best association p per tag among members
        best = {}
        for tag in assoc_tags:
            vals = [assoc[v][tag] for v in vids if v in assoc and tag in assoc[v]]
            vals = [x for x in vals if pd.notna(x)]
            if vals:
                best[tag] = min(vals)
        rows.append({
            "block_id": h,
            "n_variants": len(vids),
            "n_carriers": len(carriers),
            "carriers": ";".join(carriers),
            "n_genes": len(genes),
            "genes": ";".join(sorted(genes)) if genes else "",
            "yj016_loci": ";".join(sorted(loci)) if loci else "",
            "regions": ";".join(sorted(regions)),
            "effects": ";".join(sorted(effects)),
            "products": " | ".join(sorted(products))[:300],
            **{f"best_{t}": best.get(t, np.nan) for t in assoc_tags},
        })

    blocks = pd.DataFrame(rows).sort_values(
        ["n_carriers", "n_variants"], ascending=False)
    blocks_out = outdir / "linkage_blocks.tsv"
    blocks.to_csv(blocks_out, sep="\t", index=False)

    # 4. variant -> block map
    vmap = pd.DataFrame(
        {"variant_id": list(block_of), "block_id": list(block_of.values())})
    vmap = vmap.merge(
        annot.reset_index()[["variant_id", "region", "gene", "yj016_locus"]],
        on="variant_id", how="left")
    vmap_out = outdir / "variant_to_block.tsv"
    vmap.to_csv(vmap_out, sep="\t", index=False)

    # 5. Summary stats for the manuscript
    multi = blocks[blocks["n_variants"] >= 2]
    print("\n  ── Block size distribution ──")
    print(f"    singleton blocks (1 variant) : {(blocks['n_variants']==1).sum()}")
    print(f"    multi-variant blocks (>=2)   : {len(multi)}")
    if len(multi):
        big = blocks.iloc[0]
        print(f"    largest block: {int(big['n_variants'])} variants, "
              f"{int(big['n_carriers'])} carriers")
        print(f"      genes: {big['genes'][:80]}")
        print(f"      -> this is the kind of block the old method counted "
              f"as {int(big['n_variants'])} separate hits")

    print(f"\n  wrote {blocks_out}")
    print(f"  wrote {vmap_out}")
    print("""
  FOR THE MANUSCRIPT:
    Report the number of independent linkage BLOCKS, not the number of
    variants or genes, in the Abstract / Results / Table 1. Within a block,
    all variants share the same presence/absence pattern and therefore the
    same association statistic; they are statistically indistinguishable and
    must not be counted as independent findings (Fran's comment #55).
""")


if __name__ == "__main__":
    main()

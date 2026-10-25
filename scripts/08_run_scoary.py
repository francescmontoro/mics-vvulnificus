#!/usr/bin/env python3
"""
08_run_scoary.py
================
Complementary association analysis with Scoary (Brynildsrud et al., Genome
Biol 2016). Scoary implements the Read & Nee pairwise comparisons method
CORRECTLY over the phylogeny (using non-overlapping pairs), which is exactly
what the old home-made CP>DP filter failed to do (Fran's comment #55).

We run Scoary as an independent, well-established method to cross-check the
pyseer LMM results. It is the standard tool reviewers expect, so agreement
between pyseer-LMM and Scoary strengthens the conclusions.

Scoary needs:
  * a gene_presence_absence.csv in Roary format (first 14 columns are metadata
    that Scoary expects even if unused; then one column per strain, 0/1)
  * a traits file (strains in first column, traits in first row, 1/0)
  * optionally a newick tree (-n) for the pairwise comparisons — WE PROVIDE IT,
    because that is the whole point of doing this correctly.

Strain names must match between the gene file, the traits file and the tree.

Usage:
    mamba run -n vv_assoc python3 08_run_scoary.py \
        --var-matrix results/summary/variants_matrix.tsv \
        --var-annot  results/summary/variants_variants.tsv \
        --mics       data/mic/MICs_converted.tsv \
        --tree       results/pyseer/core.treefile \
        --outdir     results/scoary \
        --permutations 10000 \
        --env-assoc  vv_assoc

Requirements: scoary (env vv_assoc), pandas.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

# Same mapping as the other scripts
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

REFERENCE = "CECT_5769"
MIN_DELTA = 0.5
ANTIBIOTICS = ["GEN2", "STZ"]

# The 14 metadata columns Roary/Scoary expect before the strain columns
ROARY_META_COLS = [
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
    "No. sequences", "Avg sequences per isolate", "Genome Fragment",
    "Order within Fragment", "Accessory Fragment",
    "Accessory Order with Fragment", "QC", "Min group size nuc",
    "Max group size nuc", "Avg group size nuc",
]


def safe(name: str) -> str:
    """Consistent safe strain name (match tree / traits / genes)."""
    return "S_" + str(name).replace("/", "_").replace("-", "_")


def run(cmd, env=None, **kw):
    if env:
        cmd = ["mamba", "run", "-n", env] + cmd
    print("  $ " + " ".join(str(c) for c in cmd), file=sys.stderr)
    return subprocess.run(cmd, check=True, **kw)


def build_roary_csv(var_matrix, var_annot, outdir):
    """
    Convert our variant matrix (variant x strain, 0/1) into a Roary-style
    gene_presence_absence.csv that Scoary can read.
    """
    mat = pd.read_csv(var_matrix, sep="\t", index_col=0)
    annot = pd.read_csv(var_annot, sep="\t").set_index("variant_id")

    strains_safe = [safe(c) for c in mat.columns]

    rows = []
    for vid, row in mat.iterrows():
        an = annot.loc[vid] if vid in annot.index else None
        gene = vid  # unique id per variant
        gene_name = "" if an is None else (an.get("gene") or "")
        annotation = "" if an is None else (
            f"{an.get('region','')}|{an.get('effect','')}|"
            f"{an.get('yj016_locus','')}|{an.get('product','')}")
        n_iso = int(row.sum())
        meta = [gene, gene_name, annotation, n_iso, n_iso, 1.0,
                "", "", "", "", "", "", "", ""]
        presence = [1 if v else 0 for v in row.values]
        rows.append(meta + presence)

    cols = ROARY_META_COLS + strains_safe
    df = pd.DataFrame(rows, columns=cols)
    out = outdir / "gene_presence_absence.csv"
    df.to_csv(out, index=False)
    print(f"  wrote {out} ({len(df)} variants x {len(strains_safe)} strains)",
          file=sys.stderr)
    return out, strains_safe


def build_traits(mics_path, strains_safe, outdir):
    """
    Traits file: rows = strains, columns = GEN2, STZ (binary 1/0).
    Scoary wants strains in first column, traits across the top.
    """
    mics = pd.read_csv(mics_path, sep="\t", index_col=0)
    mics.index = mics.index.astype(str).str.strip()
    ref = mics.loc[REFERENCE]
    safe_to_mic = {safe(snp): mic for snp, mic in SNP_TO_MIC.items()}

    data = {}
    for ab in ANTIBIOTICS:
        col = {}
        for s in strains_safe:
            mic_name = safe_to_mic.get(s)
            if mic_name is None or mic_name not in mics.index:
                continue
            val = mics.loc[mic_name, ab]
            if pd.isna(val):
                continue
            col[s] = 1 if float(val) >= ref[ab] + MIN_DELTA else 0
        data[ab] = col

    traits = pd.DataFrame(data)
    traits.index.name = ""
    out = outdir / "traits.csv"
    traits.to_csv(out)
    for ab in ANTIBIOTICS:
        print(f"  {ab}: {int(traits[ab].sum())} resistant / {traits[ab].notna().sum()} strains",
              file=sys.stderr)
    return out


def prepare_tree(tree_in, strains_safe, outdir):
    """
    Rename the tree tips to the safe names used in the gene/traits files, and
    prune tips that are not present in the variant matrix (e.g. the Snippy
    reference, which has no variants of its own).

    Renaming is done NODE BY NODE with ete3, never by string substitution on
    the newick: strain names here are substrings of one another ('12' inside
    '94-9-112', '32' inside '32'...), so text replacement silently corrupts
    labels (it produced 'S_94_9_1S_12' and 'S_S_A14' in an earlier version).
    """
    if not tree_in:
        return None
    from ete3 import Tree

    t = Tree(str(tree_in), format=1)
    wanted = set(strains_safe)

    renamed, unmatched = 0, []
    for leaf in t.get_leaves():
        raw = leaf.name.strip().strip("'\"")
        cand = safe(raw)
        if cand in wanted:
            leaf.name = cand
            renamed += 1
        else:
            # tip is not in the variant matrix (reference, or name mismatch)
            leaf.name = cand
            unmatched.append(raw)

    # prune to exactly the strains present in the matrix
    keep = [l for l in t.get_leaves() if l.name in wanted]
    if len(keep) < len(wanted):
        missing = wanted - {l.name for l in keep}
        print(f"    WARNING: {len(missing)} strains in the matrix are absent "
              f"from the tree: {sorted(missing)}", file=sys.stderr)
    t.prune(keep, preserve_branch_length=True)

    if unmatched:
        print(f"    pruned {len(unmatched)} tip(s) not in the matrix: "
              f"{unmatched}", file=sys.stderr)
    print(f"    tree: {len(t.get_leaves())} tips after pruning "
          f"({renamed} renamed)", file=sys.stderr)

    out = outdir / "tree_safe.nwk"
    t.write(outfile=str(out), format=1)

    # sanity check: tip set must equal the strain set
    final = {l.name for l in Tree(str(out), format=1).get_leaves()}
    if final != wanted:
        print(f"    ERROR: tip set does not match strain set\n"
              f"      only in tree  : {sorted(final - wanted)}\n"
              f"      only in matrix: {sorted(wanted - final)}", file=sys.stderr)
    else:
        print(f"    tip names match the variant matrix exactly", file=sys.stderr)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--var-matrix", required=True)
    ap.add_argument("--var-annot", required=True)
    ap.add_argument("--mics", required=True)
    ap.add_argument("--tree", default="", help="newick tree (optional but recommended)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--permutations", type=int, default=10000)
    ap.add_argument("--pvalue-cutoff", default="0.05")
    ap.add_argument("--env-assoc", default="vv_assoc")
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 68)
    print("  Scoary association (correct pairwise comparisons over the tree)")
    print("=" * 68)

    print("\n[1] Roary-format gene presence/absence")
    genes_csv, strains_safe = build_roary_csv(a.var_matrix, a.var_annot, outdir)

    print("\n[2] Traits file")
    traits_csv = build_traits(a.mics, strains_safe, outdir)

    print("\n[3] Tree (renamed tips)")
    tree = prepare_tree(a.tree, strains_safe, outdir) if a.tree else None

    print("\n[4] Run Scoary")
    cmd = ["scoary",
           "-g", str(genes_csv),
           "-t", str(traits_csv),
           "-p", str(a.pvalue_cutoff),
           "-c", "BH",                       # Benjamini-Hochberg correction
           "-e", str(a.permutations),        # empirical p by permutation
           "--no-time",
           "-o", str(outdir)]
    if tree:
        cmd += ["-n", str(tree)]
    try:
        run(cmd, env=a.env_assoc)
    except subprocess.CalledProcessError as e:
        print(f"  Scoary failed: {e}", file=sys.stderr)
        sys.exit(1)

    print("\nDone. Scoary results (TRAIT.results.csv) in", outdir)
    print("  Key columns: Bonferroni_p, Benjamini_H_p, "
          "Empirical_p (permutation), Best_pairwise_comp_p")


if __name__ == "__main__":
    main()

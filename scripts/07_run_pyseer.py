#!/usr/bin/env python3
"""
07_run_pyseer.py
================
Association of variants with MIC using pyseer's linear mixed model (LMM),
correcting for population structure via a kinship matrix derived from the
core-genome phylogeny.

This REPLACES the old CP > DP pairwise-concordance filter, which did not
actually use the tree and rewarded lineage confounding (Fran's comment #55).
The LMM is the standard, defensible method: it models population structure
as a random effect and tests each variant as a fixed effect.

Two phenotype modes (run both, compare):
  * continuous : the normalised MIC score [0-5]  (preferred; no MIN_DELTA)
  * binary     : resistant/sensitive at MIN_DELTA (for comparison / Scoary)

STEPS (this script orchestrates them; heavy tools run in their conda envs):
  1. Build the tree from Panaroo's core alignment with IQ-TREE (if not given).
  2. Compute the kinship/similarity matrix from the tree
     (pyseer's phylogeny_distance.py --lmm).
  3. Build the pyseer variant input (Rtab) from build_snp_matrix's matrix.
  4. Build the phenotype files (continuous and binary) with SAFE sample names.
  5. Run pyseer --lmm for GEN2 and STZ.
  6. Apply Benjamini-Hochberg FDR SEPARATELY for coding and noncoding
     variants (Fran: pooling raises the threshold and kills the coding hits).

IMPORTANT (from pyseer docs): purely numeric sample names (12, 32, 162, ...)
break pyseer. We therefore prefix every strain with 'S_' consistently across
the kinship matrix, the variant Rtab and the phenotype file.

Usage (after run_vv_variants.sh has produced the Panaroo alignment & matrix):
    mamba run -n vv_assoc python3 07_run_pyseer.py \
        --core-aln   results/panaroo/core_gene_alignment.aln \
        --var-matrix results/summary/variants_matrix.tsv \
        --var-annot  results/summary/variants_variants.tsv \
        --mics       data/mic/MICs_converted.tsv \
        --outdir     results/pyseer \
        --env-core   vv_core \
        --threads    12

Requirements: pyseer (env vv_assoc), iqtree (env vv_core), pandas, numpy.
"""

import os
for _v in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS",
           "MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Strain name mapping (SNP/tree names -> MIC-table names), same as signal check
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


def safe(name: str) -> str:
    """Prefix to avoid purely-numeric sample names that break pyseer."""
    return "S_" + str(name).replace("/", "_").replace("-", "_")


def run(cmd, env=None, **kw):
    """Run a shell command, optionally inside a conda env via mamba run."""
    if env:
        cmd = ["mamba", "run", "-n", env] + cmd
    print("  $ " + " ".join(cmd), file=sys.stderr)
    return subprocess.run(cmd, check=True, **kw)


def build_tree(core_aln, outdir, env_core, threads):
    """Step 1: IQ-TREE from the core alignment (if tree not present)."""
    tree = outdir / "core.treefile"
    if tree.exists():
        print(f"  tree exists: {tree}", file=sys.stderr)
        return tree
    prefix = outdir / "core"
    run(["iqtree3", "-s", str(core_aln), "-m", "GTR+ASC+G",
         "-B", "1000", "-T", "AUTO", "-ntmax", str(threads),
         "--prefix", str(prefix)], env=env_core)
    return tree


def build_kinship(tree, outdir):
    """
    Step 2: kinship / similarity matrix from the tree.

    pyseer ships a helper (phylogeny_distance.py --lmm) but bioconda does not
    package it, so we implement it here. The LMM similarity between two tips
    is the SHARED branch length from the root down to their most recent common
    ancestor: tips that share more evolutionary history are more similar, and
    that is exactly the population structure the LMM must absorb as a random
    effect.

        K[i][j] = distance(root -> MRCA(i, j))
        K[i][i] = distance(root -> i)

    The matrix is written with the SAME safe sample names used for the variant
    Rtab and the phenotype files (pyseer requires them to match, and purely
    numeric names such as 12 / 32 / 162 break it).
    """
    from Bio import Phylo

    K = outdir / "phylogeny_K.tsv"
    if K.exists():
        print(f"  kinship exists: {K}", file=sys.stderr)
        return K

    t = Phylo.read(str(tree), "newick")
    tips = t.get_terminals()
    names = [tip.name for tip in tips]

    # depth (cumulative branch length) from the root for every clade
    depths = t.depths(unit_branch_lengths=False)

    n = len(tips)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            if i == j:
                shared = depths.get(tips[i], 0.0)
            else:
                mrca = t.common_ancestor(tips[i], tips[j])
                shared = depths.get(mrca, 0.0)
            mat[i, j] = shared
            mat[j, i] = shared

    safe_names = [safe(nm) for nm in names]
    df = pd.DataFrame(mat, index=safe_names, columns=safe_names)
    df.to_csv(K, sep="\t")
    print(f"  wrote kinship {K} ({n}x{n})", file=sys.stderr)
    print(f"    mean shared branch length: {mat[np.triu_indices(n,1)].mean():.4f}",
          file=sys.stderr)
    return K


def build_rtab(var_matrix, outdir):
    """
    Step 3: convert build_snp_matrix's matrix (variant x strain, 0/1) into a
    pyseer 'Rtab' presence/absence file (variant  strain1 strain2 ...), with
    SAFE strain names.
    """
    rtab = outdir / "variants.Rtab"
    df = pd.read_csv(var_matrix, sep="\t", index_col=0)
    # columns are strain (SNP) names -> safe names
    df.columns = [safe(c) for c in df.columns]
    df.index.name = "Variant"
    df.to_csv(rtab, sep="\t")
    print(f"  wrote {rtab} ({df.shape[0]} variants x {df.shape[1]} strains)",
          file=sys.stderr)
    return rtab, list(df.columns)


def build_phenotypes(mics_path, safe_strains, outdir):
    """
    Step 4: phenotype files (continuous MIC and binary R/S) with safe names.
    Returns dict {ab: {'cont': path, 'bin': path}}.
    """
    mics = pd.read_csv(mics_path, sep="\t", index_col=0)
    mics.index = mics.index.astype(str).str.strip()
    ref = mics.loc[REFERENCE]

    # inverse map: safe -> MIC name, via SNP name
    snp_to_safe = {snp: safe(snp) for snp in SNP_TO_MIC}
    safe_to_mic = {safe(snp): mic for snp, mic in SNP_TO_MIC.items()}

    out = {}
    for ab in ANTIBIOTICS:
        cont = outdir / f"pheno_{ab}_continuous.tsv"
        binp = outdir / f"pheno_{ab}_binary.tsv"
        rows_c, rows_b = [], []
        for s in safe_strains:
            mic_name = safe_to_mic.get(s)
            if mic_name is None or mic_name not in mics.index:
                continue
            val = mics.loc[mic_name, ab]
            if pd.isna(val):
                continue
            rows_c.append((s, float(val)))
            is_r = 1 if float(val) >= ref[ab] + MIN_DELTA else 0
            rows_b.append((s, is_r))
        pd.DataFrame(rows_c, columns=["samples", ab]).to_csv(
            cont, sep="\t", index=False)
        pd.DataFrame(rows_b, columns=["samples", ab]).to_csv(
            binp, sep="\t", index=False)
        n_r = sum(v for _, v in rows_b)
        print(f"  {ab}: {len(rows_c)} strains, {n_r} resistant (binary)",
              file=sys.stderr)
        out[ab] = {"cont": cont, "bin": binp}
    return out


def run_pyseer(rtab, K, pheno, ab, continuous, outdir, env_assoc, threads):
    """Step 5: run pyseer --lmm for one antibiotic / one phenotype mode."""
    mode = "continuous" if continuous else "binary"
    out = outdir / f"pyseer_{ab}_{mode}.tsv"
    patt = outdir / f"patterns_{ab}_{mode}.txt"
    cmd = ["pyseer", "--lmm",
           "--phenotypes", str(pheno),
           "--phenotype-column", ab,
           "--pres", str(rtab),
           "--similarity", str(K),
           "--output-patterns", str(patt),
           "--cpu", "1"]
    if continuous:
        cmd.append("--continuous")
    with open(out, "w") as fh:
        run(cmd, env=env_assoc, stdout=fh)
    print(f"  wrote {out}", file=sys.stderr)
    return out, patt


def apply_fdr_split(pyseer_out, var_annot, ab, mode, outdir):
    """
    Step 6: BH-FDR applied SEPARATELY to coding and noncoding variants.
    (Fran: pooling coding+noncoding raises the threshold and kills coding hits.)
    """
    from statsmodels.stats.multitest import multipletests

    res = pd.read_csv(pyseer_out, sep="\t")
    # pyseer's p-value column is 'lrt-pvalue' for --lmm
    pcol = "lrt-pvalue" if "lrt-pvalue" in res.columns else res.columns[-1]

    # map variant -> region via annotation
    annot = pd.read_csv(var_annot, sep="\t")
    # variant ids in pyseer come from the Rtab index (our variant_id)
    region = dict(zip(annot["variant_id"], annot["region"]))
    res["region"] = res["variant"].map(region).fillna("unknown")

    out_rows = []
    for reg, sub in res.groupby("region"):
        sub = sub.dropna(subset=[pcol]).copy()
        if sub.empty:
            continue
        rej, q, _, _ = multipletests(sub[pcol], method="fdr_bh")
        sub["fdr_bh"] = q
        sub["significant"] = rej
        out_rows.append(sub)
    if not out_rows:
        print(f"  {ab} {mode}: no testable variants", file=sys.stderr)
        return None
    final = pd.concat(out_rows).sort_values(pcol)
    outp = outdir / f"pyseer_{ab}_{mode}_fdr.tsv"
    final.to_csv(outp, sep="\t", index=False)
    n_sig = int(final["significant"].sum())
    print(f"  {ab} {mode}: {n_sig} significant after per-region FDR", file=sys.stderr)
    return outp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core-aln", required=True)
    ap.add_argument("--var-matrix", required=True)
    ap.add_argument("--var-annot", required=True)
    ap.add_argument("--mics", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tree", default="", help="precomputed tree (optional)")
    ap.add_argument("--env-core", default="vv_core")
    ap.add_argument("--env-assoc", default="vv_assoc")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 68)
    print("  pyseer LMM association (replaces the CP>DP filter)")
    print("=" * 68)

    # 1. tree
    print("\n[1] Tree from core alignment")
    tree = Path(a.tree) if a.tree else build_tree(
        Path(a.core_aln), outdir, a.env_core, a.threads)

    # 2. kinship
    print("\n[2] Kinship matrix from tree")
    K = build_kinship(tree, outdir)

    # 3. Rtab
    print("\n[3] Variant Rtab")
    rtab, safe_strains = build_rtab(a.var_matrix, outdir)

    # 4. phenotypes
    print("\n[4] Phenotype files (continuous + binary)")
    phenos = build_phenotypes(a.mics, safe_strains, outdir)

    # 5-6. pyseer + FDR per antibiotic and mode
    for ab in ANTIBIOTICS:
        for continuous in (True, False):
            mode = "continuous" if continuous else "binary"
            print(f"\n[5] pyseer LMM: {ab} ({mode})")
            pheno = phenos[ab]["cont" if continuous else "bin"]
            try:
                out, _ = run_pyseer(rtab, K, pheno, ab, continuous,
                                    outdir, a.env_assoc, a.threads)
            except subprocess.CalledProcessError as e:
                print(f"  pyseer failed for {ab} {mode}: {e}", file=sys.stderr)
                continue
            print(f"[6] Per-region BH-FDR: {ab} ({mode})")
            apply_fdr_split(out, a.var_annot, ab, mode, outdir)

    print("\nDone. Results in", outdir)


if __name__ == "__main__":
    main()

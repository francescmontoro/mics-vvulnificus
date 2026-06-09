#!/usr/bin/env python3
"""
09_phylo_correction.py
----------------------
Phylogenetic correction for SNP-MIC associations.
Implements Scoary's pairwise comparisons method (Kurtz et al. 2016)
using vectorized numpy operations for efficiency.

For each (SNP, antibiotic) pair significant in Fisher test, counts:
  CP (concordant pairs): pairs where SNP and resistance co-change
  DP (discordant pairs): pairs where they change independently

A SNP is phylogenetically robust if CP > DP.

Optimized to run efficiently on 30 strains without crashing WSL.

Usage:
    python3 09_phylo_correction.py

Inputs:
    ../data/processed/SNP_significant.tsv
    ../data/processed/MICs_converted.tsv
    ../data/raw/phylo_tree_renamed.nwk
    ../data/raw/allnsSNP.tsv

Outputs:
    ../data/processed/SNP_phylo_corrected.tsv
    ../results/gene_level_phylo_corrected.tsv

Author: Francesc Montoro (2025)
"""

import pandas as pd
import numpy as np
from Bio import Phylo
from pathlib import Path
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
RAW  = BASE / "data" / "raw"
PROC = BASE / "data" / "processed"
RES  = BASE / "results"
RES.mkdir(exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
REFERENCE = "CECT_5769"
MIN_DELTA = 0.5

SNP_TO_MIC = {
    "11028": "BT3_11028", "12": "BT3_12", "162": "BT3_162", "32": "BT3_32",
    "94385": "94-3-85",   "960426-1_4C": "960426_1/4C",
    "CIP8190": "C1P8190", "V4": "V_4",
    "CECT4604": "CECT_4604", "CECT4606": "CECT_4606", "CECT4608": "CECT_4608",
    "CECT4865": "CECT_4865", "CECT4866": "CECT_4866", "CECT4999": "CECT_4999",
    "CECT5198": "CECT_5198", "CECT5763": "CECT_5763", "CECT5769": "CECT_5769",
    "CECT7030": "CECT_7030", "CECT898":  "CECT_898",  "CG100":    "CG_100",
    "94-9-112": "94-9-112", "94-9-119": "94-9-119",   "95-8-161": "95-8-161",
    "95-8-6":   "95-8-6",   "95-8-7":   "95-8-7",     "A14":      "A14",
    "PD-2-51":  "PD-2-51",  "Ra_3":     "Ra_3",        "RIU1":     "RIU1",
    "YJ016":    "YJ016",
}


def compute_cp_dp_vectorized(snp_matrix: np.ndarray,
                              resistant_vec: np.ndarray) -> tuple:
    """
    Vectorized CP/DP computation for all SNPs simultaneously.

    snp_matrix: (n_snps, n_strains) binary int8 array
    resistant_vec: (n_strains,) binary int8 array

    Returns cp_array, dp_array each of shape (n_snps,)
    """
    n_snps, n_strains = snp_matrix.shape
    pair_indices = list(combinations(range(n_strains), 2))

    if not pair_indices:
        return np.zeros(n_snps, dtype=int), np.zeros(n_snps, dtype=int)

    idx_a = np.array([p[0] for p in pair_indices])
    idx_b = np.array([p[1] for p in pair_indices])

    # SNP differences: (n_snps, n_pairs)
    snp_diff = snp_matrix[:, idx_a] != snp_matrix[:, idx_b]

    # Resistance differences: (n_pairs,)
    res_diff = resistant_vec[idx_a] != resistant_vec[idx_b]

    # CP: both differ (concordant)
    cp_array = (snp_diff & res_diff[np.newaxis, :]).sum(axis=1)

    # DP: only one differs (discordant)
    dp_array = (snp_diff & ~res_diff[np.newaxis, :]).sum(axis=1) + \
               (~snp_diff & res_diff[np.newaxis, :]).sum(axis=1)

    return cp_array.astype(int), dp_array.astype(int)


if __name__ == "__main__":
    print("=" * 65)
    print("  Phylogenetic correction — Pairwise Comparisons (vectorized)")
    print("=" * 65)

    # Load tree
    tree_file = RAW / "phylo_tree_renamed.nwk"
    print(f"\nLoading tree: {tree_file}")
    tree = Phylo.read(str(tree_file), "newick")
    tree_strains   = [tip.name for tip in tree.get_terminals()]
    genome_strains = [s for s in tree_strains if s != REFERENCE]
    n_strains      = len(genome_strains)
    strain_idx     = {s: i for i, s in enumerate(genome_strains)}
    print(f"  Strains: {n_strains} (excl. reference)")

    # Load MICs
    mics = pd.read_csv(PROC / "MICs_converted.tsv", sep="\t", index_col=0)
    mics.index = mics.index.astype(str).str.strip()
    ref = mics.loc[REFERENCE]

    # Load significant SNPs
    sig = pd.read_csv(PROC / "SNP_significant.tsv", sep="\t")
    print(f"Significant pairs: {len(sig):,} | "
          f"Unique SNPs: {sig['snp_id'].nunique():,} | "
          f"Antibiotics: {sig['antibiotic'].nunique()}")

    # Load raw SNPs (only significant)
    print("\nLoading raw SNPs...")
    sig_snp_ids = set(sig["snp_id"].unique())
    snp_raw = pd.read_csv(RAW / "allnsSNP.tsv", sep="\t")
    snp_raw.columns = snp_raw.columns.str.strip()
    snp_raw["Strain"] = snp_raw["Strain"].astype(str).str.strip()
    snp_raw["strain_mic"] = snp_raw["Strain"].map(SNP_TO_MIC)
    snp_raw = snp_raw[snp_raw["strain_mic"].notna()].copy()
    snp_raw["snp_id"] = (
        snp_raw["Gene"].str.strip() + "__" +
        snp_raw["Position"].astype(str) + "__" +
        snp_raw["Reference"].str.strip() + "__" +
        snp_raw["SNP"].str.strip()
    )
    snp_raw = snp_raw[snp_raw["snp_id"].isin(sig_snp_ids)]

    # Build SNP presence matrix
    print("Building SNP matrix...")
    unique_snps = sig["snp_id"].unique()
    snp_to_idx  = {s: i for i, s in enumerate(unique_snps)}
    n_snps      = len(unique_snps)

    snp_matrix = np.zeros((n_snps, n_strains), dtype=np.int8)
    for _, row in snp_raw.iterrows():
        si = snp_to_idx.get(row["snp_id"])
        st = strain_idx.get(row["strain_mic"])
        if si is not None and st is not None:
            snp_matrix[si, st] = 1

    print(f"  SNP matrix: {snp_matrix.shape[0]:,} × {snp_matrix.shape[1]}")
    print(f"  Memory: {snp_matrix.nbytes / 1e6:.1f} MB")

    # Run per antibiotic
    all_records = []
    antibiotics = sig["antibiotic"].unique()

    for ab_idx, ab in enumerate(antibiotics, 1):
        sig_ab = sig[sig["antibiotic"] == ab]
        snp_ids_ab     = sig_ab["snp_id"].values
        ab_snp_indices = np.array([snp_to_idx[s] for s in snp_ids_ab])

        resistant_vec = np.array([
            int(s in mics.index and
                not pd.isna(mics.loc[s, ab]) and
                mics.loc[s, ab] >= ref[ab] + MIN_DELTA)
            for s in genome_strains
        ], dtype=np.int8)

        n_res = resistant_vec.sum()
        print(f"  [{ab_idx:02d}/{len(antibiotics)}] {ab}: "
              f"{n_res} resistant, {len(snp_ids_ab)} SNPs...", end="\r")

        ab_matrix = snp_matrix[ab_snp_indices, :]
        cp_arr, dp_arr = compute_cp_dp_vectorized(ab_matrix, resistant_vec)

        for i, (_, row) in enumerate(sig_ab.iterrows()):
            all_records.append({
                "antibiotic":   ab,
                "snp_id":       row["snp_id"],
                "gene":         row["gene"],
                "position":     row["position"],
                "Aa_ref":       row.get("Aa_ref", ""),
                "Aa_mod":       row.get("Aa_mod", ""),
                "fisher_fdr":   row["fdr"],
                "fisher_OR":    row["odds_ratio"],
                "CP":           int(cp_arr[i]),
                "DP":           int(dp_arr[i]),
                "CP_minus_DP":  int(cp_arr[i]) - int(dp_arr[i]),
                "phylo_robust": bool(cp_arr[i] > dp_arr[i]),
            })

        print(f"  [{ab_idx:02d}/{len(antibiotics)}] {ab}: done.              ")

    print(f"\nTotal pairs processed: {len(all_records):,}")
    results = pd.DataFrame(all_records)

    out_all = PROC / "SNP_phylo_corrected.tsv"
    results.to_csv(out_all, sep="\t", index=False)
    print(f"All results → {out_all}")

    robust = results[results["phylo_robust"]]
    print(f"Phylogenetically robust (CP > DP): {len(robust):,}")

    # Gene-level summary
    gene_records = []
    for (ab, gene), grp in robust.groupby(["antibiotic", "gene"]):
        gene_records.append({
            "antibiotic":      ab,
            "gene":            gene,
            "n_snps":          len(grp),
            "best_fisher_fdr": grp["fisher_fdr"].min(),
            "best_OR":         grp["fisher_OR"].max(),
            "max_CP":          grp["CP"].max(),
            "max_DP":          grp["DP"].max(),
            "max_CP_minus_DP": grp["CP_minus_DP"].max(),
        })

    gene_df = pd.DataFrame(gene_records).sort_values(
        ["antibiotic", "max_CP_minus_DP"], ascending=[True, False]
    )
    out_gene = RES / "gene_level_phylo_corrected.tsv"
    gene_df.to_csv(out_gene, sep="\t", index=False)
    print(f"Gene-level robust → {out_gene}")

    print(f"\n── Summary per antibiotic ──")
    for ab in sorted(robust["antibiotic"].unique()):
        sub   = robust[robust["antibiotic"] == ab]
        sub_g = gene_df[gene_df["antibiotic"] == ab]
        n_genes = sub["gene"].nunique()
        if len(sub_g) > 0:
            top = sub_g.iloc[0]
            print(f"  {ab:8s}: {len(sub):5d} robust SNPs, "
                  f"{n_genes:4d} genes | "
                  f"top: {top['gene'][:40]}  CP-DP={top['max_CP_minus_DP']}")

    print("\nDone.")

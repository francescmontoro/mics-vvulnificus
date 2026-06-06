#!/usr/bin/env python3
"""
04_gene_level_summary.py
------------------------
Two complementary analyses on top of the Fisher's exact test results:

1. GENE-LEVEL COLLAPSE (per antibiotic)
   For each (gene, antibiotic) pair, collapse all significant SNPs into
   a single row reporting:
     - number of significant SNPs in that gene
     - best FDR (minimum)
     - best odds ratio (maximum)
     - amino acid changes observed

2. ANTIBIOTIC FAMILY ANALYSIS
   Re-run Fisher's exact test grouping antibiotics by class.
   A strain is "resistant" if it is resistant to ANY antibiotic in the class.

Usage:
    python3 04_gene_level_summary.py

Inputs:
    ../data/processed/SNP_significant.tsv
    ../data/processed/MICs_converted.tsv
    ../data/raw/allnsSNP.tsv

Outputs:
    ../results/gene_level_per_ab.tsv
    ../results/gene_level_per_family.tsv
    ../data/processed/SNP_family_associations.tsv
    ../data/processed/SNP_family_significant.tsv

Author: Francesc Montoro (2025)
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
RAW  = BASE / "data" / "raw"
PROC = BASE / "data" / "processed"
RES  = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
REFERENCE     = "CECT_5769"
FDR_THRESHOLD = 0.05
MIN_STRAINS   = 2

# ── Antibiotic families ───────────────────────────────────────────────────────
AB_FAMILIES = {
    "Aminoglycoside":  ["AMI", "GEN", "GEN2", "NEO", "TOB", "STR", "SPE"],
    "Beta-lactam":     ["AMP", "FAZ", "FEP", "FOX", "AXO", "TAZ", "POD",
                        "PIP", "AMOX", "PEN", "TIO", "AZT",
                        "A/S2", "P/T4", "TIM2"],
    "Carbapenem":      ["IMI", "MERO"],
    "Fluoroquinolone": ["CIP", "ENRO", "GAT", "NOV"],
    "Macrolide":       ["ERY", "TYLT"],
    "Lincosamide":     ["CLI"],
    "Tetracycline":    ["TET", "OXY"],
    "Sulfonamide":     ["SDM", "STZ", "SXT", "SXT2", "TANS"],
    "Phenicol":        ["FFN"],
    "Other":           ["FUR", "NIT"],
}

# ── Strain name mapping (SNP → MIC) ──────────────────────────────────────────
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


def bh_correction(pvalues: pd.Series) -> pd.Series:
    n = len(pvalues)
    if n == 0:
        return pvalues.copy()
    sorted_idx = pvalues.argsort()
    sorted_p   = pvalues.iloc[sorted_idx].values
    ranks      = np.arange(1, n + 1)
    adjusted   = np.minimum(1.0, sorted_p * n / ranks)
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    result = pd.Series(index=pvalues.index, dtype=float)
    result.iloc[sorted_idx] = adjusted
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1 — Gene-level collapse
# ═══════════════════════════════════════════════════════════════════════════════

def gene_level_collapse(sig: pd.DataFrame,
                        group_col: str = "antibiotic") -> pd.DataFrame:
    """
    Collapse significant SNPs to gene level.
    Works for both per-antibiotic and per-family dataframes.
    Handles both column naming conventions (Aa_ref/Aa_mod and ref_codon/alt_codon).
    """
    # Detect which column names are present
    aa_ref_col = "Aa_ref"    if "Aa_ref"    in sig.columns else "ref_codon"
    aa_mod_col = "Aa_mod"    if "Aa_mod"    in sig.columns else "alt_codon"

    records = []
    for (group, gene), grp in sig.groupby([group_col, "gene"]):
        dedup = grp.drop_duplicates(["position"])
        aa_changes = "; ".join(
            f"{row[aa_ref_col]}{int(row['position'])}{row[aa_mod_col]}"
            for _, row in dedup.iterrows()
        )
        records.append({
            group_col:               group,
            "gene":                  gene,
            "n_snps":                len(grp),
            "best_fdr":              round(grp["fdr"].min(), 6),
            "best_OR":               grp["odds_ratio"].replace([np.inf], 999).max(),
            "best_OR_raw":           grp["odds_ratio"].max(),
            "n_resistant_with_snp":  int(grp["a_snp_res"].max()),
            "n_sensitive_with_snp":  int(grp["b_snp_sen"].min()),
            "aa_changes":            aa_changes[:300],
        })

    df = pd.DataFrame(records)
    if len(df) == 0:
        return df
    return df.sort_values([group_col, "best_fdr", "best_OR"],
                          ascending=[True, True, False])


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2 — Family-level Fisher test
# ═══════════════════════════════════════════════════════════════════════════════

def build_snp_matrix(snps: pd.DataFrame, genome_strains: list) -> pd.DataFrame:
    snps = snps.copy()
    snps["snp_id"] = (
        snps["Gene"].str.strip() + "__" +
        snps["Position"].astype(str) + "__" +
        snps["Reference"].str.strip() + "__" +
        snps["SNP"].str.strip()
    )
    mat = (snps
           .groupby(["snp_id", "strain_mic"])
           .size()
           .unstack(fill_value=0)
           .clip(upper=1))
    for s in genome_strains:
        if s not in mat.columns:
            mat[s] = 0
    return mat[genome_strains]


def run_family_fisher(snp_mat: pd.DataFrame,
                      mics: pd.DataFrame,
                      genome_strains: list) -> pd.DataFrame:
    ref_scores = mics.loc[REFERENCE]
    records    = []
    families   = list(AB_FAMILIES.keys())

    for fam_idx, family in enumerate(families, 1):
        abs_in_fam = [ab for ab in AB_FAMILIES[family] if ab in mics.columns]
        if not abs_in_fam:
            continue

        resistant = set()
        for ab in abs_in_fam:
            ref_score = ref_scores.get(ab, np.nan)
            if pd.isna(ref_score):
                continue
            resistant |= {s for s in genome_strains
                          if s in mics.index and mics.loc[s, ab] > ref_score}
        sensitive = set(genome_strains) - resistant

        n_res = len(resistant)
        n_tot = len(genome_strains)

        print(f"  [{fam_idx:02d}/{len(families)}] {family}: "
              f"{n_res} resistant / {n_tot} strains, "
              f"{len(snp_mat):,} SNPs...", end="\r")

        if n_res < MIN_STRAINS:
            print(f"\n  [SKIP] {family}: only {n_res} resistant strains")
            continue

        for snp_id, row in snp_mat.iterrows():
            with_snp    = set(row[row == 1].index)
            without_snp = set(row[row == 0].index)

            if len(with_snp) < MIN_STRAINS:
                continue

            a = len(with_snp    & resistant)
            b = len(with_snp    & sensitive)
            c = len(without_snp & resistant)
            d = len(without_snp & sensitive)

            if a + b == 0 or c + d == 0:
                continue

            odds_ratio, p_value = fisher_exact([[a, b], [c, d]],
                                               alternative="greater")
            records.append({
                "family":        family,
                "snp_id":        snp_id,
                "gene":          snp_id.split("__")[0],
                "position":      snp_id.split("__")[1],
                "ref_codon":     snp_id.split("__")[2],
                "alt_codon":     snp_id.split("__")[3],
                "n_strains_snp": len(with_snp),
                "n_resistant":   n_res,
                "n_total":       n_tot,
                "a_snp_res":     a,
                "b_snp_sen":     b,
                "c_nosnp_res":   c,
                "d_nosnp_sen":   d,
                "odds_ratio":    round(odds_ratio, 4),
                "p_value":       p_value,
            })

    print()
    return pd.DataFrame(records)


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    # ── Part 1 ────────────────────────────────────────────────────────────────
    print("=" * 65)
    print("  Part 1 — Gene-level collapse per antibiotic")
    print("=" * 65)

    sig = pd.read_csv(PROC / "SNP_significant.tsv", sep="\t")
    print(f"Loaded {len(sig):,} significant SNP associations")

    gene_ab = gene_level_collapse(sig, group_col="antibiotic")
    out1 = RES / "gene_level_per_ab.tsv"
    gene_ab.to_csv(out1, sep="\t", index=False)
    print(f"Gene-level table: {len(gene_ab)} (gene, antibiotic) pairs → {out1}")

    print(f"\n── Top 10 associations by FDR ──")
    top = gene_ab.sort_values(["best_fdr", "best_OR"], ascending=[True, False]).head(10)
    print(top[["antibiotic", "gene", "n_snps", "best_fdr", "best_OR_raw",
               "n_resistant_with_snp", "n_sensitive_with_snp"]].to_string(index=False))

    # ── Part 2 ────────────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("  Part 2 — Antibiotic family analysis")
    print("=" * 65)

    snps_raw = pd.read_csv(RAW / "allnsSNP.tsv", sep="\t")
    snps_raw.columns = snps_raw.columns.str.strip()
    snps_raw["Strain"] = snps_raw["Strain"].astype(str).str.strip()
    snps_raw["strain_mic"] = snps_raw["Strain"].map(SNP_TO_MIC)
    snps_raw = snps_raw[snps_raw["strain_mic"].notna()].copy()

    mics = pd.read_csv(PROC / "MICs_converted.tsv", sep="\t", index_col=0)
    mics.index = mics.index.astype(str).str.strip()

    genome_strains = sorted(snps_raw["strain_mic"].unique())
    snp_mat = build_snp_matrix(snps_raw, genome_strains)
    print(f"SNP matrix: {snp_mat.shape[0]:,} SNPs × {snp_mat.shape[1]} strains")
    print(f"\nRunning family-level Fisher tests...\n")

    fam_results = run_family_fisher(snp_mat, mics, genome_strains)
    print(f"Total family tests: {len(fam_results):,}")

    if len(fam_results) == 0:
        print("No family tests performed.")
    else:
        # BH correction per family
        fdr_chunks = []
        for fam, grp in fam_results.groupby("family"):
            grp = grp.copy()
            grp["fdr"] = bh_correction(grp["p_value"])
            fdr_chunks.append(grp)
        fam_results = pd.concat(fdr_chunks).sort_values(["family", "fdr"])
        fam_results["fdr"] = fam_results["fdr"].round(6)

        out_fam_all = PROC / "SNP_family_associations.tsv"
        fam_results.to_csv(out_fam_all, sep="\t", index=False)
        print(f"All family associations → {out_fam_all}")

        fam_sig = fam_results[(fam_results["fdr"] < FDR_THRESHOLD) &
                               (fam_results["odds_ratio"] > 1)]
        out_fam_sig = PROC / "SNP_family_significant.tsv"
        fam_sig.to_csv(out_fam_sig, sep="\t", index=False)
        print(f"Significant family associations: {len(fam_sig)}")

        if len(fam_sig) > 0:
            gene_fam = gene_level_collapse(fam_sig, group_col="family")
            out2 = RES / "gene_level_per_family.tsv"
            gene_fam.to_csv(out2, sep="\t", index=False)
            print(f"Gene-level family table → {out2}")

            print(f"\n── Family summary ──")
            for fam in sorted(fam_sig["family"].unique()):
                sub  = fam_sig[fam_sig["family"] == fam]
                n_genes   = sub["gene"].nunique()
                by_gene   = sub.groupby("gene")
                top_gene  = by_gene["odds_ratio"].max().sort_values(ascending=False).index[0]
                top_or    = by_gene["odds_ratio"].max().max()
                top_fdr   = by_gene["fdr"].min().min()
                print(f"  {fam:<20}: {len(sub):5d} SNPs  {n_genes:4d} genes | "
                      f"top: {top_gene[:40]}  OR={top_or:.1f}  FDR={top_fdr:.4f}")
        else:
            print("No significant family associations found.")

    print("\nDone.")

#!/usr/bin/env python3
"""
03_snp_association.py
---------------------
Statistical association between non-synonymous SNPs and antibiotic resistance.

For each (SNP, antibiotic) pair we run a Fisher's exact test:

              | MIC increased | MIC not increased |
  SNP present |      a        |        b          |
  SNP absent  |      c        |        d          |

A strain is "MIC increased" if its converted MIC score is at least
MIN_DELTA above the reference strain (CECT_5769). MIN_DELTA=0.5 avoids
spurious associations from sub-scale noise while remaining conservative.

In the absence of species-specific clinical breakpoints for V. vulnificus,
resistance was defined relative to the most susceptible reference strain,
following the approach used for species lacking established breakpoints
(Roig et al., 2009; Jones & Oliver, 2009).

Multiple testing correction: Benjamini-Hochberg FDR within each antibiotic.
Only SNPs with FDR < FDR_THRESHOLD and odds_ratio > 1 are reported.

Usage:
    python3 03_snp_association.py

Inputs:
    ../data/raw/allnsSNP.tsv
    ../data/processed/MICs_converted.tsv

Outputs:
    ../data/processed/SNP_associations.tsv
    ../data/processed/SNP_significant.tsv
    ../results/SNP_association_summary.tsv

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
MIN_DELTA     = 0.5   # minimum score difference above reference to call resistant
                      # justification: avoids sub-scale noise; conservative threshold
                      # consistent with absence of species-specific breakpoints

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


def load_snps(snp_file: Path) -> pd.DataFrame:
    print(f"  Loading SNPs: {snp_file}")
    df = pd.read_csv(snp_file, sep="\t")
    df.columns = df.columns.str.strip()
    df["Strain"] = df["Strain"].astype(str).str.strip()
    df["strain_mic"] = df["Strain"].map(SNP_TO_MIC)
    unmapped = df[df["strain_mic"].isna()]["Strain"].unique()
    if len(unmapped) > 0:
        print(f"  [WARNING] Unmapped strains: {list(unmapped)}")
    df = df[df["strain_mic"].notna()].copy()
    df["snp_id"] = (
        df["Gene"].str.strip() + "__" +
        df["Position"].astype(str) + "__" +
        df["Reference"].str.strip() + "__" +
        df["SNP"].str.strip()
    )
    return df


def load_mics(mic_file: Path) -> pd.DataFrame:
    print(f"  Loading MICs: {mic_file}")
    df = pd.read_csv(mic_file, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip()
    return df


def build_snp_matrix(snps: pd.DataFrame, all_mic_strains: list) -> pd.DataFrame:
    print(f"  Building SNP matrix...")
    mat = (snps
           .groupby(["snp_id", "strain_mic"])
           .size()
           .unstack(fill_value=0)
           .clip(upper=1))
    for s in all_mic_strains:
        if s not in mat.columns:
            mat[s] = 0
    return mat[all_mic_strains]


def run_fisher_tests(snp_mat: pd.DataFrame,
                     mics: pd.DataFrame,
                     antibiotics: list,
                     genome_strains: set) -> pd.DataFrame:
    ref_scores = mics.loc[REFERENCE]
    records    = []
    total_abs  = len(antibiotics)
    test_strains = [s for s in snp_mat.columns if s in genome_strains]
    snp_test     = snp_mat[test_strains]

    print(f"  Strains used for Fisher tests: {len(test_strains)}")
    print(f"  Resistance threshold: ref_score + {MIN_DELTA} (MIN_DELTA)")

    for ab_idx, ab in enumerate(antibiotics, 1):
        if ab not in mics.columns:
            continue

        ref_score = ref_scores[ab]

        # Resistant = score >= ref + MIN_DELTA
        resistant = set(s for s in test_strains
                        if s in mics.index
                        and mics.loc[s, ab] >= ref_score + MIN_DELTA)
        sensitive = set(test_strains) - resistant

        n_res = len(resistant)
        n_tot = len(test_strains)

        if n_res < MIN_STRAINS:
            continue

        print(f"  [{ab_idx:02d}/{total_abs}] {ab}: "
              f"{n_res} resistant / {n_tot} strains, "
              f"{len(snp_test):,} SNPs...", end="\r")

        for snp_id, row in snp_test.iterrows():
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
                "antibiotic":    ab,
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
    print("=" * 65)
    print("  SNP–MIC Association Analysis — Fisher's Exact Test + BH FDR")
    print(f"  Resistance threshold: ref_score + {MIN_DELTA}")
    print("=" * 65)

    snps = load_snps(RAW  / "allnsSNP.tsv")
    mics = load_mics(PROC / "MICs_converted.tsv")

    mapped_mic_names = set(snps["strain_mic"].dropna().unique())
    missing = mapped_mic_names - set(mics.index)
    if missing:
        print(f"\n  [WARNING] Mapped strains not in MIC table: {missing}")
    else:
        print(f"\n  All {len(mapped_mic_names)} mapped strains found in MIC table. OK")

    antibiotics = [ab for ab in mics.columns
                   if REFERENCE in mics.index
                   and not pd.isna(mics.loc[REFERENCE, ab])]
    print(f"  Antibiotics to test: {len(antibiotics)}")

    # Show resistant strain counts with new threshold
    print(f"\n  Resistant strains per antibiotic (threshold = ref + {MIN_DELTA}):")
    genome_strains_list = sorted(snps["strain_mic"].unique())
    ref_scores = mics.loc[REFERENCE]
    for ab in antibiotics:
        n_res = sum(1 for s in genome_strains_list
                    if s in mics.index
                    and mics.loc[s, ab] >= ref_scores[ab] + MIN_DELTA)
        if n_res >= MIN_STRAINS:
            print(f"    {ab:8s}: {n_res:2d} resistant")

    all_mic_strains = mics.index.tolist()
    snp_mat = build_snp_matrix(snps, all_mic_strains)
    print(f"\n  SNP matrix: {snp_mat.shape[0]:,} SNPs × {snp_mat.shape[1]} strains")
    print(f"\nRunning Fisher's exact tests...\n")

    genome_strains = set(snps["strain_mic"].unique())
    results = run_fisher_tests(snp_mat, mics, antibiotics, genome_strains)
    print(f"\nTotal tests: {len(results):,}")

    if len(results) == 0:
        print("No tests performed — check strain mapping.")
        exit(1)

    # Add amino acid info
    aa_info = (snps.drop_duplicates("snp_id")
               .set_index("snp_id")[["Aa_ref", "Aa_mod"]])
    results = results.join(aa_info, on="snp_id", how="left")

    # BH FDR per antibiotic
    print("Applying Benjamini-Hochberg FDR correction...")
    fdr_chunks = []
    for ab, grp in results.groupby("antibiotic"):
        grp = grp.copy()
        grp["fdr"] = bh_correction(grp["p_value"])
        fdr_chunks.append(grp)
    results = pd.concat(fdr_chunks).sort_values(["antibiotic", "fdr"])
    results["fdr"] = results["fdr"].round(6)

    out_all = PROC / "SNP_associations.tsv"
    results.to_csv(out_all, sep="\t", index=False)
    print(f"All associations: {out_all}")

    sig = results[(results["fdr"] < FDR_THRESHOLD) & (results["odds_ratio"] > 1)]
    out_sig = PROC / "SNP_significant.tsv"
    sig.to_csv(out_sig, sep="\t", index=False)
    print(f"Significant (FDR<{FDR_THRESHOLD}, OR>1): {len(sig)}")

    rows = []
    for ab, grp in results.groupby("antibiotic"):
        s = grp[(grp["fdr"] < FDR_THRESHOLD) & (grp["odds_ratio"] > 1)]
        rows.append({
            "antibiotic":            ab,
            "n_tests":               len(grp),
            "n_significant":         len(s),
            "n_resistant_strains":   grp.iloc[0]["n_resistant"] if len(grp) else 0,
            "top_gene":              s.iloc[0]["gene"]       if len(s) else "—",
            "top_fdr":               s.iloc[0]["fdr"]        if len(s) else None,
            "top_OR":                s.iloc[0]["odds_ratio"] if len(s) else None,
        })
    summary = pd.DataFrame(rows).sort_values("n_significant", ascending=False)
    out_sum = RES / "SNP_association_summary.tsv"
    summary.to_csv(out_sum, sep="\t", index=False)

    print(f"\n── Per-antibiotic summary ──")
    print(summary[summary["n_significant"] > 0].to_string(index=False))
    print(f"\nSummary: {out_sum}")
    print("\nDone.")

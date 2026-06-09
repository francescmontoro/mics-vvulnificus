#!/usr/bin/env python3
"""
08_scoary_prep.py
-----------------
Prepares input files for Scoary phylogenetic association analysis.

Scoary corrects for population structure by using the phylogenetic tree
to distinguish true genotype-phenotype associations from spurious ones
caused by lineage effects.

Inputs:
    ../data/raw/allnsSNP.tsv
    ../data/raw/phylo_tree.nwk
    ../data/processed/MICs_converted.tsv

Outputs:
    ../results/scoary_genes.csv    (SNP presence/absence matrix)
    ../results/scoary_traits.csv   (resistance phenotype per strain)

Then run Scoary manually:
    scoary -g ../results/scoary_genes.csv \
           -t ../results/scoary_traits.csv \
           -n ../data/raw/phylo_tree.nwk \
           -o ../results/scoary_output/ \
           --correction I --p_value_cutoff 0.05 -e 0.05

Author: Francesc Montoro (2025)
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
RAW  = BASE / "data" / "raw"
PROC = BASE / "data" / "processed"
RES  = BASE / "results"
RES.mkdir(exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
REFERENCE  = "CECT_5769"
MIN_DELTA  = 0.5
MIN_RES    = 5    # minimum resistant strains to test
MIN_SEN    = 5    # minimum sensitive strains to test

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

# ── Load SNPs ─────────────────────────────────────────────────────────────────
print("Loading SNPs...")
snps = pd.read_csv(RAW / "allnsSNP.tsv", sep="\t")
snps.columns = snps.columns.str.strip()
snps["Strain"] = snps["Strain"].astype(str).str.strip()
snps["strain_mic"] = snps["Strain"].map(SNP_TO_MIC)
snps = snps[snps["strain_mic"].notna()].copy()

snps["snp_id"] = (
    snps["Gene"].str.strip() + "__" +
    snps["Position"].astype(str) + "__" +
    snps["Reference"].str.strip() + "__" +
    snps["SNP"].str.strip()
)

# ── Build SNP presence/absence matrix ─────────────────────────────────────────
print("Building SNP matrix...")
mat = (snps
       .groupby(["snp_id", "strain_mic"])
       .size()
       .unstack(fill_value=0)
       .clip(upper=1))

genome_strains = sorted(mat.columns.tolist())
print(f"  SNPs: {len(mat):,}, Strains: {len(genome_strains)}")

# ── Load MICs ─────────────────────────────────────────────────────────────────
print("Loading MICs...")
mics = pd.read_csv(PROC / "MICs_converted.tsv", sep="\t", index_col=0)
mics.index = mics.index.astype(str).str.strip()
ref = mics.loc[REFERENCE]

# ── Build traits table ────────────────────────────────────────────────────────
# 1 = resistant (MIC >= ref + MIN_DELTA), 0 = sensitive
print("Building traits table...")
traits = {}
skipped = []

for ab in mics.columns:
    ab_data = {}
    n_res = 0
    n_sen = 0
    for s in genome_strains:
        if s in mics.index and not pd.isna(mics.loc[s, ab]):
            resistant = int(mics.loc[s, ab] >= ref[ab] + MIN_DELTA)
            ab_data[s] = resistant
            if resistant:
                n_res += 1
            else:
                n_sen += 1
        else:
            ab_data[s] = 0

    if n_res >= MIN_RES and n_sen >= MIN_SEN:
        traits[ab] = ab_data
    else:
        skipped.append(f"{ab} (res={n_res}, sen={n_sen})")

traits_df = pd.DataFrame(traits, index=genome_strains)

print(f"\nAntibiotics included in Scoary ({len(traits)}):")
for ab in traits_df.columns:
    n_res = traits_df[ab].sum()
    print(f"  {ab:8s}: {n_res:2d} resistant / {len(genome_strains)} strains")

print(f"\nAntibiotics skipped (insufficient resistant or sensitive strains):")
for s in skipped:
    print(f"  {s}")

# ── Check tree tip names match ────────────────────────────────────────────────
tree_file = RAW / "phylo_tree.nwk"
if tree_file.exists():
    tree_text = tree_file.read_text()
    print(f"\nChecking tree tip names...")
    missing_in_tree = [s for s in genome_strains if s not in tree_text]
    if missing_in_tree:
        print(f"  [WARNING] These strains are NOT in the tree: {missing_in_tree}")
        print(f"  Scoary requires all strains to be in the tree.")
        print(f"  The tree uses original SNP names — needs renaming.")
    else:
        print(f"  All {len(genome_strains)} strains found in tree. OK")
else:
    print(f"\n[WARNING] Tree file not found: {tree_file}")

# ── Save ──────────────────────────────────────────────────────────────────────
out_genes  = RES / "scoary_genes.csv"
out_traits = RES / "scoary_traits.csv"

# Scoary format: genes file has Gene as first column, then strains
mat_out = mat.T.copy()
mat_out.index.name = "Gene"
mat_out.to_csv(out_genes)

traits_df.index.name = "Strain"
traits_df.to_csv(out_traits)

print(f"\nSaved:")
print(f"  {out_genes}  ({mat_out.shape[0]:,} SNPs × {mat_out.shape[1]} strains)")
print(f"  {out_traits}  ({traits_df.shape[0]} strains × {traits_df.shape[1]} antibiotics)")

print(f"""
Next step — run Scoary:

    scoary -g {out_genes} \\
           -t {out_traits} \\
           -n {tree_file} \\
           -o {RES}/scoary_output/ \\
           --correction I --p_value_cutoff 0.05 -e 0.05

Note: The tree tip names must match the strain names in the CSV files.
      If the tree uses original SNP names (e.g. 'CECT5769' instead of
      'CECT_5769'), rename them first with 08b_rename_tree.py.
""")

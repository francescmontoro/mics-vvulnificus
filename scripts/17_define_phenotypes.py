#!/usr/bin/env python3
"""
17_define_phenotypes.py
=======================
Builds the phenotype files for the association analysis from the RAW MIC
values, using distribution-based cut-offs instead of a fixed offset from a
reference strain.

WHY THIS REPLACES THE PREVIOUS DEFINITION
-----------------------------------------
The earlier definition called a strain resistant when its normalised score was
>= reference + 0.5. Examining the raw MIC distributions across the full
116-isolate collection showed two problems:

  GENTAMICIN (GEN.1). The distribution is unimodal and continuous
  (<0.5: 51 strains, 0.5: 1, 1: 50, 2: 10, 4: 4). Twelve of the fourteen
  strains previously called resistant have MIC = 1 ug/mL, which is the modal
  value of the collection and a single dilution above the reference. A single
  dilution is within the reproducibility of broth microdilution. There is no
  wild-type/non-wild-type separation to cut. Consistently, none of these
  strains carries any known aminoglycoside determinant.
  -> GEN is retained ONLY as a worked example of a poorly defined phenotype.

  SULPHATHIAZOLE (STZ). The distribution IS bimodal with a clear valley
  (<32: 89, 64: 8, 128: 3, 256: 1, >256: 15). The three strains with a
  confirmed sul gene all sit at >256. The eight strains at 64 are unrelated by
  origin and by lineage, and one of them (BT3_12) is fully sensitive while
  sitting in the same clade as strains at 64-256; CG100 (64) is the sister
  taxon of YJ016 (>256). The variation at 64 is within clades, not between
  them, which is what a wild-type tail looks like.
  -> primary cut-off >=128, sensitivity analysis at >=64.

No formal ECOFF is claimed. With 89/116 isolates censored at <32 the wild-type
mode lies below the tested range and cannot be fitted; and the multi-centre
ECOFF study for V. harveyi explicitly discarded sulphonamide data as unsuitable
for ECOFF estimation. What is reported here is an operational cut-off derived
from the observed distribution, with the alternative shown alongside.

Outputs (pyseer format, tab-separated, header = samples <ab>):
  pheno_STZ_bin128.tsv     primary
  pheno_STZ_bin64.tsv      sensitivity
  pheno_STZ_ordinal.tsv    censored-ordinal, treated as continuous
  pheno_GEN_bin.tsv        retained as the poorly-defined example
  pheno_GEN_ordinal.tsv
  phenotype_summary.tsv    per-strain table for the manuscript

Usage:
  python3 17_define_phenotypes.py \
      --mic-xlsx  data/mic/MICs_all_Vv.xlsx \
      --outdir    results/pyseer \
      --strains   data/genomes
"""

import argparse
from pathlib import Path

import pandas as pd

# MIC-table strain name -> assembly/SNP file stem
MIC_TO_SNP = {
    "BT3_11028": "11028", "BT3_12": "12", "BT3_162": "162", "BT3_32": "32",
    "94-3-85": "94385", "960426_1/4C": "960426-1_4C", "C1P8190": "CIP8190",
    "V_4": "V4", "CECT_4604": "CECT4604", "CECT_4606": "CECT4606",
    "CECT_4608": "CECT4608", "CECT_4865": "CECT4865", "CECT_4866": "CECT4866",
    "CECT_4999": "CECT4999", "CECT_5198": "CECT5198", "CECT_5763": "CECT5763",
    "CECT_5769": "CECT5769", "CECT_7030": "CECT7030", "CECT_898": "CECT898",
    "CG_100": "CG100", "94-9-112": "94-9-112", "94-9-119": "94-9-119",
    "95-8-161": "95-8-161", "95-8-6": "95-8-6", "95-8-7": "95-8-7",
    "A14": "A14", "PD-2-51": "PD-2-51", "Ra_3": "Ra_3", "RIU1": "RIU1",
    "YJ016": "YJ016",
}

# ordinal encoding of the tested dilution series
STZ_LEVELS = {"<32": 0, "64": 1, "128": 2, "256": 3, ">256": 4}
GEN_LEVELS = {"<0.5": 0, "0.5": 1, "1": 2, "2": 3, "4": 4}

# pyseer needs non-numeric sample names
def safe(name):
    s = str(name).strip().replace("/", "_").replace("-", "_").replace(" ", "_")
    return "S_" + s


def build(series, levels, label):
    """Map raw MIC strings to ordinal levels; report anything unmapped."""
    v = series.astype(str).str.strip()
    out = v.map(levels)
    bad = sorted(set(v[out.isna()]))
    if bad:
        print(f"    WARNING [{label}]: unmapped MIC values {bad}")
    return out


def write_pheno(path, values, col):
    df = pd.DataFrame({"samples": list(values.keys()), col: list(values.values())})
    df.to_csv(path, sep="\t", index=False)
    print(f"    wrote {path}  (n={len(df)}, resistant={int(sum(df[col]))}"
          if set(df[col]) <= {0, 1} else f"    wrote {path}  (n={len(df)})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mic-xlsx", required=True)
    ap.add_argument("--sheet", default="MICs")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--strains", default="",
                    help="directory of .fna files; restricts output to the "
                         "sequenced strains")
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    d = pd.read_excel(a.mic_xlsx, sheet_name=a.sheet)
    # duplicated headers read by pandas as GEN.1/SXT.1, or plate-suffixed names
    d = d.rename(columns={"GEN_GN2F": "GEN", "GEN_AVIAN1F": "GEN.1",
                          "SXT_GN2F": "SXT", "SXT_AVIAN1F": "SXT.1"})
    d["STRAIN"] = d["STRAIN"].astype(str).str.strip()

    print("=" * 74)
    print("  Phenotype definition from raw MIC distributions")
    print("=" * 74)

    # ── full-collection distributions (this is what justifies the cut-offs) ──
    print("\n  Distribution across the whole collection (n = %d)" % len(d))
    for col, levels, name in (("STZ", STZ_LEVELS, "sulphathiazole"),
                              ("GEN.1", GEN_LEVELS, "gentamicin")):
        if col not in d.columns:
            continue
        vc = d[col].astype(str).str.strip().value_counts()
        ordered = sorted(vc.index, key=lambda k: levels.get(k, 99))
        line = "  ".join(f"{k}:{vc[k]}" for k in ordered)
        print(f"    {name:<16} {line}")

    # ── restrict to sequenced strains ───────────────────────────────────────
    wanted = None
    if a.strains:
        wanted = {p.stem for p in Path(a.strains).glob("*.fna")}
        print(f"\n  sequenced strains found: {len(wanted)}")

    rows = []
    for _, r in d.iterrows():
        snp = MIC_TO_SNP.get(r["STRAIN"])
        if not snp:
            continue
        if wanted is not None and snp not in wanted:
            continue
        stz_raw = str(r.get("STZ", "")).strip()
        gen_raw = str(r.get("GEN.1", "")).strip()
        rows.append({
            "strain": snp, "safe": safe(snp), "mic_name": r["STRAIN"],
            "STZ_raw": stz_raw, "STZ_level": STZ_LEVELS.get(stz_raw),
            "GEN_raw": gen_raw, "GEN_level": GEN_LEVELS.get(gen_raw),
        })
    t = pd.DataFrame(rows)
    print(f"  strains with both genome and MIC data: {len(t)}")

    # ── phenotypes ──────────────────────────────────────────────────────────
    print("\n  Building phenotype files")

    stz128 = {r["safe"]: int(r["STZ_level"] >= 2) for _, r in t.iterrows()
              if pd.notna(r["STZ_level"])}
    stz64 = {r["safe"]: int(r["STZ_level"] >= 1) for _, r in t.iterrows()
             if pd.notna(r["STZ_level"])}
    stzord = {r["safe"]: float(r["STZ_level"]) for _, r in t.iterrows()
              if pd.notna(r["STZ_level"])}
    genbin = {r["safe"]: int(r["GEN_level"] >= 2) for _, r in t.iterrows()
              if pd.notna(r["GEN_level"])}
    genord = {r["safe"]: float(r["GEN_level"]) for _, r in t.iterrows()
              if pd.notna(r["GEN_level"])}

    write_pheno(outdir / "pheno_STZ_bin128.tsv", stz128, "STZ")
    write_pheno(outdir / "pheno_STZ_bin64.tsv", stz64, "STZ")
    write_pheno(outdir / "pheno_STZ_ordinal.tsv", stzord, "STZ")
    write_pheno(outdir / "pheno_GEN_bin.tsv", genbin, "GEN")
    write_pheno(outdir / "pheno_GEN_ordinal.tsv", genord, "GEN")

    t.to_csv(outdir / "phenotype_summary.tsv", sep="\t", index=False)

    # ── who is resistant, for the record ────────────────────────────────────
    print("\n  STZ resistant at the primary cut-off (>=128):")
    for _, r in t[t["STZ_level"] >= 2].sort_values("STZ_level",
                                                   ascending=False).iterrows():
        print(f"    {r['strain']:<12} MIC {r['STZ_raw']}")
    extra = t[(t["STZ_level"] == 1)]
    if len(extra):
        print("\n  additionally called resistant at the >=64 sensitivity cut-off:")
        for _, r in extra.iterrows():
            print(f"    {r['strain']:<12} MIC {r['STZ_raw']}")

    print(f"""
  NOTE ON GENTAMICIN
    {int(sum(genbin.values()))} strains are called resistant at MIC >= 1 ug/mL,
    but that is the modal value of the collection and one dilution above the
    reference. This phenotype is carried through the analysis only to
    demonstrate what a poorly defined phenotype does to an association study;
    it should not be presented as resistance.
""")


if __name__ == "__main__":
    main()

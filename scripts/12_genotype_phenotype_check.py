#!/usr/bin/env python3
"""
12_genotype_phenotype_check.py   (v2 - corrected mappings)
==========================================================
Cross-validates acquired resistance genes (ResFinder / AMRFinderPlus) against
the measured MIC values.

CORRECTIONS vs v1 (they change the concordance a lot):

  1. cat genes NO LONGER mapped to florfenicol (FFN).
     Classic chloramphenicol acetyltransferases (catA, catB, catA1, catA2,
     catB7) do NOT confer florfenicol resistance: florfenicol carries a
     fluorine at C-3 that cannot be acetylated. Chloramphenicol itself is not
     on the Sensititre panels used here, so cat genes have NO testable
     antibiotic and are excluded from the concordance calculation rather than
     counted as discordant. (floR would be the florfenicol determinant.)

  2. Chromosomal / intrinsic elements excluded from "acquired":
        varG      - chromosomal carbapenemase, found in 30/30 strains
        tet(35)   - chromosomal efflux, 29/30 strains
        catB      - 17/30 strains, chromosomal background
     A determinant present in essentially every strain cannot explain MIC
     differences between strains.

  3. De-duplication: a gene detected by BOTH ResFinder and AMRFinder is one
     row, not two (v1 printed every antibiotic twice).

  4. Contamination flagging: strains whose assembly is much larger than the
     expected V. vulnificus genome (~5.26 Mb) are flagged. Resistance genes in
     a contaminated assembly may belong to the contaminant rather than to the
     strain that was phenotyped. Needs results/qc/assembly_qc.tsv (script 10).

Usage:
  python3 12_genotype_phenotype_check.py \
      --resfinder-dir ../PAPER_V2/results/resfinder_raw \
      --amrfinder-summary results/amrfinder/amrfinder_summary.tsv \
      --mics data/mic/MICs_converted.tsv \
      --qc results/qc/assembly_qc.tsv \
      --outdir results/genotype_phenotype
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

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
VV_GENOME_MB = 5.26
CONTAM_RATIO = 1.15

# Prefix -> MIC columns present on the Sensititre panels.
# []   = drug not on the panel  -> not testable, excluded from concordance
# None = gene not in this table -> flagged as unmapped
GENE_TO_MIC_COLS = {
    "sul":     ["STZ", "SDM"],
    "dfr":     ["SXT", "SXT2"],
    "tet":     ["TET", "OXY"],
    "erm":     ["ERY", "TYLT", "CLI"],
    "vat":     ["ERY"],
    "msr":     ["ERY"],
    "aada":    ["STR", "SPE"],
    "ant":     ["STR", "SPE"],
    "aph":     ["NEO", "AMI", "TOB"],
    "aac":     ["GEN", "GEN2", "TOB", "AMI"],
    "blatem":  ["AMP", "AMOX", "PIP", "A/S2"],
    "blaoxa":  ["AMP", "AMOX", "PIP"],
    "blapdc":  ["TAZ", "FEP", "AXO"],
    "blaz":    ["AMP", "AMOX"],
    "bla":     ["AMP", "AMOX", "PIP"],
    "gyra":    ["CIP", "ENRO", "GAT"],
    "parc":    ["CIP", "ENRO", "GAT"],
    "flor":    ["FFN"],
    "cat":     [],
    "fosa":    [],
    "sat":     [],
    "mex":     [],
    "mer":     [],
    "qac":     [],
    "asr":     [],
}

INTRINSIC = {"varg", "tet(35)", "catb"}

CONTAMINANT_MARKERS = {
    "blapdc": "Pseudomonas aeruginosa (intrinsic AmpC)",
    "blaoxa-50": "Pseudomonas aeruginosa (intrinsic oxacillinase)",
    "mexa": "Pseudomonas aeruginosa (MexAB efflux)",
    "aph(3')-iib": "Pseudomonas aeruginosa (intrinsic)",
    "catb7": "Pseudomonas aeruginosa (intrinsic)",
    "tet(k)": "Staphylococcus / Gram-positive",
    "erm(b)": "Gram-positive (mobile, can be genuine)",
    "tet(m)": "Gram-positive (mobile, can be genuine)",
}


def is_intrinsic(gene):
    g = gene.lower()
    return any(g == i or g.startswith(i) for i in INTRINSIC)


def map_gene_to_cols(gene):
    g = gene.lower()
    for prefix in sorted(GENE_TO_MIC_COLS, key=len, reverse=True):
        if g.startswith(prefix):
            return GENE_TO_MIC_COLS[prefix]
    return None


def load_resfinder(resdir):
    rows = []
    if not resdir or not os.path.isdir(resdir):
        return rows
    for d in sorted(os.listdir(resdir)):
        tab = os.path.join(resdir, d, "ResFinder_results_tab.txt")
        if not os.path.exists(tab):
            continue
        try:
            df = pd.read_csv(tab, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue
        gcol = df.columns[0]
        for _, r in df.iterrows():
            gene = str(r[gcol]).strip()
            if not gene or gene.lower().startswith("resistance gene"):
                continue
            rows.append({"strain": d, "gene": gene, "source": "ResFinder"})
    return rows


def load_amrfinder(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    df = pd.read_csv(path, sep="\t")
    for _, r in df.iterrows():
        for g in str(r.get("genes", "") or "").split(";"):
            g = g.strip()
            if g:
                rows.append({"strain": str(r["strain"]), "gene": g,
                             "source": "AMRFinder"})
    return rows


def load_contamination(qc_path):
    flags = {}
    if not qc_path or not os.path.exists(qc_path):
        return flags
    qc = pd.read_csv(qc_path, sep="\t")
    for _, r in qc.iterrows():
        mb = float(r["total_bp"]) / 1e6
        ratio = mb / VV_GENOME_MB
        flags[str(r["strain"])] = (mb, ratio, ratio > CONTAM_RATIO)
    return flags


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resfinder-dir", default="")
    ap.add_argument("--amrfinder-summary", default="")
    ap.add_argument("--mics", required=True)
    ap.add_argument("--qc", default="")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    mics = pd.read_csv(a.mics, sep="\t", index_col=0)
    mics.index = mics.index.astype(str).str.strip()
    ref = mics.loc[REFERENCE]

    contam = load_contamination(a.qc)
    raw = load_resfinder(a.resfinder_dir) + load_amrfinder(a.amrfinder_summary)
    if not raw:
        sys.exit("ERROR: no resistance-gene hits loaded. Check input paths.")

    merged = {}
    for h in raw:
        key = (h["strain"], h["gene"])
        if key not in merged:
            merged[key] = {**h, "sources": set()}
        merged[key]["sources"].add(h["source"])

    print("=" * 74)
    print("  Genotype-phenotype cross-check (v2, corrected mappings)")
    print("=" * 74)

    rows = []
    for (strain, gene), h in sorted(merged.items()):
        mic_name = SNP_TO_MIC.get(strain, strain)
        intrinsic = is_intrinsic(gene)
        cols = map_gene_to_cols(gene)
        mb, ratio, flagged = contam.get(strain, (float("nan"),
                                                 float("nan"), False))
        marker = CONTAMINANT_MARKERS.get(gene.lower(), "")

        base = {
            "strain": strain, "mic_strain": mic_name, "gene": gene,
            "sources": ",".join(sorted(h["sources"])),
            "intrinsic": intrinsic,
            "assembly_mb": round(mb, 2) if mb == mb else "",
            "assembly_ratio": round(ratio, 2) if ratio == ratio else "",
            "contamination_flag": flagged,
            "contaminant_marker": marker,
        }

        if intrinsic:
            rows.append({**base, "antibiotic": "", "mic_score": "",
                         "resistant": "", "testable": False,
                         "note": "chromosomal/intrinsic"})
            continue
        if cols is None:
            rows.append({**base, "antibiotic": "", "mic_score": "",
                         "resistant": "", "testable": False,
                         "note": "gene not in mapping table"})
            continue
        if not cols:
            rows.append({**base, "antibiotic": "", "mic_score": "",
                         "resistant": "", "testable": False,
                         "note": "drug not on Sensititre panel"})
            continue

        for ab in cols:
            if ab not in mics.columns or mic_name not in mics.index:
                continue
            val = mics.loc[mic_name, ab]
            if pd.isna(val):
                continue
            rows.append({**base, "antibiotic": ab, "mic_score": float(val),
                         "resistant": bool(float(val) >= ref[ab] + MIN_DELTA),
                         "testable": True, "note": ""})

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "genotype_phenotype_check.tsv", sep="\t", index=False)

    t = df[(df["testable"]) & (~df["intrinsic"])]
    clean = t[~t["contamination_flag"]]

    def summarise(sub, label):
        if not len(sub):
            print(f"\n  {label}: no testable pairs")
            return
        conc = sub.groupby(["strain", "gene"])["resistant"].any()
        print(f"\n  {label}")
        print(f"    gene-strain pairs tested : {len(conc)}")
        print(f"    concordant (strain is R) : {int(conc.sum())}"
              f"  ({100.0*conc.sum()/len(conc):.0f}%)")

    summarise(t, "ALL strains")
    summarise(clean, "EXCLUDING contamination-flagged assemblies")

    print("\n  -- Acquired genes with a testable antibiotic --")
    for (strain, gene), sub in t.groupby(["strain", "gene"]):
        mic_name = sub["mic_strain"].iloc[0]
        parts = [f"{r['antibiotic']}={r['mic_score']:.2f}"
                 f"[{'R' if r['resistant'] else 's'}]"
                 for _, r in sub.iterrows()]
        src = sub["sources"].iloc[0]
        warn = ""
        if sub["contamination_flag"].iloc[0]:
            warn = f"  <-- assembly {sub['assembly_ratio'].iloc[0]}x expected"
        if sub["contaminant_marker"].iloc[0]:
            warn += f"  [{sub['contaminant_marker'].iloc[0]}]"
        print(f"    {mic_name:<12} {gene:<14} {' '.join(parts):<44} ({src}){warn}")

    excluded = df[(~df["testable"]) & (~df["intrinsic"])]
    if len(excluded):
        print("\n  -- Excluded (drug not on panel / unmapped) --")
        for (strain, gene), sub in excluded.groupby(["strain", "gene"]):
            print(f"    {sub['mic_strain'].iloc[0]:<12} {gene:<14} "
                  f"{sub['note'].iloc[0]}")

    if contam:
        bad = {s: v for s, v in contam.items() if v[2]}
        if bad:
            print("\n  -- Assemblies flagged as possibly contaminated --")
            print(f"     (expected V. vulnificus genome ~{VV_GENOME_MB} Mb)")
            for s, (mb, ratio, _) in sorted(bad.items(), key=lambda x: -x[1][1]):
                print(f"    {SNP_TO_MIC.get(s,s):<12} {mb:>6.2f} Mb  "
                      f"({ratio:.2f}x expected)")
            print("\n    Resistance genes in these assemblies may belong to the"
                  "\n    contaminant rather than to the phenotyped strain.")

    print(f"\n  wrote {outdir/'genotype_phenotype_check.tsv'}")


if __name__ == "__main__":
    main()

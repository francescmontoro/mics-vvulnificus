#!/usr/bin/env python3
"""
19_map_resistance_to_pangenome.py
=================================
Locates the acquired resistance determinants inside the Panaroo pan-genome,
without relying on gene names.

WHY THIS IS NEEDED
------------------
Bakta was run with the LIGHT database (the full 32 GB database was not usable
on 13 GB of RAM). The light database does not recognise several acquired
resistance genes: sul1 is not annotated in CECT898, nor blaTEM-116 / catA1 in
V4, even though AMRFinderPlus and ResFinder both detect them at 100% identity.
Those coding sequences are still predicted by Bakta and still enter the
pan-genome, but as unnamed group_XXXXX clusters. Searching the pan-genome by
gene symbol therefore misses exactly the determinants that matter most.

THE LINK THAT SOLVES IT
-----------------------
AMRFinderPlus was run on Bakta's protein output, so its "Protein id" column IS
the Bakta locus tag. Panaroo's gene_presence_absence.csv lists, for every gene
cluster and every strain, the locus tags it contains. Chaining the two gives:

    resistance gene -> Bakta locus tag -> Panaroo cluster id

which is what the association analysis actually tests.

Outputs:
  <outdir>/resistance_to_pangenome.tsv   one row per (strain, gene, cluster)
  <outdir>/positive_controls.txt         cluster ids to watch in the pan-genome
                                         association, one per line

Usage:
  python3 19_map_resistance_to_pangenome.py \
      --amrfinder-dir results/amrfinder \
      --pangenome     results/panaroo/gene_presence_absence.csv \
      --bakta-dir     results/bakta \
      --outdir        results/pangenome_assoc
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import pandas as pd

# determinants worth tracking as positive controls: acquired, mobile, and with
# a phenotype we can check. varG/tet(35)/catB are chromosomal background.
INTRINSIC = {"varg", "tet(35)", "tet35", "catb", "perm", "tetd", "tetr"}


def is_intrinsic(gene):
    g = str(gene).lower()
    return any(g == i or g.startswith(i) for i in INTRINSIC)


def load_amrfinder(amr_dir):
    """strain -> [(locus_tag, gene, class)] from AMRFinderPlus output."""
    out = {}
    for f in sorted(glob.glob(os.path.join(amr_dir, "*.amr.tsv"))):
        strain = os.path.basename(f)[:-len(".amr.tsv")]
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue
        cols = {c.lower().strip(): c for c in df.columns}
        pid = cols.get("protein id") or cols.get("protein identifier")
        sym = cols.get("element symbol") or cols.get("gene symbol")
        typ = cols.get("type")
        cls = cols.get("class")
        if not pid or not sym:
            continue
        rows = []
        for _, r in df.iterrows():
            if typ and str(r[typ]).strip().upper() != "AMR":
                continue
            tag = str(r[pid]).strip()
            gene = str(r[sym]).strip()
            if not tag or tag.upper() == "NA":
                continue
            rows.append((tag, gene, str(r[cls]).strip() if cls else ""))
        if rows:
            out[strain] = rows
    return out


def load_panaroo(csv_path):
    """locus_tag -> cluster id, from gene_presence_absence.csv."""
    tag_to_cluster = {}
    df = pd.read_csv(csv_path, low_memory=False)
    gene_col = df.columns[0]
    # first 3 columns are metadata in Panaroo/Roary format; rest are strains
    strain_cols = [c for c in df.columns[3:]]
    for _, r in df.iterrows():
        cluster = str(r[gene_col])
        for c in strain_cols:
            v = r[c]
            if pd.isna(v):
                continue
            for tag in str(v).replace("\t", ";").split(";"):
                tag = tag.strip()
                if tag:
                    tag_to_cluster[tag] = cluster
    return tag_to_cluster, strain_cols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amrfinder-dir", required=True)
    ap.add_argument("--pangenome", required=True)
    ap.add_argument("--bakta-dir", default="")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 76)
    print("  Mapping resistance determinants onto Panaroo clusters")
    print("=" * 76)

    amr = load_amrfinder(a.amrfinder_dir)
    if not amr:
        sys.exit(f"ERROR: no AMRFinder results parsed from {a.amrfinder_dir}")
    print(f"  strains with AMR hits : {len(amr)}")

    tag2cl, strain_cols = load_panaroo(a.pangenome)
    print(f"  locus tags in pan-genome : {len(tag2cl)}")
    print(f"  strains in pan-genome    : {len(strain_cols)}\n")

    rows, unmapped = [], []
    for strain, hits in sorted(amr.items()):
        for tag, gene, cls in hits:
            cluster = tag2cl.get(tag)
            rows.append({
                "strain": strain, "gene": gene, "class": cls,
                "bakta_locus": tag,
                "panaroo_cluster": cluster if cluster else "",
                "intrinsic": is_intrinsic(gene),
                "mapped": bool(cluster),
            })
            if not cluster:
                unmapped.append((strain, gene, tag))

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "resistance_to_pangenome.tsv", sep="\t", index=False)

    acq = df[(~df["intrinsic"]) & (df["mapped"])]
    print("  -- acquired determinants and their pan-genome cluster --")
    for _, r in acq.sort_values(["gene", "strain"]).iterrows():
        named = "" if r["panaroo_cluster"].startswith("group_") else "  (named)"
        print(f"    {r['strain']:<12} {r['gene']:<14} -> "
              f"{r['panaroo_cluster']:<18}{named}")

    if unmapped:
        print(f"\n  -- {len(unmapped)} hits whose locus tag is NOT in the "
              f"pan-genome --")
        for s, g, t in unmapped[:20]:
            print(f"    {s:<12} {g:<14} {t}")
        print("    (these CDS were dropped by Panaroo, most likely as "
              "fragmented\n     or invalid genes; they cannot be tested)")

    # ── clusters to use as positive controls ────────────────────────────────
    controls = sorted(set(acq["panaroo_cluster"]))
    with open(outdir / "positive_controls.txt", "w") as fh:
        for c in controls:
            fh.write(c + "\n")

    # summarise how many strains carry each control cluster
    pan = pd.read_csv(a.pangenome, low_memory=False)
    gene_col = pan.columns[0]
    pan = pan.set_index(gene_col)
    print(f"\n  -- {len(controls)} clusters to watch in the association --")
    for c in controls:
        genes = sorted(set(acq[acq["panaroo_cluster"] == c]["gene"]))
        if c in pan.index:
            row = pan.loc[c, pan.columns[2:]] if len(pan.columns) > 3 else None
            n = int(row.notna().sum()) if row is not None else -1
        else:
            n = -1
        print(f"    {c:<18} n_strains={n:<3} genes: {','.join(genes)}")

    print(f"""
  wrote {outdir/'resistance_to_pangenome.tsv'}
  wrote {outdir/'positive_controls.txt'}

  NEXT
    Pass positive_controls.txt to the pan-genome association step and check
    where those clusters rank. A cluster carried by strains that are all
    phenotypically resistant is the strongest available internal validation
    of the pipeline; if such a cluster still fails to reach significance, the
    power ceiling is demonstrated on a determinant whose effect is not in
    doubt.
""")


if __name__ == "__main__":
    main()

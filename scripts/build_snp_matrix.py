#!/usr/bin/env python3
"""
build_snp_matrix.py   (reviewed edition)

Takes Snippy's snps.tab (one per strain, vs the SAME reference) and builds
the variant presence/absence matrix ready for association testing with MICs.

INCLUSION CRITERION (the study's):
  - CODING region    : NON-SYNONYMOUS variants only
                       (missense, stop_gained, stop_lost, start_lost, ...)
                       -> synonymous_variant and stop_retained are DROPPED.
  - NONCODING region : ALL variants
                       (intergenic_region, upstream/downstream_gene_variant,
                        UTR, ...) -> where regulatory mechanism lives.

CHANGES vs Fran's original (all documented; defaults chosen conservatively):
  * Removed from the default NONSYN set (they inflate counts and are low-
    confidence generic effects; can be re-added with --loose-coding):
        coding_sequence_variant, protein_protein_contact,
        structural_interaction_variant, start_retained_variant
  * intragenic_variant moved OUT of the default noncoding set (SnpEff emits it
    for "inside a gene but effect undetermined"; ambiguous). Re-add with
    --keep-intragenic.
  * --min-strains and --prefix names unchanged (these are what the pipeline
    now calls, after fixing the shell script).

Snippy snps.tab columns:
    CHROM POS TYPE REF ALT EVIDENCE FTYPE STRAND NT_POS AA_POS
    EFFECT LOCUS_TAG GENE PRODUCT

Outputs:
    <prefix>_matrix.tsv    rows = variants, cols = strains (0/1)
    <prefix>_variants.tsv  per-variant annotation (region, effect, gene, ...)
    <prefix>_summary.tsv   per-category counts (for the Results funnel)

Usage:
    python build_snp_matrix.py --snps-dir results/snps \
        --rename-map results/summary/gene_rename_map.tsv \
        --min-strains 5 --prefix results/summary/variants
"""
import argparse
import csv
import glob
import os
import sys
from collections import Counter, defaultdict

# --- SnpEff effect classification (the ones Snippy emits) -------------------

# Coding AND changes the protein -> KEPT (conservative core set)
NONSYN_STRICT = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "initiator_codon_variant",
}

# Extra generic/structural effects: only added with --loose-coding
NONSYN_LOOSE = {
    "coding_sequence_variant",
    "protein_protein_contact",
    "structural_interaction_variant",
    "start_retained_variant",
}

# Coding AND does NOT change the protein -> DROPPED
SYN = {
    "synonymous_variant",
    "stop_retained_variant",
}

# Noncoding -> ALL KEPT (conservative core set)
NONCODING_STRICT = {
    "intergenic_region",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "5_prime_UTR_premature_start_codon_gain_variant",
    "non_coding_transcript_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
}

# Ambiguous: only treated as noncoding with --keep-intragenic
NONCODING_LOOSE = {
    "intragenic_variant",
}


def build_classifier(loose_coding, keep_intragenic):
    nonsyn = set(NONSYN_STRICT)
    if loose_coding:
        nonsyn |= NONSYN_LOOSE
    noncoding = set(NONCODING_STRICT)
    if keep_intragenic:
        noncoding |= NONCODING_LOOSE

    def classify(effect, ftype):
        """Return (region, keep). region in {coding, noncoding, other}."""
        if not effect:
            # Snippy leaves EFFECT empty when there is no annotation there.
            # No FTYPE -> outside any feature -> noncoding.
            return ("noncoding", True) if not ftype else ("other", False)

        parts = {e.strip() for e in effect.split("&")}
        parts |= {effect.split()[0].strip()}

        if parts & nonsyn:
            return ("coding", True)
        if parts & SYN:
            return ("coding", False)
        if parts & noncoding:
            return ("noncoding", True)
        # frameshift / inframe indels change the protein
        if any("frameshift" in p or "inframe" in p or "disruptive" in p
               for p in parts):
            return ("coding", True)
        return ("other", False)

    return classify


def load_rename(path):
    """bakta_locus (or prokka_locus) -> (yj016_locus, gene)"""
    m = {}
    if not path or not os.path.exists(path):
        return m
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = row.get("bakta_locus") or row.get("prokka_locus") or ""
            m[key] = (row.get("yj016_locus", ""), row.get("gene", ""))
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snps-dir", required=True,
                    help="directory with <strain>.snps.tab from Snippy")
    ap.add_argument("--rename-map", default="",
                    help="gene_rename_map.tsv (optional; adds YJ016 locus)")
    ap.add_argument("--min-strains", type=int, default=5,
                    help="variant present in >= N strains (default 5)")
    ap.add_argument("--max-strains", type=int, default=0,
                    help="and in <= N strains; 0 = n_strains - min_strains")
    ap.add_argument("--keep-indels", action="store_true",
                    help="include indels in addition to substitutions")
    ap.add_argument("--loose-coding", action="store_true",
                    help="also count generic/structural coding effects "
                         "(coding_sequence_variant, protein_protein_contact, "
                         "structural_interaction_variant, start_retained)")
    ap.add_argument("--keep-intragenic", action="store_true",
                    help="treat intragenic_variant as noncoding (ambiguous)")
    ap.add_argument("--prefix", default="variants")
    a = ap.parse_args()

    classify = build_classifier(a.loose_coding, a.keep_intragenic)

    files = sorted(glob.glob(os.path.join(a.snps_dir, "*.snps.tab")))
    if not files:
        sys.exit(f"ERROR: no *.snps.tab in {a.snps_dir}")
    strains = [os.path.basename(f).replace(".snps.tab", "") for f in files]
    print(f"strains: {len(strains)}", file=sys.stderr)
    print(f"coding mode : {'loose' if a.loose_coding else 'strict'}",
          file=sys.stderr)
    print(f"intragenic  : {'kept as noncoding' if a.keep_intragenic else 'dropped'}",
          file=sys.stderr)

    rename = load_rename(a.rename_map)
    if rename:
        print(f"rename map: {len(rename)} entries", file=sys.stderr)

    presence = defaultdict(set)
    annot = {}
    stats = Counter()

    for f, s in zip(files, strains):
        with open(f) as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                vtype = (row.get("TYPE") or "").strip()
                stats[f"raw_{vtype}"] += 1
                if not a.keep_indels and vtype != "snp":
                    stats["drop_indel"] += 1
                    continue

                effect = (row.get("EFFECT") or "").strip()
                ftype = (row.get("FTYPE") or "").strip()
                region, keep = classify(effect, ftype)
                stats[f"region_{region}"] += 1
                if not keep:
                    stats[f"drop_{region}"] += 1
                    continue

                key = (row["CHROM"], int(row["POS"]),
                       row.get("REF", ""), row.get("ALT", ""))
                presence[key].add(s)
                if key not in annot:
                    plocus = (row.get("LOCUS_TAG") or "").strip()
                    yj, gname = rename.get(plocus, ("", ""))
                    annot[key] = dict(
                        region=region,
                        effect=effect.split()[0] if effect else "intergenic",
                        bakta_locus=plocus,
                        yj016_locus=yj,
                        gene=gname or (row.get("GENE") or "").strip(),
                        product=(row.get("PRODUCT") or "").strip(),
                    )

    n = len(strains)
    lo = a.min_strains
    hi = a.max_strains if a.max_strains > 0 else n - a.min_strains
    if hi < lo:
        hi = n

    kept = [k for k, v in presence.items() if lo <= len(v) <= hi]
    kept.sort()
    print(f"\nvariants after effect filter : {len(presence):,}", file=sys.stderr)
    print(f"variants with {lo} <= n_strains <= {hi} : {len(kept):,}",
          file=sys.stderr)

    os.makedirs(os.path.dirname(os.path.abspath(a.prefix)) or ".", exist_ok=True)

    with open(f"{a.prefix}_matrix.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_id"] + strains)
        for k in kept:
            vid = f"{k[0]}:{k[1]}:{k[2]}>{k[3]}"
            w.writerow([vid] + [1 if s in presence[k] else 0 for s in strains])

    with open(f"{a.prefix}_variants.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant_id", "chrom", "pos", "ref", "alt", "region",
                    "effect", "yj016_locus", "gene", "bakta_locus",
                    "product", "n_strains"])
        for k in kept:
            an = annot[k]
            vid = f"{k[0]}:{k[1]}:{k[2]}>{k[3]}"
            w.writerow([vid, k[0], k[1], k[2], k[3], an["region"],
                        an["effect"], an["yj016_locus"], an["gene"],
                        an["bakta_locus"], an["product"], len(presence[k])])

    by_region = Counter(annot[k]["region"] for k in kept)
    with open(f"{a.prefix}_summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "count"])
        for k, v in sorted(stats.items()):
            w.writerow([k, v])
        w.writerow(["--", "--"])
        w.writerow(["kept_total", len(kept)])
        for r, c in sorted(by_region.items()):
            w.writerow([f"kept_{r}", c])

    print("\nvariants kept by region:", file=sys.stderr)
    for r, c in sorted(by_region.items()):
        print(f"  {r:<10} {c:>8,}", file=sys.stderr)
    print(f"\nwritten: {a.prefix}_matrix.tsv", file=sys.stderr)
    print(f"written: {a.prefix}_variants.tsv", file=sys.stderr)
    print(f"written: {a.prefix}_summary.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
make_rename_map.py   (Bakta edition)

Map  locus_tag_Bakta -> official YJ016 locus_tag (VV####/VVA####)
by protein-protein BLAST between the proteins Bakta predicted for the
reference and the official NCBI proteins.

Panaroo does NOT see the official annotation: every strain is annotated the
same way with Bakta (homogeneous annotation, required for clustering). The
VV/VVA nomenclature is applied AFTERWARDS, as a label over the families.

Output TSV: bakta_locus  yj016_locus  gene  product  pident  qcov  replicon
Prints how many Bakta CDS have no official equivalent: that number is the
QC control (Pyrodigal vs the original 2003 annotation).

CHANGES vs Fran's original:
  * --prokka-faa renamed to --bakta-faa (kept --prokka-faa as hidden alias
    for backwards compatibility).
  * BLAST query-id handling unchanged (Bakta FAA headers are '>locus_tag ...',
    SeqIO takes the first token as id, which is the locus_tag — same as Prokka).
"""
import argparse, csv, os, subprocess, sys, tempfile
from Bio import SeqIO


def parse_official_gff(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            attrs = {}
            for kv in p[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v.replace("%2C", ",").replace("%3B", ";")
            lt = attrs.get("locus_tag")
            if not lt:
                continue
            rec = (lt, attrs.get("gene", ""), attrs.get("product", ""), p[0])
            pid = attrs.get("protein_id") or attrs.get("Name")
            if pid:
                out[pid] = rec
            out[lt] = rec
    return out


def main():
    ap = argparse.ArgumentParser()
    # Bakta faa (accept old --prokka-faa as alias so nothing breaks)
    ap.add_argument("--bakta-faa", "--prokka-faa", dest="bakta_faa",
                    required=True, help="Bakta .faa of the reference")
    ap.add_argument("--official-faa", required=True)
    ap.add_argument("--official-gff", required=True)
    ap.add_argument("--min-identity", type=float, default=90.0)
    ap.add_argument("--min-coverage", type=float, default=90.0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    for f in (a.bakta_faa, a.official_faa, a.official_gff):
        if not os.path.exists(f):
            sys.exit(f"ERROR: does not exist {f}")

    gffmap = parse_official_gff(a.official_gff)
    qlen = {r.id: len(r.seq) for r in SeqIO.parse(a.bakta_faa, "fasta")}
    print(f"Bakta CDS (reference): {len(qlen)}", file=sys.stderr)

    tmp = tempfile.mkdtemp(prefix="renamemap_")
    db = os.path.join(tmp, "official")
    subprocess.run(["makeblastdb", "-in", a.official_faa, "-dbtype", "prot",
                    "-out", db], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hits = os.path.join(tmp, "hits.tsv")
    subprocess.run(["blastp", "-query", a.bakta_faa, "-db", db,
                    "-outfmt", "6 qseqid sseqid pident length qlen slen bitscore",
                    "-evalue", "1e-10", "-max_target_seqs", "5",
                    "-num_threads", str(a.threads), "-out", hits], check=True)

    best = {}
    with open(hits) as fh:
        for line in fh:
            q, s, pid, ln, ql, sl, bits = line.rstrip("\n").split("\t")
            pid, ln, ql, sl, bits = float(pid), int(ln), int(ql), int(sl), float(bits)
            if pid < a.min_identity:
                continue
            if 100.0 * ln / ql < a.min_coverage or 100.0 * ln / sl < a.min_coverage:
                continue
            if q not in best or bits > best[q][2]:
                best[q] = (s, pid, bits, 100.0 * ln / ql)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    n_ok = 0
    with open(a.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["bakta_locus", "yj016_locus", "gene", "product",
                    "pident", "qcov", "replicon"])
        for q in sorted(qlen):
            hit = best.get(q)
            if not hit:
                w.writerow([q, "", "", "", "", "", ""])
                continue
            s, pid, bits, qcov = hit
            rec = gffmap.get(s) or gffmap.get(s.split("|")[-1])
            if not rec:
                w.writerow([q, "", "", "", f"{pid:.1f}", f"{qcov:.1f}", ""])
                continue
            lt, gene, prod, seqid = rec
            w.writerow([q, lt, gene, prod, f"{pid:.1f}", f"{qcov:.1f}", seqid])
            n_ok += 1

    unmapped = len(qlen) - n_ok
    pct = 100.0 * unmapped / max(1, len(qlen))
    print(f"mapped to YJ016 locus  : {n_ok}", file=sys.stderr)
    print(f"no official equivalent : {unmapped} ({pct:.1f}%)", file=sys.stderr)
    if pct > 5:
        print("  WARNING: >5% unmapped -> Pyrodigal and the original annotation "
              "disagree on CDS boundaries. Report it in Methods.",
              file=sys.stderr)
    print(f"written: {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

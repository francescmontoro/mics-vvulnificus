#!/usr/bin/env bash
# =============================================================================
# run_vv_variants.sh   (v2 — Bakta edition)
#
# Variant pipeline for the MIC-SNP study of Vibrio vulnificus.
#
# CHANGES vs Fran's original (documented for transparency with co-authors):
#   * Prokka  ->  Bakta 1.12   (Prokka is unmaintained since 2019 and has an
#                               unsolvable Perl dependency conflict in current
#                               bioconda; Bakta is its modern successor, gives
#                               a versioned reproducible DB, and Panaroo
#                               accepts its GFF3 output identically).
#   * Fixed:  step 5 called build_variant_matrix.py -> build_snp_matrix.py
#   * Fixed:  --out-prefix/--min-depth  ->  --prefix/--min-strains
#             (these did not match build_snp_matrix.py's actual arguments)
#   * Added:  --min-strains is now a real CLI option (was silently dropped)
#   * Added:  AMRFinderPlus run per genome (independent acquired-resistance
#             check, complements ResFinder for the paper)
#   * MIN_COV is now actually passed to snippy (--mincov); before it was dead.
#
# Produces, from N genomes and an annotated reference (YJ016):
#   - homogeneous annotation (Bakta, all strains incl. reference)
#   - pangenome / core genome (Panaroo) -> core_gene_alignment.aln
#   - name map bakta_locus -> VV####/VVA#### (BLAST vs official YJ016)
#   - per-strain variants vs reference (Snippy), genome-wide
#   - filtered variant matrix:
#         CODING     region -> non-synonymous SNPs only
#         NONCODING  region -> all SNPs (promoters, intergenic)
#
# Usage:
#   bash run_vv_variants.sh -r data/reference/YJ016.fna -d data/genomes/ \
#        --bakta-db data/bakta_db/db -t 12 -p 3 -o results
#
# Next to the reference (.fna) there must be its official .gff and .faa
# (same basename) to recover the VV/VVA nomenclature.
#
# Requirements (conda envs): vv_annot (bakta, amrfinderplus, blast),
#                            vv_core (panaroo, iqtree), snippy (own env)
# =============================================================================
set -euo pipefail

RED='\033[0;31m'; GRN='\033[0;32m'; YEL='\033[1;33m'; CYN='\033[0;36m'
BLD='\033[1m'; RST='\033[0m'
log()  { echo -e "${GRN}[✔]${RST} $*"; }
warn() { echo -e "${YEL}[!]${RST} $*"; }
err()  { echo -e "${RED}[✘]${RST} $*" >&2; exit 1; }
step() { echo -e "\n${BLD}${CYN}══ $* ══${RST}"; }

# ── Defaults ──────────────────────────────────────────────────────────────────
REF=""; GDIR=""; OUT="results"; THREADS=8; PARALLEL=1
GENUS="Vibrio"; SPECIES="vulnificus"
CORE_THRESHOLD="0.98"
MIN_COV=10           # min coverage to accept a Snippy variant (now used!)
MIN_STRAINS=5        # variant must be present in >= N strains (now a real flag)
BAKTA_DB=""          # path to the Bakta database (required)
# conda environment names (override if yours differ)
ENV_ANNOT="vv_annot"     # bakta + amrfinderplus + blast
ENV_CORE="vv_core"       # panaroo + iqtree
ENV_SNIPPY="snippy"      # snippy (its own env)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Arguments ─────────────────────────────────────────────────────────────────
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -r) REF="$2"; shift ;;
        -d) GDIR="$2"; shift ;;
        -o) OUT="$2"; shift ;;
        -t) THREADS="$2"; shift ;;
        -p) PARALLEL="$2"; shift ;;
        -G) GENUS="$2"; shift ;;
        -S) SPECIES="$2"; shift ;;
        --bakta-db)       BAKTA_DB="$2"; shift ;;
        --min-cov)        MIN_COV="$2"; shift ;;
        --min-strains)    MIN_STRAINS="$2"; shift ;;
        --core-threshold) CORE_THRESHOLD="$2"; shift ;;
        --env-annot)      ENV_ANNOT="$2"; shift ;;
        --env-core)       ENV_CORE="$2"; shift ;;
        --env-snippy)     ENV_SNIPPY="$2"; shift ;;
        -h|--help) sed -n '2,55p' "$0"; exit 0 ;;
        *) err "Unknown argument: $1" ;;
    esac
    shift
done

[[ -z "$REF"  ]] && err "Missing -r reference.fna"
[[ -z "$GDIR" ]] && err "Missing -d genomes_directory/"
[[ -z "$BAKTA_DB" ]] && err "Missing --bakta-db path/to/bakta/db"
[[ -f "$REF"  ]] || err "Reference not found: $REF"
[[ -d "$GDIR" ]] || err "Directory not found: $GDIR"
[[ -d "$BAKTA_DB" ]] || err "Bakta DB not found: $BAKTA_DB"

# helper: run a command inside a conda env
crun() { local env="$1"; shift; mamba run -n "$env" "$@"; }

CPUS_PER_JOB=$(( THREADS / PARALLEL )); [[ "$CPUS_PER_JOB" -lt 1 ]] && CPUS_PER_JOB=1

REF=$(realpath "$REF")
REF_ID="$(basename "$REF" | sed 's/\.[^.]*$//')"
REF_BASE="${REF%.*}"

# official annotation of the reference (to recover VV####/VVA####)
REF_GFF=""
for ext in gff gff3 gtf; do
    [[ -f "${REF_BASE}.${ext}" ]] && { REF_GFF="${REF_BASE}.${ext}"; break; }
done
REF_FAA="${REF_BASE}.faa"

if [[ -n "$REF_GFF" && -f "$REF_FAA" ]]; then
    log "Official reference annotation detected:"
    log "    GFF: $REF_GFF"
    log "    FAA: $REF_FAA"
    log "    -> VV####/VVA#### nomenclature will be recovered after Bakta"
else
    warn "Did not find BOTH ${REF_BASE}.faa and ${REF_BASE}.gff"
    warn "Genes will keep the Bakta locus_tag (no VV/VVA nomenclature)"
fi

mkdir -p "$OUT"/{bakta,panaroo,snps,amrfinder,summary,logs,tmp}

# genome list: those in the directory + the reference
declare -A GPATH
shopt -s nullglob
for f in "$GDIR"/*.fna "$GDIR"/*.fa "$GDIR"/*.fasta; do
    gid="$(basename "$f" | sed 's/\.[^.]*$//')"
    GPATH[$gid]="$(realpath "$f")"
done
shopt -u nullglob
GPATH[$REF_ID]="$REF"
ALL=("${!GPATH[@]}")
log "Genomes: ${#ALL[@]} (reference: $REF_ID)"
log "Threads: $THREADS | Parallel jobs: $PARALLEL | Threads/job: $CPUS_PER_JOB"
log "min-cov=$MIN_COV | min-strains=$MIN_STRAINS | core-threshold=$CORE_THRESHOLD"

# =============================================================================
# 1. Homogeneous annotation with Bakta (ALL strains, incl. reference)
#    Homogeneous on purpose: Panaroo needs all GFFs from the same gene caller.
#    Official nomenclature is recovered afterwards (step 2).
# =============================================================================
step "1/6 · Bakta annotation (${#ALL[@]} genomes)"

: > "$OUT/tmp/bakta_jobs.tsv"
for gid in "${ALL[@]}"; do
    echo -e "${gid}\t${GPATH[$gid]}" >> "$OUT/tmp/bakta_jobs.tsv"
done

run_bakta() {
    IFS=$'\t' read -r gid fpath <<< "$1"
    if [[ -f "$OUT/bakta/$gid/$gid.gff3" ]]; then
        echo "[!] Bakta $gid already exists — skipping"; return 0
    fi
    mkdir -p "$OUT/bakta/$gid"
    if mamba run -n "$ENV_ANNOT" bakta "$fpath" \
            --db "$BAKTA_DB" \
            --output "$OUT/bakta/$gid" --prefix "$gid" \
            --locus-tag "$gid" --genus "$GENUS" --species "$SPECIES" \
            --threads "$CPUS_PER_JOB" --force --skip-plot --skip-crispr \
            > "$OUT/logs/bakta_${gid}.log" 2>&1; then
        echo "[✔] Bakta $gid OK"
    else
        echo "[✘] Bakta $gid FAILED — see $OUT/logs/bakta_${gid}.log"
    fi
}
export -f run_bakta
export OUT ENV_ANNOT BAKTA_DB GENUS SPECIES CPUS_PER_JOB

xargs -d '\n' -P "$PARALLEL" -I{} bash -c 'run_bakta "$@"' _ {} < "$OUT/tmp/bakta_jobs.tsv"

# Bakta emits GFF3 (with FASTA appended), GBFF, FAA. Snippy needs GenBank:
REF_GBK="$OUT/bakta/$REF_ID/$REF_ID.gbff"
[[ -f "$REF_GBK" ]] || err "Bakta did not produce $REF_GBK; Snippy needs it."

# =============================================================================
# 2. Name map: Bakta locus_tag -> VV####/VVA#### + gene symbol
#    protein-protein BLAST between Bakta's reference .faa and the official
#    YJ016 .faa. Panaroo never sees the official annotation.
# =============================================================================
step "2/6 · Name map (Bakta -> official $REF_ID nomenclature)"

RENAME_MAP="$OUT/summary/gene_rename_map.tsv"
if [[ -n "$REF_GFF" && -f "$REF_FAA" ]]; then
    if [[ ! -f "$RENAME_MAP" ]]; then
        crun "$ENV_ANNOT" python3 "$SCRIPT_DIR/make_rename_map.py" \
            --bakta-faa   "$OUT/bakta/$REF_ID/$REF_ID.faa" \
            --official-faa "$REF_FAA" \
            --official-gff "$REF_GFF" \
            --out "$RENAME_MAP" \
            --min-identity 90 --min-coverage 90 --threads "$THREADS" \
            2>&1 | tee "$OUT/logs/rename_map.log"
    else
        warn "Name map already exists — skipping"
    fi
else
    : > "$RENAME_MAP"
fi

# =============================================================================
# 3. Pangenome / core genome with Panaroo
#    Used for: core genome definition + core_gene_alignment.aln
#    (-> IQ-TREE tree + pyseer kinship). NOT used to call SNPs.
# =============================================================================
step "3/6 · Panaroo (core genome + core alignment)"

if [[ ! -f "$OUT/panaroo/gene_presence_absence.csv" ]]; then
    : > "$OUT/tmp/gff_list.txt"
    for gid in "${ALL[@]}"; do
        echo "$OUT/bakta/$gid/$gid.gff3" >> "$OUT/tmp/gff_list.txt"
    done
    crun "$ENV_CORE" panaroo -i "$OUT/tmp/gff_list.txt" -o "$OUT/panaroo" \
        --clean-mode strict --remove-invalid-genes \
        -a core --aligner mafft --core_threshold "$CORE_THRESHOLD" \
        -t "$THREADS" 2>&1 | tee "$OUT/logs/panaroo.log"
else
    warn "Panaroo already exists — skipping"
fi

# Total gene families in the pangenome (NOT the core count; core is smaller)
N_FAM=$(awk -F',' 'NR>1' "$OUT/panaroo/gene_presence_absence.csv" 2>/dev/null | wc -l || echo 0)
log "Gene families in the pangenome: $N_FAM"
[[ -f "$OUT/panaroo/core_gene_alignment.aln" ]] && \
    log "Core alignment: $OUT/panaroo/core_gene_alignment.aln (for tree + kinship)"

# =============================================================================
# 4. Variants vs reference with Snippy (genome-wide, effect-annotated)
# =============================================================================
step "4/6 · Snippy (reference: $REF_ID)"

: > "$OUT/tmp/snippy_jobs.tsv"
for gid in "${ALL[@]}"; do
    [[ "$gid" == "$REF_ID" ]] && continue
    echo -e "${gid}\t${GPATH[$gid]}" >> "$OUT/tmp/snippy_jobs.tsv"
done

run_snippy() {
    IFS=$'\t' read -r gid fpath <<< "$1"
    local DEST="$OUT/snps/${gid}.snps.tab"
    if [[ -f "$DEST" ]]; then
        echo "[!] Snippy $gid already exists — skipping"; return 0
    fi
    if mamba run -n "$ENV_SNIPPY" snippy \
            --outdir "$OUT/tmp/snippy_$gid" --ref "$REF_GBK" --ctgs "$fpath" \
            --mincov "$MIN_COV" --cpus "$CPUS_PER_JOB" --force \
            > "$OUT/logs/snippy_${gid}.log" 2>&1 && \
       [[ -f "$OUT/tmp/snippy_$gid/snps.tab" ]]; then
        cp "$OUT/tmp/snippy_$gid/snps.tab" "$DEST"
        echo "[✔] Snippy $gid OK ($(( $(wc -l < "$DEST") - 1 )) variants)"
    else
        echo "[✘] Snippy $gid FAILED — see $OUT/logs/snippy_${gid}.log"
    fi
}
export -f run_snippy
export OUT ENV_SNIPPY REF_GBK MIN_COV CPUS_PER_JOB

xargs -d '\n' -P "$PARALLEL" -I{} bash -c 'run_snippy "$@"' _ {} < "$OUT/tmp/snippy_jobs.tsv"

N_OK=$(ls -1 "$OUT/snps"/*.snps.tab 2>/dev/null | wc -l)
log "Snippy: $N_OK genomes with variants called"
[[ "$N_OK" -eq 0 ]] && err "Snippy produced no tables. Check $OUT/logs/snippy_*.log"

# =============================================================================
# 5. AMRFinderPlus per genome (independent acquired-resistance check)
#    Complements ResFinder for the manuscript (expected: 0 acquired genes).
# =============================================================================
step "5/6 · AMRFinderPlus (acquired resistance genes)"

AMR_SUMMARY="$OUT/amrfinder/amrfinder_summary.tsv"
if [[ ! -f "$AMR_SUMMARY" ]]; then
    echo -e "strain\tn_amr_genes\tgenes" > "$AMR_SUMMARY"
    for gid in "${ALL[@]}"; do
        faa="$OUT/bakta/$gid/$gid.faa"
        [[ -f "$faa" ]] || { warn "no faa for $gid, skipping AMR"; continue; }
        out="$OUT/amrfinder/${gid}.amr.tsv"
        if [[ ! -f "$out" ]]; then
            crun "$ENV_ANNOT" amrfinder -p "$faa" --plus \
                > "$out" 2> "$OUT/logs/amrfinder_${gid}.log" || \
                warn "AMRFinder failed for $gid"
        fi
        n=$(($(wc -l < "$out") - 1)); [[ "$n" -lt 0 ]] && n=0
        genes=$(awk -F'\t' 'NR>1{print $6}' "$out" | paste -sd';' -)
        echo -e "${gid}\t${n}\t${genes}" >> "$AMR_SUMMARY"
    done
    log "AMRFinderPlus summary: $AMR_SUMMARY"
else
    warn "AMRFinder summary already exists — skipping"
fi

# =============================================================================
# 6. Filtered variant matrix + linkage blocks
#      CODING     -> non-synonymous only
#      NONCODING  -> all
# =============================================================================
step "6/6 · Variant matrix (coding: non-syn | noncoding: all)"

crun "$ENV_CORE" python3 "$SCRIPT_DIR/build_snp_matrix.py" \
    --snps-dir "$OUT/snps" \
    --rename-map "$RENAME_MAP" \
    --prefix "$OUT/summary/variants" \
    --min-strains "$MIN_STRAINS" \
    2>&1 | tee "$OUT/logs/variant_matrix.log"

echo ""
echo -e "${BLD}${GRN}══════════════════════════════════════════════${RST}"
echo -e "${BLD}  PIPELINE COMPLETE${RST}"
echo -e "${BLD}${GRN}══════════════════════════════════════════════${RST}"
echo -e "  Variant matrix   : ${BLD}$OUT/summary/variants_matrix.tsv${RST}"
echo -e "  Variant annot    : ${BLD}$OUT/summary/variants_variants.tsv${RST}"
echo -e "  Category summary : ${BLD}$OUT/summary/variants_summary.tsv${RST}"
echo -e "  Core alignment   : ${BLD}$OUT/panaroo/core_gene_alignment.aln${RST}"
echo -e "  VV/VVA map       : ${BLD}$RENAME_MAP${RST}"
echo -e "  AMRFinder summary: ${BLD}$AMR_SUMMARY${RST}"
echo ""
echo -e "  Next step: tree + kinship for pyseer"
echo -e "    ${CYN}mamba run -n $ENV_CORE iqtree3 -s $OUT/panaroo/core_gene_alignment.aln -m GTR+G -B 1000 -T AUTO${RST}"
echo -e "    ${CYN}mamba run -n vv_assoc python scripts/07_run_pyseer.py ...${RST}"
echo ""

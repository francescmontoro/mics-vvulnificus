# SNP-MIC pipeline for antibiotic resistance prediction in *Vibrio vulnificus*

This repository contains the bioinformatic pipeline described in:

> Montoro *et al.* (2025). *New pipeline for antibiotic resistance gene prediction
> using MIC and WGS: a practical approach for Vibrio vulnificus.*

## Overview

The pipeline identifies non-synonymous SNPs in the core genome of
*V. vulnificus* that are statistically associated with changes in minimum
inhibitory concentrations (MICs), using Fisher's exact test with
Benjamini-Hochberg FDR correction.

## Pipeline scripts

| Script | Description |
|--------|-------------|
| `01_convert_mics.py` | Converts raw MIC values (mg/L) to a normalised scale [0–5] |
| `02_figure_heatmap.py` | Publication-quality heatmap of MIC profiles (Fig. 1) |
| `03_snp_association.py` | Fisher's exact test for each SNP–antibiotic pair + BH FDR |
| `04_gene_level_summary.py` | Collapses SNP associations to gene level; family-level analysis |
| `05_figure_volcano.py` | Volcano plot of gene-level associations (Fig. 2, supplementary) |
| `05_figure_dotplot.py` | Dot plot of top significant genes (Fig. 2) |
| `06_figure_phylo.py` | Phylogenetic tree + resistance heatmap (Fig. 3) |
| `07_resfinder_comparison.py` | ResFinder 4.x comparison (Fig. 4) |

## Requirements

```
Python >= 3.9
pandas >= 1.4
scipy >= 1.7
numpy >= 1.20
matplotlib >= 3.5
biopython >= 1.79
adjustText
resfinder >= 4.7
BLAST+ >= 2.12 (blastn must be in PATH)
MAFFT >= 7
PhyML >= 3.0
```

Install Python dependencies:
```bash
pip install pandas scipy numpy matplotlib biopython adjustText resfinder
```

## Input data

| File | Description |
|------|-------------|
| `data/raw/MICs_all_Vv.xlsx` | Raw MIC data (116 strains × 41 antibiotics) |
| `data/raw/allnsSNP.tsv` | Non-synonymous SNPs from core-genome MSA |
| `data/raw/genomes/*.fasta` | 30 *V. vulnificus* genome assemblies |
| `data/raw/lineage_manual_mapping.tsv` | Lineage assignments (L1–L4) per strain |
| `data/raw/phylo_tree.nwk` | Newick phylogenetic tree |
| `data/raw/resfinder_db/` | ResFinder acquired resistance database |
| `data/raw/pointfinder_db/` | PointFinder chromosomal mutation database |

> MIC data belong to the Department of Microbiology, University of Valencia,
> and are available from the corresponding author upon reasonable request.

## Usage

Run scripts sequentially from the `scripts/` directory:

```bash
cd scripts/
python3 01_convert_mics.py
python3 02_figure_heatmap.py
python3 03_snp_association.py    # ~10 min
python3 04_gene_level_summary.py # ~10 min
python3 05_figure_dotplot.py
python3 06_figure_phylo.py
python3 07_resfinder_comparison.py
```

## Output structure

```
PAPER_V2/
├── data/
│   ├── raw/          # Input data
│   └── processed/    # Intermediate results (TSV files)
├── results/          # Final tables (gene-level associations, summaries)
├── figures/          # Publication figures (PDF + PNG)
└── scripts/          # Pipeline scripts (this repository)
```

## Key results

- **ResFinder**: 0/30 strains with detected acquired resistance genes
- **Our pipeline**: 4,316 significant SNP–antibiotic associations (FDR < 0.05)
  - GEN2 (gentamicin): 1,062 significant genes
  - STZ (sulphathiazole): 547 significant genes
  - Sulfonamide family: 1,528 significant genes (top hit: *zwf*, *vexH*)

## Authors

- Francesc Josep Montoro Salvador (pipeline development)
- Francisco José Roig Molina (supervision)
- Carmen Amaro (MIC data, supervision)

## License

MIT License

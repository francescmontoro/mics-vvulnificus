# mics-vvulnificus

Code and derived data for:

> Montoro Salvador FJ, Amaro C, Roig Molina FJ. **Evaluating the interpretability of bacterial resistance GWAS in small collections: a case study in *Vibrio vulnificus*.** Submitted to *Microbial Genomics*.

The repository contains every script used for the analyses, tables and figures of the article, the conda environment specifications, the patches needed to run two packages under current SciPy, and the small result tables from which the numbers in the text are taken. Large inputs (assemblies, databases) are not stored here; [Data](#data) says where to obtain them.

## Contents

```
environments/   conda environment specifications (one .yml per environment)
patches/        compatibility patch for Scoary 1.6.16 and the script that applies it
scripts/        analysis pipeline, numbered in the order in which it is run
data/           strain list with accessions; pointers to the external datasets
results/        result tables that support the tables and figures of the article
```

## Data

| Dataset | Where |
|---|---|
| MIC values, 116 isolates (Sensititre GN2F and AVIAN1F) | Zenodo, [doi:10.5281/zenodo.20679806](https://doi.org/10.5281/zenodo.20679806) |
| Filtered assemblies of the 29 strains analysed | Zenodo, [doi:10.5281/zenodo.23639878](https://doi.org/10.5281/zenodo.23639878) |
| Deposited assemblies and reference genome YJ016 | NCBI; accessions in `data/strains.tsv` (Table S1 of the article) |

Place the MIC workbook at `data/mic/MICs_all_Vv.xlsx` and the assemblies in `data/genomes/` before running the pipeline. File names in the analysis use the identifiers in the `file_id` column of `data/strains.tsv`; three of them differ from the published strain designations (94-9-112 = 94-8-112, Ra_3 = Ra3, RIU1 = Riu1), and figures display the published names.

## Software

Analyses were run under Windows Subsystem for Linux 2 (Ubuntu 24.04) on a 12-core machine with 13 GB of RAM, in conda environments kept separate to avoid dependency conflicts:

| Environment | Used for |
|---|---|
| `vv_annot` | Bakta, AMRFinderPlus, BLAST+, contamination screening |
| `vv_core` | Panaroo, IQ-TREE, SNP-sites |
| `vv_assoc` | pyseer, Scoary, statistics and figures |
| `snippy` | Snippy and its dependencies |
| `resfinder` | ResFinder 4.7.2 (minimal environment; see `environments/README.md`) |

```bash
for f in environments/*.yml; do mamba env create -f "$f"; done
bash patches/apply_patches.sh          # patches Scoary inside vv_assoc
```

Exact versions, and the commit of every package installed from a development branch, are listed in Table S2 of the article. Two points are worth knowing in advance:

* **pyseer** and **Panaroo** are installed from GitHub at fixed commits (pyseer `43cfb30fdd0cf553086c4bb492a80bb8db830bdf`, Panaroo `e96b4978327536256b5c663f741138fc7aabed64`); for pyseer this is because the bioconda release is incompatible with SciPy ≥ 1.12. Its bioconda package also lacks the helper `phylogeny_distance.py`; the kinship calculation it performs is implemented directly in `scripts/07_run_pyseer.py`, so no patch is required. pyseer is run single-threaded (its multiprocessing pool deadlocks under WSL2 with a threaded BLAS; results are identical).
* **Scoary 1.6.16** calls `scipy.stats.binom_test`, removed in SciPy 1.12. `patches/` contains a shim that maps it to `scipy.stats.binomtest`; the p-values returned are identical. Note that Scoary reports the two-sided binomial p-value for the pairwise comparisons.

Bakta was run with the light database v6.0 ([doi:10.5281/zenodo.14916843](https://doi.org/10.5281/zenodo.14916843)) because the full database does not fit in 13 GB of memory.

## Pipeline

Each script documents its inputs and options in its header (`python3 scripts/<script> --help`). Steps are listed in the order in which they were run.

| Step | Script | What it does | Article |
|---|---|---|---|
| 1 | `13_contamination_check.py` | Length-weighted GC profile and 16S assignment of every assembly | Table 1 |
| 2 | `14_decontaminate.py` | Uniform contig filter (no homology to YJ016 **and** atypical GC; < 200 bp) | Methods, Table 1 |
| 3 | `16_compare_assemblies.py` | Filtered assemblies against the deposited NCBI versions | Table S1, note b |
| 4 | `run_vv_variants.sh` (calls `build_snp_matrix.py`, `make_rename_map.py`) | Bakta annotation, Panaroo pan-genome, Snippy variants against YJ016, variant matrix | Methods |
| 5 | `10_assembly_qc.py` | Assembly statistics | Table S1 |
| 6 | `11_run_amrfinder.py`, `21_run_resfinder.py` | Acquired resistance determinants, with explicit failure handling | Table 2 |
| 7 | `12_genotype_phenotype_check.py` | Determinants against MICs | Table 2, Fig. 3 |
| 8 | `19_map_resistance_to_pangenome.py` | Determinants located in the pan-genome through locus tags | Results, short contigs |
| 9 | `00_signal_check.py` | Independent groups of resistant strains on the tree | Results |
| 10 | `17_define_phenotypes.py` | Distribution-based phenotypes and sensitivity cut-off | Methods, Fig. 2 |
| 11 | `07_run_pyseer.py` | IQ-TREE phylogeny, kinship matrix and pyseer LMM on core-genome variants | Results |
| 12 | `08_run_scoary.py` | Fisher + BH and pairwise comparisons over the tree | Fig. 4 |
| 13 | `09_linkage_blocks.py`, `15_block_level_fdr.py` | Variants collapsed to presence/absence patterns; correction per pattern | Results |
| 14 | `18_pangenome_association.py` | Gene presence/absence association and exact re-testing | Table 3 |
| 15 | `22_calibration_simulation.py` | 10 000 label permutations at r = 4, 6, 8; LMM, Fisher control and Fitch emergences | Results |
| 16 | `23_figure_calibration.py` | Calibration summary and figure | Fig. 7 |
| 17 | `20_make_figures.py` | Remaining data figures | Figs 2–6 |

Fig. 1 (diagnostic workflow) is a schematic and is not generated by code.

Example for the calibration step, as run for the article:

```bash
mamba run -n vv_assoc python3 scripts/22_calibration_simulation.py \
    --rtab     results/pangenome_assoc/pangenome.Rtab \
    --kinship  results/pyseer/phylogeny_K.tsv \
    --tree     results/pyseer/core.treefile \
    --outdir   results/calibration \
    --n-perm   10000 --r-values 4 6 8 --env vv_assoc
python3 scripts/23_figure_calibration.py --calib results/calibration --outdir figures
```

## Results included

`results/` holds the tables behind the numbers reported in the article (see `results/README.md`). Intermediate files that are large or fully regenerated by the pipeline (annotations, alignments, per-variant pyseer output) are not included.

## Citation

If you use this code, please cite the article above and the archived release of this repository ([doi:10.5281/zenodo.23639674](https://doi.org/10.5281/zenodo.23639674), which always resolves to the latest version).

## Licence

Code: MIT (see `LICENSE`). Data files in `results/` and `data/`: CC BY 4.0.

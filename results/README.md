# Results

Tables behind the numbers, tables and figures of the article. Strain identifiers follow the `file_id` column of `data/strains.tsv`.

| Folder | Content | Article |
|---|---|---|
| `contamination/` | 16S assignment and length-weighted GC profile of the unfiltered assemblies | Table 1 |
| `decontamination/` | Per-contig decision of the uniform filter (kept/removed, with reason) for every assembly | Methods, Table 1, Results |
| `contamination_postclean/` | The same screen repeated on the filtered assemblies | Methods |
| `assembly_comparison/` | Filtered assemblies against the deposited NCBI versions | Table S1, note b |
| `qc/` | Assembly statistics | Table S1 |
| `amrfinder/`, `resfinder/` | Determinants detected by each tool on the filtered assemblies | Table 2 |
| `genotype_phenotype/` | Determinants against MICs | Table 2, Fig. 3 |
| `00_signal_check/` | Independent groups of resistant strains on the tree | Results |
| `pyseer/` | Core-genome phylogeny (IQ-TREE), kinship matrix and phenotype files | Methods, Fig. 6 |
| `summary/` | Filtered variant matrix (19 525 variants × 29 strains) and summary | Results, Fig. 6 |
| `pyseer_blocks/` | LMM results collapsed to presence/absence patterns, with per-pattern correction | Results |
| `scoary/` | Fisher + BH and pairwise-comparison results | Fig. 4 |
| `pangenome_assoc/` | Gene presence/absence matrix, LMM results and exact re-tests | Table 3 |
| `calibration/` | 10 000 permutations at r = 4, 6 and 8: smallest LMM p-value, smallest Fisher p-value, genomic inflation, heritability and Fitch emergences per permutation | Results, Fig. 7 |

A note on `calibration/calibration_summary.tsv`: at r = 4 its `pct_fisher_below_floor` column reads 0.14. These are 14 permutations in which a variant matches the permuted phenotype perfectly, so that Fisher's p-value equals the floor 1/C(29, 4); the tables store p-values to seven significant digits and the rounded value falls marginally below the floor. Fisher's exact test never returns a value below the floor; `scripts/23_figure_calibration.py` applies a tolerance for this rounding.

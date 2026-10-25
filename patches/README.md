# Patches

## Scoary 1.6.16 — `scoary-1.6.16_scipy-binomtest.patch`

`scoary/methods.py` calls `scipy.stats.binom_test`, which was removed in SciPy 1.12. The patch adds a shim that maps the call to `scipy.stats.binomtest(k, n, p).pvalue`. The returned p-value is identical (two-sided, as in the original code).

Apply it with

```bash
bash patches/apply_patches.sh            # default environment: vv_assoc
bash patches/apply_patches.sh my_env     # or name another environment
```

The script locates the installed Scoary package inside the environment, checks that it is version 1.6.16 and applies the patch once (it does nothing if the patch is already in place).

## pyseer — no patch required

The bioconda package of pyseer does not include the helper script `phylogeny_distance.py`, used to derive a kinship matrix from a tree. The calculation (shared branch length from the root to the most recent common ancestor of each pair of tips, with the root-to-tip distance on the diagonal) is implemented directly in `scripts/07_run_pyseer.py`.

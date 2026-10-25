# Conda environments

One specification per environment used in the analysis. Create them with

```bash
for f in environments/*.yml; do mamba env create -f "$f"; done
```

The files were exported from the environments in which the analysis was run:

```bash
for e in vv_annot vv_core vv_assoc snippy; do   # resfinder.yml was written by hand, see below
    conda env export -n "$e" --no-builds | grep -v '^prefix:' > environments/"$e".yml
done
```

Packages installed with pip appear under the `pip:` section of the corresponding file. Two of them were installed from GitHub and are pinned to the exact commit used: pyseer 1.4.2.dev0 (`43cfb30fdd0cf553086c4bb492a80bb8db830bdf`) and Panaroo 1.8.0 (`e96b4978327536256b5c663f741138fc7aabed64`).

ResFinder 4.7.2 was installed with pip in the base environment of the analysis machine, which also holds unrelated software. Rather than export that environment, `resfinder.yml` defines a minimal one with the same ResFinder version and BLAST+ 2.17.0. Snippy needs `gsl=2.6` to resolve a shared-library conflict in BCFtools; this is pinned in `snippy.yml`.

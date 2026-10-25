#!/usr/bin/env python3
"""
22_calibration_simulation.py
============================
Assesses whether the linear mixed model is calibrated at this sample size, by
running it against phenotypes that carry no association by construction.

RATIONALE
---------
The manuscript observes that one pan-genome cluster received a mixed-model
p-value (2.05e-05) below the smallest value Fisher's exact test can return for
the same design (1/C(29,4) = 4.2e-05). Taken alone that observation is open to
an obvious objection: the two tests have different null distributions, so the
mixed model is not logically required to respect the exact-test floor.

This script closes that gap empirically. The resistance labels are permuted
among the 29 strains, holding the number of resistant strains fixed, so that no
variant can be associated with the phenotype other than by chance. Genotypes,
kinship and every analysis parameter are left untouched. Under a calibrated
test, the probability that the smallest p-value across the genome falls below
any threshold t is at most (number of independent tests) x t; for t equal to
the exact-test floor this is a small number. A materially higher observed rate
is direct evidence that the asymptotic reference distribution is
anticonservative at this scale.

Three phenotype prevalences are examined (r = 4, 6, 8 resistant of 29),
spanning the primary cut-off, the sensitivity cut-off and a more balanced
design, since the severity of the miscalibration is expected to depend on the
number of cases.

WHAT IS REPORTED
----------------
  * the proportion of permutations in which the minimum p-value across the
    genome falls below the exact-test floor 1/C(29, r)
  * the same proportion for Fisher's exact test on the identical permutations,
    which is zero by construction and serves as a control on the pipeline
  * the genomic inflation factor lambda for each permutation, the conventional
    calibration statistic
  * the distribution of minimum p-values, for a QQ-style figure
  * results stratified by the heritability the model estimates in each
    permutation, since the collapse of the random effect (h2 -> 0) is the
    proposed mechanism

RUNTIME
-------
Each permutation is one pyseer invocation. On the deduplicated pan-genome
(636 patterns) a run takes roughly 5-10 s including interpreter start-up, so
1000 permutations per prevalence is about 5 h for the three prevalences. Use
--n-perm to scale. Results are appended after every permutation, so the script
can be interrupted and resumed with --resume.

Usage:
  python3 22_calibration_simulation.py \\
      --rtab     results/pangenome_assoc/pangenome.Rtab \\
      --kinship  results/pyseer/phylogeny_K.tsv \\
      --outdir   results/calibration \\
      --n-perm   1000 \\
      --r-values 4 6 8 \\
      --env      vv_assoc
"""

import argparse
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact


def floor_p(n, r):
    """Smallest p-value Fisher's exact test can return for these margins."""
    return 1.0 / math.comb(n, r)


def lambda_gc(pvals):
    """
    Genomic inflation factor: ratio of the median observed chi-square statistic
    to its expectation under the null. Values above 1 indicate inflation.
    """
    p = np.asarray([x for x in pvals if 0 < x <= 1])
    if len(p) < 10:
        return float("nan")
    obs = chi2.isf(p, df=1)
    return float(np.median(obs) / chi2.ppf(0.5, df=1))


def fitch_changes(tree, tips, labels):
    """
    Number of state changes required on the tree under Fitch parsimony, that is
    the minimum number of independent emergences of the trait.

    This is deliberately independent of the association model. The heritability
    the model reports is an output of the very procedure under scrutiny, so
    relating the failure rate to it would be circular; the parsimony count is a
    property of the tree and the labelling alone, and is the same quantity used
    for the a priori signal check.
    """
    state = {t: {labels[i]} for i, t in enumerate(tips)}
    changes = 0

    def post(node):
        nonlocal changes
        if node.is_terminal():
            return state[node.name]
        sets = [post(c) for c in node.clades]
        inter = set.intersection(*sets)
        if inter:
            return inter
        changes += 1
        return set.union(*sets)

    post(tree.root)
    return changes


def load_rtab(path, samples=None):
    df = pd.read_csv(path, sep="\t", index_col=0)
    if samples is not None:
        missing = [s for s in samples if s not in df.columns]
        if missing:
            sys.exit(f"ERROR: samples absent from the Rtab: {missing}")
        df = df[samples]
    # one representative per presence/absence pattern: variants sharing a
    # pattern give identical statistics and would only slow the loop.
    # Patterns that are constant across the retained samples carry no
    # information and are dropped.
    df = df.loc[(df.sum(axis=1) > 0) & (df.sum(axis=1) < df.shape[1])]
    key = df.apply(lambda row: tuple(row.values), axis=1)
    return df.loc[~key.duplicated()]


def run_pyseer(rtab, pheno, kinship, env, workdir):
    out = Path(workdir) / "res.tsv"
    cmd = ["pyseer", "--lmm",
           "--phenotypes", str(pheno), "--phenotype-column", "phen",
           "--pres", str(rtab), "--similarity", str(kinship),
           "--cpu", "1"]
    if env:
        cmd = ["mamba", "run", "-n", env] + cmd
    with open(out, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        return None, None
    # heritability is printed to stderr
    h2 = float("nan")
    for line in proc.stderr.decode(errors="replace").splitlines():
        if line.strip().startswith("h^2"):
            try:
                h2 = float(line.split("=")[1])
            except (IndexError, ValueError):
                pass
    try:
        res = pd.read_csv(out, sep="\t")
    except Exception:
        return None, None
    pcol = "lrt-pvalue" if "lrt-pvalue" in res.columns else res.columns[3]
    return pd.to_numeric(res[pcol], errors="coerce").dropna().values, h2


def fisher_min(mat, labels):
    """Smallest Fisher exact p across all variants for this labelling."""
    best = 1.0
    y = labels.astype(bool)
    for i in range(mat.shape[0]):
        g = mat[i].astype(bool)
        a = int((g & y).sum()); b = int((~g & y).sum())
        c = int((g & ~y).sum()); d = int((~g & ~y).sum())
        if a + c == 0 or b + d == 0:
            continue
        _, p = fisher_exact([[a, c], [b, d]])
        if p < best:
            best = p
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rtab", required=True)
    ap.add_argument("--kinship", required=True)
    ap.add_argument("--tree", default="",
                    help="newick tree; enables the parsimony count of "
                         "independent emergences for each permuted "
                         "phenotype, an external measure of how far the "
                         "labelling is confounded with lineage")
    ap.add_argument("--phenotype", default="",
                    help="an existing phenotype file; its sample column fixes "
                         "which strains enter the simulation, so that the "
                         "design matches the association analysis exactly. "
                         "Without it the Rtab may contribute strains, such as "
                         "the mapping reference, that the real analysis "
                         "excluded, and the floor would be computed for the "
                         "wrong n.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--r-values", type=int, nargs="+", default=[4, 6, 8])
    ap.add_argument("--env", default="vv_assoc")
    ap.add_argument("--seed", type=int, default=20260829)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--skip-fisher", action="store_true",
                    help="omit the Fisher control, which is the slow part")
    a = ap.parse_args()

    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
              "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[v] = "1"

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)

    keep = None
    if a.phenotype:
        ph0 = pd.read_csv(a.phenotype, sep="\t")
        keep = [str(x) for x in ph0[ph0.columns[0]]]

    rtab = load_rtab(a.rtab, keep)
    samples = list(rtab.columns)

    tree = None
    if a.tree:
        from Bio import Phylo
        tree = Phylo.read(a.tree, "newick")
        tips = {t.name for t in tree.get_terminals()}
        missing = [s for s in samples if s not in tips]
        if missing:
            sys.exit(f"ERROR: tips absent from the tree: {missing}")
    n = len(samples)
    mat = rtab.values
    print("=" * 74)
    print("  Calibration of the mixed model under permuted phenotypes")
    print("=" * 74)
    print(f"  strains              : {n}"
          + ("  (from the phenotype file)" if keep else
             "  (all columns of the Rtab)"))
    print(f"  distinct patterns    : {mat.shape[0]}")
    print(f"  permutations each    : {a.n_perm}")
    print(f"  prevalences examined : {a.r_values}\n")

    # deduplicated Rtab written once
    tmp = Path(tempfile.mkdtemp(prefix="calib_"))
    rtab_path = tmp / "patterns.Rtab"
    rtab.to_csv(rtab_path, sep="\t")

    rng = np.random.default_rng(a.seed)

    for r in a.r_values:
        fl = floor_p(n, r)
        rec_path = outdir / f"perm_r{r}.tsv"
        done = 0
        if a.resume and rec_path.exists():
            done = sum(1 for _ in open(rec_path)) - 1
            print(f"  r={r}: resuming after {done} permutations")
        else:
            with open(rec_path, "w") as fh:
                fh.write("perm\tmin_p_lmm\tn_below_floor\tlambda\th2\t"
                         "min_p_fisher\tfitch\tresistant\n")

        print(f"\n  ── r = {r} of {n}   floor 1/C({n},{r}) = {fl:.2e} ──")

        for i in range(done, a.n_perm):
            y = np.zeros(n, dtype=int)
            y[rng.choice(n, size=r, replace=False)] = 1

            ph = tmp / "pheno.tsv"
            pd.DataFrame({"samples": samples, "phen": y}).to_csv(
                ph, sep="\t", index=False)

            pv, h2 = run_pyseer(rtab_path, ph, a.kinship, a.env, tmp)
            if pv is None or len(pv) == 0:
                continue

            minp = float(np.min(pv))
            below = int((pv < fl).sum())
            lam = lambda_gc(pv)
            fmin = float("nan") if a.skip_fisher else fisher_min(mat, y)

            fit = (fitch_changes(tree, samples, y) if tree is not None
                   else -1)
            who = ";".join(s for s, v in zip(samples, y) if v)
            with open(rec_path, "a") as fh:
                fh.write(f"{i}\t{minp:.6e}\t{below}\t{lam:.4f}\t{h2:.4f}\t"
                         f"{fmin:.6e}\t{fit}\t{who}\n")

            if (i + 1) % 25 == 0 or i + 1 == a.n_perm:
                d = pd.read_csv(rec_path, sep="\t")
                rate = (d["min_p_lmm"] < fl).mean()
                print(f"    {i+1:>5}/{a.n_perm}   below floor so far: "
                      f"{rate*100:5.1f}%   median lambda "
                      f"{d['lambda'].median():.2f}", flush=True)

    # ── summary ─────────────────────────────────────────────────────────────
    print("\n" + "=" * 74)
    rows = []
    for r in a.r_values:
        p = outdir / f"perm_r{r}.tsv"
        if not p.exists():
            continue
        d = pd.read_csv(p, sep="\t")
        if not len(d):
            continue
        fl = floor_p(n, r)
        rows.append({
            "r": r,
            "permutations": len(d),
            "floor": fl,
            "expected_pct": round(100 * min(1.0, mat.shape[0] * fl), 2),
            "pct_lmm_below_floor": round(100 * (d["min_p_lmm"] < fl).mean(), 2),
            "pct_fisher_below_floor": round(
                100 * (d["min_p_fisher"] < fl).mean(), 2)
                if d["min_p_fisher"].notna().any() else None,
            "median_lambda": round(d["lambda"].median(), 3),
            "median_h2": round(d["h2"].median(), 3),
            "pct_h2_zero": round(100 * (d["h2"] <= 1e-6).mean(), 1),
            "median_fitch": (round(d["fitch"].median(), 1)
                             if "fitch" in d and (d["fitch"] >= 0).any()
                             else None),
        })
    s = pd.DataFrame(rows)
    s.to_csv(outdir / "calibration_summary.tsv", sep="\t", index=False)
    print(s.to_string(index=False))
    print(f"\n  wrote {outdir/'calibration_summary.tsv'}")

    print("""
  HOW TO READ THIS
    pct_lmm_below_floor is the proportion of permutations, each carrying no
    association by construction, in which the mixed model returned at least one
    p-value below the smallest value the exact test can produce for the design.

    expected_pct is the upper bound on that proportion for any calibrated test:
    the union bound m x floor, where m is the number of distinct patterns. It
    holds whatever the correlation between patterns, so an observed rate above
    it cannot be attributed to linkage.

    pct_fisher_below_floor is the same quantity for the exact test and is zero
    by construction; it confirms that the permutation machinery introduces no
    signal of its own.

    median_lambda is the genomic inflation factor. A value near 1 alongside a
    high floor-crossing rate means the bulk of the null distribution is well
    behaved and only the extreme tail is not, so the conventional inflation
    diagnostic does not detect this failure mode.

    pct_h2_zero records how often the random effect collapsed. If the rate of
    floor crossings tracks this column, the collapse is the mechanism.
""")


if __name__ == "__main__":
    main()

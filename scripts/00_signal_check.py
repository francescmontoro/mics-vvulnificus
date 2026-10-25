#!/usr/bin/env python3
"""
00_signal_check.py
==================
Fran's step 1: before rebuilding anything, determine whether there is any
recoverable genotype-phenotype signal given the phylogenetic structure.

The logic (from Fran's comment #55):
    "Look at how many INDEPENDENT clades the 8 STZ-R and 14 GEN2-R strains
     fall into. If they fall in 1-2 clades, there is no recoverable signal
     with any correction. If they fall in 4+, there is."

Why this matters:
    If all resistant strains belong to a single clade, then resistance is
    perfectly confounded with lineage: every variant that differs between
    that clade and the rest will look 'associated', and no correction can
    separate a true causal variant from a passenger. Independent emergence
    of resistance in multiple clades is what makes association testable.

This script quantifies "independent clades" in three complementary ways so
the result does not depend on eyeballing the tree:

  (A) Monophyletic resistant clusters: greedily find maximal clades whose
      tips are ALL resistant. The number of such clusters needed to cover
      all resistant strains = number of independent emergences (lower bound).

  (B) Number of inferred state changes (parsimony): using Fitch's algorithm,
      count how many times resistance is gained/lost along the tree. More
      changes = more independent evolutionary events = more testable signal.

  (C) Phylogenetic clustering test: compare the mean pairwise patristic
      distance among resistant strains vs. a null of randomly labelled tips
      (1000 permutations). If resistant strains are NOT significantly more
      clustered than random, signal is recoverable.

Inputs:
    ../data/metadata/phylo_tree_renamed.nwk
    ../data/mic/MICs_converted.tsv

Outputs:
    ../results/00_signal_check/signal_check_summary.tsv
    ../results/00_signal_check/resistant_strains_per_antibiotic.tsv
    console verdict

Usage:
    python3 00_signal_check.py
"""

import pandas as pd
import numpy as np
from Bio import Phylo
from pathlib import Path
from itertools import combinations
import sys

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
TREE = BASE / "data" / "metadata" / "phylo_tree_renamed.nwk"
MICS = BASE / "data" / "mic" / "MICs_converted.tsv"
OUT  = BASE / "results" / "00_signal_check"
OUT.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
REFERENCE   = "CECT_5769"
MIN_DELTA   = 0.5
ANTIBIOTICS = ["GEN2", "STZ"]
N_PERM      = 1000
SEED        = 42

# Strain name mapping (tree/SNP names -> MIC names), from the v2 pipeline
SNP_TO_MIC = {
    "11028":"BT3_11028","12":"BT3_12","162":"BT3_162","32":"BT3_32",
    "94385":"94-3-85","960426-1_4C":"960426_1/4C","CIP8190":"C1P8190",
    "V4":"V_4","CECT4604":"CECT_4604","CECT4606":"CECT_4606",
    "CECT4608":"CECT_4608","CECT4865":"CECT_4865","CECT4866":"CECT_4866",
    "CECT4999":"CECT_4999","CECT5198":"CECT_5198","CECT5763":"CECT_5763",
    "CECT5769":"CECT_5769","CECT7030":"CECT_7030","CECT898":"CECT_898",
    "CG100":"CG_100","94-9-112":"94-9-112","94-9-119":"94-9-119",
    "95-8-161":"95-8-161","95-8-6":"95-8-6","95-8-7":"95-8-7","A14":"A14",
    "PD-2-51":"PD-2-51","Ra_3":"Ra_3","RIU1":"RIU1","YJ016":"YJ016",
}


def load_resistant_labels(mics, tree_tip_names):
    """Return {antibiotic: set(resistant tip names)} using tree tip naming."""
    # The tree tips are already MIC-named (phylo_tree_renamed.nwk), so map
    # MIC table index to the tip names present in the tree.
    mics = mics.copy()
    mics.index = mics.index.astype(str).str.strip()
    ref = mics.loc[REFERENCE]

    labels = {}
    for ab in ANTIBIOTICS:
        resistant = set()
        for tip in tree_tip_names:
            if tip == REFERENCE:
                continue
            if tip in mics.index and not pd.isna(mics.loc[tip, ab]):
                if mics.loc[tip, ab] >= ref[ab] + MIN_DELTA:
                    resistant.add(tip)
        labels[ab] = resistant
    return labels


def count_monophyletic_clusters(tree, resistant):
    """
    (A) Greedily cover the resistant set with maximal pure-resistant clades.
    Returns the number of independent monophyletic resistant clusters and
    the singletons (resistant tips not grouped with any other).
    """
    if not resistant:
        return 0, []

    # For every internal clade, check if ALL its tips are resistant
    pure_clades = []
    for clade in tree.get_nonterminals():
        tips = {t.name for t in clade.get_terminals()}
        if tips and tips.issubset(resistant) and len(tips) >= 2:
            pure_clades.append(tips)

    # Sort by size descending, greedily pick non-overlapping maximal clades
    pure_clades.sort(key=len, reverse=True)
    covered = set()
    clusters = []
    for tips in pure_clades:
        if not (tips & covered):     # non-overlapping
            clusters.append(tips)
            covered |= tips

    # Remaining resistant tips not in any multi-tip pure clade = singletons
    singletons = resistant - covered
    # Each singleton is its own independent emergence
    total_clusters = len(clusters) + len(singletons)
    return total_clusters, clusters, singletons


def fitch_changes(tree, resistant):
    """
    (B) Fitch parsimony: count minimum number of state changes (0/1) needed
    to explain the resistant/sensitive labelling on the tree.
    More changes => more independent evolutionary events.
    """
    # Post-order: assign state sets
    state = {}
    changes = [0]

    def postorder(clade):
        if clade.is_terminal():
            state[clade] = {1} if clade.name in resistant else {0}
            return state[clade]
        child_states = [postorder(c) for c in clade.clades]
        inter = set.intersection(*child_states)
        if inter:
            state[clade] = inter
        else:
            state[clade] = set.union(*child_states)
            changes[0] += 1
        return state[clade]
    postorder(tree.root)
    return changes[0]


def build_topological_distance_matrix(tree, all_tips):
    """
    Build a pairwise TOPOLOGICAL distance matrix (number of edges on the path
    between two tips). Works even when the tree has no branch lengths, unlike
    patristic distance. Uses the depths-to-MRCA method for efficiency.
    """
    tip_objs = {t.name: t for t in tree.get_terminals()}

    # Depth (number of edges from root) for every clade, via BFS
    depth = {}
    def set_depths(clade, d):
        depth[clade] = d
        for c in clade.clades:
            set_depths(c, d + 1)
    set_depths(tree.root, 0)

    # For topological distance between two tips:
    #   dist(a,b) = depth(a) + depth(b) - 2*depth(MRCA(a,b))
    names = list(all_tips)
    dmat = {}
    for a, b in combinations(names, 2):
        mrca = tree.common_ancestor(tip_objs[a], tip_objs[b])
        d = depth[tip_objs[a]] + depth[tip_objs[b]] - 2 * depth[mrca]
        dmat[(a, b)] = d
        dmat[(b, a)] = d
    return dmat


def topological_clustering_test(dmat, resistant, all_tips, n_perm=1000, seed=42):
    """
    (C) Permutation test using TOPOLOGICAL distances.
    Is the mean pairwise topological distance among resistant strains smaller
    than expected by chance (i.e. are they phylogenetically clustered)?

    Returns (observed_mean, null_mean, p_value).
      p small  => resistant strains strongly clustered => lineage confounding
      p large  => resistant strains spread across tree => signal recoverable
    """
    if len(resistant) < 2:
        return np.nan, np.nan, np.nan

    def mean_pairwise(group):
        gs = list(group)
        ds = [dmat[(a, b)] for a, b in combinations(gs, 2)]
        return np.mean(ds) if ds else np.nan

    obs = mean_pairwise(resistant)
    names = list(all_tips)
    rng = np.random.default_rng(seed)
    k = len(resistant)
    null = np.empty(n_perm)
    for i in range(n_perm):
        samp = set(rng.choice(names, size=k, replace=False))
        null[i] = mean_pairwise(samp)

    # one-sided: are resistant strains MORE clustered (smaller dist) than random?
    p = (np.sum(null <= obs) + 1) / (n_perm + 1)
    return obs, np.mean(null), p


def parsimony_permutation_test(tree, resistant, all_tips, n_perm=1000, seed=42):
    """
    (D) Permutation test on Fitch parsimony changes.
    Compares observed number of state changes against a null where the same
    number of resistant labels are randomly assigned to tips.

    A HIGH observed number of changes (relative to null) means resistance
    arose many times independently => strong, recoverable, non-lineage signal.
    Returns (observed_changes, null_mean_changes, p_value).
      p here = fraction of null datasets with >= observed changes.
      p LARGE  => observed changes are typical or fewer than random
                  (resistant strains clustered => confounded)
      p SMALL  => observed changes MORE than random is impossible by design;
                  instead we report how clustered: see note below.

    Interpretation we use: we report the observed vs null. If observed is
    close to the MAXIMUM possible (n_resistant) and >> null is not the point;
    the key is whether observed is far ABOVE the minimum. We compute the
    fraction of null with <= observed changes (clustered fraction).
    """
    if len(resistant) < 2:
        return np.nan, np.nan, np.nan

    obs = fitch_changes(tree, resistant)
    names = list(all_tips)
    rng = np.random.default_rng(seed)
    k = len(resistant)
    null = np.empty(n_perm, dtype=int)
    for i in range(n_perm):
        samp = set(rng.choice(names, size=k, replace=False))
        null[i] = fitch_changes(tree, samp)

    # Fraction of null with changes <= observed.
    # High value => observed is among the MOST dispersed (many changes) =>
    # strong evidence of independent emergence.
    frac_le = np.mean(null <= obs)
    return obs, float(np.mean(null)), float(frac_le)


def main():
    print("=" * 68)
    print("  SIGNAL CHECK (Fran's step 1) — is there recoverable signal?")
    print("=" * 68)

    # Load tree
    tree = Phylo.read(str(TREE), "newick")
    tip_names = [t.name for t in tree.get_terminals()]
    print(f"\n  Tree tips: {len(tip_names)}")

    # Detect whether the tree has branch lengths
    has_bl = any(
        (c.branch_length not in (None, 0))
        for c in tree.find_clades()
    )
    print(f"  Tree has branch lengths: {has_bl}")
    if not has_bl:
        print("  -> using TOPOLOGICAL distances (edge counts) for test (C).")

    # Precompute topological distance matrix once (used by test C)
    print("  Building topological distance matrix...")
    dmat = build_topological_distance_matrix(tree, tip_names)

    # Load MICs & resistant labels
    mics = pd.read_csv(MICS, sep="\t", index_col=0)
    labels = load_resistant_labels(mics, tip_names)

    rows = []
    detail_rows = []

    for ab in ANTIBIOTICS:
        resistant = labels[ab]
        n_res = len(resistant)
        print("\n" + "─" * 68)
        print(f"  {ab}: {n_res} resistant strains")
        print(f"    {sorted(resistant)}")

        for s in sorted(resistant):
            detail_rows.append({"antibiotic": ab, "strain": s, "status": "resistant"})

        # (A) Monophyletic clusters
        res_A = count_monophyletic_clusters(tree, resistant)
        n_clusters = res_A[0]
        clusters = res_A[1] if len(res_A) > 1 else []
        singletons = res_A[2] if len(res_A) > 2 else set()
        print(f"\n  (A) Independent monophyletic resistant clusters: {n_clusters}")
        for i, c in enumerate(clusters, 1):
            print(f"        cluster {i} ({len(c)} strains): {sorted(c)}")
        if singletons:
            print(f"        singletons ({len(singletons)}): {sorted(singletons)}")

        # (B) Fitch parsimony changes
        n_changes = fitch_changes(tree, resistant)
        print(f"\n  (B) Minimum state changes (Fitch parsimony): {n_changes}")

        # (C) Topological clustering permutation test
        obs, null_mean, p = topological_clustering_test(
            dmat, resistant, tip_names, n_perm=N_PERM, seed=SEED)
        print(f"\n  (C) Topological clustering test ({N_PERM} permutations):")
        print(f"        observed mean pairwise topo-distance : {obs:.3f}")
        print(f"        null mean (random labels)            : {null_mean:.3f}")
        print(f"        p-value (clustered)                  : {p:.4f}")
        if p < 0.05:
            print(f"        -> resistant strains ARE significantly clustered")
        else:
            print(f"        -> resistant strains NOT more clustered than random")

        # (D) Parsimony permutation test
        obs_ch, null_ch, frac_le = parsimony_permutation_test(
            tree, resistant, tip_names, n_perm=N_PERM, seed=SEED)
        print(f"\n  (D) Parsimony permutation test ({N_PERM} permutations):")
        print(f"        observed changes            : {obs_ch}")
        print(f"        null mean changes (random)  : {null_ch:.2f}")
        print(f"        percentile of observed      : {frac_le*100:.1f}%")

        # Verdict per antibiotic
        independent = max(n_clusters, n_changes)
        if independent <= 2:
            verdict = "NO recoverable signal (confounded with lineage)"
        elif independent == 3:
            verdict = "WEAK / borderline signal"
        else:
            verdict = "Signal potentially recoverable (>=4 independent events)"
        print(f"\n  ==> VERDICT [{ab}]: {verdict}")

        rows.append({
            "antibiotic": ab,
            "n_resistant": n_res,
            "monophyletic_clusters": n_clusters,
            "fitch_changes": n_changes,
            "independent_events": independent,
            "topo_obs": round(obs, 3) if not np.isnan(obs) else np.nan,
            "topo_null": round(null_mean, 3) if not np.isnan(null_mean) else np.nan,
            "topo_p_clustered": round(p, 4) if not np.isnan(p) else np.nan,
            "parsimony_null_mean": round(null_ch, 2) if not np.isnan(null_ch) else np.nan,
            "parsimony_percentile": round(frac_le*100, 1) if not np.isnan(frac_le) else np.nan,
            "verdict": verdict,
        })

    # Save
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "signal_check_summary.tsv", sep="\t", index=False)
    pd.DataFrame(detail_rows).to_csv(
        OUT / "resistant_strains_per_antibiotic.tsv", sep="\t", index=False)

    print("\n" + "=" * 68)
    print("  SUMMARY")
    print("=" * 68)
    print(summary.to_string(index=False))
    print(f"\n  Saved: {OUT/'signal_check_summary.tsv'}")

    print("""
  INTERPRETATION GUIDE (Fran's criterion):
    independent_events <= 2  ->  no correction can recover signal;
                                 paper becomes methodological.
    independent_events >= 4  ->  signal is testable; proceed with the full
                                 Prokka/Panaroo/Snippy/pyseer pipeline.
    == 3                     ->  borderline; proceed but interpret with care.

  Test (C) topological clustering: p >= 0.05 means resistant strains are NOT
  significantly more clustered than random tips -> supports recoverable signal.
  A very low p means strong lineage confounding.

  Test (D) parsimony percentile: the percentile of the observed number of
  state changes within the null distribution. A high percentile means
  resistance arose independently more times than a random labelling would
  -> strong, non-lineage signal.

  NOTE: this tree has no branch lengths, so (C) uses topological (edge-count)
  distances. When the new IQ-TREE tree (with branch lengths) from the Panaroo
  core alignment is available, re-run this script on it to get patristic
  distances as an additional robustness check.
""")


if __name__ == "__main__":
    main()

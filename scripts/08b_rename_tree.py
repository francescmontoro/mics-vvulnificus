#!/usr/bin/env python3
"""
08b_rename_tree.py
------------------
Renames tips in the Newick tree from SNP pipeline names to MIC names,
so that Scoary can match the tree tips to the traits/genes CSV files.

Usage:
    python3 08b_rename_tree.py

Inputs:
    ../data/raw/phylo_tree.nwk

Outputs:
    ../data/raw/phylo_tree_renamed.nwk

Author: Francesc Montoro (2025)
"""

from pathlib import Path

BASE = Path(__file__).parent.parent
RAW  = BASE / "data" / "raw"

# Mapping: original tree tip name → MIC name
RENAME = {
    "11028":        "BT3_11028",
    "12":           "BT3_12",
    "162":          "BT3_162",
    "32":           "BT3_32",
    "94385":        "94-3-85",
    "960426-1_4C":  "960426_1/4C",
    "CIP8190":      "C1P8190",
    "V4":           "V_4",
    "CECT4604":     "CECT_4604",
    "CECT4606":     "CECT_4606",
    "CECT4608":     "CECT_4608",
    "CECT4865":     "CECT_4865",
    "CECT4866":     "CECT_4866",
    "CECT4999":     "CECT_4999",
    "CECT5198":     "CECT_5198",
    "CECT5763":     "CECT_5763",
    "CECT5769":     "CECT_5769",
    "CECT7030":     "CECT_7030",
    "CECT898":      "CECT_898",
    "CG100":        "CG_100",
    # These already match — listed for completeness
    "94-9-112":     "94-9-112",
    "94-9-119":     "94-9-119",
    "95-8-161":     "95-8-161",
    "95-8-6":       "95-8-6",
    "95-8-7":       "95-8-7",
    "A14":          "A14",
    "PD-2-51":      "PD-2-51",
    "Ra_3":         "Ra_3",
    "RIU1":         "RIU1",
    "YJ016":        "YJ016",
}

if __name__ == "__main__":
    in_file  = RAW / "phylo_tree.nwk"
    out_file = RAW / "phylo_tree_renamed.nwk"

    tree_str = in_file.read_text().strip()
    print(f"Original tree: {in_file}")
    print(f"Original: {tree_str[:80]}...")

    renamed = tree_str
    # Replace longest names first to avoid partial matches
    for old, new in sorted(RENAME.items(), key=lambda x: -len(x[0])):
        if old != new:
            renamed = renamed.replace(old, new)

    out_file.write_text(renamed + "\n")
    print(f"\nRenamed tree: {out_file}")
    print(f"Renamed: {renamed[:80]}...")

    # Verify all MIC names are now in the tree
    from pathlib import Path
    mic_names = list(RENAME.values())
    missing = [n for n in mic_names if n not in renamed]
    if missing:
        print(f"\n[WARNING] Still missing in tree: {missing}")
    else:
        print(f"\nAll {len(mic_names)} strain names found in renamed tree. OK")

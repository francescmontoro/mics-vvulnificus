#!/usr/bin/env python3
"""
01_convert_mics.py
------------------
Converts raw MIC values (mg/L) to a normalised semi-quantitative scale [0-5].
Each antibiotic has its own concentration ladder; values are mapped linearly
so that 0 = most sensitive and 5 = most resistant.

Handles censored values:
  - "<X"  → mapped to the lowest scale step (0)
  - ">X"  → mapped to the highest scale step (5)
  - Exact values → linearly interpolated within the ladder

Usage:
    python3 01_convert_mics.py

Inputs:
    ../data/raw/MICs_all_Vv.xlsx

Outputs:
    ../data/processed/MICs_converted.tsv   (semi-quantitative scale)
    ../data/processed/MICs_summary.tsv     (per-strain resistance score)

Author: Francesc Montoro (refactored 2025)
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
RAW  = BASE / "data" / "raw"
PROC = BASE / "data" / "processed"
PROC.mkdir(parents=True, exist_ok=True)

# ── Antibiotic concentration ladders ──────────────────────────────────────────
# Format: antibiotic_code → ordered list of concentration values (mg/L)
# The script maps these linearly to [0, 5].
# "<first" → 0,  ">last" → 5,  exact values → interpolated.

LADDERS = {
    "AMI":  [8, 16, 32, 64],                       # amikacin
    "SPE":  [8, 16, 32, 64],                       # spectinomycin
    "AMP":  [4, 8, 16, 32],                        # ampicillin
    "FAZ":  [4, 8, 16, 32],                        # cefazolin
    "FEP":  [4, 8, 16, 32],                        # cefepime
    "FUR":  [4, 8, 16, 32],                        # nitrofurantoin (old code)
    "FOX":  [4, 8, 16, 32],                        # cefoxitin
    "CIP":  [0.5, 1, 2, 4],                        # ciprofloxacin
    "NOV":  [0.5, 1, 2, 4],                        # novobiocin
    "CLI":  [0.5, 1, 2, 4],                        # clindamycin
    "GAT":  [1, 2, 4, 8],                          # gatifloxacin
    "MERO": [1, 2, 4, 8],                          # meropenem
    "FFN":  [1, 2, 4, 8],                          # florfenicol
    "GEN":  [2, 4, 8, 16],                         # gentamicin
    "IMI":  [2, 4, 8, 16],                         # imipenem
    "POD":  [2, 4, 8, 16],                         # cefpodoxime
    "NIT":  [16, 32, 64, 128],                     # nitrofurantoin
    "PIP":  [16, 32, 64, 128],                     # piperacillin
    "SDM":  [32, 64, 128, 256],                    # sulphadimethoxine
    "STZ":  [32, 64, 128, 256],                    # sulphathiazole
    "TOB":  [4, 8],                                # tobramycin
    "ENRO": [0.12, 0.25, 0.5, 1, 2],              # enrofloxacin
    "GEN2": [0.5, 1, 2, 4, 8],                    # gentamicin (low range)
    "TIO":  [0.25, 0.5, 1, 2, 4],                 # ceftiofur
    "NEO":  [2, 4, 8, 16, 32],                    # neomycin
    "ERY":  [0.12, 0.25, 0.5, 1, 2, 4],           # erythromycin
    "OXY":  [0.25, 0.5, 1, 2, 4, 8],             # oxytetracycline
    "TET":  [0.25, 0.5, 1, 2, 4, 8],             # tetracycline
    "AMOX": [0.25, 0.5, 1, 2, 4, 8, 16],         # amoxicillin
    "PEN":  [0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8], # penicillin
    "STR":  [8, 16, 32, 64, 128, 256, 512, 1024], # streptomycin
    "AXO":  [1, 2, 4, 8, 16, 32, 64],            # ceftriaxone
    "TAZ":  [1, 2, 4, 8, 16, 32],                # ceftazidime
    "AZT":  [8, 16, 32],                          # aztreonam
    "TANS": [8, 16, 32],                          # trimethoprim
    "TYLT": [2.5, 5, 10, 20],                     # tylosin tartrate
    # Combination antibiotics — primary drug concentration
    "A/S2": [4, 8, 16, 32],                       # amoxicillin/sulbactam
    "P/T4": [16, 32, 64, 128],                    # piperacillin/tazobactam
    "TIM2": [16, 32, 64],                         # ticarcillin/clavulanate
    "SXT":  [0.6, 1, 2, 4],                      # trimethoprim/sulfamethoxazole
    "SXT2": [0.5, 1, 2],                         # trimethoprim/sulfamethoxazole v2
}


def parse_mic_value(raw_val):
    """
    Parse a raw MIC string into (kind, float_value).
    kind: 'low' (<X), 'high' (>X), 'exact' (plain number)
    For combination antibiotics (e.g. "16/4"), takes the first value.
    Returns None if unparseable.
    """
    if pd.isna(raw_val):
        return None
    s = str(raw_val).strip()
    censored_low  = s.startswith("<")
    censored_high = s.startswith(">")
    s_clean = s.lstrip("<>").split("/")[0]  # take primary drug for combos
    try:
        val = float(s_clean)
        if censored_low:
            return ("low", val)
        if censored_high:
            return ("high", val)
        return ("exact", val)
    except ValueError:
        return None


def mic_to_scale(raw_val, ladder):
    """
    Map a raw MIC value to [0, 5] using the antibiotic's concentration ladder.

    Scale logic:
      - '<lowest_rung'  → 0.0  (fully sensitive)
      - '>highest_rung' → 5.0  (fully resistant)
      - Exact values    → linearly interpolated between rungs
    """
    parsed = parse_mic_value(raw_val)
    if parsed is None:
        return np.nan

    kind, val = parsed

    if kind == "low":
        return 0.0
    if kind == "high":
        return 5.0

    # Exact value
    n = len(ladder)
    if val <= ladder[0]:
        return 0.0
    if val >= ladder[-1]:
        return 5.0

    step = 5.0 / n
    for i in range(n - 1):
        if ladder[i] <= val <= ladder[i + 1]:
            frac = (val - ladder[i]) / (ladder[i + 1] - ladder[i])
            return round((i + frac) * step, 4)
    return np.nan


def convert_mics(input_path: Path) -> pd.DataFrame:
    """Load raw MICs Excel and convert all antibiotic columns to [0-5] scale."""
    df = pd.read_excel(input_path, sheet_name="MICs")

    # Clean up
    # Fix duplicate column names created by Excel when two columns share the same name
    df = df.rename(columns={"GEN.1": "GEN2", "SXT.1": "SXT2"})

    if "PDF" in df.columns:
        df = df.drop(columns=["PDF"])
    df = df.dropna(subset=["STRAIN"])
    df["STRAIN"] = df["STRAIN"].astype(str).str.strip()
    df = df.set_index("STRAIN")

    converted = pd.DataFrame(index=df.index)
    unknown_abs = []

    for col in df.columns:
        ab = col.strip()
        if ab in LADDERS:
            converted[ab] = df[col].apply(
                lambda x: mic_to_scale(x, LADDERS[ab])
            )
        else:
            unknown_abs.append(ab)
            converted[ab] = np.nan

    if unknown_abs:
        print(f"\n  [WARNING] No ladder defined for: {unknown_abs}")
        print(f"  These columns are NaN — add them to LADDERS dict if needed.")

    return converted


def compute_resistance_score(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-strain resistance summary:
      - resistance_score : mean of all [0-5] values (higher = more resistant)
      - n_antibiotics    : number of antibiotics tested (non-NaN)
      - n_fully_sensitive: antibiotics where score == 0
      - n_fully_resistant: antibiotics where score == 5
      - pct_sensitive    : % of antibiotics fully sensitive
    Sorted from most sensitive to most resistant.
    """
    summary = pd.DataFrame(index=df.index)
    summary["resistance_score"]    = df.mean(axis=1).round(4)
    summary["n_antibiotics"]       = df.notna().sum(axis=1)
    summary["n_fully_sensitive"]   = (df == 0).sum(axis=1)
    summary["n_fully_resistant"]   = (df == 5).sum(axis=1)
    summary["pct_sensitive"]       = (
        summary["n_fully_sensitive"] / summary["n_antibiotics"] * 100
    ).round(2)
    return summary.sort_values("resistance_score")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    input_file = RAW / "MICs_all_Vv.xlsx"

    print("=" * 60)
    print("  MIC Conversion — V. vulnificus")
    print("=" * 60)
    print(f"\nReading: {input_file}")

    converted = convert_mics(input_file)
    summary   = compute_resistance_score(converted)

    # Save
    out_mics = PROC / "MICs_converted.tsv"
    out_summ = PROC / "MICs_summary.tsv"
    converted.to_csv(out_mics, sep="\t")
    summary.to_csv(out_summ,   sep="\t")

    print(f"\nOutputs:")
    print(f"  {out_mics}")
    print(f"  {out_summ}")
    print(f"\nDataset: {converted.shape[0]} strains × {converted.shape[1]} antibiotics")

    print(f"\n── Top 10 most SENSITIVE strains ──")
    print(summary.head(10).to_string())

    print(f"\n── Top 5 most RESISTANT strains ──")
    print(summary.tail(5).to_string())

    # Sanity check
    nan_pct = converted.isna().sum().sum() / converted.size * 100
    print(f"\nMissing values: {nan_pct:.1f}%")
    if nan_pct > 5:
        print("  [WARNING] High NaN rate — check LADDERS config or input data.")
    else:
        print("  [OK] NaN rate acceptable.")
    print("\nDone.")

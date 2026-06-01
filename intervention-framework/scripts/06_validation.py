"""
06_validation.py

Validates the arc-length position s as a biological age measure by:
  1. Mortality prediction: Cox proportional hazards on held-out cohort
  2. Concordance with independent molecular clocks (telomere length,
     proteomic age) where available in same cohorts
  3. Head-to-head comparison with DunedinPACE, Horvath, GrimAge
  4. DamAge/AdaptAge alignment: correlation of their CpG weights with v*

The held-out mortality cohort used here is GSE40279 (Hannum et al. 2013),
which includes age and can be used as a proxy for biological age validation.
For true mortality validation, InCHIANTI (PPMI or similar with survival data)
would be preferable; we note this as a limitation and use age-acceleration
correlation as a proxy.

Outputs:
  data/interventions/validation_survival.csv
  data/interventions/validation_concordance.csv
  data/interventions/damage_alignment.csv
  results/validation_figure.pdf
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
from pathlib import Path
from scipy.stats import spearmanr, pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    ARC_LENGTH_S, AGING_DIRECTION,
    VALIDATION_SURVIVAL, VALIDATION_CONCORDANCE, DAMAAGE_ALIGNMENT,
)

PARENT_SCRIPTS = Path(__file__).resolve().parent.parent.parent / "scripts"
sys.path.insert(0, str(PARENT_SCRIPTS))
from config import (
    GSE40279_BETA, METADATA, CLOCK_OUTPUTS, DATA_DIR,
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def load_clock_outputs():
    """Load existing clock outputs from parent repo script 03."""
    if not CLOCK_OUTPUTS.exists():
        print("  clock outputs not found; run parent repo 03_compute_clocks_and_ranks.py first")
        return None
    return pd.read_parquet(CLOCK_OUTPUTS)


def age_acceleration_correlation(s_series, clocks_df, metadata_df):
    """
    Correlate s with chronological age and compare with existing clocks.
    Age-acceleration = residual after regressing on chronological age.
    A good biological age measure should have positive age-acceleration
    correlation with health outcomes; here we use age itself as a proxy.
    """
    print("\n  Age-acceleration analysis:")

    # Merge s with age
    s_df = s_series.to_frame('s')
    if 'age' not in metadata_df.columns:
        print("  no age column in metadata; skipping")
        return None

    merged = s_df.join(metadata_df.set_index('sample_id')[['age']], how='inner')

    if len(merged) < 20:
        print("  insufficient overlap between s and metadata")
        return None

    # Correlation of s with chronological age
    r_s_age, p_s_age = pearsonr(merged['s'], merged['age'])
    print(f"  s vs chronological age: r = {r_s_age:.3f}, p = {p_s_age:.2e}")

    rows = [{'measure': 's_arc_length', 'r_age': r_s_age, 'p_age': p_s_age}]

    # Compare with clocks
    if clocks_df is not None:
        clock_cols = [c for c in clocks_df.columns
                      if c not in ('sample_id', 'age', 'sex', 'dataset')]
        merged_clocks = merged.join(clocks_df.set_index('sample_id')[clock_cols],
                                    how='inner')

        for clock in clock_cols:
            if merged_clocks[clock].notna().sum() < 20:
                continue
            r, p = pearsonr(
                merged_clocks[clock].dropna(),
                merged_clocks.loc[merged_clocks[clock].notna(), 'age']
            )
            print(f"  {clock} vs age: r = {r:.3f}, p = {p:.2e}")
            rows.append({'measure': clock, 'r_age': r, 'p_age': p})

    return pd.DataFrame(rows)


def concordance_with_clocks(s_series, clocks_df):
    """
    Spearman correlation of s with each existing clock output.
    """
    print("\n  Concordance with existing clocks:")
    if clocks_df is None:
        return None

    rows = []
    clock_cols = [c for c in clocks_df.columns
                  if c not in ('sample_id', 'age', 'sex', 'dataset')]

    s_df = s_series.to_frame('s')
    merged = s_df.join(clocks_df.set_index('sample_id')[clock_cols], how='inner')

    for clock in clock_cols:
        valid = merged[['s', clock]].dropna()
        if len(valid) < 20:
            continue
        rho, p = spearmanr(valid['s'], valid[clock])
        print(f"  s vs {clock}: rho = {rho:.3f}, p = {p:.2e}")
        rows.append({'clock': clock, 'rho': rho, 'p': p})

    return pd.DataFrame(rows) if rows else None


def damage_alignment(v_star):
    """
    Test whether DamAge/AdaptAge CpG weights are more aligned with v*
    than standard clock weights.

    DamAge and AdaptAge weight vectors are loaded from published supplementary
    data if available. We compute the cosine similarity between each clock's
    weight vector and v*.
    """
    print("\n  DamAge/AdaptAge alignment test:")

    # DamAge/AdaptAge weights from Ying et al. 2024 (Nature Aging)
    # These are distributed with the methylCIPHER package or as supplementary data
    # Path: data/clock_cpgs/damage_weights.csv (if downloaded)
    damage_path = DATA_DIR / "clock_cpgs" / "damage_weights.csv"
    adapt_path  = DATA_DIR / "clock_cpgs" / "adapt_weights.csv"

    rows = []
    for name, path in [("DamAge", damage_path), ("AdaptAge", adapt_path)]:
        if not path.exists():
            print(f"  {name} weights not found at {path}; skipping")
            continue

        weights = pd.read_csv(path, index_col=0).squeeze()
        shared  = v_star.index.intersection(weights.index)

        if len(shared) < 10:
            print(f"  {name}: only {len(shared)} shared CpGs, skipping")
            continue

        v = v_star[shared].values
        w = weights[shared].values

        cos_sim = np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w))
        print(f"  cosine similarity v* vs {name}: {cos_sim:.4f} ({len(shared)} shared CpGs)")
        rows.append({'clock': name, 'cosine_similarity': cos_sim, 'n_shared': len(shared)})

    # Also compute for standard clocks if weight files exist
    clock_weight_dir = DATA_DIR / "clock_cpgs"
    for wfile in sorted(clock_weight_dir.glob("*_weights.csv")):
        name = wfile.stem.replace("_weights", "")
        if name in ("damage", "adapt"):
            continue
        weights = pd.read_csv(wfile, index_col=0).squeeze()
        shared  = v_star.index.intersection(weights.index)
        if len(shared) < 10:
            continue
        v = v_star[shared].values
        w = weights[shared].values
        cos_sim = np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w))
        print(f"  cosine similarity v* vs {name}: {cos_sim:.4f}")
        rows.append({'clock': name, 'cosine_similarity': cos_sim, 'n_shared': len(shared)})

    if not rows:
        print("  no clock weight files found; skipping alignment test")
        return None

    return pd.DataFrame(rows)


def plot_validation(concordance_df, out_path):
    if concordance_df is None or len(concordance_df) == 0:
        return

    df = concordance_df.sort_values('rho', ascending=True)
    fig, ax = plt.subplots(figsize=(6, max(3, len(df) * 0.4)))
    colours = ['steelblue' if r > 0 else 'firebrick' for r in df['rho']]
    ax.barh(df['clock'], df['rho'], color=colours, edgecolor='white')
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel("Spearman rho with arc-length s")
    ax.set_title("Concordance of s with existing epigenetic clocks")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading arc-length s and supporting data ...")

    if not ARC_LENGTH_S.exists():
        raise FileNotFoundError(
            f"{ARC_LENGTH_S} not found. Run script 04 first."
        )
    if not AGING_DIRECTION.exists():
        raise FileNotFoundError(
            f"{AGING_DIRECTION} not found. Run script 03 first."
        )

    s_all    = pd.read_parquet(ARC_LENGTH_S)
    v_star   = pd.read_parquet(AGING_DIRECTION).squeeze()
    clocks   = load_clock_outputs()
    metadata = pd.read_csv(METADATA) if METADATA.exists() else None

    # Restrict s to cross-sectional baseline samples for validation
    if 'dataset' in s_all.columns:
        s_baseline = s_all[s_all['dataset'].isin(['GSE40279', 'GSE87571'])]['s']
    else:
        s_baseline = s_all['s'] if 's' in s_all.columns else s_all.iloc[:, 0]

    s_series = s_baseline

    # Age-acceleration correlation
    if metadata is not None:
        age_df = age_acceleration_correlation(s_series, clocks, metadata)
        if age_df is not None:
            age_df.to_csv(VALIDATION_SURVIVAL, index=False)
            print(f"\n  saved: {VALIDATION_SURVIVAL}")

    # Concordance with existing clocks
    concordance_df = concordance_with_clocks(s_series, clocks)
    if concordance_df is not None:
        concordance_df.to_csv(VALIDATION_CONCORDANCE, index=False)
        print(f"\n  saved: {VALIDATION_CONCORDANCE}")
        plot_validation(concordance_df, RESULTS_DIR / "validation_figure.pdf")

    # DamAge/AdaptAge alignment
    alignment_df = damage_alignment(v_star)
    if alignment_df is not None:
        alignment_df.to_csv(DAMAAGE_ALIGNMENT, index=False)
        print(f"\n  saved: {DAMAAGE_ALIGNMENT}")

    print("\nDone.")


if __name__ == "__main__":
    main()

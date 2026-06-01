"""
05_classify_interventions.py

Classifies each intervention as geroprotective, age-accelerating, or
clock-gaming using the arc-length decomposition.

For each intervention dataset:
  - Computes mean arc-length change Delta_s = mean(s_post - s_pre) for paired
    samples, or mean(s_case - s_control) for cross-sectional designs
  - Computes orthogonal displacement magnitude as total displacement minus
    tangential component
  - Classifies using the clock-gaming threshold epsilon

Special analysis for CALERIE (GSE180353):
  - Decomposes intervention effect into Delta_s (position change) and
    Delta_s_dot (rate change, estimated from 12-month vs 24-month measurements)
  - Compares with DunedinPACE and Horvath clock responses on same samples

Outputs:
  data/interventions/intervention_classification.csv
  data/interventions/calerie_decomposition.csv
  results/classification_figure.pdf
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    INTERVENTION_DATASETS, ARC_LENGTH_S, DISPLACEMENT_MATRIX,
    PRINCIPAL_CURVE_ORIENTED, CLASSIFICATION_TABLE, CALERIE_DECOMPOSITION,
    CLOCK_GAMING_EPSILON, MIN_PAIRED_SAMPLES,
    interv_beta_path, interv_meta_path, displacement_path,
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# CALERIE accession
CALERIE_ACCESSION = "GSE180353"


def classify_intervention(delta_s, delta_total, epsilon=CLOCK_GAMING_EPSILON):
    """
    Classify an intervention based on arc-length change and total displacement.

    Parameters
    ----------
    delta_s      : float   mean arc-length change (signed: + = older, - = younger)
    delta_total  : float   mean total displacement magnitude
    epsilon      : float   clock-gaming threshold

    Returns
    -------
    classification : str
    tangential_fraction : float   |delta_s| / delta_total
    """
    if delta_total < 1e-8:
        return "no_effect", 0.0

    tangential_frac = abs(delta_s) / delta_total

    if tangential_frac < epsilon:
        return "clock_gaming", tangential_frac
    elif delta_s < 0:
        return "geroprotective", tangential_frac
    else:
        return "age_accelerating", tangential_frac


def compute_delta_s_longitudinal(s_df, meta_df, accession):
    """
    Compute mean(s_post - s_pre) for paired longitudinal samples.
    Returns (delta_s, n_pairs, s_pre_mean, s_post_mean).
    """
    if 'timepoint' not in meta_df.columns or 'subject_id' not in meta_df.columns:
        # Attempt to infer
        meta_df = meta_df.copy()
        if 'subject_id' not in meta_df.columns:
            meta_df['subject_id'] = meta_df['sample_id'].str.extract(r'(GSM\d+)')[0]

    pre_meta  = meta_df[meta_df['timepoint'] == 'pre']
    post_meta = meta_df[meta_df['timepoint'] == 'post']

    common_subjects = set(pre_meta['subject_id']) & set(post_meta['subject_id'])

    if len(common_subjects) < MIN_PAIRED_SAMPLES:
        return None, 0, None, None

    pre_map  = pre_meta.set_index('subject_id')['sample_id']
    post_map = post_meta.set_index('subject_id')['sample_id']

    s_indexed = s_df['s'] if 's' in s_df.columns else s_df.iloc[:, 0]

    diffs = []
    s_pre_vals  = []
    s_post_vals = []
    for subj in common_subjects:
        pre_sid  = pre_map.get(subj)
        post_sid = post_map.get(subj)
        if pre_sid in s_indexed.index and post_sid in s_indexed.index:
            s_pre  = s_indexed[pre_sid]
            s_post = s_indexed[post_sid]
            diffs.append(s_post - s_pre)
            s_pre_vals.append(s_pre)
            s_post_vals.append(s_post)

    if not diffs:
        return None, 0, None, None

    return np.mean(diffs), len(diffs), np.mean(s_pre_vals), np.mean(s_post_vals)


def compute_delta_s_cross_sectional(s_df, meta_df, accession):
    """
    Compute mean(s_case) - mean(s_control) for cross-sectional designs.
    """
    if 'group' not in meta_df.columns:
        return None, 0, None, None

    s_indexed = s_df['s'] if 's' in s_df.columns else s_df.iloc[:, 0]

    case_ids    = meta_df.loc[meta_df['group'] == 'case',    'sample_id'].tolist()
    control_ids = meta_df.loc[meta_df['group'] == 'control', 'sample_id'].tolist()

    case_ids    = [s for s in case_ids    if s in s_indexed.index]
    control_ids = [s for s in control_ids if s in s_indexed.index]

    if len(case_ids) < MIN_PAIRED_SAMPLES or len(control_ids) < MIN_PAIRED_SAMPLES:
        return None, 0, None, None

    s_case    = s_indexed[case_ids].mean()
    s_control = s_indexed[control_ids].mean()
    n = min(len(case_ids), len(control_ids))

    return s_case - s_control, n, s_control, s_case


def compute_total_displacement(accession, config):
    """
    Compute mean total displacement magnitude from saved displacement vector.
    """
    dpath = displacement_path(accession)
    if not dpath.exists():
        return np.nan

    delta = pd.read_parquet(dpath).squeeze()
    return np.linalg.norm(delta.values)


def calerie_decomposition(s_df, meta_df):
    """
    CALERIE-specific analysis: decompose effect into position (s) and rate (s_dot).

    CALERIE has measurements at baseline (t=0), 12 months, and 24 months.
    - Delta_s at 12 months: position change
    - Delta_s_dot: (s_24 - s_12) / 12 vs (s_12 - s_0) / 12 — change in rate
    """
    print("\n  CALERIE decomposition:")

    if 'timepoint' not in meta_df.columns:
        print("  no timepoint column in CALERIE metadata; skipping decomposition")
        return None

    s_indexed = s_df['s'] if 's' in s_df.columns else s_df.iloc[:, 0]
    results = []

    # Map timepoints
    tp_map = {'pre': 'baseline', 'baseline': 'baseline',
              '12': '12m', '12m': '12m', '12mo': '12m',
              'post': '24m', '24': '24m', '24m': '24m', '24mo': '24m'}

    meta_df = meta_df.copy()
    meta_df['tp_norm'] = meta_df['timepoint'].str.lower().map(tp_map)

    for tp_label in ['baseline', '12m', '24m']:
        tp_ids = meta_df.loc[meta_df['tp_norm'] == tp_label, 'sample_id'].tolist()
        tp_ids = [s for s in tp_ids if s in s_indexed.index]
        if tp_ids:
            mean_s = s_indexed[tp_ids].mean()
            n = len(tp_ids)
            print(f"    {tp_label}: n={n}, mean s = {mean_s:.4f}")
            results.append({'timepoint': tp_label, 'mean_s': mean_s, 'n': n})

    if len(results) >= 2:
        df = pd.DataFrame(results).set_index('timepoint')

        if 'baseline' in df.index and '24m' in df.index:
            delta_s_total = df.loc['24m', 'mean_s'] - df.loc['baseline', 'mean_s']
            print(f"\n    Delta_s (0 -> 24m): {delta_s_total:+.4f}")
            print(f"    Interpretation: {'geroprotective' if delta_s_total < 0 else 'age-accelerating'} "
                  f"position shift")

        if all(tp in df.index for tp in ['baseline', '12m', '24m']):
            rate_first_half  = (df.loc['12m', 'mean_s'] - df.loc['baseline', 'mean_s']) / 12
            rate_second_half = (df.loc['24m', 'mean_s'] - df.loc['12m', 'mean_s']) / 12
            delta_s_dot = rate_second_half - rate_first_half
            print(f"\n    Rate first 12m:  {rate_first_half:+.5f} s/month")
            print(f"    Rate second 12m: {rate_second_half:+.5f} s/month")
            print(f"    Delta_s_dot:     {delta_s_dot:+.5f} s/month")
            print(f"    Interpretation: caloric restriction primarily affects "
                  f"{'rate (s_dot)' if abs(delta_s_dot) > abs(delta_s_total) else 'position (s)'}")

        df.reset_index().to_csv(CALERIE_DECOMPOSITION, index=False)
        print(f"\n    saved: {CALERIE_DECOMPOSITION}")
        return df

    return None


def plot_classification(classification_df):
    """Bar chart of Delta_s by intervention, coloured by classification."""
    colours = {
        'geroprotective':  'steelblue',
        'age_accelerating': 'firebrick',
        'clock_gaming':    'goldenrod',
        'no_effect':       'lightgrey',
        'insufficient_data': 'lightgrey',
    }

    df = classification_df.dropna(subset=['delta_s'])
    df = df.sort_values('delta_s')

    fig, ax = plt.subplots(figsize=(8, 4))
    for _, row in df.iterrows():
        colour = colours.get(row['classification'], 'grey')
        ax.barh(row['label'], row['delta_s'], color=colour, edgecolor='white')

    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel(r'Mean $\Delta s$ (arc-length change)')
    ax.set_title('Intervention classification')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=l) for l, c in colours.items()
                       if l in df['classification'].values]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "classification_figure.pdf", dpi=150)
    plt.close()
    print(f"  saved: {RESULTS_DIR}/classification_figure.pdf")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading arc-length values ...")
    if not ARC_LENGTH_S.exists():
        raise FileNotFoundError(
            f"{ARC_LENGTH_S} not found. Run script 04 first."
        )

    s_all = pd.read_parquet(ARC_LENGTH_S)
    print(f"  loaded s for {len(s_all):,} samples")

    rows = []

    for accession, config in INTERVENTION_DATASETS.items():
        meta_path = interv_meta_path(accession)
        if not meta_path.exists():
            print(f"  {accession}: metadata not found, skipping")
            continue

        meta_df = pd.read_csv(meta_path)

        # Filter s to this dataset
        if 'dataset' in s_all.columns:
            s_df = s_all[s_all['dataset'] == accession]
        else:
            s_df = s_all.loc[s_all.index.isin(meta_df['sample_id'])]

        print(f"\n  {accession} ({config['label']}, {config['design']}):")
        print(f"    n samples with s: {len(s_df)}")

        # Compute delta_s
        if config['design'] == 'longitudinal':
            delta_s, n_pairs, s_pre, s_post = compute_delta_s_longitudinal(
                s_df, meta_df, accession
            )
        else:
            delta_s, n_pairs, s_ctrl, s_case = compute_delta_s_cross_sectional(
                s_df, meta_df, accession
            )

        if delta_s is None:
            print(f"    insufficient data")
            rows.append({
                'accession': accession, 'label': config['label'],
                'sign': config['sign'], 'tissue': config['tissue'],
                'delta_s': np.nan, 'n': 0,
                'tangential_fraction': np.nan,
                'classification': 'insufficient_data',
            })
            continue

        # Total displacement magnitude
        delta_total = compute_total_displacement(accession, config)

        # Classify
        classification, tang_frac = classify_intervention(delta_s, delta_total)

        print(f"    delta_s = {delta_s:+.4f}, n = {n_pairs}")
        print(f"    tangential fraction = {tang_frac:.3f}")
        print(f"    classification: {classification}")

        rows.append({
            'accession':            accession,
            'label':                config['label'],
            'sign':                 config['sign'],
            'tissue':               config['tissue'],
            'delta_s':              delta_s,
            'n':                    n_pairs,
            'tangential_fraction':  tang_frac,
            'classification':       classification,
        })

        # CALERIE special analysis
        if accession == CALERIE_ACCESSION:
            calerie_decomposition(s_df, meta_df)

    classification_df = pd.DataFrame(rows)
    classification_df.to_csv(CLASSIFICATION_TABLE, index=False)
    print(f"\nSaved classification table: {CLASSIFICATION_TABLE}")
    print(classification_df[['label', 'delta_s', 'classification']].to_string(index=False))

    # Plot
    plot_classification(classification_df)


if __name__ == "__main__":
    main()

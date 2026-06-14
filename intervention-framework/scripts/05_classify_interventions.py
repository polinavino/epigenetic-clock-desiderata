"""
05_classify_interventions.py

Classifies each intervention by projecting its displacement vector onto the
aging direction v*, rather than differencing absolute arc-length positions.

This avoids the cross-dataset normalisation problem: within-subject (post-pre)
or within-dataset (case-control) displacements cancel batch effects, whereas
absolute projection onto the reference curve does not.

Classification (from paper Section 3.5):
  Let d = displacement vector (signed: post-pre or case-control, NOT sign-flipped)
  Let v* = aging direction (unit vector)
  tangential component  = <d, v*>           (signed scalar: + = toward old)
  total magnitude       = ||d||
  tangential fraction   = |<d, v*>| / ||d||  (in [0,1])

  geroprotective   : <d, v*> < 0  and tangential fraction >= epsilon
  age_accelerating : <d, v*> > 0  and tangential fraction >= epsilon
  clock_gaming     : tangential fraction < epsilon

Outputs:
  data/interventions/intervention_classification.csv
  results/classification_figure.pdf
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    INTERVENTION_DATASETS, AGING_DIRECTION, DISPLACEMENT_MATRIX,
    DISPLACEMENT_METADATA, CLASSIFICATION_TABLE,
    CLOCK_GAMING_EPSILON, displacement_path,
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def classify(d_proj, d_norm, epsilon=CLOCK_GAMING_EPSILON):
    """
    d_proj : signed projection <d, v*>  (+ = toward old end)
    d_norm : ||d||
    """
    if d_norm < 1e-12:
        return "no_effect", 0.0
    tang_frac = abs(d_proj) / d_norm
    if tang_frac < epsilon:
        return "clock_gaming", tang_frac
    elif d_proj < 0:
        return "geroprotective", tang_frac
    else:
        return "age_accelerating", tang_frac


def main():
    print("Loading aging direction v* ...")
    v_star = pd.read_parquet(AGING_DIRECTION).squeeze()
    v_cpgs = v_star.index.tolist()
    v_unit = v_star.values / np.linalg.norm(v_star.values)
    print(f"  v* over {len(v_cpgs):,} CpGs")

    print("\nLoading displacement metadata ...")
    meta = pd.read_csv(DISPLACEMENT_METADATA)

    rows = []
    for _, m in meta.iterrows():
        accession = m['accession']
        label     = m['label']
        sign      = m['sign']

        dpath = displacement_path(accession)
        if not dpath.exists():
            print(f"  {accession}: displacement not found, skipping")
            continue

        # Load RAW (unsigned) displacement: post-pre or case-control
        delta = pd.read_parquet(dpath).squeeze()

        # Restrict to v* CpGs
        shared = [c for c in v_cpgs if c in delta.index]
        if len(shared) < 100:
            print(f"  {accession}: only {len(shared)} shared CpGs, skipping")
            continue

        d = delta[shared].values
        v = v_star[shared].values
        v = v / np.linalg.norm(v)

        d_proj = float(np.dot(d, v))      # signed: + toward old
        d_norm = float(np.linalg.norm(d))

        classification, tang_frac = classify(d_proj, d_norm)

        print(f"\n  {accession} ({label}, expected {sign}):")
        print(f"    <d, v*>          = {d_proj:+.4f}")
        print(f"    ||d||            = {d_norm:.4f}")
        print(f"    tangential frac  = {tang_frac:.3f}")
        print(f"    classification   = {classification}")

        # Sanity: does classification match expected sign?
        expected = 'age_accelerating' if sign == 'plus' else 'geroprotective'
        match = "OK" if classification == expected else "MISMATCH"
        print(f"    expected {expected}: {match}")

        rows.append({
            'accession':           accession,
            'label':               label,
            'expected_sign':       sign,
            'projection':          d_proj,
            'magnitude':           d_norm,
            'tangential_fraction': tang_frac,
            'classification':      classification,
            'matches_expected':    classification == expected,
        })

    df = pd.DataFrame(rows)
    df.to_csv(CLASSIFICATION_TABLE, index=False)
    print(f"\nSaved: {CLASSIFICATION_TABLE}")
    print(df[['label', 'projection', 'tangential_fraction',
              'classification', 'matches_expected']].to_string(index=False))

    # Plot: signed projection, coloured by classification
    colours = {
        'geroprotective':   'steelblue',
        'age_accelerating': 'firebrick',
        'clock_gaming':     'goldenrod',
        'no_effect':        'lightgrey',
    }
    dfp = df.sort_values('projection')
    fig, ax = plt.subplots(figsize=(8, max(3, len(dfp) * 0.6)))
    for _, r in dfp.iterrows():
        ax.barh(r['label'], r['projection'],
                color=colours.get(r['classification'], 'grey'),
                edgecolor='white')
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel(r'Projection onto aging direction  $\langle d, v^* \rangle$'
                  '\n(negative = toward younger, positive = toward older)')
    ax.set_title('Intervention classification (displacement projection)')
    from matplotlib.patches import Patch
    legend = [Patch(facecolor=c, label=l) for l, c in colours.items()
              if l in dfp['classification'].values]
    ax.legend(handles=legend, loc='best', fontsize=8)
    plt.tight_layout()
    out = RESULTS_DIR / "classification_figure.pdf"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  saved: {out}")


if __name__ == "__main__":
    main()

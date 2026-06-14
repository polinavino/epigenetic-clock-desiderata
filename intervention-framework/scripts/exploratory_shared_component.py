"""
exploratory_shared_component.py

Decomposes the smoking and weight-loss displacements into:
  - a SHARED component (the direction they agree on)
  - intervention-SPECIFIC components (orthogonal residuals)

Then tests whether the shared component aligns with the chronological-age
principal curve. This addresses the central question: is the orthogonality
between interventions "noise around a shared aging signal", or is there no
shared aging signal at all?

Key outputs:
  - cosine of shared component with the age trajectory direction
  - enrichment of known aging CpGs in shared vs specific components
  - which CpGs dominate each component
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import INTERV_DIR, displacement_path

PARENT_SCRIPTS = Path(__file__).resolve().parent.parent.parent / "scripts"
sys.path.insert(0, str(PARENT_SCRIPTS))
from config import PRINCIPAL_CURVE, TOP_CPGS, WEIGHTS

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def load_displacements():
    ds = pd.read_parquet(displacement_path("GSE50660")).squeeze()
    dw = pd.read_parquet(displacement_path("GSE272137")).squeeze()
    shared = ds.index.intersection(dw.index)
    return ds[shared], dw[shared]


def decompose(ds, dw):
    """
    Decompose into shared (sum) and antisymmetric (difference) components.
    For an accelerator and a geroprotector that act on a common aging axis,
    we expect them to OPPOSE on the aging axis. So the aging-relevant shared
    structure is in the DIFFERENCE (smoking - weightloss), and the
    intervention-common artefacts are in the SUM.

    But we don't know orientation a priori, so we compute both and test
    each against the age trajectory.
    """
    s = ds.values
    w = dw.values

    # Normalise each to unit length first (direction only)
    s_unit = s / np.linalg.norm(s)
    w_unit = w / np.linalg.norm(w)

    # Symmetric (sum) and antisymmetric (difference) directions
    sum_dir  = s_unit + w_unit
    diff_dir = s_unit - w_unit

    sum_dir  = sum_dir  / (np.linalg.norm(sum_dir)  + 1e-12)
    diff_dir = diff_dir / (np.linalg.norm(diff_dir) + 1e-12)

    return pd.Series(s_unit, index=ds.index), \
           pd.Series(w_unit, index=ds.index), \
           pd.Series(sum_dir, index=ds.index), \
           pd.Series(diff_dir, index=ds.index)


def age_trajectory_direction(cpg_index):
    """
    Compute the local aging direction from the principal curve, restricted
    to the given CpG index. Uses the tangent at the curve midpoint.
    """
    if not PRINCIPAL_CURVE.exists():
        print("  principal curve not found")
        return None

    curve = pd.read_parquet(PRINCIPAL_CURVE)
    shared = [c for c in cpg_index if c in curve.columns]
    if len(shared) < 10:
        print(f"  only {len(shared)} CpGs shared with curve")
        return None

    mid = len(curve) // 2
    # Tangent: average direction over several steps for stability
    span = max(1, len(curve) // 10)
    tangent = curve.iloc[mid + span][shared].values - \
              curve.iloc[mid - span][shared].values
    tangent = tangent / np.linalg.norm(tangent)

    return pd.Series(tangent, index=shared)


def cos_between(a, b):
    shared = a.index.intersection(b.index)
    av = a[shared].values
    bv = b[shared].values
    return np.dot(av, bv) / (np.linalg.norm(av) * np.linalg.norm(bv) + 1e-12)


def load_aging_cpgs():
    if TOP_CPGS.exists():
        with open(TOP_CPGS) as f:
            return set(l.strip() for l in f if l.strip())
    return set()


def main():
    print("Loading displacements ...")
    ds, dw = load_displacements()
    print(f"  shared CpGs: {len(ds):,}")

    cos_sw = cos_between(ds, dw)
    print(f"\n  cosine(smoking, weightloss) = {cos_sw:+.4f} "
          f"({np.degrees(np.arccos(np.clip(cos_sw,-1,1))):.1f} deg)")

    print("\nDecomposing into sum and difference directions ...")
    s_unit, w_unit, sum_dir, diff_dir = decompose(ds, dw)

    print("\nComputing age trajectory direction ...")
    age_dir = age_trajectory_direction(ds.index)
    if age_dir is None:
        print("  cannot compute age direction; aborting")
        return

    print(f"  age direction over {len(age_dir):,} CpGs")

    # Test each component against the age trajectory
    print("\n--- Alignment with chronological-age trajectory ---")
    cos_smoking_age = cos_between(s_unit, age_dir)
    cos_weight_age  = cos_between(w_unit, age_dir)
    cos_sum_age     = cos_between(sum_dir, age_dir)
    cos_diff_age    = cos_between(diff_dir, age_dir)

    print(f"  smoking      vs age: {cos_smoking_age:+.4f} "
          f"({np.degrees(np.arccos(np.clip(abs(cos_smoking_age),0,1))):.1f} deg)")
    print(f"  weight loss  vs age: {cos_weight_age:+.4f} "
          f"({np.degrees(np.arccos(np.clip(abs(cos_weight_age),0,1))):.1f} deg)")
    print(f"  SUM (shared) vs age: {cos_sum_age:+.4f} "
          f"({np.degrees(np.arccos(np.clip(abs(cos_sum_age),0,1))):.1f} deg)")
    print(f"  DIFF         vs age: {cos_diff_age:+.4f} "
          f"({np.degrees(np.arccos(np.clip(abs(cos_diff_age),0,1))):.1f} deg)")

    print("\n  Interpretation:")
    best = max([('smoking', abs(cos_smoking_age)),
                ('weight loss', abs(cos_weight_age)),
                ('shared sum', abs(cos_sum_age)),
                ('difference', abs(cos_diff_age))], key=lambda x: x[1])
    print(f"  Component most aligned with age trajectory: {best[0]} "
          f"(|cos| = {best[1]:.3f})")

    # Enrichment of aging CpGs in each component
    aging_cpgs = load_aging_cpgs()
    if aging_cpgs:
        print(f"\n--- Aging CpG enrichment ({len(aging_cpgs)} known age CpGs) ---")
        for name, comp in [('smoking', s_unit), ('weight loss', w_unit),
                            ('sum', sum_dir), ('difference', diff_dir)]:
            comp_aging = comp[[c for c in comp.index if c in aging_cpgs]]
            comp_other = comp[[c for c in comp.index if c not in aging_cpgs]]
            if len(comp_aging) > 0:
                mean_aging = comp_aging.abs().mean()
                mean_other = comp_other.abs().mean()
                ratio = mean_aging / (mean_other + 1e-12)
                print(f"  {name:12s}: |weight| on aging CpGs / other = {ratio:.2f}")

    # Plot alignment summary
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = ['smoking', 'weight loss', 'sum\n(shared)', 'difference']
    vals   = [abs(cos_smoking_age), abs(cos_weight_age),
              abs(cos_sum_age), abs(cos_diff_age)]
    colours = ['firebrick', 'steelblue', 'purple', 'darkorange']
    ax.bar(labels, vals, color=colours, edgecolor='white')
    ax.set_ylabel('|cosine| with age trajectory')
    ax.set_title('Which component aligns with the chronological-age trajectory?')
    ax.axhline(0.1, color='grey', linestyle='--', alpha=0.5,
               label='~orthogonal threshold')
    ax.legend()
    plt.tight_layout()
    out = RESULTS_DIR / "shared_component_age_alignment.pdf"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"\n  saved: {out}")

    # Save summary
    summary = pd.DataFrame([
        {'component': 'smoking',     'cos_with_age': cos_smoking_age},
        {'component': 'weight_loss', 'cos_with_age': cos_weight_age},
        {'component': 'sum_shared',  'cos_with_age': cos_sum_age},
        {'component': 'difference',  'cos_with_age': cos_diff_age},
    ])
    summary.to_csv(INTERV_DIR / "shared_component_analysis.csv", index=False)
    print(f"  saved: {INTERV_DIR / 'shared_component_analysis.csv'}")


if __name__ == "__main__":
    main()

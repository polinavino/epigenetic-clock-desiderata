"""
03_svd_and_aging_direction.py

Computes the aging direction v* via SVD of the displacement matrix A,
runs the singular value dimensionality test, and runs the curvature test.

Outputs:
  data/interventions/aging_direction_v_star.parquet  -- v* (CpGs x 1)
  data/interventions/singular_values.parquet         -- all singular values
  data/interventions/curvature_test.csv              -- theta between age tertiles
  results/svd_spectrum.pdf                           -- eigenvalue spectrum plot
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    DISPLACEMENT_MATRIX, DISPLACEMENT_METADATA,
    AGING_DIRECTION, SINGULAR_VALUES, CURVATURE_TEST,
    INTERVENTION_DATASETS, INTERV_DIR,
    AGE_YOUNG_MAX, AGE_OLD_MIN,
    interv_beta_path, interv_meta_path,
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def run_svd(A):
    """
    Thin SVD of A (rows = intervention-context pairs, cols = CpGs).
    Returns U, s, Vt.
    v* is the first row of Vt (first right singular vector).
    """
    U, s, Vt = np.linalg.svd(A.values, full_matrices=False)
    return U, s, Vt


def dimensionality_test(s):
    """
    Print singular value spectrum and compute rho = lambda1 / lambda2.
    """
    rho = s[0] / s[1] if len(s) > 1 else np.inf
    print(f"\nSingular value spectrum:")
    for i, sv in enumerate(s[:10]):
        print(f"  lambda_{i+1}: {sv:.4f}")
    print(f"\n  rho = lambda_1 / lambda_2 = {rho:.3f}")
    if rho > 3:
        print("  >> One-dimensionality supported (rho > 3)")
    elif rho > 1.5:
        print("  >> Marginal one-dimensionality (1.5 < rho < 3); interpret with caution")
    else:
        print("  >> One-dimensionality NOT supported; aging signal is multi-dimensional")
    return rho


def plot_spectrum(s, rho, out_path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(range(1, len(s) + 1), s, color='steelblue', edgecolor='white')
    ax.axvline(x=1.5, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel("Singular value rank")
    ax.set_ylabel("Singular value")
    ax.set_title(f"Displacement matrix A — singular value spectrum\n"
                 f"$\\rho = \\lambda_1/\\lambda_2 = {rho:.2f}$")
    ax.set_xticks(range(1, len(s) + 1))
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  saved spectrum plot: {out_path}")


def curvature_test(v_star, cpg_ids):
    """
    Test whether aging direction rotates between young and old baseline states.

    For each intervention dataset with age metadata, split samples into
    young (age < AGE_YOUNG_MAX) and old (age > AGE_OLD_MIN) tertiles.
    Compute v* separately for each subset and measure the angle between them.
    """
    print(f"\nCurvature test (young: age < {AGE_YOUNG_MAX}, old: age > {AGE_OLD_MIN})")

    results = []

    for accession, config in INTERVENTION_DATASETS.items():
        meta_path = interv_meta_path(accession)
        beta_path = interv_beta_path(accession)

        if not meta_path.exists() or not beta_path.exists():
            continue

        meta_df = pd.read_csv(meta_path)
        if 'age' not in meta_df.columns:
            print(f"  {accession}: no age column, skipping curvature test")
            continue

        young_ids = meta_df.loc[meta_df['age'] < AGE_YOUNG_MAX, 'sample_id'].tolist()
        old_ids   = meta_df.loc[meta_df['age'] > AGE_OLD_MIN,   'sample_id'].tolist()

        if len(young_ids) < 5 or len(old_ids) < 5:
            print(f"  {accession}: insufficient age-stratified samples, skipping")
            continue

        # Import here to avoid circular dependency
        import h5py
        with h5py.File(beta_path, 'r') as f:
            sample_ids_all = [s.decode() for s in f['sample_ids'][:]]
            cpg_ids_all    = [c.decode() for c in f['cpg_ids'][:]]
            beta_all       = f['beta'][:]

        beta_df = pd.DataFrame(beta_all, index=sample_ids_all, columns=cpg_ids_all)

        # Restrict to shared CpGs
        shared = [c for c in cpg_ids if c in beta_df.columns]
        if len(shared) < 100:
            continue

        young_beta = beta_df.loc[[s for s in young_ids if s in beta_df.index], shared]
        old_beta   = beta_df.loc[[s for s in old_ids   if s in beta_df.index], shared]

        if len(young_beta) < 5 or len(old_beta) < 5:
            continue

        # Approximate v* for each age group as first PC of their beta matrices
        from sklearn.decomposition import PCA
        v_young = PCA(n_components=1).fit(young_beta.values).components_[0]
        v_old   = PCA(n_components=1).fit(old_beta.values).components_[0]

        # Angle between directions (take minimum of angle and pi - angle)
        cos_theta = np.dot(v_young, v_old) / (
            np.linalg.norm(v_young) * np.linalg.norm(v_old)
        )
        cos_theta = np.clip(cos_theta, -1, 1)
        theta = np.degrees(np.arccos(abs(cos_theta)))  # abs: directions are unsigned

        print(f"  {accession} ({config['tissue']}): theta = {theta:.1f} degrees")
        results.append({
            'accession': accession,
            'tissue':    config['tissue'],
            'n_young':   len(young_beta),
            'n_old':     len(old_beta),
            'theta_deg': theta,
        })

    if results:
        df = pd.DataFrame(results)
        df.to_csv(CURVATURE_TEST, index=False)
        print(f"\n  mean theta across datasets: {df['theta_deg'].mean():.1f} degrees")
        print(f"  saved: {CURVATURE_TEST}")
    else:
        print("  no datasets with age metadata available for curvature test")

    return results


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading displacement matrix A ...")
    if not DISPLACEMENT_MATRIX.exists():
        raise FileNotFoundError(
            f"{DISPLACEMENT_MATRIX} not found. Run script 02 first."
        )

    A = pd.read_parquet(DISPLACEMENT_MATRIX)
    meta_df = pd.read_csv(DISPLACEMENT_METADATA)
    cpg_ids = A.columns.tolist()

    print(f"  A shape: {A.shape} "
          f"({A.shape[0]} intervention-context pairs x {A.shape[1]:,} CpGs)")

    # SVD
    print("\nRunning SVD ...")
    U, s, Vt = run_svd(A)

    # Dimensionality test
    rho = dimensionality_test(s)

    # Plot spectrum
    plot_spectrum(s, rho, RESULTS_DIR / "svd_spectrum.pdf")

    # Save singular values
    sv_df = pd.DataFrame({'rank': range(1, len(s)+1), 'singular_value': s})
    sv_df.to_parquet(SINGULAR_VALUES)

    # Save v* (first right singular vector)
    v_star = pd.Series(Vt[0], index=cpg_ids, name='v_star')
    v_star.to_frame().to_parquet(AGING_DIRECTION)
    print(f"\n  v* saved: {AGING_DIRECTION}")
    print(f"  non-zero components: {(v_star.abs() > 1e-6).sum():,} / {len(v_star):,}")

    # Curvature test
    curvature_test(v_star.values, cpg_ids)

    print("\nDone.")


if __name__ == "__main__":
    main()

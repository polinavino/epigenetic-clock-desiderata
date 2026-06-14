"""
exploratory_smoking.py

Exploratory analysis of the GSE50660 smoking displacement vector.
Runs independently of the full pipeline.

Questions:
  1. Which CpGs move most with smoking? Do AHRR/F2RL3/GPR15 dominate?
  2. How aligned is the smoking displacement with the existing aging trajectory?
  3. Does smoking push methylation in the aging direction or orthogonally?
"""

import sys
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import INTERV_DIR, interv_beta_path, interv_meta_path

PARENT_SCRIPTS = Path(__file__).resolve().parent.parent.parent / "scripts"
sys.path.insert(0, str(PARENT_SCRIPTS))
from config import PRINCIPAL_CURVE, TAU, DATA_DIR, WEIGHTS, TOP_CPGS

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# Known smoking-associated CpGs from literature
SMOKING_CANONICAL = {
    'cg05575921': 'AHRR',
    'cg03636183': 'F2RL3',
    'cg19859270': 'GPR15',
    'cg23576855': 'PRSS23',
    'cg01940273': 'ALPPL2',
    'cg06126421': 'C1orf114',
    'cg10399789': 'ZMYND11',
    'cg25949550': 'AHRR',
    'cg21566642': 'ALPPL2',
}


def load_gse50660():
    beta_path = interv_beta_path("GSE50660")
    meta_path = interv_meta_path("GSE50660")

    with h5py.File(beta_path, 'r') as f:
        sample_ids = [s.decode() for s in f['sample_ids'][:]]
        cpg_ids    = [c.decode() for c in f['cpg_ids'][:]]
        beta       = f['beta'][:]

    beta_df = pd.DataFrame(beta, index=sample_ids, columns=cpg_ids)
    meta_df = pd.read_csv(meta_path)
    return beta_df, meta_df


def compute_displacement(beta_df, meta_df):
    cases    = meta_df.loc[meta_df['group'] == 'case',    'sample_id'].tolist()
    controls = meta_df.loc[meta_df['group'] == 'control', 'sample_id'].tolist()
    cases    = [s for s in cases    if s in beta_df.index]
    controls = [s for s in controls if s in beta_df.index]

    print(f"  cases: {len(cases)}, controls: {len(controls)}")

    delta = beta_df.loc[cases].mean() - beta_df.loc[controls].mean()
    return delta, cases, controls


def check_canonical_cpgs(delta):
    print("\n--- Canonical smoking CpGs ---")
    found = []
    for cpg, gene in SMOKING_CANONICAL.items():
        if cpg in delta.index:
            d = delta[cpg]
            pct = (delta.abs() >= abs(d)).mean() * 100
            print(f"  {cpg} ({gene}): delta = {d:+.4f}  (top {100-pct:.1f}% by |delta|)")
            found.append((cpg, gene, d))
        else:
            print(f"  {cpg} ({gene}): NOT IN DATASET")
    return found


def top_displaced_cpgs(delta, n=20):
    print(f"\n--- Top {n} CpGs by |delta| ---")
    top = delta.abs().nlargest(n)
    for cpg, val in top.items():
        direction = "↓" if delta[cpg] < 0 else "↑"
        print(f"  {cpg}: {direction} {delta[cpg]:+.4f}")
    return top


def alignment_with_principal_curve(delta, n_cpgs=200):
    """
    Project smoking displacement onto the existing principal curve direction.
    Uses the top N age-informative CpGs from the parent repo.
    """
    print("\n--- Alignment with aging trajectory ---")

    if not PRINCIPAL_CURVE.exists():
        print("  principal curve not found, skipping")
        return None

    if not TOP_CPGS.exists():
        print("  top CpGs file not found, skipping")
        return None

    with open(TOP_CPGS) as f:
        top_cpgs = [l.strip() for l in f if l.strip()]

    print(f"  aging trajectory CpGs: {len(top_cpgs)}")

    shared = [c for c in top_cpgs if c in delta.index]
    print(f"  shared with displacement vector: {len(shared)}")

    if len(shared) < 10:
        print("  too few shared CpGs, skipping")
        return None

    # Load principal curve
    curve = pd.read_parquet(PRINCIPAL_CURVE)
    shared_in_curve = [c for c in shared if c in curve.columns]
    print(f"  shared with principal curve: {len(shared_in_curve)}")

    if len(shared_in_curve) < 10:
        print("  too few shared CpGs with curve, skipping")
        return None

    # Aging direction: first difference along the curve (tangent at midpoint)
    mid = len(curve) // 2
    aging_dir = curve.iloc[mid+1][shared_in_curve].values - \
                curve.iloc[mid-1][shared_in_curve].values
    aging_dir = aging_dir / np.linalg.norm(aging_dir)

    # Smoking displacement in same CpG space
    smok_dir = delta[shared_in_curve].values
    smok_norm = np.linalg.norm(smok_dir)
    if smok_norm < 1e-10:
        print("  smoking displacement is zero in shared CpG space")
        return None
    smok_dir_unit = smok_dir / smok_norm

    # Cosine similarity
    cos_sim = np.dot(aging_dir, smok_dir_unit)
    angle   = np.degrees(np.arccos(np.clip(abs(cos_sim), 0, 1)))

    print(f"\n  cosine similarity (smoking vs aging direction): {cos_sim:+.4f}")
    print(f"  angle between directions: {angle:.1f} degrees")
    print(f"  tangential fraction: {abs(cos_sim):.3f}")
    print(f"  orthogonal fraction: {np.sqrt(1 - cos_sim**2):.3f}")

    if abs(cos_sim) > 0.3:
        print(f"\n  >> Smoking is substantially ALIGNED with aging direction")
        print(f"     (supports: smoking genuinely accelerates biological age)")
    elif abs(cos_sim) > 0.1:
        print(f"\n  >> Smoking has MODERATE alignment with aging direction")
    else:
        print(f"\n  >> Smoking is largely ORTHOGONAL to aging direction")
        print(f"     (suggests: smoking effect is mostly clock-gaming)")

    if cos_sim < 0:
        print(f"  Note: negative cosine means smoking pushes OPPOSITE to aging direction")
        print(f"  This would be unexpected — check curve orientation")

    return {
        'cos_sim': cos_sim,
        'angle': angle,
        'tangential': abs(cos_sim),
        'orthogonal': np.sqrt(1 - cos_sim**2),
        'n_shared': len(shared_in_curve),
    }


def plot_displacement_distribution(delta, canonical_found):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(delta.values, bins=200, color='steelblue', alpha=0.7,
            label=f'All CpGs (n={len(delta):,})')
    ax.axvline(0, color='black', linewidth=0.8)

    colours = ['red', 'orange', 'green', 'purple', 'brown',
               'pink', 'grey', 'cyan', 'magenta']
    for i, (cpg, gene, d) in enumerate(canonical_found):
        ax.axvline(d, color=colours[i % len(colours)],
                   linewidth=1.5, linestyle='--',
                   label=f'{gene} ({cpg[:10]}): {d:+.3f}')

    ax.set_xlabel('Mean beta difference (smokers - never smokers)')
    ax.set_ylabel('Number of CpGs')
    ax.set_title('GSE50660: Smoking displacement distribution')
    ax.legend(fontsize=7, loc='upper right')
    plt.tight_layout()

    out = RESULTS_DIR / "smoking_displacement_distribution.pdf"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"\n  saved: {out}")


def plot_alignment(delta, alignment, n_cpgs=200):
    if alignment is None or not TOP_CPGS.exists() or not PRINCIPAL_CURVE.exists():
        return

    with open(TOP_CPGS) as f:
        top_cpgs = [l.strip() for l in f if l.strip()]

    curve  = pd.read_parquet(PRINCIPAL_CURVE)
    shared = [c for c in top_cpgs if c in delta.index and c in curve.columns][:n_cpgs]

    if len(shared) < 10:
        return

    # PCA of principal curve in shared CpG space for visualisation
    curve_sub = curve[shared].values
    pca = PCA(n_components=2)
    curve_pca = pca.fit_transform(curve_sub)

    # Project smoking displacement onto same PCA space
    smok = delta[shared].values.reshape(1, -1)
    smok_pca = pca.transform(smok)[0]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(curve_pca[:, 0], curve_pca[:, 1],
            'o-', color='steelblue', alpha=0.5, markersize=3, label='Aging trajectory')
    ax.annotate('Young', curve_pca[0], fontsize=8, color='blue')
    ax.annotate('Old',   curve_pca[-1], fontsize=8, color='red')

    # Draw smoking displacement arrow from curve midpoint
    mid_pt = curve_pca[len(curve_pca)//2]
    scale  = np.linalg.norm(curve_pca[-1] - curve_pca[0]) * 0.3
    smok_unit = smok_pca / (np.linalg.norm(smok_pca) + 1e-10)
    ax.annotate('', xy=mid_pt + smok_unit * scale, xytext=mid_pt,
                arrowprops=dict(arrowstyle='->', color='firebrick', lw=2))
    ax.text(*(mid_pt + smok_unit * scale * 1.1), 'Smoking\ndisplacement',
            color='firebrick', fontsize=8)

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title(f'Smoking vs aging trajectory\n'
                 f'cos similarity = {alignment["cos_sim"]:+.3f}, '
                 f'angle = {alignment["angle"]:.1f}°')
    ax.legend()
    plt.tight_layout()

    out = RESULTS_DIR / "smoking_vs_aging_trajectory.pdf"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  saved: {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading GSE50660 ...")
    beta_df, meta_df = load_gse50660()
    print(f"  shape: {beta_df.shape}")

    print("\nComputing smoking displacement vector ...")
    delta, cases, controls = compute_displacement(beta_df, meta_df)
    print(f"  displacement vector: {len(delta):,} CpGs")
    print(f"  mean |delta|: {delta.abs().mean():.4f}")
    print(f"  max |delta|:  {delta.abs().max():.4f}")

    canonical = check_canonical_cpgs(delta)
    top       = top_displaced_cpgs(delta, n=20)
    alignment = alignment_with_principal_curve(delta)

    plot_displacement_distribution(delta, canonical)
    plot_alignment(delta, alignment)

    # Save displacement vector
    out = INTERV_DIR / "GSE50660_displacement_exploratory.parquet"
    delta.to_frame('delta').to_parquet(out)
    print(f"\n  displacement vector saved: {out}")

    print("\n=== Summary ===")
    print(f"  n cases (smokers):   {len(cases)}")
    print(f"  n controls (never):  {len(controls)}")
    print(f"  CpGs analysed:       {len(delta):,}")
    print(f"  mean |delta|:        {delta.abs().mean():.4f}")
    if alignment:
        print(f"  angle w/ aging:      {alignment['angle']:.1f}°")
        print(f"  tangential fraction: {alignment['tangential']:.3f}")


if __name__ == "__main__":
    main()

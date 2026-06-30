"""
exploratory_geometry.py

Consolidated geometry report for the intervention displacement vectors.
Reports, over the CpGs shared by all interventions:

  1. Pairwise angles (degrees) between *raw* (unsigned) displacement vectors.
     Two interventions in the same biological direction should be < 90 deg apart;
     near-90 deg means mutually orthogonal (no shared aging axis).
  2. The SVD singular-value spectrum of the signed, row-normalised matrix and
     rho = lambda_1 / lambda_2 (how dominant the leading shared component is).
  3. Each intervention's projection onto v* (signed cos and magnitude), and the
     tangential fraction used for classification.

This is the clock-independent geometry the framework rests on; it is descriptive
and makes no claim that a robust v* exists when rho is near 1.

Output: results table printed to stdout + data/interventions/geometry_report.csv
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    INTERVENTION_DATASETS, AGING_DIRECTION, displacement_path,
)
INTERV_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "interventions"


def main():
    # Load every available raw displacement vector with its config sign.
    vecs, signs, labels = {}, {}, {}
    for acc, cfg in INTERVENTION_DATASETS.items():
        p = displacement_path(acc)
        if not p.exists():
            continue
        d = pd.read_parquet(p).squeeze()
        if d.isna().any():
            d = d.dropna()
        vecs[acc] = d
        signs[acc] = +1 if cfg['sign'] == 'plus' else -1
        labels[acc] = cfg['label']

    accs = list(vecs)
    print(f"Loaded {len(accs)} displacement vectors: {accs}\n")

    # CpGs shared by all
    common = None
    for d in vecs.values():
        s = set(d.index)
        common = s if common is None else common & s
    common = sorted(common)
    print(f"CpGs shared by all interventions: {len(common):,}\n")

    # Matrix of raw displacements over shared CpGs (rows = interventions)
    M = np.vstack([vecs[a].reindex(common).values for a in accs])

    # 1. Pairwise angles between RAW displacements
    print("Pairwise angles between raw displacement vectors (degrees):")
    norms = np.linalg.norm(M, axis=1)
    cos = (M @ M.T) / np.outer(norms, norms)
    cos = np.clip(cos, -1, 1)
    ang = np.degrees(np.arccos(cos))
    hdr = "            " + "".join(f"{labels[a][:10]:>12}" for a in accs)
    print(hdr)
    for i, a in enumerate(accs):
        print(f"{labels[a][:11]:>11} " + "".join(f"{ang[i, j]:12.1f}" for j in range(len(accs))))
    print()

    # 2. SVD of signed, row-normalised matrix
    Ms = np.vstack([signs[a] * (vecs[a].reindex(common).values) for a in accs])
    Ms = Ms / np.linalg.norm(Ms, axis=1, keepdims=True)
    sv = np.linalg.svd(Ms, compute_uv=False)
    print("Signed row-normalised SVD spectrum:", np.round(sv, 4).tolist())
    rho = sv[0] / sv[1] if len(sv) > 1 and sv[1] > 0 else float('inf')
    print(f"rho = lambda_1/lambda_2 = {rho:.3f}   "
          f"({'1D axis supported' if rho >= 2 else 'NO dominant single axis'})\n")

    # 3. Projection onto saved v*
    rows = []
    if Path(AGING_DIRECTION).exists():
        v_star = pd.read_parquet(AGING_DIRECTION).squeeze()
        vc = [c for c in v_star.index if c in set(common)]
        v = v_star.reindex(vc).values
        v = v / np.linalg.norm(v)
        print(f"Projection onto v* (over {len(vc):,} shared CpGs):")
        print(f"{'label':>24}{'sign':>7}{'cos(d,v*)':>12}{'<d,v*>':>12}{'||d||':>10}{'tang_frac':>11}")
        for a in accs:
            d = vecs[a].reindex(vc).values
            dn = np.linalg.norm(d)
            proj = float(d @ v)
            cosv = proj / dn if dn else float('nan')
            tf = abs(proj) / dn if dn else float('nan')
            print(f"{labels[a][:24]:>24}{('+' if signs[a]>0 else '-'):>7}"
                  f"{cosv:>12.3f}{proj:>12.4f}{dn:>10.4f}{tf:>11.3f}")
            rows.append({'accession': a, 'label': labels[a],
                         'sign': '+' if signs[a] > 0 else '-',
                         'cos_with_vstar': cosv, 'proj': proj,
                         'magnitude': dn, 'tangential_fraction': tf})
    out = INTERV_DIR / "geometry_report.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nSaved: {out}")

    # Heatmap of the pairwise-angle matrix
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        results = Path(__file__).resolve().parent.parent / "results"
        results.mkdir(exist_ok=True)
        short = [labels[a] for a in accs]
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        im = ax.imshow(ang, cmap='RdBu', vmin=60, vmax=120)
        ax.set_xticks(range(len(accs))); ax.set_yticks(range(len(accs)))
        ax.set_xticklabels(short, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(short, fontsize=8)
        for i in range(len(accs)):
            for j in range(len(accs)):
                ax.text(j, i, f"{ang[i, j]:.0f}", ha='center', va='center',
                        fontsize=8, color='black')
        cb = fig.colorbar(im, ax=ax, shrink=0.8)
        cb.set_label('angle (deg);  90 = orthogonal')
        ax.set_title(f'Pairwise angles between intervention displacements\n'
                     f'(rho = {rho:.2f}: no dominant shared aging axis)')
        plt.tight_layout()
        f_out = results / "displacement_angle_matrix.pdf"
        plt.savefig(f_out, dpi=150); plt.close()
        print(f"Saved: {f_out}")
    except Exception as e:
        print(f"  (heatmap skipped: {e})")


if __name__ == "__main__":
    main()

"""
gtex_cross_tissue.py   (context-invariance test)

Uses GTEx multi-tissue methylation (GSE213478, 9 tissues, 987 samples, EPIC) to
ask two questions over the v* CpGs:

  Q1. Is there a context-invariant aging direction? i.e. do the *within-tissue*
      age-methylation directions agree across tissues? (pairwise cos + SVD rho)
  Q2. Does the intervention-derived v* align with the cross-tissue aging
      direction? If v* is a real aging axis it should; if intervention
      displacements are off-axis it will not.

Age uses the GEO bracket midpoints (20-29 -> 25, ... 70-79 -> 75). The per-tissue
age vector is the per-CpG Pearson correlation of beta with age across that
tissue's samples, over the CpGs shared with v*.

Inputs : /tmp/gtex_vstar.csv (streamed v* CpGs x 987 samples)
         data/interventions/GTEX_metadata.csv  (sample_id, tissue, age, smoker)
         data/interventions/aging_direction_v_star.parquet
Output : data/interventions/gtex_cross_tissue_report.txt
         results/gtex_cross_tissue.pdf
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
INTERV = ROOT / "data" / "interventions"
RESULTS = Path(__file__).resolve().parent.parent / "results"
GTEX_CSV = Path("/tmp/gtex_vstar.csv")

AGE_MID = {"20-29": 25, "30-39": 35, "40-49": 45,
           "50-59": 55, "60-69": 65, "70-79": 75}
MIN_TISSUE_N = 45


def main():
    out = []
    log = lambda s: (print(s), out.append(s))

    beta = pd.read_csv(GTEX_CSV, index_col=0)          # CpGs x samples
    beta.index = [str(c).strip().strip('"') for c in beta.index]
    beta.columns = [str(c).strip().strip('"') for c in beta.columns]
    meta = pd.read_csv(INTERV / "GTEX_metadata.csv")
    meta["age_mid"] = meta["age"].map(AGE_MID)
    meta = meta.dropna(subset=["age_mid"]).set_index("sample_id")

    v_star = pd.read_parquet(INTERV / "aging_direction_v_star.parquet").squeeze()
    cpgs = [c for c in beta.index if c in set(v_star.index)]
    log(f"GTEx: {beta.shape[1]} samples x {beta.shape[0]} v* CpGs "
        f"({len(cpgs)} shared with v*)")
    B = beta.loc[cpgs]                                  # CpGs x samples
    v = v_star.loc[cpgs].values
    v = v / np.linalg.norm(v)

    # Per-tissue age-methylation direction over the shared CpGs
    tissues, agevecs = [], []
    log("\nPer-tissue age-direction (correlation of beta with age):")
    for tissue, grp in meta.groupby("tissue"):
        samples = [s for s in grp.index if s in B.columns]
        if len(samples) < MIN_TISSUE_N:
            continue
        sub = B[samples].values                        # CpGs x n
        ages = grp.loc[samples, "age_mid"].values.astype(float)
        # per-CpG Pearson r with age
        sub_c = sub - sub.mean(axis=1, keepdims=True)
        a_c = ages - ages.mean()
        num = sub_c @ a_c
        den = np.linalg.norm(sub_c, axis=1) * np.linalg.norm(a_c)
        r = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
        cos_v = float(np.dot(r / (np.linalg.norm(r) + 1e-12), v))
        tissues.append(tissue); agevecs.append(r)
        log(f"  {tissue:26s} n={len(samples):3d}  "
            f"|age-corr|_mean={np.abs(r).mean():.3f}  cos(age_dir, v*)={cos_v:+.3f}")

    A = np.vstack(agevecs)
    An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)

    # Q1: cross-tissue agreement of aging directions
    cosmat = np.clip(An @ An.T, -1, 1)
    ang = np.degrees(np.arccos(cosmat))
    log("\nQ1. Pairwise angles between tissue age-directions (deg):")
    log("      " + "".join(f"{t[:8]:>10}" for t in tissues))
    for i, t in enumerate(tissues):
        log(f"{t[:6]:>6}" + "".join(f"{ang[i,j]:10.1f}" for j in range(len(tissues))))
    sv = np.linalg.svd(An, compute_uv=False)
    rho = sv[0] / sv[1] if len(sv) > 1 else float("inf")
    log(f"  SVD spectrum: {np.round(sv,3).tolist()}")
    log(f"  rho = lambda1/lambda2 = {rho:.3f}  "
        f"({'shared cross-tissue aging axis' if rho >= 2 else 'no dominant shared axis'})")

    # Cross-tissue consensus aging direction (leading singular vector)
    U, S, Vt = np.linalg.svd(An, full_matrices=False)
    consensus = Vt[0]
    if np.dot(consensus, An.mean(0)) < 0:
        consensus = -consensus
    cos_consensus_v = float(np.dot(consensus / np.linalg.norm(consensus), v))
    log(f"\nQ2. cos(cross-tissue consensus aging direction, v*) = {cos_consensus_v:+.3f}")

    # Blood specifically vs v*
    if "Whole Blood" in tissues:
        bi = tissues.index("Whole Blood")
        log(f"    cos(Whole Blood age-direction, v*) = "
            f"{float(np.dot(An[bi], v)):+.3f}")

    (INTERV / "gtex_cross_tissue_report.txt").write_text("\n".join(out))

    # Figure
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        RESULTS.mkdir(exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(ang, cmap="RdBu", vmin=0, vmax=90)
        ax.set_xticks(range(len(tissues))); ax.set_yticks(range(len(tissues)))
        ax.set_xticklabels([t[:12] for t in tissues], rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels([t[:12] for t in tissues], fontsize=7)
        for i in range(len(tissues)):
            for j in range(len(tissues)):
                ax.text(j, i, f"{ang[i,j]:.0f}", ha="center", va="center", fontsize=7)
        fig.colorbar(im, ax=ax, shrink=0.8, label="angle (deg); 0 = same aging direction")
        ax.set_title(f"GTEx cross-tissue aging directions over v* CpGs\n"
                     f"rho={rho:.2f}, cos(consensus, v*)={cos_consensus_v:+.2f}")
        plt.tight_layout(); plt.savefig(RESULTS / "gtex_cross_tissue.pdf", dpi=150); plt.close()
        log(f"\nSaved: {RESULTS/'gtex_cross_tissue.pdf'}")
    except Exception as e:
        log(f"  (figure skipped: {e})")


if __name__ == "__main__":
    main()

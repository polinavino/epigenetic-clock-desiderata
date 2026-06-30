"""
cell_composition_analysis.py   (addresses Vulnerability 1)

The core worry: the near-orthogonality of intervention displacement vectors in
blood could be an artifact of intervention-induced shifts in leukocyte
composition (chemo depletes lymphocytes, weight loss/exercise change
inflammation, smoking shifts granulocytes) rather than a statement about an
aging axis.

This script, for each blood intervention:
  1. estimates 7 leukocyte fractions per sample by constrained projection onto
     the EpiDISH centDHSbloodDMC reference (NNLS, normalised to sum 1),
  2. reports the intervention-induced fraction change (post-pre or case-control)
     -> the magnitude of the confound,
  3. recomputes the displacement vector after residualising beta on the
     estimated fractions (within dataset),
  4. compares the pairwise-angle matrix and signed-SVD rho BEFORE vs AFTER
     cell-composition adjustment.

If orthogonality (angles ~90, rho ~1) survives adjustment, it is not a pure
cell-mix artifact.

Output: data/interventions/cell_fractions_{ACC}.csv
        data/interventions/cell_adjustment_summary.txt
"""

import sys
import numpy as np
import pandas as pd
import h5py
from pathlib import Path
from scipy.optimize import nnls

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import INTERVENTION_DATASETS, interv_beta_path, interv_meta_path

SCRIPT_DIR = Path(__file__).resolve().parent
INTERV_DIR = SCRIPT_DIR.parent.parent / "data" / "interventions"
REF_CSV = SCRIPT_DIR / "centDHSbloodDMC_ref.csv"
BLOOD = ["GSE50660", "GSE272137", "GSE328810", "GSE240184", "GSE140038"]


def load_h5(path):
    with h5py.File(path, "r") as f:
        beta = f["beta"][:]
        sids = [s.decode() for s in f["sample_ids"][:]]
        cids = [c.decode() for c in f["cpg_ids"][:]]
    return pd.DataFrame(beta, index=sids, columns=cids)


def estimate_fractions(beta_df, ref):
    """Constrained projection: per sample solve nnls(ref_shared, beta_shared)."""
    shared = [c for c in ref.index if c in beta_df.columns]
    R = ref.loc[shared].values            # CpGs x 7
    B = beta_df[shared].values            # samples x CpGs
    fracs = np.zeros((B.shape[0], R.shape[1]))
    for i in range(B.shape[0]):
        w, _ = nnls(R, B[i])
        s = w.sum()
        fracs[i] = w / s if s > 0 else w
    return pd.DataFrame(fracs, index=beta_df.index, columns=ref.columns), len(shared)


def residualise(beta_df, fracs):
    """Remove the linear component of beta explained by (centered) fractions."""
    Fc = fracs.loc[beta_df.index].values
    Fc = Fc - Fc.mean(axis=0, keepdims=True)
    B = beta_df.values
    coef, *_ = np.linalg.lstsq(Fc, B, rcond=None)   # 7 x CpGs
    B_adj = B - Fc @ coef
    return pd.DataFrame(B_adj, index=beta_df.index, columns=beta_df.columns)


def displacement(beta_df, meta, design):
    """Return (raw displacement Series, group/timepoint fraction-change helper)."""
    if design == "longitudinal":
        pre = meta[meta["timepoint"] == "pre"].set_index("subject_id")["sample_id"]
        post = meta[meta["timepoint"] == "post"].set_index("subject_id")["sample_id"]
        subj = [s for s in pre.index if s in post.index]
        pre_s = [pre[s] for s in subj if pre[s] in beta_df.index and post[s] in beta_df.index]
        post_s = [post[s] for s in subj if pre[s] in beta_df.index and post[s] in beta_df.index]
        d = (beta_df.loc[post_s].values - beta_df.loc[pre_s].values).mean(axis=0)
        return pd.Series(d, index=beta_df.columns), (pre_s, post_s)
    else:  # cross_sectional
        case = meta[meta["group"] == "case"]["sample_id"]
        ctrl = meta[meta["group"] == "control"]["sample_id"]
        case = [s for s in case if s in beta_df.index]
        ctrl = [s for s in ctrl if s in beta_df.index]
        d = beta_df.loc[case].mean(axis=0).values - beta_df.loc[ctrl].mean(axis=0).values
        return pd.Series(d, index=beta_df.columns), (ctrl, case)


def angle_matrix(vectors, labels):
    accs = list(vectors)
    common = None
    for d in vectors.values():
        common = set(d.index) if common is None else common & set(d.index)
    common = sorted(common)
    M = np.vstack([vectors[a].reindex(common).values for a in accs])
    norms = np.linalg.norm(M, axis=1)
    cos = np.clip((M @ M.T) / np.outer(norms, norms), -1, 1)
    return np.degrees(np.arccos(cos)), accs, len(common)


def signed_rho(vectors, signs):
    accs = list(vectors)
    common = None
    for d in vectors.values():
        common = set(d.index) if common is None else common & set(d.index)
    common = sorted(common)
    Ms = np.vstack([signs[a] * vectors[a].reindex(common).values for a in accs])
    Ms = Ms / np.linalg.norm(Ms, axis=1, keepdims=True)
    sv = np.linalg.svd(Ms, compute_uv=False)
    return sv, (sv[0] / sv[1] if len(sv) > 1 else float("inf"))


def main():
    ref = pd.read_csv(REF_CSV, index_col=0)
    raw_disp, adj_disp, signs, labels = {}, {}, {}, {}
    out = []

    def log(s):
        print(s); out.append(s)

    for acc in BLOOD:
        cfg = INTERVENTION_DATASETS[acc]
        beta = load_h5(interv_beta_path(acc))
        meta = pd.read_csv(interv_meta_path(acc))
        meta["sample_id"] = meta["sample_id"].astype(str)
        labels[acc] = cfg["label"]
        signs[acc] = +1 if cfg["sign"] == "plus" else -1

        fracs, n_ref = estimate_fractions(beta, ref)
        fracs.to_csv(INTERV_DIR / f"cell_fractions_{acc}.csv")

        d_raw, (grp0, grp1) = displacement(beta, meta, cfg["design"])
        beta_adj = residualise(beta, fracs)
        d_adj, _ = displacement(beta_adj, meta, cfg["design"])
        raw_disp[acc], adj_disp[acc] = d_raw, d_adj

        # fraction change (intervention - baseline / case - control)
        fc = fracs.loc[grp1].mean() - fracs.loc[grp0].mean()
        log(f"\n{acc} ({cfg['label']}, {cfg['design']}, ref CpGs used={n_ref}/333)")
        log("  mean leukocyte fractions: " +
            ", ".join(f"{c}={fracs[c].mean():.3f}" for c in ref.columns))
        log("  intervention-induced fraction change (grp1-grp0): " +
            ", ".join(f"{c}={fc[c]:+.3f}" for c in ref.columns))
        log(f"  |raw displacement|={d_raw.abs().mean():.5f}  "
            f"|cell-adjusted|={d_adj.abs().mean():.5f}  "
            f"cos(raw,adj)={np.dot(d_raw.values, d_adj.reindex(d_raw.index).values)/(np.linalg.norm(d_raw)*np.linalg.norm(d_adj)):.3f}")

    # Compare angle matrices and rho
    for tag, vecs in [("RAW", raw_disp), ("CELL-ADJUSTED", adj_disp)]:
        ang, accs, ncpg = angle_matrix(vecs, labels)
        sv, rho = signed_rho(vecs, signs)
        log(f"\n===== {tag} pairwise angles (deg), shared CpGs={ncpg:,} =====")
        log("            " + "".join(f"{labels[a][:9]:>11}" for a in accs))
        for i, a in enumerate(accs):
            log(f"{labels[a][:11]:>11} " + "".join(f"{ang[i,j]:11.1f}" for j in range(len(accs))))
        log(f"  signed-SVD spectrum: {np.round(sv,3).tolist()}")
        log(f"  rho = lambda1/lambda2 = {rho:.3f}")

    (INTERV_DIR / "cell_adjustment_summary.txt").write_text("\n".join(out))
    print(f"\nSaved fractions + summary to {INTERV_DIR}")


if __name__ == "__main__":
    main()

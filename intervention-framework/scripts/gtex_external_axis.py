"""
gtex_external_axis.py

Classify interventions against an EXTERNAL, clock-independent aging reference:
the GTEx whole-blood age-methylation direction (per-CpG correlation of beta with
chronological age in GTEx whole blood), instead of the framework-internal v*
(which is non-robust, rho ~ 1.3). This anchors the geroprotective / accelerating
/ clock-gaming call outside the intervention set and outside any clock.

Also runs two placebo controls to check the GTEx findings are not artifacts of
v*-CpG selection:
  P1  cross-tissue aging tissue-specificity (rho) and blood-lung angle on a random
      CpG set vs the v* CpGs
  P2  intervention displacements projected onto a *random-CpG* GTEx "aging"
      direction (a null axis) — should not separate accelerators from geroprotectors

Inputs : /tmp/gtex_ext.csv  (S union random CpGs x 987 samples, streamed)
         /tmp/S_cpgs.txt (intervention-shared), /tmp/V_cpgs.txt (v*), /tmp/rand5k_cpgs.txt
         data/interventions/GTEX_metadata.csv, aging_direction_v_star.parquet
         data/interventions/{ACC}_displacement.parquet (raw, unsigned)
Output : data/interventions/gtex_external_axis_report.txt
         data/interventions/intervention_classification_external.csv
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import INTERVENTION_DATASETS, displacement_path

INTERV = Path(__file__).resolve().parent.parent.parent / "data" / "interventions"
AGE_MID = {"20-29":25,"30-39":35,"40-49":45,"50-59":55,"60-69":65,"70-79":75}
ACCS = ["GSE50660","GSE53045","GSE140038","GSE89218",
        "GSE272137","GSE240184","GSE328810","GSE193730"]
CLOCK_GAMING_EPS = 0.30   # |cos| below this -> off-axis / clock-gaming


def age_direction(beta, meta, cpgs, tissue):
    """Per-CpG Pearson corr(beta, age) over samples of `tissue`. Returns Series."""
    samples = [s for s in meta.index[meta["tissue"] == tissue] if s in beta.columns]
    sub = beta.loc[cpgs, samples].values
    ages = meta.loc[samples, "age_mid"].values.astype(float)
    sc = sub - sub.mean(1, keepdims=True)
    ac = ages - ages.mean()
    den = np.linalg.norm(sc, axis=1) * np.linalg.norm(ac)
    r = np.divide(sc @ ac, den, out=np.zeros(sub.shape[0]), where=den > 0)
    return pd.Series(r, index=cpgs)


def main():
    out = []
    log = lambda s: (print(s), out.append(s))

    beta = pd.read_csv("/tmp/gtex_ext.csv", index_col=0)
    beta.index = [str(c).strip().strip('"') for c in beta.index]
    beta.columns = [str(c).strip().strip('"') for c in beta.columns]
    meta = pd.read_csv(INTERV / "GTEX_metadata.csv")
    meta["age_mid"] = meta["age"].map(AGE_MID)
    meta = meta.dropna(subset=["age_mid"]).set_index("sample_id")

    S = [c for c in (l.strip() for l in open("/tmp/S_cpgs.txt")) if c in beta.index]
    V = [c for c in (l.strip() for l in open("/tmp/V_cpgs.txt")) if c in beta.index]
    R = [c for c in (l.strip() for l in open("/tmp/rand5k_cpgs.txt")) if c in beta.index]
    log(f"GTEx betas: {beta.shape[0]} CpGs x {beta.shape[1]} samples")
    log(f"  intervention-shared S={len(S)}, v*={len(V)}, random={len(R)}\n")

    # ---- External aging reference: GTEx whole-blood age-direction over S ----
    aging = age_direction(beta, meta, S, "Whole Blood")
    a = aging.values
    a_unit = a / (np.linalg.norm(a) + 1e-12)
    log("== External classification vs GTEx whole-blood aging direction (over S) ==")
    log(f"{'label':>24}{'sign':>6}{'cos':>9}{'proj':>10}{'class':>16}{'match':>7}")
    rows = []
    for acc in ACCS:
        cfg = INTERVENTION_DATASETS[acc]
        d = pd.read_parquet(displacement_path(acc)).squeeze()      # raw, unsigned
        sh = [c for c in S if c in d.index]
        dv = d.reindex(sh).values
        av = aging.reindex(sh).values
        av = av / (np.linalg.norm(av) + 1e-12)
        proj = float(dv @ av)                       # + = toward OLDER (accelerating)
        cos = proj / (np.linalg.norm(dv) + 1e-12)
        if abs(cos) < CLOCK_GAMING_EPS:
            cls = "clock_gaming"
        else:
            cls = "age_accelerating" if proj > 0 else "geroprotective"
        expected = "age_accelerating" if cfg["sign"] == "plus" else "geroprotective"
        match = cls == expected
        rows.append({"accession": acc, "label": cfg["label"], "expected_sign": cfg["sign"],
                     "cos_with_gtex_aging": cos, "projection": proj,
                     "classification": cls, "matches_expected": match})
        log(f"{cfg['label'][:24]:>24}{('+' if cfg['sign']=='plus' else '-'):>6}"
            f"{cos:>9.3f}{proj:>10.4f}{cls:>16}{('OK' if match else 'X'):>7}")
    pd.DataFrame(rows).to_csv(INTERV / "intervention_classification_external.csv", index=False)
    n_ok = sum(r["matches_expected"] for r in rows)
    log(f"\n  {n_ok}/8 match expected sign against the EXTERNAL GTEx blood-aging axis "
        f"(vs 3/8 against internal v*).")

    # consistency: v* vs GTEx blood aging on v* CpGs
    vstar = pd.read_parquet(INTERV / "aging_direction_v_star.parquet").squeeze()
    vv = [c for c in V if c in vstar.index]
    aging_V = age_direction(beta, meta, vv, "Whole Blood").values
    vs = vstar.reindex(vv).values
    cos_v = float(vs @ aging_V / (np.linalg.norm(vs)*np.linalg.norm(aging_V) + 1e-12))
    log(f"  (consistency) cos(v*, GTEx blood aging) over v* CpGs = {cos_v:+.3f}")

    # ---- Placebo P1: cross-tissue tissue-specificity, random vs v* CpGs ----
    log("\n== P1. Cross-tissue aging tissue-specificity: v* CpGs vs random CpGs ==")
    tissues = [t for t, g in meta.groupby("tissue") if (meta["tissue"] == t).sum() >= 45]
    for name, cpgs in [("v* CpGs", V), ("random CpGs", R)]:
        A = np.vstack([age_direction(beta, meta, cpgs, t).values for t in tissues])
        An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
        sv = np.linalg.svd(An, compute_uv=False)
        rho = sv[0] / sv[1]
        bi, li = tissues.index("Whole Blood"), tissues.index("Lung")
        bl = np.degrees(np.arccos(np.clip(An[bi] @ An[li], -1, 1)))
        log(f"  {name:12s}: cross-tissue rho={rho:.3f}  blood-lung angle={bl:.1f} deg  "
            f"|blood age-corr|_mean={np.abs(A[bi]).mean():.3f}")

    # ---- Placebo P2: interventions vs a random-CpG null aging axis ----
    log("\n== P2. Interventions projected onto a RANDOM-CpG GTEx aging axis (null) ==")
    aging_R = age_direction(beta, meta, R, "Whole Blood")
    seps = []
    for acc in ACCS:
        cfg = INTERVENTION_DATASETS[acc]
        d = pd.read_parquet(displacement_path(acc)).squeeze()
        sh = [c for c in R if c in d.index]
        if len(sh) < 100:
            continue
        dv = d.reindex(sh).values
        av = aging_R.reindex(sh).values; av = av/(np.linalg.norm(av)+1e-12)
        cos = float(dv @ av / (np.linalg.norm(dv)+1e-12))
        seps.append((cfg["sign"], cos))
    plus = [c for s, c in seps if s == "plus"]
    minus = [c for s, c in seps if s == "minus"]
    log(f"  mean cos onto random axis: accelerators={np.mean(plus):+.3f}, "
        f"geroprotectors={np.mean(minus):+.3f} (expect both ~0 if null)")

    (INTERV / "gtex_external_axis_report.txt").write_text("\n".join(out))
    log(f"\nSaved: {INTERV/'gtex_external_axis_report.txt'}")


if __name__ == "__main__":
    main()

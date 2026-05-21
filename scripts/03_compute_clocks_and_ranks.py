"""
03_compute_clocks_and_ranks.py

Computes all five clock outputs for every sample, then performs the
rank consistency analysis that is the core empirical contribution of the paper.

Steps:
    1. Load beta matrices and metadata
    2. Compute all clock outputs using biolearn
    3. Compute age acceleration residuals (clock output - expected for age)
    4. Compute rank consistency matrix kappa(k, k') for all clock pairs
    5. Decompose rank reversals into tau-driven vs r-driven
    6. Identify high-instability individuals
    7. Compute directional stability (D2) for each clock
    8. Save all outputs and figures

Outputs (all in data/):
    clock_outputs.parquet     -- raw clock output per sample per clock
    age_acceleration.parquet  -- age-acceleration residuals per sample per clock
    rank_consistency.parquet  -- kappa matrix (clocks x clocks)
    rank_reversals.parquet    -- per-pair rank reversal analysis
    instability_scores.parquet -- per-sample instability score

Figures (in paper/figures/):
    03_rank_consistency_matrix.png
    03_tau_r_decomposition.png
    03_high_instability_individuals.png
    03_directional_stability.png
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr, kendalltau, linregress
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings("ignore")

# ── Config ────────────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    DATA_DIR, GSE40279_BETA, GSE87571_BETA, COMMON_CPGS,
    METADATA, TAU, RESIDUALS, WEIGHTS, CLOCK_OUTPUTS, FIGURES_DIR,
    RANDOM_SEED,
)

np.random.seed(RANDOM_SEED)
FIGURES_DIR.mkdir(exist_ok=True)

# Output paths (not in config yet -- add here)
AGE_ACCEL       = DATA_DIR / "age_acceleration.parquet"
RANK_CONSISTENCY = DATA_DIR / "rank_consistency.parquet"
RANK_REVERSALS  = DATA_DIR / "rank_reversals.parquet"
INSTABILITY     = DATA_DIR / "instability_scores.parquet"

# ── Clock definitions ─────────────────────────────────────────────────────────
# Each clock has a name, biolearn model ID, and declared output type.
# Types follow our framework:
#   T_tau:   position clock (estimates tau directly)
#   T_delta: deviation clock (estimates tau + f(r))
#   T_rate:  rate clock (estimates d_tau/dt)

CLOCKS = {
    "Horvath":    {"model": "Horvathv1",  "type": "T_tau"},
    "Hannum":     {"model": "Hannum",     "type": "T_tau"},
    "PhenoAge":   {"model": "PhenoAge",   "type": "T_delta"},
    "GrimAge":    {"model": "GrimAgeV1",  "type": "T_delta"},
    "DunedinPACE":{"model": "DunedinPACE","type": "T_rate"},
}

# ── Step 1: Load data ─────────────────────────────────────────────────────────
print("=== Step 1: Loading data ===")

metadata = pd.read_csv(METADATA)
tau_df   = pd.read_parquet(TAU)
print(f"  Samples: {len(metadata)}")
print(f"  Tau range: [{tau_df['tau'].min():.1f}, {tau_df['tau'].max():.1f}]")

# Load beta matrices -- biolearn needs CpGs x samples orientation
print("\n  Loading GSE40279...")
beta_40 = pd.read_hdf(str(GSE40279_BETA), key="beta")   # CpGs x samples
print(f"    Shape: {beta_40.shape}")

print("  Loading GSE87571...")
beta_87 = pd.read_hdf(str(GSE87571_BETA), key="beta")   # CpGs x samples
print(f"    Shape: {beta_87.shape}")

# Load clock CpG lists from biolearn to restrict memory usage
# We only need the union of CpGs used by our four clocks
from biolearn.model_gallery import ModelGallery as _MG
_gallery = _MG()
clock_cpgs = set()
for _cname, _cinfo in CLOCKS.items():
    try:
        _m = _gallery.get(_cinfo["model"])
        _sites = _m.clock.methylation_sites()
        clock_cpgs.update(_sites)
        print(f"  {_cname}: {len(_sites)} CpG sites")
    except Exception as _e:
        print(f"  {_cname}: could not get sites ({_e})")

print(f"  Union of clock CpGs: {len(clock_cpgs)}")

# Restrict to clock CpGs present in both datasets
common_cpgs_set = (set(beta_40.index) & set(beta_87.index)) & clock_cpgs
print(f"  Clock CpGs present in both datasets: {len(common_cpgs_set)}")
beta_40 = beta_40.loc[sorted(common_cpgs_set)]
beta_87 = beta_87.loc[sorted(common_cpgs_set)]

# Align GSE87571 columns to metadata sample IDs
# (positional alignment as in script 1)
import re
meta_87 = metadata[metadata["dataset"] == "GSE87571"].reset_index(drop=True)
beta_87_cols = beta_87.columns.tolist()

def sort_key(c):
    m = re.match(r"X(\d+)$", c)
    return int(m.group(1)) if m else 0

beta_87_sorted = sorted(beta_87_cols, key=sort_key)
n = min(len(beta_87_sorted), len(meta_87))
col_to_id = dict(zip(beta_87_sorted[:n], meta_87["sample_id"].tolist()[:n]))
beta_87 = beta_87.rename(columns=col_to_id)

# Keep only samples in metadata
meta_40_ids = metadata[metadata["dataset"] == "GSE40279"]["sample_id"].tolist()
meta_87_ids = metadata[metadata["dataset"] == "GSE87571"]["sample_id"].tolist()

beta_40 = beta_40[[c for c in meta_40_ids if c in beta_40.columns]]
beta_87 = beta_87[[c for c in meta_87_ids if c in beta_87.columns]]

# Combine
beta_all = pd.concat([beta_40, beta_87], axis=1)   # CpGs x samples
print(f"\n  Combined beta: {beta_all.shape} (CpGs x samples)")

# Align metadata to beta columns
sample_ids = beta_all.columns.tolist()
meta_aligned = metadata.set_index("sample_id").loc[sample_ids]
ages = meta_aligned["age"].values
datasets = meta_aligned["dataset"].values

# ── Step 2: Compute clock outputs ─────────────────────────────────────────────
print("\n=== Step 2: Computing clock outputs ===")

from biolearn.model_gallery import ModelGallery
from biolearn.data_library import GeoData

gallery = ModelGallery()

# Harmonize sex encoding to biolearn format: 0=Female, 1=Male
def harmonize_sex(val):
    if str(val).strip().lower() in ("female", "f"):
        return 0
    elif str(val).strip().lower() in ("male", "m"):
        return 1
    else:
        return np.nan

sex_raw = meta_aligned["sex"].values
sex_encoded = np.array([harmonize_sex(v) for v in sex_raw], dtype=float)
n_missing_sex = np.isnan(sex_encoded).sum()
print(f"  Sex encoding: {(sex_encoded==0).sum()} female, {(sex_encoded==1).sum()} male, {n_missing_sex} missing")

# Build GeoData object with age and sex
geo_metadata = pd.DataFrame({
    "age": ages,
    "sex": sex_encoded,
}, index=sample_ids)
geo_data = GeoData(metadata=geo_metadata, dnam=beta_all)

clock_results = {}
for clock_name, clock_info in CLOCKS.items():
    print(f"  Computing {clock_name}...")
    try:
        model = gallery.get(clock_info["model"])
        result = model.predict(geo_data)
        # result is a DataFrame with index=sample_ids
        if isinstance(result, pd.DataFrame):
            # biolearn returns DataFrame with predicted age in first column
            vals = result.iloc[:, 0].values
        else:
            vals = result.values
        clock_results[clock_name] = pd.Series(vals, index=sample_ids)
        valid = ~np.isnan(vals)
        print(f"    {clock_name}: {valid.sum()} valid predictions, "
              f"mean={np.nanmean(vals):.1f}, range=[{np.nanmin(vals):.1f}, {np.nanmax(vals):.1f}]")
    except Exception as e:
        print(f"    {clock_name}: ERROR -- {e}")
        clock_results[clock_name] = pd.Series(np.nan, index=sample_ids)

# Combine into DataFrame: samples x clocks
clock_df = pd.DataFrame(clock_results)
clock_df.index.name = "sample_id"
clock_df.to_parquet(CLOCK_OUTPUTS)
print(f"\n  Saved: {CLOCK_OUTPUTS}")
print(f"  Clock output shape: {clock_df.shape}")
print(f"\n  Correlations with chronological age:")
for col in clock_df.columns:
    valid = clock_df[col].notna() & ~np.isnan(ages)
    if valid.sum() < 10:
        print(f"    {col}: insufficient valid predictions")
        continue
    r, _ = pearsonr(clock_df.loc[valid, col], ages[valid])
    print(f"    {col}: r={r:.3f}")

# ── Step 3: Age acceleration residuals ────────────────────────────────────────
# Age acceleration = clock output - expected clock output for this age
# Computed by regressing clock output on chronological age and taking residuals.
# This removes the trivial age correlation and focuses on inter-individual
# variation at the same calendar age.
# DunedinPACE is excluded -- it measures rate, not position, so residuals
# relative to chronological age are not meaningful.

print("\n=== Step 3: Computing age acceleration residuals ===")

accel_df = pd.DataFrame(index=sample_ids)
for clock_name in CLOCKS:
    if CLOCKS[clock_name]["type"] == "T_rate":
        # For rate clocks, center around 1.0 (the expected rate)
        vals = clock_df[clock_name].values
        accel_df[clock_name] = vals - np.nanmean(vals)
        print(f"  {clock_name} (rate): centered around mean")
        continue

    vals = clock_df[clock_name].values
    valid = ~np.isnan(vals) & ~np.isnan(ages)
    if valid.sum() < 10:
        accel_df[clock_name] = np.nan
        continue

    # Regress clock output on chronological age
    slope, intercept, _, _, _ = linregress(ages[valid], vals[valid])
    expected = slope * ages + intercept
    residuals = vals - expected
    accel_df[clock_name] = residuals
    print(f"  {clock_name}: slope={slope:.3f}, mean accel={residuals[valid].mean():.3f}")

accel_df.index.name = "sample_id"
accel_df.to_parquet(AGE_ACCEL)
print(f"  Saved: {AGE_ACCEL}")

# ── Step 4: Rank consistency matrix kappa(k, k') ──────────────────────────────
# kappa(k, k') = P[sign(k(A) - k(B)) == sign(k'(A) - k'(B))]
# for pairs (A, B) drawn from the population.
# This is (Kendall's tau + 1) / 2, rescaled to [0, 1].
# We compute within each dataset separately, then combined.
#
# Only compare clocks of the same type (per D0).
# We report all pairs but flag cross-type comparisons.

print("\n=== Step 4: Computing rank consistency matrix ===")

clock_names = list(CLOCKS.keys())
n_clocks = len(clock_names)

def kendall_kappa(x, y):
    """
    Kendall's tau rescaled to [0,1].
    kappa = (tau + 1) / 2
    1.0 = perfect agreement, 0.5 = random, 0.0 = perfect disagreement.
    """
    valid = ~np.isnan(x) & ~np.isnan(y)
    if valid.sum() < 10:
        return np.nan
    tau, _ = kendalltau(x[valid], y[valid])
    return (tau + 1) / 2

# Compute for all dataset subsets
subsets = {
    "combined": np.ones(len(sample_ids), dtype=bool),
    "GSE40279": datasets == "GSE40279",
    "GSE87571": datasets == "GSE87571",
}

kappa_results = {}
for subset_name, mask in subsets.items():
    kappa_mat = np.full((n_clocks, n_clocks), np.nan)
    for i, k1 in enumerate(clock_names):
        for j, k2 in enumerate(clock_names):
            x = accel_df.loc[mask, k1].values if k1 in accel_df.columns else np.full(mask.sum(), np.nan)
            y = accel_df.loc[mask, k2].values if k2 in accel_df.columns else np.full(mask.sum(), np.nan)
            kappa_mat[i, j] = kendall_kappa(x, y)
    kappa_results[subset_name] = pd.DataFrame(
        kappa_mat, index=clock_names, columns=clock_names
    )

print("\n  Rank consistency matrix (combined dataset):")
print(kappa_results["combined"].round(3).to_string())

# Save combined kappa matrix
kappa_results["combined"].to_parquet(RANK_CONSISTENCY)
print(f"\n  Saved: {RANK_CONSISTENCY}")

# ── Step 5: Decompose rank reversals ─────────────────────────────────────────
# For each pair of same-type clocks, identify rank reversals and decompose
# into tau-driven vs r-driven.
#
# A reversal between clocks k and k' for pair (A, B) is:
#   sign(k(A) - k(B)) != sign(k'(A) - k'(B))
#
# It is tau-driven if sign(tau(A) - tau(B)) agrees with one clock but not
# the other -- meaning the clocks disagree on the canonical age ordering.
#
# It is r-driven if both clocks agree on the tau ordering but disagree
# on the total output because f(r) shifts the ranking.

print("\n=== Step 5: Decomposing rank reversals ===")

tau_vals = tau_df.loc[sample_ids, "tau"].values
residual_norms = tau_df.loc[sample_ids, "residual_norm"].values

# Same-type clock pairs to analyse
same_type_pairs = [
    ("Horvath", "Hannum",   "T_tau"),
    ("Horvath", "PhenoAge", "cross"),
    ("Horvath",  "GrimAge",  "cross"),
    ("PhenoAge", "GrimAge",  "T_delta"),
]

reversal_records = []
for k1, k2, pair_type in same_type_pairs:
    if k1 not in accel_df.columns or k2 not in accel_df.columns:
        continue

    a1 = accel_df[k1].values
    a2 = accel_df[k2].values
    valid = ~np.isnan(a1) & ~np.isnan(a2) & ~np.isnan(tau_vals)
    n_valid = valid.sum()

    if n_valid < 10:
        continue

    a1v, a2v = a1[valid], a2[valid]
    tv = tau_vals[valid]
    rv = residual_norms[valid]

    # For each pair of individuals, check for rank reversals
    n_pairs = 0
    n_reversals = 0
    n_tau_driven = 0
    n_r_driven = 0

    # Sample random pairs for efficiency (n_pairs ~ 10000)
    rng = np.random.default_rng(RANDOM_SEED)
    idx = rng.choice(n_valid, size=(min(10000, n_valid*(n_valid-1)//2), 2),
                     replace=True)
    # Remove self-pairs
    idx = idx[idx[:, 0] != idx[:, 1]]

    for i, j in idx:
        n_pairs += 1
        sign_k1  = np.sign(a1v[i] - a1v[j])
        sign_k2  = np.sign(a2v[i] - a2v[j])
        sign_tau = np.sign(tv[i]  - tv[j])

        if sign_k1 == 0 or sign_k2 == 0:
            continue
        if sign_k1 != sign_k2:
            n_reversals += 1
            # Tau-driven: tau disagrees with at least one clock
            if sign_tau != sign_k1 or sign_tau != sign_k2:
                n_tau_driven += 1
            else:
                # Both clocks agree with tau ordering but still disagree
                # with each other -- r-driven
                n_r_driven += 1

    reversal_rate = n_reversals / n_pairs if n_pairs > 0 else np.nan
    tau_frac = n_tau_driven / n_reversals if n_reversals > 0 else np.nan
    r_frac   = n_r_driven   / n_reversals if n_reversals > 0 else np.nan

    kappa_val = kappa_results["combined"].loc[k1, k2]

    print(f"\n  {k1} vs {k2} (type: {pair_type})")
    print(f"    kappa:          {kappa_val:.3f}")
    print(f"    reversal rate:  {reversal_rate:.3f}")
    print(f"    tau-driven:     {tau_frac:.3f} of reversals")
    print(f"    r-driven:       {r_frac:.3f} of reversals")

    reversal_records.append({
        "clock1": k1, "clock2": k2, "pair_type": pair_type,
        "kappa": kappa_val,
        "reversal_rate": reversal_rate,
        "tau_driven_frac": tau_frac,
        "r_driven_frac": r_frac,
        "n_pairs": n_pairs,
        "n_reversals": n_reversals,
    })

reversals_df = pd.DataFrame(reversal_records)
reversals_df.to_parquet(RANK_REVERSALS)
print(f"\n  Saved: {RANK_REVERSALS}")

# ── Step 6: High-instability individuals ──────────────────────────────────────
# Instability score per individual = std of their rank percentile across clocks.
# High score = their biological age rank varies widely depending on which
# clock you use.

print("\n=== Step 6: Identifying high-instability individuals ===")

# Compute rank percentile for each clock
rank_pcts = pd.DataFrame(index=sample_ids)
for col in clock_names:
    vals = accel_df[col].values
    valid_mask = ~np.isnan(vals)
    ranks = np.full(len(vals), np.nan)
    ranks[valid_mask] = pd.Series(vals[valid_mask]).rank(pct=True).values
    rank_pcts[col] = ranks

# Instability = std of rank percentiles across clocks (excluding DunedinPACE
# since it's a different type)
position_clocks = [k for k, v in CLOCKS.items() if v["type"] != "T_rate"]
instability = rank_pcts[position_clocks].std(axis=1)
instability.name = "instability_score"

instability_df = pd.DataFrame({
    "instability_score": instability,
    "tau": tau_df.loc[sample_ids, "tau"],
    "residual_norm": tau_df.loc[sample_ids, "residual_norm"],
    "age": ages,
    "dataset": datasets,
})
instability_df.index.name = "sample_id"

# Add clock ranks
for col in position_clocks:
    instability_df[f"rank_{col}"] = rank_pcts[col]

instability_df.to_parquet(INSTABILITY)
print(f"  Saved: {INSTABILITY}")

# Top 10 most unstable individuals
top10 = instability_df.nlargest(10, "instability_score")
print(f"\n  Top 10 most unstable individuals:")
print(top10[["instability_score", "age", "tau", "residual_norm", "dataset"]].round(3).to_string())

# ── Step 7: Directional stability (D2) ────────────────────────────────────────
# For each clock, compute the directional sensitivity ratio:
#   rho(k) = sensitivity in off-manifold directions /
#             sensitivity in on-manifold direction
#
# Since all clocks are linear (weighted sums), the sensitivity in a
# direction d is |w^T d| where w is the clock's coefficient vector.
# The on-manifold direction at each sample is the tangent to gamma at
# the projected point.
#
# We approximate the tangent direction using the principal curve output:
# for each sample, the on-manifold direction ~ the local gradient of tau
# in methylation space, which we estimate from the top-CpG weights.

print("\n=== Step 7: Directional stability ===")

# Load top CpGs and weights
with open(DATA_DIR / "top_cpgs.txt") as f:
    top_cpgs = [line.strip() for line in f if line.strip()]

weights_df = pd.read_parquet(WEIGHTS)
w_top = weights_df.loc[top_cpgs, "weight"].values

# Load residuals (samples x top_cpgs)
residuals_df = pd.read_parquet(RESIDUALS)
residuals_df = residuals_df.loc[sample_ids]

# Estimate on-manifold direction from the correlation between each top CpG
# and tau: the manifold direction is approximately the vector of correlations
# between CpG values and tau (the direction of maximum age-related change).
# Load top CpGs directly from HDF5 files (not from clock-restricted beta_all)
print("  Loading top CpGs from HDF5 files for D2 analysis...")

# Load each dataset, restrict to top CpGs, rename GSE87571 columns
h5_40  = pd.read_hdf(str(GSE40279_BETA), key="beta")
h5_87  = pd.read_hdf(str(GSE87571_BETA), key="beta").rename(columns=col_to_id)

top_cpgs_available = [c for c in top_cpgs
                      if c in h5_40.index and c in h5_87.index]

# Restrict to top CpGs and combine: result is CpGs x samples
beta_top_combined = pd.concat(
    [h5_40.loc[top_cpgs_available],
     h5_87.loc[top_cpgs_available]], axis=1
)  # shape: len(top_cpgs_available) x n_samples

# Transpose to samples x CpGs, deduplicate sample index
beta_top = beta_top_combined.T   # samples x CpGs
beta_top = beta_top[~beta_top.index.duplicated(keep="first")]

# Align strictly to tau_df index (1385 samples with valid age)
common_s = [s for s in tau_df.index if s in beta_top.index]
beta_top_vals = beta_top.loc[common_s, top_cpgs_available].values
tau_vals_top  = tau_df.loc[common_s, "tau"].values
top_cpgs      = top_cpgs_available

print(f"  top_cpgs_available: {len(top_cpgs_available)}")
print(f"  common_s: {len(common_s)}")
print(f"  beta_top_vals: {beta_top_vals.shape}, tau_vals_top: {tau_vals_top.shape}")
assert beta_top_vals.shape[0] == tau_vals_top.shape[0],     f"Shape mismatch after dedup: {beta_top_vals.shape[0]} vs {tau_vals_top.shape[0]}"
manifold_direction = np.array([
    pearsonr(beta_top_vals[:, j], tau_vals_top)[0]
    for j in range(len(top_cpgs))
])
# Normalize
manifold_direction /= (np.linalg.norm(manifold_direction) + 1e-10)

# For each clock, get its coefficient vector over the top CpGs
# We estimate this by regressing clock output on each top CpG's beta values
stability_records = []
for clock_name in CLOCKS:
    # Align clock values to common_s
    vals = clock_df.loc[common_s, clock_name].values
    valid = ~np.isnan(vals)
    if valid.sum() < 100:
        continue

    from sklearn.linear_model import Ridge
    ridge = Ridge(alpha=1.0)
    ridge.fit(beta_top_vals[valid], vals[valid])
    coef = ridge.coef_   # (n_top_cpgs,)

    # Sensitivity in manifold direction (on-manifold)
    sens_on  = abs(coef @ manifold_direction)

    # Sensitivity in off-manifold directions (complement of manifold direction)
    coef_off = coef - (coef @ manifold_direction) * manifold_direction
    sens_off = np.linalg.norm(coef_off)

    rho = sens_off / (sens_on + 1e-10)

    print(f"  {clock_name}: rho={rho:.3f} "
          f"(on={sens_on:.4f}, off={sens_off:.4f})")

    stability_records.append({
        "clock": clock_name,
        "type": CLOCKS[clock_name]["type"],
        "rho": rho,
        "sens_on": sens_on,
        "sens_off": sens_off,
    })

stability_df = pd.DataFrame(stability_records)

# ── Step 8: Figures ───────────────────────────────────────────────────────────
print("\n=== Step 8: Generating figures ===")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Figure 1: Rank consistency heatmap
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, (subset_name, kappa_mat) in zip(axes, kappa_results.items()):
    # Annotate cross-type pairs
    annot = kappa_mat.round(3).astype(str)
    mask = np.zeros_like(kappa_mat.values, dtype=bool)

    im = ax.imshow(kappa_mat.values, vmin=0.5, vmax=1.0, cmap="RdYlGn")
    ax.set_xticks(range(n_clocks))
    ax.set_yticks(range(n_clocks))
    ax.set_xticklabels(clock_names, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(clock_names, fontsize=9)
    for i in range(n_clocks):
        for j in range(n_clocks):
            val = kappa_mat.values[i, j]
            if not np.isnan(val):
                # Mark cross-type pairs with asterisk
                t1 = CLOCKS[clock_names[i]]["type"]
                t2 = CLOCKS[clock_names[j]]["type"]
                label = f"{val:.3f}"
                if t1 != t2:
                    label += "*"
                ax.text(j, i, label, ha="center", va="center",
                       fontsize=7, color="black")
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title(f"Rank consistency kappa\n({subset_name})\n* = cross-type pair")

plt.tight_layout()
fig_path = FIGURES_DIR / "03_rank_consistency_matrix.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
print(f"  Saved: {fig_path}")
plt.close()

# Figure 2: Tau-r decomposition of rank reversals
if len(reversal_records) > 0:
    fig, ax = plt.subplots(figsize=(8, 5))
    pairs_labels = [f"{r['clock1']}\nvs\n{r['clock2']}" for r in reversal_records]
    tau_fracs = [r["tau_driven_frac"] for r in reversal_records]
    r_fracs   = [r["r_driven_frac"]   for r in reversal_records]
    x = np.arange(len(pairs_labels))
    width = 0.35
    ax.bar(x - width/2, tau_fracs, width, label="tau-driven", color="steelblue")
    ax.bar(x + width/2, r_fracs,   width, label="r-driven",   color="coral")
    ax.set_xticks(x)
    ax.set_xticklabels(pairs_labels, fontsize=8)
    ax.set_ylabel("Fraction of rank reversals")
    ax.set_title("Decomposition of rank reversals\ninto tau-driven vs r-driven")
    ax.legend()
    ax.set_ylim(0, 1)
    plt.tight_layout()
    fig_path = FIGURES_DIR / "03_tau_r_decomposition.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {fig_path}")
    plt.close()

# Figure 3: High-instability individuals
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.hist(instability_df["instability_score"].dropna(), bins=50,
        color="steelblue", alpha=0.7, edgecolor="white")
ax.set_xlabel("Instability score (std of rank percentiles)")
ax.set_ylabel("Count")
ax.set_title("Distribution of clock rank instability\nacross individuals")

ax = axes[1]
sc = ax.scatter(
    instability_df["age"],
    instability_df["instability_score"],
    c=instability_df["residual_norm"],
    cmap="YlOrRd", alpha=0.4, s=15,
)
plt.colorbar(sc, ax=ax, label="Residual norm ||r||_*")
ax.set_xlabel("Chronological age")
ax.set_ylabel("Instability score")
ax.set_title("Instability vs age\n(coloured by off-manifold distance)")

plt.tight_layout()
fig_path = FIGURES_DIR / "03_high_instability_individuals.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
print(f"  Saved: {fig_path}")
plt.close()

# Figure 4: Directional stability
if len(stability_records) > 0:
    fig, ax = plt.subplots(figsize=(7, 4))
    colors_by_type = {"T_tau": "steelblue", "T_delta": "coral", "T_rate": "seagreen"}
    clk_names = [r["clock"] for r in stability_records]
    rhos = [r["rho"] for r in stability_records]
    colors = [colors_by_type[r["type"]] for r in stability_records]
    bars = ax.bar(clk_names, rhos, color=colors, alpha=0.8, edgecolor="white")
    ax.set_ylabel("Directional sensitivity ratio rho\n(off-manifold / on-manifold)")
    ax.set_title("D2: Directional stability per clock\n(lower = more stable)")
    ax.axhline(1.0, color="black", linestyle="--", lw=1, alpha=0.5)
    legend_patches = [
        mpatches.Patch(color="steelblue", label="T_tau (position)"),
        mpatches.Patch(color="coral",     label="T_delta (deviation)"),
        mpatches.Patch(color="seagreen",  label="T_rate (rate)"),
    ]
    ax.legend(handles=legend_patches)
    plt.tight_layout()
    fig_path = FIGURES_DIR / "03_directional_stability.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {fig_path}")
    plt.close()

print("\n=== Script 3 complete ===")
print(f"\nSummary:")
print(f"  Clocks computed: {list(clock_df.columns)}")
print(f"  Samples: {len(sample_ids)}")
print(f"\nRank consistency matrix (combined):")
print(kappa_results["combined"].round(3).to_string())
print(f"\nRank reversal decomposition:")
print(reversals_df[["clock1","clock2","kappa","reversal_rate",
                     "tau_driven_frac","r_driven_frac"]].round(3).to_string())

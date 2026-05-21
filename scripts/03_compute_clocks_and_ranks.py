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
    "GrimAge":    {"model": "GrimAgeV2",  "type": "T_delta"},
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

# Deduplicate beta_all columns before anything else
beta_all = beta_all.loc[:, ~beta_all.columns.duplicated(keep="first")]
print(f"  After dedup: {beta_all.shape} (CpGs x samples)")

# Align metadata to beta columns
sample_ids   = beta_all.columns.tolist()
meta_aligned = metadata.set_index("sample_id").loc[sample_ids]
ages         = meta_aligned["age"].values
datasets     = meta_aligned["dataset"].values

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
        if isinstance(result, pd.DataFrame):
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
# Deduplicate sample index (can occur from GSE87571 column alignment)
clock_df = clock_df[~clock_df.index.duplicated(keep="first")]
clock_df.to_parquet(CLOCK_OUTPUTS)
# Realign sample_ids and ages to deduped clock_df
sample_ids = clock_df.index.tolist()
ages       = meta_aligned.loc[sample_ids, "age"].values
datasets   = meta_aligned.loc[sample_ids, "dataset"].values

# Calibrate GrimAge to chronological age units (years)
# biolearn returns raw composite score; rescale linearly to years
# using regression on chronological age (standard practice)
from scipy.stats import linregress as _lr
_grim = clock_df["GrimAge"].values
_valid = ~np.isnan(_grim)
_slope, _intercept, _, _, _ = _lr(ages[_valid], _grim[_valid])
# Invert: predicted_years = (raw - intercept) / slope
clock_df["GrimAge"] = (clock_df["GrimAge"] - _intercept) / _slope
print(f"  GrimAge recalibrated: mean={clock_df['GrimAge'].mean():.1f}, "
      f"range=[{clock_df['GrimAge'].min():.1f}, {clock_df['GrimAge'].max():.1f}]")
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

    # Compute tau differences and clock differences for all pairs at once
    i_idx = idx[:, 0]
    j_idx = idx[:, 1]

    diff_k1  = a1v[i_idx] - a1v[j_idx]
    diff_k2  = a2v[i_idx] - a2v[j_idx]
    diff_tau = tv[i_idx]  - tv[j_idx]
    diff_r   = rv[i_idx]  - rv[j_idx]

    sign_k1  = np.sign(diff_k1)
    sign_k2  = np.sign(diff_k2)

    # Only consider pairs where both clocks give non-tied rankings
    valid_pairs = (sign_k1 != 0) & (sign_k2 != 0)
    n_pairs = valid_pairs.sum()

    # Reversals: clocks disagree on ranking
    reversal_mask = valid_pairs & (sign_k1 != sign_k2)
    n_reversals = reversal_mask.sum()

    if n_reversals > 0:
        rev_diff_tau = np.abs(diff_tau[reversal_mask])
        rev_diff_r   = np.abs(diff_r[reversal_mask])

        # Decomposition:
        # A reversal is TAU-DOMINATED if the tau difference between the two
        # individuals is SMALL relative to the median -- meaning they are
        # close on the canonical trajectory, so the residual functions f(r)
        # are determining the ranking.
        # A reversal is RESIDUAL-DOMINATED if tau difference is LARGE --
        # meaning the clocks disagree despite the individuals being clearly
        # separated on the trajectory, suggesting they weight tau differently.
        tau_median = np.median(np.abs(diff_tau[valid_pairs]))
        r_median   = np.median(np.abs(diff_r[valid_pairs]))

        # Tau-dominated: small tau difference (below median) -- residual
        # variation is driving the disagreement
        tau_dominated = (rev_diff_tau < tau_median)
        n_tau_driven  = tau_dominated.sum()
        n_r_driven    = (~tau_dominated).sum()

        # Also compute correlation between |diff_tau| and reversal probability
        # as a continuous measure
        for i, j in zip([], []):  # placeholder to keep loop structure
            pass

    reversal_rate = int(n_reversals) / int(n_pairs) if n_pairs > 0 else np.nan
    tau_frac = int(n_tau_driven) / int(n_reversals) if n_reversals > 0 else np.nan
    r_frac   = int(n_r_driven)   / int(n_reversals) if n_reversals > 0 else np.nan

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
    # Align clock values to common_s using reindex (safe with any index)
    vals = clock_df[clock_name].reindex(common_s).values
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
    ax.set_xticklabels(pairs_labels, fontsize=9, linespacing=0.85)
    ax.set_ylabel("Fraction of rank reversals")
    ax.set_title("Decomposition of rank reversals\ninto tau-driven vs r-driven")
    ax.legend()
    ax.set_ylim(0, 1.05)
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
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
"""
Additional analyses to append to 03_compute_clocks_and_ranks.py
Steps 9-12:
    9.  D1 monotonicity test
    10. Coherence test (how well does tau/r explain each clock?)
    11. Cell type composition stability (D3 empirical)
    12. Age acceleration correlation heatmap
"""

# ── Step 9: D1 — Monotonicity test ───────────────────────────────────────────
# D1 requires that E[k(s(t))] is non-decreasing in t for position/deviation
# clocks, and E[k(s(t))] > 0 for rate clocks.
# We test this by binning samples into age deciles and checking whether
# mean clock output increases monotonically across bins.

print("\n=== Step 9: D1 Monotonicity test ===")

from scipy.stats import spearmanr

# Bin by chronological age into deciles
n_bins = 10
age_bins = pd.qcut(ages, q=n_bins, labels=False)

monotonicity_records = []
for clock_name in CLOCKS:
    clock_type = CLOCKS[clock_name]["type"]
    vals = clock_df[clock_name].values

    # Mean clock output per age bin
    bin_means = []
    bin_ages  = []
    for b in range(n_bins):
        mask_b = (age_bins == b) & ~np.isnan(vals)
        if mask_b.sum() > 0:
            bin_means.append(vals[mask_b].mean())
            bin_ages.append(ages[mask_b].mean())

    bin_means = np.array(bin_means)
    bin_ages  = np.array(bin_ages)

    # Spearman correlation between bin rank and mean clock output
    rho, pval = spearmanr(bin_ages, bin_means)

    # Count monotone violations: bins where mean decreases
    diffs = np.diff(bin_means)
    n_violations = (diffs < 0).sum()
    n_possible   = len(diffs)

    if clock_type == "T_rate":
        # For rate clocks check mean > 0 (already known for DunedinPACE)
        mean_rate = vals[~np.isnan(vals)].mean()
        passes = mean_rate > 0
        print(f"  {clock_name} (T_rate): mean={mean_rate:.3f} "
              f"({'PASS' if passes else 'FAIL'}: expected > 0)")
    else:
        passes = n_violations == 0
        print(f"  {clock_name} ({clock_type}): rho={rho:.3f}, "
              f"violations={n_violations}/{n_possible} "
              f"({'PASS' if passes else f'FAIL: {n_violations} non-monotone bins'})")

    monotonicity_records.append({
        "clock": clock_name,
        "type":  clock_type,
        "spearman_rho": rho,
        "spearman_p": pval,
        "n_violations": n_violations,
        "n_bins": n_possible,
        "passes_D1": passes,
        "bin_means": bin_means.tolist(),
        "bin_ages":  bin_ages.tolist(),
    })

monotonicity_df = pd.DataFrame([
    {k: v for k, v in r.items() if k not in ("bin_means", "bin_ages")}
    for r in monotonicity_records
])

# Figure: mean clock output per age bin
fig, axes = plt.subplots(1, len(CLOCKS), figsize=(18, 4))
colors_by_type = {"T_tau": "steelblue", "T_delta": "coral", "T_rate": "seagreen"}
for ax, rec in zip(axes, monotonicity_records):
    color = colors_by_type[rec["type"]]
    ax.plot(rec["bin_ages"], rec["bin_means"], "o-", color=color, lw=2, ms=6)
    ax.set_xlabel("Mean chronological age in bin")
    ax.set_ylabel("Mean clock output")
    ax.set_title(f"{rec['clock']}\nrho={rec['spearman_rho']:.3f}, "
                 f"violations={rec['n_violations']}")
    ax.grid(True, alpha=0.3)

plt.suptitle("D1: Monotonicity of clock outputs across age bins", y=1.02)
plt.tight_layout()
fig_path = FIGURES_DIR / "03_d1_monotonicity.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
print(f"  Saved: {fig_path}")
plt.close()

# ── Step 10: Coherence test ───────────────────────────────────────────────────
# A coherent clock should be well-explained by (tau, ||r||_*).
# We regress each clock's age acceleration on tau and residual_norm,
# and report R^2. A high R^2 means the clock is coherent with our framework.
# A low R^2 means the clock is measuring something not captured by tau or r.

print("\n=== Step 10: Coherence test ===")

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

# Align tau and residuals to sample_ids
tau_aligned      = tau_df.loc[sample_ids, "tau"].values
resid_aligned    = tau_df.loc[sample_ids, "residual_norm"].values

# Feature matrix: [tau, residual_norm, tau^2] (allow mild nonlinearity)
X_coh = np.column_stack([
    tau_aligned,
    resid_aligned,
    tau_aligned ** 2,
    resid_aligned ** 2,
    tau_aligned * resid_aligned,
])

scaler = StandardScaler()
X_coh_scaled = scaler.fit_transform(X_coh)

coherence_records = []
for clock_name in CLOCKS:
    clock_type = CLOCKS[clock_name]["type"]
    vals = accel_df[clock_name].values
    valid = ~np.isnan(vals) & ~np.isnan(tau_aligned) & ~np.isnan(resid_aligned)

    if valid.sum() < 50:
        continue

    lr = LinearRegression()
    lr.fit(X_coh_scaled[valid], vals[valid])
    preds = lr.predict(X_coh_scaled[valid])
    r2 = r2_score(vals[valid], preds)

    # Also compute partial R^2 for tau alone and r alone
    lr_tau = LinearRegression()
    lr_tau.fit(tau_aligned[valid].reshape(-1,1), vals[valid])
    r2_tau = r2_score(vals[valid], lr_tau.predict(tau_aligned[valid].reshape(-1,1)))

    lr_r = LinearRegression()
    lr_r.fit(resid_aligned[valid].reshape(-1,1), vals[valid])
    r2_r = r2_score(vals[valid], lr_r.predict(resid_aligned[valid].reshape(-1,1)))

    print(f"  {clock_name} ({clock_type}): R2(tau+r)={r2:.3f}, "
          f"R2(tau only)={r2_tau:.3f}, R2(r only)={r2_r:.3f}")

    coherence_records.append({
        "clock": clock_name,
        "type":  clock_type,
        "R2_tau_and_r": r2,
        "R2_tau_only":  r2_tau,
        "R2_r_only":    r2_r,
        "R2_unexplained": 1 - r2,
    })

coherence_df = pd.DataFrame(coherence_records)

# Figure: coherence bar chart
fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(coherence_records))
width = 0.25
colors_t = [colors_by_type[r["type"]] for r in coherence_records]
clk_labels = [r["clock"] for r in coherence_records]

ax.bar(x - width, [r["R2_tau_only"]  for r in coherence_records],
       width, label="R² (tau only)",    color="steelblue",  alpha=0.8)
ax.bar(x,         [r["R2_r_only"]    for r in coherence_records],
       width, label="R² (||r|| only)", color="coral",       alpha=0.8)
ax.bar(x + width, [r["R2_tau_and_r"] for r in coherence_records],
       width, label="R² (tau + ||r||)", color="seagreen",   alpha=0.8)

ax.set_xticks(x)
ax.set_xticklabels(clk_labels)
ax.set_ylabel("R² (fraction of age acceleration variance explained)")
ax.set_title("Coherence test: how well does (tau, ||r||) explain each clock?\n"
             "Higher = more coherent with canonical trajectory framework")
ax.legend()
ax.set_ylim(0, 1)
ax.yaxis.grid(True, alpha=0.3)
ax.set_axisbelow(True)
plt.tight_layout()
fig_path = FIGURES_DIR / "03_coherence_test.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
print(f"  Saved: {fig_path}")
plt.close()

# ── Step 11: Cell type stability (D3 empirical) ───────────────────────────────
# D3 requires clocks to be stable under cell type composition shifts.
# We estimate blood cell type proportions using biolearn's DeconvoluteBlood450K,
# then measure how much each clock's age acceleration correlates with cell
# composition after controlling for chronological age.
# A high correlation = the clock is confounded by cell type = D3 violated.

print("\n=== Step 11: Cell type stability (D3) ===")

try:
    print("  Running blood cell type deconvolution...")
    deconv_model  = gallery.get("DeconvoluteBlood450K")
    deconv_result = deconv_model.predict(geo_data)

    print(f"  Deconvolution output shape: {deconv_result.shape}")
    print(f"  Cell types: {deconv_result.columns.tolist()}")

    # Align to sample_ids
    deconv_aligned = deconv_result.reindex(sample_ids)
    cell_types = deconv_result.columns.tolist()

    # For each clock, partial correlation of age acceleration with each
    # cell type proportion, after controlling for chronological age
    from scipy.stats import pearsonr as _pr

    d3_records = []
    for clock_name in CLOCKS:
        if CLOCKS[clock_name]["type"] == "T_rate":
            continue
        accel_vals = accel_df[clock_name].values
        valid = ~np.isnan(accel_vals)

        cell_cors = {}
        for ct in cell_types:
            ct_vals = deconv_aligned[ct].values
            ct_valid = valid & ~np.isnan(ct_vals)
            if ct_valid.sum() < 50:
                continue
            r, p = _pr(accel_vals[ct_valid], ct_vals[ct_valid])
            cell_cors[ct] = r

        # Summary: max absolute correlation with any cell type
        max_cor = max(abs(v) for v in cell_cors.values()) if cell_cors else np.nan
        print(f"  {clock_name}: max |r| with cell types = {max_cor:.3f} "
              f"({'PASS D3' if max_cor < 0.2 else 'FAIL D3: cell type confound'})")

        d3_records.append({
            "clock": clock_name,
            "type":  CLOCKS[clock_name]["type"],
            "max_cell_cor": max_cor,
            "cell_cors": cell_cors,
            "passes_D3": max_cor < 0.2,
        })

    # Figure: cell type correlation heatmap
    cell_cor_mat = pd.DataFrame(
        {r["clock"]: r["cell_cors"] for r in d3_records}
    ).T  # clocks x cell_types

    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(cell_cor_mat.values, cmap="RdBu_r", vmin=-0.4, vmax=0.4,
                   aspect="auto")
    ax.set_xticks(range(len(cell_cor_mat.columns)))
    ax.set_xticklabels(cell_cor_mat.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(cell_cor_mat.index)))
    ax.set_yticklabels(cell_cor_mat.index, fontsize=10)
    for i in range(len(cell_cor_mat.index)):
        for j in range(len(cell_cor_mat.columns)):
            val = cell_cor_mat.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                       fontsize=8, color="black")
    plt.colorbar(im, ax=ax, label="Pearson r with age acceleration")
    ax.set_title("D3: Cell type composition confounding\n"
                 "(correlation of age acceleration with cell type proportions)")
    plt.tight_layout()
    fig_path = FIGURES_DIR / "03_d3_cell_type.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {fig_path}")
    plt.close()

except Exception as e:
    print(f"  Cell type deconvolution failed: {e}")
    print("  Skipping D3 cell type analysis.")
    d3_records = []

# ── Step 12: Age acceleration correlation heatmap ────────────────────────────
# Shows pairwise Pearson correlations between age acceleration residuals.
# Complements the rank consistency (Kendall kappa) matrix.
# High Pearson r but low kappa = clocks agree on average but disagree on
# individual rankings (non-linear relationship).

print("\n=== Step 12: Age acceleration correlation heatmap ===")

from scipy.stats import pearsonr as _pr2

corr_mat = np.full((n_clocks, n_clocks), np.nan)
for i, k1 in enumerate(clock_names):
    for j, k2 in enumerate(clock_names):
        v1 = accel_df[k1].values
        v2 = accel_df[k2].values
        valid = ~np.isnan(v1) & ~np.isnan(v2)
        if valid.sum() < 10:
            continue
        r, _ = _pr2(v1[valid], v2[valid])
        corr_mat[i, j] = r

corr_df = pd.DataFrame(corr_mat, index=clock_names, columns=clock_names)
print("\n  Pearson correlation of age accelerations:")
print(corr_df.round(3).to_string())

# Figure: side-by-side Pearson r and Kendall kappa
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (mat, title) in zip(axes, [
    (corr_mat,                        "Pearson r\n(age acceleration)"),
    (kappa_results["combined"].values, "Kendall kappa\n(rank consistency)"),
]):
    im = ax.imshow(mat, vmin=0.5, vmax=1.0, cmap="RdYlGn")
    ax.set_xticks(range(n_clocks))
    ax.set_yticks(range(n_clocks))
    ax.set_xticklabels(clock_names, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(clock_names, fontsize=9)
    for i in range(n_clocks):
        for j in range(n_clocks):
            val = mat[i, j]
            if not np.isnan(val):
                t1 = CLOCKS[clock_names[i]]["type"]
                t2 = CLOCKS[clock_names[j]]["type"]
                label = f"{val:.3f}" + ("*" if t1 != t2 else "")
                ax.text(j, i, label, ha="center", va="center",
                       fontsize=7.5, color="black")
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title(title + "\n(* = cross-type pair)")

plt.suptitle("Agreement between clocks: Pearson r vs Kendall kappa",
             fontsize=12, y=1.02)
plt.tight_layout()
fig_path = FIGURES_DIR / "03_pearson_vs_kappa.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
print(f"  Saved: {fig_path}")
plt.close()

# ── Save all additional analysis outputs ──────────────────────────────────────
print("\n=== Saving additional analysis outputs ===")

monotonicity_df.to_parquet(DATA_DIR / "monotonicity.parquet")
coherence_df.to_parquet(DATA_DIR / "coherence.parquet")
if d3_records:
    pd.DataFrame([{k: v for k, v in r.items() if k != "cell_cors"}
                  for r in d3_records]).to_parquet(DATA_DIR / "d3_cell_type.parquet")
corr_df.to_parquet(DATA_DIR / "age_accel_correlations.parquet")

print("  All outputs saved.")
print("\n=== All analyses complete ===")
print(f"\nD1 Monotonicity summary:")
print(monotonicity_df[["clock","type","spearman_rho","n_violations","passes_D1"]].to_string())
print(f"\nCoherence summary:")
print(coherence_df[["clock","type","R2_tau_only","R2_r_only","R2_tau_and_r"]].round(3).to_string())

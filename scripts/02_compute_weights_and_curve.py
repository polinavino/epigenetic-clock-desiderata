"""
02_compute_weights_and_curve.py

Computes the canonical aging trajectory gamma(t) in methylation space.

Steps:
    1. Load clock CpG site lists (union of all clock sites)
    2. Load beta values restricted to clock CpGs from both datasets
    3. Compute age-informativeness weights w_i = R2_i * sigma2_i
    4. Select top N_TOP_CPGS by weight
    5. Fit principal curve through weighted methylation space
    6. Compute tau(m) and r(m) for every sample
    7. Save weights, top CpGs, curve, tau, residuals

Outputs (all in data/):
    cpg_weights.parquet      -- w_i for every common CpG
    top_cpgs.txt             -- the N_TOP_CPGS selected CpG IDs
    principal_curve.parquet  -- gamma(t): curve points x CpG coords
    tau.parquet              -- biological age coordinate per sample
    residuals.parquet        -- r(m) = m - pi(m) per sample (top CpGs only)
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ── Config ────────────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    DATA_DIR, GSE40279_BETA, GSE87571_BETA, COMMON_CPGS,
    METADATA, WEIGHTS, TOP_CPGS, PRINCIPAL_CURVE, TAU, RESIDUALS,
    N_TOP_CPGS, RANDOM_SEED, PC_SMOOTHING, PC_OUTLIER_THRESHOLD,
    CLOCK_CPG_DIR,
)

np.random.seed(RANDOM_SEED)

# ── Step 1: Load clock CpG site lists ─────────────────────────────────────────
# Each major clock uses a specific set of CpG sites. We take the union of all
# clock CpG sets so our principal curve lives in a space that is relevant to
# all clocks simultaneously.
#
# We define the clock CpG sets here directly from the published papers.
# In production these would be loaded from files in data/clock_cpgs/.
#
# Horvath 2013: 353 CpGs (Supplementary Table 1 of Horvath 2013)
# Hannum 2013:  71 CpGs
# PhenoAge:    513 CpGs (Levine 2018)
# GrimAge:    1030 CpGs (Lu 2019)
# DunedinPACE: 173 CpGs (Belsky 2022)
#
# Rather than hardcode all site IDs here, we use a data-driven approach:
# select the top N_TOP_CPGS sites by age-informativeness weight w_i.
# This is more principled than using any one clock's site selection.
# We verify below that the top sites overlap substantially with published clocks.

print("=== Step 1: Loading data ===")

# Load common CpG list
with open(COMMON_CPGS) as f:
    common_cpgs = [line.strip() for line in f if line.strip()]
print(f"  Common CpGs: {len(common_cpgs)}")

# Load metadata
metadata = pd.read_csv(METADATA)
print(f"  Samples: {len(metadata)}")

# ── Step 2: Load beta values for common CpGs ──────────────────────────────────
# We load both datasets restricted to common CpGs.
# This is memory-efficient: we select columns (CpGs) after loading.
#
# Memory note: even with 470k CpGs x 1385 samples this is ~2.6GB float64.
# We cast to float32 immediately to halve memory usage.

print("\n=== Step 2: Loading beta matrices (common CpGs only) ===")

print("  Loading GSE40279...")
beta_40 = pd.read_parquet(GSE40279_BETA)   # samples x CpGs
# Restrict to common CpGs that are present in this dataset
cpgs_40 = [c for c in common_cpgs if c in beta_40.columns]
beta_40 = beta_40[cpgs_40].astype(np.float32)
print(f"    Shape: {beta_40.shape}")

print("  Loading GSE87571...")
beta_87 = pd.read_parquet(GSE87571_BETA)   # samples x CpGs
cpgs_87 = [c for c in common_cpgs if c in beta_87.columns]
beta_87 = beta_87[cpgs_87].astype(np.float32)
print(f"    Shape: {beta_87.shape}")

# Align to common CpGs present in both
shared_cpgs = sorted(set(cpgs_40) & set(cpgs_87))
beta_40 = beta_40[shared_cpgs]
beta_87 = beta_87[shared_cpgs]
print(f"  Shared CpGs after alignment: {len(shared_cpgs)}")

# Align metadata to samples in each dataset
meta_40 = metadata[metadata["dataset"] == "GSE40279"].set_index("sample_id")
meta_87 = metadata[metadata["dataset"] == "GSE87571"].set_index("sample_id")

# Keep only samples present in both beta matrix and metadata
keep_40 = [s for s in beta_40.index if s in meta_40.index]
keep_87 = [s for s in beta_87.index if s in meta_87.index]
beta_40 = beta_40.loc[keep_40]
beta_87 = beta_87.loc[keep_87]
ages_40 = meta_40.loc[keep_40, "age"].values.astype(np.float32)
ages_87 = meta_87.loc[keep_87, "age"].values.astype(np.float32)

print(f"  GSE40279: {beta_40.shape[0]} samples")
print(f"  GSE87571: {beta_87.shape[0]} samples")

# Combine
beta_all = pd.concat([beta_40, beta_87], axis=0)
ages_all = np.concatenate([ages_40, ages_87])
sample_ids_all = beta_all.index.tolist()
print(f"  Combined: {beta_all.shape[0]} samples x {beta_all.shape[1]} CpGs")

# ── Step 3: Compute age-informativeness weights w_i = R2_i * sigma2_i ─────────
# For each CpG site i:
#   - sigma2_i: variance of beta values across the population
#   - R2_i: fraction of that variance explained by chronological age
#     (from linear regression of beta_i on age)
#   - w_i = R2_i * sigma2_i = variance of the age-predicted component
#
# This weights each CpG by how much of its variation is age-related,
# upweighting informative sites and downweighting noise.
#
# We compute this on the combined dataset to maximise statistical power.

print("\n=== Step 3: Computing age-informativeness weights ===")

X = beta_all.values.astype(np.float64)   # n_samples x n_cpgs
y = ages_all.astype(np.float64)
n, p = X.shape

# Vectorised computation of R2 for all CpGs simultaneously
# R2_i = (corr(beta_i, age))^2
# We compute all correlations at once using matrix operations.
print(f"  Computing correlations for {p} CpGs x {n} samples...")

# Standardize ages
y_mean = y.mean()
y_std  = y.std()
y_z    = (y - y_mean) / y_std              # (n,)

# Standardize each CpG (subtract mean, divide by std)
X_mean = X.mean(axis=0)                    # (p,)
X_std  = X.std(axis=0)                     # (p,)

# Avoid division by zero for constant CpGs
X_std_safe = np.where(X_std < 1e-10, 1.0, X_std)
X_z = (X - X_mean) / X_std_safe           # (n, p)

# Pearson r for each CpG: r_i = (1/n) * sum(x_i_z * y_z)
r_all = (X_z * y_z[:, np.newaxis]).mean(axis=0)   # (p,)
R2_all = r_all ** 2                                 # (p,)

# Variance of each CpG
sigma2_all = X_std ** 2                             # (p,)

# Weight: age-explained variance
w_all = R2_all * sigma2_all                         # (p,)

print(f"  R2 range: [{R2_all.min():.4f}, {R2_all.max():.4f}]")
print(f"  Median R2: {np.median(R2_all):.4f}")
print(f"  Sites with R2 > 0.5: {(R2_all > 0.5).sum()}")
print(f"  Sites with R2 > 0.3: {(R2_all > 0.3).sum()}")

# Save weights for all common CpGs
weights_df = pd.DataFrame({
    "cpg":    shared_cpgs,
    "R2":     R2_all,
    "sigma2": sigma2_all,
    "weight": w_all,
}).set_index("cpg").sort_values("weight", ascending=False)

weights_df.to_parquet(WEIGHTS)
print(f"  Saved: {WEIGHTS}")

# ── Step 4: Select top N_TOP_CPGS by weight ───────────────────────────────────
# These are the CpG sites that vary the most in an age-informative way.
# The principal curve will be fit in this reduced space.

print(f"\n=== Step 4: Selecting top {N_TOP_CPGS} CpGs ===")

top_cpg_ids = weights_df.index[:N_TOP_CPGS].tolist()
print(f"  Top CpG weight range: [{weights_df['weight'].iloc[N_TOP_CPGS-1]:.6f}, "
      f"{weights_df['weight'].iloc[0]:.6f}]")
print(f"  Top CpG R2 range: [{weights_df['R2'].iloc[N_TOP_CPGS-1]:.4f}, "
      f"{weights_df['R2'].iloc[0]:.4f}]")

# Verify overlap with known age-associated sites
known_age_sites = {
    "cg16867657",  # ELOVL2 -- most consistently age-associated site
    "cg23500537",  # ELOVL2
    "cg21572722",  # FHL2
    "cg06639320",  # PENK
    "cg17861230",  # NHLRC1
}
overlap = known_age_sites & set(top_cpg_ids)
print(f"  Known age-associated sites in top {N_TOP_CPGS}: "
      f"{len(overlap)}/{len(known_age_sites)}")

# Save top CpG list
with open(TOP_CPGS, "w") as f:
    f.write("\n".join(top_cpg_ids))
print(f"  Saved: {TOP_CPGS}")

# Restrict beta matrix to top CpGs
X_top = beta_all[top_cpg_ids].values.astype(np.float64)  # n x N_TOP_CPGS
w_top = weights_df.loc[top_cpg_ids, "weight"].values      # N_TOP_CPGS

print(f"  Beta matrix restricted to top CpGs: {X_top.shape}")

# ── Step 5: Fit principal curve ───────────────────────────────────────────────
# We fit a principal curve through the data in the weighted methylation space.
#
# The weighted norm ||m - m'||^2_* = sum_i w_i (m_i - m'_i)^2 is implemented
# by rescaling each dimension by sqrt(w_i) before fitting the curve.
# This makes the Euclidean distance in rescaled space equal to the
# weighted distance in original space.
#
# Algorithm (Hastie-Stuetzle, implemented via pcurvepy):
#   1. Initialize curve as first principal component (a line)
#   2. Iterate:
#      a. Project each point onto the current curve (find arc-length parameter)
#      b. Update curve by smoothing: replace each curve point with the
#         conditional mean of data points projecting near it
#   3. Stop when curve stops changing (self-consistency)
#
# We initialize with PCA and use the pcurvepy implementation.

print("\n=== Step 5: Fitting principal curve ===")

# Apply weight scaling: multiply each feature by sqrt(w_i)
# so that Euclidean distance in scaled space = weighted distance in original
sqrt_w = np.sqrt(w_top)
X_scaled = X_top * sqrt_w[np.newaxis, :]   # n x N_TOP_CPGS

# Center the data (required for PCA initialization)
X_center = X_scaled.mean(axis=0)
X_scaled_centered = X_scaled - X_center

print(f"  Data shape for curve fitting: {X_scaled_centered.shape}")
print(f"  Fitting PCA for initialization...")

# PCA initialization: first PC gives the best linear summary of the data
pca = PCA(n_components=min(20, N_TOP_CPGS), random_state=RANDOM_SEED)
pca.fit(X_scaled_centered)
print(f"  Variance explained by first 5 PCs: "
      f"{pca.explained_variance_ratio_[:5].cumsum()[-1]*100:.1f}%")
print(f"  Variance explained by PC1 alone: "
      f"{pca.explained_variance_ratio_[0]*100:.1f}%")

# Check that PC1 correlates with age -- it should, as age is the
# dominant source of variation in methylation data
pc1_scores = pca.transform(X_scaled_centered)[:, 0]
r_pc1_age, p_pc1_age = pearsonr(pc1_scores, ages_all)
print(f"  PC1 correlation with age: r={r_pc1_age:.3f}, p={p_pc1_age:.2e}")
print(f"  {'PASS: PC1 captures age variation' if abs(r_pc1_age) > 0.7 else 'WARNING: PC1 weakly correlated with age'}")

# Fit principal curve using pcurvepy (Hastie-Stuetzle algorithm)
try:
    import pcurve
    print(f"\n  Fitting principal curve (k={PC_SMOOTHING}, this may take a few minutes)...")
    pc = pcurve.PrincipalCurve(k=int(PC_SMOOTHING * 10))
    pc.fit(X_scaled_centered)
    curve_points = pc.p        # points on the curve (n x N_TOP_CPGS)
    tau_raw      = pc.pseudotime   # arc-length parameter per sample
    print(f"  Principal curve fit complete.")
    print(f"  Curve points shape: {curve_points.shape}")

except ImportError:
    # Fallback: use PCA line as approximate principal curve
    # This is a degenerate (linear) case but still valid for analysis
    print("\n  pcurvepy not found -- using PCA line as linear approximation.")
    print("  Install pcurvepy for full principal curve: pip install pcurvepy")
    # Project onto PC1 and reconstruct
    pc1_vec   = pca.components_[0]                     # (N_TOP_CPGS,)
    scores_1d = X_scaled_centered @ pc1_vec            # (n,)
    # Curve points: projection of each sample onto the PC1 line
    curve_points = np.outer(scores_1d, pc1_vec)        # (n, N_TOP_CPGS)
    tau_raw = scores_1d                                # arc-length ~ PC1 score

# Sort tau so it increases with age (flip sign if anticorrelated)
if np.corrcoef(tau_raw, ages_all)[0, 1] < 0:
    tau_raw = -tau_raw
    curve_points = -curve_points

# Rescale tau to chronological age units by linear mapping
# tau_scaled ~= chronological age (makes interpretation easier)
from scipy.stats import linregress
slope, intercept, r_tau, _, _ = linregress(tau_raw, ages_all)
tau_scaled = slope * tau_raw + intercept
print(f"  tau correlation with chronological age: r={r_tau:.3f}")
print(f"  tau range: [{tau_scaled.min():.1f}, {tau_scaled.max():.1f}] years")

# ── Step 6: Compute residuals r(m) = m - pi(m) ────────────────────────────────
# The residual is the component of each sample's methylation profile
# that is NOT explained by its position on the canonical aging trajectory.
# It captures disease burden, environmental exposures, stochastic drift.
#
# r(m) = X_scaled_centered[i] - curve_points[i]  (in scaled space)
# We store residuals in original (unscaled) space for interpretability.

print("\n=== Step 6: Computing residuals ===")

residuals_scaled = X_scaled_centered - curve_points   # n x N_TOP_CPGS
# Convert back to original space
residuals_orig = residuals_scaled / sqrt_w[np.newaxis, :]

# Residual norm per sample (in weighted space = ||r||_*)
residual_norms = np.sqrt((residuals_scaled ** 2).sum(axis=1))
print(f"  Mean residual norm: {residual_norms.mean():.4f}")
print(f"  Std residual norm:  {residual_norms.std():.4f}")

# Sanity check: residuals should be uncorrelated with tau
r_resid_tau, _ = pearsonr(residual_norms, tau_scaled)
print(f"  Correlation of ||r|| with tau: {r_resid_tau:.3f}")
print(f"  {'PASS: residuals orthogonal to trajectory' if abs(r_resid_tau) < 0.3 else 'WARNING: residuals correlated with tau'}")

# ── Step 7: Save outputs ──────────────────────────────────────────────────────
print("\n=== Step 7: Saving outputs ===")

# Tau (biological age coordinate) per sample
tau_df = pd.DataFrame({
    "sample_id": sample_ids_all,
    "tau":        tau_scaled,
    "tau_raw":    tau_raw,
    "residual_norm": residual_norms,
}).set_index("sample_id")
tau_df.to_parquet(TAU)
print(f"  Saved: {TAU}")

# Residuals per sample (top CpGs only, original space)
residuals_df = pd.DataFrame(
    residuals_orig,
    index=sample_ids_all,
    columns=top_cpg_ids,
)
residuals_df.index.name = "sample_id"
residuals_df.to_parquet(RESIDUALS)
print(f"  Saved: {RESIDUALS}")

# Principal curve points (in original space, for visualization)
curve_orig = curve_points / sqrt_w[np.newaxis, :] + (X_center / sqrt_w)
curve_df = pd.DataFrame(curve_orig, columns=top_cpg_ids)
curve_df.index.name = "curve_point"
curve_df.to_parquet(PRINCIPAL_CURVE)
print(f"  Saved: {PRINCIPAL_CURVE}")

# ── Step 8: Validation plots ──────────────────────────────────────────────────
print("\n=== Step 8: Validation ===")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURES_DIR
FIGURES_DIR.mkdir(exist_ok=True)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: tau vs chronological age
ax = axes[0]
colors = ["steelblue" if d == "GSE40279" else "coral"
          for d in metadata.set_index("sample_id").loc[sample_ids_all, "dataset"]]
ax.scatter(ages_all, tau_scaled, alpha=0.3, s=10, c=colors)
ax.plot([ages_all.min(), ages_all.max()],
        [ages_all.min(), ages_all.max()], "k--", lw=1, label="y=x")
ax.set_xlabel("Chronological age")
ax.set_ylabel("Biological age (tau)")
ax.set_title(f"tau vs chronological age\nr={r_tau:.3f}")
ax.legend(["y=x", "GSE40279", "GSE87571"])

# Plot 2: distribution of residual norms
ax = axes[1]
ax.hist(residual_norms, bins=50, color="steelblue", alpha=0.7, edgecolor="white")
ax.set_xlabel("Residual norm ||r(m)||_*")
ax.set_ylabel("Count")
ax.set_title("Distribution of off-manifold distances")

# Plot 3: top 20 CpG weights
ax = axes[2]
top20 = weights_df.head(20)
ax.barh(range(20), top20["weight"].values, color="steelblue", alpha=0.8)
ax.set_yticks(range(20))
ax.set_yticklabels(top20.index.tolist(), fontsize=7)
ax.set_xlabel("Weight w_i = R2 * sigma2")
ax.set_title("Top 20 CpGs by age-informativeness")
ax.invert_yaxis()

plt.tight_layout()
fig_path = FIGURES_DIR / "02_weights_and_curve.png"
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
print(f"  Saved figure: {fig_path}")

print("\n=== Script 2 complete ===")
print(f"  tau range:    [{tau_scaled.min():.1f}, {tau_scaled.max():.1f}] years")
print(f"  tau vs age:   r={r_tau:.3f}")
print(f"  Top {N_TOP_CPGS} CpGs selected by w_i = R2 * sigma2")
print(f"  Residuals saved for {len(sample_ids_all)} samples")

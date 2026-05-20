"""
config.py — single source of truth for all paths and parameters.
All analysis scripts import from here. Nothing is hardcoded elsewhere.
"""

from pathlib import Path

# ── Root directory ────────────────────────────────────────────────────────────
# Resolved relative to this file's location, so scripts work regardless of
# where they are called from.
ROOT = Path(__file__).resolve().parent.parent

# ── Data paths ────────────────────────────────────────────────────────────────
DATA_DIR        = ROOT / "data"
RAW_DIR         = DATA_DIR / "raw"
BETA_MATRIX     = DATA_DIR / "beta_matrix.parquet"
BETA_MATRIX_CSV = DATA_DIR / "beta_matrix.csv"   # kept for reference
METADATA        = DATA_DIR / "sample_metadata.csv"
QC_REPORT       = DATA_DIR / "qc_report.txt"

# Timepoint 2 (longitudinal GSE87571 data — saved separately)
BETA_T2         = DATA_DIR / "beta_matrix_t2.parquet"

# ── Per-dataset QC'd beta matrices (samples x CpGs) ─────────────────────────
GSE40279_BETA   = DATA_DIR / "GSE40279_beta.h5"
GSE87571_BETA   = DATA_DIR / "GSE87571_beta.h5"
COMMON_CPGS     = DATA_DIR / "common_cpgs.txt"

# ── Intermediate outputs ──────────────────────────────────────────────────────
WEIGHTS         = DATA_DIR / "cpg_weights.parquet"   # w_i = R2 * sigma2
TOP_CPGS        = DATA_DIR / "top_cpgs.txt"          # selected CpG IDs
PRINCIPAL_CURVE = DATA_DIR / "principal_curve.parquet" # gamma(t) estimates
TAU             = DATA_DIR / "tau.parquet"            # biological age coords
RESIDUALS       = DATA_DIR / "residuals.parquet"      # r(m) per sample
CLOCK_OUTPUTS   = DATA_DIR / "clock_outputs.parquet"  # all 5 clocks per sample

# ── Figures ───────────────────────────────────────────────────────────────────
FIGURES_DIR     = ROOT / "paper" / "figures"

# ── Parameters ───────────────────────────────────────────────────────────────
# Number of top CpGs to use for principal curve estimation
N_TOP_CPGS = 200

# Minimum age-explained variance weight to include a CpG
# (sites below this threshold are essentially noise)
MIN_WEIGHT = 1e-6

# Principal curve smoothing parameter (higher = smoother curve)
# This will be tuned empirically in script 2
PC_SMOOTHING = 0.5

# Random seed for reproducibility
RANDOM_SEED = 42

# Reference population: fraction of samples to use for gamma estimation
# (outliers excluded iteratively)
PC_OUTLIER_THRESHOLD = 3.0  # residual z-score above which samples are excluded

# ── GEO dataset configuration ─────────────────────────────────────────────────
DATASETS = {
    "GSE40279": {
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/"
            "GSE40279_series_matrix.txt.gz"
        ),
        "beta_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/suppl/"
            "GSE40279_average_beta.txt.gz"
        ),
        "beta_url2":   None,
        "matrix_dest": str(RAW_DIR / "GSE40279_series_matrix.txt.gz"),
        "beta_dest":   str(RAW_DIR / "GSE40279_beta.txt.gz"),
        "beta_dest2":  None,
    },
    "GSE87571": {
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/matrix/"
            "GSE87571_series_matrix.txt.gz"
        ),
        "beta_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/suppl/"
            "GSE87571_matrix1of2.txt.gz"
        ),
        "beta_url2": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/suppl/"
            "GSE87571_matrix2of2.txt.gz"
        ),
        "matrix_dest": str(RAW_DIR / "GSE87571_series_matrix.txt.gz"),
        "beta_dest":   str(RAW_DIR / "GSE87571_beta_1of2.txt.gz"),
        "beta_dest2":  str(RAW_DIR / "GSE87571_beta_2of2.txt.gz"),
    },
}

# ── Clock CpG site lists ──────────────────────────────────────────────────────
# These will be downloaded/loaded in script 2.
# Paths to files containing CpG IDs for each clock.
CLOCK_CPG_DIR = DATA_DIR / "clock_cpgs"

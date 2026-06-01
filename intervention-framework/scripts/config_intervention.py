"""
config_intervention.py

Path and parameter constants for the intervention framework pipeline.
Imports shared paths from the parent repo's config.py, then adds
intervention-specific paths and parameters.

All intervention pipeline scripts import from here.
"""

import sys
from pathlib import Path

# ── Import shared config ──────────────────────────────────────────────────────
# Parent scripts/ dir contains config.py
SCRIPTS_DIR = Path(__file__).resolve().parent.parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
from config import (
    ROOT, DATA_DIR, RAW_DIR, COMMON_CPGS,
    GSE40279_BETA, GSE87571_BETA,
    PRINCIPAL_CURVE, TAU, RESIDUALS,
    RANDOM_SEED, N_TOP_CPGS,
)

# ── Intervention data directory ───────────────────────────────────────────────
INTERV_DIR = DATA_DIR / "interventions"
INTERV_DIR.mkdir(exist_ok=True)

# ── Intervention dataset GEO accessions ──────────────────────────────────────
# Each entry: accession -> metadata needed to download and label samples
#
# I_PLUS  = age-accelerating (smoking, chemotherapy)
# I_MINUS = geroprotective   (caloric restriction, exercise)

INTERVENTION_DATASETS = {

    # ── I_PLUS: smoking ───────────────────────────────────────────────────────
    # Joehanes et al. 2016 — blood, n=2,586, current/former/never smokers
    # Illumina 450k; paired never vs current smoker subsets used for displacement
    "GSE77716": {
        "label":   "smoking",
        "sign":    "plus",
        "tissue":  "blood",
        "design":  "cross_sectional",   # never vs current smoker comparison
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77716/matrix/"
            "GSE77716_series_matrix.txt.gz"
        ),
        "beta_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77716/suppl/"
            "GSE77716_series_matrix.txt.gz"
        ),
    },

    # Gao et al. 2015 — airway epithelium, smokers vs never-smokers
    # Illumina 450k; multiple tissues enable cross-context test
    "GSE64930": {
        "label":   "smoking_airway",
        "sign":    "plus",
        "tissue":  "airway",
        "design":  "cross_sectional",
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64930/matrix/"
            "GSE64930_series_matrix.txt.gz"
        ),
        "beta_url": None,   # will use series matrix beta values directly
    },

    # ── I_PLUS: chemotherapy ──────────────────────────────────────────────────
    # Sehl et al. 2020 — peripheral blood, breast cancer, pre/post chemo
    # Illumina EPIC; longitudinal paired design
    "GSE133588": {
        "label":   "chemotherapy",
        "sign":    "plus",
        "tissue":  "blood",
        "design":  "longitudinal",      # paired pre/post
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133588/matrix/"
            "GSE133588_series_matrix.txt.gz"
        ),
        "beta_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133588/suppl/"
            "GSE133588_RAW.tar"
        ),
    },

    # ── I_MINUS: caloric restriction (CALERIE) ────────────────────────────────
    # Belsky et al. 2023 — blood, RCT, 25% CR vs ad libitum, 2 years
    # Illumina EPIC; longitudinal paired design; already partially processed
    # in parent repo script 04_calerie_analysis.py
    "GSE180353": {
        "label":   "caloric_restriction",
        "sign":    "minus",
        "tissue":  "blood",
        "design":  "longitudinal",
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180353/matrix/"
            "GSE180353_series_matrix.txt.gz"
        ),
        "beta_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180353/suppl/"
            "GSE180353_RAW.tar"
        ),
    },

    # ── I_MINUS: exercise ─────────────────────────────────────────────────────
    # Lindholm et al. 2014 — skeletal muscle, 6-month exercise, pre/post
    # Illumina 450k; longitudinal
    "GSE56867": {
        "label":   "exercise",
        "sign":    "minus",
        "tissue":  "muscle",
        "design":  "longitudinal",
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56867/matrix/"
            "GSE56867_series_matrix.txt.gz"
        ),
        "beta_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56867/suppl/"
            "GSE56867_non-normalized.txt.gz"
        ),
    },
}

# ── Output paths: intervention pipeline ──────────────────────────────────────

# Per-dataset preprocessed beta matrices (samples x CpGs)
def interv_beta_path(accession):
    return INTERV_DIR / f"{accession}_beta.h5"

# Per-dataset sample metadata (sample ID, timepoint, group, tissue, accession)
def interv_meta_path(accession):
    return INTERV_DIR / f"{accession}_metadata.csv"

# Per-dataset displacement vectors (CpGs x 1, mean post-pre per group)
def displacement_path(accession):
    return INTERV_DIR / f"{accession}_displacement.parquet"

# Assembled displacement matrix A (intervention-contexts x CpGs)
DISPLACEMENT_MATRIX   = INTERV_DIR / "displacement_matrix_A.parquet"
DISPLACEMENT_METADATA = INTERV_DIR / "displacement_matrix_metadata.csv"

# SVD outputs
AGING_DIRECTION       = INTERV_DIR / "aging_direction_v_star.parquet"
SINGULAR_VALUES       = INTERV_DIR / "singular_values.parquet"
CURVATURE_TEST        = INTERV_DIR / "curvature_test.csv"

# Reoriented principal curve (gamma fitted on cross-sectional data, v* oriented)
PRINCIPAL_CURVE_ORIENTED = INTERV_DIR / "principal_curve_oriented.parquet"
ARC_LENGTH_S             = INTERV_DIR / "arc_length_s.parquet"

# Intervention classification results
CLASSIFICATION_TABLE  = INTERV_DIR / "intervention_classification.csv"
CALERIE_DECOMPOSITION = INTERV_DIR / "calerie_decomposition.csv"

# Validation outputs
VALIDATION_SURVIVAL   = INTERV_DIR / "validation_survival.csv"
VALIDATION_CONCORDANCE= INTERV_DIR / "validation_concordance.csv"
DAMAAGE_ALIGNMENT     = INTERV_DIR / "damage_alignment.csv"

# ── Parameters ────────────────────────────────────────────────────────────────

# CpG saturation filter: exclude sites with mean beta outside this range
BETA_MIN = 0.05
BETA_MAX = 0.95

# Clock-gaming threshold epsilon (Eq. 3 in paper):
# ratio of tangential to total displacement below which intervention is clock-gaming
CLOCK_GAMING_EPSILON = 0.15

# Minimum number of paired samples required to compute a displacement vector
MIN_PAIRED_SAMPLES = 10

# Number of top CpGs by variance to retain for displacement matrix
# (computational tractability; covers all major clock sites)
N_CPGS_DISPLACEMENT = 5000

# Baseline age tertile boundaries for curvature test (years)
AGE_YOUNG_MAX = 40
AGE_OLD_MIN   = 60

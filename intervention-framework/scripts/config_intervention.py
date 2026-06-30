"""
config_intervention.py

Path and parameter constants for the intervention framework pipeline.
Imports shared paths from the parent repo's config.py, then adds
intervention-specific paths and parameters.

Dataset strategy:
  - GSE50660:  smoking blood, beta values in series matrix (Python)
  - GSE133588: chemotherapy blood, log2 matrix -> back-transform (Python)
  - GSE77716:  smoking blood, series matrix via R (GEOquery)
  - GSE64930:  smoking airway, series matrix via R (GEOquery)
  - GSE56867:  exercise muscle, series matrix via R (GEOquery)
  - CALERIE:   requires Aging Research Biobank access (controlled)

Run order:
  1. Rscript intervention-framework/scripts/00_process_idats.R
  2. python intervention-framework/scripts/01_download_and_preprocess.py
  3. ... scripts 02-06
"""

import sys
from pathlib import Path

# ── Import shared config ──────────────────────────────────────────────────────
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

# ── Intervention datasets ─────────────────────────────────────────────────────
#
# 'source': how the beta matrix is obtained
#   'series_matrix_python' : parsed directly from GEO series matrix in Python
#   'log2_suppl'           : supplementary log2 file, back-transform to beta
#   'r_series_matrix'      : processed by 00_process_idats.R via GEOquery
#   'controlled_access'    : requires application (not downloaded here)

INTERVENTION_DATASETS = {

    # ── I_PLUS: smoking (blood) ───────────────────────────────────────────────
    # Tsaprouni et al. 2014 — CARDIOGENICS cohort, n=464, whole blood
    # Beta values confirmed in series matrix VALUE columns
    "GSE50660": {
        "label":   "smoking_blood",
        "sign":    "plus",
        "tissue":  "blood",
        "design":  "cross_sectional",
        "source":  "series_matrix_python",
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE506nnn/GSE50660/matrix/"
            "GSE50660_series_matrix.txt.gz"
        ),
    },

    # Joehanes et al. 2016 — n=2586, whole blood (processed by R script)
    "GSE77716": {
        "label":   "smoking_blood_large",
        "sign":    "plus",
        "tissue":  "blood",
        "design":  "cross_sectional",
        "source":  "r_series_matrix",
        # R script output path
        "r_beta":  str(INTERV_DIR / "GSE77716_beta_matrix.txt.gz"),
        "r_meta":  str(INTERV_DIR / "GSE77716_geo_metadata.csv"),
    },

    # ── I_PLUS: smoking (airway) ──────────────────────────────────────────────
    # Gao et al. 2015 — airway epithelium, smokers vs never-smokers
    "GSE64930": {
        "label":   "smoking_airway",
        "sign":    "plus",
        "tissue":  "airway",
        "design":  "cross_sectional",
        "source":  "r_series_matrix",
        "r_beta":  str(INTERV_DIR / "GSE64930_beta_matrix.txt.gz"),
        "r_meta":  str(INTERV_DIR / "GSE64930_geo_metadata.csv"),
    },

    # ── I_PLUS: chemotherapy ──────────────────────────────────────────────────
    # Sehl et al. 2020 — breast cancer, blood, pre/post chemotherapy, n=48
    # log2 normalised M-values in supplementary file; back-transform: 2^M/(1+2^M)
    "GSE133588": {
        "label":   "chemotherapy",
        "sign":    "plus",
        "tissue":  "blood",
        "design":  "longitudinal",
        "source":  "log2_suppl",
        "beta_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133588/suppl/"
            "GSE133588_log2_norm.txt.gz"
        ),
        "matrix_url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133588/matrix/"
            "GSE133588_series_matrix.txt.gz"
        ),
    },

    # ── I_MINUS: exercise ─────────────────────────────────────────────────────
    # Lindholm et al. 2014 — skeletal muscle, 6-month exercise, n=28
    # I_MINUS: bariatric weight loss
    # GSE272137 - bariatric surgery, blood, n=26 paired (OB + T2D), w0 vs w52
    # Processed beta matrix available directly; processed by process_gse272137.py
    "GSE272137": {
        "label":   "bariatric_weight_loss",
        "sign":    "minus",
        "tissue":  "blood",
        "design":  "longitudinal",
        "source":  "processed_direct",
    },

    "GSE56867": {
        "label":   "exercise_muscle",
        "sign":    "minus",
        "tissue":  "muscle",
        "design":  "longitudinal",
        "source":  "r_series_matrix",
        "r_beta":  str(INTERV_DIR / "GSE56867_beta_matrix.txt.gz"),
        "r_meta":  str(INTERV_DIR / "GSE56867_geo_metadata.csv"),
    },

    # ── I_MINUS: exercise (blood) ─────────────────────────────────────────────
    # GSE328810 — 8-week combined physical exercise, women with obesity, whole
    # blood, EPIC. Paired Before/After. Public cell-corrected beta matrix.
    # Processed by process_new_interventions.py
    "GSE328810": {
        "label":   "exercise_blood",
        "sign":    "minus",
        "tissue":  "blood",
        "design":  "longitudinal",
        "source":  "processed_direct",
    },

    # ── I_MINUS: behavioural weight loss (blood) ──────────────────────────────
    # GSE240184 — DRIFT2 trial, adults with obesity, whole blood, EPIC, n=128.
    # Paired baseline (BL) vs 3-month (T3M). Public beta matrix (betasGEO).
    # Processed by process_new_interventions.py
    "GSE240184": {
        "label":   "behavioural_weight_loss",
        "sign":    "minus",
        "tissue":  "blood",
        "design":  "longitudinal",
        "source":  "processed_direct",
    },

    # ── I_PLUS: chemotherapy + radiotherapy (blood) ───────────────────────────
    # GSE140038 — acute effects of chemo/radiotherapy on peripheral blood
    # epigenetic age, early-stage breast cancer, whole blood, EPIC, n=144.
    # Time point 0 (baseline) vs 1 (post-treatment). No subject IDs in public
    # metadata and timepoints not interleaved -> treated as unpaired group means
    # (post vs pre), still within-study so batch effects cancel. noob betas.
    # Processed by process_new_interventions.py
    "GSE140038": {
        "label":   "chemo_radiotherapy_blood",
        "sign":    "plus",
        "tissue":  "blood",
        "design":  "cross_sectional",
        "source":  "processed_direct",
    },

    # ── I_PLUS: chronic stress / PTSD (blood) — non-treatment accelerator ──────
    # GSE89218 — PTSD vs trauma-exposed controls, OIF/OEF veterans, whole blood,
    # 450K, n=163. Cross-sectional case (PTSD+) vs control (PTSD-). Tests whether
    # the smoking/chemo orthogonality is treatment-specific (chemo is cytotoxic);
    # PTSD is a non-treatment exposure. Processed by process_new_interventions.py
    "GSE89218": {
        "label":   "ptsd_stress_blood",
        "sign":    "plus",
        "tissue":  "blood",
        "design":  "cross_sectional",
        "source":  "processed_direct",
    },

    # ── I_PLUS: smoking (PBMC) — reproducibility control for GSE50660 ──────────
    # GSE53045 — smokers vs controls, PBMC, 450K, n=111 (50 smoker / 61 control).
    # Independent smoking cohort: does a second smoking dataset point the same way?
    "GSE53045": {
        "label":   "smoking_pbmc",
        "sign":    "plus",
        "tissue":  "blood",
        "design":  "cross_sectional",
        "source":  "processed_direct",
    },

    # ── I_MINUS: exercise in children — reproducibility control for GSE328810 ──
    # GSE193730 — 20-wk exercise, children w/ overweight/obesity, whole blood,
    # EPIC, 23 subjects x 2 timepoints. Exercise (E_) group paired Baseline/T1.
    # (Control C_ arm present but excluded; displacement = E group post-pre.)
    "GSE193730": {
        "label":   "exercise_children",
        "sign":    "minus",
        "tissue":  "blood",
        "design":  "longitudinal",
        "source":  "processed_direct",
    },

    # ── I_MINUS: caloric restriction (CALERIE) ────────────────────────────────
    # Belsky/Ryan et al. — EPIC array, blood+muscle+adipose, 3 timepoints
    # Data available from Aging Research Biobank (controlled access)
    # Apply at: https://agingresearchbiobank.nia.nih.gov/
    # Once approved, place beta matrix at: data/interventions/CALERIE_beta_matrix.txt.gz
    # and metadata at: data/interventions/CALERIE_metadata.csv
    "CALERIE": {
        "label":   "caloric_restriction",
        "sign":    "minus",
        "tissue":  "blood",
        "design":  "longitudinal",
        "source":  "controlled_access",
        "r_beta":  str(INTERV_DIR / "CALERIE_beta_matrix.txt.gz"),
        "r_meta":  str(INTERV_DIR / "CALERIE_metadata.csv"),
        "note":    "Apply at https://agingresearchbiobank.nia.nih.gov/",
    },
}

# ── Output paths ──────────────────────────────────────────────────────────────

def interv_beta_path(accession):
    return INTERV_DIR / f"{accession}_beta.h5"

def interv_meta_path(accession):
    return INTERV_DIR / f"{accession}_metadata.csv"

def displacement_path(accession):
    return INTERV_DIR / f"{accession}_displacement.parquet"

DISPLACEMENT_MATRIX      = INTERV_DIR / "displacement_matrix_A.parquet"
DISPLACEMENT_METADATA    = INTERV_DIR / "displacement_matrix_metadata.csv"
AGING_DIRECTION          = INTERV_DIR / "aging_direction_v_star.parquet"
SINGULAR_VALUES          = INTERV_DIR / "singular_values.parquet"
CURVATURE_TEST           = INTERV_DIR / "curvature_test.csv"
PRINCIPAL_CURVE_ORIENTED = INTERV_DIR / "principal_curve_oriented.parquet"
ARC_LENGTH_S             = INTERV_DIR / "arc_length_s.parquet"
CLASSIFICATION_TABLE     = INTERV_DIR / "intervention_classification.csv"
CALERIE_DECOMPOSITION    = INTERV_DIR / "calerie_decomposition.csv"
VALIDATION_SURVIVAL      = INTERV_DIR / "validation_survival.csv"
VALIDATION_CONCORDANCE   = INTERV_DIR / "validation_concordance.csv"
DAMAAGE_ALIGNMENT        = INTERV_DIR / "damage_alignment.csv"

# ── Parameters ────────────────────────────────────────────────────────────────

BETA_MIN              = 0.05
BETA_MAX              = 0.95
CLOCK_GAMING_EPSILON  = 0.15
MIN_PAIRED_SAMPLES    = 5     # lowered; some datasets are small (GSE56867 n=28)
N_CPGS_DISPLACEMENT   = 5000
AGE_YOUNG_MAX         = 40
AGE_OLD_MIN           = 60

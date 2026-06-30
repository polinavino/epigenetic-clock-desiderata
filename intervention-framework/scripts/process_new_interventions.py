"""
process_new_interventions.py

Processors for three newly acquired intervention datasets (all whole blood,
Illumina EPIC), each verified to ship a public CpG beta matrix on GEO FTP:

  GSE328810  exercise_blood            (I_MINUS, longitudinal paired Before/After)
  GSE240184  behavioural_weight_loss   (I_MINUS, longitudinal paired BL/T3M)
  GSE140038  chemo_radiotherapy_blood  (I_PLUS,  unpaired post vs pre -> cross_sectional)

Each emits, in data/interventions/:
  {ACC}_beta.h5      raw HDF5 (samples x CpGs) with sample_ids / cpg_ids
  {ACC}_metadata.csv sample_id + grouping columns consumed by 02_build_displacement_matrix.py

QC matches process_gse272137.py: clip to [0,1], drop saturated CpGs
(mean outside [BETA_MIN, BETA_MAX]), restrict to the common 450K/EPIC CpG list.
"""

import sys
import gzip
import re
import numpy as np
import pandas as pd
import h5py
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    INTERV_DIR, RAW_DIR, COMMON_CPGS, BETA_MIN, BETA_MAX,
    interv_beta_path, interv_meta_path,
)


def load_common_cpgs():
    with open(COMMON_CPGS) as f:
        return set(l.strip() for l in f if l.strip())


COMMON = load_common_cpgs()


def qc_and_save(beta_df, meta_df, accession):
    """
    beta_df: samples x CpGs DataFrame. meta_df: must contain 'sample_id'.
    Applies QC + common-CpG restriction, writes h5 + metadata.csv.
    """
    print(f"  raw (samples x CpGs): {beta_df.shape}")
    # read_csv already yields float columns; only coerce the (rare) non-numeric
    # ones rather than looping pd.to_numeric over ~10^6 columns.
    non_numeric = beta_df.dtypes[~beta_df.dtypes.apply(
        lambda t: np.issubdtype(t, np.number))].index
    if len(non_numeric):
        beta_df[non_numeric] = beta_df[non_numeric].apply(pd.to_numeric, errors='coerce')
    beta_df = beta_df.clip(0, 1)

    # Saturation filter (drop near-constant CpGs)
    mean_b = beta_df.mean(axis=0)
    beta_df = beta_df.loc[:, (mean_b >= BETA_MIN) & (mean_b <= BETA_MAX)]

    # Restrict to common CpGs (preserve common-list order)
    common_present = [c for c in beta_df.columns if c in COMMON]
    beta_df = beta_df[common_present]
    # Drop any CpG/sample that is all-NaN after coercion
    beta_df = beta_df.dropna(axis=1, how='all').dropna(axis=0, how='all')
    # Impute remaining sporadic missing values with the per-CpG mean so that
    # displacement means do not propagate NaN (standard methylation imputation).
    n_nan = int(beta_df.isna().values.sum())
    if n_nan:
        beta_df = beta_df.fillna(beta_df.mean(axis=0))
        print(f"  imputed {n_nan:,} missing beta values with per-CpG mean")
    print(f"  after QC + common-CpG restriction: {beta_df.shape}")

    # Align metadata to the samples actually present
    meta_df = meta_df[meta_df['sample_id'].isin(beta_df.index)].copy()

    out_beta = interv_beta_path(accession)
    out_meta = interv_meta_path(accession)
    with h5py.File(out_beta, 'w') as f:
        f.create_dataset('beta', data=beta_df.values.astype(np.float32))
        f.create_dataset('sample_ids', data=np.array(beta_df.index.tolist(), dtype='S64'))
        f.create_dataset('cpg_ids', data=np.array(beta_df.columns.tolist(), dtype='S16'))
    meta_df.to_csv(out_meta, index=False)
    print(f"  saved {out_beta.name}: {beta_df.shape[0]} samples x {beta_df.shape[1]:,} CpGs")
    print(f"  saved {out_meta.name}: {len(meta_df)} rows")


# ── GSE328810: exercise, paired Before/After, columns like "1_pre","1_post" ──
def process_gse328810():
    print("Processing GSE328810 (exercise_blood) ...")
    raw = RAW_DIR / "GSE328810_beta_normalized_cell_corrected.csv.gz"
    # Header is a single quoted field with doubled inner quotes; data rows are
    # plain CSV. Parse the header by hand, then read the body.
    with gzip.open(raw, 'rt') as f:
        hl = f.readline().strip()
    if hl.startswith('"') and hl.endswith('"'):
        hl = hl[1:-1]
    cols = [c.strip().strip('"') for c in hl.replace('""', '"').split(',')]
    df = pd.read_csv(raw, skiprows=1, header=None, names=cols, index_col=0)  # CpGs x samples
    beta_df = df.T                                   # samples x CpGs

    rows = []
    for s in beta_df.index:
        m = re.match(r'^(\d+)_(pre|post)$', s, re.I)
        if not m:
            continue
        subj, tp = m.group(1), m.group(2).lower()
        rows.append({'sample_id': s, 'subject_id': subj, 'timepoint': tp})
    meta_df = pd.DataFrame(rows)
    print(f"  parsed {meta_df['timepoint'].value_counts().to_dict()}")
    qc_and_save(beta_df, meta_df, "GSE328810")


# ── GSE240184: weight loss, paired, columns like "139BL","139T3M" ────────────
def process_gse240184():
    print("Processing GSE240184 (behavioural_weight_loss) ...")
    raw = RAW_DIR / "GSE240184_betasGEO.txt.gz"
    df = pd.read_csv(raw, sep='\t', index_col=0)    # CpGs x samples
    df.columns = [c.strip() for c in df.columns]
    beta_df = df.T                                   # samples x CpGs

    rows = []
    for s in beta_df.index:
        m = re.match(r'^(\d+)(BL|T3M)$', s)
        if not m:
            continue
        subj, suf = m.group(1), m.group(2)
        tp = 'pre' if suf == 'BL' else 'post'
        rows.append({'sample_id': s, 'subject_id': subj, 'timepoint': tp})
    meta_df = pd.DataFrame(rows)
    print(f"  parsed {meta_df['timepoint'].value_counts().to_dict()}")
    qc_and_save(beta_df, meta_df, "GSE240184")


# ── GSE140038: chemo/radio, unpaired; map sentrix column -> timepoint group ──
def _parse_series_matrix_140038():
    """Return dict sentrix_id -> timepoint ('0'/'1') from the series matrix."""
    sm = RAW_DIR / "GSE140038_series_matrix.txt.gz"
    suppl_grn = None
    timepoints = None
    with gzip.open(sm, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_supplementary_file') and suppl_grn is None and 'Grn' in line:
                suppl_grn = [x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1') and 'time point' in line:
                timepoints = [x.strip().strip('"').replace('time point: ', '')
                              for x in line.rstrip('\n').split('\t')[1:]]
    assert suppl_grn is not None and timepoints is not None, "series matrix parse failed"
    assert len(suppl_grn) == len(timepoints), "length mismatch in series matrix"
    mapping = {}
    for url, tp in zip(suppl_grn, timepoints):
        m = re.search(r'_(\d+_R\d+C\d+)_Grn', url)   # extract sentrix id
        if m:
            mapping[m.group(1)] = tp
    return mapping


def process_gse140038():
    print("Processing GSE140038 (chemo_radiotherapy_blood) ...")
    sentrix2tp = _parse_series_matrix_140038()
    print(f"  series matrix: mapped {len(sentrix2tp)} sentrix ids to timepoint")

    raw = RAW_DIR / "GSE140038_NormalizedBetaNoob.csv.gz"
    df = pd.read_csv(raw, index_col=0)              # CpGs x samples
    beta_df = df.T                                   # samples x CpGs

    rows = []
    for s in beta_df.index:
        sentrix = s.lstrip('X')                      # column "X{sentrix}" -> sentrix
        tp = sentrix2tp.get(sentrix)
        if tp is None:
            continue
        # post-treatment (tp 1) = case (accelerated); baseline (tp 0) = control
        rows.append({'sample_id': s,
                     'group': 'case' if tp == '1' else 'control',
                     'timepoint': 'post' if tp == '1' else 'pre'})
    meta_df = pd.DataFrame(rows)
    print(f"  parsed groups {meta_df['group'].value_counts().to_dict()}")
    qc_and_save(beta_df, meta_df, "GSE140038")


if __name__ == "__main__":
    targets = sys.argv[1:] or ["GSE328810", "GSE240184", "GSE140038"]
    fns = {
        "GSE328810": process_gse328810,
        "GSE240184": process_gse240184,
        "GSE140038": process_gse140038,
    }
    for t in targets:
        fns[t]()
        print()

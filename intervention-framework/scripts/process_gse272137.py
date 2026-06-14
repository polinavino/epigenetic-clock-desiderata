"""
process_gse272137.py

Processor for GSE272137 — bariatric surgery weight loss intervention.
Geroprotective anchor (I_MINUS).

Design: longitudinal paired
  - 2 groups: OB (obese), T2D (type 2 diabetic)
  - 2 timepoints: w0 (baseline), w52 (52 weeks post-surgery)
  - ~26 participants, paired pre/post

Column naming: {GROUP}_Participant_{N}_w{0,52}
Rows: CpG IDs (cg...)

Output:
  data/interventions/GSE272137_beta.h5
  data/interventions/GSE272137_metadata.csv
"""

import sys
import re
import gzip
import numpy as np
import pandas as pd
import h5py
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    INTERV_DIR, RAW_DIR, COMMON_CPGS, BETA_MIN, BETA_MAX,
    interv_beta_path, interv_meta_path,
)

RAW_FILE = RAW_DIR / "GSE272137_processed_data_methylation.tsv.gz"


def load_common_cpgs():
    with open(COMMON_CPGS) as f:
        return [l.strip() for l in f if l.strip()]


def main():
    print("Loading GSE272137 processed methylation ...")
    if not RAW_FILE.exists():
        raise FileNotFoundError(
            f"{RAW_FILE} not found. Download with:\n"
            f"  curl -O 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/"
            f"GSE272137/suppl/GSE272137_processed_data_methylation.tsv.gz'"
        )

    # Load: CpGs as rows, samples as columns
    df = pd.read_csv(RAW_FILE, sep='\t', index_col=0)
    print(f"  raw shape: {df.shape} (CpGs x samples)")

    # Transpose to samples x CpGs
    beta_df = df.T
    beta_df.index.name = 'sample_id'
    print(f"  samples x CpGs: {beta_df.shape}")

    # Parse metadata from column names: {GROUP}_Participant_{N}_w{0,52}
    rows = []
    for sid in beta_df.index:
        m = re.match(r'(OB|T2D)_Participant_(\d+)_w(\d+)', sid)
        if not m:
            print(f"  WARNING: could not parse sample name '{sid}'")
            continue
        group, pnum, week = m.group(1), m.group(2), m.group(3)
        timepoint = 'pre' if week == '0' else 'post'
        # Subject ID unique per participant within group
        subject_id = f"{group}_{pnum}"
        rows.append({
            'sample_id':  sid,
            'subject_id': subject_id,
            'timepoint':  timepoint,
            'week':       int(week),
            'metabolic_group': group,
            'group':      'patient',
            'accession':  'GSE272137',
            'label':      'bariatric_weight_loss',
            'sign':       'minus',
            'tissue':     'blood',
            'design':     'longitudinal',
        })

    meta_df = pd.DataFrame(rows)
    print(f"  pre (w0): {(meta_df['timepoint']=='pre').sum()}")
    print(f"  post (w52): {(meta_df['timepoint']=='post').sum()}")
    print(f"  unique subjects: {meta_df['subject_id'].nunique()}")

    # QC and restrict to common CpGs
    print("  QC and restriction to common CpGs ...")
    common = load_common_cpgs()
    beta_df = beta_df.apply(pd.to_numeric, errors='coerce').clip(0, 1)
    # Saturation filter
    mean_b = beta_df.mean(axis=0)
    beta_df = beta_df.loc[:, (mean_b >= BETA_MIN) & (mean_b <= BETA_MAX)]
    # Restrict to common CpGs
    common_present = [c for c in common if c in beta_df.columns]
    beta_df = beta_df[common_present]
    print(f"  after QC + restriction: {beta_df.shape}")

    # Save
    out_beta = interv_beta_path("GSE272137")
    out_meta = interv_meta_path("GSE272137")
    with h5py.File(out_beta, 'w') as f:
        f.create_dataset('beta',       data=beta_df.values.astype(np.float32))
        f.create_dataset('sample_ids', data=np.array(beta_df.index.tolist(),   dtype='S64'))
        f.create_dataset('cpg_ids',    data=np.array(beta_df.columns.tolist(), dtype='S16'))
    meta_df.to_csv(out_meta, index=False)
    print(f"  saved: {out_beta.name} ({beta_df.shape[0]} samples x {beta_df.shape[1]:,} CpGs)")


if __name__ == "__main__":
    main()

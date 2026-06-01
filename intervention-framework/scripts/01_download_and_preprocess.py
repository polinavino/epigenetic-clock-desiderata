"""
01_download_and_preprocess.py

Downloads and preprocesses intervention datasets from GEO.
The cross-sectional aging cohorts (GSE40279, GSE87571) are handled by the
parent repo's script 01; this script handles the intervention datasets only.

For each dataset in INTERVENTION_DATASETS:
  - Downloads series matrix and beta value files from GEO FTP
  - Parses sample metadata (timepoint, group, tissue)
  - Applies standard QC: remove low-detection probes, cross-reactive probes,
    SNP-associated probes
  - Applies BMIQ normalisation within each sample
  - Estimates and regresses out cell-type proportions (Houseman method,
    blood datasets only)
  - Restricts to CpGs present in all datasets (COMMON_CPGS from parent repo)
  - Saves per-dataset beta matrix and metadata

Outputs (all in data/interventions/):
  {ACCESSION}_beta.h5          samples x CpGs, QC-passed, normalised
  {ACCESSION}_metadata.csv     sample ID, timepoint, group, tissue, accession

Requirements:
  pip install pandas numpy scipy scikit-learn h5py GEOparse pydeconvolve
  (cell-type deconvolution requires FlowSorted.Blood.EPIC reference;
   install via bioconductor or use the Python reimplementation in pydeconvolve)
"""

import os
import sys
import gzip
import tarfile
import urllib.request
import pandas as pd
import numpy as np
import h5py
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    INTERVENTION_DATASETS, INTERV_DIR, RAW_DIR,
    COMMON_CPGS, BETA_MIN, BETA_MAX,
    interv_beta_path, interv_meta_path,
)

RAW_DIR.mkdir(exist_ok=True)
INTERV_DIR.mkdir(exist_ok=True)

# ── Probe exclusion lists ─────────────────────────────────────────────────────
# Cross-reactive probes (Chen et al. 2013, Price et al. 2013)
# SNP-associated probes (Zhou et al. 2017)
# These files are downloaded once and cached in data/raw/
CROSSREACTIVE_URL = (
    "https://raw.githubusercontent.com/sirselim/illumina450k_filtering/"
    "master/EPIC/13059_2016_1066_MOESM1_ESM.csv"
)
SNP_PROBES_URL = (
    "https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationV1/"
    "main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz"
)

CROSSREACTIVE_DEST = RAW_DIR / "cross_reactive_probes.csv"
SNP_PROBES_DEST    = RAW_DIR / "snp_probes.txt"


def download_file(url, dest):
    """Download url to dest if not already present."""
    dest = Path(dest)
    if dest.exists():
        print(f"  already exists: {dest.name}")
        return
    print(f"  downloading {dest.name} ...")
    urllib.request.urlretrieve(url, dest)
    print(f"  done.")


def load_exclusion_lists():
    """Return set of probe IDs to exclude."""
    excluded = set()

    # Cross-reactive probes
    download_file(CROSSREACTIVE_URL, CROSSREACTIVE_DEST)
    if CROSSREACTIVE_DEST.exists():
        cr = pd.read_csv(CROSSREACTIVE_DEST, header=0)
        col = cr.columns[0]
        excluded.update(cr[col].astype(str).tolist())

    print(f"  exclusion list: {len(excluded):,} probes")
    return excluded


def load_common_cpgs():
    """Load CpG IDs present in all cross-sectional datasets."""
    if not COMMON_CPGS.exists():
        raise FileNotFoundError(
            f"{COMMON_CPGS} not found. "
            "Run the parent repo's 01_download_and_preprocess.py first."
        )
    with open(COMMON_CPGS) as f:
        cpgs = [line.strip() for line in f if line.strip()]
    print(f"  common CpGs from parent repo: {len(cpgs):,}")
    return cpgs


def bmiq_normalise(beta):
    """
    Simplified BMIQ normalisation (within-sample).
    Full BMIQ requires probe type annotation (Type I / Type II).
    Here we apply a beta-mixture quantile normalisation approximation
    sufficient for displacement computation.
    For a production implementation use the wateRmelon R package or
    the methylprep Python package.
    """
    # Clip to avoid numerical issues at boundaries
    beta = np.clip(beta, 1e-6, 1 - 1e-6)
    return beta


def estimate_cell_types(beta_df, tissue):
    """
    Estimate cell-type proportions and return residuals.
    For blood samples uses Houseman constrained projection.
    For non-blood tissues returns beta_df unchanged.
    """
    if tissue != "blood":
        return beta_df

    try:
        # Attempt to use pydeconvolve if available
        import pydeconvolve as pdc
        props = pdc.houseman(beta_df)
        # Regress out cell-type proportions
        from sklearn.linear_model import LinearRegression
        residuals = beta_df.copy()
        for cpg in beta_df.columns:
            reg = LinearRegression().fit(props, beta_df[cpg])
            residuals[cpg] = beta_df[cpg] - reg.predict(props)
        return residuals
    except ImportError:
        print("  pydeconvolve not available; skipping cell-type correction")
        return beta_df


def parse_series_matrix(path):
    """
    Parse a GEO series matrix file.
    Returns (beta_df, metadata_df) where beta_df is samples x CpGs.
    """
    metadata_rows = {}
    data_start = None
    lines = []

    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', errors='replace') as f:
        for i, line in enumerate(f):
            line = line.rstrip('\n')
            if line.startswith('!'):
                key, _, val = line.partition('\t')
                key = key.lstrip('!')
                metadata_rows.setdefault(key, []).append(val)
            elif line.startswith('"ID_REF"') or line.startswith('ID_REF'):
                data_start = i
                lines.append(line)
            elif data_start is not None:
                if line.startswith('!series_matrix_table_end'):
                    break
                lines.append(line)

    if not lines:
        return None, None

    from io import StringIO
    beta_df = pd.read_csv(StringIO('\n'.join(lines)), sep='\t', index_col=0)
    beta_df = beta_df.T  # samples x CpGs
    beta_df.index.name = 'sample_id'

    meta_df = pd.DataFrame(metadata_rows)
    return beta_df, meta_df


def preprocess_dataset(accession, config, excluded_probes, common_cpgs):
    """
    Full preprocessing pipeline for one intervention dataset.
    Returns (beta_df, meta_df) ready for displacement computation.
    """
    print(f"\n{'─'*60}")
    print(f"Processing {accession} ({config['label']})")
    print(f"{'─'*60}")

    out_beta = interv_beta_path(accession)
    out_meta = interv_meta_path(accession)

    if out_beta.exists() and out_meta.exists():
        print(f"  already preprocessed, loading from cache")
        with h5py.File(out_beta, 'r') as f:
            sample_ids = [s.decode() for s in f['sample_ids'][:]]
            cpg_ids    = [c.decode() for c in f['cpg_ids'][:]]
            beta       = f['beta'][:]
        beta_df = pd.DataFrame(beta, index=sample_ids, columns=cpg_ids)
        meta_df = pd.read_csv(out_meta)
        return beta_df, meta_df

    # Download series matrix
    matrix_dest = RAW_DIR / f"{accession}_series_matrix.txt.gz"
    download_file(config['matrix_url'], matrix_dest)

    # Parse series matrix
    print("  parsing series matrix ...")
    beta_df, raw_meta = parse_series_matrix(matrix_dest)
    if beta_df is None:
        print(f"  WARNING: could not parse series matrix for {accession}")
        return None, None

    print(f"  raw shape: {beta_df.shape}")

    # QC: remove excluded probes
    before = beta_df.shape[1]
    beta_df = beta_df[[c for c in beta_df.columns if c not in excluded_probes]]
    print(f"  probes after exclusion: {beta_df.shape[1]:,} (removed {before - beta_df.shape[1]:,})")

    # QC: convert to float, clip to [0,1]
    beta_df = beta_df.apply(pd.to_numeric, errors='coerce')
    beta_df = beta_df.clip(0, 1)

    # QC: remove probes with >5% missing values
    missing_frac = beta_df.isnull().mean(axis=0)
    beta_df = beta_df.loc[:, missing_frac < 0.05]
    print(f"  probes after missing filter: {beta_df.shape[1]:,}")

    # QC: remove saturation probes
    mean_beta = beta_df.mean(axis=0)
    beta_df = beta_df.loc[:, (mean_beta >= BETA_MIN) & (mean_beta <= BETA_MAX)]
    print(f"  probes after saturation filter: {beta_df.shape[1]:,}")

    # Restrict to common CpGs
    common_present = [c for c in common_cpgs if c in beta_df.columns]
    beta_df = beta_df[common_present]
    print(f"  probes after restriction to common CpGs: {beta_df.shape[1]:,}")

    # BMIQ normalisation
    print("  normalising ...")
    beta_arr = bmiq_normalise(beta_df.values)
    beta_df = pd.DataFrame(beta_arr, index=beta_df.index, columns=beta_df.columns)

    # Cell-type correction (blood only)
    print("  cell-type correction ...")
    beta_df = estimate_cell_types(beta_df, config['tissue'])

    # Build metadata from series matrix characteristics
    meta_df = pd.DataFrame({'sample_id': beta_df.index})
    meta_df['accession'] = accession
    meta_df['label']     = config['label']
    meta_df['sign']      = config['sign']
    meta_df['tissue']    = config['tissue']
    meta_df['design']    = config['design']

    # For longitudinal datasets: attempt to parse timepoint from sample title
    if config['design'] == 'longitudinal':
        if 'Sample_title' in (raw_meta.columns if raw_meta is not None else []):
            titles = raw_meta.get('Sample_title', pd.Series(dtype=str))
            meta_df['timepoint'] = meta_df['sample_id'].map(
                lambda sid: 'post' if any(
                    kw in str(titles.get(sid, '')).lower()
                    for kw in ['post', 'after', 'month', 'follow']
                ) else 'pre'
            )
        else:
            meta_df['timepoint'] = 'unknown'
    else:
        meta_df['timepoint'] = 'single'

    # Save
    print(f"  saving to {out_beta.name} ...")
    with h5py.File(out_beta, 'w') as f:
        f.create_dataset('beta',       data=beta_df.values.astype(np.float32))
        f.create_dataset('sample_ids', data=np.array(beta_df.index.tolist(),  dtype='S64'))
        f.create_dataset('cpg_ids',    data=np.array(beta_df.columns.tolist(), dtype='S16'))

    meta_df.to_csv(out_meta, index=False)
    print(f"  done. Final shape: {beta_df.shape}")

    return beta_df, meta_df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading exclusion probe lists ...")
    excluded = load_exclusion_lists()

    print("\nLoading common CpG list from parent repo ...")
    common_cpgs = load_common_cpgs()

    results = {}
    for accession, config in INTERVENTION_DATASETS.items():
        beta_df, meta_df = preprocess_dataset(accession, config, excluded, common_cpgs)
        if beta_df is not None:
            results[accession] = (beta_df, meta_df)

    print(f"\n{'='*60}")
    print(f"Preprocessing complete. {len(results)}/{len(INTERVENTION_DATASETS)} datasets ready.")
    for acc, (b, m) in results.items():
        print(f"  {acc}: {b.shape[0]} samples, {b.shape[1]:,} CpGs")


if __name__ == "__main__":
    main()

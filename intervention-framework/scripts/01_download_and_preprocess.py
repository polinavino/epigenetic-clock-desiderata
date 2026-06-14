"""
01_download_and_preprocess.py

Downloads and preprocesses intervention datasets with confirmed public beta matrices.

Datasets:
  GSE50660  — smoking, blood, beta values in series matrix (Tsaprouni 2014)
  GSE133588 — chemotherapy, blood, log2 matrix in suppl file (Sehl 2020)

Datasets pending data access or author response:
  CALERIE   — caloric restriction (Aging Research Biobank application)
  Ronn 2013 — exercise, adipose (author data request)
  Lindholm  — exercise, muscle (author data request)

Once additional data is received, add processors below and update config.
"""

import sys
import re
import gzip
import urllib.request
import numpy as np
import pandas as pd
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

CROSSREACTIVE_URL  = (
    "https://raw.githubusercontent.com/sirselim/illumina450k_filtering/"
    "master/EPIC/13059_2016_1066_MOESM1_ESM.csv"
)
CROSSREACTIVE_DEST = RAW_DIR / "cross_reactive_probes.csv"


# ── Utilities ─────────────────────────────────────────────────────────────────

def load_exclusion_list():
    if not CROSSREACTIVE_DEST.exists():
        print("  downloading cross-reactive probe list ...")
        urllib.request.urlretrieve(CROSSREACTIVE_URL, CROSSREACTIVE_DEST)
    cr  = pd.read_csv(CROSSREACTIVE_DEST, header=0)
    exc = set(cr.iloc[:, 0].astype(str).tolist())
    print(f"  exclusion list: {len(exc):,} probes")
    return exc


def load_common_cpgs():
    if not COMMON_CPGS.exists():
        raise FileNotFoundError(
            f"{COMMON_CPGS} not found. "
            "Run parent repo 01_download_and_preprocess.py first."
        )
    with open(COMMON_CPGS) as f:
        cpgs = [l.strip() for l in f if l.strip()]
    print(f"  common CpGs: {len(cpgs):,}")
    return cpgs


def download(url, dest):
    dest = Path(dest)
    if dest.exists():
        print(f"  already downloaded: {dest.name}")
        return
    print(f"  downloading {dest.name} ...")
    urllib.request.urlretrieve(url, dest)
    print(f"  done.")


def qc_and_restrict(beta_df, common_cpgs, excluded, label):
    print(f"  [{label}] raw shape: {beta_df.shape}")
    beta_df = beta_df[[c for c in beta_df.columns if c not in excluded]]
    beta_df = beta_df.apply(pd.to_numeric, errors='coerce').clip(0, 1)
    beta_df = beta_df.loc[:, beta_df.isnull().mean(axis=0) < 0.05]
    mean_b  = beta_df.mean(axis=0)
    beta_df = beta_df.loc[:, (mean_b >= BETA_MIN) & (mean_b <= BETA_MAX)]
    common  = [c for c in common_cpgs if c in beta_df.columns]
    beta_df = beta_df[common]
    print(f"  [{label}] after QC + restriction: {beta_df.shape}")
    return beta_df


def save_h5(beta_df, meta_df, accession):
    out_beta = interv_beta_path(accession)
    out_meta = interv_meta_path(accession)
    with h5py.File(out_beta, 'w') as f:
        f.create_dataset('beta',       data=beta_df.values.astype(np.float32))
        f.create_dataset('sample_ids', data=np.array(beta_df.index.tolist(),   dtype='S64'))
        f.create_dataset('cpg_ids',    data=np.array(beta_df.columns.tolist(), dtype='S16'))
    meta_df.to_csv(out_meta, index=False)
    print(f"  saved: {out_beta.name} ({beta_df.shape[0]} samples x {beta_df.shape[1]:,} CpGs)")


def already_done(accession):
    return interv_beta_path(accession).exists() and interv_meta_path(accession).exists()


# ── GSE50660: smoking, blood ──────────────────────────────────────────────────

def process_GSE50660(common_cpgs, excluded):
    """
    Tsaprouni et al. 2014 — CARDIOGENICS cohort, n=464, whole blood.
    Beta values stored in series matrix VALUE columns.
    Cross-sectional: current vs never smokers.
    """
    print(f"\n{'─'*60}\nGSE50660 — smoking, blood\n{'─'*60}")

    if already_done("GSE50660"):
        print("  already preprocessed, skipping")
        return

    matrix_dest = RAW_DIR / "GSE50660_series_matrix.txt.gz"
    download(INTERVENTION_DATASETS["GSE50660"]["matrix_url"], matrix_dest)

    print("  parsing series matrix (fast path) ...")

    # First pass: collect metadata and locate data table header row
    opener2 = gzip.open if str(matrix_dest).endswith('.gz') else open
    meta_lines  = {}
    skip_rows   = 0
    sample_ids  = []

    with opener2(matrix_dest, 'rt', errors='replace') as fh:
        for i, line in enumerate(fh):
            line = line.rstrip('\n')
            if line.startswith('!'):
                key, _, val = line.partition('\t')
                meta_lines.setdefault(key.lstrip('!'), []).append(val)
            elif line.startswith('"ID_REF"') or line.startswith('ID_REF'):
                parts      = line.split('\t')
                sample_ids = [p.strip('"') for p in parts[1:]]
                skip_rows  = i
                break

    if not sample_ids:
        print("  ERROR: could not find ID_REF header")
        return

    print(f"  {len(sample_ids)} samples, data at line {skip_rows}")
    print("  loading with pandas ...")

    beta_df = pd.read_csv(
        matrix_dest,
        sep='\t',
        skiprows=skip_rows,
        index_col=0,
        na_values=['null', 'NULL', 'NA', ''],
        low_memory=False,
    )

    # Keep only cg* rows and drop end sentinel
    beta_df = beta_df[beta_df.index.astype(str).str.startswith('cg')]
    beta_df.columns = sample_ids[:beta_df.shape[1]]
    beta_df = beta_df.T
    beta_df.index.name = 'sample_id'
    print(f"  parsed {beta_df.shape[1]:,} CpGs x {beta_df.shape[0]} samples")

    # char_meta is a list of full lines, one per characteristics field.
    # Each line contains tab-separated values for all samples.
    # Find the smoking line and extract per-sample values.
    char_lines = meta_lines.get('Sample_characteristics_ch1', [])
    smoking_vals = {}   # sid -> '0', '1', or '2'
    age_vals     = {}
    for line in char_lines:
        parts = [p.strip().strip('"') for p in line.split('\t')]
        if parts and 'smoking' in parts[0].lower():
            for j, sid in enumerate(sample_ids):
                val = parts[j] if j < len(parts) else ''
                # Extract trailing digit: "smoking (...): 2" -> "2"
                m = re.search(r':\s*(\d+)\s*$', val)
                smoking_vals[sid] = m.group(1) if m else 'unknown'
        elif parts and 'age' in parts[0].lower():
            for j, sid in enumerate(sample_ids):
                val = parts[j] if j < len(parts) else ''
                m = re.search(r':\s*(\d+)', val)
                age_vals[sid] = int(m.group(1)) if m else None

    # Map smoking code to group
    # 0 = never (control), 1 = former (excluded), 2 = current (case)
    code_to_group = {'0': 'control', '1': 'former', '2': 'case'}

    rows = []
    for sid in sample_ids:
        code  = smoking_vals.get(sid, 'unknown')
        group = code_to_group.get(code, 'unknown')
        rows.append({
            'sample_id':  sid,
            'group':      group,
            'age':        age_vals.get(sid),
            'timepoint':  'single',
            'subject_id': sid,
            'accession':  'GSE50660',
            'label':      'smoking_blood',
            'sign':       'plus',
            'tissue':     'blood',
            'design':     'cross_sectional',
        })

    meta_df = pd.DataFrame(rows)
    print(f"  cases (current, code=2): {(meta_df['group']=='case').sum()}")
    print(f"  controls (never, code=0): {(meta_df['group']=='control').sum()}")
    print(f"  former (code=1, excluded): {(meta_df['group']=='former').sum()}")

    beta_df = qc_and_restrict(beta_df, common_cpgs, excluded, "GSE50660")
    save_h5(beta_df, meta_df, "GSE50660")


# ── GSE133588: chemotherapy, blood ────────────────────────────────────────────

def process_GSE133588(common_cpgs, excluded):
    """
    Sehl et al. 2020 — breast cancer, blood, pre/post cytotoxic chemotherapy.
    EPIC array, n=48 paired samples.
    Supplementary file: log2 normalised M-values.
    Back-transform: beta = 2^M / (1 + 2^M)
    """
    print(f"\n{'─'*60}\nGSE133588 — chemotherapy, blood\n{'─'*60}")

    if already_done("GSE133588"):
        print("  already preprocessed, skipping")
        return

    config = INTERVENTION_DATASETS["GSE133588"]

    # Download log2 matrix
    beta_dest   = RAW_DIR / "GSE133588_log2_norm.txt.gz"
    matrix_dest = RAW_DIR / "GSE133588_series_matrix.txt.gz"
    download(config["beta_url"],   beta_dest)
    download(config["matrix_url"], matrix_dest)

    # Load log2 matrix and back-transform
    print("  loading log2 matrix ...")
    log2_df = pd.read_csv(beta_dest, sep='\t', index_col=0, low_memory=False)
    print(f"  log2 matrix shape: {log2_df.shape} (CpGs x samples)")

    # Back-transform M-values to beta: beta = 2^M / (1 + 2^M)
    print("  back-transforming M-values to beta ...")
    beta_arr = np.power(2, log2_df.values) / (1 + np.power(2, log2_df.values))
    beta_df  = pd.DataFrame(beta_arr, index=log2_df.index, columns=log2_df.columns).T
    beta_df.index.name = 'sample_id'

    val_range = (beta_df.values.min(), beta_df.values.max())
    print(f"  beta range after transform: [{val_range[0]:.3f}, {val_range[1]:.3f}]")

    # Parse metadata from series matrix
    print("  parsing series matrix for metadata ...")
    opener = gzip.open if str(matrix_dest).endswith('.gz') else open
    meta_lines = {}
    with opener(matrix_dest, 'rt', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!Sample_'):
                key, _, val = line.partition('\t')
                key = key.lstrip('!')
                meta_lines.setdefault(key, []).append(val.strip('"'))
            elif line.startswith('"ID_REF"') or line.startswith('ID_REF'):
                break

    sample_geo    = meta_lines.get('Sample_geo_accession', [])
    sample_titles = meta_lines.get('Sample_title', [])
    sample_chars  = meta_lines.get('Sample_characteristics_ch1', [])

    # GSE133588: samples are labeled pre/post chemo
    # title format typically: "Patient X pre" or "Patient X post"
    rows = []
    for i, sid in enumerate(beta_df.index):
        title = sample_titles[i] if i < len(sample_titles) else sid
        chars = sample_chars[i]  if i < len(sample_chars)  else ''
        combined = (title + ' ' + chars).lower()

        if any(kw in combined for kw in ['post', 'after', 'cycle', 'follow']):
            timepoint = 'post'
        elif any(kw in combined for kw in ['pre', 'before', 'baseline', 'prior']):
            timepoint = 'pre'
        else:
            timepoint = 'unknown'

        # Subject ID from title — extract number
        m = re.search(r'patient\s*(\d+)|subject\s*(\d+)|\bp(\d+)\b', combined)
        subject_id = (m.group(1) or m.group(2) or m.group(3)) if m else sid

        rows.append({
            'sample_id':  sid,
            'timepoint':  timepoint,
            'subject_id': subject_id,
            'group':      'patient',
            'accession':  'GSE133588',
            'label':      'chemotherapy',
            'sign':       'plus',
            'tissue':     'blood',
            'design':     'longitudinal',
        })

    meta_df = pd.DataFrame(rows)
    print(f"  pre: {(meta_df['timepoint']=='pre').sum()}, "
          f"post: {(meta_df['timepoint']=='post').sum()}, "
          f"unknown: {(meta_df['timepoint']=='unknown').sum()}")

    beta_df = qc_and_restrict(beta_df, common_cpgs, excluded, "GSE133588")
    save_h5(beta_df, meta_df, "GSE133588")


# ── Placeholder loader for author-supplied / access-controlled data ───────────

def load_if_available(accession):
    """
    For datasets pending data access (CALERIE, exercise):
    check if beta matrix has been manually placed in interventions dir.
    """
    config   = INTERVENTION_DATASETS.get(accession, {})
    r_beta   = Path(config.get('r_beta', ''))
    r_meta   = Path(config.get('r_meta', ''))

    if r_beta.exists() and r_meta.exists():
        print(f"  found manually supplied data for {accession}, loading ...")
        # These are expected as gzipped TSV (CpGs x samples) from R or author
        beta_df = pd.read_csv(r_beta, sep='\t', index_col=0).T
        meta_df = pd.read_csv(r_meta)
        print(f"  shape: {beta_df.shape}")
        return beta_df, meta_df
    else:
        note = config.get('note', '')
        print(f"  {accession} not yet available ({config.get('source','')})")
        if note:
            print(f"  note: {note}")
        return None, None


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading exclusion probe list ...")
    excluded = load_exclusion_list()

    print("\nLoading common CpG list ...")
    common_cpgs = load_common_cpgs()

    # Datasets with confirmed public access
    process_GSE50660(common_cpgs, excluded)
    process_GSE133588(common_cpgs, excluded)

    # Datasets pending access — load if manually placed
    print(f"\n{'─'*60}\nChecking for pending datasets\n{'─'*60}")
    for accession in ['CALERIE', 'GSE77716', 'GSE56867', 'GSE64930']:
        if accession not in INTERVENTION_DATASETS:
            continue
        config = INTERVENTION_DATASETS[accession]
        if config.get('source') in ('controlled_access', 'r_series_matrix'):
            print(f"\n  {accession} ({config['label']}):")
            beta_df, meta_df = load_if_available(accession)
            if beta_df is not None:
                beta_df = qc_and_restrict(beta_df, common_cpgs, excluded, accession)
                meta_df['accession'] = accession
                meta_df['label']     = config['label']
                meta_df['sign']      = config['sign']
                meta_df['tissue']    = config['tissue']
                meta_df['design']    = config['design']
                save_h5(beta_df, meta_df, accession)

    # Summary
    print(f"\n{'='*60}")
    print("Summary:")
    for accession in INTERVENTION_DATASETS:
        status = "ready" if already_done(accession) else "pending"
        print(f"  {accession}: {status}")


if __name__ == "__main__":
    main()

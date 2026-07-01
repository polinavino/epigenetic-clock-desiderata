"""
02_build_displacement_matrix.py

Computes per-intervention mean displacement vectors and assembles the
signed displacement matrix A used to identify the aging direction v*.

For cross-sectional datasets (smoking): displacement = mean(smokers) - mean(never)
For longitudinal datasets (chemo, CR, exercise): displacement = mean(post - pre)

Each row of A is one intervention-context pair, signed:
  +displacement  if the intervention is in I_PLUS  (age-accelerating)
  -displacement  if the intervention is in I_MINUS (geroprotective)

Each entry is further normalised by per-CpG baseline standard deviation
(the weighted inner product from Eq. 1 in the paper).

Outputs:
  data/interventions/displacement_matrix_A.parquet
  data/interventions/displacement_matrix_metadata.csv
  data/interventions/{ACCESSION}_displacement.parquet   (per-dataset)
"""

import sys
import numpy as np
import pandas as pd
import h5py
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    INTERVENTION_DATASETS, INTERV_DIR,
    GSE40279_BETA, GSE87571_BETA,
    DISPLACEMENT_MATRIX, DISPLACEMENT_METADATA,
    N_CPGS_DISPLACEMENT, MIN_PAIRED_SAMPLES,
    interv_beta_path, interv_meta_path, displacement_path,
)


def load_h5_beta(path):
    """
    Load beta matrix from HDF5, return DataFrame (samples x CpGs).

    Handles two formats:
      - Our format: datasets 'beta', 'sample_ids', 'cpg_ids' (samples x CpGs)
      - Parent repo format: pandas HDF5 with key 'beta' (CpGs x samples)
    """
    with h5py.File(path, 'r') as f:
        keys = list(f.keys())
        has_our_format = 'sample_ids' in keys and 'cpg_ids' in keys

    if has_our_format:
        with h5py.File(path, 'r') as f:
            sample_ids = [s.decode() for s in f['sample_ids'][:]]
            cpg_ids    = [c.decode() for c in f['cpg_ids'][:]]
            beta       = f['beta'][:]
        return pd.DataFrame(beta, index=sample_ids, columns=cpg_ids)
    else:
        # Parent repo pandas HDF5: CpGs x samples -> transpose to samples x CpGs
        df = pd.read_hdf(path, key='beta')
        return df.T


def compute_baseline_std(common_cpgs):
    """
    Compute per-CpG standard deviation across control individuals
    from the cross-sectional aging cohorts (parent repo datasets).
    Used for the weighted inner product normalisation.
    Returns Series indexed by CpG ID.
    """
    print("  computing baseline std from cross-sectional cohorts ...")
    dfs = []
    for path in [GSE40279_BETA, GSE87571_BETA]:
        if path.exists():
            df = load_h5_beta(path)
            common = [c for c in common_cpgs if c in df.columns]
            dfs.append(df[common])

    if not dfs:
        raise FileNotFoundError(
            "Cross-sectional beta matrices not found. "
            "Run parent repo 01_download_and_preprocess.py first."
        )

    combined = pd.concat(dfs, axis=0)[common_cpgs]
    sigma = combined.std(axis=0)
    # Floor at small value to avoid division by zero
    sigma = sigma.clip(lower=1e-4)
    return sigma


def displacement_cross_sectional(beta_df, meta_df, accession):
    """
    Compute displacement for cross-sectional designs.
    Returns mean(case) - mean(control) as a Series indexed by CpG.
    Expects meta_df to have a 'group' column with values 'case'/'control'
    or inferred from sample characteristics.
    """
    # For smoking datasets: infer group from metadata if not explicit
    if 'group' not in meta_df.columns:
        # Heuristic: samples with 'current' in title -> case, 'never' -> control
        meta_df = meta_df.copy()
        meta_df['group'] = 'unknown'

    cases    = meta_df.loc[meta_df['group'] == 'case',    'sample_id'].tolist()
    controls = meta_df.loc[meta_df['group'] == 'control', 'sample_id'].tolist()

    cases    = [s for s in cases    if s in beta_df.index]
    controls = [s for s in controls if s in beta_df.index]

    if len(cases) < MIN_PAIRED_SAMPLES or len(controls) < MIN_PAIRED_SAMPLES:
        print(f"  WARNING {accession}: insufficient samples "
              f"(cases={len(cases)}, controls={len(controls)}), skipping")
        return None

    delta = beta_df.loc[cases].mean(axis=0) - beta_df.loc[controls].mean(axis=0)
    print(f"  {accession}: {len(cases)} cases, {len(controls)} controls, "
          f"mean |delta| = {delta.abs().mean():.4f}")
    return delta


def displacement_longitudinal(beta_df, meta_df, accession):
    """
    Compute displacement for longitudinal paired designs.
    Returns mean(post - pre) as a Series indexed by CpG.
    Expects meta_df to have 'timepoint' column with values 'pre'/'post'
    and a 'subject_id' column for pairing.
    """
    if 'subject_id' not in meta_df.columns:
        # Attempt to infer subject ID from sample ID prefix
        meta_df = meta_df.copy()
        meta_df['subject_id'] = meta_df['sample_id'].str.extract(r'(GSM\d+)')[0]

    pre_meta  = meta_df[meta_df['timepoint'] == 'pre']
    post_meta = meta_df[meta_df['timepoint'] == 'post']

    # Match on subject_id
    common_subjects = set(pre_meta['subject_id']) & set(post_meta['subject_id'])

    if len(common_subjects) < MIN_PAIRED_SAMPLES:
        print(f"  WARNING {accession}: only {len(common_subjects)} paired subjects, "
              f"skipping")
        return None

    pre_samples  = pre_meta.set_index('subject_id').loc[list(common_subjects), 'sample_id']
    post_samples = post_meta.set_index('subject_id').loc[list(common_subjects), 'sample_id']

    pre_valid  = [s for s in pre_samples  if s in beta_df.index]
    post_valid = [s for s in post_samples if s in beta_df.index]

    if len(pre_valid) < MIN_PAIRED_SAMPLES:
        print(f"  WARNING {accession}: insufficient paired samples in beta matrix, skipping")
        return None

    delta = (beta_df.loc[post_valid].values - beta_df.loc[pre_valid].values).mean(axis=0)
    delta = pd.Series(delta, index=beta_df.columns)

    print(f"  {accession}: {len(pre_valid)} pairs, "
          f"mean |delta| = {delta.abs().mean():.4f}")
    return delta


def compute_displacement(accession, config):
    """
    Load preprocessed data and compute displacement vector.
    Returns signed displacement Series (+ for I_PLUS, - for I_MINUS).
    """
    beta_path = interv_beta_path(accession)
    meta_path = interv_meta_path(accession)

    if not beta_path.exists():
        print(f"  {accession}: beta matrix not found, skipping")
        return None, None

    beta_df = load_h5_beta(beta_path)
    meta_df = pd.read_csv(meta_path)
    # sample_id must be str to match the (string) h5 sample index; purely
    # numeric IDs (e.g. GSE89218 "201227") would otherwise be read as int.
    meta_df['sample_id'] = meta_df['sample_id'].astype(str)

    if config['design'] == 'cross_sectional':
        delta = displacement_cross_sectional(beta_df, meta_df, accession)
    elif config['design'] == 'longitudinal':
        delta = displacement_longitudinal(beta_df, meta_df, accession)
    else:
        print(f"  {accession}: unknown design '{config['design']}', skipping")
        return None, None

    if delta is None:
        return None, None

    # Apply sign: + for accelerating, - for geroprotective
    sign = +1 if config['sign'] == 'plus' else -1
    delta_signed = delta * sign

    return delta_signed, {
        'accession': accession,
        'label':     config['label'],
        'sign':      config['sign'],
        'tissue':    config['tissue'],
        'design':    config['design'],
        'n_cpgs':    len(delta),
    }


def select_top_cpgs(displacements, n=None):
    """
    Select CpGs present in all displacement vectors, restricted to
    those with highest variance across displacements.
    Returns list of CpG IDs.
    """
    if n is None:
        n = N_CPGS_DISPLACEMENT

    # Intersection of all displacement CpG sets
    common = None
    for d in displacements:
        cpgs = set(d.index)
        common = cpgs if common is None else common & cpgs

    common = list(common)
    print(f"  CpGs present in all displacement vectors: {len(common):,}")

    # Select top n by variance across displacement vectors
    disp_matrix = pd.DataFrame({i: d[common] for i, d in enumerate(displacements)}).T
    variances = disp_matrix.var(axis=0)
    top_cpgs = variances.nlargest(n).index.tolist()

    print(f"  selected top {len(top_cpgs):,} CpGs by displacement variance")
    return top_cpgs


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Computing per-intervention displacement vectors ...")

    displacements = []
    metadata_rows = []

    for accession, config in INTERVENTION_DATASETS.items():
        out_path = displacement_path(accession)

        if out_path.exists():
            print(f"  {accession}: loading cached displacement")
            delta = pd.read_parquet(out_path).squeeze()
            sign  = +1 if config['sign'] == 'plus' else -1
            delta_signed = delta * sign
            displacements.append(delta_signed)
            metadata_rows.append({
                'accession': accession,
                'label':     config['label'],
                'sign':      config['sign'],
                'tissue':    config['tissue'],
            })
            continue

        delta_signed, meta = compute_displacement(accession, config)

        if delta_signed is not None:
            # Save raw (unsigned) displacement
            delta_unsigned = delta_signed * (+1 if config['sign'] == 'plus' else -1)
            delta_unsigned.to_frame('displacement').to_parquet(out_path)

            displacements.append(delta_signed)
            metadata_rows.append(meta)

    if len(displacements) < 2:
        raise RuntimeError(
            f"Only {len(displacements)} displacement vectors computed. "
            "Need at least 2 to build matrix A. Check dataset downloads."
        )

    print(f"\nComputed {len(displacements)} displacement vectors")

    # Select top CpGs
    print("\nSelecting top CpGs ...")
    top_cpgs = select_top_cpgs(displacements)

    # Compute baseline std for normalisation
    print("\nComputing baseline std for weighted inner product ...")
    sigma = compute_baseline_std(top_cpgs)
    sigma = sigma.reindex(top_cpgs).fillna(1e-4)

    # Assemble matrix A: rows = intervention-context pairs, cols = CpGs
    # Normalise each row by sigma (weighted inner product)
    print("\nAssembling displacement matrix A ...")
    rows = []
    for d in displacements:
        row = d.reindex(top_cpgs).fillna(0)
        row_normalised = row / sigma
        rows.append(row_normalised.values)

    A = pd.DataFrame(
        np.array(rows),
        columns=top_cpgs
    )

    A.to_parquet(DISPLACEMENT_MATRIX)
    print(f"  A shape: {A.shape}")

    meta_df = pd.DataFrame(metadata_rows)
    meta_df.to_csv(DISPLACEMENT_METADATA, index=False)

    print(f"\nSaved:")
    print(f"  {DISPLACEMENT_MATRIX}")
    print(f"  {DISPLACEMENT_METADATA}")


if __name__ == "__main__":
    main()

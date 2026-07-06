"""
process_gse42861_stream.py

Smoking cohort #3 (GSE42861, Liu et al. 2013 RA EWAS; 450k whole blood/PBL, n=689,
smoking status current/ex/never/occasional). Accelerator (I_PLUS), cross-sectional
(current vs never).

The processed matrix is a 2.7 GB gz (~485k CpGs x 690 samples). On this
memory-constrained machine, building the full samples x CpGs matrix thrashes for
hours, so we do NOT store a beta h5. Instead we STREAM the file once (O(n_cpg)
memory) and write the displacement vector directly to
data/interventions/GSE42861_displacement.parquet, which 02_build_displacement_matrix.py
loads from cache (bypassing the h5 entirely).

Displacement = mean(beta | current) - mean(beta | never), per CpG, restricted to the
common 450k/EPIC CpG list, with the same saturation filter used elsewhere
(all-sample mean beta in [BETA_MIN, BETA_MAX]). Methodologically identical to the
cross-sectional path in 02, just computed by streaming.

Note: because there is no h5, GSE42861 participates in 02/03/05 + geometry (SVD,
angles, classification) but NOT in cell_composition_analysis.py (which needs the
full matrix). It is a smoking-reproducibility check; the 8-dataset cell-composition
result is unaffected. RA disease status is a confound (kept, not adjusted).
"""

import sys
import gzip
import re
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_intervention import (
    RAW_DIR, COMMON_CPGS, BETA_MIN, BETA_MAX, interv_meta_path, displacement_path,
)

MATRIX = RAW_DIR / "GSE42861_processed_methylation_matrix.txt.gz"
SERIES = RAW_DIR / "GSE42861_series_matrix.txt.gz"


def sentrix_to_group():
    """Map sentrix id -> 'case' (current) / 'control' (never) from series matrix."""
    suppl_grn, smoking = None, None
    with gzip.open(SERIES, 'rt') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!Sample_supplementary_file') and suppl_grn is None and 'Grn' in line:
                suppl_grn = [x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1') and 'smoking status' in line:
                smoking = [x.strip().strip('"').split(': ', 1)[-1]
                           for x in line.rstrip('\n').split('\t')[1:]]
    assert suppl_grn and smoking and len(suppl_grn) == len(smoking), "series parse failed"
    mp = {}
    for url, sm in zip(suppl_grn, smoking):
        m = re.search(r'_(\d+_R\d+C\d+)_Grn', url)
        if not m:
            continue
        if sm == 'current':
            mp[m.group(1)] = 'case'
        elif sm == 'never':
            mp[m.group(1)] = 'control'
    return mp


def main():
    common = set(l.strip() for l in open(COMMON_CPGS) if l.strip())
    s2g = sentrix_to_group()
    n_case = sum(v == 'case' for v in s2g.values())
    n_ctrl = sum(v == 'control' for v in s2g.values())
    print(f"mapped {n_case} current (case) / {n_ctrl} never (control) sentrix ids")

    with gzip.open(MATRIX, 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        cols = [c.strip().strip('"') for c in header[1:]]
        case_mask = np.array([s2g.get(c) == 'case' for c in cols])
        ctrl_mask = np.array([s2g.get(c) == 'control' for c in cols])
        print(f"matrix has {len(cols)} sample columns; "
              f"{case_mask.sum()} case / {ctrl_mask.sum()} control matched")

        cpgs, case_sum, case_n, ctrl_sum, ctrl_n, all_sum, all_n = [], [], [], [], [], [], []
        seen = 0
        for line in f:
            tab = line.find('\t')
            cpg = line[:tab].strip().strip('"')
            if cpg not in common:
                continue
            vals = np.fromstring(line[tab + 1:], sep='\t')
            if vals.size != len(cols):
                continue
            c = vals[case_mask]; k = vals[ctrl_mask]
            cpgs.append(cpg)
            case_sum.append(np.nansum(c)); case_n.append(np.count_nonzero(~np.isnan(c)))
            ctrl_sum.append(np.nansum(k)); ctrl_n.append(np.count_nonzero(~np.isnan(k)))
            all_sum.append(np.nansum(vals)); all_n.append(np.count_nonzero(~np.isnan(vals)))
            seen += 1
            if seen % 100000 == 0:
                print(f"  ...{seen:,} common CpGs streamed")

    case_mean = np.array(case_sum) / np.maximum(np.array(case_n), 1)
    ctrl_mean = np.array(ctrl_sum) / np.maximum(np.array(ctrl_n), 1)
    all_mean = np.array(all_sum) / np.maximum(np.array(all_n), 1)
    disp = pd.Series(case_mean - ctrl_mean, index=cpgs)

    # saturation filter (match qc_and_save): all-sample mean in [BETA_MIN, BETA_MAX]
    keep = (all_mean >= BETA_MIN) & (all_mean <= BETA_MAX)
    disp = disp[keep]
    print(f"streamed {len(cpgs):,} common CpGs; kept {int(keep.sum()):,} after saturation filter")

    # AHRR sanity check
    if 'cg05575921' in disp.index:
        print(f"  AHRR cg05575921 displacement (current-never) = {disp['cg05575921']:+.4f} "
              f"(expect strongly negative = hypomethylation in smokers)")
    print(f"  most negative CpG: {disp.idxmin()} ({disp.min():+.4f}); "
          f"rank of AHRR by most-negative: "
          f"{int((disp < disp.get('cg05575921', 0)).sum()) + 1 if 'cg05575921' in disp.index else 'NA'}")

    out = displacement_path("GSE42861")
    disp.to_frame("displacement").to_parquet(out)
    # minimal metadata for provenance
    meta = pd.DataFrame({"sentrix": list(s2g), "group": list(s2g.values())})
    meta.rename(columns={"sentrix": "sample_id"}).to_csv(interv_meta_path("GSE42861"), index=False)
    print(f"saved displacement: {out.name} ({len(disp):,} CpGs)")


if __name__ == "__main__":
    main()

# Epigenetic Clock Desiderata

**Paper:** Towards a Formal Definition of Biological Age: Empirical Characterization of Clock Disagreement and Proposed Desiderata  
**Author:** Polina Vinogradova  
**Preprint:** (forthcoming)

---

## Plain-language explanation

(to be written after empirical results)

---

## Reproducing the results

### Requirements
pip install -r requirements.txt

### Data

Downloaded automatically by script 1. Requires ~2GB disk space and ~8GB RAM.

### Scripts

| Script | Description |
|--------|-------------|
| `scripts/01_download_and_preprocess.py` | Download GEO datasets, QC, output clean beta matrix |
| `scripts/02_compute_weights_and_curve.py` | Compute w_i weights, select top CpGs, fit principal curve |
| `scripts/03_compute_clocks_and_ranks.py` | Run all clocks, compute tau/r decomposition, rank analysis |
| `scripts/04_calerie_analysis.py` | Intervention type consistency test on CALERIE data |

### Replication steps

```bash
python scripts/01_download_and_preprocess.py
python scripts/02_compute_weights_and_curve.py
python scripts/03_compute_clocks_and_ranks.py
python scripts/04_calerie_analysis.py
```

---

## Repository structure
.
├── README.md
├── requirements.txt
├── scripts/
│   ├── 01_download_and_preprocess.py
│   ├── 02_compute_weights_and_curve.py
│   ├── 03_compute_clocks_and_ranks.py
│   └── 04_calerie_analysis.py
├── data/                        # populated by scripts, not tracked by git
└── paper/
├── main.tex
├── references.bib
├── figures/                 # populated by scripts
└── sections/
├── abstract.tex
├── introduction.tex
├── related_work.tex
├── methods.tex
├── results.tex
├── desiderata.tex
├── discussion.tex
└── conclusion.tex

## Citation
Vinogradova, P. (2026). Towards a Formal Definition of Biological Age:
Empirical Characterization of Clock Disagreement and Proposed Desiderata.

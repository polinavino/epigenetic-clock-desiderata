# A Formal Classification of Interventions on Biological Age

**Paper:** A Formal Classification of Interventions on Biological Age: An Invariance-Based Framework for Epigenetic Aging
**Author:** Polina Vinogradova
**Status:** In preparation

## Relationship to companion paper

This is the second paper in this repository. The first paper (`paper/`) establishes
that epigenetic clocks disagree systematically on intervention response and proposes
desiderata for a valid biological age measure. This paper provides the formal
resolution: a clock-independent definition of what an intervention does to biological
age, derived from a geometric framework in methylation space.

The two papers are intended to be read together. The clock disagreement on CALERIE
documented in the first paper is the primary motivating finding for the framework
developed here.

Shared data infrastructure (downloaded GEO datasets, preprocessing) lives in the
top-level `data/` directory of this repository and is populated by script 01.

## Structure

```
intervention-framework/
├── paper/
│   ├── main.tex
│   ├── references.bib
│   └── sections/
│       ├── abstract.tex
│       ├── introduction.tex
│       ├── related_work.tex
│       ├── framework.tex
│       ├── identification.tex
│       ├── results.tex
│       ├── discussion.tex
│       └── conclusion.tex
├── scripts/
│   ├── 01_download_and_preprocess.py
│   ├── 02_build_displacement_matrix.py
│   ├── 03_svd_and_aging_direction.py
│   ├── 04_fit_principal_curve.py
│   ├── 05_classify_interventions.py
│   └── 06_validation.py
├── results/        (populated by scripts, not tracked)
└── requirements.txt
```

## Pipeline

| Script | Description |
|--------|-------------|
| 01 | GEO download, QC, normalization, batch correction |
| 02 | Per-intervention displacement vectors, assemble matrix A |
| 03 | SVD of A, aging direction v*, dimensionality and curvature tests |
| 04 | Principal curve fitting, arc-length computation |
| 05 | Intervention classification, CALERIE position/rate decomposition |
| 06 | Mortality validation, concordance with independent clocks |

## Citation

Vinogradova, P. (2026). A Formal Classification of Interventions on Biological Age:
An Invariance-Based Framework for Epigenetic Aging. In preparation.

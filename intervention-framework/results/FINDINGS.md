# Findings (current) — 8 interventions, cell-composition control, GTEx cross-tissue

Supersedes earlier notes. Pipeline runs on **8 blood interventions** (balanced
4 accelerating / 4 geroprotective), with a leukocyte-composition control and a
GTEx multi-tissue cross-tissue test.

| Intervention | Accession | Sign | Design | Tissue |
|--------------|-----------|------|--------|--------|
| Smoking | GSE50660 | accel (+) | current vs never | whole blood |
| Smoking (replication) | GSE53045 | accel (+) | smoker vs control | PBMC |
| Chemo + radiotherapy | GSE140038 | accel (+) | post vs pre | whole blood |
| Chronic stress / PTSD | GSE89218 | accel (+) | PTSD+ vs PTSD− | whole blood |
| Bariatric weight loss | GSE272137 | gero (−) | paired w0/w52 | whole blood |
| Behavioural weight loss | GSE240184 | gero (−) | paired BL/3m | whole blood |
| Exercise (8 wk) | GSE328810 | gero (−) | paired Before/After | whole blood |
| Exercise (children, 20 wk) | GSE193730 | gero (−) | paired Baseline/T1 | whole blood |

Shared CpG support: 30,214 CpGs present in all eight displacement vectors.

---

## 1. No dominant shared aging axis — robust to n and to cell adjustment

Signed row-normalised SVD (raw displacements): ρ = λ₁/λ₂ = **1.30**
(σ-weighted A in `03`: 1.37). After EpiDISH leukocyte adjustment: **ρ = 1.12**,
even flatter. Growing 2 → 5 → 6 → 8 interventions never concentrates the signal
onto one axis. The framework declines to certify a v\*.

## 2. Near-orthogonality is real, not a cell-mix artifact

Raw pairwise angles cluster near 90°; after removing the leukocyte-fraction
component (within-dataset residualisation on EpiDISH fractions) every pair lands
in **75–99°** and ρ drops to 1.12. So "no shared axis" is not cell mix.

Fraction of each raw displacement that is *not* cell composition, cos(raw, adj):
smoking-blood 0.93, smoking-PBMC 0.93, exercise 0.97, behavioural-WL 0.91 (mostly
biology) vs **chemo 0.36, bariatric 0.57** (substantially leukocyte shifts —
chemo's raw vector is largely lymphodepletion, CD4T −0.03 / Neutro +0.04).
**Retraction of the earlier note:** the dramatic raw "smoking vs chemo 114.6°,
opposite accelerators" was mostly that chemo artifact — cell-adjusted they are
87° (orthogonal, not anti-parallel).

## 3. Accelerators do not share a direction (answers "is it treatment-specific?")

The new non-treatment accelerator, **PTSD (GSE89218), is orthogonal to
everything** — 87–92° to all seven other interventions (and ~90° to v\*). So the
earlier smoking↔chemo split is not just chemo being cytotoxic: even a
non-treatment stress exposure fails to align with smoking or chemo. The four
"accelerators" occupy four different directions in methylation space. There is no
common acceleration axis.

## 4. Smoking reproduces across cohort and tissue (positive control)

The two independent smoking datasets (GSE50660 whole blood, GSE53045 PBMC) are the
**most-aligned pair** both raw (78.6°) and cell-adjusted (82.8°), and project the
same way onto v\* (cos −0.67, −0.60). The geometry reflects biology, not
single-cohort noise. (Caveat: even same-exposure cohorts are only ~80° aligned.)

## 5. Classification against the provisional v\* (3/8 match sign)

| Intervention | proj⟨d,v\*⟩ | classified | expected | match |
|---|---:|---|---|:--:|
| chemo/radio | +1.03 | accelerating | accel | ✓ |
| bariatric WL | −1.10 | geroprotective | gero | ✓ |
| behavioural WL | −0.21 | geroprotective | gero | ✓ |
| smoking (blood) | −0.77 | geroprotective | accel | ✗ |
| smoking (PBMC) | −1.34 | geroprotective | accel | ✗ |
| PTSD | −0.10 | ~perpendicular | accel | ✗ |
| exercise (adult) | −0.10 | clock-gaming | gero | ✗ |
| exercise (children) | +0.30 | accelerating | gero | ✗ |

v\* is defined by the chemo + weight-loss axis; the other six interventions sit
off it. The mismatches are the finding, not a failure — a single-axis label is
not meaningful when the interventions don't share an axis.

## 6. GTEx cross-tissue test (context-invariance) — v\* points *away* from aging

Run locally by streaming just the v\* CpGs out of the 6 GB GTEx methylation matrix
(GSE213478, 9 tissues, 987 samples, EPIC). **Age brackets ARE in the GEO metadata**
(20-29 … 70-79; earlier "no age" was a truncated-download error). Per-tissue aging
direction = per-CpG correlation of beta with age, over the v\* CpGs.

- **Aging is itself largely tissue-specific.** Pairwise angles between tissue
  age-directions are mostly 60–100°, ρ = **1.25** — no dominant cross-tissue aging
  axis either, mirroring the intervention result. The exception is a cluster:
  **whole blood and lung age-directions are only 37° apart**, with colon / kidney /
  prostate somewhat aligned to them.
- **The intervention v\* is anti-aligned with real aging.** cos(cross-tissue
  consensus aging direction, v\*) = **−0.52**; cos(**whole-blood** aging direction,
  v\*) = **−0.72**; lung −0.67. The direction interventions collectively define
  points *opposite* to how methylation actually changes with chronological age in
  blood. This is the deepest form of the dissociation the framework is built
  around: intervention effects and chronological aging are not the same axis, and
  can be anti-correlated.

*(Caveat: the age-direction test is computed over the v\*-selected CpGs, which are
biased toward intervention-displacement variance; an unbiased CpG panel would
strengthen it. Robust across v\* versions: blood cos was −0.75 / −0.75 / −0.72 for
the 6/7/8-dataset v\*.)*

---

## Status of the publishability vulnerabilities

- **V1 (cell composition): addressed.** Orthogonality survives EpiDISH adjustment
  (ρ 1.30 → 1.12). Chemo/bariatric partly cell-driven; 115° sub-claim retracted.
- **V2 (n of interventions): addressed** 2 → 8; ρ stays flat as n grows.
- **V3 (accelerator diversity): addressed.** Added a second smoking cohort
  (reproducibility) and PTSD, a *non-treatment* accelerator — which is orthogonal
  to smoking and chemo, showing the accelerator disagreement is not a treatment
  artifact.
- **New:** the GTEx cross-tissue test now runs and shows v\* is anti-aligned with
  blood aging — a strong, independent statement of intervention/aging dissociation.

## Still open
- CALERIE (caloric restriction) — the motivating dataset, still requested.
- An unbiased (non-v\*-selected) CpG panel for the GTEx age-direction test.
- Mortality validation; DamAge/AdaptAge alignment.

## Artifacts (`data/interventions/`, `results/`)
`geometry_report.csv`, `cell_fractions_{ACC}.csv`, `cell_adjustment_summary.txt`,
`intervention_classification.csv`, `gtex_cross_tissue_report.txt`, `GTEX_metadata.csv`;
figures `displacement_angle_matrix.pdf`, `svd_spectrum.pdf`, `gtex_cross_tissue.pdf`.
Repro: `02 → 03 → 04 → 05`, then `exploratory_geometry.py`,
`cell_composition_analysis.py`, `gtex_cross_tissue.py`
(after `/tmp/gtex_extract.py` streams the v\* CpGs out of the GTEx matrix).

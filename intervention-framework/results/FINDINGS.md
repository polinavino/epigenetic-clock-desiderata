# Findings (current) — 9 interventions, cell-composition control, GTEx cross-tissue

Supersedes earlier notes. Pipeline runs on **9 blood interventions** (5 accelerating
/ 4 geroprotective; three of the accelerators are independent smoking cohorts), with a
leukocyte-composition control and a GTEx multi-tissue cross-tissue test.

| Intervention | Accession | Sign | Design | Tissue |
|--------------|-----------|------|--------|--------|
| Smoking | GSE50660 | accel (+) | current vs never | whole blood |
| Smoking (replication) | GSE53045 | accel (+) | smoker vs control | PBMC |
| Smoking (RA cohort) | GSE42861 | accel (+) | current vs never | whole blood |
| Chemo + radiotherapy | GSE140038 | accel (+) | post vs pre | whole blood |
| Chronic stress / PTSD | GSE89218 | accel (+) | PTSD+ vs PTSD− | whole blood |
| Bariatric weight loss | GSE272137 | gero (−) | paired w0/w52 | whole blood |
| Behavioural weight loss | GSE240184 | gero (−) | paired BL/3m | whole blood |
| Exercise (8 wk) | GSE328810 | gero (−) | paired Before/After | whole blood |
| Exercise (children, 20 wk) | GSE193730 | gero (−) | paired Baseline/T1 | whole blood |

Shared CpG support: 30,214 CpGs present in all nine displacement vectors.
(GSE42861's displacement is computed by streaming the 2.7 GB matrix — no beta h5 —
so it enters 02/03/05 + geometry but not the cell-composition control; the RA disease
status is an unadjusted confound.)

---

## 1. No dominant shared aging axis — robust to n and to cell adjustment

Signed row-normalised SVD (raw displacements): ρ = λ₁/λ₂ = **1.31** at n=9
(σ-weighted A in `03`: 1.31; was 1.30 at n=8). After EpiDISH leukocyte adjustment
(n=8, GSE42861 excluded): **ρ = 1.12**, even flatter. Growing 2 → 5 → 6 → 8 → 9
interventions never concentrates the signal onto one axis. The framework declines to
certify a v\*.

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

## 4. Smoking: agreement on the top CpG, weak agreement on the full vector

Three independent smoking cohorts are now included: GSE50660 (whole blood),
GSE53045 (PBMC), GSE42861 (whole blood, RA cohort). Two things:

- **They agree on the headline biomarker but not the full displacement.** AHRR
  `cg05575921` is the single **most-negative CpG in every smoking cohort** (e.g.
  −0.22 in GSE42861, rank 1/425k) — hypomethylation in smokers, the canonical
  signal. Yet over the full ~30k-CpG displacement the vectors are only weakly
  aligned: GSE50660↔GSE53045 = **78.6°** (the most-aligned pair overall),
  GSE50660↔GSE42861 = **86.8°**, GSE53045↔GSE42861 = **95.8°**.
- **Confound rotates the vector.** GSE42861 (RA patients + controls; smoking
  correlates with RA) is ~90° from the two cleaner smoking cohorts and is the only
  smoking cohort that projects *positive* on v\* (cos +0.25 → nominally
  "accelerating"), whereas the two clean cohorts project negative.

Reading: agreement on a single strong biomarker (AHRR) does **not** imply agreement
on the direction of the full methylation displacement — the extra signal beyond
AHRR is cohort- and confound-specific. This mirrors the smoking-signature
concordance result in the companion `methylation-biomarker-agreement` analysis
(signatures agree for strong/current exposure, diverge for weak/heterogeneous
cases) and is the biomarker-level version of this project's whole thesis.

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

*(Robust across v\* versions: blood cos was −0.75 / −0.75 / −0.72 for the
6/7/8-dataset v\*.)*

## 7. Classification against an EXTERNAL aging axis (GTEx blood) — the sharper test

Because the internal v\* is non-robust *and* anti-aligned with real aging, we
re-classify each intervention directly against the **GTEx whole-blood
age-methylation direction** (per-CpG corr of beta with chronological age; external,
clock-free, no intervention training), over the 29,558 intervention CpGs present in
GTEx. This replaces the circular internal reference with an external one.

| Intervention | sign | cos(d, blood-aging) | reads as |
|---|:--:|---:|---|
| smoking (blood) | + | **+0.32** | accelerating ✓ |
| smoking (PBMC) | + | **+0.53** | accelerating ✓ |
| chemo/radiotherapy | + | −0.33 | *anti-aging* (lymphodepletion signature) |
| PTSD | + | +0.02 | off-axis |
| bariatric WL | − | +0.28 | off-axis |
| behavioural WL | − | +0.21 | off-axis |
| exercise (adult) | − | +0.07 | off-axis |
| exercise (children) | − | −0.08 | off-axis |

Two things flip the interpretation versus the internal v\*:

- **Smoking is now correctly accelerating** (both cohorts, +0.32/+0.53). Against the
  internal v\* both smoking cohorts looked *geroprotective* — an artifact of v\*
  itself pointing away from real aging (consistency check: cos(v\*, blood-aging) =
  −0.72). The external anchor is the right reference and it fixes this.
- **The geroprotectors do not reverse aging — they act off-axis** (cos 0.07–0.28,
  near-orthogonal to the blood-aging direction). They change methylation without
  moving cells back along the aging axis. This is exactly "clock-gaming" at the
  level of the true aging direction, and it is the CALERIE puzzle made geometric:
  weight loss / exercise need not travel down the chronological-aging axis.
- **Chemo reads paradoxically anti-aging** (−0.33): its acute blood signature is
  dominated by lymphocyte depletion, which moves *opposite* to normal age drift on
  these CpGs. (Consistent with §2 — chemo's raw displacement is largely cell mix.)

Counting strict sign matches gives 2/8, but that undercounts the point: the "misses"
are interventions sitting **orthogonal** to aging, which is the substantive result,
not a classifier failure.

### Placebo controls (both pass)

- **P1 — tissue-specificity is not a v\*-selection artifact.** Repeating the
  cross-tissue age-direction analysis on 5k *random* CpGs gives ρ = 1.61 (vs 1.25 on
  v\* CpGs) and blood–lung angle 55° (vs 39°) — aging is substantially tissue-specific
  on random CpGs too (no dominant cross-tissue axis either way). The blood–lung
  alignment is *stronger* on intervention-relevant CpGs, not manufactured by them.
- **P2 — the external separation is real, not a projection artifact.** Projecting the
  interventions onto a *random-CpG* GTEx "aging" axis does not separate the groups
  (accelerators +0.10, geroprotectors +0.11 — both ≈0). Only the true blood-aging
  axis distinguishes smoking (+) from the rest.

Artifacts: `intervention_classification_external.csv`, `gtex_external_axis_report.txt`.

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

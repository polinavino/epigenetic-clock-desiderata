# Findings — five-intervention run

First run of the full pipeline with **5 interventions** spanning both directions
(previously only 2). All whole blood, 450K/EPIC, displacement computed
within-dataset (case−control or post−pre) so batch effects cancel.

| Intervention | Accession | Sign | Design |
|--------------|-----------|------|--------|
| Smoking | GSE50660 | accel (+) | cross-sectional (current vs never) |
| Chemo + radiotherapy | GSE140038 | accel (+) | post vs pre (unpaired) |
| Bariatric weight loss | GSE272137 | gero (−) | paired w0/w52 |
| Behavioural weight loss | GSE240184 | gero (−) | paired BL/3m |
| Exercise (8 wk) | GSE328810 | gero (−) | paired Before/After |

Shared CpG support: **41,626** CpGs present in all five displacement vectors;
top 5,000 by displacement variance used for SVD.

## 1. No dominant shared aging axis emerges (still)

Signed, row-normalised SVD spectrum: `[1.354, 1.041, 0.983, 0.835, 0.646]`

> **ρ = λ₁/λ₂ = 1.30** — essentially flat. One-dimensionality is **not**
> supported. Going from 2 → 5 interventions did not concentrate the signal onto
> a single axis (ρ was 1.65 with two near-orthogonal vectors). The framework
> again **declines to manufacture a v\***; what we report below uses the leading
> singular vector as a *provisional* v\*, not a validated aging direction.

## 2. The interventions are nearly mutually orthogonal

Pairwise angles between raw (unsigned) displacement vectors (degrees):

|                | smoking | bariatric | exercise | beh.WL | chemo |
|----------------|--------:|----------:|---------:|-------:|------:|
| smoking        |    0.0  |   80.9    |   88.0   |  81.8  | 114.6 |
| bariatric      |         |    0.0    |   86.6   |  87.4  |  97.1 |
| exercise       |         |           |    0.0   |  91.5  |  88.0 |
| beh.weightloss |         |           |          |   0.0  | 108.9 |
| chemo          |         |           |          |        |   0.0 |

Almost every pair is within a few degrees of 90° (orthogonal). Two consequences
worth highlighting:

- **The two accelerators disagree.** Smoking and chemo/radiotherapy are
  **114.6° apart** — they push blood methylation into nearly opposite
  half-spaces despite both being "age-accelerating" exposures. There is no
  common "acceleration direction."
- **The geroprotectors don't align either** (bariatric/exercise/behavioural
  weight loss are 87–91° apart). Even three nominally similar "lose weight / move
  more" interventions do not share a direction in methylation space.

The only mild structure is that chemo sits obtuse to both weight-loss
interventions (97°, 109°) — the one place the expected accel-vs-gero opposition
appears, and weakly.

## 3. Forcing a v\* misclassifies smoking and exercise

Projecting each raw displacement onto the provisional v\* (3/5 match expected sign):

| Intervention | cos(d, v\*) | ‖d‖ | tang. frac | classified | expected | match |
|--------------|-----------:|----:|-----------:|------------|----------|:-----:|
| chemo/radio | +0.81 | 1.57 | 0.81 | accelerating | accelerating | ✓ |
| bariatric WL | −0.32 | 2.48 | 0.32 | geroprotective | geroprotective | ✓ |
| behavioural WL | −0.64 | 0.32 | 0.64 | geroprotective | geroprotective | ✓ |
| **smoking** | **−0.72** | 1.19 | 0.72 | geroprotective | accelerating | ✗ |
| **exercise** | +0.20 | 1.84 | 0.20 | accelerating | geroprotective | ✗ |

The provisional v\* is effectively defined by the chemo + weight-loss axis;
smoking and exercise lie off it (smoking actually **anti-aligned**, cos −0.72).
Against the chronological-age principal curve the earlier smoking cos was −0.34;
against this intervention-consensus v\* it is −0.72 — i.e. smoking is *even more*
off the interventions' shared direction than off the age curve. This is the
sharpest version yet of the clock-gaming / off-axis observation: smoking has a
large, aging-clock-loaded displacement that does not point the way the other
interventions move.

## Interpretation

With 5 interventions the central empirical claim of the two-dataset draft is
**reinforced, not overturned**: intervention-induced methylation displacements in
blood are largely orthogonal, and no robust low-dimensional aging axis v\* falls
out of them. The framework's value here is diagnostic — it quantifies the absence
of a shared axis rather than papering over it. Two readings remain open:

1. **Substantive** — "intervention" is genuinely heterogeneous; accelerators and
   geroprotectors do not act along one axis, so any single-clock adjudication is
   ill-posed (the CALERIE puzzle generalises).
2. **Methodological** — blood cell-composition shifts, tissue/platform, and
   per-study technical structure dominate the displacement and bury a real but
   smaller shared component. Cell-composition adjustment and the cross-tissue
   (GTEx) test are the next levers.

CALERIE remains the key missing geroprotector; a second *non-treatment*
accelerator (e.g. a clean smoking-free exposure) would test whether the
smoking/chemo 115° split is biology or treatment-specific confounding.

## Artifacts
- `data/interventions/geometry_report.csv` — angles + v\* projections
- `data/interventions/intervention_classification.csv` — classification table
- `results/displacement_angle_matrix.pdf` — angle heatmap
- `results/svd_spectrum.pdf` — singular-value spectrum
- `results/classification_figure.pdf` — projection bar chart

Reproduce: `02 → 03 → 04 → 05`, then `exploratory_geometry.py`.

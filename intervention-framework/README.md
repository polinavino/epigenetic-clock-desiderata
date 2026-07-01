# A Formal Classification of Interventions on Biological Age

**An invariance-based framework for epigenetic aging**

Author: Polina Vinogradova
Status: Work in progress — framework and pipeline implemented; empirical validation
in progress pending additional intervention datasets.

This is the second project in this repository. The first (`../paper/`) establishes that
epigenetic clocks disagree systematically on intervention response and proposes desiderata
for a valid biological age measure. This project provides the formal resolution: a
clock-independent definition of what an intervention does to biological age.

---

## 1. The problem

The word "intervention" is load-bearing throughout the epigenetic aging literature but is
used without a formal definition, and almost always carries an unstated assumption that an
intervention is geroprotective. This causes two concrete problems:

1. **Circularity.** Epigenetic clocks are validated by their response to interventions, and
   interventions are evaluated by their effect on clock readings. Each is used to justify the
   other. The circularity has been noted in the literature (e.g. Moqri et al. 2023, Cell,
   Figure 4) but not resolved.

2. **No way to separate genuine effects from artefacts.** If an intervention reduces a clock
   reading, current methods cannot tell whether it slowed the aging process or merely perturbed
   the specific CpG sites the clock happens to read. The CALERIE trial is the motivating case:
   caloric restriction slows DunedinPACE but does not move Horvath age, and neither clock-based
   framework can adjudicate which is correct.

---

## 2. The core idea

Represent a sample's epigenetic state as a point in methylation space (a vector of CpG beta
values). Aging is hypothesised to follow a trajectory through this space. The key move is to
identify the **aging direction** not by correlation with chronological age, and not by any
clock's training objective, but by **what is invariant across interventions and contexts**:

> The aging direction is the direction in methylation space that is consistently induced by
> interventions known to accelerate biological deterioration and consistently opposed by
> interventions known to slow it, across independent biological contexts (tissues, cohorts).

This grounds the framework in mortality and healthspan evidence (why we label smoking as
accelerating and caloric restriction as protective) rather than in clock readings. The
residual circularity — that biological age is ultimately anchored to health outcomes — is
unavoidable and appropriate; what is eliminated is clock-dependence.

From the trajectory, every intervention can be classified by **projecting its displacement
vector onto the aging direction**:

- **Geroprotective** — displaces a sample toward younger positions on the trajectory.
- **Age-accelerating** — displaces toward older positions.
- **Clock-gaming** — produces a large methylation displacement that is mostly *orthogonal* to
  the trajectory; it changes clock readings without moving the sample along the aging axis.
- **Curve-deforming** — alters the shape of the trajectory itself (e.g. cellular reprogramming);
  noted as a category requiring an extended framework, not treated here.

Biological age decomposes into two independent quantities: **position** on the trajectory (what
static clocks like Horvath estimate) and **rate** of traversal (what pace-of-aging measures like
DunedinPACE estimate). An intervention can affect either or both. This decomposition is the
framework's proposed resolution of the CALERIE puzzle: caloric restriction may affect rate
without affecting position.

---

## 3. How the aging trajectory is currently computed

Two distinct constructions appear in this work, and it is important not to confuse them.

### 3a. The chronological-age trajectory (from the companion paper)

Computed in `../scripts/02_compute_weights_and_curve.py`:

1. Two cross-sectional cohorts (GSE40279 + GSE87571), 1,385 individuals aged 14–101, all with
   known chronological age.
2. For each CpG, an age-informativeness weight = R² (correlation with chronological age) ×
   variance. Top 200 CpGs selected.
3. A principal curve is fit through the 1,385 samples in this 200-CpG space. Because the CpGs
   were selected for age-correlation, the curve runs from young-typical to old-typical
   methylation states.
4. Arc-length position along the curve = biological age coordinate (correlates with
   chronological age at r = 0.84).

This is a **cross-sectional, correlational** construct. It is *not* the intervention-invariant
direction the framework calls for. It is used in current analyses only as a fallback reference,
and any comparison against it is labelled accordingly.

### 3b. The intervention-invariant direction v* (the framework's actual target)

Computed in the intervention pipeline (`03_svd_and_aging_direction.py`):

1. For each intervention, compute a displacement vector (case − control for cross-sectional
   designs; post − pre for longitudinal designs).
2. Sign each displacement by its known direction (+ for accelerators, − for geroprotectors),
   normalise each to unit length so direction rather than effect magnitude drives the estimate,
   and assemble into a matrix A.
3. v* = first right singular vector of A; the singular value ratio ρ = λ₁/λ₂ tests whether a
   single shared direction dominates (one-dimensionality of aging).
4. Orient v* toward the accelerating direction.

**v* is the scientifically meaningful object. It currently cannot be robustly estimated because
we have only two interventions (see Section 5).**

---

## 4. Pipeline

| Script | Purpose |
|--------|---------|
| `00_process_idats.R` | Process GEO IDAT/series-matrix datasets via R/GEOquery (for datasets without public beta matrices) |
| `01_download_and_preprocess.py` | Download and QC intervention datasets with public beta matrices |
| `process_gse272137.py` | Dataset-specific processor for the bariatric weight-loss dataset |
| `02_build_displacement_matrix.py` | Per-intervention displacement vectors; assemble matrix A |
| `03_svd_and_aging_direction.py` | SVD of A → aging direction v*; singular value / dimensionality test |
| `04_fit_principal_curve.py` | Orient the principal curve along v*; compute arc-length positions |
| `05_classify_interventions.py` | Classify each intervention by projecting its displacement onto v* |
| `06_validation.py` | Mortality / concordance / DamAge alignment validation (pending data) |
| `exploratory_smoking.py` | Single-dataset analysis of the smoking displacement |
| `exploratory_shared_component.py` | Decompose two interventions into shared vs specific components |

Shared infrastructure (cross-sectional cohorts, common CpG list, principal curve) is produced by
the parent repo's scripts and lives in the top-level `data/` directory.

---

## 5. Data situation

Acquiring usable intervention methylation data is the central practical bottleneck. The large
majority of GEO methylation datasets store only IDAT files (requiring R/minfi processing) or are
under controlled access. Many candidate accessions turned out to be RNA-seq, the wrong study, or
empty of beta values.

### Currently usable (public beta matrices)

All blood, on 450K or EPIC — directly comparable to the cross-sectional anchors.
The three EPIC datasets were added via a systematic FTP-verified GEO sweep
(see `CANDIDATE_DATASETS.md`); processed by `scripts/process_new_interventions.py`.

| Dataset | Intervention | Direction | Design | n |
|---------|-------------|-----------|--------|---|
| GSE50660 | Smoking (blood) | accelerating | cross-sectional | 22 current vs 179 never |
| GSE53045 | Smoking (PBMC) | accelerating | cross-sectional | 50 smoker vs 61 control |
| GSE140038 | Chemo + radiotherapy (blood) | accelerating | post vs pre, unpaired | 72 pre / 72 post |
| GSE89218 | Chronic stress / PTSD (blood) | accelerating | cross-sectional | 81 PTSD+ vs 82 PTSD− |
| GSE272137 | Bariatric weight loss (blood) | geroprotective | longitudinal (w0 vs w52) | 26 paired |
| GSE240184 | Behavioural weight loss, DRIFT2 (blood) | geroprotective | longitudinal (BL vs 3m) | 64 paired |
| GSE328810 | 8-week combined exercise (blood) | geroprotective | longitudinal (Before/After) | 13 paired |
| GSE193730 | Exercise, children 20 wk (blood) | geroprotective | longitudinal (Baseline/T1) | 10 paired (exercise arm) |

This is **8 interventions, balanced 4 accelerating / 4 geroprotective** — past the
4–6 needed for a first genuine v\*. The two smoking cohorts and two exercise
cohorts also serve as within-direction reproducibility controls.

### Requested from authors / pending access

| Source | Intervention | Direction | Status |
|--------|-------------|-----------|--------|
| CALERIE (Aging Research Biobank; Belsky) | Caloric restriction | geroprotective | access application + direct request |
| Fiorito (DAMA study) | Diet + physical activity RCT | geroprotective | requested |
| Janelsins (Yao et al. 2019) | Chemotherapy (450k, n=93) | accelerating | requested |
| Sehl (GSE133588) | Chemotherapy (EPIC, n=48) | accelerating | requested |
| Rönn | Exercise (adipose) | geroprotective | requested |
| Lindholm | Exercise (muscle) | geroprotective | declined — full beta matrix no longer available; DMP list usable for validation |

CALERIE is the highest priority: it is the central motivating dataset for the framework.

### Checked and rejected (do not re-pursue as methylation interventions)

- **GSE133588** — *not methylation.* Its supplementary `log2_norm.txt.gz` is
  Agilent **gene-expression** data (control probes `GE_BrightCorner`/`DarkCorner`,
  `CUST_*`), ~10 samples. Sehl's chemo methylation data is elsewhere / controlled.
- **GSE77716** — has a real ~2500-sample 450K beta matrix on FTP, but **no smoking
  phenotype** in the public metadata (only sex, cell fractions, ethnicity), so it
  cannot be used as a smoking intervention without labels from the authors.

See `CANDIDATE_DATASETS.md` for the full FTP-verified candidate list and backups.

---

## 6. Preliminary findings

> **Update (8 interventions + cell control + GTEx cross-tissue).** The pipeline now
> runs on eight blood interventions, balanced 4 accelerating (smoking ×2, chemo,
> PTSD) / 4 geroprotective (bariatric WL, behavioural WL, exercise ×2), plus a
> leukocyte-deconvolution control and a GTEx multi-tissue test. The two-dataset
> findings below (6a–6d) still hold. Current results are canonical in
> **[`results/FINDINGS.md`](results/FINDINGS.md)**. Headlines:
> - Still **no dominant shared aging axis** (ρ = λ₁/λ₂ = 1.30 raw), and it is
>   **not a cell-mix artifact** — orthogonality survives EpiDISH cell adjustment
>   (ρ → 1.12, all pairwise angles 75–99°).
> - **The four "accelerators" share no direction:** the non-treatment accelerator
>   PTSD is ~90° from smoking and chemo, so the disagreement isn't a chemo/
>   treatment artifact — accelerators occupy four different directions.
> - **The two independent smoking cohorts agree** (79–83° apart, both anti-aligned
>   to v\*) — reproducibility control showing the orthogonality is biological.
> - **GTEx cross-tissue:** aging is itself largely tissue-specific (ρ = 1.25;
>   exception: blood–lung age-directions 37° apart), and the intervention v\* is
>   **anti-aligned with real blood aging (cos = −0.72)** — interventions move cells
>   off, even against, the chronological-aging axis.
> - **Retraction:** the earlier "smoking vs chemo = 114.6°, opposite accelerators"
>   was largely a chemo cell-composition (lymphodepletion) artifact; adjusted they
>   are 87° (orthogonal, not anti-parallel).

### Two-dataset baseline (smoking + bariatric weight loss)

With only smoking and bariatric weight loss available, results are necessarily preliminary, but
they are coherent and point somewhere specific.

### 6a. Sanity check passed
The smoking displacement is dominated by canonical smoking CpGs — AHRR (cg05575921) is the single
largest-displaced CpG in the dataset (Δβ = −0.24), with F2RL3, PRSS23, ALPPL2, C1orf114 all in the
top percentile. The pipeline is measuring real biology.

### 6b. The two interventions are nearly orthogonal
The smoking and weight-loss displacement vectors are **81.8° apart** (cosine = +0.14) over 237,974
shared CpGs. They do not push in opposite directions along a shared axis. Consequently the
normalised singular value ratio is low (ρ = 1.65), and a robust one-dimensional aging direction
cannot be identified from these two interventions alone. This is the framework behaving correctly:
it declines to manufacture a clean aging axis from inputs that do not share one. Identifying v*
genuinely requires more interventions, so that a shared aging component can be separated from
large intervention-specific effects.

### 6c. Smoking appears to be a concrete case of clock-gaming
Against the chronological-age trajectory (156 shared CpGs):

- Smoking's displacement is **2.5× more concentrated on known aging-clock CpGs** than on other CpGs.
  This is why smokers register as epigenetically "older" on clocks.
- But the *direction* of smoking's displacement is **not aligned with the aging trajectory**
  (cosine = −0.34, ≈70°) — it is orthogonal-to-slightly-antiparallel.

In other words: smoking changes methylation heavily at exactly the sites clocks read, but in a
direction that is not the direction methylation moves with age. A clock computes a weighted sum
over its CpGs and reports "age acceleration"; it cannot see that the change is off-axis. The
geometric framework can. This is the clearest articulation so far of what clock-gaming means in
practice, and it is a hypothesis the literature has gestured at (some clock CpGs associate more
with smoking than age) but never tested geometrically.

### 6d. Bariatric weight loss is near-orthogonal to the age trajectory
Weight loss displacement aligns with the age trajectory at only cosine = +0.04 (≈87°) and shows
weak aging-CpG enrichment (1.2×). Over 52 weeks it does not appear to move cells substantially
along the chronological-age trajectory, consistent with the literature's mixed findings on whether
weight loss reverses epigenetic age versus changing metabolic CpGs.

### 6e. Six interventions + cell control: orthogonality is real, not cell mix

Over the ~41,600 CpGs shared by all six interventions, the displacement vectors are **nearly
mutually orthogonal** and the signed SVD spectrum is flat (ρ = λ₁/λ₂ = **1.29** raw). The decisive
test — is this just intervention-induced leukocyte-composition shifts? — is answered by EpiDISH
deconvolution and within-dataset residualisation on the estimated fractions:

- **Orthogonality survives cell adjustment**: ρ drops to **1.16** and every pairwise angle lands in
  **77–99°**. So "no shared aging axis" is not a cell-mix artifact.
- **But cell composition drove several raw vectors** — chemo most of all (cos(raw,adj) = 0.36; its
  raw displacement was largely lymphodepletion), bariatric partly (0.57). The earlier dramatic
  "smoking vs chemo 114.6°" was mostly that artifact; adjusted it is **88.7°**. *(Retracted.)*
- **Smoking reproduces**: the two independent smoking cohorts (whole blood vs PBMC) are the
  most-aligned pair raw (77.2°) and adjusted (82.2°) and project the same way onto v\* — evidence
  the geometry reflects biology, not single-cohort noise.
- **Forcing a provisional v\* still misclassifies both smoking cohorts (anti-aligned) and exercise
  (perpendicular).** v\* is defined by the chemo + weight-loss axis; the mismatches *are* the
  finding — these interventions don't share one axis.

Full tables and caveats in **[`results/FINDINGS.md`](results/FINDINGS.md)**.

---

## 7. Important caveats

- **Five interventions still do not define a robust axis.** ρ = 1.30 means v* is provisional, not
  validated; the orthogonality finding is robust but the v*-based classification (which flips
  smoking and exercise) should be read as "these interventions don't share an axis," not as a
  trustworthy per-intervention label. More interventions — and cell-composition adjustment — are
  needed before v* is meaningful.
- **Cell composition is an unaddressed confound.** All five are whole blood; intervention-induced
  shifts in leukocyte fractions can dominate the displacement and masquerade as (or bury) a shared
  aging axis. GSE328810 is cell-corrected; the others are not. Harmonising this is a priority.
- **The reference for the age-curve statements is the chronological-age curve, not v*.** Statements
  like "smoking is orthogonal to aging" against the age curve are a weaker claim than orthogonality
  to a validated aging direction; the v* projections in 6e are against the provisional v*.
- **The negative sign on smoking (−0.34)** could be a real phenomenon (smoking hypomethylation
  opposing age-related hypermethylation at shared sites) or an artefact of the small shared-CpG
  set (156). These cannot yet be distinguished.
- **Cross-dataset normalisation is unsolved.** Absolute arc-length projection across separately
  normalised datasets is dominated by batch effects; this is why classification uses within-dataset
  displacement projection rather than absolute position differences.
- **No mortality validation yet.** The ultimate anchor (does position on v* predict mortality in
  held-out data?) requires data not yet in hand.

---

## 8. What's needed next

1. **More interventions, especially geroprotective ones.** A robust v* needs at least 4–6
   interventions spanning both directions, so the shared aging component separates from
   intervention-specific noise. CALERIE is the priority.
2. **A within-study cross-tissue dataset** to test the context-invariance criterion directly
   (the GTEx multi-tissue methylation resource, GSE213478, is a candidate reference).
3. **Mortality-linked cohort** to test whether position on v* predicts all-cause mortality
   independent of any clock.
4. **DamAge/AdaptAge alignment test** — do the causally-enriched CpGs of Ying et al. 2024 align
   with v* more than standard clock CpGs? (Prediction: yes.)

---

## 9. The bigger picture

The framework gives a way to ask a question the field currently cannot: *does an intervention move
cells along the aging trajectory, or merely perturb methylation in other directions?* The
preliminary smoking result suggests the answer for smoking is "the latter" — its epigenetic
signature is concentrated on clock CpGs but off the aging axis. If this pattern holds for other
canonical accelerators once v* is robustly estimated, it would imply that much of the literature on
"X accelerates epigenetic aging" is measuring clock-specific artefacts rather than genuine
trajectory displacement — with direct consequences for how longevity interventions are evaluated.
Conversely, if a geroprotector like caloric restriction turns out to move cells genuinely along the
aging axis, that would be the first rigorous evidence distinguishing real geroprotection from
clock perturbation.

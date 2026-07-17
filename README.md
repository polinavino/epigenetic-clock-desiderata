# Epigenetic Clock Desiderata

**Paper:** *Towards a Formal Definition of Biological Age: Empirical Characterization of Clock Disagreement and Proposed Desiderata*  
**Author:** Polina Vinogradova  
**Preprint:** (forthcoming on bioRxiv)

---

## Plain-language summary

Your doctor could measure your "biological age" using a blood test — not your birthday, but how old your cells look. Several such tests exist, each producing a number in years. The problem: run all of them on the same person and they often disagree, sometimes by decades. A 55-year-old might score 48 on one test and 67 on another. Which one is right?

This paper argues that the question "which one is right?" is malformed — the tests are not measuring the same thing, and comparing them directly is like comparing your car's odometer reading to its speedometer reading. Both tell you something real about your car, but they answer different questions.

We make this precise. Using data from 1,385 blood samples across two large public studies, we:

1. **Show that the disagreement is structured, not random.** Tests designed to measure the same type of quantity agree with each other much more than tests designed to measure different types. When two tests of the same type disagree about who is biologically older, we can diagnose whether the disagreement is due to them tracking different trajectories through the data, or weighting biological deviations differently.

2. **Find a mathematical relationship between the test types.** Position-based tests (how much has your biology aged?) are linearly related to each other — they're essentially the same measurement on different scales. But rate-based tests (how fast are you aging right now?) are related to position tests via a *logarithm* — the current rate of aging increases sublinearly with accumulated biological age. This is new.

3. **Show that all tests are confounded by blood cell composition.** When your blood has more of one cell type than another (due to infection, stress, or normal variation), every test shifts — by up to r=0.37 correlation with cell type proportions. This is a technical artefact, not biology.

4. **Propose five formal criteria** (desiderata) that any biological age test should satisfy, and show which tests satisfy which criteria. No existing test satisfies all five.

5. **Show that "GrimAge" has a fundamentally different mathematical structure** from all other major tests — it is a two-stage pipeline rather than a single linear function, making it structurally incomparable to the others. This is a D0 (type declaration) violation.

---

## Scientific summary

DNA methylation — chemical marks on the genome that accumulate and drift over a lifetime — can be used to predict chronological age from blood with surprising accuracy. Over the past decade, multiple "epigenetic clocks" have been developed using different training objectives:

| Clock | Generation | Training objective | Our type |
|-------|-----------|-------------------|----------|
| Horvath (2013) | 1st | Chronological age | $\mathcal{T}_\tau$ (position) |
| Hannum (2013) | 1st | Chronological age | $\mathcal{T}_\tau$ (position) |
| PhenoAge (2018) | 2nd | Clinical phenotype composite | $\mathcal{T}_\delta$ (deviation) |
| GrimAge (2019) | 2nd | Time-to-death | $\mathcal{T}_\delta$† |
| DunedinPACE (2022) | 3rd | Rate of physiological decline | $\mathcal{T}_{\dot\tau}$ (rate) |

†GrimAge is structurally a two-stage composite model, not a linear clock.

### The canonical aging trajectory

We model the methylation profiles of a healthy population as tracing a path through methylation space. Formally, we define a **canonical aging trajectory** $\gamma: [0, T_{\max}] \to [0,1]^{|C|}$ as the principal curve through the data under an age-informativeness weighted norm:

$$\|m - m'\|^2_* = \sum_{i \in C} w_i (m_i - m'_i)^2, \quad w_i = R^2_i \cdot \sigma^2_i$$

where $R^2_i$ is the fraction of variance in CpG site $i$ explained by chronological age, and $\sigma^2_i$ is its population variance. This weights each genomic position by how much of its variation is age-related rather than noise.

Each sample's methylation profile $m$ is then decomposed into:
- $\tau(m)$ — biological age coordinate (arc-length projection onto $\gamma$)
- $r(m) = m - \pi(m)$ — residual (off-manifold component)

Position clocks estimate $\tau$, deviation clocks estimate $\tau + f(r)$, and rate clocks estimate $\dot{\tau}$.

### The five desiderata

**D0 — Type declaration:** Every clock must declare whether it measures position ($\mathcal{T}_\tau$), deviation ($\mathcal{T}_\delta$), or rate ($\mathcal{T}_{\dot\tau}$). Comparisons between clocks of different types are not well-formed.

**D1 — Monotonicity:** The expected clock output should increase with chronological age across the population.

**D2 — Directional stability:** A clock should be sensitive to motion *along* the canonical aging trajectory, and insensitive to motion *off* the trajectory (measurement noise, cell type shifts).

**D3 — Rank consistency:** Two clocks of the same declared type should agree on the relative biological age ordering of individuals.

**D4 — Intervention type consistency:** Clocks should respond to interventions in a manner consistent with their declared type — position clocks to accumulated changes, rate clocks to changes in current aging rate.

### Key empirical results

**Result 1 — Rank consistency matrix**

We compute Kendall's $\kappa$ (rank correlation rescaled to [0,1]) between all clock pairs on 1,385 blood samples:

| | Horvath | Hannum | PhenoAge | GrimAge | DunedinPACE |
|--|--|--|--|--|--|
| **Horvath** | 1.000 | 0.886 | 0.858 | 0.682 | 0.643 |
| **Hannum** | 0.886 | 1.000 | 0.867 | 0.691 | 0.662 |
| **PhenoAge** | 0.858 | 0.867 | 1.000 | 0.715 | 0.674 |
| **GrimAge** | 0.682 | 0.691 | 0.715 | 1.000 | 0.679 |
| **DunedinPACE** | 0.643 | 0.662 | 0.674 | 0.679 | 1.000 |

First-generation position clocks agree with each other ($\kappa$=0.886). DunedinPACE agrees least with everything ($\kappa$=0.643–0.679). The block structure confirms the type classification: same-type clocks agree more. Clock disagreement is higher in older cohorts (GSE40279, mean age 64: $\kappa_{HH}$=0.800) than younger ones (GSE87571, mean age 47: $\kappa_{HH}$=0.929).

**Result 2 — Functional relationships between clocks**

Position clock pairs are linearly related (R²=0.991–0.995, slopes ~1):

$$k_{\text{Hannum}} \approx 0.952 \cdot k_{\text{Horvath}} + b$$

The rate clock (DunedinPACE) is related to position clocks via a signed logarithm (R²=0.933–0.951 vs R²=0.748–0.793 for linear):

$$k_{\text{DunedinPACE}} \approx 0.0135 \cdot \log(|k_{\text{position}} + c| + 1) \cdot \text{sign}(k_{\text{position}} + c) + b$$

The coefficient $a \approx 0.0135$ is consistent across all three position clocks; only the shift $c$ varies (reflecting their different scales). This saturation means the current aging rate increases sublinearly with accumulated biological age — a Weber-Fechner-like relationship.

**Result 3 — Rank reversal decomposition**

For each clock pair, we sample random pairs of individuals and classify rank reversals as:
- *Tau-dominated*: individuals are close on the canonical trajectory (small $|\Delta\tau|$), so residual functions $f(r)$ drive the disagreement
- *Residual-dominated*: individuals are clearly separated on the trajectory but clocks weight their off-manifold profiles differently

GrimAge pairs show ~45–47% residual-dominated reversals vs ~35% for same-type position clock pairs, confirming that deviation clocks disagree more due to different residual weighting.

**Result 4 — Directional stability (D2)**

Estimating each clock's sensitivity to on-manifold vs off-manifold perturbations via Ridge regression on the top 200 age-informative CpG sites:

| Clock | Type | $\rho$ (off/on ratio) |
|-------|------|----------------------|
| Horvath | $\mathcal{T}_\tau$ | 3.47 |
| Hannum | $\mathcal{T}_\tau$ | 3.51 |
| PhenoAge | $\mathcal{T}_\delta$ | 3.51 |
| GrimAge | $\mathcal{T}_\delta$ | 6.53 |
| DunedinPACE | $\mathcal{T}_{\dot\tau}$ | 10.04 |

Higher $\rho$ = more off-manifold sensitivity. DunedinPACE is most sensitive to off-manifold variation (expected for a rate clock under violations of the autonomy assumption T3). GrimAge is intermediate. Position clocks are most stable.

**Result 5 — Cell type confounding (D3)**

Blood cell type proportions (estimated via reference-based deconvolution) correlate significantly with age-acceleration residuals for all clocks:

| Clock | Max |r| with cell types | D3 |
|-------|--------------------------|-----|
| Horvath | 0.293 | FAIL |
| Hannum | 0.372 | FAIL |
| PhenoAge | 0.350 | FAIL |
| GrimAge | 0.278 | FAIL |

All clocks fail D3. Hannum is most confounded; GrimAge least.

**Result 6 — GrimAge structural finding**

GrimAge is a two-stage composite model: it first predicts protein biomarker levels (GDF15, B2M, cystatin C, etc.) from methylation using separate linear sub-models, then combines these predictions. This is structurally incomparable to the other clocks, which are single linear functions. Coefficient vector analysis confirms near-orthogonality. The off-diagonal cosine similarities range from −0.171 (PhenoAge–DunedinPACE) to 0.115 (Horvath–PhenoAge), with Horvath–DunedinPACE = 0.000 (verified from `data/cosine_similarity.parquet`). All |cosine| ≤ 0.17, so the clocks measure genuinely independent biological signals. (An earlier version of this section reported the range as "0.03–0.11", which omitted the negative PhenoAge–DunedinPACE value; the near-orthogonal conclusion is unchanged.)

**Result 7 — Coherence**

Regressing each clock's age acceleration on $(\tau, \|r\|_*)$ yields R²=0.09–0.28 — only 9–28% of clock variance is explained by the 1D canonical trajectory. DunedinPACE is most coherent (R²=0.279, driven by $\tau$), confirming its nature as a rate clock. The low coherence for position clocks suggests the 1D manifold is an approximation, motivating a higher-dimensional or multi-modal extension.

---

## Reproducing the results

### Requirements

```bash
pip install -r requirements.txt
# Also required:
pip install biolearn tables pyarrow
```

All dependencies: Python 3.10+, numpy, pandas, scipy, scikit-learn, matplotlib, seaborn, statsmodels, biolearn, tables, pyarrow.

### Data

Downloaded automatically by script 1 from NCBI GEO. Requires ~4GB disk space and ~8GB RAM for preprocessing. All data files are excluded from git and regenerated by the scripts.

**Datasets used:**
- [GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) — Hannum et al. 2013, n=656, ages 19–101, whole blood 450k array
- [GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571) — Johansson et al. 2013, n=729, ages 14–94, whole blood 450k array

### Scripts

| Script | Description | Runtime |
|--------|-------------|---------|
| `scripts/01_download_and_preprocess.py` | Download GEO data, QC, save HDF5 | ~70 min |
| `scripts/02_compute_weights_and_curve.py` | Compute $w_i$ weights, fit principal curve, compute $\tau$ and $r$ | ~30 min |
| `scripts/03_compute_clocks_and_ranks.py` | All clock computation and analyses (Steps 1–14) | ~20 min |
| `scripts/principal_curve.py` | Hastie-Stuetzle principal curve implementation (utility, not run directly) | — |
| `scripts/config.py` | All paths and parameters (no hardcoding elsewhere) | — |

### Replication

```bash
python scripts/01_download_and_preprocess.py
python scripts/02_compute_weights_and_curve.py
python scripts/03_compute_clocks_and_ranks.py
```

All figures are saved to `paper/figures/`. All intermediate data to `data/`. Nothing is hardcoded — all paths and parameters are in `scripts/config.py`.

### Key parameters (in `config.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N_TOP_CPGS` | 200 | CpG sites used for principal curve |
| `PC_SMOOTHING` | 0.5 | Principal curve smoothing (× n_samples) |
| `RANDOM_SEED` | 42 | Random seed for reproducibility |

---

## Repository structure

```
.
├── README.md
├── requirements.txt
├── scripts/                             # main methylation-clock pipeline (heavy)
│   ├── config.py                        # all paths and parameters
│   ├── principal_curve.py               # Hastie-Stuetzle implementation
│   ├── 01_download_and_preprocess.py
│   ├── 02_compute_weights_and_curve.py
│   └── 03_compute_clocks_and_ranks.py
├── data/                                # populated by scripts, not tracked (large)
├── outputs/                             # tracked console outputs of the light analysis scripts
├── intervention-framework/             # sub-project: geometric intervention classification
│   ├── scripts/                         # its own pipeline + light analysis scripts
│   └── results/                         # its figures
├── eye_aging/                           # extension of the clock instance to retina (transcriptomic)
│   ├── analysis/eye_clocks.py           # self-teeing to eye_aging/analysis/outputs/
│   └── data/                            # rat retinal aging data (GSE314970) + senescence gene sets
└── paper/
    ├── figures/                         # all output figures
    └── sections/                        # LaTeX sections (forthcoming)
```

## Sub-projects and the `outputs/` convention

- **`intervention-framework/`** is a distinct sub-project that classifies interventions geometrically
  (geroprotective / age-accelerating / clock-gaming) against the aging trajectory. It has its own
  desiderata and should not be conflated with the clock-comparison analyses above.
- **`eye_aging/`** extends the clock instance to a new tissue and modality: a rat retinal bulk RNA-seq
  aging series (GSE314970), on which several competing transcriptomic aging measures (SenMayo and Fridman
  senescence signatures, an inflammaging proxy, and a fitted elastic-net clock) are compared with the same
  framework. The fitted clock tracks age at Spearman 0.98 while the curated senescence signatures track it
  weakly (down to 0.13), they split into a senescence/inflammaging family and a distinct clock direction,
  and the consensus tracks age. It is a single-cohort, rat, single-sex extension, so it is qualitative.
  See `eye_aging/README.md`.
- **`outputs/` convention.** The light, self-contained analysis scripts tee their console output to a
  tracked text file, so the reported numbers always have a stored source of truth. The repo-level
  `outputs/` holds the retrofitted light `intervention-framework` analysis scripts
  (`03_svd_and_aging_direction`, `04_fit_principal_curve`, `05_classify_interventions`,
  `cell_composition_analysis`, `exploratory_shared_component`, `exploratory_smoking`); `eye_aging` keeps
  its own `analysis/outputs/`. The heavy main pipeline (`scripts/01–03`, which downloads GEO data and runs
  biolearn) is not retrofitted — its results are the derived parquet tables in `data/`.

---

## Figures

| Figure | Description |
|--------|-------------|
| `03_rank_consistency_matrix.png` | Kappa matrix across all clock pairs and datasets |
| `03_functional_relationships.png` | Linear within position clocks; logarithmic vs DunedinPACE |
| `03_tau_r_decomposition.png` | Rank reversal decomposition: tau-dominated vs residual-dominated |
| `03_directional_stability.png` | D2: off-manifold sensitivity ratio by clock type |
| `03_d3_cell_type.png` | D3: cell type confounding for all clocks |
| `03_cosine_similarity.png` | Geometric relationships between clock coefficient vectors |
| `03_coherence_test.png` | How well does (τ, ‖r‖) explain each clock? |
| `03_d1_monotonicity.png` | D1: monotonicity test across age deciles |
| `03_pearson_vs_kappa.png` | Pearson r vs Kendall κ: why r overstates agreement |
| `03_high_instability_individuals.png` | Individuals whose biological age rank varies most across clocks |
| `02_weights_and_curve.png` | Principal curve validation: τ vs chronological age |

---

## Relation to prior work

The closest precursor is Klemera and Doubal (2006), who made the same complaint — that biological age lacks an exact definition, making comparisons between methods meaningless — but predated methylation clocks entirely and worked only with linear clinical biomarkers.

This paper applies an analogous programme to the kinase inhibitor selectivity problem in [Vinogradova (2025)](https://github.com/polinavino/kinase-selectivity-definitions), where the same issue — multiple competing scalar summaries of a complex profile, treated as interchangeable — was addressed by formal desiderata and empirical instability analysis.

---

## Citation

```bibtex
@article{vinogradova2026clocks,
  title   = {Towards a Formal Definition of Biological Age: 
             Empirical Characterization of Clock Disagreement 
             and Proposed Desiderata},
  author  = {Vinogradova, Polina},
  year    = {2026},
  note    = {Preprint, forthcoming on bioRxiv}
}
```

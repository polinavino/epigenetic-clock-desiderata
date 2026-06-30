# Intervention dataset hunt — FTP-verified candidates

Result of a systematic GEO sweep (NCBI eutils `db=gds`,
`"methylation profiling by array"[DataSet Type] AND Homo sapiens[Organism]`)
across geroprotective and accelerating intervention terms, FTP-verifying each
hit for **(a)** a public CpG beta matrix on the FTP `suppl/` directory and
**(b)** usable phenotype + intervention design in the series matrix.

Last run: 2026-06-29.

---

## Ingested into the pipeline (verified: public betas + phenotype + design)

All three are whole blood on Illumina EPIC, directly comparable to the existing
blood datasets (GSE50660 smoking, GSE272137 bariatric, and the GSE40279/GSE87571
cross-sectional anchors). Processed by `scripts/process_new_interventions.py`.

| Accession | Intervention | Direction | Design | n | Beta file (FTP `suppl/`) |
|-----------|-------------|-----------|--------|---|--------------------------|
| **GSE140038** | Chemo + radiotherapy, early breast cancer | accelerating (+) | post vs pre, **unpaired** (no subject IDs public; timepoints not interleaved) | 144 (72 tp0 / 72 tp1) | `GSE140038_NormalizedBetaNoob.csv.gz` |
| **GSE240184** | Behavioural weight loss (DRIFT2 trial) | geroprotective (−) | paired baseline (BL) vs 3-month (T3M) | 128 (64 pairs) | `GSE240184_betasGEO.txt.gz` |
| **GSE328810** | 8-week combined exercise, women w/ obesity | geroprotective (−) | paired Before/After | 26 (13 pairs) | `GSE328810_beta_normalized_cell_corrected.csv.gz` |

With these, the usable set goes from 2 → **5 interventions** spanning both
directions: **2 accelerating** (smoking GSE50660, chemo/radio GSE140038) and
**3 geroprotective** (bariatric GSE272137, behavioural weight loss GSE240184,
exercise GSE328810) — enough for a first genuine v* estimate (target was 4–6).

> Note on GSE140038: the public metadata has no patient/subject IDs and the
> tp0/tp1 ordering is not interleaved, so pre/post cannot be paired. It is
> processed as an unpaired group contrast `mean(post) − mean(pre)`. This is
> still within-study, so batch effects cancel as with the paired designs.

---

## Verified-available backups (public betas confirmed; not yet ingested)

| Accession | Intervention | Direction | Tissue | n | Notes |
|-----------|-------------|-----------|--------|---|-------|
| GSE213363 | Resistance+aerobic exercise, PCOS women | gero (−) | blood | 112 | EPIC; `matrix_processed.csv.gz` (betas + detection p interleaved). PCOS confound. |
| GSE193730 | 20-wk exercise, children w/ overweight (ActiveBrains) | gero (−) | blood | 46 | `BETA_VALUES_ACTIVEBRAIN.csv.gz`; European decimal commas; paediatric. |
| GSE53045 | Smoking (PBMC) | accel (+) | PBMC | 111 | 450K; `matrix_processed_GEO.txt.gz`. Second smoking cohort (redundant w/ GSE50660). |
| GSE268211 | HIIT skeletal-muscle "memory" | gero (−) | muscle | 20 | `..._matrix_processed_for_GEO.txt.gz`. Muscle (tissue arm). |
| GSE213029 | Aerobic exercise rejuvenates muscle methylome | gero (−) | muscle | 32 | only `signal_intensities` → must back-calc beta = M/(M+U). |
| GSE171140 | Exercise, skeletal muscle (4/8/12 wk) | gero (−) | muscle | 195 | only `signal_intensities` → back-calc needed. |

Muscle datasets are kept separate because cross-tissue projection onto a
blood-derived axis is confounded; useful later for the context-invariance test.

---

## Checked and rejected

| Accession | Why rejected |
|-----------|--------------|
| **GSE133588** | **Not methylation.** Supplementary `log2_norm.txt.gz` is Agilent gene-expression (probe IDs `GE_BrightCorner`, `DarkCorner`, `CUST_*`), ~10 samples. The chemo *methylation* data sought from Sehl is elsewhere / controlled. Remove from config. |
| **GSE77716** | Beta matrix exists (`GSE77716_Matrix_processed.tsv.gz`, ~2500 samples, real `cg` betas) **but the public metadata has no smoking phenotype** — only sex, blood cell fractions, ethnicity (Latino cohort). Unusable as a smoking intervention without phenotype from authors. |
| GSE109914 | Arsenic exposure, blood — `suppl/` has **RAW IDATs only**, no public beta matrix. |

---

## How to add a verified candidate

1. Add an entry to `INTERVENTION_DATASETS` in `scripts/config_intervention.py`
   (`label`, `sign` plus/minus, `tissue`, `design` longitudinal/cross_sectional,
   `source: processed_direct`).
2. Download the beta file into `data/raw/`.
3. Add a processor branch in `scripts/process_new_interventions.py` that emits
   `{ACC}_beta.h5` (samples × CpGs) + `{ACC}_metadata.csv` with `sample_id`
   plus either `subject_id`+`timepoint` (longitudinal) or `group` (cross_sectional).
4. Re-run `02` → `05`.

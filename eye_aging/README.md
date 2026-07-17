# Eye aging — an extension of the epigenetic-clock instance to a new tissue

This folder extends the clock instance (this repository) from blood DNA-methylation clocks to the retina,
and from methylation to transcriptomic aging measures. It is folded into the meta-paper's clock
demonstration (Section 4.2 of the measurement-protocol manuscript). It is not a standalone domain. Every
number below is copied from the stored output `analysis/outputs/eye_clocks.txt` (the single source of
truth). Rerun `python analysis/eye_clocks.py` after any edit.

## Why this exists
Recent human eye-cell aging work is single-cell with too few donors (6–18) to compare aging clocks, so it
cannot power a measure-agreement analysis. The tractable substrate is a well-powered bulk time-course,
which exists for rat retina (the mouse companion atlases contain no eye tissue). So this is a rat,
transcriptomic, single-sex demonstration, treated honestly as a qualitative extension.

## Data (`data/`)
- **GSE314970** (Shavlakadze et al. 2026, Calico/Regeneron multi-tissue aging atlas), rat neural-retina
  bulk RNA-seq. 85 samples across 7 ages (6, 9, 12, 18, 21, 24, 27 months, ~12 each), male, one tissue.
  30,454 genes, raw counts. Age is continuous, a better anchor than the discrete anchors elsewhere.
- Competing aging measures (all oriented so higher = older): **SenMayo** (SAUL_SEN_MAYO, 113/124 genes
  matched), **Fridman senescence-up** (74/77), **Fridman senescence-down** negated (12/13), a **Hallmark
  inflammatory** inflammaging proxy (199/200), and a **fitted elastic-net transcriptomic clock**
  (out-of-fold predictions). Anchor = chronological age.

## Results (from `analysis/outputs/eye_clocks.txt`)
- **The measures disagree strongly on the aging axis.** Against true age, the fitted clock scores
  Spearman 0.979, the inflammatory proxy 0.616, SenMayo 0.448, the down-signature 0.302, and Fridman-up
  only 0.128. Curated senescence signatures are weak individual age-trackers while the fitted clock is
  near-perfect.
- **Families.** SenMayo, Fridman-up, and the inflammatory proxy cluster (r 0.72–0.90), the down-signature
  is the outlier, and the fitted clock sits partly apart (r 0.15–0.66).
- **Consensus tracks age.** The mean-rank consensus correlates with age at 0.755 and rises monotonically
  across tertiles (mean age 9.7 → 18.3 → 22.7 months).
- **Near-tie law holds.** Discordance falls with consensus separation (0.99 → 0.54), though the overall
  rate is high (0.80) and the separation signal is weaker (R² 0.185) than in the well-behaved domains.
- **Canonical aggregate.** Consensus poset 80.2% incomparable; Bubley–Dyer chains converge to Spearman
  0.94 (approximate at this incomparability). The inflammatory proxy (0.787) and the fitted clock (0.741)
  best proxy the canonical average-rank, the down-signature worst (0.280). The canonical tracks age 0.699.

## Honest limits
Rat, not human. One sex, one tissue, one dataset, so no cross-cohort reproducibility was run (the study's
RPE/choroid subset, 75 samples, is the obvious second cohort and a clean next step). The fitted clock is
trained on these samples (out-of-fold, so not circular, but still in-distribution). A qualitative
extension on par with the thin inflammaging-clock topic, not a strong standalone domain.

## Files
`analysis/eye_clocks.py` self-tees to `analysis/outputs/eye_clocks.txt`. `data/` holds the counts matrix
(the 5.9 MB parquet is tracked so the analysis reruns; the redundant 13 MB CSV is git-ignored), the
metadata, the gene-symbol map, the MSigDB senescence gene sets, the Hallmark inflammatory set, and the
derived `eye_scores.csv`.

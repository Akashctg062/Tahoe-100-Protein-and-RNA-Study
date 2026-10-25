# Consolidated results — 16 Oct 2026 (all analyses complete)

Study: *What does a 100-million-cell drug-perturbation atlas predict? A pre-specified evaluation of Tahoe-100M against clinical cardiotoxicity and cancer-cell drug sensitivity, with RNA and protein baselines.*

## Data
Tahoe-100M pseudobulk DE (379 drugs × 50 measured lines × 3 doses; 15 lung lines) · FDA DICTrank · FDA FAERS (openFDA) · PRISM 20Q2 secondary (MTS010) · GDSC2 · CCLE RNA-seq (Ghandi 2019) · CCLE MS proteomics (Nusinow 2020) · CCLE RPPA. All public, de-identified.

## Discovery analyses (exploratory)
| ID | Question | Result |
|---|---|---|
| D2 | Are Tahoe signatures reproducible and mechanism-specific? (positive control) | **Yes.** Mid vs high dose, same drug ρ = 0.68 vs 0.23 between drugs; anticancer MoA 1-NN balanced accuracy 0.32 vs chance 0.09 (perm. 95th pct 0.15) |
| D1 | Do tumour signatures predict DICTrank cardiotoxicity (100 compounds)? | **No.** Chemistry AUROC 0.67; Tahoe 0.46–0.54 across 5 representations; chemistry + Tahoe ≤ chemistry; same for lung-only, Most+Less, anticancer subset, DICTrank subtypes |
| D3 | Do 24-h pooled cell counts measure killing? | **No.** ρ ≈ 0 with PRISM/GDSC (raw and pool-share); poor dose consistency; salt-form replicates disagree |
| D4 | Does response magnitude predict which lines are killed? | **No.** Within-drug median ρ ≈ −0.02 (PRISM 51 drugs, p ≈ 0.3; GDSC 43 drugs, one nominal p = 0.01 not surviving correction) |

## Confirmatory analyses (pre-specified, frozen 15 Oct 2026)
| ID | Question | Frozen result | Amendment |
|---|---|---|---|
| C1 | Does Tahoe add to chemistry for FAERS cardiac signal? (128 drugs) | AUROC chem 0.459, Tahoe 0.444, chem+Tahoe 0.456; Δ −0.003 (−0.012, +0.006) → **does not add** | A2: excluding CV-indication drugs, unchanged (0.487 / 0.432 / 0.484) |
| C2 | Proliferation-program drop predicts sensitivity? | **Not supported** (PRISM median ρ +0.015, p 0.34; GDSC +0.055, p 0.035; threshold 0.025) | Post hoc: per-drug ρ agrees across screens (ρ = 0.44, 22 drugs); targeted cytostatic agents strongest |
| C3 | Which data predict which lines a drug kills? (≈27 lines/drug; PRISM 43, GDSC 37 drugs) | No set passes frozen "median ρ > 0.10" rule; frozen pairwise "Response worse" = LOO artefact | A1 (null-calibrated): RPPA +0.17 GDSC (FDR 0.039), +0.10 PRISM; MS protein +0.07/+0.13; RNA ≈ 0; Response +0.10 GDSC / −0.05 PRISM; RPPA > RNA in GDSC (FDR 0.039) |
| Robustness | k = 3, 10; α = 0.1, 10 | — | GDSC RPPA positive in 5/5 settings (p ≤ 0.004); RPPA > RNA 4/5; RNA ≈ 0; Response not reproducible across screens |

## One-sentence conclusion
Tahoe-100M signatures robustly encode drug mechanism but carry no reproducible information about cardiotoxicity, early cell killing or which cancer lines a drug kills; baseline protein (especially RPPA phospho-/total protein) carries a modest, reproducible signal for cell-line sensitivity that baseline RNA does not.

## Caveats to state in the manuscript
Small effect sizes; ≈27 lines per drug; cell lines only; single 24-h time point; label noise (DICTrank is label-based; FAERS confounded by indication); PRISM vs GDSC platform differences; amendments are post hoc and labelled.

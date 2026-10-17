---
title: "Amendment 3 — Pre-specification of extension analyses (Tier 1 and Tier 2)"
study: "What does a 100-million-cell drug-perturbation atlas predict? (extended journal version, v2)"
author: "Asaduzzaman (GIBMTE, Taipei Medical University) — prepared with Claude"
date: "Frozen 18 October 2026, before any Tier 1 analysis was run and before any Tier 2 dataset was downloaded"
---

# 0. Purpose and status

This amendment fixes the hypotheses, data, models, metrics and decision rules of the extension analyses **before** they are run. It does not change the frozen plan (15 Oct 2026) or Amendments 1–2; the conference manuscript (v1, 16 Oct 2026) is left unchanged.

* **Tier 1** re-uses data already analysed (Tahoe-100M, DICTrank, FAERS, PRISM, GDSC, CCLE). Its results are **secondary / robustness** evidence, except T1-C, which tests a new hypothesis (antibody–target concordance) at a level of the data not yet examined.
* **Tier 2** uses **independent data not yet downloaded** (CTRP v2, ProCan-DepMapSanger proteomics, LINCS L1000, iPSC-cardiomyocyte signatures, ChEMBL hERG). Its tests are **confirmatory replications**.
* Every result is reported whatever its direction. Any later change is recorded as a dated note in Section 9.
* Seeds: 2026 unless stated. Multiple testing: Benjamini–Hochberg (BH) within each family named below; one-sided tests where a direction is pre-specified.

# 1. T1-A — "Best-effort" models (are the null results a model limitation?)

## T1-A1 Cardiotoxicity (D1 task "most vs no concern", n = 100; C1 FAERS, n = 128)
* **Tahoe representations:** (i) PCA-50 of the concatenated three-dose mean signature (fitted inside each training fold); (ii) 14 gene programs × 3 doses + 6 potency scalars (as in D1).
* **Chemistry:** RDKit 2-D descriptors + Morgan (r = 2, 1,024 bits), as in D1/C1.
* **Model families (hyper-parameters tuned in an inner stratified 5-fold CV, AUROC):**
  1. Elastic-net logistic regression: C ∈ {0.001, 0.01, 0.1, 1}, l1_ratio ∈ {0.1, 0.5, 0.9}
  2. Histogram gradient boosting: learning_rate ∈ {0.03, 0.1}, max_depth ∈ {2, 3}, max_iter = 200
  3. RBF support-vector machine: C ∈ {0.1, 1, 10}, gamma = "scale"
  4. Multilayer perceptron: hidden units ∈ {16, 64}, alpha ∈ {0.001, 0.1}, early stopping
* **Outer evaluation:** stratified 5-fold CV repeated 10 times; out-of-fold probabilities averaged over repeats.
* **Primary metric (per family):** ΔAUROC = AUROC(chemistry + Tahoe) − AUROC(chemistry), with 2,000 bootstrap resamples of drugs.
* **Decision rule:** Tahoe is "rescued" only if, for at least one family, the Holm-adjusted lower confidence limit of ΔAUROC (98.75% CI for 4 families) is > 0, for C1. D1 results are reported as secondary.
* Also reported: Tahoe-alone AUROC versus a 1,000-permutation null.

## T1-A2 Sensitivity (C3 common-line sets, PRISM 43 drugs, GDSC 37 drugs)
* **Feature sets:** Response, RNA, Protein_MS, RPPA, Response+RNA, Response+Protein_MS (as in C3); each block standardised and reduced to 10 principal components inside every training fold.
* **Model families:** (1) ridge with internal generalised CV over α ∈ 10^[−2, 3] (12 values); (2) elastic net (internal 5-fold CV, l1_ratio ∈ {0.2, 0.5, 0.8}); (3) random forest (200 trees, min_samples_leaf = 3, max_features = 0.33); (4) RBF kernel ridge (α = 1, γ = 1/number of features).
* **Evaluation:** within-drug 5-fold CV (fixed seed) instead of leave-one-out, for computational tractability; per-drug Spearman ρ between predicted and observed AUC, **null-calibrated** by the mean of 50 within-drug label permutations passed through the identical procedure (same permutations for every set and model).
* **Decision rules:** (a) Response is "rescued" if, for at least one model family, Response+RNA > RNA or Response+Protein_MS > Protein_MS (two-sided Wilcoxon on calibrated ρ, BH within the model family) at FDR < 0.05 **in both screens**; (b) the protein-over-RNA result is called "model-robust" if RPPA − RNA has a positive median calibrated difference in GDSC for ≥ 3 of 4 families.

# 2. T1-B — Leave-one-out (LOO) bias: simulation and analytic description (descriptive)
* **Analytic case:** for an intercept-only predictor, the LOO prediction is (S − y~i~)/(n − 1), so corr(ŷ, y) = −1 exactly for any data.
* **Simulation grid (null: y independent of X):** n ∈ {10, 15, 20, 27, 40, 60}; p = 5,000 Gaussian features reduced in-fold to k ∈ {1, 3, 5, 10} components; ridge α ∈ {0.1, 1, 10, 100}; 500 replicates per cell. Report mean LOO Spearman ρ (raw), the permutation-calibrated ρ, and 5-fold-CV ρ.
* **With signal:** y = βz + ε with true correlation ∈ {0.2, 0.4} on a latent component; report bias of raw and calibrated estimates.
* No hypothesis test; the output is a methods figure and a short recommendation.

# 3. T1-C — RPPA antibody–drug-target concordance (new hypothesis, existing data)
* **Hypothesis:** for drugs whose targets lie in a pathway measured by RPPA, the antibodies of that pathway are more strongly associated with the drug's within-drug sensitivity than other antibodies.
* **Units:** drug × screen pairs from the C3 common-line sets (PRISM, GDSC) whose Tahoe target annotation maps to ≥ 1 pathway group below.
* **Pre-specified pathway groups (genes; antibodies are mapped to genes by name, e.g. "MAPK_pT202_Y204" → MAPK1/MAPK3):**
  * ERBB: EGFR, ERBB2, ERBB3
  * MAPK: ARAF, BRAF, RAF1, MAP2K1, MAP2K2, MAPK1, MAPK3, RPS6KA1–3
  * PI3K/AKT/mTOR: PIK3CA, PIK3R1, AKT1–3, MTOR, RPS6KB1, RPS6, EIF4EBP1, PDPK1, PTEN, TSC1, TSC2, AKT1S1, RICTOR, RPTOR, GSK3A, GSK3B
  * Cell cycle/CDK: CDK1, CDK2, CDK4, CDK6, RB1, CCND1, CCNE1, CCNE2, CCNB1, CDKN1A, CDKN1B
  * RTK (MET/KIT/VEGFR/PDGFR/FLT3): MET, KIT, KDR, PDGFRA, PDGFRB, FLT3
  * SRC family: SRC, LCK, YES1, LYN
  * JAK/STAT: JAK1, JAK2, STAT3, STAT5A
  * Apoptosis/BCL2: BCL2, BCL2L1, BAX, BAK1, BAD, BCL2L11
  * DNA damage: ATM, CHEK1, CHEK2, PARP1, XRCC1, ERCC1, MRE11, RAD50, RAD51, BRCA2, TP53BP1
  * Nuclear receptor: ESR1, PGR, AR
* **Statistic:** per unit, Spearman ρ of each of the 214 antibodies with AUC across common lines; concordance = mean percentile rank (by |ρ|) of the target-matched antibodies. Summary = median over units.
* **Null:** 10,000 permutations re-assigning pathway target sets among units within screen.
* **Decision rule:** concordance supported if one-sided permutation p < 0.05 (pooled over screens); per-screen results reported. Secondary: sign of ρ for phospho-antibodies of the target pathway (expected negative = higher pathway activity in sensitive lines).

# 4. T1-D — Middle-dose robustness (from the original frozen plan)
* Re-run Amendment 1 (calibrated C3, k = 5, α = 1, LOO, 500 permutations) with the Tahoe **middle** dose instead of the top dose. Baseline sets are unchanged, so only Response-containing sets can differ. Reported as robustness.

# 5. T2-E — Independent screen replication: CTRP v2
* **Data:** CTRP v2 per-drug, per-line AUC (original CTD^2 release if reachable; otherwise the re-processed CTRPv2 table of the drevalpy dataset, Zenodo 10.5281/zenodo.17600096, cited as such).
* **Analysis:** Amendment-1 calibrated C3 unchanged (k = 5, α = 1, LOO, 500 permutations), drugs with ≥ 20 common lines.
* **Hypotheses (one-sided Wilcoxon across drugs):** E1 RPPA calibrated ρ > 0; E2 RPPA − RNA > 0. Secondary: all other sets and pairs (BH).
* **Decision rule:** the GDSC RPPA finding is "replicated in CTRP" if E1 and E2 both have p < 0.05.
* Overlap: drugs also present in PRISM/GDSC are kept (different assay); a sensitivity analysis restricted to drugs absent from GDSC is reported.

# 6. T2-F — Second mass-spectrometry proteome: ProCan-DepMapSanger (949 lines)
* **Data:** Gonçalves et al. 2022 protein matrix (Zenodo 10.5281/zenodo.6563157).
* **Analysis:** calibrated C3 with "Protein_ProCan" replacing Protein_MS, for PRISM, GDSC and CTRP; the common-line set requires ProCan, RNA and RPPA, and all sets are recomputed on that set for paired comparison.
* **Hypotheses:** F1 ProCan protein calibrated ρ > 0; F2 ProCan − RNA > 0 (one-sided Wilcoxon).
* **Decision rule:** "protein > RNA" is supported by ProCan if F2 has p < 0.05 in ≥ 2 of the 3 screens.

# 7. T2-G — Bulk L1000 versus pooled single-cell Tahoe, same drugs and lines (conditional)
* **Condition:** carried out only if LINCS L1000 Level-5 signatures can be obtained for ≥ 4 cell lines and ≥ 20 drugs shared with Tahoe at 24 h; otherwise reported as not feasible.
* **Signatures:** L1000 consensus (MODZ) signature at 24 h at the dose closest to the Tahoe top dose; landmark genes plus well-inferred genes; Tahoe top-dose log~2~FC on the same genes.
* **Tests (per shared line, across shared drugs; label = PRISM, GDSC or CTRP AUC):** G1 Spearman ρ between the 20-gene proliferation score (C2 list, genes available on both platforms) and AUC; G2 the same for response magnitude (root-mean-square signature).
* **Comparison:** paired one-sided Wilcoxon over (line, screen) units, hypothesis ρ~L1000~ > ρ~Tahoe~.
* **Decision rule:** "bulk L1000 carries viability information that Tahoe lacks" if G1 p < 0.05.

# 8. T2-H — Cardiac context
## T2-H1 iPSC-cardiomyocyte signatures (Hansen et al. 2024, supplementary differential-expression tables; 6 donors, 48 h)
* (a) **Cross-context identity:** for drugs present in both resources, correlation of the Tahoe signature (all-line mean, top dose) with the cardiomyocyte signature of the same drug versus different drugs (shared genes); one-sided permutation test (10,000 label shuffles), p < 0.05.
* (b) **Context matters for toxicity:** if ≥ 30 shared drugs have DICTrank or FAERS labels, compare leave-one-out logistic AUROC of cardiomyocyte versus Tahoe signatures (bootstrap ΔAUROC, 2,000); otherwise descriptive only.
## T2-H2 hERG inhibition (ChEMBL target CHEMBL240)
* **Labels:** human hERG IC~50~/K~i~ records (standard units nM, exact relations "="), median pIC~50~ per compound; matched to Tahoe drugs by InChIKey connectivity block computed from SMILES.
* **Tasks:** regression (pIC~50~; Spearman ρ of out-of-fold predictions) and classification (pIC~50~ ≥ 5, i.e. IC~50~ ≤ 10 µM; AUROC), chemistry versus chemistry + Tahoe PCA-10, models and CV as in C1.
* **Decision rule:** Tahoe adds information if the lower 95% bootstrap confidence limit of Δρ (regression, primary) is > 0.

# 9. Change log
* (none yet)

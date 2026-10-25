---
title: "Pre-specified Analysis Plan and Manuscript Outline"
subtitle: "What does a 100-million-cell drug-perturbation atlas predict? A pre-specified evaluation of Tahoe-100M against clinical cardiotoxicity and cancer-cell drug sensitivity, with RNA and protein baselines"
author: "Asaduzzaman (GIBMTE, Taipei Medical University) — draft prepared with Claude"
date: "Frozen: 15 October 2026 (before running confirmatory analyses C1–C3)"
---

# 1. Why this plan is frozen now

Exploratory (discovery) analyses were run on 14–15 Oct 2026 (Section 3). To protect the paper against "testing until something is significant", the **confirmatory analyses C1–C3, their metrics, thresholds and decision rules are fixed here before they are run**. Any analysis added later will be labelled *post hoc* in the manuscript.

Recommended: upload this file (unchanged) to OSF (https://osf.io, free) or a dated GitHub release to obtain a public timestamp, and cite it in the Methods.

# 2. Novelty statement (checked 15 Oct 2026)

| Prior work | What it did | How this study differs |
|---|---|---|
| Tahoe-100M (Zhang et al., Cell 2026; bioRxiv 2025) | Built the atlas; E-distance benchmarks vs sci-Plex/Replogle; MoA used only as a grouping variable (no accuracy reported); drug-response links mentioned only as future work | First quantitative test of Tahoe signatures against **clinical cardiotoxicity**, **external viability (PRISM/GDSC)**, and **baseline RNA/protein**; first MoA-recognition accuracy with permutation nulls |
| Tahoe-x1 (bioRxiv 2025) | Foundation model; essentiality, hallmarks, cell types, perturbation prediction | No drug sensitivity, toxicity or MoA evaluation |
| State (Adduri et al., Cell 2026), Virtual Cell Challenge | Predict perturbation *transcriptomes* | We test what the transcriptomes predict about **phenotypic outcomes** |
| Szalai et al., NAR 2019 (L1000 + CTRP) | Bulk L1000 signatures carry a strong cell-death/viability component | We test whether this holds in **pooled single-cell pseudobulk at 24 h**, within-drug across lines, with proliferation programs pre-specified (C2) — discovery data suggest it does **not** |
| Seal et al., JCIM 2024; DICTrank (Qu 2023) | Cardiotoxicity from chemistry, targets, Cell Painting | We add **tumour-cell transcriptomic responses at atlas scale** and show they add no signal beyond chemistry |
| Gonçalves et al., Cancer Cell 2022; Nusinow et al., Cell 2020 | Baseline proteomics predicts drug response | We compare **post-treatment response signatures** against **baseline RNA, MS proteomics and RPPA (incl. phospho-proteins)** on the same lines and drugs |

Novelty is never absolute: re-run the searches (Google Scholar, bioRxiv: "Tahoe-100M" + toxicity/viability/sensitivity/benchmark) monthly and immediately before submission.

# 3. Discovery results already obtained (exploratory)

* **D1 Cardiotoxicity (FDA DICTrank, 100 clear-label compounds):** chemistry AUROC 0.67; Tahoe tumour response 0.46–0.54 across raw genes, PCA, 14 gene programs, potency; adding Tahoe to chemistry lowered AUROC (0.57–0.64). Same for lung-only lines, Most+Less vs No (n = 167), anticancer subset (n = 43), and DICTrank subtypes (heart damage / arrhythmia / mixed).
* **D2 Positive controls (pipeline valid):** same-drug mid vs high dose signature correlation 0.68 vs 0.23 between drugs; anticancer MoA recognition 1-NN balanced accuracy 0.32 vs chance 0.09 (permutation 95th percentile 0.15).
* **D3 Tahoe cell counts are not a viability readout at 24 h:** raw and pool-share measures do not agree with PRISM or GDSC (ρ ≈ 0), poor dose consistency, salt-form replicates disagree.
* **D4 Response magnitude does not predict sensitivity:** within-drug median ρ −0.02 (PRISM, 51 drugs; p ≈ 0.3) and −0.02 to −0.08 (GDSC, 43 drugs; one nominal p = 0.01 not surviving correction, not replicated in PRISM).

# 4. Confirmatory analyses (frozen)

## C1 — Independent cardiotoxicity labels (replication of D1)
* **Labels:** FAERS cardiac reporting-odds-ratio signal (openFDA; lower 95% CI > 1 and ≥ 3 reports; MedDRA terms fixed in script 04), restricted to FDA-approved Tahoe drugs **not used for model selection**.
* **Features / models:** identical to D1 (chemistry; Tahoe PCA-10; chemistry + Tahoe), logistic regression, repeated 5-fold CV (5 × 10).
* **Primary metric:** AUROC difference (chemistry + Tahoe) − (chemistry), 95% CI by bootstrap over drugs (2,000 resamples).
* **Decision rule:** Tahoe "adds information" only if the lower 95% CI of the difference is > 0.

## C2 — Proliferation / cell-death program test (Szalai hypothesis)
* **Program (fixed list, 20 genes):** MKI67, TOP2A, CCNB1, CCNB2, CDK1, CCNA2, BUB1, PLK1, AURKA, AURKB, CENPF, E2F1, MCM2, MCM3, MCM4, MCM5, MCM6, MCM7, PCNA, TYMS.
* **Score:** mean top-dose log2FC of the program genes per drug × line.
* **Test:** within each drug (≥ 20 labelled lines), Spearman ρ between score and AUC (PRISM and GDSC separately). Expected under the hypothesis: **ρ > 0** (lines whose proliferation genes fall most are the most sensitive, i.e., low AUC).
* **Statistic:** median ρ across drugs; one-sided permutation p (1,000 permutations of lines within drug).
* **Decision rule:** supported if p < 0.025 in **both** PRISM and GDSC (Bonferroni for 2 sources).

## C3 — What predicts which cell lines a drug kills? Response vs baseline RNA vs baseline protein
* **Unit:** drug with ≥ 20 cell lines having the label (PRISM or GDSC) **and** all modalities below ("common-line set").
* **Feature sets (each reduced by PCA to 5 components fitted inside every training fold):**
  1. **Response:** Tahoe top-dose log2FC (5,000 most-expressed genes, fixed list from the data, not from labels)
  2. **Baseline RNA:** CCLE RNA-seq TPM (log2(TPM+1)), same gene universe
  3. **Baseline protein (MS):** CCLE proteomics (Nusinow 2020)
  4. **Baseline protein (RPPA):** CCLE RPPA (≈ 200 antibodies incl. phospho)
  5. **Response + baseline RNA**, 6. **Response + baseline protein (MS)**
* **Model:** ridge regression (standardised features, α = 1.0, fixed — no tuning), leave-one-line-out within drug.
* **Metric:** per-drug Spearman ρ (predicted vs observed AUC); summary = median across drugs.
* **Tests:** (a) each set vs null — permutation of AUC within drug (500 permutations), one-sided p; (b) pairwise comparisons Response vs RNA, Response vs MS protein, MS protein vs RNA, RNA+Response vs RNA — two-sided Wilcoxon signed-rank across drugs.
* **Multiple testing:** Benjamini–Hochberg across all C3 tests, FDR 5%.
* **Decision rules:** a set "predicts sensitivity" if FDR < 0.05 and median ρ > 0.10. Response "adds to" a baseline if the combined set beats that baseline at FDR < 0.05.
* **Lung focus (descriptive only, too few lines for inference):** for the 15 lung lines, report predicted vs observed rankings for EGFR/MEK/KRAS-relevant drugs.

## Robustness (pre-specified, reported regardless of outcome)
* PCA components 3 and 10 instead of 5; α = 0.1 and 10; middle dose instead of top dose; PRISM vs GDSC agreement on shared pairs.

# 5. Manuscript outline

**Title (working):** *Mechanism without consequence: a pre-specified evaluation of the Tahoe-100M single-cell perturbation atlas for predicting drug cardiotoxicity and cancer-cell sensitivity*

**Abstract** (structured, ≤ 250 words): Background (atlases + virtual cells) → Methods (Tahoe-100M, DICTrank/FAERS, PRISM/GDSC, CCLE RNA/MS/RPPA; pre-specified plan) → Results (D1–D4, C1–C3) → Conclusion (what atlas signatures can/cannot be used for).

**Introduction:** rise of perturbation atlases and virtual-cell models; untested assumption that perturbation transcriptomes inform clinical and phenotypic outcomes; gap (no systematic external evaluation of Tahoe-100M); aims.

**Results**
1. Tahoe signatures are reproducible and encode mechanism of action (D2; Fig. 1).
2. Tumour-cell signatures do not predict clinical cardiotoxicity beyond chemistry (D1 + C1; Fig. 2).
3. 24-h pooled single-cell counts do not measure drug killing (D3; Fig. 3a).
4. Neither response magnitude nor proliferation programs predict which lines are killed (D4 + C2; Fig. 3b–c) — contrast with L1000 (Szalai 2019).
5. Response signatures vs baseline RNA vs baseline protein for sensitivity prediction (C3; Fig. 4).
6. Lung-cancer lines: descriptive case studies (Fig. 5).

**Discussion:** what atlases are good for (MoA, target deconvolution) vs not (organ toxicity, viability at 24 h); implications for virtual-cell evaluation (phenotype-anchored benchmarks); why bulk L1000 differs (time point, pooling, cell-number normalisation); limitations (cell lines only, label noise, 24-h single time point, label sources' platforms); recommendations for future atlases (include viability readouts, later time points, cardiac contexts).

**Methods:** data sources and versions (dates accessed); name harmonisation (salt/synonym rules); feature construction; models; CV; statistics; pre-registration; code and data availability (GitHub + Zenodo).

**Target journals (computational):** Genome Biology; Briefings in Bioinformatics; NAR Genomics and Bioinformatics; Bioinformatics; Patterns. Ethics: public de-identified data only (TMU-JIRB exemption letter).

# 6. Data for C3 (all open, no account)

| Data | File | URL |
|---|---|---|
| CCLE RNA-seq (Ghandi 2019) | CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz (63 MB) | https://data.broadinstitute.org/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz |
| CCLE RPPA | CCLE_RPPA_20181003.csv (1.9 MB) | https://data.broadinstitute.org/ccle/CCLE_RPPA_20181003.csv |
| CCLE MS proteomics (Nusinow 2020) | protein_quant_current_normalized.csv.gz | already downloaded |
| PRISM 20Q2 secondary, GDSC2 | — | already downloaded |

# 7. Amendment 1 (added 15 Oct 2026, AFTER the confirmatory C2/C3 results were seen)

**What happened.** The frozen C3 used leave-one-line-out (LOO) ridge regression with a training-mean intercept. Under the null hypothesis this procedure is negatively biased (leaving one line out shifts the training mean in the opposite direction), and the size of the bias differs between feature sets (median null ρ: Response −0.38; RNA −0.15 to −0.16; MS protein −0.15 to −0.16; RPPA −0.16 to −0.17).

**Consequences for the frozen analysis.** (i) The permutation p-values remain valid, because every permutation passes through the identical LOO procedure. (ii) The effect-size rule "median ρ > 0.10" is not interpretable on raw ρ. (iii) The frozen pairwise Wilcoxon comparisons on raw ρ are confounded by set-specific bias; e.g., the frozen result "Response worse than RNA / protein in PRISM (FDR < 0.001)" is an artefact.

**Amended analysis (script 09b).** Identical data, model, seed and 500 permutations; per drug and set we report calibrated ρ = observed ρ − mean permutation ρ; sets are tested for calibrated ρ > 0 (one-sided Wilcoxon across drugs) and compared pairwise on calibrated ρ (two-sided Wilcoxon), with BH-FDR over all amended tests. **Frozen results are reported unchanged in the manuscript, with the amendment clearly labelled as post hoc.**

**Amendment 1 results (15 Oct 2026).** After calibration, the frozen "Response worse than RNA/protein" difference disappears (PRISM p = 0.78 and 0.37). Baseline protein shows small, consistent positive calibrated effects in both screens (MS: +0.07 GDSC, +0.13 PRISM; RPPA: +0.17 GDSC, FDR 0.039, +0.10 PRISM, p = 0.059), baseline RNA shows none (+0.003, +0.04), and RPPA outperforms RNA in GDSC (FDR 0.039). The Tahoe 24-h response adds no reproducible information (GDSC +0.10, PRISM −0.05).

# v2 extension — results summary (18 Oct 2026)

All tests pre-specified in Amendment 3 (frozen 18 Oct 2026, before running). Calibrated ρ = observed − mean permutation ρ.

## Tier 2 (independent data; confirmatory)

| Test | Result | Frozen decision |
|---|---|---|
| **T2-E CTRP replication** (46 drugs, median 27 lines) | RPPA +0.180 (p < 1e-4); RPPA − RNA −0.073 (one-sided p 0.999); RNA +0.225, MS protein +0.182, Response +0.102; Response < RNA (FDR 0.040); RNA > RPPA (FDR 0.0035) | GDSC RPPA finding **NOT replicated** (E2 fails) |
| **T2-F ProCan proteome** (3,594 proteins; PRISM 44, GDSC 37, CTRP 47 drugs) | ProCan: PRISM +0.046 (p 0.21), GDSC +0.293 (p < 1e-4), CTRP +0.282 (p 1e-4); ProCan − RNA: −0.085 (p 0.86), **+0.180 (p 0.0002)**, −0.012 (p 0.63) | "protein > RNA" **NOT supported** (1 of 3 screens) |
| **T2-H1a cardiomyocyte identity** (25 drugs, 5,518 genes) | same-drug r +0.062 vs different +0.013; permutation p = 0.0001 (Spearman and drug-selective DEGs: p ≤ 0.0002) | Identity **preserved** across contexts |
| T2-H1b | 25 shared labelled drugs (< 30) | descriptive only |
| **T2-H2 hERG** (ChEMBL; 43 drugs, 28 with pIC50 ≥ 5) | Spearman chem 0.027, chem+Tahoe 0.002; Δ −0.025 (95% CI −0.087 to +0.018); AUROC 0.512 vs 0.507 | Tahoe **does not add** (weak test: chemistry also ≈ chance at n = 43) |
| T2-G L1000 | not performed (data volume; dropped before download) | — |

## Tier 1 (existing data; secondary)

| Test | Result |
|---|---|
| **T1-D middle dose** (Amendment-1 procedure, same 43/37 drugs) | Baselines identical (check). Response: PRISM +0.074 (FDR 0.019; top dose −0.046), GDSC +0.070 (top +0.102). Response never adds to RNA or protein (all pairwise FDR ≥ 0.33). |
| **T1-C RPPA–target concordance** (46 drug × screen units) | median 0.511 vs null 0.496; p = 0.148 → **not supported**. Descriptively EGFR/ERBB inhibitors top (gefitinib 0.89, afatinib 0.83, osimertinib 0.77, erlotinib 0.74 in GDSC; tucatinib 0.72 in PRISM). |
| **T1-A2 best-effort sensitivity** (4 families: ridge-GCV, RBF kernel ridge, elastic net, random forest; 5-fold CV, PCA-10, 50 perms) | Response rescued: **NO** (no family improves on RNA or protein in either screen). RPPA − RNA in GDSC: −0.026, −0.007, −0.028, −0.101 → positive in **0/4** families → protein advantage **not model-robust**. Baseline RNA is positive in GDSC for every family (+0.07 to +0.12). |
| **T1-A1 best-effort cardiotoxicity** (4 families × 2 Tahoe representations, nested CV, Holm over 8) | C1 FAERS: Tahoe **not rescued** (all Holm p ≥ 0.16; every model AUROC 0.35–0.55). D1 DICTrank: chemistry 0.60–0.69, Tahoe 0.29–0.56 (none above permutation null), adding Tahoe lowers AUROC (Δ down to −0.32). |
| **T1-B LOO-bias simulation** | Done — see section below. |

## What the v2 paper can now claim (all pre-registered analyses complete, 19 Oct 2026)
1. **Baseline cell state predicts within-drug sensitivity** across three independent screens once the leave-one-out bias is removed: RNA (CTRP +0.22), MS protein (CCLE/ProCan; GDSC +0.29, CTRP +0.28), RPPA (CTRP +0.18, GDSC +0.17).
2. **No baseline modality is consistently best.** Protein > RNA appears only in GDSC (with both RPPA and ProCan), not in PRISM or CTRP. Hypothesis for the Discussion (not tested): GDSC and ProCan were both generated at the Sanger Institute on shared cell stocks.
3. **The 24-h Tahoe response is weaker than baseline state and never adds to it**: top and middle dose, PRISM/GDSC/CTRP, multiple model families.
4. **Cardiotoxicity:** Tahoe adds nothing to chemistry (DICTrank, FAERS, hERG), yet drug identity is preserved between tumour-line and iPSC-cardiomyocyte responses (H1a). Toxicity-relevant biology appears context-specific rather than absent.
5. Methods contribution: leave-one-out within-drug evaluation is strongly negatively biased under the null (−0.77 at 27 lines); permutation calibration gives a valid test statistic but not a correlation-scale effect size.
6. Best-effort models (4 families) confirm both nulls (cardiotoxicity, response added value) and show the GDSC protein advantage is model-dependent.

## Files
`v2/results/T2EF/` (D, E, F), `T2H1/`, `T2H2/`, `T1C/`, `T1A2/`; scripts in `v2/scripts/`; Amendment 3 with change log in `v2/`.

## Sanity checks (19 Oct 2026)
* CTRP AUC direction: Spearman with GDSC AUC +0.526 (2,738 drug–line pairs, 89 drugs) and with PRISM +0.553 (1,579 pairs) → same direction (lower = more sensitive).
* Cross-screen agreement within drugs (across lines): median ρ 0.267 (CTRP–GDSC) and 0.205 (CTRP–PRISM) → practical ceiling for within-drug prediction. NOTE (T1-B): calibrated ρ is a valid test statistic but not on the correlation scale (LOO calibration can overshoot the true ρ), so it should not be compared numerically with this ceiling.
* ProCan vs CCLE MS proteome: median per-protein Spearman across 29 shared lines +0.448 (3,249 proteins).
* Key Wilcoxon tests recomputed independently from per-drug files: identical to the summaries.

## T1-B leave-one-out bias simulation (19 Oct 2026)
* Null bias of the frozen procedure (k = 5, α = 1): mean LOO ρ −0.94 (n = 10), −0.84 (20), −0.77 (27), −0.64 (40), −0.51 (60); intercept-only predictor exactly −1.
* Permutation calibration centres the null at 0 for LOO and 5-fold CV at every n (|mean| ≤ 0.04), i.e. calibrated ρ is a valid test statistic.
* With signal, calibrated ρ increases monotonically with the true correlation but is not on the correlation scale: LOO overshoots at larger n (n = 40, true 0.37 → 0.54; n = 60, 0.38 → 0.59) and undershoots at small n; 5-fold CV is closer (n = 40: 0.29; n = 60: 0.37).
* Recommendation: report calibrated ρ for testing and ranking feature sets; do not interpret its magnitude as a correlation; prefer K-fold over LOO for effect sizes.

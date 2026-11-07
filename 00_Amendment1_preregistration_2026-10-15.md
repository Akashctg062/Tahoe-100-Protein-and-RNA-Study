---
title: "Amendment 1 — Null-calibrated C3 (post hoc)"
study: "What does a 100-million-cell drug-perturbation atlas predict? Pre-specified evaluation of Tahoe-100M"
---

This file reproduces, unchanged, Section 7 of the analysis plan as amended
(`preregistration/03_Prespecified_Analysis_Plan_and_Manuscript_Outline_v3_with_Amendments1-2.md`). The frozen plan is `preregistration/03_Prespecified_Analysis_Plan_and_Manuscript_Outline.md`;
the plan with this amendment alone is `preregistration/03_Prespecified_Analysis_Plan_and_Manuscript_Outline_v2_with_Amendment1.md`.

# 7. Amendment 1 (AFTER the confirmatory C2/C3 results were seen)

**What happened.** The frozen C3 used leave-one-line-out (LOO) ridge regression with a training-mean intercept. Under the null hypothesis this procedure is negatively biased (leaving one line out shifts the training mean in the opposite direction), and the size of the bias differs between feature sets (median null ρ: Response −0.38; RNA −0.15 to −0.16; MS protein −0.15 to −0.16; RPPA −0.16 to −0.17).

**Consequences for the frozen analysis.** (i) The permutation p-values remain valid, because every permutation passes through the identical LOO procedure. (ii) The effect-size rule "median ρ > 0.10" is not interpretable on raw ρ. (iii) The frozen pairwise Wilcoxon comparisons on raw ρ are confounded by set-specific bias; e.g., the frozen result "Response worse than RNA / protein in PRISM (FDR < 0.001)" is an artefact.

**Amended analysis (script 09b).** Identical data, model, seed and 500 permutations; per drug and set we report calibrated ρ = observed ρ − mean permutation ρ; sets are tested for calibrated ρ > 0 (one-sided Wilcoxon across drugs) and compared pairwise on calibrated ρ (two-sided Wilcoxon), with BH-FDR over all amended tests. **Frozen results are reported unchanged in the manuscript, with the amendment clearly labelled as post hoc.**

**Amendment 1 results (15 Oct 2026).** After calibration, the frozen "Response worse than RNA/protein" difference disappears (PRISM p = 0.78 and 0.37). Baseline protein shows small, consistent positive calibrated effects in both screens (MS: +0.07 GDSC, +0.13 PRISM; RPPA: +0.17 GDSC, FDR 0.039, +0.10 PRISM, p = 0.059), baseline RNA shows none (+0.003, +0.04), and RPPA outperforms RNA in GDSC (FDR 0.039). The Tahoe 24-h response adds no reproducible information (GDSC +0.10, PRISM −0.05).

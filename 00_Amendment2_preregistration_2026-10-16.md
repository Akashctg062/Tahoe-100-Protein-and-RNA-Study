---
title: "Amendment 2 — C1 implementation note and sensitivity analysis excluding cardiovascular-indication drugs (post hoc)"
study: "What does a 100-million-cell drug-perturbation atlas predict? Pre-specified evaluation of Tahoe-100M"
---

This file reproduces, unchanged, Section 8 of the analysis plan as amended
(`preregistration/03_Prespecified_Analysis_Plan_and_Manuscript_Outline_v3_with_Amendments1-2.md`). The frozen plan is `preregistration/03_Prespecified_Analysis_Plan_and_Manuscript_Outline.md`.

# 8. C1 implementation note and Amendment 2

**Implementation correction made BEFORE any C1 result was computed.** The first openFDA search used only `patient.drug.openfda.generic_name`, which returned 0 reports for 37 of 148 drugs, including FDA-approved abemaciclib, darolutamide and belzutifan (missing harmonisation for recent drugs). The search was changed to `generic_name OR activesubstance.activesubstancename`; zero-report drugs fell to 3. Only this v2 search was used for C1.

**C1 frozen result.** 128 drugs (66 FAERS cardiac signal). AUROC chemistry 0.459 (95% CI 0.360–0.558); Tahoe 0.444 (0.346–0.543); chemistry + Tahoe 0.456 (0.357–0.554). Primary difference −0.003 (95% CI −0.012 to +0.006) → **Tahoe does not add information** (frozen decision rule).

**Amendment 2 (post hoc, labelled as such).** The strongest FAERS signals were cardiovascular drugs (eplerenone, sacubitril/valsartan, furosemide, captopril, clopidogrel …), i.e. confounding by indication. Excluding 13 cardiovascular-indication drugs (signal rate 0.92 vs 0.47) and recomputing AUROC from the frozen out-of-fold predictions: chemistry 0.487 (0.382–0.596), Tahoe 0.432 (0.327–0.538), chemistry + Tahoe 0.484 (0.372–0.590). Conclusion unchanged. Interpretation: FAERS cardiac signal is a weak, indication-confounded label; C1 is consistent with D1 but is a weak replication on its own.

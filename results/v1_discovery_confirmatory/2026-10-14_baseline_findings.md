# Baseline findings — 14 Oct 2026

## Data
- Tahoe-100M pseudobulk DE, 196 drugs × 3 doses; summarised over all 50 lines and over 27 lung lines (script 03b).
- FDA DICTrank labels: 179 compounds matched with features; 100 clear (65 most-concern / 35 no-concern); 43 of the 100 are anticancer (28 / 15).

## Main result: tumour responses do NOT predict DICTrank cardiotoxicity
Repeated 5-fold CV, AUROC mean ± SD (logistic regression unless stated)

| Task | Structure | Tumour raw 2000 genes | Tumour PCA-10 | Tumour 14 gene programs | Structure + programs |
|---|---|---|---|---|---|
| Most vs No (n=100) | **0.67 ± 0.11** | 0.48 ± 0.09 | 0.46 ± 0.11 | 0.53 ± 0.14 | 0.64 ± 0.11 |
| Most+Less vs No (n=167) | **0.59 ± 0.09** | 0.44 ± 0.10 | 0.48 ± 0.09 | 0.54 ± 0.11 | 0.57 ± 0.10 |
| Anticancer only, Most vs No (n=43) | 0.45 ± 0.24 | 0.43 ± 0.21 | 0.50 ± 0.17 | 0.34 ± 0.15 | 0.43 ± 0.22 |
| Lung lines only, Most vs No (n=100) | 0.67 | 0.46 | 0.46 | 0.53 | 0.65 |

Subtype test (pre-specified): heart damage vs no — structure 0.53, tumour 0.33; arrhythmia vs no — 0.51 vs 0.53; mixed vs no — 0.63 vs 0.44. No subtype is predicted by tumour responses.

## Positive controls: the pipeline is sound
- Same drug, mid vs high dose signature correlation: median 0.68 (different drugs: 0.23).
- Anticancer MoA recognition (1-NN, correlation, 11 classes): balanced accuracy 0.32 vs chance 0.09 (permutation 95th pct 0.15). MEK 3/3, HDAC 2/3, JAK 2/3 correct.

## Interpretation
Tahoe tumour-cell transcriptional responses encode drug mechanism of action robustly, but this information does not transfer to clinical cardiotoxicity labels. Cardiac liability appears to require cardiac-context data (cardiomyocyte responses, cardiac ion-channel/target pharmacology), not tumour-cell context.

Scripts: scripts/exploratory/ (prep.py, exp.py, posctrl2.py, subtype.py) — run in the cloud workspace on a copy of drug_gene_features.parquet.

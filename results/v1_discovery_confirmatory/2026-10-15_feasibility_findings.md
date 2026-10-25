# Feasibility findings — 15 Oct 2026 (lung drug-response study, Option 2)

## Corrected cell-line numbers
- Tahoe cell_line_metadata.csv lists 102 lines, but only **50 lines were measured** in Tahoe-100M.
- **15 measured lung lines**: A-427, A549, HOP62, NCI-H1573, NCI-H1792, NCI-H2030, NCI-H2122, NCI-H23, NCI-H2347, NCI-H460, NCI-H596, NCI-H661, SHP-77, SW 1271, SW 900.
- (Earlier "27 lung lines" counted metadata, not measured lines. Baseline results are unaffected.)

## Proteomics coverage (CCLE, Nusinow 2020)
- 33 / 50 measured lines; **10 / 15 lung lines** (A549, NCI-H1573, NCI-H1792, NCI-H2030, NCI-H2122, NCI-H23, NCI-H460, NCI-H661, SHP-77, SW 1271).

## Label overlap (198-drug subset; all-drug extraction now complete, numbers will grow)
- PRISM secondary (MTS010): 1,639 drug-line pairs, 45 drugs, 43 lines (lung 543 pairs, 14 lines).
- GDSC2: 1,650 pairs, 39 drugs, 44 lines (lung 513 pairs, 14 lines).

## Tahoe cell-count "survival" is NOT a valid viability label (idea dropped)
- Raw ratio (treated/DMSO): agreement with PRISM rho 0.03, GDSC rho 0.09; no dose monotonicity (51%).
- Pool-share corrected measure (line's share of the well vs DMSO share): dose consistency median rho 0.31,
  salt-variant replicates inconsistent, agreement with PRISM 0.01 / GDSC 0.00.
- Interpretation: 24 h is too early and pooled-well cell counts too noisy to measure killing.

## Next decisive test
- Script 08: does each line's 24-h transcriptional response strength predict its PRISM/GDSC sensitivity
  (within drug, across lines; permutation test; cell-count adjusted)?

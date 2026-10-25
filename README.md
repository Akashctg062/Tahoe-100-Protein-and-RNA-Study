# What does a 100-million-cell drug-perturbation atlas predict?

Code, pre-registration documents and result files for a pre-registered, fully computational evaluation of the
**Tahoe-100M** single-cell drug-perturbation atlas against clinical cardiotoxicity, hERG potency,
iPSC-cardiomyocyte transcriptomes, three cancer-cell drug-sensitivity screens (PRISM, GDSC, CTRP) and baseline
RNA, mass-spectrometry protein (CCLE, ProCan) and RPPA profiles.

Author: Asaduzzaman, Graduate Institute of Biomedical Materials and Tissue Engineering, Taipei Medical University.
Manuscript in preparation — please cite the article once it is published (citation will be added here).

---

## Main findings

| Question | Answer |
|---|---|
| Are Tahoe signatures reproducible and mechanism-specific? | Yes. Same drug at two doses r = 0.68 vs 0.22 between drugs; mechanism of action recognised (balanced accuracy 0.32 vs chance 0.09). |
| Is drug identity preserved between tumour cells and iPSC-cardiomyocytes? | Yes. Same-drug r = 0.062 vs 0.013 between drugs; permutation p = 0.0001 (25 drugs). |
| Do Tahoe signatures add to chemistry for cardiotoxicity (DICTrank, FAERS) or hERG potency? | No. This holds for the frozen models and for four flexible model families (all Holm p ≥ 0.16). |
| Is leave-one-out within-drug evaluation biased? | Yes, strongly negatively (mean null ρ −0.77 at 27 lines). Permutation calibration restores a valid test, but the calibrated value is not an effect size. |
| Does baseline cell state predict which lines a drug kills? | Yes, in PRISM, GDSC and CTRP. No single data type is consistently best. |
| Does the 24-h Tahoe response add to baseline state? | No. Tested at two doses, in three screens and with five model families. |

---

## Pre-registration timeline

| Date | Document | Status |
|---|---|---|
| 15 Oct 2026 | `preregistration/03_Prespecified_Analysis_Plan_and_Manuscript_Outline.md` | Frozen plan (before confirmatory analyses C1–C3) |
| 15 Oct 2026 | `preregistration/..._v2_with_Amendment1.md` | Amendment 1 (post hoc): null-calibrated C3 |
| 16 Oct 2026 | `preregistration/..._v3_with_Amendments1-2.md` | Amendment 2 (post hoc): C1 without cardiovascular-indication drugs |
| 18 Oct 2026 | `00_Amendment3_preregistration_2026-10-18.md` | Amendment 3: extension analyses T1-A–T1-D and T2-E–T2-H, frozen **before** they were run and before the new data were downloaded; includes a dated change log |

Analysis labels: **D** = discovery, **C** = confirmatory (frozen), **T** = pre-registered extension (Amendment 3).
Exploratory analyses are labelled as such.

---

## Repository structure

```
preregistration/              frozen plan and Amendments 1-2 (Amendment 3 is in the repository root)
scripts/
  01_core_analysis/           steps 01-11: data download, name harmonisation, D1-D4, C1-C3, Amendments 1-2, lung subset
  02_exploratory/             discovery helpers (D1 models and subtypes, D2 positive controls)
  03_extension_v2/            Amendment 3 analyses (T1-A to T1-D, T2-E, T2-F, T2-H) and data preparation
  04_figures_v1/              figure scripts, original analyses
  05_figures_v2/              figure scripts, extended analyses
  06_checks/                  independent recomputation of all numbers shown in the v2 figures
results/
  v1_discovery_confirmatory/  D1-D4, C1-C3, Amendments 1-2, robustness, lung subset (per-drug tables and summaries)
  v2_extension/               T1A1, T1A2, T1B, T1C, T2EF (T1-D, T2-E, T2-F), T2H1, T2H2 + results summary
  checks/                     output of the independent number check (40 checks, all matching)
figures/
  v1/                         Figures 1-5 and S1 of the original analysis (PNG, PDF) + source data
  v2/                         Figures 1-7 and S1-S3 of the extended analysis (PNG, PDF) + source data
requirements.txt              Python packages used
environment.yml               original conda environment file (includes some packages not used in the final analyses)
```

### Script order

**Original analyses (`scripts/01_core_analysis/`)**

1. `01_explore_tahoe_metadata.py`: download Tahoe-100M metadata over plain HTTPS
2. `02_drug_overlap.py`: drug-name harmonisation and label overlap (DICTrank)
3. `03_extract_tahoe_pseudobulk.py`, `03b_build_drug_features.py`: filter the 1,026 pseudobulk differential-expression files and build drug-level features
4. `04_faers_cardiac_ror.py`, `06_temporal_faers_validation.py`: FAERS reporting odds ratios (discovery)
5. `05_baseline_models.py` (+ `02_exploratory/`): D1 cardiotoxicity models, D2 positive controls
6. `07_feasibility_lung_response.py`, `08_response_vs_sensitivity.py`: feasibility, D3 and D4
7. `09_confirmatory_C2_C3.py`: frozen C2 and C3
8. `09b_amendment_calibrated_C3.py`, `09c_run_robustness.py`: Amendment 1 and the robustness grid
9. `10_confirmatory_C1_faers.py`, `10b_C1_sensitivity_cv_indication.py`: frozen C1 and Amendment 2
10. `11_exploratory_lung_subset.py`: exploratory lung-subset evaluation

**Extension (`scripts/03_extension_v2/`)**

1. `download_tier2.ps1`, `T2H2_download_chembl_herg.py`: download CTRP v2 / ProCan (DrEval release) and ChEMBL hERG data
2. `T2_prepare_ctrp_procan.py`, `T2_build_responses.py`, `T2_build_meanstat_signatures.py`, `T2_build_alldrug_signatures.py`, `T1_build_c3_cache.py`: data preparation
3. `T1A1_best_effort_cardiotox.py`, `T1A2_best_effort_sensitivity.py`, `T1B_loo_bias_simulation.py`, `T1C_rppa_target_concordance.py`
4. `T2EF_ctrp_procan.py`: T1-D (middle dose), T2-E (CTRP) and T2-F (ProCan). `T1D_middle_dose.py` is an earlier version superseded by part D of this script.
5. `T2H1_cardiomyocyte_context.py`, `T2H2_herg.py`

**Note on paths.** The scripts were run partly on a Windows computer and partly in a Linux cloud workspace. Some
scripts therefore contain absolute paths at the top (for example `U = ...`, `OUT = ...`). Edit these to point
to your copy of the data before running.

---

## Data sources

All input data are public and de-identified. They are **not** redistributed here; the scripts download or read them.

| Data | Source |
|---|---|
| Tahoe-100M (metadata, pseudobulk differential expression) | https://huggingface.co/datasets/tahoebio/Tahoe-100M (downloaded over plain HTTPS; no Hugging Face library) |
| FDA DICTrank | https://www.fda.gov/science-research/bioinformatics-tools/drug-induced-cardiotoxicity-rank-dictrank-dataset |
| FDA FAERS (openFDA drug-event API) | https://open.fda.gov/apis/drug/event/ (queried 16 Oct 2026) |
| ChEMBL hERG (CHEMBL240) | https://www.ebi.ac.uk/chembl/ |
| iPSC-cardiomyocyte signatures (Hansen et al. 2024, Suppl. Data 05 and 10) | https://doi.org/10.1038/s41467-024-52145-4 |
| PRISM 20Q2, CCLE RNA-seq (2018-09-29), CCLE MS proteomics (Nusinow 2020), CCLE RPPA (2018-10-03) | https://depmap.org/portal/ |
| GDSC2 (release 24 Jul 2022) | https://www.cancerrxgene.org/downloads/bulk_download |
| CTRP v2 and ProCan-DepMapSanger (DrEval release) | https://doi.org/10.5281/zenodo.17600096 |

**Not included because of file size** (GitHub limit 100 MB): the derived Tahoe response tables and signature
arrays (`response_dosetop.parquet`, `response_top_80drugs.parquet`, `response_mid_80drugs.parquet`,
`tahoe_features_C1.npz`, `tahoe_signatures_remaining.npz`). They are regenerated from the public data by steps
03/03b, `10_confirmatory_C1_faers.py` and the `T2_build_*` scripts.

**openFDA API key.** The FAERS scripts read an optional key from the environment variable `OPENFDA_KEY`.
No key is stored in this repository. Never commit a key.

---

## Installation

```bash
python -m pip install -r requirements.txt
```

The analyses were run with Python 3.11–3.13. Random seeds are fixed in every script (2026 unless stated).

## License

[Choose a license, for example MIT for the code, and add a LICENSE file.]

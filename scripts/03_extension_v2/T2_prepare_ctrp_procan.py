"""
v2 / Tier 2 data preparation (Amendment 3, change-log note of 18 Oct 2026).
Extracts from the drevalpy CTRPv2 package (Zenodo 10.5281/zenodo.17600096):
  (1) CTRP v2 drug response (CurveCurator re-fit; AUC_curvecurator) -> ctrp_auc.csv (cell line, drug, AUC)
  (2) ProCan-DepMapSanger replicate-level protein intensities (Goncalves et al. 2022) -> procan_log2_mean.csv
      processed as published: log2(intensity), replicates of a cell line combined by geometric mean (= mean of log2).
Run (any Python with pandas):  python T2_prepare_ctrp_procan.py
Writes to  E:\\Protein and RNA\\data\\raw\\v2\\derived\\   (raw zip is read only)
"""
import os, zipfile
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.abspath(os.path.join(HERE, "..", "..", "data", "raw", "v2"))
OUT = os.path.join(RAW, "derived"); os.makedirs(OUT, exist_ok=True)
z = zipfile.ZipFile(os.path.join(RAW, "ctrp", "drevalpy_CTRPv2.zip"))

with z.open("CTRPv2/CTRPv2.csv") as f:
    d = pd.read_csv(f, usecols=["cellosaurus_id", "cell_line_name", "pubchem_id", "drug_name", "AUC_curvecurator",
                                "IC50_curvecurator", "R2"], low_memory=False)
d.to_csv(os.path.join(OUT, "ctrp_auc.csv"), index=False)
print("CTRP rows:", len(d), "| drugs:", d.drug_name.nunique(), "| lines:", d.cellosaurus_id.nunique())
for nm in ["cell_line_names.csv", "drug_names.csv", "drug_smiles.csv"]:
    with z.open("CTRPv2/" + nm) as f: pd.read_csv(f).to_csv(os.path.join(OUT, "ctrp_" + nm), index=False)

with z.open("CTRPv2/proteomics.csv") as f:
    p = pd.read_csv(f, low_memory=False)
meta = p[["cellosaurus_id", "cell_line_name"]]
vals = p.drop(columns=["cellosaurus_id", "cell_line_name"]).apply(pd.to_numeric, errors="coerce")
neg = int((vals <= 0).sum().sum())
L = np.log2(vals.where(vals > 0))
L.insert(0, "cellosaurus_id", meta.cellosaurus_id.values); L.insert(1, "cell_line_name", meta.cell_line_name.values)
A = L.groupby(["cellosaurus_id", "cell_line_name"], as_index=False).mean(numeric_only=True)
A.to_csv(os.path.join(OUT, "procan_log2_mean.csv"), index=False)
print(f"ProCan replicate rows: {len(p)} -> cell lines: {len(A)} | proteins: {A.shape[1]-2} | non-positive values set missing: {neg}")
print("replicates per line (median):", int(L.groupby("cellosaurus_id").size().median()))

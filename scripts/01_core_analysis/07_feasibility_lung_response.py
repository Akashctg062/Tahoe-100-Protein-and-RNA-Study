"""
Step 7 - FEASIBILITY CHECK for the lung-cancer drug-response study (Option 2).

Answers four questions:
  Q1. How many Tahoe cell lines (all 50 / 27 lung) have CCLE proteomics (Nusinow 2020)?
  Q2. How many Tahoe drug x cell-line pairs also have PRISM (secondary, 20Q2) or GDSC2 viability?
  Q3. Does Tahoe's built-in survival measure (treated vs DMSO cell counts) agree with PRISM/GDSC?
  Q4. Summary table saved for the manuscript's data-overview figure.

Inputs (download first - see the guide / chat for copy-paste commands):
  ../data/processed/pseudobulk_parts/*.parquet                 (step 3; already on your disk)
  ../data/raw/tahoe_metadata/cell_line_metadata.csv, drug_metadata.csv   (step 1)
  ../data/raw/depmap/GDSC2_fitted_dose_response_24Jul22.csv
  ../data/raw/depmap/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv
  ../data/raw/depmap/protein_quant_current_normalized.csv.gz

Run:  python 07_feasibility_lung_response.py
Outputs (../results/feasibility/):
  tahoe_viability.csv          drug x dose x cell line survival measure
  pairs_with_labels.csv        Tahoe pairs matched to PRISM / GDSC AUC
  proteomics_coverage.csv      which Tahoe lines have CCLE proteomics
  feasibility_summary.txt      all numbers printed below
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
norm_drug = importlib.import_module("02_drug_overlap").norm          # same drug-name matching as step 2

def norm_cell(x):
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())

def spearman(a, b):
    a, b = pd.Series(a).rank(), pd.Series(b).rank()
    return np.corrcoef(a, b)[0, 1] if len(a) > 2 else np.nan

def main(a):
    os.makedirs(a.out, exist_ok=True)
    log = open(os.path.join(a.out, "feasibility_summary.txt"), "w", encoding="utf-8")
    def say(*x):
        s = " ".join(str(i) for i in x); print(s); log.write(s + "\n")

    cl = pd.read_csv(a.cells)
    cl["lung"] = cl["Organ"].astype(str).str.contains("Lung", case=False)
    lines = cl.drop_duplicates("cell_name")[["cell_name", "Cell_ID_DepMap", "Organ", "lung"]]
    lines["cell_key"] = lines.cell_name.map(norm_cell)
    say(f"Tahoe cell lines: {len(lines)} (lung: {int(lines.lung.sum())})")

    # ---------- Tahoe built-in survival measure ----------
    cols = ["drug", "concentration", "Cell_Name_Vevo", "Cell_ID_DepMap", "n_cells_trt", "n_cells_ctrl"]
    parts = sorted(glob.glob(os.path.join(a.parts, "*.parquet")))
    say(f"Reading cell counts from {len(parts)} filtered Tahoe files ...")
    v = pd.concat([pd.read_parquet(p, columns=cols).drop_duplicates() for p in parts]).drop_duplicates()
    v = v.groupby(["drug", "concentration", "Cell_Name_Vevo", "Cell_ID_DepMap"], as_index=False)[
        ["n_cells_trt", "n_cells_ctrl"]].max()
    # survival = log2 ratio of treated to DMSO cells, centred per cell line on the median drug
    # (most drugs barely affect growth, so the median acts as the 'no effect' reference)
    v["log2_ratio"] = np.log2((v.n_cells_trt + 1) / (v.n_cells_ctrl + 1))
    v["dose_rank"] = v.groupby("drug").concentration.rank(method="dense").astype(int)
    v["survival"] = v.log2_ratio - v.groupby(["Cell_Name_Vevo", "dose_rank"]).log2_ratio.transform("median")
    v["drug_key"] = v.drug.map(norm_drug)
    v.to_csv(os.path.join(a.out, "tahoe_viability.csv"), index=False)
    say(f"Tahoe survival measure: {len(v):,} drug-dose-line rows, {v.drug.nunique()} drugs, "
        f"{v.Cell_Name_Vevo.nunique()} lines")
    top = v[v.dose_rank == v.groupby("drug").dose_rank.transform("max")]
    say(f"  strongly killed (survival < -1, i.e. <50% of typical) at top dose: "
        f"{(top.survival < -1).mean():.1%} of drug-line pairs")

    # ---------- Q1 proteomics coverage ----------
    hdr = pd.read_csv(a.ccle_prot, nrows=0, compression="infer").columns
    samples = [c for c in hdr if re.search(r"_TenPx\d+$", c)]
    prot_lines = {norm_cell(c.split("_")[0]): c for c in samples}
    lines["has_ccle_proteomics"] = lines.cell_key.isin(prot_lines)
    lines.to_csv(os.path.join(a.out, "proteomics_coverage.csv"), index=False)
    say(f"\nQ1. CCLE proteomics: {int(lines.has_ccle_proteomics.sum())}/{len(lines)} Tahoe lines "
        f"(lung: {int(lines[lines.lung].has_ccle_proteomics.sum())}/{int(lines.lung.sum())})")
    say("    lung lines WITH proteomics:", ", ".join(sorted(lines[lines.lung & lines.has_ccle_proteomics].cell_name)))
    say("    lung lines WITHOUT:", ", ".join(sorted(lines[lines.lung & ~lines.has_ccle_proteomics].cell_name)))

    # ---------- Q2 PRISM / GDSC labels ----------
    pr = pd.read_csv(a.prism, usecols=lambda c: c in {"depmap_id", "name", "auc", "ic50", "screen_id", "ccle_name"},
                     low_memory=False)
    pr = pr.dropna(subset=["auc"])
    if "screen_id" in pr and pr.screen_id.astype(str).str.contains("MTS010").any():
        pr = pr[pr.screen_id.astype(str).str.contains("MTS010")]          # recommended secondary screen
    pr["drug_key"] = pr["name"].map(norm_drug)
    pr = pr.groupby(["depmap_id", "drug_key"], as_index=False).auc.mean().rename(
        columns={"depmap_id": "Cell_ID_DepMap", "auc": "prism_auc"})
    gd = pd.read_csv(a.gdsc, usecols=["CELL_LINE_NAME", "DRUG_NAME", "AUC", "LN_IC50"])
    gd["cell_key"] = gd.CELL_LINE_NAME.map(norm_cell); gd["drug_key"] = gd.DRUG_NAME.map(norm_drug)
    gd = gd.groupby(["cell_key", "drug_key"], as_index=False).AUC.mean().rename(columns={"AUC": "gdsc_auc"})
    t = top.merge(lines[["cell_name", "cell_key", "lung"]], left_on="Cell_Name_Vevo", right_on="cell_name")
    t = t.merge(pr, on=["Cell_ID_DepMap", "drug_key"], how="left").merge(gd, on=["cell_key", "drug_key"], how="left")
    t.to_csv(os.path.join(a.out, "pairs_with_labels.csv"), index=False)
    for lab in ["prism_auc", "gdsc_auc"]:
        s = t.dropna(subset=[lab])
        say(f"\nQ2. {lab.split('_')[0].upper()}: {len(s):,} Tahoe drug-line pairs with labels | "
            f"{s.drug.nunique()} drugs | {s.cell_name.nunique()} lines | "
            f"lung: {int(s.lung.sum()):,} pairs in {s[s.lung].cell_name.nunique()} lines")

        # ---------- Q3 agreement ----------
        if len(s) > 20:
            rho_all = spearman(s.survival, s[lab])
            per_line = s.groupby("cell_name").apply(lambda g: spearman(g.survival, g[lab]) if len(g) >= 8 else np.nan).dropna()
            per_drug = s.groupby("drug").apply(lambda g: spearman(g.survival, g[lab]) if len(g) >= 8 else np.nan).dropna()
            say(f"Q3. Agreement Tahoe survival vs {lab}: pooled Spearman rho = {rho_all:.2f}")
            say(f"    within cell line (across drugs): median rho = {per_line.median():.2f} over {len(per_line)} lines")
            say(f"    within drug (across lines):      median rho = {per_drug.median():.2f} over {len(per_drug)} drugs")
            say("    (positive rho = agreement: both low when the line is sensitive)")
    log.close()
    print(f"\nSaved results to {a.out}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="../data/processed/pseudobulk_parts")
    ap.add_argument("--cells", default="../data/raw/tahoe_metadata/cell_line_metadata.csv")
    ap.add_argument("--gdsc", default="../data/raw/depmap/GDSC2_fitted_dose_response_24Jul22.csv")
    ap.add_argument("--prism", default="../data/raw/depmap/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv")
    ap.add_argument("--ccle-prot", dest="ccle_prot", default="../data/raw/depmap/protein_quant_current_normalized.csv.gz")
    ap.add_argument("--out", default="../results/feasibility")
    main(ap.parse_args())

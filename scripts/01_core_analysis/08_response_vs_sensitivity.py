"""
Step 8 - DECISIVE TEST: does a cell line's EARLY (24 h) transcriptional response to a drug
(Tahoe-100M) predict its LATER drug sensitivity (PRISM 5-day / GDSC 72-h viability AUC)?

For every drug x dose x cell line it computes, from the per-line pseudobulk DE results:
  n_de        number of genes with padj < 0.05 and |log2FC| > 0.5
  rms_lfc     root-mean-square log2 fold change over well-expressed genes (baseMean >= 10)
  rms_stat    root-mean-square DE statistic
Then, WITHIN each drug and ACROSS cell lines, it correlates response strength with PRISM/GDSC AUC.
Expected if the hypothesis holds: NEGATIVE correlation (stronger response <-> lower AUC = more sensitive).
A permutation test (shuffling cell lines within each drug) gives a p-value.
Cell-number confounding is checked with a partial correlation that removes log(n_cells_trt).

Run:  python 08_response_vs_sensitivity.py
Inputs:  ../data/processed/pseudobulk_parts_all/  (step 3 --all-drugs)
         ../data/raw/depmap/ PRISM + GDSC files (step 7)
Outputs: ../results/response_vs_sensitivity/
"""
import argparse, glob, os, sys, importlib
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
norm_drug = importlib.import_module("02_drug_overlap").norm

def norm_cell(x):
    import re
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())

def rank(x): return pd.Series(x).rank().values

def partial_spearman(x, y, z):
    """Spearman correlation of x and y after removing the linear effect of z (on ranks)."""
    rx, ry, rz = rank(x), rank(y), rank(z)
    if np.std(rz) == 0:
        return np.corrcoef(rx, ry)[0, 1]
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return np.corrcoef(ex, ey)[0, 1]

def main(a):
    os.makedirs(a.out, exist_ok=True)
    log = open(os.path.join(a.out, "summary.txt"), "w", encoding="utf-8")
    def say(*x):
        s = " ".join(str(i) for i in x); print(s); log.write(s + "\n")

    # ---------- 1. response strength per drug x dose x line ----------
    cache = os.path.join(a.out, "line_response_strength.csv")
    if os.path.exists(cache):
        r = pd.read_csv(cache); say(f"Loaded cached response table ({len(r):,} rows)")
    else:
        cols = ["drug", "concentration", "Cell_Name_Vevo", "Cell_ID_DepMap", "baseMean",
                "log2FoldChange", "stat", "padj", "n_cells_trt", "n_cells_ctrl"]
        keys = ["drug", "concentration", "Cell_Name_Vevo", "Cell_ID_DepMap"]
        acc = []
        parts = sorted(glob.glob(os.path.join(a.parts, "*.parquet")))
        for i, p in enumerate(parts, 1):
            d = pd.read_parquet(p, columns=cols)
            d["de"] = ((d.padj < 0.05) & (d.log2FoldChange.abs() > 0.5)).astype("int32")
            d["lfc2"] = np.where(d.baseMean >= 10, d.log2FoldChange.astype("float64") ** 2, np.nan)
            d["stat2"] = d.stat.astype("float64") ** 2
            d["expr"] = (d.baseMean >= 10).astype("int32")
            g = d.groupby(keys, observed=True).agg(n_de=("de", "sum"), sum_lfc2=("lfc2", "sum"),
                                                   n_expr=("expr", "sum"), sum_stat2=("stat2", "sum"),
                                                   n_genes=("stat2", "size"),
                                                   n_cells_trt=("n_cells_trt", "max"),
                                                   n_cells_ctrl=("n_cells_ctrl", "max")).reset_index()
            acc.append(g)
            if i % 100 == 0 or i == len(parts): print(f"  {i}/{len(parts)} files")
        r = pd.concat(acc)
        # a drug x dose x line can be split across files -> combine
        r = r.groupby(keys, observed=True).agg(n_de=("n_de", "sum"), sum_lfc2=("sum_lfc2", "sum"),
                                               n_expr=("n_expr", "sum"), sum_stat2=("sum_stat2", "sum"),
                                               n_genes=("n_genes", "sum"), n_cells_trt=("n_cells_trt", "max"),
                                               n_cells_ctrl=("n_cells_ctrl", "max")).reset_index()
        r["rms_lfc"] = np.sqrt(r.sum_lfc2 / r.n_expr.clip(lower=1))
        r["rms_stat"] = np.sqrt(r.sum_stat2 / r.n_genes.clip(lower=1))
        r["dose_rank"] = r.groupby("drug").concentration.rank(method="dense").astype(int)
        r["drug_key"] = r.drug.map(norm_drug)
        r = r.drop(columns=["sum_lfc2", "sum_stat2"])
        r.to_csv(cache, index=False)
    say(f"Response table: {r.drug.nunique()} drugs x {r.Cell_Name_Vevo.nunique()} lines x {r.dose_rank.max()} doses")

    # ---------- 2. labels ----------
    pr = pd.read_csv(a.prism, usecols=lambda c: c in {"depmap_id", "name", "auc", "screen_id"}, low_memory=False).dropna(subset=["auc"])
    if pr.screen_id.astype(str).str.contains("MTS010").any():
        pr = pr[pr.screen_id.astype(str).str.contains("MTS010")]
    pr["drug_key"] = pr["name"].map(norm_drug)
    pr = pr.groupby(["depmap_id", "drug_key"], as_index=False).auc.mean().rename(
        columns={"depmap_id": "Cell_ID_DepMap", "auc": "prism_auc"})
    gd = pd.read_csv(a.gdsc, usecols=["CELL_LINE_NAME", "DRUG_NAME", "AUC"])
    gd["cell_key"] = gd.CELL_LINE_NAME.map(norm_cell); gd["drug_key"] = gd.DRUG_NAME.map(norm_drug)
    gd = gd.groupby(["cell_key", "drug_key"], as_index=False).AUC.mean().rename(columns={"AUC": "gdsc_auc"})
    r["cell_key"] = r.Cell_Name_Vevo.map(norm_cell)
    cl = pd.read_csv(a.cells).drop_duplicates("cell_name")
    lung = set(cl.loc[cl.Organ.astype(str).str.contains("Lung", case=False), "cell_name"])
    r["lung"] = r.Cell_Name_Vevo.isin(lung)

    # ---------- 3. within-drug association + permutation test ----------
    rng = np.random.default_rng(0)
    rows = []
    for lab, tab, on in [("prism_auc", pr, ["Cell_ID_DepMap", "drug_key"]), ("gdsc_auc", gd, ["cell_key", "drug_key"])]:
        m = r.merge(tab, on=on, how="inner")
        m = m[m.dose_rank == m.groupby("drug").dose_rank.transform("max")]       # top dose
        m.to_csv(os.path.join(a.out, f"matched_{lab}.csv"), index=False)
        say(f"\n=== {lab}: {m.drug.nunique()} drugs, {m.Cell_Name_Vevo.nunique()} lines, {len(m):,} pairs "
            f"(lung: {int(m.lung.sum())} pairs) ===")
        for feat in ["n_de", "rms_lfc", "rms_stat"]:
            per = []
            for dname, g in m.groupby("drug"):
                if len(g) < a.min_lines or g[feat].nunique() < 3: continue
                rho = np.corrcoef(rank(g[feat]), rank(g[lab]))[0, 1]
                prho = partial_spearman(g[feat], g[lab], np.log1p(g.n_cells_trt))
                null = [np.corrcoef(rank(g[feat]), rank(rng.permutation(g[lab].values)))[0, 1] for _ in range(a.n_perm)]
                per.append(dict(label=lab, feature=feat, drug=dname, n_lines=len(g), rho=rho,
                                partial_rho=prho, null_mean=np.mean(null), null_sd=np.std(null)))
            if not per: continue
            P = pd.DataFrame(per); rows.append(P)
            # overall test: is the median within-drug rho below the permutation null?
            obs = P.rho.median()
            null_medians = []
            for _ in range(a.n_perm):
                null_medians.append(np.median([np.corrcoef(rank(g[feat]), rank(rng.permutation(g[lab].values)))[0, 1]
                                               for _, g in m.groupby("drug") if len(g) >= a.min_lines and g[feat].nunique() >= 3]))
            pval = (np.sum(np.array(null_medians) <= obs) + 1) / (a.n_perm + 1)
            say(f"  {feat:9s} median within-drug rho = {obs:+.3f} | partial (cell-count adjusted) = "
                f"{P.partial_rho.median():+.3f} | drugs with rho<0: {(P.rho < 0).mean():.0%} of {len(P)} | "
                f"one-sided permutation p = {pval:.4f}")
            top = P.sort_values("rho").head(6)
            say("     most negative drugs:", ", ".join(f"{d} ({v:+.2f})" for d, v in zip(top.drug, top.rho)))
    if rows:
        pd.concat(rows).to_csv(os.path.join(a.out, "per_drug_correlations.csv"), index=False)
    log.close(); print(f"\nSaved results to {a.out}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="../data/processed/pseudobulk_parts_all")
    ap.add_argument("--cells", default="../data/raw/tahoe_metadata/cell_line_metadata.csv")
    ap.add_argument("--gdsc", default="../data/raw/depmap/GDSC2_fitted_dose_response_24Jul22.csv")
    ap.add_argument("--prism", default="../data/raw/depmap/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv")
    ap.add_argument("--out", default="../results/response_vs_sensitivity")
    ap.add_argument("--min-lines", dest="min_lines", type=int, default=10)
    ap.add_argument("--n-perm", dest="n_perm", type=int, default=500)
    main(ap.parse_args())

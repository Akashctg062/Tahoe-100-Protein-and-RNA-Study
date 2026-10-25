"""
Step 9 - CONFIRMATORY analyses C2 and C3, exactly as frozen in
         03_Prespecified_Analysis_Plan_and_Manuscript_Outline (15 Oct 2026).

C2  Proliferation program (20 fixed genes): within-drug Spearman(score, AUC) > 0 ?   (PRISM & GDSC)
C3  Which data predict the cell lines a drug kills?  Leave-one-line-out ridge (alpha=1, PCA-5 in fold)
      1 Response (Tahoe top-dose log2FC)  2 Baseline RNA (CCLE)  3 Baseline protein MS (CCLE)
      4 Baseline protein RPPA (CCLE)      5 Response+RNA         6 Response+protein MS
    Permutation test per set (500), Wilcoxon pairwise comparisons, Benjamini-Hochberg FDR.

Inputs (all open; download commands in the chat / plan):
  ../data/processed/pseudobulk_parts_all/              (step 3 --all-drugs)
  ../data/raw/tahoe_metadata/gene_metadata.csv, cell_line_metadata.csv
  ../data/raw/depmap/  PRISM secondary, GDSC2, protein_quant_current_normalized.csv.gz,
                       CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz, CCLE_RPPA_20181003.csv
Run:
  python 09_confirmatory_C2_C3.py                 (primary analysis, ~30-60 min first time)
  python 09_confirmatory_C2_C3.py --robustness    (pre-specified robustness grid, after primary)
Outputs: ../results/confirmatory/
"""
import argparse, glob, os, re, sys, importlib, time, warnings
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
norm_drug = importlib.import_module("02_drug_overlap").norm

PROLIF = ["MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "CCNA2", "BUB1", "PLK1", "AURKA", "AURKB",
          "CENPF", "E2F1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "TYMS"]
R = "../data/raw/"

def norm_cell(x): return re.sub(r"[^A-Z0-9]", "", str(x).upper())
def rank(x): return pd.Series(x).rank().values
def spear(a, b):
    if np.std(a) == 0 or np.std(b) == 0: return np.nan
    return np.corrcoef(rank(a), rank(b))[0, 1]
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)

# ---------------------------------------------------------------- labels
def load_labels(lines):
    pr = pd.read_csv(R + "depmap/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv",
                     usecols=lambda c: c in {"depmap_id", "name", "auc", "screen_id"}, low_memory=False).dropna(subset=["auc"])
    if pr.screen_id.astype(str).str.contains("MTS010").any():
        pr = pr[pr.screen_id.astype(str).str.contains("MTS010")]
    pr["drug_key"] = pr["name"].map(norm_drug)
    d2c = dict(zip(lines.Cell_ID_DepMap, lines.cell_name))
    pr["cell_name"] = pr.depmap_id.map(d2c)
    pr = pr.dropna(subset=["cell_name"]).groupby(["drug_key", "cell_name"], as_index=False).auc.mean()
    gd = pd.read_csv(R + "depmap/GDSC2_fitted_dose_response_24Jul22.csv", usecols=["CELL_LINE_NAME", "DRUG_NAME", "AUC"])
    k2c = dict(zip(lines.cell_key, lines.cell_name))
    gd["cell_name"] = gd.CELL_LINE_NAME.map(norm_cell).map(k2c); gd["drug_key"] = gd.DRUG_NAME.map(norm_drug)
    gd = gd.dropna(subset=["cell_name"]).groupby(["drug_key", "cell_name"], as_index=False).AUC.mean().rename(columns={"AUC": "auc"})
    return {"PRISM": pr, "GDSC": gd}

# ---------------------------------------------------------------- Tahoe response cache
def build_response(a, drug_keys, out):
    cache = os.path.join(out, f"response_dose{a.dose}.parquet")
    if os.path.exists(cache): return pd.read_parquet(cache)
    parts = sorted(glob.glob(os.path.join(a.parts, "*.parquet")))
    # fixed, label-independent gene universe: 5,000 genes with highest mean baseMean in the first 20 files
    bm = pd.concat([pd.read_parquet(p, columns=["gene_name", "baseMean"]) for p in parts[:20]])
    universe = set(bm.groupby("gene_name").baseMean.mean().sort_values(ascending=False).index[:a.n_genes]) | set(PROLIF)
    print(f"  gene universe: {len(universe)} genes")
    acc = []
    for i, p in enumerate(parts, 1):
        d = pd.read_parquet(p, columns=["drug", "concentration", "Cell_Name_Vevo", "gene_name", "log2FoldChange"])
        d = d[d.drug.map(norm_drug).isin(drug_keys) & d.gene_name.isin(universe)]
        if len(d): acc.append(d)
        if i % 100 == 0 or i == len(parts): print(f"  {i}/{len(parts)} files scanned")
    d = pd.concat(acc)
    d["dose_rank"] = d.groupby("drug").concentration.rank(method="dense").astype(int)
    want = d.groupby("drug").dose_rank.transform("max") if a.dose == "top" else 2
    d = d[d.dose_rank == want]
    d["drug_key"] = d.drug.map(norm_drug)
    d = d.groupby(["drug_key", "Cell_Name_Vevo", "gene_name"], as_index=False).log2FoldChange.mean()
    d.to_parquet(cache, index=False)
    return d

# ---------------------------------------------------------------- baselines
def load_baselines(lines, gm):
    k2c = dict(zip(lines.cell_key, lines.cell_name))
    # RNA
    hdr = pd.read_csv(R + "depmap/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz", sep="\t", nrows=0).columns
    keep = [c for c in hdr if norm_cell(c.split("_")[0]) in k2c]
    rna = pd.read_csv(R + "depmap/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz", sep="\t", usecols=["gene_id"] + keep)
    ens2sym = dict(zip(gm.ensembl_id.astype(str).str.split(".").str[0], gm.gene_symbol))
    rna["sym"] = rna.gene_id.str.split(".").str[0].map(ens2sym)
    rna = rna.dropna(subset=["sym"]).groupby("sym")[keep].mean()
    rna.columns = [k2c[norm_cell(c.split("_")[0])] for c in rna.columns]
    rna = np.log2(rna.T.groupby(level=0).mean() + 1)
    # MS proteomics
    ph = pd.read_csv(R + "depmap/protein_quant_current_normalized.csv.gz", nrows=0).columns
    pcols = [c for c in ph if re.search(r"_TenPx\d+$", c) and norm_cell(c.split("_")[0]) in k2c]
    ms = pd.read_csv(R + "depmap/protein_quant_current_normalized.csv.gz", usecols=["Gene_Symbol"] + pcols)
    ms = ms.dropna(subset=["Gene_Symbol"]).groupby("Gene_Symbol")[pcols].mean().T
    ms.index = [k2c[norm_cell(c.split("_")[0])] for c in ms.index]
    ms = ms.groupby(level=0).mean()
    ms = ms.loc[:, ms.isna().mean() < 0.2]; ms = ms.fillna(ms.mean())
    # RPPA
    rp = pd.read_csv(R + "depmap/CCLE_RPPA_20181003.csv")
    rp = rp.rename(columns={rp.columns[0]: "ccle"}); rp["cell"] = rp.ccle.map(lambda c: k2c.get(norm_cell(str(c).split("_")[0])))
    rp = rp.dropna(subset=["cell"]).drop(columns="ccle").groupby("cell").mean()
    print(f"  baselines: RNA {rna.shape}, MS protein {ms.shape}, RPPA {rp.shape} (lines x features)")
    return rna, ms, rp

# ---------------------------------------------------------------- LOO ridge on in-fold PCA
def pca_fold_scores(X, k):
    """For each left-out i: fit standardise+PCA on the rest, return (train_scores, test_score) per fold."""
    n = X.shape[0]; folds = []
    for i in range(n):
        tr = np.arange(n) != i
        mu, sd = X[tr].mean(0), X[tr].std(0); sd[sd == 0] = 1
        Z = (X[tr] - mu) / sd; z = (X[i] - mu) / sd
        U, S, Vt = np.linalg.svd(Z, full_matrices=False)
        kk = min(k, Vt.shape[0]); V = Vt[:kk].T
        folds.append((Z @ V, z @ V))
    return folds

def loo_predict(folds_list, y, alpha):
    n = len(y); pred = np.empty(n)
    for i in range(n):
        Ztr = np.hstack([f[i][0] for f in folds_list]); ztest = np.hstack([f[i][1] for f in folds_list])
        m, s = Ztr.mean(0), Ztr.std(0); s[s == 0] = 1
        A = (Ztr - m) / s; b = (ztest - m) / s
        ytr = np.delete(y, i); ym = ytr.mean()
        w = np.linalg.solve(A.T @ A + alpha * np.eye(A.shape[1]), A.T @ (ytr - ym))
        pred[i] = ym + b @ w
    return pred

# ---------------------------------------------------------------- main
def main(a):
    out = a.out; os.makedirs(out, exist_ok=True)
    log = open(os.path.join(out, "summary_robustness.txt" if a.robustness else "summary.txt"), "w", encoding="utf-8")
    def say(*x):
        s = " ".join(str(i) for i in x); print(s); log.write(s + "\n"); log.flush()
    t0 = time.time(); rng = np.random.default_rng(2026)

    cl = pd.read_csv(R + "tahoe_metadata/cell_line_metadata.csv").drop_duplicates("cell_name")
    gm = pd.read_csv(R + "tahoe_metadata/gene_metadata.csv")
    labels = None
    # measured lines only (from the Tahoe response table of step 8)
    measured = set(pd.read_csv("../results/response_vs_sensitivity/line_response_strength.csv",
                               usecols=["Cell_Name_Vevo"]).Cell_Name_Vevo)
    lines = cl[cl.cell_name.isin(measured)].copy(); lines["cell_key"] = lines.cell_name.map(norm_cell)
    lines["lung"] = lines.Organ.astype(str).str.contains("Lung", case=False)
    labels = load_labels(lines)
    tahoe_keys = set(pd.read_csv("../results/response_vs_sensitivity/line_response_strength.csv", usecols=["drug_key"]).drug_key)
    keys = set()
    for L in labels.values():
        c = L[L.drug_key.isin(tahoe_keys)].groupby("drug_key").cell_name.nunique()
        keys |= set(c[c >= a.min_lines].index)
    say(f"Drugs with >= {a.min_lines} labelled Tahoe lines (PRISM or GDSC): {len(keys)}")
    resp = build_response(a, keys, out)
    say(f"Response table: {resp.drug_key.nunique()} drugs, {resp.Cell_Name_Vevo.nunique()} lines, {resp.gene_name.nunique()} genes")
    rna, ms, rp = load_baselines(lines, gm)

    # ======================= C2 =======================
    if not a.robustness:
        say("\n=== C2  Proliferation program (20 genes) vs AUC, within drug ===")
        pr_score = resp[resp.gene_name.isin(PROLIF)].groupby(["drug_key", "Cell_Name_Vevo"]).log2FoldChange.mean().rename("prolif").reset_index()
        c2rows = []
        for src, L in labels.items():
            m = pr_score.merge(L, left_on=["drug_key", "Cell_Name_Vevo"], right_on=["drug_key", "cell_name"])
            groups = [g for _, g in m.groupby("drug_key") if len(g) >= a.min_lines]
            rhos = np.array([spear(g.prolif.values, g.auc.values) for g in groups]); ok = ~np.isnan(rhos)
            obs = np.median(rhos[ok])
            null = np.array([np.nanmedian([spear(g.prolif.values, rng.permutation(g.auc.values)) for g in groups])
                             for _ in range(a.n_perm_c2)])
            p = (np.sum(null >= obs) + 1) / (a.n_perm_c2 + 1)
            say(f"  {src}: {ok.sum()} drugs | median rho = {obs:+.3f} | drugs rho>0: {(rhos[ok] > 0).mean():.0%} | one-sided p = {p:.4f} "
                f"| {'SUPPORTED' if p < 0.025 else 'not supported'} at Bonferroni 0.025")
            for g, r_ in zip(groups, rhos): c2rows.append(dict(source=src, drug_key=g.drug_key.iloc[0], n_lines=len(g), rho=r_))
        pd.DataFrame(c2rows).to_csv(os.path.join(out, "C2_per_drug.csv"), index=False)

    # ======================= C3 =======================
    grid = [(a.k, a.alpha)] if not a.robustness else [(3, 1.0), (10, 1.0), (5, 0.1), (5, 10.0)]
    for k, alpha in grid:
        say(f"\n=== C3  LOO ridge | PCA k={k} | alpha={alpha} | dose={a.dose} ===")
        allrows, preds = [], []
        for src, L in labels.items():
            drugs = sorted(set(L.drug_key) & set(resp.drug_key))
            per = {s: [] for s in ["Response", "RNA", "Protein_MS", "RPPA", "Response+RNA", "Response+Protein_MS"]}
            nulls = {s: [] for s in per}; used = []
            for dk in drugs:
                Ld = L[L.drug_key == dk].set_index("cell_name").auc
                R_ = resp[resp.drug_key == dk].pivot_table(index="Cell_Name_Vevo", columns="gene_name", values="log2FoldChange")
                common = sorted(set(Ld.index) & set(R_.index) & set(rna.index) & set(ms.index) & set(rp.index))
                if len(common) < a.min_lines: continue
                y = Ld.loc[common].values.astype(float)
                Xr = R_.loc[common].fillna(0).values
                g_shared = [g for g in R_.columns if g in rna.columns]
                feats = {"Response": Xr, "RNA": rna.loc[common, g_shared].values,
                         "Protein_MS": ms.loc[common].values, "RPPA": rp.loc[common].values}
                F = {s: pca_fold_scores(X, k) for s, X in feats.items()}
                combos = {"Response": ["Response"], "RNA": ["RNA"], "Protein_MS": ["Protein_MS"], "RPPA": ["RPPA"],
                          "Response+RNA": ["Response", "RNA"], "Response+Protein_MS": ["Response", "Protein_MS"]}
                used.append(dk)
                for s, comp in combos.items():
                    fl = [F[c] for c in comp]
                    pred = loo_predict(fl, y, alpha); per[s].append(spear(pred, y))
                    if not a.robustness:
                        nulls[s].append([spear(loo_predict(fl, yp, alpha), yp)
                                         for yp in (rng.permutation(y) for _ in range(a.n_perm_c3))])
                    preds += [dict(source=src, drug_key=dk, set=s, cell_name=c, observed=o, predicted=pp,
                                   lung=c in set(lines[lines.lung].cell_name)) for c, o, pp in zip(common, y, pred)]
            if not used: say(f"  {src}: no drug with >= {a.min_lines} common lines"); continue
            say(f"  {src}: {len(used)} drugs (common lines with all modalities, median n = "
                f"{int(np.median([ (L[L.drug_key==d].cell_name.isin(set(rna.index)&set(ms.index)&set(rp.index))).sum() for d in used]))})")
            for s in per:
                rh = np.array(per[s]); med = np.nanmedian(rh); row = dict(source=src, set=s, k=k, alpha=alpha, n_drugs=len(rh), median_rho=med)
                if not a.robustness:
                    null_med = np.nanmedian(np.array(nulls[s]), axis=0)          # median across drugs, per permutation
                    row["perm_p"] = (np.sum(null_med >= med) + 1) / (len(null_med) + 1)
                allrows.append(row)
            if not a.robustness:
                for s1, s2 in [("Response", "RNA"), ("Response", "Protein_MS"), ("Protein_MS", "RNA"),
                               ("Response+RNA", "RNA"), ("Response+Protein_MS", "Protein_MS"), ("RPPA", "Protein_MS")]:
                    d1, d2 = np.array(per[s1]), np.array(per[s2]); ok = ~np.isnan(d1) & ~np.isnan(d2)
                    p = wilcoxon(d1[ok], d2[ok]).pvalue if ok.sum() >= 5 and np.any(d1[ok] != d2[ok]) else np.nan
                    allrows.append(dict(source=src, set=f"{s1} vs {s2}", k=k, alpha=alpha, n_drugs=int(ok.sum()),
                                        median_rho=np.median(d1[ok] - d2[ok]), wilcoxon_p=p))
        T = pd.DataFrame(allrows)
        if not a.robustness:
            T["p"] = T["perm_p"].fillna(T["wilcoxon_p"]) if "wilcoxon_p" in T else T["perm_p"]
            T["fdr_bh"] = bh(T.p.fillna(1).values)
            T["verdict"] = np.where(T.set.str.contains(" vs "),
                                    np.where(T.fdr_bh < 0.05, "different", "no difference"),
                                    np.where((T.fdr_bh < 0.05) & (T.median_rho > 0.10), "PREDICTS", "does not predict"))
        say(T.round(4).to_string(index=False))
        T.to_csv(os.path.join(out, f"C3_summary_k{k}_a{alpha}_{a.dose}.csv"), index=False)
        pd.DataFrame(preds).to_csv(os.path.join(out, f"C3_predictions_k{k}_a{alpha}_{a.dose}.csv"), index=False)
    say(f"\nDone in {(time.time() - t0) / 60:.1f} min. Results in {out}")
    log.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="../data/processed/pseudobulk_parts_all")
    ap.add_argument("--out", default="../results/confirmatory")
    ap.add_argument("--min-lines", dest="min_lines", type=int, default=20)
    ap.add_argument("--n-genes", dest="n_genes", type=int, default=5000)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--dose", default="top", choices=["top", "mid"])
    ap.add_argument("--n-perm-c2", dest="n_perm_c2", type=int, default=1000)
    ap.add_argument("--n-perm-c3", dest="n_perm_c3", type=int, default=500)
    ap.add_argument("--robustness", action="store_true")
    main(ap.parse_args())

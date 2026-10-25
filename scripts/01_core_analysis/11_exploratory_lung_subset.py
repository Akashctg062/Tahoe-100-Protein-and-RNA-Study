"""
Step 11 - EXPLORATORY (post hoc, not in the frozen plan; added 16 Oct 2026). Label as exploratory in the paper.

Question: within the C3 leave-one-line-out models (trained on all lines, as frozen), how well are LUNG
lines ranked compared with non-lung lines?  Each drug's LOO predictions are evaluated separately on its
lung lines (>= 6 required) and its non-lung lines, and each subset rho is null-calibrated with the same
within-drug label permutations (200) passed through the identical LOO procedure.

Run:  python 11_exploratory_lung_subset.py        Output: ../results/exploratory/lung_subset/
"""
import os, sys, time, importlib, warnings
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_confirmatory_C2_C3")
spear, bh, R = m09.spear, m09.bh, m09.R

def main(n_perm=200, k=5, alpha=1.0, min_lines=20, min_lung=6, out="../results/exploratory/lung_subset"):
    os.makedirs(out, exist_ok=True); t0 = time.time(); rng = np.random.default_rng(2026)
    cl = pd.read_csv(R + "tahoe_metadata/cell_line_metadata.csv").drop_duplicates("cell_name")
    gm = pd.read_csv(R + "tahoe_metadata/gene_metadata.csv")
    lrs = pd.read_csv("../results/response_vs_sensitivity/line_response_strength.csv", usecols=["Cell_Name_Vevo", "drug_key"])
    lines = cl[cl.cell_name.isin(set(lrs.Cell_Name_Vevo))].copy(); lines["cell_key"] = lines.cell_name.map(m09.norm_cell)
    lines["lung"] = lines.Organ.astype(str).str.contains("Lung", case=False)
    is_lung = lines.set_index("cell_name").lung
    labels = m09.load_labels(lines)
    resp = pd.read_parquet("../results/confirmatory/response_dosetop.parquet")
    rna, ms, rp = m09.load_baselines(lines, gm)
    sets = ["Response", "RNA", "Protein_MS", "RPPA"]
    rows = []
    for src, L in labels.items():
        for dk in sorted(set(L.drug_key) & set(resp.drug_key)):
            Ld = L[L.drug_key == dk].set_index("cell_name").auc
            R_ = resp[resp.drug_key == dk].pivot_table(index="Cell_Name_Vevo", columns="gene_name", values="log2FoldChange")
            common = sorted(set(Ld.index) & set(R_.index) & set(rna.index) & set(ms.index) & set(rp.index))
            if len(common) < min_lines: continue
            lung = is_lung.reindex(common).fillna(False).values.astype(bool)
            if lung.sum() < min_lung or (~lung).sum() < min_lung: continue
            y = Ld.loc[common].values.astype(float)
            g_shared = [g for g in R_.columns if g in rna.columns]
            feats = {"Response": R_.loc[common].fillna(0).values, "RNA": rna.loc[common, g_shared].values,
                     "Protein_MS": ms.loc[common].values, "RPPA": rp.loc[common].values}
            F = {s: m09.pca_fold_scores(X, k) for s, X in feats.items()}
            perms = [rng.permutation(y) for _ in range(n_perm)]
            for s in sets:
                pr = m09.loo_predict([F[s]], y, alpha)
                nl = np.array([[spear(p_[lung], yp[lung]), spear(p_[~lung], yp[~lung])]
                               for yp in perms for p_ in [m09.loo_predict([F[s]], yp, alpha)]])
                for j, (sub, msk) in enumerate([("lung", lung), ("non-lung", ~lung)]):
                    o = spear(pr[msk], y[msk])
                    rows.append(dict(source=src, drug_key=dk, set=s, subset=sub, n_lines=int(msk.sum()), rho=o,
                                     null_mean=np.nanmean(nl[:, j]), calibrated_rho=o - np.nanmean(nl[:, j])))
            print(f"  {src} {dk}: lung {lung.sum()} / other {(~lung).sum()} ({(time.time()-t0)/60:.1f} min)", flush=True)
    D = pd.DataFrame(rows); D.to_csv(os.path.join(out, "per_drug_subset_calibrated.csv"), index=False)
    S = D.groupby(["source", "set", "subset"]).agg(n_drugs=("drug_key", "nunique"), median_lines=("n_lines", "median"),
                                                   median_raw=("rho", "median"), median_null=("null_mean", "median"),
                                                   median_calibrated=("calibrated_rho", "median")).reset_index()
    S["p_gt0"] = [wilcoxon(D[(D.source == a) & (D.set == b) & (D.subset == c)].calibrated_rho.dropna(), alternative="greater").pvalue
                  for a, b, c in zip(S.source, S.set, S.subset)]
    S["fdr_bh"] = bh(S.p_gt0.values)
    S.to_csv(os.path.join(out, "summary.csv"), index=False)
    txt = ("EXPLORATORY lung-subset evaluation of C3 (post hoc)\n\n" + S.round(4).to_string(index=False)
           + f"\n\nDone in {(time.time()-t0)/60:.1f} min")
    open(os.path.join(out, "summary.txt"), "w", encoding="utf-8").write(txt); print("\n" + txt)

if __name__ == "__main__":
    main()

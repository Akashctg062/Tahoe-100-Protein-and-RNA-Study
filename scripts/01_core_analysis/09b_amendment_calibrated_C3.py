"""
Step 9b - AMENDMENT 1 to the frozen plan (added 15 Oct 2026, AFTER seeing C3 results; label as such).

Reason: leave-one-line-out (LOO) prediction with a training-mean intercept is negatively biased under
the null (the left-out value shifts the training mean in the opposite direction), and the size of the
bias differs between feature sets. The frozen permutation p-values remain valid (each permutation goes
through the same LOO procedure), but (i) the "median rho > 0.10" effect-size rule and (ii) the pairwise
Wilcoxon comparisons on RAW rho are not interpretable.

This script re-runs C3 exactly as frozen (same data, PCA k=5, ridge alpha=1, LOO, 500 permutations,
same seed) and additionally reports, per drug and feature set:
    calibrated_rho = observed rho - mean(rho under 500 within-drug permutations)
Sets are then compared on calibrated rho (two-sided Wilcoxon signed-rank across drugs, BH-FDR).
The frozen results are reported unchanged alongside.

Run:  python 09b_amendment_calibrated_C3.py           Output: ../results/confirmatory/amendment1/
"""
import os, sys, time, importlib, warnings
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m09 = importlib.import_module("09_confirmatory_C2_C3")
spear, bh = m09.spear, m09.bh
R = m09.R

def main(n_perm=500, k=5, alpha=1.0, min_lines=20, out="../results/confirmatory/amendment1"):
    os.makedirs(out, exist_ok=True); t0 = time.time(); rng = np.random.default_rng(2026)
    cl = pd.read_csv(R + "tahoe_metadata/cell_line_metadata.csv").drop_duplicates("cell_name")
    gm = pd.read_csv(R + "tahoe_metadata/gene_metadata.csv")
    lrs = pd.read_csv("../results/response_vs_sensitivity/line_response_strength.csv", usecols=["Cell_Name_Vevo", "drug_key"])
    lines = cl[cl.cell_name.isin(set(lrs.Cell_Name_Vevo))].copy(); lines["cell_key"] = lines.cell_name.map(m09.norm_cell)
    lines["lung"] = lines.Organ.astype(str).str.contains("Lung", case=False)
    labels = m09.load_labels(lines)
    resp = pd.read_parquet("../results/confirmatory/response_dosetop.parquet")
    rna, ms, rp = m09.load_baselines(lines, gm)
    combos = {"Response": ["Response"], "RNA": ["RNA"], "Protein_MS": ["Protein_MS"], "RPPA": ["RPPA"],
              "Response+RNA": ["Response", "RNA"], "Response+Protein_MS": ["Response", "Protein_MS"]}
    rows = []
    for src, L in labels.items():
        for dk in sorted(set(L.drug_key) & set(resp.drug_key)):
            Ld = L[L.drug_key == dk].set_index("cell_name").auc
            R_ = resp[resp.drug_key == dk].pivot_table(index="Cell_Name_Vevo", columns="gene_name", values="log2FoldChange")
            common = sorted(set(Ld.index) & set(R_.index) & set(rna.index) & set(ms.index) & set(rp.index))
            if len(common) < min_lines: continue
            y = Ld.loc[common].values.astype(float)
            g_shared = [g for g in R_.columns if g in rna.columns]
            feats = {"Response": R_.loc[common].fillna(0).values, "RNA": rna.loc[common, g_shared].values,
                     "Protein_MS": ms.loc[common].values, "RPPA": rp.loc[common].values}
            F = {s: m09.pca_fold_scores(X, k) for s, X in feats.items()}
            perms = [rng.permutation(y) for _ in range(n_perm)]           # same permutations for every set
            for s, comp in combos.items():
                fl = [F[c] for c in comp]
                obs = spear(m09.loo_predict(fl, y, alpha), y)
                null = np.array([spear(m09.loo_predict(fl, yp, alpha), yp) for yp in perms])
                rows.append(dict(source=src, drug_key=dk, set=s, n_lines=len(common), rho=obs,
                                 null_mean=np.nanmean(null), null_sd=np.nanstd(null),
                                 calibrated_rho=obs - np.nanmean(null),
                                 z=(obs - np.nanmean(null)) / (np.nanstd(null) + 1e-12),
                                 p_drug=(np.sum(null >= obs) + 1) / (n_perm + 1)))
            print(f"  {src} {dk}: done ({(time.time()-t0)/60:.1f} min)", flush=True)
    D = pd.DataFrame(rows); D.to_csv(os.path.join(out, "per_drug_calibrated.csv"), index=False)

    # summaries
    S = D.groupby(["source", "set"]).agg(n_drugs=("drug_key", "nunique"), median_raw_rho=("rho", "median"),
                                         median_null_rho=("null_mean", "median"),
                                         median_calibrated_rho=("calibrated_rho", "median"),
                                         frac_drugs_calibrated_gt0=("calibrated_rho", lambda x: (x > 0).mean()),
                                         frac_drugs_p_lt_0_05=("p_drug", lambda x: (x < 0.05).mean())).reset_index()
    # one-sample test of calibrated rho > 0 across drugs (Wilcoxon, one-sided)
    S["p_calibrated_gt0"] = [wilcoxon(D[(D.source == a) & (D.set == b)].calibrated_rho.dropna(), alternative="greater").pvalue
                             for a, b in zip(S.source, S.set)]
    pairs = [("Response", "RNA"), ("Response", "Protein_MS"), ("Protein_MS", "RNA"), ("RPPA", "RNA"),
             ("Response+RNA", "RNA"), ("Response+Protein_MS", "Protein_MS"), ("RPPA", "Protein_MS")]
    C = []
    for src in D.source.unique():
        W = D[D.source == src].pivot(index="drug_key", columns="set", values="calibrated_rho")
        for a, b in pairs:
            d = (W[a] - W[b]).dropna()
            C.append(dict(source=src, comparison=f"{a} vs {b}", n_drugs=len(d), median_difference=d.median(),
                          p=wilcoxon(d).pvalue if (d != 0).sum() >= 5 else np.nan))
    C = pd.DataFrame(C)
    allp = np.r_[S.p_calibrated_gt0.values, C.p.fillna(1).values]; q = bh(allp)
    S["fdr_bh"] = q[:len(S)]; C["fdr_bh"] = q[len(S):]
    S.to_csv(os.path.join(out, "summary_sets.csv"), index=False); C.to_csv(os.path.join(out, "summary_pairs.csv"), index=False)
    pd.set_option("display.width", 220)
    txt = ("AMENDMENT 1 - null-calibrated C3 (post hoc; frozen results unchanged)\n\n"
           + S.round(4).to_string(index=False) + "\n\n" + C.round(4).to_string(index=False)
           + f"\n\nDone in {(time.time()-t0)/60:.1f} min")
    open(os.path.join(out, "summary.txt"), "w", encoding="utf-8").write(txt); print("\n" + txt)

if __name__ == "__main__":
    main()

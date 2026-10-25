"""
v2 / T1-D (pre-specified in the frozen plan of 15 Oct 2026; Amendment 3 Section 4) - middle-dose robustness of C3.
Re-runs Amendment 1 (calibrated C3: PCA-5 in fold, ridge alpha = 1, leave-one-line-out, 500 permutations, seed 2026)
with the Tahoe MIDDLE dose instead of the top dose.  Baseline sets are identical, so only Response sets can change.

Run on your computer (PowerShell; about 20-30 minutes):
    conda activate oncocardio
    cd "E:\\Protein and RNA\\v2_extended\\scripts"
    python T1D_middle_dose.py
It reads the v1 scripts and data (read-only) and writes ONLY to  E:\\Protein and RNA\\v2_extended\\results\\T1D\\
"""
import os, sys, time, importlib, warnings, types
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
V1 = os.path.abspath(os.path.join(HERE, "..", "..", "scripts"))           # E:\Protein and RNA\scripts (read-only)
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "T1D"))
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, V1); os.chdir(V1)                                        # v1 relative paths (../data, ../results)
m09 = importlib.import_module("09_confirmatory_C2_C3")
spear, bh, R = m09.spear, m09.bh, m09.R

def main(n_perm=500, k=5, alpha=1.0, min_lines=20):
    t0 = time.time(); rng = np.random.default_rng(2026)
    cl = pd.read_csv(R + "tahoe_metadata/cell_line_metadata.csv").drop_duplicates("cell_name")
    gm = pd.read_csv(R + "tahoe_metadata/gene_metadata.csv")
    lrs = pd.read_csv("../results/response_vs_sensitivity/line_response_strength.csv", usecols=["Cell_Name_Vevo", "drug_key"])
    lines = cl[cl.cell_name.isin(set(lrs.Cell_Name_Vevo))].copy(); lines["cell_key"] = lines.cell_name.map(m09.norm_cell)
    labels = m09.load_labels(lines)
    keys = set(labels["PRISM"].drug_key) | set(labels["GDSC"].drug_key)
    a = types.SimpleNamespace(parts="../data/processed/pseudobulk_parts_all", dose="mid", n_genes=5000)
    print("Building middle-dose response table (first run scans all Tahoe files, ~10-20 min) ...", flush=True)
    resp = m09.build_response(a, keys, OUT)                                 # cached at v2_extended/results/T1D/response_dosemid.parquet
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
            perms = [rng.permutation(y) for _ in range(n_perm)]
            for s, comp in combos.items():
                fl = [F[c] for c in comp]
                obs = spear(m09.loo_predict(fl, y, alpha), y)
                null = np.array([spear(m09.loo_predict(fl, yp, alpha), yp) for yp in perms])
                rows.append(dict(source=src, drug_key=dk, set=s, n_lines=len(common), rho=obs, null_mean=np.nanmean(null),
                                 calibrated_rho=obs - np.nanmean(null), p_drug=(np.sum(null >= obs) + 1) / (n_perm + 1)))
            print(f"  {src} {dk}: done ({(time.time()-t0)/60:.1f} min)", flush=True)
    D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, "per_drug_calibrated_middose.csv"), index=False)
    S = D.groupby(["source", "set"]).agg(n_drugs=("drug_key", "nunique"), median_raw_rho=("rho", "median"),
                                         median_null_rho=("null_mean", "median"), median_calibrated_rho=("calibrated_rho", "median")).reset_index()
    S["p_calibrated_gt0"] = [wilcoxon(D[(D.source == x) & (D.set == z)].calibrated_rho.dropna(), alternative="greater").pvalue
                             for x, z in zip(S.source, S.set)]
    pairs = [("Response", "RNA"), ("Response", "Protein_MS"), ("Protein_MS", "RNA"), ("RPPA", "RNA"),
             ("Response+RNA", "RNA"), ("Response+Protein_MS", "Protein_MS"), ("RPPA", "Protein_MS")]
    C = []
    for src in D.source.unique():
        W = D[D.source == src].pivot(index="drug_key", columns="set", values="calibrated_rho")
        for x, z in pairs:
            d = (W[x] - W[z]).dropna()
            C.append(dict(source=src, comparison=f"{x} vs {z}", n_drugs=len(d), median_difference=d.median(),
                          p=wilcoxon(d).pvalue if (d != 0).sum() >= 5 else np.nan))
    C = pd.DataFrame(C); q = bh(np.r_[S.p_calibrated_gt0.values, C.p.fillna(1).values])
    S["fdr_bh"] = q[:len(S)]; C["fdr_bh"] = q[len(S):]
    S.to_csv(os.path.join(OUT, "summary_sets_middose.csv"), index=False); C.to_csv(os.path.join(OUT, "summary_pairs_middose.csv"), index=False)
    txt = ("T1-D middle-dose robustness of C3 (Amendment 1 procedure)\n\n" + S.round(4).to_string(index=False)
           + "\n\n" + C.round(4).to_string(index=False) + f"\n\nDone in {(time.time()-t0)/60:.1f} min")
    open(os.path.join(OUT, "summary_middose.txt"), "w", encoding="utf-8").write(txt); print("\n" + txt)

if __name__ == "__main__":
    main()

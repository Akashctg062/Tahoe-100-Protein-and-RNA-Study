"""
v2 / T2-E and T2-F (Amendment 3, frozen 18 Oct 2026).
E: CTRP v2 replication of the Amendment-1 calibrated C3 (PCA-5 in fold, ridge alpha 1, leave-one-line-out, 500 permutations).
   E1 RPPA calibrated rho > 0; E2 RPPA - RNA > 0 (one-sided Wilcoxon across drugs); replicated if both p < 0.05.
F: ProCan-DepMapSanger proteome replacing CCLE MS proteomics, in PRISM, GDSC and CTRP (common lines need ProCan, RNA, RPPA).
   F1 ProCan calibrated rho > 0; F2 ProCan - RNA > 0 (one-sided Wilcoxon); supported if F2 p < 0.05 in >= 2 of 3 screens.
Inputs: response_top_80drugs.parquet (built by T2_build_responses.py), v1 baselines (cache), CTRP and ProCan derived tables.
Run:  python T2EF_ctrp_procan.py --part E   |  --part F
"""
import os, sys, re, pickle, importlib, argparse, time, warnings
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/claude/mirror/scripts")
U = "/mnt/user-data/uploads/Protein and RNA/"; D = U + "data/raw/v2/derived/"
OUT = "/home/claude/v2/results/T2EF"; os.makedirs(OUT, exist_ok=True)
PAIRS = [("Response", "RNA"), ("Response", "Protein"), ("Protein", "RNA"), ("RPPA", "RNA"),
         ("Response+RNA", "RNA"), ("Response+Protein", "Protein"), ("RPPA", "Protein")]

def setup():
    cwd = os.getcwd(); os.chdir("/home/claude/mirror/scripts")
    m09 = importlib.import_module("09_confirmatory_C2_C3"); norm = importlib.import_module("02_drug_overlap").norm
    B = pickle.load(open("/home/claude/v2/results/cache/baselines.pkl", "rb"))
    labels = m09.load_labels(B["lines"]); os.chdir(cwd)
    return m09, norm, B, labels

def ctrp_labels(norm, lines):
    c = pd.read_csv(D + "ctrp_auc.csv", low_memory=False)
    dm = pd.read_csv(U + "data/raw/tahoe_metadata/drug_metadata.csv").drop_duplicates("drug"); dm["key"] = dm.drug.map(norm)
    pc = pd.to_numeric(dm.pubchem_cid, errors="coerce"); pub = pd.Series(dm.key.values, index=pc).dropna(); pub = pub[~pub.index.duplicated()]
    c["drug_key"] = pd.to_numeric(c.pubchem_id, errors="coerce").map(pub).fillna(c.drug_name.map(norm))
    k2c = dict(zip(lines.cell_key, lines.cell_name))
    c["cell_name"] = c.cell_line_name.map(lambda x: k2c.get(re.sub(r"[^A-Z0-9]", "", str(x).upper())))
    c = c.dropna(subset=["cell_name", "AUC_curvecurator"])
    return c.groupby(["drug_key", "cell_name"], as_index=False).AUC_curvecurator.mean().rename(columns={"AUC_curvecurator": "auc"})

def procan(lines):
    p = pd.read_csv(D + "procan_log2_mean.csv")
    k2c = dict(zip(lines.cell_key, lines.cell_name))
    p["cell"] = p.cell_line_name.map(lambda x: k2c.get(re.sub(r"[^A-Z0-9]", "", str(x).upper())))
    p = p.dropna(subset=["cell"]).drop(columns=["cellosaurus_id", "cell_line_name"]).groupby("cell").mean()
    p = p.loc[:, p.isna().mean() < 0.2]; return p.fillna(p.mean())

def analyse(m09, labels_by_src, resp, rna, prot, rp, tag, n_perm=500, k=5, alpha=1.0, min_lines=20):
    rng = np.random.default_rng(2026); rows = []; t0 = time.time()
    combos = {"Response": ["Response"], "RNA": ["RNA"], "Protein": ["Protein"], "RPPA": ["RPPA"],
              "Response+RNA": ["Response", "RNA"], "Response+Protein": ["Response", "Protein"]}
    have = set(resp.drug_key.astype(str).unique())
    for src, L in labels_by_src.items():
        for dk in sorted(set(L.drug_key) & have):
            Ld = L[L.drug_key == dk].set_index("cell_name").auc
            sub = resp[resp.drug_key == dk]
            R_ = sub.pivot_table(index=sub.Cell_Name_Vevo.astype(str), columns=sub.gene_name.astype(str), values="log2FoldChange", observed=True)
            common = sorted(set(Ld.index) & set(R_.index) & set(rna.index) & set(prot.index) & set(rp.index))
            if len(common) < min_lines: continue
            y = Ld.loc[common].values.astype(float); g_shared = [g for g in R_.columns if g in rna.columns]
            feats = {"Response": R_.loc[common].fillna(0).values, "RNA": rna.loc[common, g_shared].values,
                     "Protein": prot.loc[common].values, "RPPA": rp.loc[common].values}
            F = {s: m09.pca_fold_scores(X, k) for s, X in feats.items()}
            perms = [rng.permutation(y) for _ in range(n_perm)]
            for s, comp in combos.items():
                fl = [F[c] for c in comp]; obs = m09.spear(m09.loo_predict(fl, y, alpha), y)
                null = np.array([m09.spear(m09.loo_predict(fl, yp, alpha), yp) for yp in perms])
                rows.append(dict(source=src, drug_key=dk, set=s, n_lines=len(common), rho=obs, null_mean=np.nanmean(null),
                                 calibrated_rho=obs - np.nanmean(null), p_drug=(np.sum(null >= obs) + 1) / (n_perm + 1)))
            print(f"  [{tag}] {src} {dk}: {len(common)} lines ({(time.time()-t0)/60:.1f} min)", flush=True)
    Dd = pd.DataFrame(rows); Dd.to_csv(os.path.join(OUT, f"per_drug_{tag}.csv"), index=False)
    S = Dd.groupby(["source", "set"]).agg(n_drugs=("drug_key", "nunique"), median_lines=("n_lines", "median"),
                                          median_raw=("rho", "median"), median_null=("null_mean", "median"),
                                          median_calibrated=("calibrated_rho", "median")).reset_index()
    S["p_gt0"] = [wilcoxon(Dd[(Dd.source == a) & (Dd.set == b)].calibrated_rho.dropna(), alternative="greater").pvalue for a, b in zip(S.source, S.set)]
    C = []
    for src in Dd.source.unique():
        W = Dd[Dd.source == src].pivot(index="drug_key", columns="set", values="calibrated_rho")
        for a, b in PAIRS:
            d = (W[a] - W[b]).dropna()
            C.append(dict(source=src, comparison=f"{a} vs {b}", n_drugs=len(d), median_difference=d.median(),
                          p_two_sided=wilcoxon(d).pvalue, p_greater=wilcoxon(d, alternative="greater").pvalue))
    C = pd.DataFrame(C)
    q = m09.bh(np.r_[S.p_gt0.values, C.p_two_sided.values]); S["fdr_bh"] = q[:len(S)]; C["fdr_bh"] = q[len(S):]
    S.to_csv(os.path.join(OUT, f"summary_sets_{tag}.csv"), index=False); C.to_csv(os.path.join(OUT, f"summary_pairs_{tag}.csv"), index=False)
    return Dd, S, C

def main(a):
    m09, norm, B, labels = setup()
    resp = pd.read_parquet(a.response)
    for c in ["drug_key", "Cell_Name_Vevo", "gene_name"]: resp[c] = resp[c].astype("category")
    resp["log2FoldChange"] = resp.log2FoldChange.astype("float32")
    rna, ms, rp, lines = B["rna"], B["ms"], B["rppa"], B["lines"]
    if a.part == "D":      # T1-D middle-dose robustness: Amendment-1 procedure, PRISM + GDSC, CCLE MS proteome
        frozen = set(pd.read_parquet("/home/claude/mirror/results/confirmatory/response_dosetop.parquet", columns=["drug_key"]).drug_key)
        resp = resp[resp.drug_key.astype(str).isin(frozen)]                 # same drug set as Amendment 1
        Dd, S, C = analyse(m09, {"PRISM": labels["PRISM"], "GDSC": labels["GDSC"]}, resp, rna, ms, rp, "T1D_middose")
        txt = "T1-D middle-dose robustness of calibrated C3 (Protein = CCLE MS)\n\n" + S.round(4).to_string(index=False) + "\n\n" + C.round(4).to_string(index=False)
    elif a.part == "E":
        ctrp = ctrp_labels(norm, lines)
        Dd, S, C = analyse(m09, {"CTRP": ctrp}, resp, rna, ms, rp, "E_ctrp")
        e1 = S[S.set == "RPPA"].p_gt0.iloc[0]; e2 = C[C.comparison == "RPPA vs RNA"].p_greater.iloc[0]
        gd = set(labels["GDSC"].drug_key)
        sub = Dd[~Dd.drug_key.isin(gd)]
        W = sub.pivot(index="drug_key", columns="set", values="calibrated_rho") if len(sub) else None
        sens = (f"sensitivity (drugs not in GDSC, n = {sub.drug_key.nunique()}): RPPA median {W.RPPA.median():+.3f}, RPPA-RNA median {(W.RPPA - W.RNA).median():+.3f}"
                if W is not None and len(W) >= 5 else "sensitivity (drugs not in GDSC): too few drugs")
        txt = (f"T2-E CTRP replication: {S.n_drugs.max()} drugs, median {int(S.median_lines.median())} lines\n"
               f"E1 RPPA calibrated rho > 0: median {S[S.set=='RPPA'].median_calibrated.iloc[0]:+.3f}, p = {e1:.4f}\n"
               f"E2 RPPA - RNA > 0: median {C[C.comparison=='RPPA vs RNA'].median_difference.iloc[0]:+.3f}, p = {e2:.4f}\n"
               f"DECISION: GDSC RPPA finding " + ("REPLICATED in CTRP" if (e1 < 0.05 and e2 < 0.05) else "NOT replicated in CTRP") + "\n" + sens
               + "\n\n" + S.round(4).to_string(index=False) + "\n\n" + C.round(4).to_string(index=False))
    else:
        pr = procan(lines); ctrp = ctrp_labels(norm, lines)
        Dd, S, C = analyse(m09, {"PRISM": labels["PRISM"], "GDSC": labels["GDSC"], "CTRP": ctrp}, resp, rna, pr, rp, "F_procan")
        lines_ = [f"T2-F ProCan proteome ({pr.shape[0]} Tahoe lines x {pr.shape[1]} proteins)"]; n_ok = 0
        for src in ["PRISM", "GDSC", "CTRP"]:
            s1 = S[(S.source == src) & (S.set == "Protein")]; c2 = C[(C.source == src) & (C.comparison == "Protein vs RNA")]
            if not len(s1): continue
            f2 = c2.p_greater.iloc[0]; n_ok += f2 < 0.05
            lines_.append(f"  {src}: {int(s1.n_drugs.iloc[0])} drugs | F1 ProCan median {s1.median_calibrated.iloc[0]:+.3f} (p {s1.p_gt0.iloc[0]:.4f}) "
                          f"| F2 ProCan - RNA median {c2.median_difference.iloc[0]:+.3f} (p {f2:.4f})")
        lines_.append(f"DECISION: 'protein > RNA' {'SUPPORTED' if n_ok >= 2 else 'NOT supported'} by ProCan (F2 p < 0.05 in {n_ok} of 3 screens)")
        txt = "\n".join(lines_) + "\n\n" + S.round(4).to_string(index=False) + "\n\n" + C.round(4).to_string(index=False)
    open(os.path.join(OUT, f"summary_{a.part}.txt"), "w").write(txt); print(txt)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["D", "E", "F"], required=True)
    ap.add_argument("--response", default="/mnt/user-data/uploads/Protein and RNA/v2_extended/results/responses/response_top_80drugs.parquet")
    main(ap.parse_args())

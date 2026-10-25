"""
v2 / T2-H2 (Amendment 3, frozen 18 Oct 2026) - does Tahoe add to chemistry for hERG inhibition (ChEMBL CHEMBL240)?
Labels: human hERG IC50/Ki, standard_units nM, relation "=", no validity comment; pIC50 = 9 - log10(nM); median per
compound; matched to Tahoe drugs by InChIKey connectivity block (largest fragment) computed from SMILES.
Tasks: regression (pIC50; Spearman rho of out-of-fold predictions; PRIMARY) and classification (pIC50 >= 5; AUROC).
Models as in C1: chemistry (standardised; ridge alpha 100 / logistic C 0.01) vs chemistry + Tahoe PCA-10 (in fold);
stratified/plain 5-fold x 10 repeats (seed 2026), OOF averaged.  Tahoe adds if lower 95% bootstrap CI of delta-rho > 0.
"""
import os, sys, importlib, numpy as np, pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger; RDLogger.DisableLog("rdApp.*")
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import RepeatedKFold, RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
sys.path.insert(0, "/home/claude/mirror/scripts"); norm = importlib.import_module("02_drug_overlap").norm
sys.path.insert(0, "/home/claude/v2/scripts")
U = "/mnt/user-data/uploads/Protein and RNA/"; OUT = "/home/claude/v2/results/T2H2"; os.makedirs(OUT, exist_ok=True)
REM = U + "v2_extended/results/responses/tahoe_signatures_remaining.npz"
lf = rdMolStandardize.LargestFragmentChooser()

def ik(s):
    m = Chem.MolFromSmiles(str(s))
    if m is None: return None
    try: return Chem.MolToInchiKey(lf.choose(m)).split("-")[0]
    except Exception: return None

def signatures(wanted):
    """Load only the wanted drugs; return {key: (3 x genes array, genes)} using the first file that has the drug."""
    files = [U + "results/confirmatory/C1/tahoe_features_C1.npz", "/home/claude/v2/results/cache/tahoe_alldrug_signatures.npz"]
    if os.path.exists(REM): files.append(REM)
    sig = {}
    for f in files:
        z = np.load(f, allow_pickle=True); keys = list(z["keys"]); genes = list(z["genes"]); X = None
        for k in wanted:
            if k in sig or k not in keys: continue
            if X is None: X = z["X"]
            sig[k] = (np.asarray(X[keys.index(k)]), genes)
        del X
    return sig, files

def main():
    import importlib.util
    spec = importlib.util.spec_from_file_location("t1a1", "/home/claude/v2/scripts/T1A1_best_effort_cardiotox.py")
    t1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t1)
    dm = pd.read_csv(U + "data/raw/tahoe_metadata/drug_metadata.csv").drop_duplicates("drug"); dm["key"] = dm.drug.map(norm); dm["ik"] = dm.canonical_smiles.map(ik)
    h = pd.read_csv(U + "data/raw/v2/chembl_herg/chembl_CHEMBL240_IC50_Ki.csv", low_memory=False)
    f = h[(h.standard_units == "nM") & (h.standard_relation == "=") & (h.standard_value > 0) & h.data_validity_comment.isna()
          & (h.target_organism == "Homo sapiens")].copy()
    f["pIC50"] = 9 - np.log10(f.standard_value.astype(float)); f["ik"] = f.canonical_smiles.map(ik)
    med = f.dropna(subset=["ik"]).groupby("ik").pIC50.median()
    dm["pIC50"] = dm.ik.map(med)
    lab = dm.dropna(subset=["pIC50"]).groupby("key").agg(pIC50=("pIC50", "median"), smiles=("canonical_smiles", "first"))
    sig, files = signatures(list(lab.index))
    keys = [k for k in lab.index if k in sig]; missing = [k for k in lab.index if k not in sig]
    # common gene set across signature files; genes measured at top dose in >= 90% of drugs
    G = sorted(set.intersection(*[set(sig[k][1]) for k in keys]))
    T3 = []
    for k in keys:
        pos = {g: i for i, g in enumerate(sig[k][1])}; T3.append(sig[k][0][:, [pos[g] for g in G]])
    T3 = np.array(T3)
    keep = (T3[:, 2, :] != 0).mean(0) >= 0.9; T = T3[:, :, keep].reshape(len(keys), -1)
    Xc = np.nan_to_num(np.array([t1.chem(lab.loc[k, "smiles"]) for k in keys], float), posinf=0, neginf=0)
    y = lab.loc[keys, "pIC50"].values; yb = (y >= 5).astype(int)
    print(f"hERG: {len(lab)} Tahoe drugs with pIC50; analysed {len(keys)} (missing signatures: {len(missing)}); pIC50>=5: {yb.sum()}", flush=True)
    nc = Xc.shape[1]

    def oof(kind):
        cv = (RepeatedKFold(n_splits=5, n_repeats=10, random_state=2026).split(Xc) if kind == "reg"
              else RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=2026).split(Xc, yb))
        P = {"chem": np.zeros(len(y)), "chem+tahoe": np.zeros(len(y))}; n = np.zeros(len(y))
        for tr, te in cv:
            sc = StandardScaler().fit(Xc[tr]); A, B = sc.transform(Xc[tr]), sc.transform(Xc[te])
            st = StandardScaler().fit(T[tr]); pca = PCA(10, random_state=0).fit(st.transform(T[tr]))
            At = np.hstack([A, pca.transform(st.transform(T[tr]))]); Bt = np.hstack([B, pca.transform(st.transform(T[te]))])
            s2 = StandardScaler().fit(At); At, Bt = s2.transform(At), s2.transform(Bt)
            if kind == "reg":
                P["chem"][te] += Ridge(alpha=100).fit(A, y[tr]).predict(B); P["chem+tahoe"][te] += Ridge(alpha=100).fit(At, y[tr]).predict(Bt)
            else:
                P["chem"][te] += LogisticRegression(C=0.01, max_iter=5000).fit(A, yb[tr]).predict_proba(B)[:, 1]
                P["chem+tahoe"][te] += LogisticRegression(C=0.01, max_iter=5000).fit(At, yb[tr]).predict_proba(Bt)[:, 1]
            n[te] += 1
        return {k: v / n for k, v in P.items()}

    rng = np.random.default_rng(2026); res = []
    Pr = oof("reg"); Pc = oof("clf")
    boots = [rng.integers(0, len(y), len(y)) for _ in range(2000)]
    m = lambda p, i: spearmanr(p[i], y[i])[0]
    d = np.array([m(Pr["chem+tahoe"], i) - m(Pr["chem"], i) for i in boots])
    res.append(dict(task="regression_pIC50_spearman", n=len(y), chem=spearmanr(Pr["chem"], y)[0], chem_tahoe=spearmanr(Pr["chem+tahoe"], y)[0],
                    delta=spearmanr(Pr["chem+tahoe"], y)[0] - spearmanr(Pr["chem"], y)[0], ci_lo=np.nanpercentile(d, 2.5), ci_hi=np.nanpercentile(d, 97.5)))
    bb = [i for i in boots if yb[i].min() != yb[i].max()]
    d2 = np.array([roc_auc_score(yb[i], Pc["chem+tahoe"][i]) - roc_auc_score(yb[i], Pc["chem"][i]) for i in bb])
    res.append(dict(task="classification_pIC50>=5_auroc", n=len(y), chem=roc_auc_score(yb, Pc["chem"]), chem_tahoe=roc_auc_score(yb, Pc["chem+tahoe"]),
                    delta=roc_auc_score(yb, Pc["chem+tahoe"]) - roc_auc_score(yb, Pc["chem"]), ci_lo=np.percentile(d2, 2.5), ci_hi=np.percentile(d2, 97.5)))
    R = pd.DataFrame(res); R.to_csv(os.path.join(OUT, "summary.csv"), index=False)
    pd.DataFrame({"key": keys, "pIC50": y, "oof_chem": Pr["chem"], "oof_chem_tahoe": Pr["chem+tahoe"]}).to_csv(os.path.join(OUT, "oof_regression.csv"), index=False)
    prim = R.iloc[0]
    txt = (f"T2-H2 hERG (ChEMBL CHEMBL240): {len(lab)} Tahoe drugs with exact pIC50; {len(keys)} analysed; missing signatures: {', '.join(missing) or 'none'}\n"
           f"signature files used: {len(files)}\n"
           f"PRIMARY delta Spearman (chem+Tahoe - chem) = {prim.delta:+.3f} (95% CI {prim.ci_lo:+.3f} to {prim.ci_hi:+.3f})\n"
           "DECISION: Tahoe " + ("ADDS information (lower CI > 0)" if prim.ci_lo > 0 else "does NOT add information (lower CI <= 0)")
           + "\n\n" + R.round(3).to_string(index=False))
    open(os.path.join(OUT, "summary.txt"), "w").write(txt); print(txt)

if __name__ == "__main__":
    main()

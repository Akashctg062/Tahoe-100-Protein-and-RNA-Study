"""
v2 / T1-A2 (Amendment 3, frozen 18 Oct 2026) - best-effort models for within-drug sensitivity prediction.
Sets: Response, RNA, Protein_MS, RPPA, Response+RNA, Response+Protein_MS; each block -> PCA-10 inside each fold.
Models: RidgeCV (GCV, alpha 1e-2..1e3, 12 values), ElasticNetCV (5-fold, l1 0.2/0.5/0.8, 30 alphas),
        RandomForest (200 trees, leaf 3, max_features 0.33), RBF kernel ridge (alpha 1, gamma 1/p).
Evaluation: within-drug 5-fold CV (fixed seed); Spearman(pred, AUC) calibrated by 50 within-drug permutations.
Run:  python T1A2_best_effort_sensitivity.py [--models ridge,enet,rf,krr] [--n-perm 50]
"""
import os, sys, pickle, time, argparse, warnings
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
from sklearn.model_selection import KFold
from sklearn.linear_model import RidgeCV, ElasticNetCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.kernel_ridge import KernelRidge
warnings.filterwarnings("ignore")

SETS = {"Response": ["Response"], "RNA": ["RNA"], "Protein_MS": ["Protein_MS"], "RPPA": ["RPPA"],
        "Response+RNA": ["Response", "RNA"], "Response+Protein_MS": ["Response", "Protein_MS"]}
PAIRS = [("Response", "RNA"), ("Response", "Protein_MS"), ("Protein_MS", "RNA"), ("RPPA", "RNA"),
         ("Response+RNA", "RNA"), ("Response+Protein_MS", "Protein_MS"), ("RPPA", "Protein_MS")]

def rank(x): return pd.Series(x).rank().values
def spear(a, b):
    if np.std(a) == 0 or np.std(b) == 0: return np.nan
    return np.corrcoef(rank(a), rank(b))[0, 1]
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)

def fold_pca(X, folds, k=10):
    """label-free: per fold, standardise + PCA on training rows; returns list of (Ztr, Zte)."""
    out = []
    for tr, te in folds:
        mu, sd = X[tr].mean(0), X[tr].std(0); sd[sd == 0] = 1
        Z = (X[tr] - mu) / sd; z = (X[te] - mu) / sd
        U, S, Vt = np.linalg.svd(Z, full_matrices=False); V = Vt[:min(k, Vt.shape[0])].T
        out.append((Z @ V, z @ V))
    return out

def make(model, p):
    if model == "ridge": return RidgeCV(alphas=np.logspace(-2, 3, 12))
    if model == "enet": return ElasticNetCV(l1_ratio=[0.2, 0.5, 0.8], alphas=30, cv=5, max_iter=5000, random_state=0)
    if model == "rf": return RandomForestRegressor(200, min_samples_leaf=3, max_features=0.33, random_state=0, n_jobs=1)
    if model == "krr": return KernelRidge(alpha=1.0, kernel="rbf", gamma=1.0 / p)

def cv_predict(blocks, folds, y, model):
    pred = np.empty(len(y))
    for f, (tr, te) in enumerate(folds):
        Ztr = np.hstack([b[f][0] for b in blocks]); Zte = np.hstack([b[f][1] for b in blocks])
        m, s = Ztr.mean(0), Ztr.std(0); s[s == 0] = 1
        A, B = (Ztr - m) / s, (Zte - m) / s
        mdl = make(model, A.shape[1]); ym = y[tr].mean()
        mdl.fit(A, y[tr] - ym); pred[te] = mdl.predict(B) + ym
    return pred

def main(a):
    data = pickle.load(open(a.cache, "rb")); os.makedirs(a.out, exist_ok=True)
    models = a.models.split(","); t0 = time.time()
    part = os.path.join(a.out, "per_drug_" + "_".join(models) + ".csv")
    done = pd.read_csv(part) if os.path.exists(part) else pd.DataFrame()
    rows = done.to_dict("records")
    for src in ["PRISM", "GDSC"]:
        for d in data[src]:
            if len(done) and ((done.source == src) & (done.drug_key == d["drug"])).any(): continue
            y = d["y"]; n = len(y)
            folds = list(KFold(5, shuffle=True, random_state=2026).split(np.arange(n)))
            rng = np.random.default_rng(abs(hash((src, d["drug"]))) % (2**32) if a.seed_by_drug else 2026)
            perms = [rng.permutation(y) for _ in range(a.n_perm)]
            P = {b: fold_pca(d["X"][b], folds, 10) for b in ["Response", "RNA", "Protein_MS", "RPPA"]}
            for model in models:
                for s, comp in SETS.items():
                    blocks = [P[c] for c in comp]
                    obs = spear(cv_predict(blocks, folds, y, model), y)
                    null = np.array([spear(cv_predict(blocks, folds, yp, model), yp) for yp in perms])
                    rows.append(dict(source=src, drug_key=d["drug"], model=model, set=s, n_lines=n, rho=obs,
                                     null_mean=np.nanmean(null), calibrated_rho=obs - np.nanmean(null)))
            pd.DataFrame(rows).to_csv(part, index=False)
            print(f"{src} {d['drug']}: done ({(time.time()-t0)/60:.1f} min)", flush=True)
    D = pd.DataFrame(rows)
    # ---- summaries per model family (BH within family across sets + pairs, both screens)
    S, C = [], []
    for model in models:
        Dm = D[D.model == model]
        for (src, s), g in Dm.groupby(["source", "set"]):
            S.append(dict(model=model, source=src, set=s, n_drugs=len(g), median_calibrated=g.calibrated_rho.median(),
                          p_gt0=wilcoxon(g.calibrated_rho.dropna(), alternative="greater").pvalue))
        for src in Dm.source.unique():
            W = Dm[Dm.source == src].pivot(index="drug_key", columns="set", values="calibrated_rho")
            for x, z in PAIRS:
                dd = (W[x] - W[z]).dropna()
                C.append(dict(model=model, source=src, comparison=f"{x} vs {z}", n_drugs=len(dd),
                              median_difference=dd.median(), p=wilcoxon(dd).pvalue))
    S, C = pd.DataFrame(S), pd.DataFrame(C)
    S["fdr_bh"] = np.nan; C["fdr_bh"] = np.nan
    for model in models:
        i, j = S.model == model, C.model == model
        q = bh(np.r_[S.loc[i, "p_gt0"].values, C.loc[j, "p"].fillna(1).values])
        S.loc[i, "fdr_bh"] = q[:i.sum()]; C.loc[j, "fdr_bh"] = q[i.sum():]
    S.to_csv(os.path.join(a.out, "summary_sets_" + "_".join(models) + ".csv"), index=False)
    C.to_csv(os.path.join(a.out, "summary_pairs_" + "_".join(models) + ".csv"), index=False)
    # ---- frozen decision rules
    lines = ["T1-A2 decisions (Amendment 3)"]
    rescued = []
    for model in models:
        ok = []
        for comp in ["Response+RNA vs RNA", "Response+Protein_MS vs Protein_MS"]:
            r = C[(C.model == model) & (C.comparison == comp)]
            if len(r) == 2 and (r.fdr_bh < 0.05).all() and (r.median_difference > 0).all(): ok.append(comp)
        if ok: rescued.append(f"{model}: {ok}")
    lines.append("(a) Response rescued: " + ("YES - " + "; ".join(rescued) if rescued else "NO (no family met the rule in both screens)"))
    pos = [m for m in models if (C[(C.model == m) & (C.source == "GDSC") & (C.comparison == "RPPA vs RNA")].median_difference > 0).any()]
    lines.append(f"(b) RPPA - RNA median > 0 in GDSC for {len(pos)}/{len(models)} families ({', '.join(pos)}) -> "
                 + ("model-robust" if len(pos) >= 3 else "NOT model-robust (needs >= 3 of 4)"))
    txt = "\n".join(lines) + "\n\n" + S.round(4).to_string(index=False) + "\n\n" + C.round(4).to_string(index=False)
    open(os.path.join(a.out, "summary_" + "_".join(models) + ".txt"), "w").write(txt); print(txt)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default="/home/claude/v2/results/cache/c3_cache.pkl")
    ap.add_argument("--out", default="/home/claude/v2/results/T1A2")
    ap.add_argument("--models", default="ridge,enet,rf,krr")
    ap.add_argument("--n-perm", type=int, default=50)
    ap.add_argument("--seed-by-drug", action="store_true")
    main(ap.parse_args())

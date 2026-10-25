"""
v2 / T1-A1 (Amendment 3, frozen 18 Oct 2026) - best-effort models for cardiotoxicity (C1 FAERS primary; D1 secondary).
Representations: Tahoe PCA-50 (3-dose concatenated mean signature; PCA fitted in fold) and Tahoe programs+potency;
chemistry = RDKit 2-D descriptors + Morgan(r=2, 1024).  Families tuned by inner stratified 5-fold CV (AUROC):
elastic-net logistic, histogram gradient boosting, RBF-SVM, MLP.  Outer: stratified 5-fold x 10 repeats.
Primary per family x representation: dAUROC = AUROC(chem+Tahoe) - AUROC(chem), 2,000 drug bootstraps,
one-sided bootstrap p = P(d* <= 0), Holm over 8 combinations; rescued if any adjusted p < 0.025 (C1).
Run:  python T1A1_best_effort_cardiotox.py --task C1     (then --task D1)
"""
import os, re, sys, time, argparse, warnings
import numpy as np, pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold, RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")

from sklearn.base import BaseEstimator, TransformerMixin
class GramPCA(BaseEstimator, TransformerMixin):
    """Standardise + PCA computed exactly from the n x n Gram matrix (n samples << p features).
    Numerically equivalent to Pipeline(StandardScaler(), PCA(k)) up to the sign of each component."""
    def __init__(self, k=50): self.k = k
    def fit(self, X, y=None):
        X = np.asarray(X, float); self.mu_ = X.mean(0); sd = X.std(0); sd[sd == 0] = 1; self.sd_ = sd
        Z = (X - self.mu_) / self.sd_; w, U = np.linalg.eigh(Z @ Z.T)
        o = np.argsort(w)[::-1][:min(self.k, len(w))]; w, U = np.clip(w[o], 1e-12, None), U[:, o]
        self.V_ = Z.T @ (U / np.sqrt(w)); return self
    def transform(self, X): return ((np.asarray(X, float) - self.mu_) / self.sd_) @ self.V_
U = "/mnt/user-data/uploads/Protein and RNA/"

PROGS = {
 "OXPHOS": r"^(NDUF[ABSCV]\d|SDH[ABCD]$|UQCR|COX[4-8]|ATP5)", "mito_ribo": r"^MRP[LS]\d", "cyto_ribo": r"^RP[LS]\d+[AXY]?$",
 "cholesterol": r"^(HMGCR|HMGCS1|LDLR|INSIG1|SQLE|IDI1|FDFT1|MVD|MVK|LSS|CYP51A1|DHCR7|DHCR24|SC5D|MSMO1|FDPS|STARD4|NSDHL|HSD17B7)$",
 "UPR_ERstress": r"^(DDIT3|ATF4|ATF3|HSPA5|XBP1|TRIB3|CHAC1|ASNS|HERPUD1|SESN2|DNAJB9|SEL1L|PDIA4|HYOU1|CALR|VEGFA|STC2)$",
 "p53_DNAdamage": r"^(CDKN1A|MDM2|GDF15|TP53I3|FAS|BAX|BBC3|SESN1|RRM2B|TNFRSF10B|ZMAT3|PLK3|FDXR|AEN|TP53INP1|DDB2|XPC|GADD45A)$",
 "cell_cycle": r"^(MKI67|TOP2A|CCNB1|CCNB2|CDK1|CCNA2|BUB1|PLK1|AURKA|AURKB|CENPF|E2F1|MCM[2-7]|PCNA|TYMS|RRM2|CDC20|KIF11|TPX2)$",
 "NRF2_oxstress": r"^(HMOX1|NQO1|TXNRD1|GCLM|GCLC|SRXN1|SLC7A11|OSGIN1|AKR1C[123]|G6PD|ME1|FTL|FTH1)$",
 "metallothionein": r"^MT[12][A-Z]?$", "histone": r"^H[1234][A-C]?[0-9A-Z-]*$|^HIST",
 "hypoxia_glyco": r"^(SLC2A1|PGK1|LDHA|ENO1|ALDOA|PFKFB3|BNIP3|CA9|NDRG1|P4HA1|EGLN3|ANKRD37)$",
 "ion_channel": r"^(KCN|SCN\d|CACN|HCN\d|RYR|ATP2A|PLN$|SLC8A)", "sarcomere_cytoskel": r"^(MYL\d|MYH\d|TNN[TIC]|ACTN\d|DES$|TTN$|ACTC1|TPM\d)",
 "IFN": r"^(ISG15|IFIT\d|IFI\d+|MX[12]$|OAS[123L]|IRF7|STAT1|BST2|RSAD2|XAF1)$",
 "NFkB_inflam": r"^(NFKBIA|TNFAIP3|CXCL\d+|IL6$|IL1[AB]$|CCL\d+|ICAM1|RELB|NFKB2|BIRC3|IER3)$"}

def chem(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    mg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    m = Chem.MolFromSmiles(str(smiles))
    if m is None: return None
    return [f(m) for _, f in Descriptors.descList] + list(mg.GetFingerprintAsNumPy(m))

def programs(T3, genes):
    P = []
    for n, rx in PROGS.items():
        msk = np.array([bool(re.match(rx, str(g))) for g in genes])
        if msk.sum() < 3: continue
        for d in range(3): P.append(T3[:, d, msk].mean(1))
    pot = [np.sqrt((T3[:, d, :] ** 2).mean(1)) for d in range(3)] + [(np.abs(T3[:, d, :]) > 3).mean(1) for d in range(3)]
    return np.column_stack(P + pot)

def load(task):
    dm = pd.read_csv(U + "data/raw/tahoe_metadata/drug_metadata.csv").drop_duplicates("drug")
    if task == "C1":
        c1 = pd.read_csv("/home/claude/ana/c1/c1_with_indication.csv")
        z = np.load(U + "results/confirmatory/C1/tahoe_features_C1.npz", allow_pickle=True)
        X3, keys, genes = z["X"], list(z["keys"]), list(z["genes"])
        smi = dict(zip(dm.drug, dm.canonical_smiles))
        ks = list(c1.key); y = c1.signal.values.astype(int); S = [smi.get(d) for d in c1.drug_tahoe]
    else:
        lab = pd.read_csv("/home/claude/ana/labels_meta.csv").drop_duplicates("key").set_index("key")
        X3 = np.load("/home/claude/ana/X_all.npy"); keys = pd.read_csv("/home/claude/ana/keys_all.csv").iloc[:, 0].tolist()
        genes = pd.read_csv("/home/claude/ana/genes_all.csv").iloc[:, 0].astype(str).tolist()
        sub = lab[lab.dict_concern.isin(["most", "no"])]
        ks = [k for k in sub.index if k in set(keys)]; y = (sub.loc[ks, "dict_concern"] == "most").astype(int).values
        S = list(sub.loc[ks, "canonical_smiles"])
    ki = {k: i for i, k in enumerate(keys)}
    keep = (X3[:, 2, :] != 0).mean(0) >= 0.9
    T3 = X3[:, :, keep]; G = [g for g, k in zip(genes, keep) if k]
    rows = [(k, chem(s), ki.get(k)) for k, s in zip(ks, S)]
    ok = [i for i, (k, c, t) in enumerate(rows) if c is not None and t is not None]
    Xc = np.nan_to_num(np.array([rows[i][1] for i in ok], float), posinf=0, neginf=0)
    idx = [rows[i][2] for i in ok]
    Tflat = T3[idx].reshape(len(idx), -1); Tprog = programs(T3[idx], G)
    return [ks[i] for i in ok], y[ok], Xc, Tflat, Tprog

FAMILIES = {
    "enet_logistic": (LogisticRegression(solver="saga", max_iter=5000, tol=1e-3, l1_ratio=0.5, C=0.1),
                      {"m__C": [0.001, 0.01, 0.1, 1], "m__l1_ratio": [0.1, 0.5, 0.9]}),
    "grad_boost": (HistGradientBoostingClassifier(max_iter=200, random_state=0),
                   {"m__learning_rate": [0.03, 0.1], "m__max_depth": [2, 3]}),
    "rbf_svm": (SVC(kernel="rbf", gamma="scale", probability=True, random_state=0), {"m__C": [0.1, 1, 10]}),
    "mlp": (MLPClassifier(early_stopping=True, max_iter=2000, random_state=0),
            {"m__hidden_layer_sizes": [(16,), (64,)], "m__alpha": [0.001, 0.1]}),
}

def pipe(model, nc, nt, tahoe_pca, memory=None):
    blocks = []
    if nc: blocks.append(("chem", StandardScaler(), list(range(nc))))
    if nt:
        t = GramPCA(50) if tahoe_pca else StandardScaler()
        blocks.append(("tahoe", t, list(range(nc, nc + nt))))
    return Pipeline([("prep", ColumnTransformer(blocks)), ("m", model)], memory=memory)

def main(a):
    keys, y, Xc, Tflat, Tprog = load(a.task); t0 = time.time()
    print(f"{a.task}: {len(y)} drugs ({y.sum()} positive); chem {Xc.shape[1]}, Tahoe flat {Tflat.shape[1]}, programs {Tprog.shape[1]}", flush=True)
    sets = {"chemistry": (Xc, Xc.shape[1], 0, False),
            "tahoe_pca50": (Tflat, 0, Tflat.shape[1], True), "tahoe_programs": (Tprog, 0, Tprog.shape[1], False),
            "chem+tahoe_pca50": (np.hstack([Xc, Tflat]), Xc.shape[1], Tflat.shape[1], True),
            "chem+tahoe_programs": (np.hstack([Xc, Tprog]), Xc.shape[1], Tprog.shape[1], False)}
    os.makedirs(a.out, exist_ok=True)
    oof_file = os.path.join(a.out, f"oof_{a.task}.csv")
    O = pd.read_csv(oof_file) if os.path.exists(oof_file) else pd.DataFrame({"key": keys, "y": y})
    outer = list(RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=2026).split(Xc, y))
    for fam, (model, grid) in FAMILIES.items():
        for s, (X, nc, nt, tp) in sets.items():
            col = f"{fam}|{s}"
            if col in O.columns: continue
            pred = np.zeros(len(y)); cnt = np.zeros(len(y))
            from sklearn.model_selection import ParameterGrid
            from sklearn.base import clone
            params = list(ParameterGrid({k.replace("m__", ""): v for k, v in grid.items()}))
            for tr, te in outer:
                # nested CV written out: preprocessing fitted once per inner split, then every grid point (= GridSearchCV)
                scores = np.zeros(len(params)); inner = StratifiedKFold(5, shuffle=True, random_state=0).split(X[tr], y[tr])
                for itr, iva in inner:
                    prep = pipe(model, nc, nt, tp).named_steps["prep"].fit(X[tr][itr], y[tr][itr])
                    A, Bv = prep.transform(X[tr][itr]), prep.transform(X[tr][iva])
                    for j, pr in enumerate(params):
                        try:
                            mdl = clone(model).set_params(**pr).fit(A, y[tr][itr])
                            scores[j] += roc_auc_score(y[tr][iva], mdl.predict_proba(Bv)[:, 1]) / 5
                        except Exception:
                            scores[j] += 0.5 / 5
                best = params[int(np.argmax(scores))]
                prep = pipe(model, nc, nt, tp).named_steps["prep"].fit(X[tr], y[tr])
                mdl = clone(model).set_params(**best).fit(prep.transform(X[tr]), y[tr])
                pred[te] += mdl.predict_proba(prep.transform(X[te]))[:, 1]; cnt[te] += 1
            O[col] = pred / cnt; O.to_csv(oof_file, index=False)
            print(f"  {col}: AUROC {roc_auc_score(y, O[col]):.3f}  ({(time.time()-t0)/60:.1f} min)", flush=True)
    # ---- statistics
    rng = np.random.default_rng(2026); B = 2000; boots = []
    while len(boots) < B:
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max(): boots.append(i)
    res = []
    for fam in FAMILIES:
        for rep in ["pca50", "programs"]:
            c, t, ct = O[f"{fam}|chemistry"].values, O[f"{fam}|tahoe_{rep}"].values, O[f"{fam}|chem+tahoe_{rep}"].values
            au = lambda p, i=slice(None): roc_auc_score(y[i], p[i])
            d_b = np.array([au(ct, i) - au(c, i) for i in boots])
            t_b = np.array([au(t, i) for i in boots])
            perm = np.array([roc_auc_score(rng.permutation(y), t) for _ in range(1000)])
            res.append(dict(family=fam, representation=rep, auroc_chem=au(c), auroc_tahoe=au(t), auroc_chem_tahoe=au(ct),
                            tahoe_ci_lo=np.percentile(t_b, 2.5), tahoe_ci_hi=np.percentile(t_b, 97.5),
                            tahoe_perm_p=(np.sum(perm >= au(t)) + 1) / 1001,
                            delta=au(ct) - au(c), delta_ci_lo=np.percentile(d_b, 2.5), delta_ci_hi=np.percentile(d_b, 97.5),
                            p_one_sided=(np.sum(d_b <= 0) + 1) / (B + 1)))
    R = pd.DataFrame(res)
    o = np.argsort(R.p_one_sided.values); m = len(R); adj = np.empty(m); run = 0
    for r, i in enumerate(o): run = max(run, min(1, (m - r) * R.p_one_sided.values[i])); adj[i] = run
    R["p_holm"] = adj
    R.to_csv(os.path.join(a.out, f"summary_{a.task}.csv"), index=False)
    dec = "RESCUED (some Holm-adjusted p < 0.025)" if (R.p_holm < 0.025).any() else "NOT rescued (all Holm-adjusted p >= 0.025)"
    txt = f"T1-A1 {a.task}: {len(y)} drugs, {y.sum()} positive\nDecision (frozen rule applies to C1): Tahoe {dec}\n\n" + R.round(3).to_string(index=False)
    open(os.path.join(a.out, f"summary_{a.task}.txt"), "w").write(txt); print(txt)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", default="C1", choices=["C1", "D1"])
    ap.add_argument("--out", default="/home/claude/v2/results/T1A1")
    main(ap.parse_args())

"""
Step 5 - Baselines FIRST. Predict cardiotoxicity (DICTrank) for Tahoe drugs from
  (a) chemical structure only   <- the strong published baseline (Seal et al. 2024, JCIM)
  (b) tumour-cell transcriptional response from Tahoe-100M
  (c) both combined
Your scientific claim is about (b) and (c): does BIOLOGY add signal beyond chemistry?

Inputs
  ../data/processed/shared_drugs.csv               (step 2, has dict_concern)
  ../data/processed/drug_gene_features.parquet      (step 3b)
  ../data/raw/tahoe_metadata/drug_metadata.csv      (step 1, has SMILES)
Run:  python 05_baseline_models.py            (Most-concern vs No-concern)
      python 05_baseline_models.py --positive most_less
      python 05_baseline_models.py --scope lung        (lung-only tumour features)
Output: ../results/baseline_table_<positive>_<scope>.csv and oof_predictions_<positive>_<scope>.csv
"""
import argparse, os
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score

ap = argparse.ArgumentParser(); ap.add_argument("--positive", default="most", choices=["most", "most_less"])
ap.add_argument("--scope", default="all", choices=["all", "lung"])   # tumour features from all lines or lung only
a = ap.parse_args()

lab = pd.read_csv("../data/processed/shared_drugs.csv")
c = lab.dict_concern.astype(str).str.strip().str.lower()
if a.positive == "most":
    lab = lab[c.isin(["most", "no"])]; y_all = lab.dict_concern.str.strip().str.lower() == "most"
else:
    lab = lab[c.isin(["most", "less", "no"])]; y_all = lab.dict_concern.str.strip().str.lower() != "no"
lab = lab.assign(y=y_all.astype(int).values)
name2key = dict(zip(lab.drug_tahoe, lab.key))
lab = lab.drop_duplicates("key").set_index("key")       # one row per COMPOUND (salt variants merged)

# (b) tumour response: mean signed DE statistic across cell lines at each drug's highest dose
feat = pd.read_parquet("../data/processed/drug_gene_features.parquet",
                       filters=[("scope", "==", a.scope)])        # from step 3b
feat = feat[feat.drug.isin(name2key)].assign(key=lambda x: x.drug.map(name2key))
top = feat.loc[feat.groupby("drug").concentration.transform("max") == feat.concentration]
Xb = top.pivot_table(index="key", columns="gene_name", values="mean_stat", aggfunc="mean").fillna(0.0)
Xb = Xb.loc[:, Xb.var().sort_values(ascending=False).index[:2000]]

# (a) structure: RDKit descriptors + Morgan fingerprint from Tahoe SMILES
dm = pd.read_csv("../data/raw/tahoe_metadata/drug_metadata.csv")
smi_col = next(col for col in dm.columns if "smiles" in col.lower())
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
def featurise(s):
    m = Chem.MolFromSmiles(str(s))
    if m is None: return None
    desc = [f(m) for _, f in Descriptors.descList]
    fp = list(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024))
    return desc + fp
dm = dm[dm.drug.isin(name2key)].assign(key=lambda x: x.drug.map(name2key)).drop_duplicates("key")
feats = {r.key: featurise(r[smi_col]) for _, r in dm.iterrows()}
Xa = pd.DataFrame({k: v for k, v in feats.items() if v is not None}).T
Xa = Xa.replace([np.inf, -np.inf], np.nan).fillna(0.0)

common = lab.index.intersection(Xb.index).intersection(Xa.index)
y = lab.loc[common, "y"].values
sets = {"structure": Xa.loc[common].values, "tumour_response": Xb.loc[common].values,
        "combined": np.hstack([Xa.loc[common].values, Xb.loc[common].values])}
print(f"{len(y)} labelled drugs ({y.sum()} positive)")

models = {"majority": lambda: DummyClassifier(strategy="prior"),
          "logistic": lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.05, max_iter=5000)),
          "random_forest": lambda: RandomForestClassifier(n_estimators=500, min_samples_leaf=2, random_state=0)}
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)
rows, oof = [], pd.DataFrame(index=common)
for fs, X in sets.items():
    for mn, mk in models.items():
        aucs, aps, pred = [], [], np.zeros(len(y)); cnt = np.zeros(len(y))
        for tr, te in cv.split(X, y):
            m = mk().fit(X[tr], y[tr]); p = m.predict_proba(X[te])[:, 1]
            aucs.append(roc_auc_score(y[te], p)); aps.append(average_precision_score(y[te], p))
            pred[te] += p; cnt[te] += 1
        oof[f"{fs}|{mn}"] = pred / cnt
        rows.append(dict(features=fs, model=mn, AUROC=np.mean(aucs), AUROC_sd=np.std(aucs),
                         AUPRC=np.mean(aps), AUPRC_sd=np.std(aps), prevalence=y.mean()))
        print(f"{fs:16s} {mn:14s} AUROC={np.mean(aucs):.3f}±{np.std(aucs):.3f}  AUPRC={np.mean(aps):.3f}")
os.makedirs("../results", exist_ok=True)
tag = f"{a.positive}_{a.scope}"
oof["label"] = y; oof.to_csv(f"../results/oof_predictions_{tag}.csv")
pd.DataFrame(rows).to_csv(f"../results/baseline_table_{tag}.csv", index=False)
print(f"saved ../results/baseline_table_{tag}.csv and ../results/oof_predictions_{tag}.csv")

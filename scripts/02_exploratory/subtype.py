import pandas as pd, numpy as np, re, warnings, importlib.util; warnings.filterwarnings("ignore")
from rdkit import Chem, RDLogger; RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
spec=importlib.util.spec_from_file_location("m","/home/claude/proj/out/scripts/02_drug_overlap.py"); mm=importlib.util.module_from_spec(spec); spec.loader.exec_module(mm)
d=pd.read_excel("/mnt/user-data/uploads/Protein and RNA/data/raw/dictrank_dataset_508.xlsx"); d["key"]=d["Generic/Proper Name(s)"].astype(str).map(mm.norm)
sev={"most":3,"less":2,"ambiguous":1,"no":0}; d["sev"]=d["DICT _ Concern"].astype(str).str.strip().str.lower().map(sev)
d=d.sort_values("sev",ascending=False).drop_duplicates("key")
lab=pd.read_csv("labels_meta.csv").drop_duplicates("key").merge(d[["key","Cardiotoxicity"]],on="key",how="left").set_index("key")
lab["ctype"]=lab.Cardiotoxicity.astype(str).str.lower().str.strip()
X3=np.load("X_all.npy"); keys=pd.read_csv("keys_all.csv").iloc[:,0].tolist(); kidx={k:i for i,k in enumerate(keys)}
keep=(X3[:,2,:]!=0).mean(0)>=0.9; T=X3[:,:,keep]
mg=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024)
def feat(s):
    m=Chem.MolFromSmiles(str(s)); return None if m is None else [f(m) for _,f in Descriptors.descList]+list(mg.GetFingerprintAsNumPy(m))
def evaluate(pos_mask, name):
    sub=lab[(pos_mask | (lab.dict_concern=="no"))]
    ks=[k for k in sub.index if k in kidx]; S=[feat(lab.loc[k,"canonical_smiles"]) for k in ks]
    ks=[k for k,s in zip(ks,S) if s is not None]; Xs=np.nan_to_num(np.array([feat(lab.loc[k,"canonical_smiles"]) for k in ks],float),posinf=0,neginf=0)
    y=(lab.loc[ks,"dict_concern"]!="no").astype(int).values; idx=[kidx[k] for k in ks]
    Xb=np.hstack([T[idx,dd,:] for dd in range(3)])
    cv=RepeatedStratifiedKFold(n_splits=5,n_repeats=10,random_state=1)
    res={}
    for nm,X,pipe in [("structure",Xs,lambda: make_pipeline(StandardScaler(),LogisticRegression(C=0.01,max_iter=5000))),
                      ("tumour PCA10",Xb,lambda: make_pipeline(StandardScaler(),PCA(10,random_state=0),LogisticRegression(C=0.1,max_iter=5000)))]:
        a=[roc_auc_score(y[te],pipe().fit(X[tr],y[tr]).predict_proba(X[te])[:,1]) for tr,te in cv.split(X,y)]
        res[nm]=(np.mean(a),np.percentile(a,2.5),np.percentile(a,97.5))
    print(f"{name:34s} n={len(y)} pos={y.sum()} | "+" | ".join(f"{k} AUROC {v[0]:.2f} [{v[1]:.2f}-{v[2]:.2f}]" for k,v in res.items()))
mostless=lab.dict_concern.isin(["most","less"])
evaluate(mostless & (lab.ctype=="heart damage"), "heart damage vs no")
evaluate(mostless & (lab.ctype=="arrhythmia"), "arrhythmia vs no")
evaluate(mostless & (lab.ctype=="mixed"), "mixed vs no")
evaluate(mostless & (lab.ctype=="heart damage") & lab.oncology, "heart damage (oncology) vs no")

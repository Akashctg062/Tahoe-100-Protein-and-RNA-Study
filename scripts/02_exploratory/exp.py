import pandas as pd, numpy as np, re, warnings, sys
warnings.filterwarnings("ignore")
from rdkit import Chem, RDLogger; RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import KNeighborsClassifier

SC=sys.argv[1] if len(sys.argv)>1 else "all"
lab=pd.read_csv("labels_meta.csv").drop_duplicates("key").set_index("key")
X3=np.load(f"X_{SC}.npy"); keys=pd.read_csv(f"keys_{SC}.csv").iloc[:,0].tolist(); genes=pd.read_csv(f"genes_{SC}.csv").iloc[:,0].astype(str).values
kidx={k:i for i,k in enumerate(keys)}
# keep genes measured (non-zero) in >=90% of drug-doses at top dose
top=X3[:,2,:]; keep=(top!=0).mean(0)>=0.9
print(SC,"genes kept:",keep.sum())
G=genes[keep]; T=X3[:,:,keep]

# ---- curated gene programs (prefix/regex) ----
progs={
 "OXPHOS":r"^(NDUF[ABSCV]\d|SDH[ABCD]$|UQCR|COX[4-8]|ATP5)",
 "mito_ribo":r"^MRP[LS]\d", "cyto_ribo":r"^RP[LS]\d+[AXY]?$",
 "cholesterol":r"^(HMGCR|HMGCS1|LDLR|INSIG1|SQLE|IDI1|FDFT1|MVD|MVK|LSS|CYP51A1|DHCR7|DHCR24|SC5D|MSMO1|FDPS|STARD4|NSDHL|HSD17B7)$",
 "UPR_ERstress":r"^(DDIT3|ATF4|ATF3|HSPA5|XBP1|TRIB3|CHAC1|ASNS|HERPUD1|SESN2|DNAJB9|SEL1L|PDIA4|HYOU1|CALR|VEGFA|STC2)$",
 "p53_DNAdamage":r"^(CDKN1A|MDM2|GDF15|TP53I3|FAS|BAX|BBC3|SESN1|RRM2B|TNFRSF10B|ZMAT3|PLK3|FDXR|AEN|TP53INP1|DDB2|XPC|GADD45A)$",
 "cell_cycle":r"^(MKI67|TOP2A|CCNB1|CCNB2|CDK1|CCNA2|BUB1|PLK1|AURKA|AURKB|CENPF|E2F1|MCM[2-7]|PCNA|TYMS|RRM2|CDC20|KIF11|TPX2)$",
 "NRF2_oxstress":r"^(HMOX1|NQO1|TXNRD1|GCLM|GCLC|SRXN1|SLC7A11|OSGIN1|AKR1C[123]|G6PD|ME1|FTL|FTH1)$",
 "metallothionein":r"^MT[12][A-Z]?$",
 "histone":r"^H[1234][A-C]?[0-9A-Z-]*$|^HIST",
 "hypoxia_glyco":r"^(SLC2A1|PGK1|LDHA|ENO1|ALDOA|PFKFB3|BNIP3|CA9|NDRG1|P4HA1|EGLN3|ANKRD37)$",
 "ion_channel":r"^(KCN|SCN\d|CACN|HCN\d|RYR|ATP2A|PLN$|SLC8A)",
 "sarcomere_cytoskel":r"^(MYL\d|MYH\d|TNN[TIC]|ACTN\d|DES$|TTN$|ACTC1|TPM\d)",
 "IFN":r"^(ISG15|IFIT\d|IFI\d+|MX[12]$|OAS[123L]|IRF7|STAT1|BST2|RSAD2|XAF1)$",
 "NFkB_inflam":r"^(NFKBIA|TNFAIP3|CXCL\d+|IL6$|IL1[AB]$|CCL\d+|ICAM1|RELB|NFKB2|BIRC3|IER3)$",
}
P=[];pn=[]
for n,rx in progs.items():
    m=np.array([bool(re.match(rx,g)) for g in G])
    if m.sum()<3: continue
    for d in range(3): P.append(T[:,d,m].mean(1)); pn.append(f"{n}_d{d}")
P=np.vstack(P).T
# potency scalars
pot=np.column_stack([np.sqrt((T[:,d,:]**2).mean(1)) for d in range(3)]+[(np.abs(T[:,d,:])>3).mean(1) for d in range(3)])

# structure
smi=lab["canonical_smiles"]
mg=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024)
def feat(s):
    m=Chem.MolFromSmiles(str(s))
    if m is None: return None
    return [f(m) for _,f in Descriptors.descList]+list(mg.GetFingerprintAsNumPy(m))
S={k:feat(smi[k]) for k in keys if k in smi.index}
def run(task):
    if task=="most_vs_no": sub=lab[lab.dict_concern.isin(["most","no"])]; y=(sub.dict_concern=="most")
    elif task=="mostless_vs_no": sub=lab[lab.dict_concern.isin(["most","less","no"])]; y=(sub.dict_concern!="no")
    elif task=="onc_most_vs_no": sub=lab[lab.dict_concern.isin(["most","no"])&lab.oncology]; y=(sub.dict_concern=="most")
    elif task=="onc_mostless_vs_no": sub=lab[lab.dict_concern.isin(["most","less","no"])&lab.oncology]; y=(sub.dict_concern!="no")
    ks=[k for k in sub.index if k in kidx and S.get(k) is not None]
    y=y.loc[ks].astype(int).values; idx=[kidx[k] for k in ks]
    Xs=np.nan_to_num(np.array([S[k] for k in ks],dtype=float),posinf=0,neginf=0)
    topv=T[idx,2,:]; 
    sets={
     "structure":Xs,
     "bio_raw2000":topv[:,np.argsort(-topv.var(0))[:2000]],
     "bio_PCA10_alldoses":np.hstack([T[idx,d,:] for d in range(3)]),   # PCA applied in-fold
     "bio_programs":P[idx],
     "bio_programs+potency":np.hstack([P[idx],pot[idx]]),
     "structure+programs":np.hstack([Xs,P[idx],pot[idx]]),
    }
    cv=RepeatedStratifiedKFold(n_splits=5,n_repeats=4,random_state=0)
    out={}
    for name,X in sets.items():
        for mname in ["logit","rf"]:
            aucs=[]
            for tr,te in cv.split(X,y):
                if mname=="logit":
                    steps=[StandardScaler()]
                    if name=="bio_PCA10_alldoses": steps.append(PCA(10,random_state=0))
                    m=make_pipeline(*steps,LogisticRegression(C=0.1 if X.shape[1]<200 else 0.01,max_iter=5000))
                else:
                    m=make_pipeline(PCA(10,random_state=0),RandomForestClassifier(200,min_samples_leaf=3,random_state=0,n_jobs=2)) if name=="bio_PCA10_alldoses" else RandomForestClassifier(200,min_samples_leaf=3,random_state=0,n_jobs=2)
                m.fit(X[tr],y[tr]); aucs.append(roc_auc_score(y[te],m.predict_proba(X[te])[:,1]))
            out[(name,mname)]=(np.mean(aucs),np.std(aucs))
    print(f"\n== {SC} | {task}: n={len(y)} pos={y.sum()}")
    for (n,mn),(a,s) in out.items(): print(f"  {n:24s} {mn:6s} AUROC {a:.3f} ± {s:.3f}")
for t in ["most_vs_no","mostless_vs_no","onc_most_vs_no","onc_mostless_vs_no"]: run(t)

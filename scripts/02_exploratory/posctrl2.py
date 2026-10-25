import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import balanced_accuracy_score
lab=pd.read_csv("labels_meta.csv").drop_duplicates("key").set_index("key")
X3=np.load("X_all.npy"); keys=pd.read_csv("keys_all.csv").iloc[:,0].tolist()
keep=(X3[:,2,:]!=0).mean(0)>=0.9; X3=X3[:,:,keep]
# 1) dose-replicate check: corr(mid, high) same drug vs other drugs
def cor(a,b): a=a-a.mean(); b=b-b.mean(); return (a@b)/np.sqrt((a@a)*(b@b))
pot=np.sqrt((X3[:,2,:]**2).mean(1))
same=[cor(X3[i,1],X3[i,2]) for i in range(len(keys))]
other=[cor(X3[i,1],X3[(i+7)%len(keys),2]) for i in range(len(keys))]
print(f"mid-vs-high dose correlation: same drug median {np.median(same):.2f} | different drug median {np.median(other):.2f}")
q=pd.qcut(pot,3,labels=["weak","medium","strong"])
print(pd.DataFrame({"same":same,"pot":q}).groupby("pot").same.median().round(2).to_dict(),"(by response strength tertile)")
# 2) MoA recognition among oncology MoA classes, strong responders, correlation kNN on top-var genes
moa=lab.loc[keys,"moa-fine"].fillna("unclear").astype(str)
onc=["DNA synthesis/repair inhibitor","EGFR/ERBB inhibitor","Other TK inhibitor","Multi-TK inhibitor","Microtubule inhibitor","MEK inhibitor","RAF inhibitor","CDK inhibitor","HDAC inhibitor","Proteasome inhibitor","JAK/STAT inhibitor","Sonic inhibitor"]
m=moa.isin(onc).values
T=X3[m,2,:]; y=np.asarray(moa[m].values,dtype=object)
vc=pd.Series(y).value_counts(); ok=np.isin(y,vc[vc>=3].index); T,y=T[ok],y[ok]
g=np.argsort(-T.var(0))[:3000]; Z=T[:,g]; Z=(Z-Z.mean(1,keepdims=True))/Z.std(1,keepdims=True)
print("classes:",dict(pd.Series(y).value_counts()))
cv=StratifiedKFold(3,shuffle=True,random_state=0)
p=cross_val_predict(KNeighborsClassifier(1,metric="correlation"),Z,y,cv=cv)
rng=np.random.default_rng(0); null=[balanced_accuracy_score(rng.permutation(y),p) for _ in range(500)]
print(f"oncology MoA 1-NN: balanced acc {balanced_accuracy_score(y,p):.2f} | chance {1/len(set(y)):.2f} | perm 95th {np.percentile(null,95):.2f}")
for c in sorted(set(y)): print(f"   {c:32s} n={int((y==c).sum())} correct={int(((y==c)&(p==c)).sum())}")

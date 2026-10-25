import pandas as pd, numpy as np, pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow as pa
U="/mnt/user-data/uploads/Protein and RNA/"
s=pd.read_csv(U+"data/processed/shared_drugs.csv")
dm=pd.read_csv(U+"data/raw/tahoe_metadata/drug_metadata.csv").drop_duplicates("drug")
onc_re=r"cancer|carcinoma|leuk|lymphoma|myeloma|tumou?r|melanoma|sarcoma|neoplas|oncolog|malignan|nsclc|glioma|myelofibrosis|myelodysplas|gist"
dm["oncology"]=dm["gpt-notes-approval"].astype(str).str.lower().str.contains(onc_re)
s=s.merge(dm[["drug","oncology","moa-fine","targets","canonical_smiles"]],left_on="drug_tahoe",right_on="drug",how="left")
s.to_csv("labels_meta.csv",index=False)
pf=pq.ParquetFile(U+"data/processed/drug_gene_features.parquet")
print(pf.metadata.num_row_groups, pf.schema_arrow)
for sc in ["all","lung"]:
    parts=[]
    for i in range(pf.metadata.num_row_groups):
        t=pf.read_row_group(i,columns=["scope","drug","concentration","gene_name","mean_stat"])
        t=t.filter(pc.equal(t["scope"],sc)).drop(["scope"])
        parts.append(t.to_pandas().astype({"drug":"category","gene_name":"category","mean_stat":"float32","concentration":"float32"}))
    d=pd.concat(parts,ignore_index=True); del parts
    d["drug"]=d.drug.astype(str)
    d=d[d.drug.isin(set(s.drug_tahoe))]
    d["key"]=d.drug.map(dict(zip(s.drug_tahoe,s.key)))
    d["dose"]=(d.groupby("drug").concentration.rank(method="dense").astype(int)-1).astype("int8")
    genes=pd.Index(sorted(d.gene_name.astype(str).unique())); keys=sorted(d.key.unique())
    gi=genes.get_indexer(d.gene_name.astype(str)); ki=pd.Index(keys).get_indexer(d.key)
    X=np.zeros((len(keys),3,len(genes)),dtype=np.float32); N=np.zeros_like(X)
    np.add.at(X,(ki,d.dose.values,gi),d.mean_stat.values); np.add.at(N,(ki,d.dose.values,gi),1)
    X=np.divide(X,N,out=np.zeros_like(X),where=N>0)
    np.save(f"X_{sc}.npy",X); pd.Series(keys).to_csv(f"keys_{sc}.csv",index=False); pd.Series(genes).to_csv(f"genes_{sc}.csv",index=False)
    print(sc,X.shape, "doses per drug:", d.groupby("drug").dose.nunique().value_counts().to_dict()); del d
u=s.drop_duplicates("key")
c=u[u.dict_concern.isin(["most","no"])]
print("clear:",len(c),"oncology among clear:",int(c.oncology.sum()), c[c.oncology].dict_concern.value_counts().to_dict())
m=u[u.dict_concern.isin(["most","less","no"])]
print("most+less+no:", m.dict_concern.value_counts().to_dict(), "oncology:", m[m.oncology].dict_concern.value_counts().to_dict())

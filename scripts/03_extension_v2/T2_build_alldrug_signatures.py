"""v2 / Tier 2 helper: all-drug Tahoe signatures (mean Wald statistic over the 50 lines, per dose) from
data/processed/drug_gene_features.parquet (scope 'all'), as used for D1/D2/C1.  Output: X (drugs x 3 doses x genes)."""
import sys, importlib, numpy as np, pandas as pd, pyarrow.parquet as pq, pyarrow.compute as pc
sys.path.insert(0, "/home/claude/mirror/scripts"); norm = importlib.import_module("02_drug_overlap").norm
U = "/mnt/user-data/uploads/Protein and RNA/"; OUT = "/home/claude/v2/results/cache/"
pf = pq.ParquetFile(U + "data/processed/drug_gene_features.parquet")
# pass 1: drug x concentration and gene lists
dc, genes = set(), set()
for i in range(pf.metadata.num_row_groups):
    t = pf.read_row_group(i, columns=["scope", "drug", "concentration", "gene_name"]); t = t.filter(pc.equal(t["scope"], "all"))
    d = t.drop(["scope"]).to_pandas(); dc |= set(zip(d.drug, d.concentration)); genes |= set(d.gene_name.astype(str))
dcf = pd.DataFrame(sorted(dc), columns=["drug", "concentration"])
dcf["dose"] = dcf.groupby("drug").concentration.rank(method="dense").astype(int) - 1
dcf["key"] = dcf.drug.map(norm)
keys = sorted(dcf.key.unique()); genes = sorted(genes)
ki = {k: i for i, k in enumerate(keys)}; gi = {g: i for i, g in enumerate(genes)}
dose = {(r.drug, r.concentration): (ki[r.key], r.dose) for r in dcf.itertuples()}
X = np.zeros((len(keys), 3, len(genes)), np.float32); N = np.zeros_like(X)
for i in range(pf.metadata.num_row_groups):
    t = pf.read_row_group(i, columns=["scope", "drug", "concentration", "gene_name", "mean_stat"]); t = t.filter(pc.equal(t["scope"], "all"))
    d = t.drop(["scope"]).to_pandas()
    if not len(d): continue
    kd = np.array([dose[(a, b)] for a, b in zip(d.drug, d.concentration)])
    g = d.gene_name.astype(str).map(gi).values
    np.add.at(X, (kd[:, 0], np.clip(kd[:, 1], 0, 2), g), d.mean_stat.values.astype(np.float32))
    np.add.at(N, (kd[:, 0], np.clip(kd[:, 1], 0, 2), g), 1)
X = np.divide(X, N, out=np.zeros_like(X), where=N > 0)
np.savez_compressed(OUT + "tahoe_alldrug_signatures.npz", X=X, keys=np.array(keys, dtype=object), genes=np.array(genes, dtype=object))
print(X.shape, "doses per drug:", dcf.groupby("key").dose.nunique().value_counts().to_dict())

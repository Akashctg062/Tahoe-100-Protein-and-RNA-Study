"""
v2 / Tier 1-D + Tier 2-E data step - build Tahoe per-line log2 fold-change tables (top AND middle dose) for the
80 drugs in T2_drug_keys.txt (PRISM/GDSC drugs of C3 + CTRP drugs), in ONE pass over the Tahoe files.
Same rules as the frozen C3 (09_confirmatory_C2_C3.build_response): gene universe = 5,000 genes with the highest
mean baseMean in the first 20 files + the 20 proliferation genes; dose rank per Tahoe drug name; salt forms averaged.

Run on your computer (PowerShell, in the window where 'python' worked for the ChEMBL script; ~15-25 min):
    cd "E:\\Protein and RNA\\v2_extended\\scripts"
    python T2_build_responses.py
Reads (read-only):  E:\\Protein and RNA\\data\\processed\\pseudobulk_parts_all\\*.parquet
Writes:             E:\\Protein and RNA\\v2_extended\\results\\responses\\response_top_80drugs.parquet
                    E:\\Protein and RNA\\v2_extended\\results\\responses\\response_mid_80drugs.parquet
"""
import os, sys, glob, time, importlib
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
V1 = os.path.abspath(os.path.join(HERE, "..", "..", "scripts"))
PARTS = os.path.abspath(os.path.join(HERE, "..", "..", "data", "processed", "pseudobulk_parts_all"))
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "responses")); os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, V1)
norm = importlib.import_module("02_drug_overlap").norm
PROLIF = ["MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "CCNA2", "BUB1", "PLK1", "AURKA", "AURKB",
          "CENPF", "E2F1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "TYMS"]

def main():
    t0 = time.time()
    keys = {l.strip() for l in open(os.path.join(HERE, "T2_drug_keys.txt"), encoding="utf-8") if l.strip()}
    parts = sorted(glob.glob(os.path.join(PARTS, "*.parquet")))
    print(f"{len(keys)} drug keys | {len(parts)} Tahoe files", flush=True)
    bm = pd.concat([pd.read_parquet(p, columns=["gene_name", "baseMean"]) for p in parts[:20]])
    universe = set(bm.groupby("gene_name").baseMean.mean().sort_values(ascending=False).index[:5000]) | set(PROLIF)
    print(f"gene universe: {len(universe)} genes", flush=True)
    acc = []
    for i, p in enumerate(parts, 1):
        d = pd.read_parquet(p, columns=["drug", "concentration", "Cell_Name_Vevo", "gene_name", "log2FoldChange"])
        d = d[d.gene_name.isin(universe)]
        d = d[d.drug.map(norm).isin(keys)]
        if len(d):
            for c in ["drug", "Cell_Name_Vevo", "gene_name"]: d[c] = d[c].astype("category")
            d["log2FoldChange"] = d.log2FoldChange.astype("float32"); acc.append(d)
        if i % 100 == 0 or i == len(parts): print(f"  {i}/{len(parts)} files scanned ({(time.time()-t0)/60:.1f} min)", flush=True)
    d = pd.concat(acc, ignore_index=True)
    for c in ["drug", "Cell_Name_Vevo", "gene_name"]: d[c] = d[c].astype(str)
    d["dose_rank"] = d.groupby("drug").concentration.rank(method="dense").astype(int)
    d["max_rank"] = d.groupby("drug").dose_rank.transform("max")
    d["drug_key"] = d.drug.map(norm)
    for name, sel in [("top", d.dose_rank == d.max_rank), ("mid", d.dose_rank == 2)]:
        x = d[sel].groupby(["drug_key", "Cell_Name_Vevo", "gene_name"], as_index=False).log2FoldChange.mean()
        f = os.path.join(OUT, f"response_{name}_80drugs.parquet"); x.to_parquet(f, index=False)
        print(f"{name} dose: {x.drug_key.nunique()} drugs, {x.Cell_Name_Vevo.nunique()} lines, {len(x):,} rows -> {f}", flush=True)
    print(f"Done in {(time.time()-t0)/60:.1f} min")

if __name__ == "__main__":
    main()

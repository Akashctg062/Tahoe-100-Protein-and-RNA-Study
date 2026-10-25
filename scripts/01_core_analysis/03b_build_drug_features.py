"""
Step 3b - Summarise the ~1,026 filtered Tahoe files into compact drug-level features,
reading ONE file at a time (needs only ~1-2 GB RAM, takes ~5-15 minutes).

For every drug x dose x gene it computes the mean differential-expression statistic over
  scope = "all"  : all cell lines
  scope = "lung" : lung cell lines only (A549, NCI-H460, ... from cell_line_metadata 'Organ')
Output: ../data/processed/drug_gene_features.parquet
        columns: scope, drug, concentration, gene_name, mean_stat, mean_lfc, n_lines

Run:  python 03b_build_drug_features.py
"""
import argparse, glob, os
import pandas as pd

def main(a):
    parts = sorted(glob.glob(os.path.join(a.parts, "*.parquet")))
    print(f"{len(parts)} filtered files found")
    cl = pd.read_csv(a.cells)
    lung = set(cl.loc[cl["Organ"].astype(str).str.contains("Lung", case=False), "cell_name"])
    print(f"{len(lung)} lung cell lines")
    keys = ["drug", "concentration", "gene_name"]
    acc = []
    for i, p in enumerate(parts, 1):
        d = pd.read_parquet(p, columns=keys + ["Cell_Name_Vevo", "stat", "log2FoldChange"])
        d = d.dropna(subset=["stat"])
        for scope, sub in (("all", d), ("lung", d[d.Cell_Name_Vevo.isin(lung)])):
            if len(sub) == 0: continue
            g = sub.groupby(keys, observed=True).agg(sum_stat=("stat", "sum"),
                                                     sum_lfc=("log2FoldChange", "sum"),
                                                     n=("stat", "size")).reset_index()
            g["scope"] = scope; acc.append(g)
        if i % 50 == 0 or i == len(parts):
            # compact the accumulator periodically to keep memory low
            acc = [pd.concat(acc).groupby(["scope"] + keys, observed=True)
                     [["sum_stat", "sum_lfc", "n"]].sum().reset_index()]
            print(f"  {i}/{len(parts)} files summarised")
    f = acc[0]
    f["mean_stat"] = (f.sum_stat / f.n).astype("float32")
    f["mean_lfc"] = (f.sum_lfc / f.n).astype("float32")
    f = f.rename(columns={"n": "n_lines"})[["scope"] + keys + ["mean_stat", "mean_lfc", "n_lines"]]
    f.to_parquet(a.out, index=False, compression="zstd")
    print(f"\nSaved {a.out}: {len(f):,} rows")
    for sc in ("all", "lung"):
        s = f[f.scope == sc]
        print(f"  scope={sc}: {s.drug.nunique()} drugs, {s.groupby('drug').concentration.nunique().median():.0f} doses/drug, "
              f"{s.gene_name.nunique()} genes, median lines per drug-dose = {s.n_lines.median():.0f}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="../data/processed/pseudobulk_parts")
    ap.add_argument("--cells", default="../data/raw/tahoe_metadata/cell_line_metadata.csv")
    ap.add_argument("--out", default="../data/processed/drug_gene_features.parquet")
    main(ap.parse_args())

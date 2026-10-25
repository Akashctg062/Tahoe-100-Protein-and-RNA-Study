"""
v2 / Tier 1 helper - build a cache of the C3 per-drug datasets (identical inputs to the frozen C3 / Amendment 1),
so that T1-A2, T1-C and Tier 2 analyses use exactly the same drugs, lines and features.
Run from the scripts folder that contains 09_confirmatory_C2_C3.py (paths ../data, ../results as in v1).
Output: <out>/c3_cache.pkl  -> {screen: [ {drug, cells, y, lung, X: {Response, RNA, Protein_MS, RPPA}} ]}
"""
import os, sys, pickle, importlib, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def main(a):
    sys.path.insert(0, a.v1_scripts); os.chdir(a.v1_scripts)
    m09 = importlib.import_module("09_confirmatory_C2_C3")
    R = m09.R
    cl = pd.read_csv(R + "tahoe_metadata/cell_line_metadata.csv").drop_duplicates("cell_name")
    gm = pd.read_csv(R + "tahoe_metadata/gene_metadata.csv")
    lrs = pd.read_csv("../results/response_vs_sensitivity/line_response_strength.csv", usecols=["Cell_Name_Vevo", "drug_key"])
    lines = cl[cl.cell_name.isin(set(lrs.Cell_Name_Vevo))].copy(); lines["cell_key"] = lines.cell_name.map(m09.norm_cell)
    lines["lung"] = lines.Organ.astype(str).str.contains("Lung", case=False)
    is_lung = lines.set_index("cell_name").lung
    labels = m09.load_labels(lines)
    resp = pd.read_parquet(a.response)
    rna, ms, rp = m09.load_baselines(lines, gm)
    out = {}
    for src, L in labels.items():
        out[src] = []
        for dk in sorted(set(L.drug_key) & set(resp.drug_key)):
            Ld = L[L.drug_key == dk].set_index("cell_name").auc
            R_ = resp[resp.drug_key == dk].pivot_table(index="Cell_Name_Vevo", columns="gene_name", values="log2FoldChange")
            common = sorted(set(Ld.index) & set(R_.index) & set(rna.index) & set(ms.index) & set(rp.index))
            if len(common) < a.min_lines: continue
            g_shared = [g for g in R_.columns if g in rna.columns]
            out[src].append(dict(drug=dk, cells=common, y=Ld.loc[common].values.astype(float),
                                 lung=is_lung.reindex(common).fillna(False).values.astype(bool),
                                 X={"Response": R_.loc[common].fillna(0).values, "RNA": rna.loc[common, g_shared].values,
                                    "Protein_MS": ms.loc[common].values, "RPPA": rp.loc[common].values},
                                 rppa_names=list(rp.columns), genes=g_shared))
        print(src, len(out[src]), "drugs")
    os.makedirs(a.out, exist_ok=True)
    pickle.dump(out, open(os.path.join(a.out, "c3_cache.pkl"), "wb"))
    # also keep the full baselines for Tier 2 (new screens / proteome)
    pickle.dump(dict(rna=rna, ms=ms, rppa=rp, lines=lines), open(os.path.join(a.out, "baselines.pkl"), "wb"))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--v1-scripts", default="/home/claude/mirror/scripts")
    ap.add_argument("--response", default="/home/claude/mirror/results/confirmatory/response_dosetop.parquet")
    ap.add_argument("--out", default="/home/claude/v2/results/cache")
    ap.add_argument("--min-lines", type=int, default=20)
    main(ap.parse_args())

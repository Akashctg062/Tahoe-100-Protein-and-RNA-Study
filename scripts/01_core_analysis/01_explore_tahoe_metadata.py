"""
Step 1 - Download the small Tahoe-100M metadata tables (~2.5 MB in total, not the 429 GB atlas)
with plain HTTPS (no account, no token, no Hugging Face library) and save them as CSV.

Run:  python 01_explore_tahoe_metadata.py --out ../data/raw/tahoe_metadata
"""
import argparse, io, os
import pandas as pd, requests

BASE = "https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/{t}.parquet"
TABLES = ["drug_metadata", "cell_line_metadata", "sample_metadata", "gene_metadata"]

def main(out):
    os.makedirs(out, exist_ok=True)
    for t in TABLES:
        r = requests.get(BASE.format(t=t), timeout=120); r.raise_for_status()
        df = pd.read_parquet(io.BytesIO(r.content))
        df.to_csv(os.path.join(out, f"{t}.csv"), index=False)
        print(f"{t:22s} rows={len(df):>7,}  cols={list(df.columns)[:8]}")
    cl = pd.read_csv(os.path.join(out, "cell_line_metadata.csv"))
    lung = cl[cl["Organ"].astype(str).str.contains("Lung", case=False)]
    print("\nLung cell lines in Tahoe-100M:", sorted(lung["cell_name"].unique()))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="../data/raw/tahoe_metadata")
    main(ap.parse_args().out)

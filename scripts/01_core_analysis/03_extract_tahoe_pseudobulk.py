"""
Step 3 - Extract ONLY the Tahoe-100M pseudobulk differential-expression rows you need.

NO account, NO token, NO Hugging Face library needed.
The public parquet files are downloaded with plain HTTPS (like a browser download),
filtered on your computer, and the raw file is deleted immediately afterwards.

What it does, one file at a time (1,026 files, ~90 MB each, ~90 GB in total):
  1. download file -> ../data/processed/_tmp_download/
  2. keep only drugs in ../data/processed/shared_drugs.csv and expressed genes (baseMean >= 1)
  3. save the small filtered file to ../data/processed/pseudobulk_parts/
  4. delete the raw download (disk use stays ~100 MB)
It is RESUMABLE: rerun the same command after any interruption; finished files are skipped.
When all files are done, run 03b_build_drug_features.py (memory-light summary; no 9 GB merge).

Requirements (already installed with pandas):  pip install pandas pyarrow requests

All 379 drugs (lung drug-response study), saved to a separate folder:
  python 03_extract_tahoe_pseudobulk.py --all-drugs --parts ../data/processed/pseudobulk_parts_all

Run (copy-paste), first a 2-file speed test:
  python 03_extract_tahoe_pseudobulk.py --max-shards 2
then the full run (can take several hours; leave it running):
  python 03_extract_tahoe_pseudobulk.py
"""
import argparse, glob, os, time
import pandas as pd, pyarrow as pa, pyarrow.parquet as pq, requests

N_SHARDS = 1026
BASE = ("https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/"
        "metadata/pseudobulk_differential_expression/train-{i:05d}-of-01026.parquet")
COLS = ["drug", "concentration", "Cell_Name_Vevo", "Cell_ID_DepMap", "gene_name", "baseMean",
        "log2FoldChange", "stat", "padj", "n_cells_trt", "n_cells_ctrl"]

def download(url, dest, tries=5):
    for k in range(tries):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(dest + ".part", "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            os.replace(dest + ".part", dest); return True
        except Exception as e:
            print(f"    download attempt {k+1} failed: {e}"); time.sleep(15 * (k + 1))
    return False

def main(a):
    if a.all_drugs:
        drugs = None; print("Keeping ALL Tahoe drugs")
    else:
        drugs = set(pd.read_csv(a.drugs)["drug_tahoe"].dropna())
        print(f"{len(drugs)} Tahoe drug names to keep")
    tmpdir = os.path.join(os.path.dirname(a.parts), "_tmp_download")
    os.makedirs(a.parts, exist_ok=True); os.makedirs(tmpdir, exist_ok=True)

    todo = range(N_SHARDS if not a.max_shards else a.max_shards)
    t0, done_now = time.time(), 0
    for i in todo:
        name = f"train-{i:05d}-of-01026.parquet"
        out = os.path.join(a.parts, name)
        if os.path.exists(out):
            continue
        raw = os.path.join(tmpdir, name)
        if not download(BASE.format(i=i), raw):
            print(f"  skipping {name} for now - rerun later"); continue
        tbl = pq.read_table(raw, columns=COLS,
                            filters=([("drug", "in", list(drugs))] if drugs else []) + [("baseMean", ">=", a.min_basemean)])
        # store floats as float32 to halve disk use
        tbl = tbl.cast(pa.schema([(f.name, pa.float32() if pa.types.is_floating(f.type) else f.type)
                                  for f in tbl.schema]))
        pq.write_table(tbl, out + ".tmp", compression="zstd"); os.replace(out + ".tmp", out)
        os.remove(raw); done_now += 1
        el = (time.time() - t0) / 60
        print(f"[{i+1}/{N_SHARDS}] {name}: kept {tbl.num_rows:,} rows | {el:.1f} min "
              f"| ~{el/done_now:.2f} min per file")

    parts = sorted(glob.glob(os.path.join(a.parts, "*.parquet")))
    print(f"\n{len(parts)} of {N_SHARDS} files finished")
    if len(parts) < N_SHARDS:
        if parts:
            df = pd.read_parquet(parts[0])
            print(f"Example file: {len(df):,} rows, {df.drug.nunique()} drugs, "
                  f"{df.Cell_Name_Vevo.nunique()} cell lines, {df.gene_name.nunique()} genes")
        print("Not finished yet -> run the script again (without --max-shards) to continue.")
        return
    print("ALL FILES DONE. Next step (summarise into drug-level features, low memory):")
    print("  python 03b_build_drug_features.py")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--drugs", default="../data/processed/shared_drugs.csv")
    ap.add_argument("--parts", default="../data/processed/pseudobulk_parts")
    ap.add_argument("--min-basemean", dest="min_basemean", type=float, default=1.0)
    ap.add_argument("--max-shards", dest="max_shards", type=int, default=0)
    ap.add_argument("--all-drugs", dest="all_drugs", action="store_true",
                    help="keep all 379 Tahoe drugs (use with --parts ../data/processed/pseudobulk_parts_all)")
    main(ap.parse_args())

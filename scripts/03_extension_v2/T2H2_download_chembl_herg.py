"""
T2-H2 data step - download human hERG (ChEMBL target CHEMBL240) IC50 / Ki activity records from the
public ChEMBL REST API.  Pre-specified in Amendment 3 (18 Oct 2026).

Run (PowerShell):
    conda activate oncocardio
    cd "E:\\Protein and RNA\\v2_extended\\scripts"
    python T2H2_download_chembl_herg.py
Output: E:\\Protein and RNA\\data\\raw\\v2\\chembl_herg\\chembl_CHEMBL240_IC50_Ki.csv  (+ a small log)
Needs only 'requests' and 'pandas'. Takes about 2-10 minutes. Safe to re-run (resumes from saved pages).
"""
import os, sys, time, json
import requests, pandas as pd

OUT = os.path.join("..", "..", "data", "raw", "v2", "chembl_herg")
os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(OUT, "pages"); os.makedirs(CACHE, exist_ok=True)
URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
FIELDS = ["activity_id", "molecule_chembl_id", "canonical_smiles", "molecule_pref_name", "standard_type",
          "standard_relation", "standard_value", "standard_units", "pchembl_value", "assay_chembl_id",
          "assay_type", "assay_description", "target_organism", "data_validity_comment", "document_year"]

def get(params, tries=6):
    for k in range(tries):
        try:
            r = requests.get(URL, params=params, timeout=120)
            if r.status_code == 200: return r.json()
            print(f"  HTTP {r.status_code}; retrying", flush=True)
        except requests.RequestException as e:
            print(f"  {type(e).__name__}; retrying", flush=True)
        time.sleep(10 * (k + 1))
    sys.exit("ChEMBL API not reachable - please send me this message.")

rows, offset, limit = [], 0, 1000
for stype in ["IC50", "Ki"]:
    offset = 0
    while True:
        f = os.path.join(CACHE, f"{stype}_{offset}.json")
        if os.path.exists(f):
            js = json.load(open(f, encoding="utf-8"))
        else:
            js = get({"target_chembl_id": "CHEMBL240", "standard_type": stype, "limit": limit, "offset": offset})
            json.dump(js, open(f, "w", encoding="utf-8"))
        acts = js.get("activities", [])
        rows += [{k: a.get(k) for k in FIELDS} for a in acts]
        total = js.get("page_meta", {}).get("total_count", 0)
        print(f"{stype}: {min(offset + limit, total)}/{total}", flush=True)
        offset += limit
        if not js.get("page_meta", {}).get("next"): break

df = pd.DataFrame(rows).drop_duplicates("activity_id")
dst = os.path.join(OUT, "chembl_CHEMBL240_IC50_Ki.csv")
df.to_csv(dst, index=False)
msg = (f"Saved {len(df):,} records ({df.molecule_chembl_id.nunique():,} molecules) to {os.path.abspath(dst)}\n"
       f"ChEMBL download date: {time.strftime('%Y-%m-%d')}")
open(os.path.join(OUT, "download_log.txt"), "w").write(msg); print(msg)

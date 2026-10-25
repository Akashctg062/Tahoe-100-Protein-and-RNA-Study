"""
v2 / Tier 2 data step - Tahoe drug signatures (mean Wald statistic over all lines, 3 doses), identical to the
D1/C1 features, for the Tahoe drugs not yet covered (list in T2_missing_signature_keys.txt).
Uses the frozen v1 function tahoe_features() from 10_confirmatory_C1_faers.py (read-only).
Run (same window as before; ~15-25 min):
    cd "E:\\Protein and RNA\\v2_extended\\scripts"
    python T2_build_meanstat_signatures.py
Writes: E:\\Protein and RNA\\v2_extended\\results\\responses\\tahoe_signatures_remaining.npz
"""
import os, sys, importlib
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
V1 = os.path.abspath(os.path.join(HERE, "..", "..", "scripts"))
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "responses")); os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, V1); os.chdir(V1)
m02 = importlib.import_module("02_drug_overlap"); m10 = importlib.import_module("10_confirmatory_C1_faers")
want = {l.strip() for l in open(os.path.join(HERE, "T2_missing_signature_keys.txt"), encoding="utf-8") if l.strip()}
dm = pd.read_csv("../data/raw/tahoe_metadata/drug_metadata.csv").drop_duplicates("drug")
names = {n: m02.norm(n) for n in dm.drug if m02.norm(n) in want}
print(f"{len(want)} drug keys requested, {len(names)} Tahoe drug names", flush=True)
X, keys, genes = m10.tahoe_features("../data/processed/pseudobulk_parts_all", names, os.path.join(OUT, "tahoe_signatures_remaining.npz"))
print("Done:", X.shape, "->", os.path.join(OUT, "tahoe_signatures_remaining.npz"))

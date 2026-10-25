"""
Step 2 - FEASIBILITY GATE: which Tahoe-100M drugs have cardiotoxicity labels?

Primary label source (dry-lab study): FDA DICTrank (1,318 drugs, Excel file from
  https://www.fda.gov/science-research/bioinformatics-tools/drug-induced-cardiotoxicity-rank-dictrank-dataset)
Optional: any CSV with a column 'drug' (e.g., predicTox drug list) via --cardiac

Run (copy-paste):
  python 02_drug_overlap.py --dictrank ../data/raw/dictrank_dataset_508.xlsx
  python 02_drug_overlap.py --cardiac ../data/raw/predictox_drugs.csv --out ../data/processed/shared_predictox.csv
Decision rule: >= 60 Tahoe drugs with a clear DICTrank label (Most- or No-concern) -> proceed.
"""
import argparse, re
import pandas as pd

# salt / hydrate / formulation words removed from BOTH sides before matching
SALT_WORDS = set("""hydrochloride monohydrochloride dihydrochloride hcl mesylate dimesylate maleate
ditosylate tosylate sulfate disulfate bisulfate citrate malate s-malate phosphate disodium sodium
potassium calcium magnesium acetate diacetate fumarate hemifumarate bitartrate tartrate anhydrous
hydrate monohydrate dihydrate trihydrate pentahydrate hexahydrate sesquihydrate hemihydrate besylate
succinate bromide hydrobromide olamine isethionate erbumine pamoate lactate gluconate nitrate
salt free base liposome liposomal injection glucuronate lactobionate stearate
ethylsuccinate estolate""".split())

# Tahoe name -> DICTrank (USAN/INN) name, for common non-matching synonyms
SYNONYMS = {"5azacytidine": "azacitidine", "5fluorouracil": "fluorouracil", "rapamycin": "sirolimus",
            "cyclosporina": "cyclosporine", "retinoicacid": "tretinoin", "lthyroxine": "levothyroxine",
            "aspirin": "acetylsalicylicacid", "leflornithine": "eflornithine", "vitamink4": "menadiol",
            "sacubitrilvalsartan": "sacubitril", "quinidine15dihydroquinidine": "quinidine",
            "torkinib": "pp242"}

def norm(name: str) -> str:
    s = str(name).lower()
    while re.search(r"\([^()]*\)", s):                  # drop "(hydrochloride)", nested too
        s = re.sub(r"\([^()]*\)", " ", s)
    allw = [w for w in re.split(r"[\s,]+", s) if w]
    words = [w for w in allw if w not in SALT_WORDS] or allw   # all-salt names (e.g. calcium lactate) keep every word
    key = re.sub(r"[^a-z0-9]", "", "".join(words))
    return SYNONYMS.get(key, key)

def pick(cols, prefs, exclude=()):
    for p in prefs:
        for c in cols:
            if p in str(c).lower() and c not in exclude:
                return c
    return None

def load_dictrank(path, name_col=None, label_col=None):
    raw = pd.read_excel(path, header=None)
    # find the header row (FDA files sometimes have title rows above the table)
    hdr = 0
    for i in range(min(15, len(raw))):
        cells = [str(v).lower() for v in raw.iloc[i].values if pd.notna(v)]
        has_name = any("name" in c for c in cells)
        has_label = any(("concern" in c or "dict" in c or "category" in c) and "name" not in c for c in cells)
        if has_name and has_label and len(cells) >= 2:
            hdr = i; break
    d = pd.read_excel(path, header=hdr)
    print("DICTrank columns found:", list(d.columns))
    name_col = name_col or pick(d.columns, ["generic", "proper", "drug", "name"])
    label_col = label_col or pick(d.columns, ["concern", "dictrank", "dict", "category", "rank"], exclude=(name_col,))
    if name_col is None or label_col is None:
        raise SystemExit("Could not detect columns - rerun with --name-col \"...\" --label-col \"...\" using names printed above")
    print(f"Using name column = '{name_col}', label column = '{label_col}'")
    d = d.rename(columns={label_col: "dict_concern"})
    d["dict_concern"] = d.dict_concern.astype(str).str.strip().str.lower()
    ai_col = pick(d.columns, ["active ingredient"], exclude=(name_col,))
    parts = [d[[name_col, "dict_concern"]].rename(columns={name_col: "drug_cardiac"})]
    if ai_col is not None:                               # also match on active ingredient
        parts.append(d[[ai_col, "dict_concern"]].rename(columns={ai_col: "drug_cardiac"}))
    d = pd.concat(parts).dropna(subset=["drug_cardiac"])
    # some rows list several names separated by ';' or ',' -> one row per name
    d = d.assign(drug_cardiac=d.drug_cardiac.astype(str).str.split(r"[;,]")).explode("drug_cardiac")
    d["drug_cardiac"] = d.drug_cardiac.str.strip()
    return d[["drug_cardiac", "dict_concern"]].drop_duplicates()

def main(a):
    t = pd.read_csv(a.tahoe)
    tcol = "drug" if "drug" in t.columns else t.columns[0]
    t = t.rename(columns={tcol: "drug_tahoe"}); t["key"] = t.drug_tahoe.map(norm)
    t = t[t.key.str.len() > 0]
    if a.dictrank:
        c = load_dictrank(a.dictrank, a.name_col, a.label_col)
    else:
        c = pd.read_csv(a.cardiac).rename(columns={"drug": "drug_cardiac"})
    c["key"] = c.drug_cardiac.map(norm); c = c[c.key.str.len() > 0]
    if "dict_concern" in c:
        # a compound can appear several times in DICTrank (different products/formulations).
        # Rule (state it in Methods): keep the HIGHEST concern; record conflicts.
        sev = {"most": 3, "less": 2, "ambiguous": 1, "no": 0}
        c["sev"] = c.dict_concern.map(sev).fillna(-1)
        g = c.groupby("key")
        lab = g.agg(drug_cardiac=("drug_cardiac", "first"), sev=("sev", "max"),
                    n_labels=("dict_concern", "nunique")).reset_index()
        inv = {v: k for k, v in sev.items()}
        lab["dict_concern"] = lab.sev.map(inv); lab["label_conflict"] = lab.n_labels > 1
        lab = lab.drop(columns=["sev", "n_labels"])
    else:
        lab = c[["key", "drug_cardiac"]].drop_duplicates("key")
    # keep EVERY Tahoe name (salt variants are separate Tahoe conditions); 'key' = compound ID
    shared = t[["key", "drug_tahoe"]].drop_duplicates().merge(lab, on="key")
    print(f"Tahoe compounds: {t.key.nunique()} | DICTrank compounds: {c.key.nunique()} | "
          f"SHARED compounds: {shared.key.nunique()} (Tahoe names: {len(shared)})")
    if "dict_concern" in shared:
        per = shared.drop_duplicates("key")
        print(per.dict_concern.value_counts().to_string())
        print("compounds with conflicting DICTrank labels:", int(per.label_conflict.sum()))
        clear = per.dict_concern.isin(["most", "no"]).sum()
        print(f"Clear labels (most/no): {clear}")
        print("DECISION:", "proceed" if clear >= 60 else "too few clear labels - add predicTox/FAERS labels")
    approved = t[t.get("human-approved", pd.Series(dtype=str)).astype(str).str.lower() == "yes"] if "human-approved" in t else None
    if approved is not None and len(approved):
        miss = approved[~approved.key.isin(shared.key)].drug_tahoe.unique()
        print(f"\nFDA-approved Tahoe drugs still unmatched ({len(miss)}): check manually ->")
        print(", ".join(sorted(miss)))
    shared.to_csv(a.out, index=False); print("saved ->", a.out)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tahoe", default="../data/raw/tahoe_metadata/drug_metadata.csv")
    ap.add_argument("--dictrank"); ap.add_argument("--cardiac")
    ap.add_argument("--name-col", dest="name_col"); ap.add_argument("--label-col", dest="label_col")
    ap.add_argument("--out", default="../data/processed/shared_drugs.csv")
    a = ap.parse_args()
    assert a.dictrank or a.cardiac, "give --dictrank or --cardiac"
    main(a)

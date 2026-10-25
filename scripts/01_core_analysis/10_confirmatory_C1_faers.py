"""
Step 10 - CONFIRMATORY analysis C1 (frozen plan, 15 Oct 2026):
Does adding Tahoe tumour-cell responses to chemistry improve prediction of an INDEPENDENT
cardiotoxicity label (FDA FAERS reporting-odds-ratio signal)?

Frozen specification
  Label      : FAERS cardiac signal = lower 95% CI of ROR > 1 AND >= 3 cardiac reports
               (MedDRA preferred terms fixed below; openFDA, all report dates)
  Drugs      : FDA-approved Tahoe compounds NOT among the 100 DICTrank compounds used in D1
  Features   : chemistry (RDKit descriptors + Morgan); Tahoe PCA-10 of the 3-dose mean DE statistic;
               chemistry + Tahoe
  Model / CV : logistic regression, repeated stratified 5-fold CV (10 repeats)
  Primary    : AUROC(chemistry + Tahoe) - AUROC(chemistry); 95% CI from 2,000 bootstrap resamples of drugs
  Decision   : Tahoe adds information only if the lower 95% CI > 0
Implementation details fixed in this script (not specified in the plan; reported in Methods):
  * a drug enters C1 only if it has >= 100 FAERS reports in total (ROR otherwise unstable)
  * drug names are searched in patient.drug.openfda.generic_name OR
    patient.drug.activesubstance.activesubstancename (tokenised, so "SUNITINIB MALATE" matches
    "sunitinib"); salt words removed with the step-2 rules. Fixed on 16 Oct 2026 BEFORE any C1 result
    was computed, because v1 (generic name only) gave 0 reports for approved drugs such as abemaciclib.
  * out-of-fold probabilities are averaged over the 10 CV repeats before bootstrapping

First test the connection (10 seconds):  python 10_confirmatory_C1_faers.py --check-only --api-key YOUR_KEY
Run (PowerShell; a free openFDA key is recommended: https://open.fda.gov/apis/authentication/):
  $env:OPENFDA_KEY="your_key"          (optional)
  python 10_confirmatory_C1_faers.py
It is resumable: FAERS counts are cached in ../results/confirmatory/C1/faers_counts_v2.csv
(v2 = generic name OR active substance; v1 file from the first run is kept for the record)
Outputs: ../results/confirmatory/C1/
"""
import argparse, glob, math, os, re, sys, time, importlib, warnings
import numpy as np, pandas as pd, requests
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m02 = importlib.import_module("02_drug_overlap")

API = "https://api.fda.gov/drug/event.json"
CARDIAC_PTS = ["Cardiac failure", "Cardiac failure congestive", "Left ventricular dysfunction",
               "Ejection fraction decreased", "Cardiomyopathy", "Cardiotoxicity", "Myocarditis",
               "Electrocardiogram QT prolonged", "Torsade de pointes", "Ventricular tachycardia",
               "Atrial fibrillation", "Myocardial infarction"]
CARD = "patient.reaction.reactionmeddrapt.exact:(" + " ".join(f'"{p}"' for p in CARDIAC_PTS) + ")"
ALL = "_exists_:safetyreportid"

# ------------------------------------------------------------------ openFDA
def make_total(key):
    sess = requests.Session()
    def total(q):
        params = {"search": q, "limit": 1}
        if key: params["api_key"] = key
        for attempt in range(8):
            try:
                r = sess.get(API, params=params, timeout=60)
                if r.status_code == 404: return 0                       # no matching reports
                if r.status_code == 429 or r.status_code >= 500:         # rate limit / server busy
                    time.sleep(5 * (attempt + 1)); continue
                r.raise_for_status(); time.sleep(0.25 if key else 0.6)
                return int(r.json()["meta"]["results"]["total"])
            except requests.RequestException as e:
                print(f"    network problem ({e}); retrying"); time.sleep(5 * (attempt + 1))
        raise RuntimeError(f"openFDA failed repeatedly for query: {q}")
    return total

def search_term(drug_tahoe):
    """Tahoe name -> plain generic name for openFDA (salt words and brackets removed)."""
    raw = {"lactate (calcium)": "calcium lactate", "l-eflornithine (monohydrochloride)": "eflornithine"}
    if str(drug_tahoe).strip().lower() in raw: return raw[str(drug_tahoe).strip().lower()]
    s = str(drug_tahoe).lower()
    while re.search(r"\([^()]*\)", s): s = re.sub(r"\([^()]*\)", " ", s)
    s = s.split("/")[0]                                   # combinations: search the first component
    words = [w for w in re.split(r"[\s,]+", s) if w and w not in m02.SALT_WORDS]
    phrase = " ".join(words).strip()
    syn = {"5-azacytidine": "azacitidine", "5-fluorouracil": "fluorouracil", "rapamycin": "sirolimus",
           "cyclosporin a": "cyclosporine", "l-thyroxine": "levothyroxine", "retinoic acid": "tretinoin"}
    phrase = syn.get(phrase, phrase)
    return phrase if phrase else None                    # None -> drug is skipped and reported

def drug_query(term):
    """Match the harmonised generic name OR the reported active substance (recent drugs often lack
    openFDA harmonisation, e.g. abemaciclib, darolutamide, belzutifan)."""
    return (f'(patient.drug.openfda.generic_name:"{term}" OR '
            f'patient.drug.activesubstance.activesubstancename:"{term}")')

def faers_counts(drugs, total, cache):
    done = pd.read_csv(cache) if os.path.exists(cache) else pd.DataFrame(columns=["key"])
    have = set(done.key)
    N = total(ALL); C = total(CARD)
    print(f"FAERS totals: {N:,} reports, {C:,} with a cardiac term")
    rows = []
    for i, (key, name) in enumerate(drugs.items(), 1):
        if key in have: continue
        term = search_term(name)
        if not term:
            print(f"  [{i}/{len(drugs)}] {name}: no usable generic name - skipped"); continue
        d = drug_query(term)
        n = total(d); a = total(f"{d} AND {CARD}") if n > 0 else 0
        rows.append(dict(key=key, drug_tahoe=name, search_term=term, n_drug=n, a=a, N=N, C=C))
        print(f"  [{i}/{len(drugs)}] {name} -> '{term}': {n:,} reports, {a:,} cardiac")
        if len(rows) % 10 == 0:
            done = pd.concat([done, pd.DataFrame(rows)]); done.to_csv(cache, index=False); rows = []
    done = pd.concat([done, pd.DataFrame(rows)]); done.to_csv(cache, index=False)
    return done

def add_ror(f):
    a = f.a.astype(float); b = f.n_drug - a; c = f.C - a; d = f.N - f.n_drug - c
    ok = (a > 0) & (b > 0) & (c > 0) & (d > 0)
    ror = np.where(ok, (a / b.where(b > 0, 1)) / (c / d.where(d > 0, 1)), np.nan)
    se = np.where(ok, np.sqrt(1 / a.where(a > 0, 1) + 1 / b.where(b > 0, 1) + 1 / c + 1 / d), np.nan)
    f["ROR"] = ror; f["ROR_lo95"] = np.exp(np.log(ror) - 1.96 * se); f["ROR_hi95"] = np.exp(np.log(ror) + 1.96 * se)
    f["signal"] = ((f.ROR_lo95 > 1) & (f.a >= 3)).astype(int)
    return f

# ------------------------------------------------------------------ Tahoe features (3 doses, mean over all lines)
def tahoe_features(parts_dir, keys_by_name, cache, dose_map_csv="../results/response_vs_sensitivity/line_response_strength.csv"):
    """Mean DE statistic over all cell lines, per drug x dose (3) x gene. Memory-light: each file is
    added straight into fixed numpy arrays (no large intermediate tables)."""
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True); return z["X"], list(z["keys"]), list(z["genes"])
    # dose rank per (drug, concentration), from step 8 (covers all 379 drugs)
    dm = pd.read_csv(dose_map_csv, usecols=["drug", "concentration", "dose_rank"]).drop_duplicates(["drug", "concentration"])
    dm["conc_r"] = dm.concentration.astype(float).round(6)
    dose_index = pd.MultiIndex.from_arrays([dm.drug.astype(str), dm.conc_r]); dose_vals = (dm.dose_rank.astype(int) - 1).to_numpy()
    keys = sorted(set(keys_by_name.values())); kidx = {k: i for i, k in enumerate(keys)}
    genes = sorted(pd.read_csv("../data/raw/tahoe_metadata/gene_metadata.csv").gene_symbol.dropna().astype(str).unique())
    gidx = {g: i for i, g in enumerate(genes)}
    S = np.zeros((len(keys), 3, len(genes)), np.float64); N = np.zeros((len(keys), 3, len(genes)), np.float32)
    parts = sorted(glob.glob(os.path.join(parts_dir, "*.parquet")))
    unknown_dose = 0
    for i, p in enumerate(parts, 1):
        d = pd.read_parquet(p, columns=["drug", "concentration", "gene_name", "stat"])
        d = d[d.drug.isin(keys_by_name.keys())].dropna(subset=["stat"])
        if len(d):
            ki = d.drug.astype(str).map({n: kidx[k] for n, k in keys_by_name.items()}).to_numpy().astype(int)
            pos = dose_index.get_indexer(pd.MultiIndex.from_arrays([d.drug.astype(str).to_numpy(),
                                                                    d.concentration.astype(float).round(6).to_numpy()]))
            di = np.where(pos >= 0, dose_vals[np.clip(pos, 0, None)], -1)
            gi = d.gene_name.astype(str).map(gidx).fillna(-1).astype(int).to_numpy()
            ok = (di >= 0) & (di <= 2) & (gi >= 0); unknown_dose += int((di < 0).sum())
            np.add.at(S, (ki[ok], di[ok], gi[ok]), d.stat.to_numpy()[ok])
            np.add.at(N, (ki[ok], di[ok], gi[ok]), 1)
        if i % 100 == 0 or i == len(parts): print(f"  {i}/{len(parts)} Tahoe files scanned")
    if unknown_dose: print(f"  note: {unknown_dose:,} rows had a dose not found in the step-8 table (ignored)")
    X = np.divide(S, N, out=np.zeros_like(S), where=N > 0).astype(np.float32)
    used = N.sum(axis=(0, 1)) > 0
    X = X[:, :, used]; genes = [g for g, u in zip(genes, used) if u]
    np.savez_compressed(cache, X=X, keys=np.array(keys, dtype=object), genes=np.array(genes, dtype=object))
    return X, keys, genes

def chem_features(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    mg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    m = Chem.MolFromSmiles(str(smiles))
    if m is None: return None
    return [f(m) for _, f in Descriptors.descList] + list(mg.GetFingerprintAsNumPy(m))

# ------------------------------------------------------------------ main
def main(a):
    from sklearn.linear_model import LogisticRegression
    from sklearn.decomposition import PCA
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import RepeatedStratifiedKFold
    from sklearn.metrics import roc_auc_score
    os.makedirs(a.out, exist_ok=True)
    log = open(os.path.join(a.out, "summary.txt"), "w", encoding="utf-8")
    def say(*x):
        s = " ".join(str(i) for i in x); print(s); log.write(s + "\n"); log.flush()

    dm = pd.read_csv("../data/raw/tahoe_metadata/drug_metadata.csv").drop_duplicates("drug")
    dm["key"] = dm.drug.map(m02.norm)
    appr = dm[dm["human-approved"].astype(str).str.lower() == "yes"]
    d1 = pd.read_csv("../data/processed/shared_drugs.csv").drop_duplicates("key")
    d1_keys = set(d1[d1.dict_concern.astype(str).str.strip().str.lower().isin(["most", "no"])].key)
    cand = appr[~appr.key.isin(d1_keys)].drop_duplicates("key")
    say(f"FDA-approved Tahoe compounds: {appr.key.nunique()} | excluded (used in D1): {len(d1_keys & set(appr.key))} "
        f"| C1 candidates: {len(cand)}")

    key = (a.api_key or os.getenv("OPENFDA_KEY") or "").strip().strip('"').strip("'")
    say("openFDA API key:", "provided" if key else "NOT provided (works, but slower and rate-limited)")
    total = make_total(key or None)
    # ---- connection self-check (stops early with a clear message if openFDA or the query format fails)
    n_all = total(ALL); n_su = total('patient.drug.openfda.generic_name:"sunitinib"')
    n_su_c = total(f'patient.drug.openfda.generic_name:"sunitinib" AND {CARD}')
    say(f"openFDA self-check: all reports = {n_all:,} | sunitinib = {n_su:,} | sunitinib + cardiac = {n_su_c:,}")
    if not (n_all > 1_000_000 and n_su > 1_000 and 0 < n_su_c < n_su):
        say("SELF-CHECK FAILED: openFDA did not return sensible counts. Check internet access / try again later.")
        return
    if a.check_only:
        say("Self-check passed (--check-only): stopping here."); return
    f = faers_counts(dict(zip(cand.key, cand.drug)), total, os.path.join(a.out, "faers_counts_v2.csv"))
    f = add_ror(f[f.key.isin(set(cand.key))].drop_duplicates("key").copy())
    f.to_csv(os.path.join(a.out, "faers_ror.csv"), index=False)
    f = f[f.n_drug >= a.min_reports]
    say(f"C1 drugs with >= {a.min_reports} FAERS reports: {len(f)} | cardiac signal: {int(f.signal.sum())} "
        f"| no signal: {int((f.signal == 0).sum())}")
    if f.signal.sum() < 8 or (f.signal == 0).sum() < 8:
        say("Too few drugs in one class for a reliable test - stopping (report descriptively)."); return

    names_by_key = dict(zip(dm.drug, dm.key))
    wanted = {n: k for n, k in names_by_key.items() if k in set(f.key)}
    X3, tkeys, genes = tahoe_features(a.parts, wanted, os.path.join(a.out, "tahoe_features_C1.npz"))
    keep = (X3[:, 2, :] != 0).mean(0) >= 0.9
    T = X3[:, :, keep].reshape(len(tkeys), -1)
    tk = {k: i for i, k in enumerate(tkeys)}
    smi = dict(zip(cand.key, cand.canonical_smiles))
    rows = [(k, chem_features(smi.get(k)), tk.get(k)) for k in f.key]
    no_smiles = [k for k, c, t in rows if c is None]; no_tahoe = [k for k, c, t in rows if t is None]
    if no_smiles: say(f"  dropped {len(no_smiles)} drug(s) with unreadable SMILES: {', '.join(no_smiles)}")
    if no_tahoe: say(f"  dropped {len(no_tahoe)} drug(s) without Tahoe features: {', '.join(no_tahoe)}")
    rows = [(k, c, t) for k, c, t in rows if c is not None and t is not None]
    keys = [k for k, _, _ in rows]
    Xc = np.nan_to_num(np.array([c for _, c, _ in rows], float), posinf=0, neginf=0)
    Xt = T[[t for _, _, t in rows]]
    y = f.set_index("key").loc[keys, "signal"].values
    say(f"Analysed: {len(y)} drugs ({y.sum()} signal) | Tahoe genes kept: {keep.sum()}")

    sets = {"chemistry": (Xc, None), "tahoe": (Xt, "pca"), "chemistry+tahoe": (np.hstack([Xc, Xt]), "pca_t")}
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=2026)
    oof = {s: np.zeros(len(y)) for s in sets}; cnt = np.zeros(len(y))
    for tr, te in cv.split(Xc, y):
        cnt[te] += 1
        for s, (X, mode) in sets.items():
            if mode is None:
                mdl = make_pipeline(StandardScaler(), LogisticRegression(C=0.01, max_iter=5000))
                mdl.fit(X[tr], y[tr]); oof[s][te] += mdl.predict_proba(X[te])[:, 1]
            elif mode == "pca":
                mdl = make_pipeline(StandardScaler(), PCA(10, random_state=0), LogisticRegression(C=0.1, max_iter=5000))
                mdl.fit(X[tr], y[tr]); oof[s][te] += mdl.predict_proba(X[te])[:, 1]
            else:   # chemistry block + Tahoe PCA-10 fitted inside the fold
                nc = Xc.shape[1]
                sc_t = StandardScaler().fit(X[tr, nc:]); pca = PCA(10, random_state=0).fit(sc_t.transform(X[tr, nc:]))
                Ztr = np.hstack([X[tr, :nc], pca.transform(sc_t.transform(X[tr, nc:]))])
                Zte = np.hstack([X[te, :nc], pca.transform(sc_t.transform(X[te, nc:]))])
                mdl = make_pipeline(StandardScaler(), LogisticRegression(C=0.01, max_iter=5000))
                mdl.fit(Ztr, y[tr]); oof[s][te] += mdl.predict_proba(Zte)[:, 1]
    for s in oof: oof[s] /= cnt
    auc = {s: roc_auc_score(y, oof[s]) for s in oof}
    rng = np.random.default_rng(2026); diffs, bs = [], {s: [] for s in oof}
    for _ in range(2000):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max(): continue
        for s in oof: bs[s].append(roc_auc_score(y[i], oof[s][i]))
        diffs.append(bs["chemistry+tahoe"][-1] - bs["chemistry"][-1])
    for s in oof:
        lo, hi = np.percentile(bs[s], [2.5, 97.5]); say(f"  AUROC {s:16s} = {auc[s]:.3f}  (95% CI {lo:.3f}-{hi:.3f})")
    lo, hi = np.percentile(diffs, [2.5, 97.5]); d = auc["chemistry+tahoe"] - auc["chemistry"]
    say(f"\nPRIMARY: AUROC(chemistry+Tahoe) - AUROC(chemistry) = {d:+.3f}  (95% CI {lo:+.3f} to {hi:+.3f})")
    say("DECISION:", "Tahoe ADDS information (lower CI > 0)" if lo > 0 else "Tahoe does NOT add information (lower CI <= 0)")
    pd.DataFrame({"key": keys, "signal": y, **{f"oof_{s}": oof[s] for s in oof}}).to_csv(os.path.join(a.out, "oof_predictions.csv"), index=False)
    log.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="../data/processed/pseudobulk_parts_all")
    ap.add_argument("--out", default="../results/confirmatory/C1")
    ap.add_argument("--min-reports", dest="min_reports", type=int, default=100)
    ap.add_argument("--check-only", dest="check_only", action="store_true", help="only test the openFDA connection")
    ap.add_argument("--api-key", dest="api_key", default=None, help="openFDA API key (alternative to $env:OPENFDA_KEY)")
    main(ap.parse_args())

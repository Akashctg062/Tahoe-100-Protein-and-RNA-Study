"""
Step 6 - TEMPORAL VALIDATION (dry-lab substitute for wet-lab validation).
Question: do model scores computed WITHOUT post-2019 information anticipate cardiac
safety signals that only appeared in FAERS from 2020 onward?

  1. FAERS ROR with reports received 2004-2019  -> 'early' signal
  2. FAERS ROR with reports received 2020-today -> 'late' signal
  3. 'Emerging signal' = no early signal but a late signal
  4. AUROC of your model's score for emerging vs. never-signalled drugs

Important: for this test, train the model only on labels that existed before 2020
(e.g., DICTrank is based on labels up to ~2023, so for a strict test use FAERS-early
labels for training). Document this choice in Methods.

Run:  python 06_temporal_faers_validation.py --scores ../results/oof_predictions_most_all.csv --col "combined|logistic"
"""
import argparse, math, os, time
import pandas as pd, requests
from sklearn.metrics import roc_auc_score

API = "https://api.fda.gov/drug/event.json"
PTS = ["Cardiac failure", "Cardiac failure congestive", "Left ventricular dysfunction",
       "Ejection fraction decreased", "Cardiomyopathy", "Cardiotoxicity", "Myocarditis",
       "Electrocardiogram QT prolonged", "Torsade de pointes", "Ventricular tachycardia",
       "Atrial fibrillation", "Myocardial infarction"]
CARD = "patient.reaction.reactionmeddrapt.exact:(" + " ".join(f'"{p}"' for p in PTS) + ")"
KEY = os.getenv("OPENFDA_KEY")

def total(q):
    p = {"search": q, "limit": 1}
    if KEY: p["api_key"] = KEY
    r = requests.get(API, params=p, timeout=60); time.sleep(0.3)
    if r.status_code == 404: return 0
    r.raise_for_status(); return r.json()["meta"]["results"]["total"]

def ror(drug, period):
    d = f'patient.drug.openfda.generic_name.exact:"{drug.upper()}"'
    N = total(period); C = total(f"{period} AND {CARD}")
    n = total(f"{period} AND {d}"); a = total(f"{period} AND {d} AND {CARD}")
    b, c = n - a, C - a; dd = N - n - c
    if min(a, b, c, dd) <= 0: return a, None
    r = (a / b) / (c / dd); se = math.sqrt(1/a + 1/b + 1/c + 1/dd)
    return a, math.exp(math.log(r) - 1.96 * se)          # lower 95% CI

def main(a):
    sc = pd.read_csv(a.scores, index_col=0)
    names = (pd.read_csv("../data/processed/shared_drugs.csv").drop_duplicates("key")
               .set_index("key").drug_cardiac)                    # scores are indexed by compound key
    early = "receivedate:[20040101 TO 20191231]"; late = f"receivedate:[20200101 TO {pd.Timestamp.today():%Y%m%d}]"
    rows = []
    for drug in sc.index:
        g = str(names.get(drug, drug))
        a1, lo1 = ror(g, early); a2, lo2 = ror(g, late)
        rows.append(dict(drug=drug, generic=g, early_n=a1, early_lo=lo1, late_n=a2, late_lo=lo2))
        print(rows[-1])
    r = pd.DataFrame(rows).set_index("drug")
    r["early_sig"] = (r.early_lo > 1) & (r.early_n >= 3)
    r["late_sig"] = (r.late_lo > 1) & (r.late_n >= 3)
    r["score"] = sc[a.col]
    test = r[~r.early_sig]                                  # drugs without an early signal
    y = test.late_sig.astype(int)
    print(f"\nDrugs without early signal: {len(test)} | emerging signals: {y.sum()}")
    if 0 < y.sum() < len(y):
        print(f"Temporal AUROC (score -> emerging signal): {roc_auc_score(y, test.score):.3f}")
    r.to_csv("../results/temporal_validation.csv"); print("saved ../results/temporal_validation.csv")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default="../results/oof_predictions_most_all.csv")
    ap.add_argument("--col", default="combined|logistic")
    main(ap.parse_args())

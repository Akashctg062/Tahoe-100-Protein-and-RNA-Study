"""
Step 4 - Clinical cardiotoxicity labels from FDA FAERS via openFDA (public, de-identified,
no IRB). Computes the Reporting Odds Ratio (ROR) and 95% CI for a cardiac event group.

  ROR = (a/b) / (c/d)
  a = reports with drug & cardiac event      b = drug & no cardiac event
  c = other drugs & cardiac event            d = other drugs & no cardiac event
Signal if lower 95% CI > 1 and a >= 3.

Run:  python 04_faers_cardiac_ror.py --drugs ../data/processed/shared_drugs.csv
Free API key (recommended, 120k requests/day): https://open.fda.gov/apis/authentication/
"""
import argparse, math, os, time
import pandas as pd, requests

API = "https://api.fda.gov/drug/event.json"
CARDIAC_PTS = ["Cardiac failure", "Cardiac failure congestive", "Left ventricular dysfunction",
               "Ejection fraction decreased", "Cardiomyopathy", "Cardiotoxicity", "Myocarditis",
               "Electrocardiogram QT prolonged", "Torsade de pointes", "Ventricular tachycardia",
               "Atrial fibrillation", "Myocardial infarction"]
KEY = os.getenv("OPENFDA_KEY")

def total(search):
    p = {"search": search, "limit": 1}
    if KEY: p["api_key"] = KEY
    r = requests.get(API, params=p, timeout=60)
    if r.status_code == 404: return 0
    r.raise_for_status(); time.sleep(0.3)
    return r.json()["meta"]["results"]["total"]

def main(drug_csv, out):
    d = pd.read_csv(drug_csv)
    col = "drug_cardiac" if "drug_cardiac" in d.columns else "drug"
    cardiac = "patient.reaction.reactionmeddrapt.exact:(" + " ".join(f'"{p}"' for p in CARDIAC_PTS) + ")"
    N = total("_exists_:safetyreportid"); C = total(cardiac)
    rows = []
    for name in sorted(set(d[col].dropna())):
        drug = f'patient.drug.openfda.generic_name.exact:"{str(name).upper()}"'
        n_drug = total(drug); a = total(f"{drug} AND {cardiac}")
        b, c = n_drug - a, C - a; dd = N - n_drug - c
        if min(a, b, c, dd) <= 0:
            rows.append(dict(drug=name, a=a, n_drug=n_drug, ROR=None, lo=None, hi=None)); continue
        ror = (a / b) / (c / dd); se = math.sqrt(1/a + 1/b + 1/c + 1/dd)
        rows.append(dict(drug=name, a=a, n_drug=n_drug, ROR=ror,
                         lo=math.exp(math.log(ror) - 1.96*se), hi=math.exp(math.log(ror) + 1.96*se)))
        print(rows[-1])
    res = pd.DataFrame(rows); res["signal"] = (res.lo > 1) & (res.a >= 3)
    res.to_csv(out, index=False); print("saved ->", out)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--drugs", default="../data/processed/shared_drugs.csv")
    ap.add_argument("--out", default="../data/processed/faers_cardiac_ror.csv")
    a = ap.parse_args(); main(a.drugs, a.out)

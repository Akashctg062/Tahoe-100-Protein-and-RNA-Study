"""Step 10b - C1 post-hoc sensitivity (Amendment 2): exclude drugs with a cardiovascular indication
(confounding by indication in FAERS) and recompute AUROCs from the frozen out-of-fold predictions."""
import pandas as pd, numpy as np, re, sys
from sklearn.metrics import roc_auc_score
U = sys.argv[1] if len(sys.argv) > 1 else "../"      # default: run from the scripts folder
o = pd.read_csv(U + "results/confirmatory/C1/oof_predictions.csv")
f = pd.read_csv(U + "results/confirmatory/C1/faers_ror.csv")[["key", "drug_tahoe"]]
dm = pd.read_csv(U + "data/raw/tahoe_metadata/drug_metadata.csv").drop_duplicates("drug")[["drug", "gpt-notes-approval"]]
x = o.merge(f, on="key").merge(dm, left_on="drug_tahoe", right_on="drug", how="left")
cv = r"hypertens|heart failure|angina|arrhythm|atrial|cardi|myocard|coronary|thromb|stroke|edema|hypotens|vasopress|blood pressure|platelet|anticoag"
x["cardiovascular_indication"] = x["gpt-notes-approval"].astype(str).str.lower().str.contains(cv)
print("cardiovascular-indication drugs:", int(x.cardiovascular_indication.sum()), "of", len(x))
print("  ", ", ".join(sorted(x[x.cardiovascular_indication].drug_tahoe)))
print("signal rate: CV-indication", round(x[x.cardiovascular_indication].signal.mean(), 2), "| other", round(x[~x.cardiovascular_indication].signal.mean(), 2))
rng = np.random.default_rng(2026)
def boot(sub, col):
    y = sub.signal.values; p = sub[col].values; bs = []
    for _ in range(2000):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max(): bs.append(roc_auc_score(y[i], p[i]))
    return roc_auc_score(y, p), *np.percentile(bs, [2.5, 97.5])
for name, sub in [("all drugs (frozen C1)", x), ("excluding cardiovascular-indication drugs", x[~x.cardiovascular_indication])]:
    print(f"\n{name}: n = {len(sub)}, signal = {int(sub.signal.sum())}")
    for col in ["oof_chemistry", "oof_tahoe", "oof_chemistry+tahoe"]:
        a, lo, hi = boot(sub, col); print(f"   AUROC {col[4:]:16s} {a:.3f} (95% CI {lo:.3f}-{hi:.3f})")
x.to_csv(U + "results/confirmatory/C1/c1_with_indication.csv", index=False)

"""Independent recomputation of the numbers displayed in v2 figures, from per-drug / raw files (not the summary tables)."""
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, spearmanr
from sklearn.metrics import roc_auc_score
R = "/home/claude/v2/results/"; A1 = "/home/claude/mirror/results/confirmatory/amendment1/"
fails = []
def chk(name, got, exp, tol=0.006):
    ok = abs(got - exp) <= tol; print(f"{'OK ' if ok else 'FAIL'} {name}: recomputed {got:+.4f} vs shown {exp:+.4f}")
    if not ok: fails.append(name)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / (np.arange(len(p)) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
# Fig 6a / b: medians of calibrated rho from per-drug
E_ = pd.read_csv(R + "T2EF/per_drug_E_ctrp.csv"); F_ = pd.read_csv(R + "T2EF/per_drug_F_procan.csv"); A = pd.read_csv(A1 + "per_drug_calibrated.csv")
med = E_.groupby("set").calibrated_rho.median()
for st, v in [("RNA", .22), ("RPPA", .18), ("Response", .10)]: chk(f"Fig6 CTRP {st} median", med[st], v)
print("   CTRP sets:", med.round(3).to_dict())
fp = F_[F_.set == "Protein"].groupby("source").calibrated_rho.median()
print("   ProCan sets present:", F_.set.unique())
for src, v in [("PRISM", .046), ("GDSC", .293), ("CTRP", .282)]: chk(f"Fig6 ProCan {src}", fp[src], v, 0.002)
am = A.groupby(["source", "set"]).calibrated_rho.median()
for (s, st), v in [(("GDSC", "RPPA"), .17), (("PRISM", "Protein_MS"), .13), (("PRISM", "Response"), -.05), (("GDSC", "RNA"), .00)]: chk(f"Fig6 C3 {s} {st}", am[(s, st)], v)
# Fig 6b pairwise FDRs in CTRP: Response vs RNA, RPPA vs RNA (two-sided Wilcoxon of paired differences, BH over the family in the file)
W = E_.pivot_table(index="drug_key", columns="set", values="calibrated_rho")
d1 = W.Response - W.RNA; d2 = W.RPPA - W.RNA
print(f"   CTRP Response-RNA median {d1.median():+.3f}, Wilcoxon two-sided p {wilcoxon(d1.dropna()).pvalue:.4f}; RPPA-RNA median {d2.median():+.3f}, p {wilcoxon(d2.dropna()).pvalue:.4f}")
sp = pd.read_csv(R + "T2EF/summary_pairs_E_ctrp.csv"); print(sp[["comparison", "median_difference", "p_two_sided", "fdr_bh"]].round(4).to_string(index=False))
ss = pd.read_csv(R + "T2EF/summary_sets_E_ctrp.csv")   # BH is joint over the 6 set tests + 7 pairwise tests (Amendment-1 procedure)
q = bh(np.r_[ss.p_gt0.values, sp.p_two_sided.values])[len(ss):]
chk("Fig6b joint BH (Response vs RNA)", q[sp.comparison.eq("Response vs RNA").values][0], 0.040, 0.0005)
chk("Fig6b joint BH (RPPA vs RNA)", q[sp.comparison.eq("RPPA vs RNA").values][0], 0.0035, 0.0002)
# GDSC ProCan - RNA paired
G = F_[F_.source == "GDSC"].pivot_table(index="drug_key", columns="set", values="calibrated_rho"); pc = "Protein"
chk("Fig6c ProCan-RNA GDSC median diff", (G[pc] - G.RNA).median(), 0.180, 0.003)
# Fig 6e top vs mid dose correlation
D = pd.read_csv(R + "T2EF/per_drug_T1D_middose.csv")
for src, v in [("PRISM", .12), ("GDSC", .08)]:
    t = A[(A.source == src) & (A.set == "Response")].set_index("drug_key").calibrated_rho
    m = D[(D.source == src) & (D.set == "Response")].set_index("drug_key").calibrated_rho
    j = t.index.intersection(m.index); chk(f"Fig6e top vs mid r {src} (n={len(j)})", np.corrcoef(t[j], m[j])[0, 1], v)
# Fig 5: null bias from raw grid; calibration from raw reps
g = pd.read_csv(R + "T1B/null_grid_raw.csv")
for n, v in [(10, -.94), (20, -.84), (27, -.77), (40, -.64), (60, -.51)]: chk(f"Fig5 null k5 a1 n={n}", g[(g.k == 5) & (g.alpha == 1) & (g.n == n)].rho.mean(), v)
chk("Fig5 intercept-only n=27", g[(g.k == 0) & (g.n == 27)].rho.mean(), -1.0, 1e-9)
c = pd.read_csv(R + "T1B/calibration_raw.csv")
print("   max |calibrated| under null:", c[c.true_r == 0].groupby(["scheme", "n"]).calibrated.mean().abs().max().round(3))
for sch, n, v in [("loo", 40, .54), ("loo", 60, .59), ("5fold", 40, .29), ("5fold", 60, .37)]:
    chk(f"Fig5 calibrated {sch} n={n} r=0.4", c[(c.scheme == sch) & (c.n == n) & (c.true_r == 0.4)].calibrated.mean(), v)
# Fig 3a-c: AUROCs from OOF predictions
for task in ["C1", "D1"]:
    o = pd.read_csv(R + f"T1A1/oof_{task}.csv"); s = pd.read_csv(R + f"T1A1/summary_{task}.csv")
    for _, r in s.iterrows():
        a = roc_auc_score(o.y, o[f"{r.family}|tahoe_{r.representation}"]); b = roc_auc_score(o.y, o[f"{r.family}|chemistry"]); cc = roc_auc_score(o.y, o[f"{r.family}|chem+tahoe_{r.representation}"])
        if abs(a - r.auroc_tahoe) > 1e-6 or abs(b - r.auroc_chem) > 1e-6 or abs(cc - b - r.delta) > 1e-6: fails.append(f"{task} {r.family} {r.representation}")
    print(f"OK  Fig3 {task}: all 8 AUROC/Δ values re-derived from OOF (n = {len(o)}, positives = {int(o.y.sum())}); min Holm p = {s.p_holm.min():.3f}")
# Fig 3d hERG
h = pd.read_csv(R + "T2H2/oof_regression.csv")
chk("Fig3d hERG rho chem", spearmanr(h.oof_chem, h.pIC50)[0], 0.027, 0.001); chk("Fig3d hERG rho chem+Tahoe", spearmanr(h.oof_chem_tahoe, h.pIC50)[0], 0.002, 0.001)
print(f"   hERG n = {len(h)}, pIC50 >= 5: {(h.pIC50 >= 5).sum()}")
# Fig 3e
M = pd.read_csv(R + "T2H1/corr_matrix_averaged_DEGs.csv", index_col=0).values; n = len(M)
chk("Fig3e same-drug mean r", np.diag(M).mean(), 0.062, 0.001); chk("Fig3e different mean r", M[~np.eye(n, dtype=bool)].mean(), 0.013, 0.001)
# S2: frozen + best-effort medians from per-drug
for m in ["enet", "rf", "ridge_krr"]:
    pdd = pd.read_csv(R + f"T1A2/per_drug_{m}.csv"); ss = pd.read_csv(R + f"T1A2/summary_sets_{m}.csv")
    rec = pdd.groupby(["model", "source", "set"]).calibrated_rho.median()
    mx = max(abs(rec[(r.model, r.source, r.set)] - r.median_calibrated) for r in ss.itertuples())
    print(f"{'OK ' if mx < 1e-9 else 'FAIL'} S2 {m}: per-drug medians match summary (max diff {mx:.1e})")
    if mx >= 1e-9: fails.append(m)
u = pd.read_csv(R + "T1C/per_unit.csv"); chk("S2c pooled median concordance", u.concordance.median(), 0.511, 0.001)
print("\nFAILS:", fails if fails else "none")

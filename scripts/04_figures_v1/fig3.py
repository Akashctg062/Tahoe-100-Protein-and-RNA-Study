"""Figure 3 - Early cell counts and response magnitude do not track drug sensitivity (D3, D4, C2)."""
import sys; sys.path.insert(0, "/home/claude/figs")
from figcommon import *
from scipy.stats import spearmanr, wilcoxon
style()
F = U + "results/feasibility/"; RS = U + "results/response_vs_sensitivity/"
v = pd.read_csv(ANA + "tahoe_viability_v2.csv")            # pool-share relative fitness (D3)
p = pd.read_csv(F + "pairs_with_labels.csv")
d4 = pd.read_csv(RS + "per_drug_correlations.csv")
c2 = pd.read_csv(U + "results/confirmatory/C2_per_drug.csv")
LAB = {"prism_auc": ("PRISM", BLUE), "gdsc_auc": ("GDSC", ORANGE)}
rng = np.random.default_rng(7)
src = []

top = v[v.dose_rank == 3][["drug", "Cell_Name_Vevo", "rel_fitness"]]
q = p.merge(top, on=["drug", "Cell_Name_Vevo"])

fig = plt.figure(figsize=(7.1, 5.6))
gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.55, left=0.08, right=0.985, top=0.93, bottom=0.12)

# a - pooled scatter, Tahoe pool-share fitness vs PRISM AUC
ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", "D3: 24-h cell counts vs PRISM")
s = q.dropna(subset=["prism_auc"])
ax.scatter(s.rel_fitness, s.prism_auc, s=3, color=BLUE, alpha=0.25, lw=0, rasterized=True)
r, pv = spearmanr(s.rel_fitness, s.prism_auc)
ax.text(0.03, 0.04, f"Spearman ρ = {r:.2f}\n{len(s):,} drug–line pairs\n{s.drug.nunique()} drugs", transform=ax.transAxes, fontsize=5.8, color=INK)
ax.set_xlabel("Tahoe relative fitness, top dose\n(log2 pool-share, treated / control)"); ax.set_ylabel("PRISM AUC (lower = sensitive)")
ax.set_xlim(np.percentile(s.rel_fitness, 0.5), np.percentile(s.rel_fitness, 99.5))

# b - within-drug correlations, both screens
ax = fig.add_subplot(gs[0, 1]); panel(ax, "b", "D3: within-drug agreement")
for i, (lab, (nm, col)) in enumerate(LAB.items()):
    s = q.dropna(subset=[lab])
    pdg = s.groupby("drug").apply(lambda x: spearmanr(x.rel_fitness, x[lab])[0] if len(x) >= 8 else np.nan).dropna()
    for d_, r_ in pdg.items(): src.append(dict(panel="3b", screen=nm, drug=d_, rho=r_))
    ax.scatter(i + rng.uniform(-0.17, 0.17, len(pdg)), pdg, s=8, color=col, alpha=0.8, lw=0.3, edgecolor="white", zorder=3)
    ax.hlines(pdg.median(), i - 0.3, i + 0.3, color=INK, lw=1.2, zorder=4)
    w = wilcoxon(pdg)[1]
    ax.text(i, 0.9, f"median {pdg.median():+.2f}\np = {w:.2f}\n{len(pdg)} drugs", ha="center", fontsize=5.6, color=INK2, va="top")
ax.axhline(0, color=GRAY, lw=0.7, ls="--"); ax.set_xticks([0, 1]); ax.set_xticklabels(["PRISM", "GDSC"])
ax.set_xlim(-0.6, 1.6); ax.set_ylim(-1, 1); ax.set_ylabel("Spearman ρ across cell lines\n(fitness vs AUC, per drug)")

# c - internal consistency: dose 2 vs dose 3, salt-form replicates
ax = fig.add_subplot(gs[0, 2]); panel(ax, "c", "D3: internal consistency")
wv = v.pivot_table(index=["drug", "Cell_Name_Vevo"], columns="dose_rank", values="rel_fitness")
dc = wv.groupby(level=0).apply(lambda x: spearmanr(x[2], x[3])[0]).dropna()
dup = v.groupby("drug_key").drug.nunique(); salt = []
for k in dup[dup > 1].index:
    a = v[(v.drug_key == k) & (v.dose_rank == 3)]; ds = sorted(a.drug.unique())
    m = a[a.drug == ds[0]].merge(a[a.drug == ds[1]], on="Cell_Name_Vevo")
    salt.append(spearmanr(m.rel_fitness_x, m.rel_fitness_y)[0])
bins = np.linspace(-0.8, 1, 28)
ax.hist(dc, bins=bins, color=GRAY, alpha=0.85, label=f"Mid vs top dose, same drug\n(median {dc.median():.2f}, {len(dc)} drugs)")
ax.hist(salt, bins=bins, color=AQUA, alpha=0.95, label=f"Two salt forms, top dose\n(median {np.median(salt):.2f}, {len(salt)} pairs)")
ax.axvline(0, color=INK2, lw=0.6, ls="--")
ax.set_xlabel("Spearman ρ of relative fitness\nacross cell lines"); ax.set_ylabel("Drugs")
ax.legend(fontsize=5.4, loc="upper left", bbox_to_anchor=(-0.02, 1.04)); ax.set_ylim(0, ax.get_ylim()[1] * 1.55)

# d - D4: response magnitude vs sensitivity, per drug
ax = fig.add_subplot(gs[1, 0]); panel(ax, "d", "D4: response magnitude vs sensitivity")
FE = [("n_de", "No. DE genes"), ("rms_lfc", "RMS log2FC"), ("rms_stat", "RMS statistic")]
xt, xl = [], []
for j, (lab, (nm, col)) in enumerate(LAB.items()):
    for k, (f, fl) in enumerate(FE):
        xx = j * 3.6 + k
        r_ = d4[(d4.label == lab) & (d4.feature == f)].rho
        ax.scatter(xx + rng.uniform(-0.2, 0.2, len(r_)), r_, s=5, color=col, alpha=0.6, lw=0, zorder=3)
        ax.hlines(r_.median(), xx - 0.33, xx + 0.33, color=INK, lw=1.2, zorder=4)
        xt.append(xx); xl.append(fl)
    ax.text(j * 3.6 + 1, 0.97, f"{nm} ({r_.size} drugs)", ha="center", fontsize=5.8, fontweight="bold", va="top")
ax.axhline(0, color=GRAY, lw=0.7, ls="--"); ax.set_xticks(xt); ax.set_xticklabels(xl, rotation=40, ha="right", fontsize=5.6)
ax.set_ylim(-1, 1); ax.set_ylabel("Spearman ρ across cell lines\n(magnitude vs AUC, per drug)")

# e - C2: proliferation-program drop vs sensitivity
ax = fig.add_subplot(gs[1, 1]); panel(ax, "e", "C2: proliferation program (pre-specified)")
for i, (nm, col) in enumerate([("PRISM", BLUE), ("GDSC", ORANGE)]):
    r_ = c2[c2.source == nm].rho
    ax.scatter(i + rng.uniform(-0.17, 0.17, len(r_)), r_, s=8, color=col, alpha=0.8, lw=0.3, edgecolor="white", zorder=3)
    ax.hlines(r_.median(), i - 0.3, i + 0.3, color=INK, lw=1.2, zorder=4)
    ax.text(i, 0.95, f"median {r_.median():+.3f}\n{len(r_)} drugs", ha="center", fontsize=5.6, color=INK2, va="top")
ax.axhline(0, color=GRAY, lw=0.7, ls="--")
ax.axhline(0.025, color=MAGENTA, lw=0.9, ls=":")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], color=MAGENTA, ls=":", lw=0.9, label="Pre-specified threshold (+0.025)")], fontsize=5.4, loc="lower right")
ax.set_xticks([0, 1]); ax.set_xticklabels(["PRISM", "GDSC"]); ax.set_xlim(-0.6, 1.7); ax.set_ylim(-0.75, 1.0)
ax.set_ylabel("Spearman ρ across cell lines\n(program drop vs sensitivity)")

# f - C2 agreement between screens (post hoc)
ax = fig.add_subplot(gs[1, 2]); panel(ax, "f", "C2 per drug: PRISM vs GDSC (post hoc)")
w = c2.pivot(index="drug_key", columns="source", values="rho").dropna()
r, pv = spearmanr(w.PRISM, w.GDSC)
ax.scatter(w.PRISM, w.GDSC, s=14, color=VIOLET, edgecolor="white", lw=0.4, zorder=3)
ax.axhline(0, color=GRAY, lw=0.6, ls="--"); ax.axvline(0, color=GRAY, lw=0.6, ls="--")
lim = (-0.55, 0.55); ax.plot(lim, lim, color=LIGHTGRAY, lw=0.7, zorder=1); ax.set_xlim(lim); ax.set_ylim(lim)
OFF = {"palbociclib": (4, -1, "left"), "crizotinib": (-4, 0, "right"), "dabrafenib": (4, 0, "left"), "gefitinib": (3, 3, "left"),
       "erlotinib": (3, 3, "left")}
for d_, (dx, dy, ha) in OFF.items():
    if d_ in w.index:
        ax.annotate(d_.capitalize(), (w.loc[d_, "PRISM"], w.loc[d_, "GDSC"]), xytext=(dx, dy), textcoords="offset points", fontsize=5, color=INK2, ha=ha, va="center")
ax.text(0.97, 0.03, f"Spearman ρ = {r:.2f} (p = {pv:.3f})\n{len(w)} drugs in both screens", transform=ax.transAxes, fontsize=5.6, va="bottom", ha="right")
ax.set_xlabel("C2 ρ in PRISM"); ax.set_ylabel("C2 ρ in GDSC"); ax.set_aspect("equal")
w.to_csv(OUT + "source_Fig3f.csv")
pd.DataFrame(src).to_csv(OUT + "source_Fig3b.csv", index=False)
save(fig, "Figure3_viability_sensitivity")
print("dose consistency", dc.median(), "salt", salt, "C2 agree", r, pv)

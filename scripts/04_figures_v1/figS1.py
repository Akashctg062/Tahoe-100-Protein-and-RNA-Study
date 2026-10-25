"""Supplementary Figure 1 - Sensitivity analyses: Amendment 2 (C1), cross-screen reproducibility (C3), pairwise robustness."""
import sys; sys.path.insert(0, "/home/claude/figs")
from figcommon import *
from sklearn.metrics import roc_auc_score
style()
ci = pd.read_csv(ANA + "c1/c1_with_indication.csv")
xs = pd.read_csv(OUT + "source_crossscreen_C3.csv")
rp = pd.read_csv(CONF + "robustness/robustness_all_pairs.csv")
rng = np.random.default_rng(2026)
def boot(y, p, n=2000):
    b = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max(): b.append(roc_auc_score(y[i], p[i]))
    return roc_auc_score(y, p), *np.percentile(b, [2.5, 97.5])

fig = plt.figure(figsize=(7.1, 2.75))
gs = fig.add_gridspec(1, 3, wspace=0.75, left=0.08, right=0.985, top=0.86, bottom=0.3, width_ratios=[1, 1, 1.25])

# a - C1 with and without cardiovascular-indication drugs
ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", "C1: Amendment 2", x=-0.3)
M = [("oof_chemistry", "Chemistry", ORANGE), ("oof_tahoe", "Tahoe", BLUE), ("oof_chemistry+tahoe", "Chem. + Tahoe", AQUA)]
for k, (sub, nm) in enumerate([(ci, f"All drugs\n(n = {len(ci)})"), (ci[~ci.cardiovascular_indication], f"Excl. CV indication\n(n = {(~ci.cardiovascular_indication).sum()})")]):
    for j, (c, lab, col) in enumerate(M):
        a, lo, hi = boot(sub.signal.values, sub[c].values)
        ax.errorbar(k + (j - 1) * 0.22, a, yerr=[[a - lo], [hi - a]], fmt="o", ms=3.5, color=col, elinewidth=0.8, label=lab if k == 0 else None)
ax.axhline(0.5, color=GRAY, lw=0.7, ls="--"); ax.set_xticks([0, 1])
ax.set_xticklabels([f"All drugs\n(n = {len(ci)})", f"Excl. CV indication\n(n = {(~ci.cardiovascular_indication).sum()})"], fontsize=5.8)
ax.set_xlim(-0.5, 1.5); ax.set_ylim(0.25, 0.75); ax.set_ylabel("AUROC (95% bootstrap CI)")
ax.legend(fontsize=5.3, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=3, columnspacing=0.6, handletextpad=0.2)

# b - cross-screen reproducibility of per-drug calibrated rho
ax = fig.add_subplot(gs[0, 1]); panel(ax, "b", "C3: PRISM–GDSC agreement", x=-0.3)
SH = {"RNA": "RNA", "Protein_MS": "Protein (MS)", "RPPA": "Protein (RPPA)", "Response": "Tahoe response",
      "Response+RNA": "Response + RNA", "Response+Protein_MS": "Response + MS"}
ax.bar(range(len(xs)), xs.rho, color=[SET_COLORS[s] for s in xs.set], width=0.65)
for i, r in xs.iterrows(): ax.text(i, r.rho + (0.03 if r.rho >= 0 else -0.03), f"p={r.p:.2f}", ha="center", va="bottom" if r.rho >= 0 else "top", fontsize=4.8, color=INK2)
ax.axhline(0, color=INK2, lw=0.6)
ax.set_xticks(range(len(xs))); ax.set_xticklabels([SH[s] for s in xs.set], rotation=45, ha="right", fontsize=5.6)
ax.set_ylim(-0.4, 0.75); ax.set_ylabel(f"Spearman ρ of per-drug calibrated ρ\n({xs.n_drugs.iloc[0]} drugs in both screens)")

# c - pairwise comparisons across settings
ax = fig.add_subplot(gs[0, 2]); panel(ax, "c", "Pairwise comparisons across settings", x=-0.42)
CMP = ["RPPA vs RNA", "Protein_MS vs RNA", "Response vs RNA", "RPPA vs Protein_MS"]
CFG = {"k5_a1.0": "Primary", "k3_a1.0": "k = 3", "k10_a1.0": "k = 10", "k5_a0.1": "α = 0.1", "k5_a10.0": "α = 10"}
A1 = pd.read_csv(CONF + "amendment1/summary_pairs.csv").assign(config="k5_a1.0")
allp = pd.concat([A1, rp[rp.config != "k5_a1.0 (primary, 500 perm)"]])
allp["config"] = allp.config.str.split(" ").str[0]
for si, (src, mk) in enumerate([("PRISM", "o"), ("GDSC", "s")]):
    for ci_, cmp_ in enumerate(CMP):
        for kj, cf in enumerate(CFG):
            r = allp[(allp.source == src) & (allp.comparison == cmp_) & (allp.config == cf)]
            if not len(r): continue
            r = r.iloc[0]; y = ci_ + (si - 0.5) * 0.42 + (kj - 2) * 0.06
            ax.scatter(r.median_difference, y, marker=mk, s=9 if kj else 16, color=BLUE if src == "PRISM" else ORANGE,
                       edgecolor=INK if r.p < 0.05 else "none", lw=0.6, alpha=1 if kj == 0 else 0.6, zorder=3)
ax.axvline(0, color=GRAY, lw=0.7, ls="--")
ax.set_yticks(range(len(CMP))); ax.set_yticklabels([c.replace("Protein_MS", "MS").replace("Response", "Tahoe resp.") for c in CMP], fontsize=5.8); ax.invert_yaxis()
ax.set_xlabel("Median difference in calibrated ρ")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=BLUE, ms=3.5, label="PRISM"), Line2D([], [], marker="s", ls="", color=ORANGE, ms=3.5, label="GDSC"),
                   Line2D([], [], marker="o", ls="", mfc="white", mec=INK, ms=3.5, label="p < 0.05")],
          fontsize=5.3, loc="upper center", bbox_to_anchor=(0.45, -0.3), ncol=3)
save(fig, "FigureS1_sensitivity_analyses")

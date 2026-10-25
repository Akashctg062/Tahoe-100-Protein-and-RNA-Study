"""v2 Supplementary Figure S2 - Best-effort sensitivity models (T1-A2) and RPPA-target concordance (T1-C)."""
import sys; sys.path.insert(0, "/home/claude/v2/figs")
from figcommon import *
from matplotlib.colors import TwoSlopeNorm
style()
R = "/home/claude/v2/results/"; A1 = "/home/claude/mirror/results/confirmatory/amendment1/"
fr = pd.read_csv(A1 + "summary_sets.csv").rename(columns={"median_calibrated_rho": "median_calibrated"}).assign(model="frozen")
be = pd.concat([pd.read_csv(R + f"T1A2/summary_sets_{m}.csv") for m in ["ridge_krr", "enet", "rf"]])
S = pd.concat([fr[["model", "source", "set", "median_calibrated", "fdr_bh"]], be[["model", "source", "set", "median_calibrated", "fdr_bh"]]])
MODELS = [("frozen", "Ridge LOO\n(frozen)"), ("ridge", "Ridge-GCV"), ("krr", "Kernel\nridge"), ("enet", "Elastic\nnet"), ("rf", "Random\nforest")]
ORDER = ["RNA", "Protein_MS", "RPPA", "Response", "Response+RNA", "Response+Protein_MS"]
SH = {"RNA": "RNA", "Protein_MS": "Protein (MS)", "RPPA": "RPPA", "Response": "Tahoe response", "Response+RNA": "Response + RNA", "Response+Protein_MS": "Response + MS"}
fig = plt.figure(figsize=(7.1, 5.6))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.55, left=0.14, right=0.93, top=0.94, bottom=0.08, height_ratios=[1, 0.95])
norm = TwoSlopeNorm(0, -0.2, 0.2)
src_rows = []
for j, (src, letter) in enumerate([("PRISM", "a"), ("GDSC", "b")]):
    ax = fig.add_subplot(gs[0, j]); panel(ax, letter, f"{src}: median calibrated ρ by model family", x=-0.32)
    M = np.full((len(ORDER), len(MODELS)), np.nan); F = np.ones_like(M)
    for a, st in enumerate(ORDER):
        for b, (m, _) in enumerate(MODELS):
            r = S[(S.model == m) & (S.source == src) & (S.set == st)]
            if len(r): M[a, b] = r.median_calibrated.iloc[0]; F[a, b] = r.fdr_bh.iloc[0]
            src_rows.append(dict(source=src, set=st, model=m, median_calibrated=M[a, b], fdr=F[a, b]))
    im = ax.imshow(M, cmap="RdBu_r", norm=norm, aspect="auto")
    for a in range(len(ORDER)):
        for b in range(len(MODELS)):
            ax.text(b, a, f"{M[a,b]:+.2f}" + ("*" if F[a, b] < 0.05 else ""), ha="center", va="center", fontsize=5.4,
                    color="white" if abs(M[a, b]) > 0.15 else INK, fontweight="bold" if F[a, b] < 0.05 else "normal")
    ax.set_xticks(range(len(MODELS))); ax.set_xticklabels([x[1] for x in MODELS], fontsize=5.4)
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([SH[x] for x in ORDER], fontsize=5.8)
    ax.tick_params(length=0); [s_.set_visible(False) for s_ in ax.spines.values()]
    ax.axvline(0.5, color="white", lw=1.5)
cb = fig.colorbar(im, ax=fig.axes[-1], fraction=0.05, pad=0.03); cb.set_label("median calibrated ρ", fontsize=5.5); cb.ax.tick_params(labelsize=5)
fig.text(0.5, 0.505, "* FDR < 0.05 (Wilcoxon > 0, BH within model and screen). Frozen = Amendment 1 ridge LOO, PCA-5, α = 1, 500 permutations; "
         "best effort = 5-fold CV within drug, PCA-10, 50 permutations.", ha="center", fontsize=5.1, color=INK2)
pd.DataFrame(src_rows).to_csv(OUT + "source_FigS2ab.csv", index=False)

# c - RPPA-target concordance per drug x screen unit
u = pd.read_csv(R + "T1C/per_unit.csv"); nl = pd.read_csv(R + "T1C/null.csv").null_median
ax = fig.add_subplot(gs[1, 0]); panel(ax, "c", "RPPA-target concordance (46 units)", x=-0.32)
rng = np.random.default_rng(3)
for i, (src, mk, c) in enumerate([("PRISM", "o", BLUE), ("GDSC", "s", ORANGE)]):
    g = u[u.source == src]
    ax.scatter(i + rng.uniform(-0.17, 0.17, len(g)), g.concordance, marker=mk, s=11, color=c, edgecolor="white", lw=0.3, zorder=3)
    ax.hlines(g.concordance.median(), i - 0.3, i + 0.3, color=INK, lw=1)
    prev = 9
    for _, r in g.nlargest(3, "concordance").iterrows():
        yy = min(r.concordance, prev - 0.045); prev = yy
        ax.text(i + 0.22, yy, r.drug, fontsize=4.8, color=INK2, va="center")
ax.axhline(0.5, color=GRAY, lw=0.7, ls="--")
ax.set_xticks([0, 1]); ax.set_xticklabels([f"PRISM (n = {(u.source=='PRISM').sum()})", f"GDSC (n = {(u.source=='GDSC').sum()})"], fontsize=6)
ax.set_ylabel("Concordance (mean percentile of\ntarget-matched antibodies)"); ax.set_xlim(-0.5, 1.75); ax.set_ylim(0, 1.02)

ax = fig.add_subplot(gs[1, 1]); panel(ax, "d", "Pooled median vs permutation null", x=-0.32)
obs = u.concordance.median(); p = (np.sum(nl >= obs) + 1) / (len(nl) + 1)
ax.hist(nl, bins=40, color=LIGHTGRAY); ax.axvline(obs, color=YELLOW, lw=1.6, label=f"observed {obs:.3f}")
ax.axvline(0.5, color=INK2, lw=0.6, ls=":", label="chance 0.5")
ax.set_xlabel("Median concordance"); ax.set_ylabel("Permutations"); ax.legend(fontsize=5.4, loc="upper left")
ax.text(0.97, 0.95, f"one-sided p = {p:.3f}\n(not supported)", transform=ax.transAxes, ha="right", va="top", fontsize=5.6, color=INK2)
save(fig, "FigureS2_bestfit_rppa_targets")

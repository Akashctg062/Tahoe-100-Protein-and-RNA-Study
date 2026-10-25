"""Figure 4 - Which data predict which cell lines a drug kills? (C3, Amendment 1, robustness)."""
import sys; sys.path.insert(0, "/home/claude/figs")
from figcommon import *
from scipy.stats import spearmanr
style()
A1 = CONF + "amendment1/"; RB = CONF + "robustness/"
ss = pd.read_csv(A1 + "summary_sets.csv"); sp = pd.read_csv(A1 + "summary_pairs.csv")
pdc = pd.read_csv(A1 + "per_drug_calibrated.csv"); rob = pd.read_csv(RB + "robustness_all_sets.csv")
ORDER = ["RNA", "Protein_MS", "RPPA", "Response", "Response+RNA", "Response+Protein_MS"]
SHORT = {"RNA": "RNA", "Protein_MS": "Protein (MS)", "RPPA": "Protein (RPPA)", "Response": "Tahoe response",
         "Response+RNA": "Response + RNA", "Response+Protein_MS": "Response + MS"}
rng = np.random.default_rng(11)
def stars(q): return "**" if q < 0.01 else "*" if q < 0.05 else ""

fig = plt.figure(figsize=(7.1, 6.6))
gs = fig.add_gridspec(2, 3, hspace=0.78, wspace=0.62, left=0.1, right=0.985, top=0.94, bottom=0.11, width_ratios=[1, 1, 1.05])

# a - leave-one-out bias: raw vs permutation-null medians
ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", "Leave-one-out bias", x=-0.5)
yy = np.arange(len(ORDER))
for k, (src, mk, off) in enumerate([("PRISM", "o", -0.17), ("GDSC", "s", 0.17)]):
    s = ss[ss.source == src].set_index("set").loc[ORDER]
    for i, st in enumerate(ORDER):
        ax.plot([s.median_null_rho[st], s.median_raw_rho[st]], [i + off] * 2, color=LIGHTGRAY, lw=1.0, zorder=1)
    ax.scatter(s.median_null_rho, yy + off, marker=mk, s=12, color=GRAY, zorder=3, label=f"{src}: null" if True else None)
    ax.scatter(s.median_raw_rho, yy + off, marker=mk, s=12, color=[SET_COLORS[x] for x in ORDER], edgecolor=INK, lw=0.3, zorder=4, label=f"{src}: observed")
ax.axvline(0, color=INK2, lw=0.6, ls="--")
ax.set_yticks(yy); ax.set_yticklabels([SHORT[x] for x in ORDER], fontsize=5.8); ax.invert_yaxis()
ax.set_xlabel("Median per-drug Spearman ρ\n(observed vs label-permuted)"); ax.set_xlim(-0.55, 0.12)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=INK2, ms=3.5, label="PRISM"), Line2D([], [], marker="s", ls="", color=INK2, ms=3.5, label="GDSC"),
                   Line2D([], [], marker="o", ls="", color=GRAY, ms=3.5, label="Permutation null"), Line2D([], [], marker="o", ls="", mfc=BLUE, mec=INK, ms=3.5, label="Observed")],
          fontsize=5.2, loc="upper center", bbox_to_anchor=(0.45, -0.3), ncol=2)

# b, c - calibrated per-drug rho per set
for letter, src, pos in [("b", "PRISM", gs[0, 1]), ("c", "GDSC", gs[0, 2])]:
    ax = fig.add_subplot(pos); n = ss[ss.source == src].n_drugs.iloc[0]
    panel(ax, letter, f"Null-calibrated ρ, {src} ({n} drugs)", x=-0.3)
    s = ss[ss.source == src].set_index("set")
    for i, st in enumerate(ORDER):
        r = pdc[(pdc.source == src) & (pdc.set == st)].calibrated_rho
        ax.scatter(i + rng.uniform(-0.2, 0.2, len(r)), r, s=4, color=SET_COLORS[st], alpha=0.65, lw=0, zorder=2)
        ax.hlines(r.median(), i - 0.34, i + 0.34, color=INK, lw=1.2, zorder=4)
        q = s.loc[st, "fdr_bh"]
        ax.text(i, 0.86, f"{r.median():+.2f}{stars(q)}", ha="center", fontsize=5.4, color=INK, fontweight="bold" if q < 0.05 else "normal")
    ax.axhline(0, color=GRAY, lw=0.7, ls="--")
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([SHORT[x] for x in ORDER], rotation=45, ha="right", fontsize=5.6)
    ax.set_ylim(-0.85, 0.95); ax.set_ylabel("Calibrated ρ (observed − null)")
    if src == "GDSC":
        q = sp[(sp.source == "GDSC") & (sp.comparison == "RPPA vs RNA")].iloc[0]
        ax.annotate("", xy=(2, 0.66), xytext=(0, 0.66), arrowprops=dict(arrowstyle="-", lw=0.6, color=INK2))
        ax.text(1, 0.68, f"RPPA > RNA\nFDR {q.fdr_bh:.3f}", ha="center", fontsize=5.2, color=INK2, va="bottom")
        ax.set_ylim(-0.85, 1.05)

# d - robustness heatmap
ax = fig.add_subplot(gs[1, 0:2]); panel(ax, "d", "Robustness across model settings (median calibrated ρ; * p < 0.05)", x=-0.2, y=1.1); ax.set_title(ax.get_title(loc="left"), pad=13, loc="left")
CFG = [("k5_a1.0 (primary, 500 perm)", "Primary\nk = 5, α = 1"), ("k3_a1.0 (200 perm)", "k = 3"), ("k10_a1.0 (200 perm)", "k = 10"),
       ("k5_a0.1 (200 perm)", "α = 0.1"), ("k5_a10.0 (200 perm)", "α = 10")]
rows = [(src, st) for src in ["PRISM", "GDSC"] for st in ORDER]
M = np.array([[rob[(rob.source == s_) & (rob.set == st) & (rob.config == c)].median_calibrated_rho.iloc[0] for c, _ in CFG] for s_, st in rows])
P = np.array([[rob[(rob.source == s_) & (rob.set == st) & (rob.config == c)].p_calibrated_gt0.iloc[0] for c, _ in CFG] for s_, st in rows])
Mt, Pt = M.T, P.T
im = ax.imshow(Mt, cmap="RdBu_r", vmin=-0.2, vmax=0.2, aspect="auto")
for i in range(Mt.shape[0]):
    for j in range(Mt.shape[1]):
        ax.text(j, i, f"{Mt[i,j]:+.2f}" + ("*" if Pt[i, j] < 0.05 else ""), ha="center", va="center", fontsize=4.7,
                color="white" if abs(Mt[i, j]) > 0.14 else INK, linespacing=0.8)
ax.set_yticks(range(len(CFG))); ax.set_yticklabels([c[1].replace("\n", " ") for c in CFG], fontsize=5.6)
ax.set_xticks(range(len(rows))); ax.set_xticklabels([SHORT[st] for _, st in rows], rotation=45, ha="right", fontsize=5.6)
ax.axvline(5.5, color="white", lw=2)
for k, src in enumerate(["PRISM", "GDSC"]): ax.text(2.5 + 6 * k, -0.62, src, ha="center", fontsize=6.3, fontweight="bold")
for sp_ in ax.spines.values(): sp_.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("Median calibrated ρ", fontsize=5.8); cb.ax.tick_params(labelsize=5.4)

# e - RPPA vs RNA per drug (GDSC + PRISM)
ax = fig.add_subplot(gs[1, 2]); panel(ax, "e", "Per drug: RPPA vs RNA", x=-0.3)
w = pdc.pivot_table(index=["source", "drug_key"], columns="set", values="calibrated_rho").reset_index()
for src, col, mk in [("PRISM", BLUE, "o"), ("GDSC", ORANGE, "s")]:
    g = w[w.source == src]
    ax.scatter(g.RNA, g.RPPA, s=10, marker=mk, color=col, edgecolor="white", lw=0.3, alpha=0.9, zorder=3,
               label=f"{src}: RPPA > RNA in {(g.RPPA > g.RNA).sum()}/{len(g)}")
lim = (-0.75, 0.85); ax.plot(lim, lim, color=GRAY, lw=0.7, ls="--"); ax.set_xlim(lim); ax.set_ylim(lim)
ax.axhline(0, color=LIGHTGRAY, lw=0.5); ax.axvline(0, color=LIGHTGRAY, lw=0.5)
ax.set_xlabel("Calibrated ρ, baseline RNA"); ax.set_ylabel("Calibrated ρ, baseline RPPA"); ax.set_aspect("equal")
ax.legend(fontsize=5.2, loc="upper center", bbox_to_anchor=(0.45, -0.3))
w.to_csv(OUT + "source_Fig4e.csv", index=False)
save(fig, "Figure4_protein_vs_rna")

# cross-screen reproducibility of per-drug calibrated rho (reported in text / Supplementary)
x = pdc.pivot_table(index=["drug_key", "set"], columns="source", values="calibrated_rho").dropna()
out = []
for st in ORDER:
    g = x.xs(st, level="set"); r, p = spearmanr(g.PRISM, g.GDSC); out.append(dict(set=st, n_drugs=len(g), rho=r, p=p))
pd.DataFrame(out).to_csv(OUT + "source_crossscreen_C3.csv", index=False); print(pd.DataFrame(out))

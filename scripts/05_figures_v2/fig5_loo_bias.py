"""v2 Figure 5 - Leave-one-out bias of within-drug evaluation and its permutation calibration (T1-B + Amendment 1)."""
import sys; sys.path.insert(0, "/home/claude/v2/figs")
from figcommon import *
from matplotlib.lines import Line2D
style()
R = "/home/claude/v2/results/T1B/"; A1 = "/home/claude/mirror/results/confirmatory/amendment1/"
G = pd.read_csv(R + "null_grid_summary.csv"); C = pd.read_csv(R + "calibration_summary.csv")
fig = plt.figure(figsize=(7.1, 5.4))
gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.55, left=0.09, right=0.985, top=0.93, bottom=0.12)

# a - null LOO rho vs n, by number of components (alpha = 1)
ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", "Null bias vs lines per drug", x=-0.3)
cols = {1: AQUA, 3: YELLOW, 5: BLUE, 10: MAGENTA}
for k, c in cols.items():
    g = G[(G.k == k) & (G.alpha == 1)].sort_values("n")
    ax.plot(g.n, g["mean"], "-o", color=c, ms=3, lw=1.1, label=f"k = {k} PCs")
g = G[G.k == 0].sort_values("n"); ax.plot(g.n, g["mean"], "--", color=INK2, lw=1, label="intercept only")
ax.axhline(0, color=GRAY, lw=0.6, ls=":"); ax.axvline(27, color=LIGHTGRAY, lw=0.8); ax.text(28, -0.22, "median\ndrug (27)", fontsize=5, color=INK2, va="top")
ax.set_xlabel("Cell lines per drug (n)"); ax.set_ylabel("Mean LOO ρ under the null"); ax.set_ylim(-1.05, 0.1)
ax.legend(fontsize=5.0, loc="upper right", ncol=2, columnspacing=0.6, handlelength=1.4)

# b - effect of ridge penalty (k = 5)
ax = fig.add_subplot(gs[0, 1]); panel(ax, "b", "Null bias vs ridge penalty (k = 5)", x=-0.3)
for a_, c in zip([0.1, 1, 10, 100], [AQUA, BLUE, ORANGE, MAGENTA]):
    g = G[(G.k == 5) & (G.alpha == a_)].sort_values("n")
    ax.plot(g.n, g["mean"], "-o", color=c, ms=3, lw=1.1, label=f"α = {a_:g}")
ax.axhline(0, color=GRAY, lw=0.6, ls=":"); ax.set_xlabel("Cell lines per drug (n)"); ax.set_ylabel("Mean LOO ρ under the null")
ax.set_ylim(-1.05, 0.1); ax.legend(fontsize=5.0, loc="upper left", ncol=2, columnspacing=0.6, handlelength=1.4)

# c - empirical: observed vs permutation medians in the C3 data
ax = fig.add_subplot(gs[0, 2]); panel(ax, "c", "Same bias in the real C3 data", x=-0.42)
s = pd.read_csv(A1 + "summary_sets.csv")
ORDER = ["RNA", "Protein_MS", "RPPA", "Response", "Response+RNA", "Response+Protein_MS"]
SH = {"RNA": "RNA", "Protein_MS": "Protein (MS)", "RPPA": "RPPA", "Response": "Tahoe response", "Response+RNA": "Resp. + RNA", "Response+Protein_MS": "Resp. + MS"}
for k, (src, mk, off) in enumerate([("PRISM", "o", -0.17), ("GDSC", "s", 0.17)]):
    g = s[s.source == src].set_index("set").loc[ORDER]
    for i, st in enumerate(ORDER):
        ax.plot([g.median_null_rho[st], g.median_raw_rho[st]], [i + off] * 2, color=LIGHTGRAY, lw=1)
    ax.scatter(g.median_null_rho, np.arange(6) + off, marker=mk, s=11, color=GRAY, zorder=3)
    ax.scatter(g.median_raw_rho, np.arange(6) + off, marker=mk, s=11, color=[SET_COLORS[x] for x in ORDER], edgecolor=INK, lw=0.3, zorder=4)
ax.axvline(0, color=INK2, lw=0.6, ls="--"); ax.set_yticks(range(6)); ax.set_yticklabels([SH[x] for x in ORDER], fontsize=5.7); ax.invert_yaxis()
ax.set_xlabel("Median per-drug ρ"); ax.set_xlim(-0.55, 0.12)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=GRAY, ms=3.5, label="permutation null"), Line2D([], [], marker="o", ls="", mfc=BLUE, mec=INK, ms=3.5, label="observed"),
                   Line2D([], [], marker="o", ls="", color=INK2, ms=3, label="PRISM"), Line2D([], [], marker="s", ls="", color=INK2, ms=3, label="GDSC")],
          fontsize=5.0, loc="upper center", bbox_to_anchor=(0.4, -0.25), ncol=2, columnspacing=0.5, handletextpad=0.2)

# d, e - calibration under null and with signal (LOO and 5-fold)
for letter, scheme, pos, ttl in [("d", "loo", gs[1, 0], "Calibration, leave-one-out"), ("e", "5fold", gs[1, 1], "Calibration, 5-fold CV")]:
    ax = fig.add_subplot(pos); panel(ax, letter, ttl, x=-0.3)
    for tr, c in zip([0.0, 0.2, 0.4], [GRAY, BLUE, ORANGE]):
        g = C[(C.scheme == scheme) & (C.true_r == tr)].sort_values("n")
        ax.plot(g.n, g.raw, ":", color=c, lw=1)
        ax.plot(g.n, g.calibrated, "-o", color=c, ms=3, lw=1.2)
        ax.plot(g.n, g.oracle, "--", color=c, lw=0.8, alpha=0.7)
    ax.axhline(0, color=LIGHTGRAY, lw=0.6); ax.set_ylim(-0.9, 0.65)
    ax.set_xlabel("Cell lines per drug (n)"); ax.set_ylabel("Mean ρ")
    if letter == "d":
        ax.legend(handles=[Line2D([], [], color=c, lw=1.4, label=f"true r = {t}") for t, c in zip([0, 0.2, 0.4], [GRAY, BLUE, ORANGE])] +
                  [Line2D([], [], color=INK2, ls=":", label="raw"), Line2D([], [], color=INK2, marker="o", ms=3, label="calibrated"), Line2D([], [], color=INK2, ls="--", label="true (sample)")],
                  fontsize=5.0, loc="upper center", bbox_to_anchor=(1.15, -0.27), ncol=6, columnspacing=0.8, handletextpad=0.3)

# f - summary text panel: practical recommendation (as a small table)
ax = fig.add_subplot(gs[1, 2]); ax.set_axis_off(); panel(ax, "f", "Take-home", x=-0.05)
lines = ["• LOO with a training-mean intercept is", "  negatively biased under the null", "  (−0.77 at 27 lines; −1 without features).",
         "• The bias depends on model size, so", "  raw ρ cannot compare feature sets.", "• Permutation calibration centres the null", "  at 0 → valid test statistic.",
         "• Calibrated ρ is not a correlation:", "  LOO overshoots at large n; K-fold is closer.", "• Report calibrated ρ for testing; avoid", "  reading its size as an effect size."]
for i, t in enumerate(lines): ax.text(0.0, 0.95 - i * 0.088, t, transform=ax.transAxes, fontsize=5.7, va="top", color=INK)
save(fig, "Figure5_loo_bias")

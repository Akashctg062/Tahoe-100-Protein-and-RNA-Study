"""v2 Figure 6 - Baseline cell state, not the 24-h Tahoe response, predicts within-drug sensitivity (3 screens)."""
import sys; sys.path.insert(0, "/home/claude/v2/figs")
from figcommon import *
from matplotlib.lines import Line2D
style()
R = "/home/claude/v2/results/"; A1 = "/home/claude/mirror/results/confirmatory/amendment1/"
a1 = pd.read_csv(A1 + "summary_sets.csv").rename(columns={"median_calibrated_rho": "med", "fdr_bh": "fdr"})
a1["set"] = a1.set.replace({"Protein_MS": "MS protein (CCLE)"})
e = pd.read_csv(R + "T2EF/summary_sets_E_ctrp.csv").rename(columns={"median_calibrated": "med", "fdr_bh": "fdr"})
e["set"] = e.set.replace({"Protein": "MS protein (CCLE)", "Response+Protein": "Response+Protein_MS"})
f = pd.read_csv(R + "T2EF/summary_sets_F_procan.csv").rename(columns={"median_calibrated": "med", "fdr_bh": "fdr"})
fp = f[f.set == "Protein"].assign(set="MS protein (ProCan)")
S = pd.concat([a1[["source", "set", "med", "fdr", "n_drugs"]], e[["source", "set", "med", "fdr", "n_drugs"]], fp[["source", "set", "med", "fdr", "n_drugs"]]])
ROWS = ["RNA", "MS protein (CCLE)", "MS protein (ProCan)", "RPPA", "Response", "Response+RNA"]
LAB = {"RNA": "Baseline RNA", "MS protein (CCLE)": "Protein, MS (CCLE)", "MS protein (ProCan)": "Protein, MS (ProCan)",
       "RPPA": "Protein, RPPA", "Response": "Tahoe response", "Response+RNA": "Response + RNA"}
COLS = ["PRISM", "GDSC", "CTRP"]

fig = plt.figure(figsize=(7.1, 6.9))
gs = fig.add_gridspec(2, 3, hspace=0.6, wspace=0.62, left=0.13, right=0.985, top=0.94, bottom=0.12, width_ratios=[1.15, 1, 1])

# a - heatmap of median calibrated rho
ax = fig.add_subplot(gs[0, 0]); panel(ax, "a", "Median calibrated ρ", x=-0.55)
M = np.full((len(ROWS), 3), np.nan); Q = np.full_like(M, np.nan)
for i, r in enumerate(ROWS):
    for j, c in enumerate(COLS):
        x = S[(S.set == r) & (S.source == c)]
        if len(x): M[i, j], Q[i, j] = x.med.iloc[0], x.fdr.iloc[0]
im = ax.imshow(M, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto")
for i in range(len(ROWS)):
    for j in range(3):
        if np.isfinite(M[i, j]):
            ax.text(j, i, f"{M[i,j]:+.2f}" + ("*" if Q[i, j] < 0.05 else ""), ha="center", va="center", fontsize=5.6,
                    color="white" if abs(M[i, j]) > 0.2 else INK, fontweight="bold" if Q[i, j] < 0.05 else "normal")
ax.set_xticks(range(3)); ax.set_xticklabels(COLS, fontsize=6.2); ax.xaxis.tick_top()
ax.set_yticks(range(len(ROWS))); ax.set_yticklabels([LAB[r] for r in ROWS], fontsize=5.9)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.tick_params(length=0); ax.axhline(3.5, color="white", lw=2)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03); cb.ax.tick_params(labelsize=5.2); cb.set_label("calibrated ρ", fontsize=5.6)
ax.text(0.5, -0.1, "* FDR < 0.05 (Wilcoxon > 0)", transform=ax.transAxes, ha="center", va="top", fontsize=5.2, color=INK2)

# b - CTRP per-drug calibrated rho by set
ax = fig.add_subplot(gs[0, 1:]); panel(ax, "b", "Independent replication in CTRP (46 drugs)", x=-0.09)
pe = pd.read_csv(R + "T2EF/per_drug_E_ctrp.csv")
ORDB = ["RNA", "Protein", "RPPA", "Response", "Response+RNA", "Response+Protein"]
LB = {"RNA": "RNA", "Protein": "Protein (MS)", "RPPA": "RPPA", "Response": "Tahoe response", "Response+RNA": "Response + RNA", "Response+Protein": "Response + MS"}
COLB = {"RNA": ORANGE, "Protein": AQUA, "RPPA": YELLOW, "Response": BLUE, "Response+RNA": MAGENTA, "Response+Protein": GREEN}
rng = np.random.default_rng(3); es = pd.read_csv(R + "T2EF/summary_sets_E_ctrp.csv").set_index("set")
for i, s in enumerate(ORDB):
    v = pe[pe.set == s].calibrated_rho
    ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=6, color=COLB[s], alpha=0.7, lw=0)
    ax.hlines(v.median(), i - 0.33, i + 0.33, color=INK, lw=1.3, zorder=4)
    ax.text(i, 1.02, f"{v.median():+.2f}" + ("*" if es.loc[s, "fdr_bh"] < 0.05 else ""), ha="center", fontsize=5.8, fontweight="bold")
ax.axhline(0, color=GRAY, lw=0.7, ls="--"); ax.set_ylim(-0.9, 1.12)
ax.set_xticks(range(len(ORDB))); ax.set_xticklabels([LB[s] for s in ORDB], rotation=30, ha="right", fontsize=5.9)
ax.set_ylabel("Calibrated ρ per drug")
ep = pd.read_csv(R + "T2EF/summary_pairs_E_ctrp.csv").set_index("comparison")
ax.text(0.99, 0.04, f"Response < RNA: FDR {ep.loc['Response vs RNA','fdr_bh']:.3f}\nRPPA < RNA: FDR {ep.loc['RPPA vs RNA','fdr_bh']:.4f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4, color=INK2)

# c - pairwise differences across screens
ax = fig.add_subplot(gs[1, 0]); panel(ax, "c", "Paired differences vs baseline RNA", x=-0.55)
pa = pd.read_csv(A1 + "summary_pairs.csv").replace({"Protein_MS": "Protein"}, regex=True)
pe2 = pd.read_csv(R + "T2EF/summary_pairs_E_ctrp.csv"); pf = pd.read_csv(R + "T2EF/summary_pairs_F_procan.csv")
rowsc = [("Response vs RNA", "Tahoe response − RNA", "A"), ("RPPA vs RNA", "RPPA − RNA", "A"), ("Protein vs RNA", "MS (CCLE) − RNA", "A"),
         ("Protein vs RNA", "MS (ProCan) − RNA", "F"), ("Response+RNA vs RNA", "Response + RNA − RNA", "A")]
mk = {"PRISM": ("o", BLUE), "GDSC": ("s", ORANGE), "CTRP": ("D", AQUA)}
for i, (comp, lab, src) in enumerate(rowsc):
    for k, scr in enumerate(COLS):
        if src == "F": d = pf[(pf.source == scr) & (pf.comparison == comp)]
        elif scr == "CTRP": d = pe2[pe2.comparison == comp.replace("Protein", "Protein")]
        else: d = pa[(pa.source == scr) & (pa.comparison == comp)]
        if not len(d): continue
        d = d.iloc[0]; y = i + (k - 1) * 0.22
        ax.scatter(d.median_difference, y, marker=mk[scr][0], s=16, color=mk[scr][1], edgecolor=INK if d.fdr_bh < 0.05 else "none", lw=0.7, zorder=3)
ax.axvline(0, color=GRAY, lw=0.7, ls="--")
ax.set_yticks(range(len(rowsc))); ax.set_yticklabels([r[1] for r in rowsc], fontsize=5.9); ax.invert_yaxis()
ax.set_xlabel("Median paired difference\nin calibrated ρ"); ax.set_xlim(-0.3, 0.3)
ax.legend(handles=[Line2D([], [], marker=mk[s][0], ls="", color=mk[s][1], ms=4, label=s) for s in COLS] +
          [Line2D([], [], marker="o", ls="", mfc="white", mec=INK, ms=4, label="FDR < 0.05")],
          fontsize=5.2, loc="upper center", bbox_to_anchor=(0.4, -0.32), ncol=4, columnspacing=0.6, handletextpad=0.2)

# d - model robustness (4 best-effort families + frozen ridge LOO)
ax = fig.add_subplot(gs[1, 1]); panel(ax, "d", "Robust to model choice", x=-0.32)
P4 = pd.concat([pd.read_csv(R + f"T1A2/summary_pairs_{m}.csv") for m in ["ridge_krr", "enet", "rf"]])
FAM = [("primary", "Ridge LOO\n(frozen)"), ("ridge", "Ridge-GCV"), ("krr", "Kernel ridge"), ("enet", "Elastic net"), ("rf", "Random\nforest")]
comps = [("Response+RNA vs RNA", "Response adds to RNA", MAGENTA), ("RPPA vs RNA", "RPPA − RNA", YELLOW)]
for j, (fam, fl) in enumerate(FAM):
    for k, (comp, cl, col) in enumerate(comps):
        for scr, m_ in [("PRISM", "o"), ("GDSC", "s")]:
            if fam == "primary":
                d = pa[(pa.source == scr) & (pa.comparison == comp)]
            else:
                d = P4[(P4.model == fam) & (P4.source == scr) & (P4.comparison == comp.replace("Protein", "Protein_MS"))]
            if not len(d): continue
            d = d.iloc[0]; x = j + (k - 0.5) * 0.36 + (0.08 if scr == "GDSC" else -0.08)
            ax.scatter(x, d.median_difference, marker=m_, s=13, color=col, edgecolor=INK if d.fdr_bh < 0.05 else "none", lw=0.7, zorder=3)
ax.axhline(0, color=GRAY, lw=0.7, ls="--")
ax.set_xticks(range(len(FAM))); ax.set_xticklabels([f[1] for f in FAM], fontsize=5.4, rotation=35, ha="right")
ax.set_ylabel("Median paired difference"); ax.set_ylim(-0.25, 0.25)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=c[2], ms=4, label=c[1]) for c in comps] +
          [Line2D([], [], marker="o", ls="", color=INK2, ms=3.5, label="PRISM"), Line2D([], [], marker="s", ls="", color=INK2, ms=3.5, label="GDSC")],
          fontsize=5.1, loc="upper center", bbox_to_anchor=(0.45, -0.42), ncol=2, columnspacing=0.6, handletextpad=0.2)

# e - top vs middle dose response
ax = fig.add_subplot(gs[1, 2]); panel(ax, "e", "Tahoe response: top vs middle dose", x=-0.32)
top = pd.read_csv(A1 + "per_drug_calibrated.csv"); mid = pd.read_csv(R + "T2EF/per_drug_T1D_middose.csv")
w = top[top.set == "Response"].merge(mid[mid.set == "Response"], on=["source", "drug_key"], suffixes=("_top", "_mid"))
for scr, m_, col in [("PRISM", "o", BLUE), ("GDSC", "s", ORANGE)]:
    g = w[w.source == scr]
    ax.scatter(g.calibrated_rho_top, g.calibrated_rho_mid, s=9, marker=m_, color=col, edgecolor="white", lw=0.3, alpha=0.9,
               label=f"{scr}: r = {np.corrcoef(g.calibrated_rho_top, g.calibrated_rho_mid)[0,1]:.2f}")
lim = (-0.8, 0.95); ax.plot(lim, lim, color=GRAY, lw=0.7, ls="--"); ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
ax.axhline(0, color=LIGHTGRAY, lw=0.5); ax.axvline(0, color=LIGHTGRAY, lw=0.5)
ax.set_xlabel("Calibrated ρ, top dose"); ax.set_ylabel("Calibrated ρ, middle dose")
ax.legend(fontsize=5.2, loc="upper center", bbox_to_anchor=(0.45, -0.32))
w.to_csv(OUT + "source_Fig6e.csv", index=False); S.to_csv(OUT + "source_Fig6a.csv", index=False)
save(fig, "Figure6_baseline_vs_response")

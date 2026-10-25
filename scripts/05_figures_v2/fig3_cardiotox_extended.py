"""v2 Figure 3 - Best-effort cardiotoxicity models, hERG potency, and cross-context drug identity (T1-A1, T2-H2, T2-H1)."""
import sys; sys.path.insert(0, "/home/claude/v2/figs")
from figcommon import *
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
style()
R = "/home/claude/v2/results/"
CHEM, TAHOE, COMB = ORANGE, BLUE, AQUA
FAM = [("enet_logistic", "Elastic-net\nlogistic"), ("grad_boost", "Gradient\nboosting"), ("rbf_svm", "RBF SVM"), ("mlp", "MLP")]
fig = plt.figure(figsize=(7.1, 6.3))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.62, left=0.1, right=0.97, top=0.94, bottom=0.12, width_ratios=[1, 1, 1.25])

def auroc_panel(ax, task, letter, title, show_legend):
    s = pd.read_csv(R + f"T1A1/summary_{task}.csv")
    panel(ax, letter, title, x=-0.42)
    for i, (f, fl) in enumerate(FAM):
        g = s[s.family == f].set_index("representation")
        ax.plot([g.auroc_chem.iloc[0]], [i - 0.24], "D", color=CHEM, ms=3.6)
        ax.plot([g.auroc_tahoe["pca50"]], [i - 0.08], "o", color=TAHOE, ms=3.6)
        ax.plot([g.auroc_tahoe["programs"]], [i + 0.08], "s", color=TAHOE, ms=3.4, mfc="white", mew=0.9)
        ax.plot([g.auroc_chem_tahoe["pca50"]], [i + 0.24], "o", color=COMB, ms=3.6)
        ax.plot([g.auroc_chem_tahoe["programs"]], [i + 0.40], "s", color=COMB, ms=3.4, mfc="white", mew=0.9)
    ax.axvline(0.5, color=GRAY, lw=0.7, ls="--")
    ax.set_yticks(np.arange(4) + 0.08); ax.set_yticklabels([x[1] for x in FAM], fontsize=5.8); ax.invert_yaxis()
    ax.set_xlim(0.2, 0.8); ax.set_xlabel("AUROC (nested CV, out-of-fold)")
    if show_legend:
        ax.legend(handles=[Line2D([], [], marker="D", ls="", color=CHEM, ms=3.5, label="Chemistry"),
                           Line2D([], [], marker="o", ls="", color=TAHOE, ms=3.5, label="Tahoe PCA-50"),
                           Line2D([], [], marker="s", ls="", color=TAHOE, mfc="white", ms=3.5, label="Tahoe programs"),
                           Line2D([], [], marker="o", ls="", color=COMB, ms=3.5, label="Chem. + PCA-50"),
                           Line2D([], [], marker="s", ls="", color=COMB, mfc="white", ms=3.5, label="Chem. + programs")],
                  fontsize=5.0, loc="upper center", bbox_to_anchor=(1.2, -0.2), ncol=5, columnspacing=0.8, handletextpad=0.2)

auroc_panel(fig.add_subplot(gs[0, 0]), "D1", "a", "DICTrank (best effort)", True)
auroc_panel(fig.add_subplot(gs[0, 1]), "C1", "b", "FAERS (best effort)", False)

# c - delta AUROC forest (chem + Tahoe - chem), 95% CI, Holm-adjusted p
ax = fig.add_subplot(gs[0, 2]); panel(ax, "c", "Added value of Tahoe (Δ AUROC)", x=-0.5)
rows = []
for task, tn in [("C1", "FAERS"), ("D1", "DICTrank")]:
    s = pd.read_csv(R + f"T1A1/summary_{task}.csv")
    for f, fl in FAM:
        for rep, rl in [("pca50", "PCA-50"), ("programs", "programs")]:
            r = s[(s.family == f) & (s.representation == rep)].iloc[0]
            rows.append(dict(task=tn, lab=f"{fl.replace(chr(10), ' ')}, {rl}", d=r.delta, lo=r.delta_ci_lo, hi=r.delta_ci_hi, ph=r.p_holm))
rows = pd.DataFrame(rows); y = np.arange(len(rows)) + np.where(rows.task == "DICTrank", 0.8, 0)
for yy, r in zip(y, rows.itertuples()):
    c = COMB if r.task == "FAERS" else MAGENTA
    ax.plot([r.lo, r.hi], [yy, yy], color=c, lw=0.9); ax.plot(r.d, yy, "o", color=c, ms=2.8)
ax.axvline(0, color=INK2, lw=0.6, ls="--")
ax.set_yticks(y); ax.set_yticklabels(rows.lab, fontsize=4.6); ax.invert_yaxis(); ax.set_xlim(-0.55, 0.3)
ax.set_xlabel("Δ AUROC, chem. + Tahoe − chem. (95% CI)")
ax.text(0.27, y[0] - 0.4, "FAERS (C1, n = 128)", fontsize=5.3, color=COMB, ha="right", fontweight="bold", va="bottom")
ax.text(0.27, y[8] - 0.4, "DICTrank (D1)", fontsize=5.3, color=MAGENTA, ha="right", fontweight="bold", va="bottom")
ax.text(-0.53, (y[7] + y[8]) / 2, f"All Holm-adjusted p ≥ {rows.ph.min():.2f}", fontsize=5.0, color=INK2, va="center"); ax.set_ylim(y[-1] + 0.7, -1.4)

# d - hERG: observed pIC50 vs out-of-fold predictions
ax = fig.add_subplot(gs[1, 0]); panel(ax, "d", "hERG potency (ChEMBL, 43 drugs)", x=-0.36)
h = pd.read_csv(R + "T2H2/oof_regression.csv"); hs = pd.read_csv(R + "T2H2/summary.csv").set_index("task")
ax.scatter(h.oof_chem, h.pIC50, s=10, color=CHEM, edgecolor="white", lw=0.3, label=f"Chemistry  ρ = {spear(h.oof_chem, h.pIC50):+.2f}", zorder=3)
ax.scatter(h.oof_chem_tahoe, h.pIC50, s=16, marker="s", facecolor="none", edgecolor=COMB, lw=0.7, label=f"Chem. + Tahoe  ρ = {spear(h.oof_chem_tahoe, h.pIC50):+.2f}", zorder=3)
ax.axhline(5, color=LIGHTGRAY, lw=0.7, ls=":")
ax.set_xlabel("Predicted pIC50 (out-of-fold)"); ax.set_ylabel("Measured hERG pIC50")
r = hs.loc["regression_pIC50_spearman"]
ax.legend(fontsize=5.0, loc="upper left", handletextpad=0.2, borderaxespad=0.2)
ax.text(0.03, 0.80, f"Δρ = {r.delta:+.3f}\n({r.ci_lo:+.3f} to {r.ci_hi:+.3f})", transform=ax.transAxes, fontsize=5.2, ha="left", va="top", color=INK2)

# e - same vs different drug correlation (tumour vs cardiomyocyte)
M = pd.read_csv(R + "T2H1/corr_matrix_averaged_DEGs.csv", index_col=0)
Rm = M.values; n = len(M); same = np.diag(Rm); diff = Rm[~np.eye(n, dtype=bool)]
ax = fig.add_subplot(gs[1, 1]); panel(ax, "e", "Drug identity across contexts", x=-0.36)
rng = np.random.default_rng(7)
parts = ax.violinplot([diff], positions=[0], widths=0.7, showextrema=False)
for b in parts["bodies"]: b.set_facecolor(LIGHTGRAY); b.set_alpha(0.8); b.set_edgecolor("none")
ax.scatter(1 + rng.uniform(-0.15, 0.15, n), same, s=10, color=AQUA, edgecolor="white", lw=0.3, zorder=3)
ax.hlines([np.mean(diff)], -0.3, 0.3, color=INK, lw=1); ax.hlines([np.mean(same)], 0.7, 1.3, color=INK, lw=1)
ax.axhline(0, color=GRAY, lw=0.6, ls=":")
ax.set_xticks([0, 1]); ax.set_xticklabels([f"Different drugs\n(n = {len(diff)})", f"Same drug\n(n = {n})"], fontsize=5.8)
ax.set_ylabel("Pearson r, tumour vs cardiomyocyte\nsignature (5,518 genes)"); ax.set_xlim(-0.6, 1.6)
prev = 9
for k in np.argsort(-same)[:3]:
    yy = min(same[k], prev - 0.026); prev = yy
    ax.text(1.2, yy, M.index[k], fontsize=4.8, color=INK2, va="center")
null = np.array([Rm[np.arange(n), rng.permutation(n)].mean() for _ in range(10000)])
p = (np.sum(null >= same.mean()) + 1) / 10001
ins = ax.inset_axes([0.05, 0.72, 0.4, 0.24]); ins.hist(null, bins=40, color=LIGHTGRAY); ins.axvline(same.mean(), color=AQUA, lw=1.2)
ins.set_yticks([]); ins.tick_params(labelsize=4.5, length=1.5); ins.set_title(f"permutation p = {p:.4f}", fontsize=4.8, pad=1.5)
for sp_ in ["top", "right", "left"]: ins.spines[sp_].set_visible(False)
ax.set_ylim(min(diff.min(), same.min()) - 0.02, max(diff.max(), same.max()) + 0.12)

# f - heatmap ordered by drug class
CLS = [("Anthracycline", ["doxorubicin", "epirubicin", "idarubicin"]),
       ("EGFR/HER2", ["afatinib", "erlotinib", "gefitinib", "lapatinib", "vandetanib"]),
       ("Multi-kinase", ["cabozantinib", "ponatinib", "regorafenib", "sunitinib", "crizotinib"]),
       ("BRAF/MEK", ["dabrafenib", "vemurafenib", "trametinib"]),
       ("Hypomethylating", ["azacitidine", "decitabine"]),
       ("Other", ["bortezomib", "cyclosporine", "tofacitinib", "pioglitazone", "rosiglitazone", "phenylephrine", "verapamil"])]
order = [d for _, ds in CLS for d in ds]; assert sorted(order) == sorted(M.index), set(M.index) ^ set(order)
Mo = M.loc[order, order]
ax = fig.add_subplot(gs[1, 2]); panel(ax, "f", "Tumour (rows) vs cardiomyocyte (columns)", x=-0.38)
v = np.nanpercentile(np.abs(Mo.values), 98)
im = ax.imshow(Mo.values, cmap="RdBu_r", norm=TwoSlopeNorm(0, -v, v), aspect="equal")
for i in range(n): ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, ec=INK, lw=0.5))
b = 0
for _, ds in CLS[:-1]:
    b += len(ds); ax.axhline(b - 0.5, color="white", lw=0.8); ax.axvline(b - 0.5, color="white", lw=0.8)
ax.set_yticks(range(n)); ax.set_yticklabels(order, fontsize=4.2); ax.set_xticks(range(n)); ax.set_xticklabels(order, fontsize=4.2, rotation=90)
ax.tick_params(length=1)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=5); cb.set_label("Pearson r", fontsize=5.5)
pd.DataFrame(dict(drug=M.index, r_same=same)).to_csv(OUT + "source_Fig3e.csv", index=False)
save(fig, "Figure3_cardiotox_extended")

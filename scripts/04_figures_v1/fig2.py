"""Figure 2 - Tumour-cell signatures do not predict clinical cardiotoxicity (D1, C1)."""
import sys, re; sys.path.insert(0, "/home/claude/figs")
from figcommon import *
from sklearn.metrics import roc_auc_score, roc_curve
style()
CHEM, TAHOE, COMB = ORANGE, BLUE, AQUA

# ---- D1 results parsed from the saved discovery outputs (logistic models) ----
def parse(path):
    rows, task = [], None
    for line in open(path):
        m = re.match(r"== (\w+) \| (\w+): n=(\d+) pos=(\d+)", line)
        if m: task = (m.group(2), int(m.group(3)), int(m.group(4))); continue
        m = re.match(r"\s+(\S+)\s+(logit|rf)\s+AUROC ([\d.]+) ± ([\d.]+)", line)
        if m and task: rows.append(dict(task=task[0], n=task[1], pos=task[2], feat=m.group(1), model=m.group(2),
                                        auc=float(m.group(3)), sd=float(m.group(4))))
    return pd.DataFrame(rows)
d1 = parse(ANA + "out_all.txt").assign(scope="all"); d1l = parse(ANA + "out_lung.txt").assign(scope="lung")
FEATS = [("structure", "Chemistry (structure)", CHEM), ("bio_raw2000", "Tahoe: 2,000 genes", TAHOE),
         ("bio_PCA10_alldoses", "Tahoe: PCA-10, 3 doses", TAHOE), ("bio_programs", "Tahoe: 14 gene programs", TAHOE),
         ("bio_programs+potency", "Tahoe: programs + potency", TAHOE), ("structure+programs", "Chemistry + Tahoe", COMB)]
TASKS = [("most_vs_no", "Most vs no concern"), ("mostless_vs_no", "Most+less vs no"),
         ("onc_most_vs_no", "Anticancer: most vs no"), ("onc_mostless_vs_no", "Anticancer: most+less vs no")]
# ---- subtypes ----
sub = []
for line in open(ANA + "subtype_out.txt"):
    m = re.match(r"(.+?)\s+n=(\d+) pos=(\d+) \| structure AUROC ([\d.]+) \[([\d.]+)-([\d.]+)\] \| tumour PCA10 AUROC ([\d.]+) \[([\d.]+)-([\d.]+)\]", line)
    if m: sub.append(dict(name=m.group(1).strip(), n=int(m.group(2)), s=float(m.group(4)), slo=float(m.group(5)), shi=float(m.group(6)),
                          t=float(m.group(7)), tlo=float(m.group(8)), thi=float(m.group(9))))
sub = pd.DataFrame(sub)
# ---- C1 ----
c1 = pd.read_csv(U + "results/confirmatory/C1/oof_predictions.csv")
ci = pd.read_csv("/home/claude/ana/c1/c1_with_indication.csv")
ror = pd.read_csv(U + "results/confirmatory/C1/faers_ror.csv")
rng = np.random.default_rng(2026)
def boot_auc(y, p, n=2000):
    bs = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max(): bs.append(roc_auc_score(y[i], p[i]))
    return roc_auc_score(y, p), np.percentile(bs, 2.5), np.percentile(bs, 97.5)

fig = plt.figure(figsize=(7.1, 6.9))
gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.75, left=0.2, right=0.985, top=0.94, bottom=0.13)

# a - D1, all lines
def d1_panel(ax, df, tasks, letter, title):
    panel(ax, letter, title, x=-1.05, y=1.08); ax.title.set_position((0, 1.04))
    ylabels, ypos, y0 = [], [], 0.6
    for tk, tname in tasks:
        sub_ = df[(df.task == tk) & (df.model == "logit")].set_index("feat")
        n_ = int(sub_.n.iloc[0]) if len(sub_) else 0
        ax.text(0.205, y0 + 0.15, f"{tname} (n = {n_})", fontsize=5.8, fontweight="bold", color=INK, va="bottom")
        for k, (f, fl, col) in enumerate(FEATS):
            yy = y0 + 0.9 + k * 0.62
            if f in sub_.index:
                r = sub_.loc[f]; ax.errorbar(r.auc, yy, xerr=r.sd, fmt="o", ms=3.4, color=col, ecolor=col, elinewidth=0.8, capsize=0)
            ylabels.append(fl); ypos.append(yy)
        y0 += 0.9 + len(FEATS) * 0.62 + 0.35
    ax.axvline(0.5, color=GRAY, lw=0.7, ls="--")
    ax.set_yticks(ypos); ax.set_yticklabels(ylabels, fontsize=5.3); ax.set_ylim(y0 - 0.2, 0.25)
    ax.set_xlim(0.2, 1.0); ax.set_xlabel("AUROC (mean ± SD, repeated 5-fold CV)")
ax = fig.add_subplot(gs[0, 0]); d1_panel(ax, d1, TASKS[:2], "a", "DICTrank, all 50 lines")
ax = fig.add_subplot(gs[1, 0]); d1_panel(ax, d1, TASKS[2:], "d", "DICTrank, anticancer drugs")

# b - lung-only Tahoe features
ax = fig.add_subplot(gs[0, 1]); panel(ax, "b", "Lung lines only (15)", x=-0.22)
tasks_l = TASKS[:2]; w = 0.36
for j, (tk, tname) in enumerate(tasks_l):
    a_all = d1[(d1.task == tk) & (d1.model == "logit")].set_index("feat").auc
    a_lng = d1l[(d1l.task == tk) & (d1l.model == "logit")].set_index("feat").auc
    feats = ["structure", "bio_PCA10_alldoses", "bio_programs", "structure+programs"]
    x = np.arange(len(feats)) + j * (len(feats) + 1)
    ax.bar(x - w / 2, [a_all[f] for f in feats], width=w, color=LIGHTGRAY, label="All lines" if j == 0 else None)
    ax.bar(x + w / 2, [a_lng[f] for f in feats], width=w, color=[CHEM, TAHOE, TAHOE, COMB], label="Lung lines" if j == 0 else None)
    for xi, f in zip(x, feats): ax.text(xi + w / 2, a_lng[f] + 0.015, f"{a_lng[f]:.2f}", ha="center", fontsize=5, color=INK2)
    ax.text(x.mean(), 0.965, tname.replace(" concern",""), ha="center", fontsize=5.8, fontweight="bold")
lab_short = ["Chem.", "PCA-10", "Programs", "Chem.+T."]
ax.set_xticks(list(np.arange(4)) + list(np.arange(4) + 5)); ax.set_xticklabels(lab_short * 2, rotation=45, ha="right", fontsize=5.6)
ax.axhline(0.5, color=GRAY, lw=0.7, ls="--"); ax.set_ylim(0.3, 1.0); ax.set_ylabel("AUROC (logistic)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=LIGHTGRAY, label="All 50 lines"), Patch(color=INK2, label="15 lung lines (coloured)")], fontsize=5.4, loc="upper right", bbox_to_anchor=(1.0, 0.93))

# c - DICTrank subtypes
ax = fig.add_subplot(gs[0, 2]); panel(ax, "c", "DICTrank cardiotoxicity subtypes", x=-0.42)
names = {"heart damage vs no": "Heart damage", "arrhythmia vs no": "Arrhythmia", "mixed vs no": "Mixed",
         "heart damage (oncology) vs no": "Heart damage\n(anticancer)"}
yy = np.arange(len(sub))
ax.errorbar(sub.s, yy - 0.15, xerr=[sub.s - sub.slo, sub.shi - sub.s], fmt="o", ms=3.4, color=CHEM, elinewidth=0.8, label="Chemistry")
ax.errorbar(sub.t, yy + 0.15, xerr=[sub.t - sub.tlo, sub.thi - sub.t], fmt="s", ms=3.2, color=TAHOE, elinewidth=0.8, label="Tahoe PCA-10")
ax.set_yticks(yy); ax.set_yticklabels([f"{names[n]}\n(n = {k})" for n, k in zip(sub.name, sub.n)], fontsize=5.6); ax.invert_yaxis()
ax.axvline(0.5, color=GRAY, lw=0.7, ls="--"); ax.set_xlim(0, 1); ax.set_xlabel("AUROC (95% CV interval)")
ax.legend(fontsize=5.6, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2)

# e - C1 ROC curves
ax = fig.add_subplot(gs[1, 1]); panel(ax, "e", "C1: FAERS cardiac signal", x=-0.22)
y = c1.signal.values
for col, nm, colr in [("oof_chemistry", "Chemistry", CHEM), ("oof_tahoe", "Tahoe", TAHOE), ("oof_chemistry+tahoe", "Chemistry + Tahoe", COMB)]:
    fpr, tpr, _ = roc_curve(y, c1[col]); a, lo, hi = boot_auc(y, c1[col].values)
    ax.plot(fpr, tpr, color=colr, lw=1.2, label=f"{nm}  {a:.2f} ({lo:.2f}–{hi:.2f})")
ax.plot([0, 1], [0, 1], color=GRAY, lw=0.7, ls="--")
ax.set_xlabel("False-positive rate"); ax.set_ylabel("True-positive rate"); ax.set_xlim(0, 1); ax.set_ylim(0, 1.0)
ax.legend(fontsize=5.3, loc="upper left", bbox_to_anchor=(-0.05, -0.2), title="AUROC (95% CI)", title_fontsize=5.5, ncol=1)
# primary difference
d = []
for _ in range(2000):
    i = rng.integers(0, len(y), len(y))
    if y[i].min() != y[i].max(): d.append(roc_auc_score(y[i], c1["oof_chemistry+tahoe"].values[i]) - roc_auc_score(y[i], c1["oof_chemistry"].values[i]))
dd = roc_auc_score(y, c1["oof_chemistry+tahoe"]) - roc_auc_score(y, c1["oof_chemistry"])
ax.text(0.03, 0.97, f"n = 128 drugs\nΔAUROC (Chem.+Tahoe − Chem.)\n= {dd:+.3f} ({np.percentile(d,2.5):+.3f} to {np.percentile(d,97.5):+.3f})",
        transform=ax.transAxes, fontsize=5.3, va="top", color=INK)

# f - confounding by indication
ax = fig.add_subplot(gs[1, 2]); panel(ax, "f", "FAERS: confounding by indication", x=-0.42)
x = ci.merge(ror[["key", "ROR", "ROR_lo95"]], on="key", how="left")
x = x[x.ROR > 0]
grp = [("Cardiovascular\nindication", x[x.cardiovascular_indication]), ("Other\nindications", x[~x.cardiovascular_indication])]
for i, (nm, g) in enumerate(grp):
    jit = rng.uniform(-0.18, 0.18, len(g))
    ax.scatter(np.full(len(g), i) + jit, np.log2(g.ROR), s=9, c=np.where(g.signal == 1, ORANGE, GRAY), edgecolor="white", linewidth=0.3, zorder=3)
    ax.boxplot(np.log2(g.ROR), positions=[i], widths=0.55, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="none", edgecolor=INK2, lw=0.6), medianprops=dict(color=INK, lw=0.9),
               whiskerprops=dict(color=INK2, lw=0.6), capprops=dict(color=INK2, lw=0.6))
    ax.text(i, ax.get_ylim()[1] if False else 3.75, f"signal {g.signal.mean():.0%}\n(n = {len(g)})", ha="center", fontsize=5.6, color=INK2)
ax.axhline(0, color=GRAY, lw=0.7, ls="--")
ax.set_xticks([0, 1]); ax.set_xticklabels([g[0] for g in grp], fontsize=6); ax.set_ylabel("log2 reporting odds ratio (cardiac)")
ax.set_ylim(-4, 4.7); ax.set_xlim(-0.6, 1.6)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=ORANGE, ms=4, label="Signal"), Line2D([], [], marker="o", ls="", color=GRAY, ms=4, label="No signal")],
          fontsize=5.6, loc="lower left")
for lbl in []:
    r = x[x.drug_tahoe.str.startswith(lbl)]
    if len(r): ax.annotate(lbl, (0, np.log2(r.ROR.iloc[0])), xytext=(-0.55, np.log2(r.ROR.iloc[0]) + 0.5), fontsize=5, color=INK2,
                           arrowprops=dict(arrowstyle="-", lw=0.4, color=GRAY), va="center")
save(fig, "Figure2_cardiotoxicity")

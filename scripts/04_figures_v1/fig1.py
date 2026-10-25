"""Figure 1 - Study design and positive controls (D2)."""
import sys; sys.path.insert(0, "/home/claude/figs")
from figcommon import *
from matplotlib.patches import FancyBboxPatch
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import balanced_accuracy_score
style()

# ---------------- data (same computations as the discovery analysis D2) ----------------
lab = pd.read_csv(ANA + "labels_meta.csv").drop_duplicates("key").set_index("key")
X3 = np.load(ANA + "X_all.npy"); keys = pd.read_csv(ANA + "keys_all.csv").iloc[:, 0].tolist()
keep = (X3[:, 2, :] != 0).mean(0) >= 0.9; X3 = X3[:, :, keep]
def cor(a, b): a = a - a.mean(); b = b - b.mean(); return (a @ b) / np.sqrt((a @ a) * (b @ b))
n = len(keys)
same = np.array([cor(X3[i, 1], X3[i, 2]) for i in range(n)])
diff = np.array([cor(X3[i, 1], X3[j, 2]) for i in range(n) for j in range(n) if i != j])
moa = lab.loc[keys, "moa-fine"].fillna("unclear").astype(str)
onc = ["DNA synthesis/repair inhibitor", "EGFR/ERBB inhibitor", "Other TK inhibitor", "Multi-TK inhibitor",
       "Microtubule inhibitor", "MEK inhibitor", "RAF inhibitor", "CDK inhibitor", "HDAC inhibitor",
       "Proteasome inhibitor", "JAK/STAT inhibitor", "Sonic inhibitor"]
m = moa.isin(onc).values; T = X3[m, 2, :]; y = np.asarray(moa[m].values, dtype=object)
vc = pd.Series(y).value_counts(); ok = np.isin(y, vc[vc >= 3].index); T, y = T[ok], y[ok]
g = np.argsort(-T.var(0))[:3000]; Z = T[:, g]; Z = (Z - Z.mean(1, keepdims=True)) / Z.std(1, keepdims=True)
p = cross_val_predict(KNeighborsClassifier(1, metric="correlation"), Z, y, cv=StratifiedKFold(3, shuffle=True, random_state=0))
obs = balanced_accuracy_score(y, p)
rng = np.random.default_rng(0); null = np.array([balanced_accuracy_score(rng.permutation(y), p) for _ in range(2000)])
classes = sorted(set(y), key=lambda c: (-np.sum(y == c), c))

# ---------------- figure ----------------
fig = plt.figure(figsize=(7.1, 8.4))
gsA = fig.add_gridspec(1, 1, top=0.985, bottom=0.785, left=0.17, right=0.98)
gs = fig.add_gridspec(2, 3, top=0.70, bottom=0.04, hspace=0.62, wspace=0.62, height_ratios=[1, 1.1], left=0.17, right=0.98)

# a - study design schematic (full width)
ax = fig.add_subplot(gsA[0, 0]); ax.set_axis_off(); panel(ax, "a", "Study design", x=-0.02, y=1.0)
ax.set_xlim(-0.8, 20.2); ax.set_ylim(0.1, 5.25)
def box(x, y, w, h, txt, fc, ec, fs=6.3, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.15", fc=fc, ec=ec, lw=0.7))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", color=INK)
ax.text(3.5, 5.0, "Open data (no IRB)", ha="center", fontsize=6.5, fontweight="bold", color=INK2)
ax.text(14.45, 5.0, "Analyses (D = discovery, C = confirmatory)", ha="center", fontsize=6.5, fontweight="bold", color=INK2)
box(-0.6, 3.35, 8.2, 1.35, "Tahoe-100M single-cell drug atlas\n379 drugs · 50 cancer cell lines (15 lung)\n3 doses · 24 h · pseudobulk differential expression", "#e8f1fc", BLUE)
box(-0.6, 1.8, 8.2, 1.2, "Clinical cardiotoxicity labels\nFDA DICTrank · FDA FAERS (openFDA)", "#fdece4", ORANGE)
box(-0.6, 0.2, 8.2, 1.25, "Cell-line drug sensitivity and baseline state\nPRISM · GDSC · CCLE RNA-seq,\nmass-spectrometry proteomics, RPPA", "#e3f6ef", AQUA)
items = [("D2   Reproducibility and mechanism of action (positive controls)", 4.45),
         ("D1 · C1   Clinical cardiotoxicity: Tahoe response vs chemistry", 3.65),
         ("D3   Pooled 24-h cell counts as a viability readout", 2.85),
         ("D4 · C2   Response magnitude and proliferation program vs sensitivity", 2.05),
         ("C3   Tahoe response vs baseline RNA vs baseline protein (MS, RPPA)", 1.25),
         ("Pre-registered plan · permutation nulls · robustness · dated amendments", 0.45)]
for t, yy in items:
    pre = "Pre-reg" in t
    box(8.9, yy - 0.31, 11.1, 0.62, t, "#f0efec" if pre else "white", GRAY if pre else INK2, fs=6.3)
for yy in (4.02, 2.4, 0.82):
    ax.annotate("", xy=(8.8, yy), xytext=(7.7, yy), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK2))

# b - analysis sizes
ax = fig.add_subplot(gs[0, 0]); panel(ax, "b", "Drugs per analysis", x=-0.62)
rows = [("D1 DICTrank most vs no", 100), ("D1 DICTrank most+less vs no", 167), ("C1 FAERS signal", 128),
        ("D4 PRISM", 51), ("D4 GDSC", 43), ("C2 PRISM", 44), ("C2 GDSC", 37), ("C3 PRISM", 43), ("C3 GDSC", 37),
        ("D2 mechanism classes", len(y))]
names, vals = zip(*rows[::-1])
ax.barh(range(len(vals)), vals, color=BLUE, height=0.62)
for i, v in enumerate(vals): ax.text(v + 4, i, str(v), va="center", fontsize=5.8, color=INK2)
ax.set_yticks(range(len(vals))); ax.set_yticklabels(names, fontsize=5.8); ax.set_xlabel("Number of drugs"); ax.set_xlim(0, 210)

# c - dose reproducibility
ax = fig.add_subplot(gs[0, 1]); panel(ax, "c", "Signature reproducibility")
bins = np.linspace(-0.4, 1, 29)
ax.hist(diff, bins=bins, density=True, color=GRAY, alpha=0.75, label=f"Different drugs (median {np.median(diff):.2f})")
ax.hist(same, bins=bins, density=True, color=BLUE, alpha=0.8, label=f"Same drug (median {np.median(same):.2f})")
hmax = max(np.histogram(same, bins, density=True)[0].max(), np.histogram(diff, bins, density=True)[0].max())
ax.vlines([np.median(same)], 0, hmax * 1.02, color=BLUE, lw=0.8, ls="--"); ax.vlines([np.median(diff)], 0, hmax * 1.02, color=INK2, lw=0.8, ls="--")
ax.set_xlabel("Correlation of mid- vs high-dose\nsignatures"); ax.set_ylabel("Density")
ax.legend(loc="upper left", fontsize=5.6, bbox_to_anchor=(0, 1.02)); ax.set_ylim(0, ax.get_ylim()[1] * 1.35)

# d - MoA recognition vs permutation null
ax = fig.add_subplot(gs[0, 2]); panel(ax, "d", "Mechanism recognition")
pv = ((null >= obs).sum() + 1) / (len(null) + 1)
ax.hist(null, bins=30, color=GRAY, alpha=0.8, label="Permutation null")
ax.vlines([obs], 0, np.histogram(null, 30)[0].max() * 1.02, color=BLUE, lw=1.4, label=f"Observed {obs:.2f} (p = {pv:.4f})")
ax.vlines([1 / len(classes)], 0, np.histogram(null, 30)[0].max() * 1.02, color=INK2, lw=0.8, ls=":", label=f"Chance {1/len(classes):.2f}")
ax.set_xlabel("Balanced accuracy\n(1-nearest neighbour, 11 classes)"); ax.set_ylabel("Permutations")
ax.legend(fontsize=5.6, loc="upper left", bbox_to_anchor=(0, 1.02)); ax.set_ylim(0, ax.get_ylim()[1] * 1.45)

# e - per-class recognition
ax = fig.add_subplot(gs[1, 0]); panel(ax, "e", "Correct by mechanism class", x=-0.62)
tot = [int((y == c).sum()) for c in classes]; cor_ = [int(((y == c) & (p == c)).sum()) for c in classes]
short = [c.replace("Sonic inhibitor", "Hedgehog inhibitor").replace(" inhibitor", "").replace("DNA synthesis/repair", "DNA synth./repair") for c in classes]
yy = np.arange(len(classes))[::-1]
ax.barh(yy, tot, color=LIGHTGRAY, height=0.62, label="Drugs in class")
ax.barh(yy, cor_, color=BLUE, height=0.62, label="Correctly recognised")
for i, (t, c) in enumerate(zip(tot, cor_)): ax.text(t + 0.4, yy[i], f"{c}/{t}", va="center", fontsize=5.6, color=INK2)
ax.set_yticks(yy); ax.set_yticklabels(short, fontsize=5.8); ax.set_xlabel("Drugs")
ax.legend(fontsize=5.6, loc="lower right", bbox_to_anchor=(1.02, -0.02)); ax.set_xlim(0, max(tot) + 13)

# f - class-averaged signature similarity
ax = fig.add_subplot(gs[1, 1:]); panel(ax, "f", "Mean signature similarity between mechanism classes", x=-0.1)
Cfull = np.corrcoef(Z); K = len(classes); M = np.zeros((K, K))
for i, ci in enumerate(classes):
    for j, cj in enumerate(classes):
        a_ = np.where(y == ci)[0]; b2 = np.where(y == cj)[0]
        blk = Cfull[np.ix_(a_, b2)]
        M[i, j] = blk[np.triu_indices(len(a_), 1)].mean() if i == j else blk.mean()
im = ax.imshow(M, cmap="RdBu_r", vmin=-0.5, vmax=0.5, interpolation="nearest")
for i in range(K):
    for j in range(K):
        ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=4.6,
                color="white" if abs(M[i, j]) > 0.33 else INK, fontweight="bold" if i == j else "normal")
ax.set_xticks(range(K)); ax.set_xticklabels(short, rotation=45, ha="right", fontsize=5.6)
ax.set_yticks(range(K)); ax.set_yticklabels(short, fontsize=5.6)
for spine in ax.spines.values(): spine.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("Mean Pearson r", fontsize=6); cb.ax.tick_params(labelsize=5.6)
diag = np.diag(M).mean(); off = M[~np.eye(K, dtype=bool)].mean()
ax.text(0.5, -0.36, f"Within-class (diagonal) mean r = {diag:.2f}   ·   between-class mean r = {off:.2f}", transform=ax.transAxes,
        ha="center", va="top", fontsize=6, color=INK2)
pd.DataFrame(M, index=classes, columns=classes).to_csv(OUT + "source_Fig1f.csv")

save(fig, "Figure1_design_positive_controls")
pd.DataFrame({"class": classes, "n": tot, "correct": cor_}).to_csv(OUT + "source_Fig1e.csv", index=False)
print("D2:", np.median(same), np.median(diff), obs, (null >= obs).mean())

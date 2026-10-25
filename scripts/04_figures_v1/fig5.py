"""Figure 5 - Lung adenocarcinoma/NSCLC focus: data coverage and exploratory lung-subset evaluation."""
import sys, gzip; sys.path.insert(0, "/home/claude/figs")
from figcommon import *
style()
EX = "/home/claude/mirror/results/exploratory/lung_subset/"
lrs = pd.read_csv(U + "results/response_vs_sensitivity/line_response_strength.csv")
cov = pd.read_csv(U + "results/feasibility/proteomics_coverage.csv")
pairs = pd.read_csv(U + "results/feasibility/pairs_with_labels.csv")
pred = pd.read_csv(U + "results/confirmatory/C3_predictions_k5_a1.0_top.csv")
ex = pd.read_csv(EX + "per_drug_subset_calibrated.csv"); exs = pd.read_csv(EX + "summary.csv")
meas = set(lrs.Cell_Name_Vevo)
lung = cov[cov.lung & cov.cell_name.isin(meas)].copy()
with gzip.open(U + "data/raw/depmap/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz", "rt") as f: hdr = f.readline().strip().split("\t")
rna = {h.split("_")[0] for h in hdr[2:]}
rpk = {str(x).split("_")[0] for x in pd.read_csv(U + "data/raw/depmap/CCLE_RPPA_20181003.csv", usecols=[0]).iloc[:, 0]}
top = lrs[lrs.dose_rank == 3]
lung["Tahoe drugs"] = lung.cell_name.map(top.groupby("Cell_Name_Vevo").drug_key.nunique())
lung["RNA-seq"] = lung.cell_key.isin(rna).astype(int)
lung["MS protein"] = lung.has_ccle_proteomics.astype(int)
lung["RPPA"] = lung.cell_key.isin(rpk).astype(int)
lung["PRISM drugs"] = lung.cell_name.map(pairs.dropna(subset=["prism_auc"]).groupby("Cell_Name_Vevo").drug_key.nunique()).fillna(0).astype(int)
lung["GDSC drugs"] = lung.cell_name.map(pairs.dropna(subset=["gdsc_auc"]).groupby("Cell_Name_Vevo").drug_key.nunique()).fillna(0).astype(int)
lung["In C3"] = lung.cell_name.isin(set(pred[pred.lung].cell_name)).astype(int)
lung = lung.sort_values(["In C3", "cell_name"], ascending=[False, True])
COLS = ["Tahoe drugs", "RNA-seq", "MS protein", "RPPA", "PRISM drugs", "GDSC drugs", "In C3"]
lung[["cell_name"] + COLS].to_csv(OUT + "source_Fig5a.csv", index=False)

fig = plt.figure(figsize=(7.1, 6.4))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.6, left=0.1, right=0.985, top=0.94, bottom=0.09, height_ratios=[1, 1.05])

# a - coverage matrix
ax = fig.add_subplot(gs[0, 0:2]); panel(ax, "a", "Data coverage of the 15 lung lines profiled in Tahoe-100M", x=-0.13)
nL, nC = len(lung), len(COLS)
for j, c in enumerate(COLS):
    for i, (_, r) in enumerate(lung.iterrows()):
        v = r[c]; binary = c in ["RNA-seq", "MS protein", "RPPA", "In C3"]
        if binary:
            fc = (AQUA if c != "In C3" else BLUE) if v else "#f3f2ef"; txt = "+" if v else "–"
            tc = "white" if v else GRAY
        else:
            mx = {"Tahoe drugs": 379, "PRISM drugs": 45, "GDSC drugs": 39}[c]
            fc = plt.cm.Greys(0.12 + 0.55 * v / mx); txt = str(int(v)); tc = "white" if v / mx > 0.75 else INK
        ax.add_patch(plt.Rectangle((i - 0.46, j - 0.44), 0.92, 0.88, fc=fc, ec="white", lw=0.8))
        ax.text(i, j, txt, ha="center", va="center", fontsize=5.4, color=tc, fontweight="bold" if binary and v else "normal")
ax.set_xlim(-0.6, nL - 0.4); ax.set_ylim(nC - 0.5, -0.5)
ax.set_xticks(range(nL)); ax.set_xticklabels(lung.cell_name, rotation=45, ha="right", fontsize=5.8)
ax.set_yticks(range(nC)); ax.set_yticklabels(COLS, fontsize=6)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.tick_params(length=0)

# b - response strength by line
ax = fig.add_subplot(gs[0, 2]); panel(ax, "b", "Tahoe response strength by line", x=-0.35)
st = top.groupby("Cell_Name_Vevo").n_de.median().sort_values()
isl = st.index.isin(lung.cell_name)
ax.barh(range(len(st)), st.values, color=np.where(isl, BLUE, LIGHTGRAY), height=0.8)
ax.set_yticks([]); ax.set_ylabel("50 cell lines (ranked)"); ax.set_xlabel("Median no. DE genes per drug\n(top dose)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=BLUE, label="Lung (15)"), Patch(color=LIGHTGRAY, label="Other (35)")], fontsize=5.6, loc="lower right")
from scipy.stats import mannwhitneyu
p_ = mannwhitneyu(st[isl], st[~isl]).pvalue
ax.text(0.97, 0.3, f"Mann–Whitney p = {p_:.2f}", transform=ax.transAxes, ha="right", fontsize=5.4, color=INK2)

# c, d - exploratory: lung vs non-lung ranking (null-calibrated)
SETS4 = ["RNA", "Protein_MS", "RPPA", "Response"]
SH = {"RNA": "RNA", "Protein_MS": "Protein (MS)", "RPPA": "Protein (RPPA)", "Response": "Tahoe response"}
rng = np.random.default_rng(5)
for letter, src, pos in [("c", "PRISM", gs[1, 0]), ("d", "GDSC", gs[1, 1])]:
    ax = fig.add_subplot(pos); nd = exs[exs.source == src].n_drugs.iloc[0]
    panel(ax, letter, f"Lung vs other lines, {src} ({nd} drugs)", x=-0.3)
    for i, s in enumerate(SETS4):
        for k, (sub, alpha_, off) in enumerate([("lung", 1.0, -0.19), ("non-lung", 0.45, 0.19)]):
            r = ex[(ex.source == src) & (ex.set == s) & (ex.subset == sub)].calibrated_rho
            col = SET_COLORS[s] if sub == "lung" else GRAY
            ax.scatter(i + off + rng.uniform(-0.09, 0.09, len(r)), r, s=3, color=col, alpha=0.6, lw=0)
            ax.hlines(r.median(), i + off - 0.15, i + off + 0.15, color=INK, lw=1.1, zorder=4)
            q = exs[(exs.source == src) & (exs.set == s) & (exs.subset == sub)].fdr_bh.iloc[0]
            ax.text(i + off, 1.05, f"{r.median():+.2f}" + ("*" if q < 0.05 else ""), ha="center", fontsize=4.8, rotation=90, va="bottom",
                    color=INK, fontweight="bold" if q < 0.05 else "normal")
    ax.axhline(0, color=GRAY, lw=0.7, ls="--")
    ax.set_xticks(range(4)); ax.set_xticklabels([SH[s] for s in SETS4], rotation=40, ha="right", fontsize=5.8)
    ax.set_ylim(-1.3, 1.45); ax.set_yticks([-1, -0.5, 0, 0.5, 1]); ax.set_ylabel("Calibrated ρ within subset")
    if src == "PRISM":
        ax.legend(handles=[Patch(color=INK2, label="Lung lines (coloured; ≈10)"), Patch(color=GRAY, alpha=0.45, label="Other lines (≈18)")],
                  fontsize=5.4, loc="lower left", ncol=1)

# e - per drug: lung-line ranking by RPPA vs Tahoe response (both screens, all drugs)
ax = fig.add_subplot(gs[1, 2]); panel(ax, "e", "Lung lines: RPPA vs Tahoe response", x=-0.3)
w = ex[ex.subset == "lung"].pivot_table(index=["source", "drug_key"], columns="set", values="calibrated_rho").reset_index()
for src, col, mk in [("PRISM", BLUE, "o"), ("GDSC", ORANGE, "s")]:
    g = w[w.source == src]
    ax.scatter(g.Response, g.RPPA, s=10, marker=mk, color=col, edgecolor="white", lw=0.3, alpha=0.9, zorder=3,
               label=f"{src}: RPPA > response in {(g.RPPA > g.Response).sum()}/{len(g)}")
lim = (-1.05, 1.15); ax.plot(lim, lim, color=GRAY, lw=0.7, ls="--"); ax.set_xlim(lim); ax.set_ylim(lim)
ax.axhline(0, color=LIGHTGRAY, lw=0.5); ax.axvline(0, color=LIGHTGRAY, lw=0.5); ax.set_aspect("equal")
ax.set_xlabel("Calibrated ρ, Tahoe response"); ax.set_ylabel("Calibrated ρ, baseline RPPA")
ax.legend(fontsize=5.2, loc="upper center", bbox_to_anchor=(0.45, -0.3))
w.to_csv(OUT + "source_Fig5e.csv", index=False)
save(fig, "Figure5_lung_focus")
print(st[isl].median(), st[~isl].median(), p_)

"""
v2 / T2-H1 (Amendment 3, frozen 18 Oct 2026) - cardiac context: Tahoe tumour-line signatures vs iPSC-cardiomyocyte
signatures (Hansen et al. 2024, Nat Commun 15:7968; Supplementary Data 05 = averaged DEGs over cell lines, 48 h).
(a) cross-context identity: Pearson r between the Tahoe signature (mean Wald statistic over 50 lines, top dose) and the
    cardiomyocyte signature (signed -log10 P of the top-600 DEGs, 0 otherwise) over shared genes, same drug vs different
    drugs; one-sided permutation test (10,000 shuffles of the drug matching), p < 0.05.  Sensitivity: Spearman r, and the
    drug-selective signatures (Supplementary Data 10).
(b) toxicity comparison only if >= 30 shared drugs with DICTrank/FAERS labels (otherwise descriptive).
"""
import os, sys, importlib, numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/mirror/scripts"); norm = importlib.import_module("02_drug_overlap").norm
U = "/mnt/user-data/uploads/Protein and RNA/"; H = U + "data/raw/v2/hansen2024/"
OUT = "/home/claude/v2/results/T2H1"; os.makedirs(OUT, exist_ok=True)

def load_tahoe():
    a = np.load("/home/claude/v2/results/cache/tahoe_alldrug_signatures.npz", allow_pickle=True)
    b = np.load(U + "results/confirmatory/C1/tahoe_features_C1.npz", allow_pickle=True)
    sig = {}
    for z in [b, a]:                                  # identical values where both exist
        for i, k in enumerate(z["keys"]): sig[k] = pd.Series(z["X"][i, 2, :], index=list(z["genes"]))
    return sig

def hansen(fname):
    d = pd.read_csv(H + fname, sep="\t")
    d["key"] = d["Drug name"].map(norm)
    M = d.pivot_table(index="key", columns="NCBI gene symbol", values="Signed -log10(pvalue)", aggfunc="mean").fillna(0)
    tox = d.drop_duplicates("key").set_index("key")["Is cardiotoxic"]
    return M, tox

def corr(a, b, method):
    if method == "pearson": return np.corrcoef(a, b)[0, 1]
    return pd.Series(a).rank().corr(pd.Series(b).rank())

def run(M, sig, label, rng):
    shared = sorted(set(M.index) & set(sig))
    genes = sorted(set(M.columns) & set.intersection(*[set(sig[k].index) for k in shared]))
    T = np.array([sig[k].reindex(genes).values for k in shared]); C = M.loc[shared, genes].values
    out = []
    for method in ["pearson", "spearman"]:
        R = np.array([[corr(T[i], C[j], method) for j in range(len(shared))] for i in range(len(shared))])
        same = np.diag(R); diff = R[~np.eye(len(shared), dtype=bool)]
        obs = same.mean(); null = np.array([R[np.arange(len(shared)), rng.permutation(len(shared))].mean() for _ in range(10000)])
        p = (np.sum(null >= obs) + 1) / 10001
        rank_pct = np.array([(R[i] < R[i, i]).sum() / (len(shared) - 1) for i in range(len(shared))])  # own CM signature ranked vs others
        out.append(dict(signature=label, method=method, n_drugs=len(shared), n_genes=len(genes), mean_same=obs,
                        median_same=np.median(same), mean_different=diff.mean(), perm_p=p,
                        mean_rank_percentile_own=rank_pct.mean(), n_own_best=int((rank_pct == 1).sum())))
        if method == "pearson":
            pd.DataFrame(R, index=shared, columns=shared).to_csv(os.path.join(OUT, f"corr_matrix_{label}.csv"))
            pd.DataFrame({"drug": shared, "r_same": same, "own_rank_percentile": rank_pct}).to_csv(os.path.join(OUT, f"per_drug_{label}.csv"), index=False)
    return out

def main():
    rng = np.random.default_rng(2026); sig = load_tahoe()
    M5, tox = hansen("Supplementary Data 05 - Averaged DEGs across all cell lines.txt")
    M10, _ = hansen("Supplementary Data 10 - Averaged drug-selective DEGs across all cell lines.txt")
    res = run(M5, sig, "averaged_DEGs", rng) + run(M10, sig, "drug_selective_DEGs", rng)
    R = pd.DataFrame(res); R.to_csv(os.path.join(OUT, "summary.csv"), index=False)
    prim = R[(R.signature == "averaged_DEGs") & (R.method == "pearson")].iloc[0]
    lines = [f"T2-H1a cross-context identity (primary: averaged DEGs, Pearson): n = {prim.n_drugs} drugs, {prim.n_genes} genes",
             f"  mean r same drug = {prim.mean_same:+.3f} vs different drugs = {prim.mean_different:+.3f}; one-sided permutation p = {prim.perm_p:.4f}",
             "  DECISION: " + ("identity preserved across contexts (p < 0.05)" if prim.perm_p < 0.05 else "identity NOT detectably preserved (p >= 0.05)")]
    # (b) feasibility
    lab = pd.read_csv("/home/claude/ana/labels_meta.csv").drop_duplicates("key").set_index("key")
    c1 = pd.read_csv("/home/claude/ana/c1/c1_with_indication.csv").set_index("key")
    shared = sorted(set(M5.index) & set(sig))
    n_lab = sum((k in lab.index and lab.loc[k, "dict_concern"] in ["most", "less", "no"]) or (k in c1.index) for k in shared)
    lines.append(f"T2-H1b: shared drugs with DICTrank (most/less/no) or FAERS labels = {n_lab} -> "
                 + ("model comparison run" if n_lab >= 30 else "below the pre-specified 30; reported descriptively only"))
    desc = pd.DataFrame({"drug": shared, "hansen_is_cardiotoxic": [tox.get(k) for k in shared],
                         "dictrank": [lab.loc[k, "dict_concern"] if k in lab.index else None for k in shared],
                         "faers_signal": [c1.loc[k, "signal"] if k in c1.index else None for k in shared]})
    desc.to_csv(os.path.join(OUT, "H1b_descriptive_labels.csv"), index=False)
    txt = "\n".join(lines) + "\n\n" + R.round(4).to_string(index=False)
    open(os.path.join(OUT, "summary.txt"), "w").write(txt); print(txt)

if __name__ == "__main__":
    main()

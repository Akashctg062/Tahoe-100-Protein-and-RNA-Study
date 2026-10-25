"""
v2 / T1-C (Amendment 3, frozen 18 Oct 2026) - RPPA antibody / drug-target concordance.
Units: drug x screen pairs of the C3 common-line sets whose Tahoe target genes fall in a pre-specified pathway group.
Statistic per unit: Spearman rho of every RPPA antibody with AUC across lines; concordance = mean percentile rank (by |rho|,
1 = strongest) of the antibodies whose gene is in the drug's target pathway group(s).  Summary = median over units.
Null: 10,000 permutations re-assigning matched-antibody sets among units within screen.  One-sided p (higher = concordant).
"""
import os, re, sys, pickle, importlib
import numpy as np, pandas as pd
OUT = "/home/claude/v2/results/T1C"; os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, "/home/claude/mirror/scripts"); norm = importlib.import_module("02_drug_overlap").norm

GROUPS = {
 "ERBB": ["EGFR", "ERBB2", "ERBB3"],
 "MAPK": ["ARAF", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RPS6KA1", "RPS6KA2", "RPS6KA3"],
 "PI3K_AKT_MTOR": ["PIK3CA", "PIK3R1", "AKT1", "AKT2", "AKT3", "MTOR", "RPS6KB1", "RPS6", "EIF4EBP1", "PDPK1", "PTEN",
                   "TSC1", "TSC2", "AKT1S1", "RICTOR", "RPTOR", "GSK3A", "GSK3B"],
 "CELL_CYCLE": ["CDK1", "CDK2", "CDK4", "CDK6", "RB1", "CCND1", "CCNE1", "CCNE2", "CCNB1", "CDKN1A", "CDKN1B"],
 "RTK": ["MET", "KIT", "KDR", "PDGFRA", "PDGFRB", "FLT3"],
 "SRC": ["SRC", "LCK", "YES1", "LYN"],
 "JAK_STAT": ["JAK1", "JAK2", "STAT3", "STAT5A"],
 "APOPTOSIS": ["BCL2", "BCL2L1", "BAX", "BAK1", "BAD", "BCL2L11"],
 "DNA_DAMAGE": ["ATM", "CHEK1", "CHEK2", "PARP1", "XRCC1", "ERCC1", "MRE11", "RAD50", "RAD51", "BRCA2", "TP53BP1"],
 "NUCLEAR_RECEPTOR": ["ESR1", "PGR", "AR"],
}
# antibody (name prefix before the first "_" or "(") -> gene symbol(s); only antibodies of genes in GROUPS need mapping
AB2GENE = {
 "EGFR": ["EGFR"], "HER2": ["ERBB2"], "HER3": ["ERBB3"],
 "A-Raf": ["ARAF"], "B-Raf": ["BRAF"], "C-Raf": ["RAF1"], "MEK1": ["MAP2K1"], "MAPK": ["MAPK1", "MAPK3"], "ERK2": ["MAPK1"],
 "p90RSK": ["RPS6KA1"], "RSK1-2-3": ["RPS6KA1", "RPS6KA2", "RPS6KA3"],
 "PI3K-p110-alpha": ["PIK3CA"], "PI3K-p85": ["PIK3R1"], "Akt": ["AKT1", "AKT2", "AKT3"], "mTOR": ["MTOR"],
 "p70S6K": ["RPS6KB1"], "S6": ["RPS6"], "4E-BP1": ["EIF4EBP1"], "PDK1": ["PDPK1"], "PTEN": ["PTEN"], "TSC1": ["TSC1"],
 "Tuberin": ["TSC2"], "PRAS40": ["AKT1S1"], "Rictor": ["RICTOR"], "Raptor": ["RPTOR"],
 "GSK3-alpha-beta": ["GSK3A", "GSK3B"], "GSK3": ["GSK3B"], "GSK-3-BETA": ["GSK3B"],
 "CDK1": ["CDK1"], "Rb": ["RB1"], "Cyclin": None, "p21": ["CDKN1A"], "p27": ["CDKN1B"],
 "c-Met": ["MET"], "c-Kit": ["KIT"], "VEGFR2": ["KDR"],
 "Src": ["SRC"], "Lck": ["LCK"], "JAK2": ["JAK2"], "STAT3": ["STAT3"], "STAT5-alpha": ["STAT5A"],
 "Bcl-2": ["BCL2"], "Bcl-xL": ["BCL2L1"], "Bax": ["BAX"], "Bak": ["BAK1"], "Bad": ["BAD"], "Bim": ["BCL2L11"],
 "ATM": ["ATM"], "Chk1": ["CHEK1"], "Chk2": ["CHEK2"], "PARP": ["PARP1"], "XRCC1": ["XRCC1"], "ERCC1": ["ERCC1"],
 "Mre11": ["MRE11"], "Rad50": ["RAD50"], "RAD51": ["RAD51"], "BRCA2": ["BRCA2"], "53BP1": ["TP53BP1"],
 "ER-alpha": ["ESR1"], "PR": ["PGR"], "AR": ["AR"],
}
CYCLIN = {"Cyclin_B1": ["CCNB1"], "Cyclin_D1": ["CCND1"], "Cyclin_E1": ["CCNE1"], "Cyclin_E2_Caution": ["CCNE2"]}

def ab_genes(name):
    if name in CYCLIN: return CYCLIN[name]
    key = re.split(r"[_(]", name)[0].strip()
    return AB2GENE.get(key) or []

def main():
    data = pickle.load(open("/home/claude/v2/results/cache/c3_cache.pkl", "rb"))
    dm = pd.read_csv("/mnt/user-data/uploads/Protein and RNA/data/raw/tahoe_metadata/drug_metadata.csv").drop_duplicates("drug")
    dm["key"] = dm.drug.map(norm)
    tg = dm.dropna(subset=["targets"]).groupby("key").targets.apply(lambda s: sorted({t.strip() for x in s for t in str(x).split(",") if t.strip()}))
    gene2grp = {g: grp for grp, gs in GROUPS.items() for g in gs}
    units = []
    for src in ["PRISM", "GDSC"]:
        for d in data[src]:
            targets = tg.get(d["drug"], [])
            grps = sorted({gene2grp[t] for t in targets if t in gene2grp})
            if not grps: continue
            names = d["rppa_names"]; X = d["X"]["RPPA"]; y = d["y"]
            rho = pd.DataFrame(X, columns=names).rank().corrwith(pd.Series(y).rank(), method="pearson").values
            matched_genes = {g for grp in grps for g in GROUPS[grp]}
            matched = np.array([bool(set(ab_genes(n)) & matched_genes) for n in names])
            if matched.sum() == 0: continue
            units.append(dict(source=src, drug=d["drug"], targets=";".join(targets), groups=";".join(grps),
                              rho=rho, matched=matched, names=names))
    print(f"units: {len(units)} (PRISM {sum(u['source']=='PRISM' for u in units)}, GDSC {sum(u['source']=='GDSC' for u in units)})")
    def score(rho, matched):
        r = np.abs(np.nan_to_num(rho)); pct = (pd.Series(r).rank(ascending=True).values - 0.5) / len(r)  # 1 = strongest
        return pct[matched].mean()
    obs_u = np.array([score(u["rho"], u["matched"]) for u in units]); obs = np.median(obs_u)
    rng = np.random.default_rng(2026); src = np.array([u["source"] for u in units]); null = np.empty(10000)
    for b in range(10000):
        sets = [None] * len(units)
        for s in ["PRISM", "GDSC"]:
            idx = np.where(src == s)[0]; p = rng.permutation(idx)
            for i, j in zip(idx, p): sets[i] = units[j]["matched"]
        null[b] = np.median([score(units[i]["rho"], sets[i]) for i in range(len(units))])
    p = (np.sum(null >= obs) + 1) / 10001
    res = [f"T1-C RPPA target concordance: {len(units)} drug x screen units",
           f"median concordance (mean percentile of matched antibodies, 0.5 = chance) = {obs:.3f}",
           f"permutation null median {np.median(null):.3f} (95% range {np.percentile(null,2.5):.3f}-{np.percentile(null,97.5):.3f}); one-sided p = {p:.4f}",
           "DECISION: " + ("concordance SUPPORTED (p < 0.05)" if p < 0.05 else "concordance NOT supported (p >= 0.05)")]
    for s in ["PRISM", "GDSC"]:
        idx = np.where(src == s)[0]; o = np.median(obs_u[idx]); nl = np.empty(5000)
        for b in range(5000):
            p_ = rng.permutation(idx); nl[b] = np.median([score(units[i]["rho"], units[j]["matched"]) for i, j in zip(idx, p_)])
        res.append(f"  {s}: n = {len(idx)}, median {o:.3f}, null median {np.median(nl):.3f}, p = {(np.sum(nl >= o)+1)/5001:.4f}")
    # secondary: sign of rho for matched phospho-antibodies (expected negative = active pathway in sensitive lines)
    ph = [u["rho"][i] for u in units for i in np.where(u["matched"])[0] if "_p" in u["names"][i]]
    res.append(f"secondary: matched phospho-antibodies, median rho = {np.median(ph):+.3f} (n = {len(ph)}), fraction negative = {np.mean(np.array(ph) < 0):.2f}")
    rows = []
    for u, sc in zip(units, obs_u):
        best = np.argsort(-np.abs(np.nan_to_num(u["rho"])))[:3]
        top_m = [u["names"][i] for i in np.where(u["matched"])[0]]
        rows.append(dict(source=u["source"], drug=u["drug"], targets=u["targets"], groups=u["groups"], n_matched=int(u["matched"].sum()),
                         concordance=sc, top3_antibodies=";".join(f"{u['names'][i]}({u['rho'][i]:+.2f})" for i in best),
                         matched_antibodies=";".join(top_m)))
    pd.DataFrame(rows).sort_values("concordance", ascending=False).to_csv(os.path.join(OUT, "per_unit.csv"), index=False)
    pd.DataFrame({"null_median": null}).to_csv(os.path.join(OUT, "null.csv"), index=False)
    open(os.path.join(OUT, "summary.txt"), "w").write("\n".join(res)); print("\n".join(res))

if __name__ == "__main__":
    main()

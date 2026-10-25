"""
v2 / T1-B (Amendment 3) - leave-one-out (LOO) bias of within-drug prediction: simulation (descriptive).
Model as in C3: centre features on the training rows, PCA to k components on the training rows, ridge (alpha) on
standardised scores with a training-mean intercept; predict the held-out row.  PCA is computed exactly from the
n x n Gram matrix (features are i.i.d. N(0,1), so per-feature scaling is omitted).
Grid (null, y independent of X): n in {10,15,20,27,40,60}, k in {1,3,5,10}, alpha in {0.1,1,10,100}, 500 reps, p=5000.
Also: intercept-only predictor; permutation-calibrated rho (100 perms) and 5-fold CV for k=5, alpha=1;
signal scenario y = b*z + e with true corr(y, z) in {0.2, 0.4}, z = a latent factor carried by 50 features.
"""
import os, numpy as np, pandas as pd, time
rng = np.random.default_rng(2026); OUT = "/home/claude/v2/results/T1B"; os.makedirs(OUT, exist_ok=True)

def rank(x): return np.argsort(np.argsort(x, axis=-1), axis=-1).astype(float)
def spear_rows(P, y):          # P: (m, n) predictions, y: (m, n)
    a, b = rank(P), rank(y); a -= a.mean(1, keepdims=True); b -= b.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))

def fold_scores(K, tr, te, k):
    """PCA scores from the full Gram matrix K = X X^T for training rows tr and test rows te (centred on training mean)."""
    Ktt, Ket = K[np.ix_(tr, tr)], K[np.ix_(te, tr)]
    n = len(tr); one = np.ones((n, n)) / n
    Kc = Ktt - one @ Ktt - Ktt @ one + one @ Ktt @ one
    Kec = Ket - Ket.mean(1, keepdims=True) - Ktt.mean(0, keepdims=True) + Ktt.mean()
    w, U = np.linalg.eigh(Kc); o = np.argsort(w)[::-1][:k]; w, U = np.clip(w[o], 1e-12, None), U[:, o]
    return U * np.sqrt(w), Kec @ U / np.sqrt(w)            # training scores, test scores

def predict(folds, Y, alpha):
    """folds: list of (tr, te, Ztr, Zte); Y: (m, n) several label vectors at once -> predictions (m, n)."""
    P = np.empty_like(Y)
    for tr, te, Ztr, Zte in folds:
        m_, s_ = Ztr.mean(0), Ztr.std(0); s_[s_ == 0] = 1
        A, B = (Ztr - m_) / s_, (Zte - m_) / s_
        H = np.linalg.solve(A.T @ A + alpha * np.eye(A.shape[1]), A.T)      # k x ntr
        ytr = Y[:, tr]; ym = ytr.mean(1, keepdims=True)
        P[:, te] = ym + (ytr - ym) @ H.T @ B.T
    return P

def make_folds(K, n, k, scheme):
    idx = np.arange(n)
    if scheme == "loo": splits = [(idx[idx != i], np.array([i])) for i in range(n)]
    else:
        perm = rng.permutation(n); splits = [(np.setdiff1d(idx, f), f) for f in np.array_split(perm, 5)]
    return [(tr, te, *fold_scores(K, tr, te, k)) for tr, te in splits]

t0 = time.time(); rows = []
GRID_DONE = os.path.exists(os.path.join(OUT, "null_grid_summary.csv"))
for n in ([] if GRID_DONE else [10, 15, 20, 27, 40, 60]):
    for rep in range(500):
        X = rng.standard_normal((n, 5000)); K = X @ X.T; y = rng.standard_normal((1, n))
        loo = {k: make_folds(K, n, k, "loo") for k in [1, 3, 5, 10] if k < n - 1}
        for k, F in loo.items():
            for alpha in [0.1, 1, 10, 100]:
                rows.append(dict(n=n, k=k, alpha=alpha, rep=rep, rho=spear_rows(predict(F, y, alpha), y)[0]))
        # intercept-only LOO
        yi = y[0]; pi = (yi.sum() - yi) / (n - 1)
        rows.append(dict(n=n, k=0, alpha=np.nan, rep=rep, rho=np.corrcoef(pi, yi)[0, 1]))
    print("null grid n =", n, f"{(time.time()-t0)/60:.1f} min", flush=True)
if not GRID_DONE:
    G = pd.DataFrame(rows); G.to_csv(os.path.join(OUT, "null_grid_raw.csv"), index=False)
    S = G.groupby(["n", "k", "alpha"], dropna=False).rho.agg(["mean", "std", "count"]).reset_index()
    S.to_csv(os.path.join(OUT, "null_grid_summary.csv"), index=False)
else:
    S = pd.read_csv(os.path.join(OUT, "null_grid_summary.csv")); print("null grid loaded from file (already complete)", flush=True)

# ---- calibration and 5-fold, k=5, alpha=1; null and signal
rows = []
for n in [15, 20, 27, 40, 60]:
    for r_true in [0.0, 0.2, 0.4]:
        for rep in range(200):
            X = rng.standard_normal((n, 5000)); z = rng.standard_normal(n)
            X[:, :50] += z[:, None]                                 # latent factor carried by 50 features
            e = rng.standard_normal(n); y = r_true * z + np.sqrt(1 - r_true ** 2) * e
            K = X @ X.T; Y = np.vstack([y] + [rng.permutation(y) for _ in range(100)])
            for scheme in ["loo", "5fold"]:
                F = make_folds(K, n, 5, scheme); rho = spear_rows(predict(F, Y, 1.0), Y)
                rows.append(dict(n=n, true_r=r_true, scheme=scheme, rep=rep, raw=rho[0], null_mean=rho[1:].mean(),
                                 calibrated=rho[0] - rho[1:].mean(), oracle=np.corrcoef(rank(z[None])[0], rank(y[None])[0])[0, 1]))
        print("calibration n =", n, f"{(time.time()-t0)/60:.1f} min", flush=True)
C = pd.DataFrame(rows); C.to_csv(os.path.join(OUT, "calibration_raw.csv"), index=False)
CS = C.groupby(["n", "true_r", "scheme"])[["raw", "null_mean", "calibrated", "oracle"]].mean().reset_index()
CS.to_csv(os.path.join(OUT, "calibration_summary.csv"), index=False)
print(S[(S.k.isin([0, 5])) & (S.alpha.isin([1]) | S.k.eq(0))].to_string(index=False))
print(CS.round(3).to_string(index=False))

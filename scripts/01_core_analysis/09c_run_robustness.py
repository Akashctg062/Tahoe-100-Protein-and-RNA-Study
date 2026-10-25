"""Step 9c - pre-specified robustness grid (PCA k = 3, 10; ridge alpha = 0.1, 10), run with the
null-calibrated method of Amendment 1 (09b) and 200 permutations each. Run from the scripts folder:
    python 09c_run_robustness.py      Output: ../results/confirmatory/robustness/
"""
import importlib, sys, os
sys.path.insert(0, os.getcwd())
m = importlib.import_module("09b_amendment_calibrated_C3")
for k, a in [(3, 1.0), (10, 1.0), (5, 0.1), (5, 10.0)]:
    print(f"##### robustness k={k} alpha={a}", flush=True)
    m.main(n_perm=200, k=k, alpha=a, out=f"../results/confirmatory/robustness/k{k}_a{a}")

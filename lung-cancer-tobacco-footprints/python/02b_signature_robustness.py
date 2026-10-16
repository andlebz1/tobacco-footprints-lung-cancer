"""
Step 2b - How robust is the de novo tobacco signature?

Checks requested for any NMF signature result:
  1. Does a SBS4-like signature appear for different numbers of signatures (k = 2..6)?
  2. Is it stable across random starting points (20 random initialisations per k)?
  3. How well does each k reconstruct the data (reconstruction error, cosine)?
  4. How well does the COSMIC refit (step 2) reconstruct each tumour?

Run:  python python/02b_signature_robustness.py
"""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF

ROOT = Path(__file__).resolve().parents[1]
PROC, FIG, TAB = ROOT / "data/processed", ROOT / "results/figures", ROOT / "results/tables"

mat = pd.read_csv(PROC / "mutation_matrix_96.csv", index_col=0)          # 96 x N
cosmic = pd.read_csv(ROOT / "data/reference/COSMIC_v3.2_SBS_GRCh37.txt", sep="\t")
cosmic.index = [f"{s[0]}[{t}]{s[2]}" for t, s in zip(cosmic["Type"], cosmic["Subtype"])]
cosmic = cosmic.drop(columns=["Type", "Subtype"]).loc[mat.index]
sbs4 = cosmic["SBS4"].values
X = mat.T.values.astype(float)                                           # N x 96


def cos(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


rows = []
for k in range(2, 7):
    sbs4_like = []
    for seed in range(20):
        m = NMF(n_components=k, init="random", max_iter=3000, random_state=seed)
        H = m.fit_transform(X)
        W = m.components_.T
        W = W / W.sum(axis=0)
        recon = H @ m.components_
        sims = [cos(W[:, j], sbs4) for j in range(k)]
        best = int(np.argmax(sims))
        sbs4_like.append(W[:, best])
        per_tumour = [cos(X[i], recon[i]) for i in range(len(X)) if X[i].sum() > 0]
        rows.append({"k": k, "seed": seed, "best_cosine_to_SBS4": sims[best],
                     "reconstruction_error": m.reconstruction_err_,
                     "median_tumour_reconstruction_cosine": float(np.median(per_tumour))})
    # stability: mean pairwise cosine between the SBS4-like signatures of different seeds
    pair = [cos(a, b) for a, b in itertools.combinations(sbs4_like, 2)]
    for r in rows:
        if r["k"] == k:
            r["stability_mean_pairwise_cosine"] = float(np.mean(pair))

res = pd.DataFrame(rows)
summary = (res.groupby("k")
             .agg(SBS4_cosine_median=("best_cosine_to_SBS4", "median"),
                  SBS4_cosine_min=("best_cosine_to_SBS4", "min"),
                  stability_across_seeds=("stability_mean_pairwise_cosine", "first"),
                  reconstruction_error=("reconstruction_error", "median"),
                  median_tumour_reconstruction_cosine=("median_tumour_reconstruction_cosine", "median"))
             .round(3).reset_index())
summary.to_csv(TAB / "signature_robustness_nmf.csv", index=False)
print("NMF robustness (20 random initialisations per k):")
print(summary.to_string(index=False))

# refit quality from step 2
expo = pd.read_csv(TAB / "signature_exposures_per_tumour.csv")
ok = expo[expo.n_snv_in_matrix >= 30]
print(f"\nCOSMIC refit (7 signatures), tumours with >= 30 SNVs: n = {len(ok)}")
print(f"  median fit cosine = {ok.fit_cosine.median():.3f}; "
      f"tumours with fit cosine >= 0.90: {(ok.fit_cosine >= 0.90).mean()*100:.0f}%")

fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
axes[0].boxplot([res.loc[res.k == k, "best_cosine_to_SBS4"] for k in range(2, 7)],
                tick_labels=list(range(2, 7)))
axes[0].set(title="SBS4-like signature found?", xlabel="number of signatures (k)",
            ylabel="cosine to COSMIC SBS4", ylim=(0.5, 1.0))
axes[1].plot(summary.k, summary.reconstruction_error, "o-", color="black")
axes[1].set(title="Reconstruction error", xlabel="k", ylabel="Frobenius error")
axes[2].hist(ok.fit_cosine, bins=25, color="grey")
axes[2].set(title="COSMIC refit quality per tumour", xlabel="cosine (observed vs refit)",
            ylabel="tumours")
fig.tight_layout()
fig.savefig(FIG / "09_signature_robustness.png", dpi=150)
print("figure saved: results/figures/09_signature_robustness.png")

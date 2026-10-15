"""
Step 2 - Mutational signatures: can we "see" tobacco in the tumour DNA?

Idea: every mutagen leaves its own pattern of mutations (a "signature").
Tobacco smoke (benzo[a]pyrene) mainly causes C>A changes = COSMIC SBS4.
We:
  a) plot the average 96-channel spectrum of never vs ever smokers
  b) learn signatures from the data WITHOUT labels (NMF, de novo)
  c) fit known COSMIC signatures to each tumour (NNLS, refitting)
  d) test whether SBS4 is higher in smokers and rises with pack-years

Run:  python python/02_mutational_signatures.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import nnls
from scipy import stats
from sklearn.decomposition import NMF

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

mat = pd.read_csv(PROC / "mutation_matrix_96.csv", index_col=0)      # 96 x N
meta = pd.read_csv(PROC / "sample_metrics.csv").set_index("PATIENT_ID")
meta = meta.loc[mat.columns]

cosmic = pd.read_csv(ROOT / "data/reference/COSMIC_v3.2_SBS_GRCh37.txt", sep="\t")
cosmic.index = [f"{s[0]}[{t}]{s[2]}" for t, s in zip(cosmic["Type"], cosmic["Subtype"])]
cosmic = cosmic.drop(columns=["Type", "Subtype"]).loc[mat.index]

SUB_COLOURS = {"C>A": "#1EBFF0", "C>G": "#050708", "C>T": "#E62725",
               "T>A": "#CBCACB", "T>C": "#A1CF64", "T>G": "#EDC8C5"}
sub_of_channel = [c[2:5] for c in mat.index]
bar_colours = [SUB_COLOURS[s] for s in sub_of_channel]


def cosine(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def plot_spectrum(ax, profile, title):
    ax.bar(range(96), profile, color=bar_colours, width=0.8)
    ax.set_xlim(-1, 96)
    ax.set_xticks([])
    ax.set_ylabel("fraction of SNVs")
    ax.set_title(title, fontsize=10, loc="left")
    for i, s in enumerate(SUB_COLOURS):
        ax.text(i * 16 + 8, ax.get_ylim()[1] * 0.92, s, ha="center",
                fontsize=8, color="black")


# ---------------------------------------------------------------------------
# a) average spectrum per group
# ---------------------------------------------------------------------------
prop = mat / mat.sum(axis=0)                        # each tumour sums to 1
never = prop.loc[:, meta.ever_smoker == "Never"].mean(axis=1)
ever = prop.loc[:, meta.ever_smoker == "Ever"].mean(axis=1)

fig, axes = plt.subplots(3, 1, figsize=(11, 7), sharex=True)
plot_spectrum(axes[0], never, f"Never-smokers (n={int((meta.ever_smoker=='Never').sum())})")
plot_spectrum(axes[1], ever, f"Ever-smokers (n={int((meta.ever_smoker=='Ever').sum())})")
plot_spectrum(axes[2], cosmic["SBS4"], "COSMIC SBS4 (tobacco smoking) reference")
fig.tight_layout()
fig.savefig(FIG / "01_mutation_spectrum_never_vs_ever.png", dpi=150)
plt.close(fig)

print(f"cosine(ever-smoker spectrum, SBS4)  = {cosine(ever, cosmic.SBS4):.3f}")
print(f"cosine(never-smoker spectrum, SBS4) = {cosine(never, cosmic.SBS4):.3f}")

# ---------------------------------------------------------------------------
# b) de novo extraction with NMF (no smoking labels used!)
# ---------------------------------------------------------------------------
# NMF splits the matrix  M (96 x N)  into  W (96 x k) x H (k x N)
# W = signatures, H = how much each tumour uses each signature.
K = 3
model = NMF(n_components=K, init="nndsvda", max_iter=2000, random_state=2026)
H = model.fit_transform(mat.T.values)               # N x k  (exposures)
W = model.components_.T                             # 96 x k (signatures)
W = W / W.sum(axis=0)

denovo_rows = []
for j in range(K):
    sims = {s: cosine(W[:, j], cosmic[s]) for s in cosmic.columns}
    best = max(sims, key=sims.get)
    denovo_rows.append({"denovo_signature": f"Sig_{chr(65 + j)}",
                        "best_COSMIC_match": best,
                        "cosine_similarity": round(sims[best], 3)})
denovo = pd.DataFrame(denovo_rows)
denovo.to_csv(TAB / "denovo_signatures_vs_cosmic.csv", index=False)
print("\nDe novo signatures (NMF, k=3):")
print(denovo.to_string(index=False))

fig, axes = plt.subplots(K, 1, figsize=(11, 2.3 * K), sharex=True)
for j in range(K):
    r = denovo.iloc[j]
    plot_spectrum(axes[j], W[:, j],
                  f"{r.denovo_signature}  ->  closest COSMIC: {r.best_COSMIC_match} "
                  f"(cosine {r.cosine_similarity})")
fig.tight_layout()
fig.savefig(FIG / "02_denovo_signatures_nmf.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------------
# c) refit known COSMIC signatures to every tumour (NNLS)
# ---------------------------------------------------------------------------
# Signatures known to be active in lung adenocarcinoma (COSMIC / PCAWG):
#   SBS1, SBS5  clock-like (age)          SBS2, SBS13  APOBEC enzymes
#   SBS4        tobacco smoking           SBS18        reactive oxygen
#   SBS40       flat, unknown cause
SIGS = ["SBS1", "SBS2", "SBS4", "SBS5", "SBS13", "SBS18", "SBS40"]
S = cosmic[SIGS].values

rows = []
for pid in mat.columns:
    counts = mat[pid].values.astype(float)
    expo, _ = nnls(S, counts)                        # non-negative fit
    recon = S @ expo
    rows.append({"PATIENT_ID": pid,
                 **{s: e for s, e in zip(SIGS, expo)},
                 "n_snv_in_matrix": counts.sum(),
                 "fit_cosine": cosine(counts, recon) if counts.sum() > 0 else np.nan})
expo = pd.DataFrame(rows).set_index("PATIENT_ID")
rel = expo[SIGS].div(expo[SIGS].sum(axis=1), axis=0).add_suffix("_rel")
expo = expo.join(rel).join(meta[["smoking_group", "ever_smoker", "pack_years",
                                  "TMB_per_Mb", "frac_C>A"]])
expo.to_csv(TAB / "signature_exposures_per_tumour.csv")

# tumours with very few mutations give unstable fits -> flag them
MIN_SNV = 30
ok = expo[expo.n_snv_in_matrix >= MIN_SNV]
print(f"\ntumours with >= {MIN_SNV} SNVs used for statistics: {len(ok)} / {len(expo)}")

# ---------------------------------------------------------------------------
# d) statistics
# ---------------------------------------------------------------------------
results = []

def mw(var, data, label):
    a = data.loc[data.ever_smoker == "Ever", var].dropna()
    b = data.loc[data.ever_smoker == "Never", var].dropna()
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    results.append({"test": f"Mann-Whitney: {label}, ever vs never",
                    "n_ever": len(a), "n_never": len(b),
                    "median_ever": round(a.median(), 3),
                    "median_never": round(b.median(), 3),
                    "statistic": round(u, 1), "p_value": p})

mw("TMB_per_Mb", expo, "TMB (mut/Mb)")
mw("frac_C>A", expo, "C>A fraction")
mw("SBS4_rel", ok, "SBS4 relative contribution")

pk = ok.dropna(subset=["pack_years"])
rho, p = stats.spearmanr(pk["pack_years"], pk["SBS4_rel"])
results.append({"test": "Spearman: pack-years vs SBS4_rel (smokers with data)",
                "n_ever": len(pk), "n_never": np.nan, "median_ever": np.nan,
                "median_never": np.nan, "statistic": round(rho, 3), "p_value": p})

# 3 groups: Never / Former / Current
groups = [ok.loc[ok.smoking_group == g, "SBS4_rel"] for g in ["Never", "Former", "Current"]]
h, p = stats.kruskal(*groups)
results.append({"test": "Kruskal-Wallis: SBS4_rel across Never/Former/Current",
                "n_ever": sum(len(g) for g in groups), "n_never": np.nan,
                "median_ever": np.nan, "median_never": np.nan,
                "statistic": round(h, 2), "p_value": p})

res = pd.DataFrame(results)
res.to_csv(TAB / "statistics_signatures.csv", index=False)
print("\n" + res.to_string(index=False))

# ---------------------------------------------------------------------------
# figures for c) and d)
# ---------------------------------------------------------------------------
order = ["Never", "Former", "Current"]
col = {"Never": "#4C9F70", "Former": "#E0A030", "Current": "#C0392B"}

fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
rng = np.random.default_rng(1)
for ax, var, lab, data in [(axes[0], "TMB_per_Mb", "TMB (mutations / Mb, log scale)", expo),
                           (axes[1], "SBS4_rel", "SBS4 (tobacco) share of SNVs", ok)]:
    vals = [data.loc[data.smoking_group == g, var].dropna() for g in order]
    ax.boxplot(vals, showfliers=False, widths=0.5)
    for i, (g, v) in enumerate(zip(order, vals), start=1):
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=9, alpha=0.6, color=col[g])
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([f"{g}\n(n={len(v)})" for g, v in zip(order, vals)])
    ax.set_ylabel(lab)
    if var == "TMB_per_Mb":
        ax.set_yscale("log")

axes[2].scatter(pk.pack_years, pk.SBS4_rel, s=12, alpha=0.7, color="#C0392B")
axes[2].set_xlabel("pack-years")
axes[2].set_ylabel("SBS4 share of SNVs")
axes[2].set_title(f"Spearman rho = {rho:.2f}", fontsize=10)
fig.suptitle("Tobacco footprint in TCGA lung adenocarcinoma", fontsize=12)
fig.tight_layout()
fig.savefig(FIG / "03_tmb_sbs4_by_smoking.png", dpi=150)
plt.close(fig)

# stacked bar of signature shares, tumours sorted by SBS4
srt = ok.sort_values("SBS4_rel")
fig, ax = plt.subplots(figsize=(12, 3.8))
bottom = np.zeros(len(srt))
palette = ["#8E8E8E", "#F4A261", "#C0392B", "#BDBDBD", "#E9C46A", "#457B9D", "#A8DADC"]
for s, c in zip(SIGS, palette):
    ax.bar(range(len(srt)), srt[f"{s}_rel"], bottom=bottom, color=c, width=1, label=s)
    bottom += srt[f"{s}_rel"].values
strip = [col[g] for g in srt.smoking_group]
ax.scatter(range(len(srt)), np.full(len(srt), -0.05), c=strip, marker="s", s=10, clip_on=False)
ax.set_xlim(-0.5, len(srt) - 0.5)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_ylabel("share of SNVs")
ax.set_xlabel("tumours (sorted by SBS4); coloured strip = smoking group "
              "(green never, orange former, red current)", labelpad=14)
ax.legend(ncol=7, fontsize=8, loc="upper left", bbox_to_anchor=(0, 1.15))
fig.tight_layout()
fig.savefig(FIG / "04_signature_contributions.png", dpi=150)
plt.close(fig)
print("\nfigures saved in results/figures/")

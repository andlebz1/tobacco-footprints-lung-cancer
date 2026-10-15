"""
Step 3 - Driver genes and genome instability.

Questions:
  a) Do smokers and never-smokers have different DRIVER mutations?
     (EGFR is the classic never-smoker driver; KRAS G12C is linked to
      smoking and is now a drug target - sotorasib / adagrasib.)
  b) Is the genome more UNSTABLE (more copy-number changes) when TP53
     is mutated, or in smokers?

Run:  python python/03_drivers_and_instability.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"

meta = pd.read_csv(PROC / "sample_metrics.csv").set_index("PATIENT_ID")

# ---------------------------------------------------------------------------
# a) driver mutations
# ---------------------------------------------------------------------------
maf = pd.read_csv(RAW / "data_mutations.txt", sep="\t", comment="#",
                  usecols=["Hugo_Symbol", "Tumor_Sample_Barcode",
                           "Variant_Classification", "HGVSp_Short"],
                  low_memory=False)
maf["PATIENT_ID"] = maf["Tumor_Sample_Barcode"].str[:12]
maf = maf[maf.PATIENT_ID.isin(meta.index)]

# protein-changing mutations only (silent ones do not change the protein)
PROTEIN_CHANGING = ["Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
                    "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del",
                    "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site"]
pc = maf[maf.Variant_Classification.isin(PROTEIN_CHANGING)]

GENES = ["TP53", "KRAS", "EGFR", "STK11", "KEAP1", "BRAF", "NF1", "SMARCA4"]
status = pd.DataFrame(False, index=meta.index, columns=GENES)
for g in GENES:
    status.loc[status.index.isin(pc.loc[pc.Hugo_Symbol == g, "PATIENT_ID"]), g] = True

# KRAS G12C: the specific variant that new KRAS drugs target
g12c = pc[(pc.Hugo_Symbol == "KRAS") & (pc.HGVSp_Short == "p.G12C")].PATIENT_ID
status["KRAS_G12C"] = status.index.isin(g12c)
status.to_csv(PROC / "driver_status.csv")

rows = []
for g in status.columns:
    ever = status.loc[meta.ever_smoker == "Ever", g]
    never = status.loc[meta.ever_smoker == "Never", g]
    table = [[ever.sum(), (~ever).sum()], [never.sum(), (~never).sum()]]
    odds, p = stats.fisher_exact(table)
    rows.append({"gene": g, "mut_ever": int(ever.sum()), "n_ever": len(ever),
                 "pct_ever": round(100 * ever.mean(), 1),
                 "mut_never": int(never.sum()), "n_never": len(never),
                 "pct_never": round(100 * never.mean(), 1),
                 "odds_ratio": round(odds, 2), "p_value": p})
drv = pd.DataFrame(rows)

# many genes tested -> correct p-values (Benjamini-Hochberg FDR)
def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out

drv["q_value_BH"] = bh(drv.p_value)
drv.to_csv(TAB / "driver_genes_ever_vs_never.csv", index=False)
print(drv.to_string(index=False))

fig, ax = plt.subplots(figsize=(9, 4))
x = np.arange(len(drv))
ax.bar(x - 0.2, drv.pct_never, 0.4, label=f"Never (n={drv.n_never[0]})", color="#4C9F70")
ax.bar(x + 0.2, drv.pct_ever, 0.4, label=f"Ever (n={drv.n_ever[0]})", color="#C0392B")
for i, (p, q) in enumerate(zip(drv.p_value, drv.q_value_BH)):
    mark = "*" if q < 0.05 else ("+" if p < 0.05 else "")
    ax.text(i, max(drv.pct_never[i], drv.pct_ever[i]) + 2, mark, ha="center", fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(drv.gene, rotation=30)
ax.set_ylabel("% tumours mutated")
ax.set_title("Driver mutations by smoking history (+ = p < 0.05 before correction; * = FDR < 0.05)", fontsize=10)
ax.legend()
fig.tight_layout()
fig.savefig(FIG / "05_driver_genes_by_smoking.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------------
# b) genome instability from copy-number data
# ---------------------------------------------------------------------------
# data_cna.txt: genes x tumours, GISTIC calls
#   -2 deep deletion, -1 loss, 0 normal, 1 gain, 2 amplification
cna = pd.read_csv(RAW / "data_cna.txt", sep="\t", index_col=0)
cna = cna.drop(columns=[c for c in cna.columns if not c.startswith("TCGA")])
cna.columns = [c[:12] for c in cna.columns]
cna = cna.loc[:, ~cna.columns.duplicated()]

# "fraction of genome altered" (gene-based approximation)
fga = (cna != 0).mean(axis=0).rename("fraction_genes_altered")
deep = ((cna == 2) | (cna == -2)).mean(axis=0).rename("fraction_genes_deep_altered")
inst = pd.concat([fga, deep], axis=1).join(status, how="inner").join(
    meta[["ever_smoker", "smoking_group", "TMB_per_Mb"]], how="inner")
inst.to_csv(TAB / "genome_instability_per_tumour.csv")
print(f"\ntumours with copy-number + mutation + smoking data: {len(inst)}")

tests = []
def mw(var, grp_col, a_lab, b_lab):
    a = inst.loc[inst[grp_col] == a_lab, var]
    b = inst.loc[inst[grp_col] == b_lab, var]
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    tests.append({"variable": var, "comparison": f"{grp_col}: {a_lab} vs {b_lab}",
                  "n_a": len(a), "n_b": len(b), "median_a": round(a.median(), 3),
                  "median_b": round(b.median(), 3), "p_value": p})

mw("fraction_genes_altered", "TP53", True, False)
mw("fraction_genes_altered", "ever_smoker", "Ever", "Never")
rho, p = stats.spearmanr(inst.fraction_genes_altered, inst.TMB_per_Mb)
tests.append({"variable": "fraction_genes_altered vs TMB", "comparison": "Spearman rho",
              "n_a": len(inst), "n_b": np.nan, "median_a": round(rho, 3),
              "median_b": np.nan, "p_value": p})
tst = pd.DataFrame(tests)
tst.to_csv(TAB / "statistics_instability.csv", index=False)
print(tst.to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
rng = np.random.default_rng(2)
for ax, col, labs, names in [
        (axes[0], "TP53", [False, True], ["TP53 wild-type", "TP53 mutated"]),
        (axes[1], "ever_smoker", ["Never", "Ever"], ["Never-smoker", "Ever-smoker"])]:
    vals = [inst.loc[inst[col] == l, "fraction_genes_altered"] for l in labs]
    ax.boxplot(vals, showfliers=False, widths=0.5)
    for i, v in enumerate(vals, start=1):
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=8, alpha=0.5, color="#457B9D")
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"{n}\n(n={len(v)})" for n, v in zip(names, vals)])
    ax.set_ylabel("fraction of genes with copy-number change")
p_tp53 = tst.p_value[0]
p_smk = tst.p_value[1]
axes[0].set_title(f"Mann-Whitney p = {p_tp53:.1e}", fontsize=10)
axes[1].set_title(f"Mann-Whitney p = {p_smk:.2g}", fontsize=10)
fig.suptitle("Genome instability (copy-number burden)", fontsize=12)
fig.tight_layout()
fig.savefig(FIG / "06_genome_instability.png", dpi=150)
plt.close(fig)
print("figures saved")

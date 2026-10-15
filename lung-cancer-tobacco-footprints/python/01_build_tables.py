"""
Step 1 - Clean the raw data and build three tidy tables.

Output (data/processed/):
  samples.csv          one row per tumour: smoking group, age, sex, pack-years
  mutation_matrix_96.csv  96 x tumours: counts of each SNV type in its
                          trinucleotide context (the input for signatures)
  sample_metrics.csv   per tumour: number of mutations, TMB, C>A fraction ...

Run:  python python/01_build_tables.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Clinical data: who smoked?
# ---------------------------------------------------------------------------
# TCGA code for TOBACCO_SMOKING_HISTORY_INDICATOR:
#   1 = lifelong non-smoker
#   2 = current smoker
#   3 = former smoker, quit more than 15 years ago
#   4 = former smoker, quit 15 years ago or less
#   5 = former smoker, duration not specified
clin = pd.read_csv(RAW / "data_clinical_patient.txt", sep="\t",
                   comment="#", low_memory=False)

smoke_map = {"1": "Never", "2": "Current",
             "3": "Former", "4": "Former", "5": "Former"}
clin["smoking_group"] = clin["TOBACCO_SMOKING_HISTORY_INDICATOR"].astype(str).map(smoke_map)
clin["ever_smoker"] = np.where(clin["smoking_group"] == "Never", "Never", "Ever")
clin.loc[clin["smoking_group"].isna(), "ever_smoker"] = np.nan
clin["pack_years"] = pd.to_numeric(clin["SMOKING_PACK_YEARS"], errors="coerce")
clin["age"] = pd.to_numeric(clin["AGE"], errors="coerce")

samples = clin[["PATIENT_ID", "smoking_group", "ever_smoker",
                "pack_years", "age", "SEX"]].copy()

# ---------------------------------------------------------------------------
# 2. Mutations (MAF file): one row = one mutation in one tumour
# ---------------------------------------------------------------------------
cols = ["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Type",
        "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele2",
        "ref_context"]
maf = pd.read_csv(RAW / "data_mutations.txt", sep="\t", comment="#",
                  usecols=cols, low_memory=False)
maf["PATIENT_ID"] = maf["Tumor_Sample_Barcode"].str[:12]   # TCGA-XX-XXXX

# keep only patients that have smoking information
maf = maf[maf["PATIENT_ID"].isin(samples.loc[samples.smoking_group.notna(), "PATIENT_ID"])]
print(f"tumours with mutation data + smoking info: {maf.PATIENT_ID.nunique()}")

# ---------------------------------------------------------------------------
# 3. Build the 96-channel mutation matrix
# ---------------------------------------------------------------------------
# A single base substitution (SNV) is written with the pyrimidine (C or T)
# as the reference base. If the reference is a purine (G or A), we take the
# reverse complement of the whole trinucleotide. So G>T becomes C>A.
COMP = str.maketrans("ACGT", "TGCA")

def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]

snv = maf[(maf["Variant_Type"] == "SNP") & maf["ref_context"].notna()].copy()
snv = snv[snv["ref_context"].str.len() % 2 == 1]          # need a clear centre
mid = (snv["ref_context"].str.len() // 2).astype(int)    # centre position
snv["tri"] = [c[m - 1:m + 2] for c, m in zip(snv["ref_context"], mid)]
snv = snv[snv["tri"].str[1] == snv["Reference_Allele"]]   # sanity check
snv = snv[snv["tri"].str.fullmatch("[ACGT]{3}")]

ref = snv["Reference_Allele"].to_numpy()
alt = snv["Tumor_Seq_Allele2"].to_numpy()
tri = snv["tri"].to_numpy()
purine = np.isin(ref, ["A", "G"])
ref = np.where(purine, [r.translate(COMP) for r in ref], ref)
alt = np.where(purine, [a.translate(COMP) for a in alt], alt)
tri = np.where(purine, [revcomp(t) for t in tri], tri)

# label in COSMIC style, e.g. A[C>A]T
snv["channel"] = [f"{t[0]}[{r}>{a}]{t[2]}" for t, r, a in zip(tri, ref, alt)]
snv["sub6"] = [f"{r}>{a}" for r, a in zip(ref, alt)]

# order of the 96 channels exactly as in the COSMIC reference file
cosmic = pd.read_csv(ROOT / "data/reference/COSMIC_v3.2_SBS_GRCh37.txt", sep="\t")
channels = [f"{s[0]}[{t}]{s[2]}" for t, s in zip(cosmic["Type"], cosmic["Subtype"])]

mat96 = (snv.groupby(["channel", "PATIENT_ID"]).size()
            .unstack(fill_value=0)
            .reindex(channels, fill_value=0))
mat96.to_csv(OUT / "mutation_matrix_96.csv")
print(f"96-channel matrix: {mat96.shape[0]} channels x {mat96.shape[1]} tumours")

# ---------------------------------------------------------------------------
# 4. Simple per-tumour metrics
# ---------------------------------------------------------------------------
EXOME_MB = 38.0   # approximate size of the exome that was sequenced (Mb)

n_all = maf.groupby("PATIENT_ID").size().rename("n_mutations")
sub6 = snv.groupby(["PATIENT_ID", "sub6"]).size().unstack(fill_value=0)
frac6 = sub6.div(sub6.sum(axis=1), axis=0).add_prefix("frac_")

metrics = pd.concat([n_all, sub6.sum(axis=1).rename("n_snv"), frac6], axis=1)
metrics["TMB_per_Mb"] = metrics["n_mutations"] / EXOME_MB
metrics = metrics.reset_index().rename(columns={"index": "PATIENT_ID"})
metrics = metrics.merge(samples, on="PATIENT_ID", how="left")
metrics.to_csv(OUT / "sample_metrics.csv", index=False)

samples.to_csv(OUT / "samples.csv", index=False)
print(metrics.groupby("smoking_group")[["n_mutations", "frac_C>A"]].median().round(3))

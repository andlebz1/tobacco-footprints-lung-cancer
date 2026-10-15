"""
Step 0 - Download the public TCGA lung adenocarcinoma (LUAD) data.

Source: cBioPortal DataHub, study "luad_tcga" (TCGA Firehose Legacy).
All files are public, open-access, de-identified summary data.

Run:  python python/00_download_data.py
"""
from pathlib import Path
import urllib.request

BASE = ("https://media.githubusercontent.com/media/cBioPortal/datahub/"
        "master/public/luad_tcga/")

FILES = [
    "data_clinical_patient.txt",   # age, sex, smoking history ...
    "data_clinical_sample.txt",    # sample IDs
    "data_mutations.txt",          # somatic mutations (MAF format)
    "data_cna.txt",                # copy-number calls per gene (GISTIC)
    "data_methylation_hm450.txt",  # DNA methylation, Illumina 450K, per gene
]

out_dir = Path(__file__).resolve().parents[1] / "data" / "raw"
out_dir.mkdir(parents=True, exist_ok=True)

for name in FILES:
    target = out_dir / name
    if target.exists():
        print(f"already present: {name}")
        continue
    print(f"downloading {name} ...")
    urllib.request.urlretrieve(BASE + name, target)
    print(f"  saved {target.stat().st_size / 1e6:.1f} MB")

# ---------------------------------------------------------------------------
# COSMIC v3.2 SBS reference signatures (GRCh37)
# A copy is already included in data/reference/. If it is missing, it is
# downloaded from the Alexandrov lab's SigProfilerAssignment repository and
# re-formatted to the classic COSMIC layout (columns Type, Subtype, SBS...).
# ---------------------------------------------------------------------------
ref = Path(__file__).resolve().parents[1] / "data" / "reference" / "COSMIC_v3.2_SBS_GRCh37.txt"
if ref.exists():
    print("already present: COSMIC_v3.2_SBS_GRCh37.txt")
else:
    import pandas as pd
    url = ("https://raw.githubusercontent.com/AlexandrovLab/SigProfilerAssignment/main/"
           "SigProfilerAssignment/data/Reference_Signatures/GRCh37/COSMIC_v3.2_SBS_GRCh37.txt")
    print("downloading COSMIC v3.2 reference signatures ...")
    cosmic = pd.read_csv(url, sep="\t")
    channel = cosmic.pop("Type")                       # e.g. "A[C>A]T"
    cosmic.insert(0, "Subtype", [c[0] + c[2] + c[6] for c in channel])   # "ACT"
    cosmic.insert(0, "Type", [c[2:5] for c in channel])                  # "C>A"
    ref.parent.mkdir(parents=True, exist_ok=True)
    cosmic.to_csv(ref, sep="\t", index=False)
    print(f"  saved {ref}")

print("done")

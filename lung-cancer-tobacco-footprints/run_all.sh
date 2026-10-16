#!/usr/bin/env bash
# Run the whole pipeline from the project folder:  bash run_all.sh
set -e
python python/00_download_data.py
python python/01_build_tables.py
python python/02_mutational_signatures.py
python python/02b_signature_robustness.py
python python/03_drivers_and_instability.py
Rscript R/04_methylation_smoking.R
echo "All steps finished. See results/figures and results/tables."

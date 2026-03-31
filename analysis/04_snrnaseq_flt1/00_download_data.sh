#!/usr/bin/env bash
# Download SEA-AD Microglia-Immune h5ad from S3 to HPC
# Run on Minerva login node (not via LSF -- just a download)
set -euo pipefail

DATA_DIR="/sc/arion/projects/vascbrain/andres/vasc-sflt1-alphafold/data/sea_ad"
S3_BUCKET="s3://sea-ad-single-cell-profiling/Microglia-and-Immune-for-AAIC"
H5AD_FILE="SEA-AD_Microglia-and-Immune_multi-regional_final-nuclei_AAIC-pre-release.2025-07-24.h5ad"
README_FILE="README.md"

mkdir -p "${DATA_DIR}"

module load awscli/2.17.37

# Download README first (small)
if [[ ! -f "${DATA_DIR}/${README_FILE}" ]]; then
    echo "Downloading README..."
    aws s3 cp --no-sign-request "${S3_BUCKET}/${README_FILE}" "${DATA_DIR}/${README_FILE}"
fi

# Download h5ad (~3 GB)
if [[ ! -f "${DATA_DIR}/${H5AD_FILE}" ]]; then
    echo "Downloading h5ad (~3 GB)..."
    aws s3 cp --no-sign-request "${S3_BUCKET}/${H5AD_FILE}" "${DATA_DIR}/${H5AD_FILE}"
    echo "Download complete: $(ls -lh "${DATA_DIR}/${H5AD_FILE}")"
else
    echo "h5ad already exists: $(ls -lh "${DATA_DIR}/${H5AD_FILE}")"
fi

echo "Done. Files in ${DATA_DIR}:"
ls -lh "${DATA_DIR}/"

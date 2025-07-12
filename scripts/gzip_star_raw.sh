#!/bin/bash

# gzip_star_raw.sh
# Compress barcodes.tsv, features.tsv, and matrix.mtx in STARsolo 'raw' folders
# Usage: ./gzip_star_raw.sh /absolute/path/to/main_work sample1 sample2 ...

# Example: ./gzip_star_raw.sh /export/storage/users/azaid/vallabh_lab SRR16502301 SRR16502302

BASE="$1"
shift
samples=("$@")

for sample in "${samples[@]}"; do
  raw_dir="${BASE}/${sample}/starsolo/Solo.out/Gene/raw"
  if [ -d "$raw_dir" ]; then
    echo "Processing $raw_dir"
    for fname in barcodes.tsv features.tsv matrix.mtx; do
      if [ -f "$raw_dir/$fname" ] && [ ! -f "$raw_dir/$fname.gz" ]; then
        echo "Compressing $raw_dir/$fname"
        gzip "$raw_dir/$fname"
      fi
    done
  else
    echo "Directory $raw_dir does not exist, skipping."
  fi
done

#!/bin/bash
#$ -N soupx_SRR16502308
#$ -o /export/storage/users/azaid/vallabh_lab/work/SRR16502308/soupx_matrix/stdout.log
#$ -e /export/storage/users/azaid/vallabh_lab/work/SRR16502308/soupx_matrix/stderr.log
#$ -l h_vmem=80G
#$ -pe smp 8

# Activate environment
source ~/.bashrc
mamba activate aav

# Run the R script for this specific sample
# The paths are passed in from the generator script
Rscript "/home/azaid/vallabh_lab/scripts/soupx_filter.R" \
  --base_path "/export/storage/users/azaid/vallabh_lab/work" \
  --samples "SRR16502308"


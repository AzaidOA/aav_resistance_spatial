#!/bin/bash

# This script generates a qsub job script for a single sample.
#
# Arguments:
#   $1: Sample ID (e.g., SRR16502301)
#   $2: Base directory for scripts (e.g., /home/azaid/vallabh_lab/scripts)
#   $3: Full path to the R script to execute (e.g., /home/azaid/vallabh_lab/scripts/soupx_filter.R)
#   $4: Base path for the work/data directories (e.g., /export/storage/users/azaid/vallabh_lab/work)

# Validate Arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <sample_id> <scripts_base_dir> <r_script_path> <work_base_path>"
    exit 1
fi

# Assign Arguments to Variables
sample_id="$1"
scripts_base_dir="$2"
r_script_path="$3"
work_base_path="$4"

# Define Output Paths
# Directory to save the generated job script for this specific sample
output_dir="${scripts_base_dir}/${sample_id}"
# Full path for the new job script
output_script="${output_dir}/run_soupx_${sample_id}.sh"
# Base path for the work directory
work_dir=${work_base_path}/${sample_id}/soupx_matrix
# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Create the work directory if it doesn't exist
mkdir -p "$work_dir"

# Create the Job Script using a Here Document (cat << EOF)
cat << EOF > "$output_script"
#!/bin/bash
#$ -N soupx_${sample_id}
#$ -o ${work_dir}/stdout.log
#$ -e ${work_dir}/stderr.log
#$ -l h_vmem=80G
#$ -pe smp 8

# Activate environment
source ~/.bashrc
mamba activate aav

# Run the R script for this specific sample
# The paths are passed in from the generator script
Rscript "${r_script_path}" \\
  --base_path "${work_base_path}" \\
  --samples "${sample_id}"

EOF

# Make the generated script executable
chmod +x "$output_script"

echo "Created and configured job script: ${output_script}"
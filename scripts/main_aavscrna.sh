# I need to correct and state how was the data structure created and set up
cd /home/azaid/vallabh_lab/work
STAR --runMode genomeGenerate \
  --genomeDir /export/storage/users/azaid/reference/GRCm38_index \
  --genomeFastaFiles /home/azaid/vallabh_lab/data/reference/GRCm38.fa \
  --sjdbGTFfile /home/azaid/vallabh_lab/data/reference/gencode.vM17.gtf \
  --runThreadN 8

ln -s /export/storage/users/azaid/reference/GRCm38_index /home/azaid/vallabh_lab/data/reference/GRCm38_index
ln -s /export/storage/users/azaid/bin/cellranger-9.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt

# Until here
# Main working and data directories
main_data="/home/azaid/vallabh_lab"
main_work="/export/storage/users/azaid/vallabh_lab"
scripts_main="${main_data}/scripts"
aavsc_metadata="${main_data}/data/aav_scrna/metadata/AAVscRNA_run_metadata.csv"

chmod +x "${scripts_main}/run_starsolo.sh"

# Extract unique samples (change $1 if run is not the first column)
samples=$(awk -F',' 'NR>1 {print $1}' "$aavsc_metadata" | sort | uniq)

for sample in $samples; do
  # Generate the job script for each sample
  "${scripts_main}/run_starsolo.sh" "$sample" "$main_data" "$main_work" "$scripts_main"
  # Submit the generated job
  qsub "${scripts_main}/${sample}/run_starsolo_${sample}.sh"
done

# Step 1: Replicate the directory structure from the source to the destination
# This ensures all subfolders exist in the destination before linking files
find /export/storage/users/azaid/vallabh_lab/work -type d | while read dir; do
  # Remove the source prefix and add the destination prefix to get the new path
  mkdir -p "/home/azaid/vallabh_lab${dir#/export/storage/users/azaid/vallabh_lab}"
done

# Step 2: Create symbolic links for all files from the source to the destination
# The link in the destination points to the original file in the source location
find /export/storage/users/azaid/vallabh_lab/work -type f | while read file; do
  # Remove the source prefix and add the destination prefix to get the new link path
  ln -sf "$file" "/home/azaid/vallabh_lab${file#/export/storage/users/azaid/vallabh_lab}"
done

# Step 3: Filter raw STARsolo data using the filter_droplets.R script.
Rscript "${scripts_main}/filter_droplets.R" \
  -i "${main_work}" \
  -o "${main_work}" \
  -s SRR16502301,SRR16502302,SRR16502308,SRR16502309 \
  --lower_prop 0.05 \
  --seed 100

# Step 4: Run SoupX filtering to remove ambient RNA from STARsolo outputs
# Compress STARsolo raw matrices (.tsv/.mtx) for specified samples using the gzip_star_raw.sh script
bash "${scripts_main}/gzip_star_raw.sh" "$main_work/work" SRR16502301 SRR16502302 SRR16502308 SRR16502309
# Run SoupX filter script
chmod +x "${scripts_main}/run_soupx.sh"
for sample in $samples; do
  echo "Processing sample: ${sample}"
  # 1. Create the specific job script by calling run_soupx.sh
  bash "${scripts_main}/run_soupx.sh" \
    "${sample}" \
    "${scripts_main}" \
    "${scripts_main}/soupx_filter.R" \
    "${main_work}/work"
  
  # 2. Submit the job script that was just created
  qsub "${scripts_main}/${sample}/run_soupx_${sample}.sh"
  
  echo "Submitted job for sample ${sample}."
done

# Step 5: Run scds filtering to remove doublets.
source ~/.bashrc
mamba activate scds_env
Rscript "${scripts_main}/scds_doublets.R" \
  --base_path "${main_work}/work" \
  --fig_dir "${main_data}/work" \
  --samples "$(echo "$samples" | head -n 2 | paste -s -d ',')"
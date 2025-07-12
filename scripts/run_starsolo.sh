#!/bin/bash

# Arguments
sample="$1"
main_data="$2"
main_work="$3"
scripts_main="$4"

# Paths
sample_dir="${scripts_main}/${sample}"
job_script="${sample_dir}/run_starsolo_${sample}.sh"
out_dir="${main_work}/work/${sample}/starsolo"
read_file="${main_data}/data/aav_scrna/${sample}/${sample}_S1_L001_R2_001.fastq"
genome_dir="${main_data}/data/reference/GRCm38_index"
whitelist="${main_data}/data/reference/3M-february-2018_TRU.txt"

# Create necessary folders
mkdir -p "$sample_dir"
mkdir -p "$out_dir"

# Generate job script
cat <<EOF > "$job_script"
#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=4G
#$ -N starsolo_${sample}
#$ -o ${out_dir}/stdout.log
#$ -e ${out_dir}/stderr.log

# Activate conda environment
source ~/.bashrc
mamba activate aav

# Run STARsolo
STAR --runThreadN 8 \
  --genomeDir ${genome_dir} \
  --readFilesIn ${main_data}/data/aav_scrna/${sample}/${sample}_S1_L001_R2_001.fastq ${main_data}/data/aav_scrna/${sample}/${sample}_S1_L001_R1_001.fastq \
  --readFilesCommand cat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist ${whitelist} \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --soloUMIstart 17 \
  --soloUMIlen 12 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes CR UR CY UY CB UB \
  --outFileNamePrefix ${out_dir}/
EOF

# Make job script executable
chmod +x "$job_script"

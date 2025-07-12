#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=4G
#$ -N starsolo_SRR16502308
#$ -o /export/storage/users/azaid/vallabh_lab/work/SRR16502308/starsolo/stdout.log
#$ -e /export/storage/users/azaid/vallabh_lab/work/SRR16502308/starsolo/stderr.log

# Activate conda environment
source ~/.bashrc
mamba activate aav

# Run STARsolo
STAR --runThreadN 8   --genomeDir /home/azaid/vallabh_lab/data/reference/GRCm38_index   --readFilesIn /home/azaid/vallabh_lab/data/aav_scrna/SRR16502308/SRR16502308_S1_L001_R2_001.fastq /home/azaid/vallabh_lab/data/aav_scrna/SRR16502308/SRR16502308_S1_L001_R1_001.fastq   --readFilesCommand cat   --soloType CB_UMI_Simple   --soloCBwhitelist /home/azaid/vallabh_lab/data/reference/3M-february-2018_TRU.txt   --soloCBstart 1   --soloCBlen 16   --soloUMIstart 17   --soloUMIlen 12   --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts   --soloUMIfiltering MultiGeneUMI_CR   --soloUMIdedup 1MM_CR   --clipAdapterType CellRanger4   --outFilterScoreMin 30   --outSAMtype BAM SortedByCoordinate   --outSAMattributes CR UR CY UY CB UB   --outFileNamePrefix /export/storage/users/azaid/vallabh_lab/work/SRR16502308/starsolo/

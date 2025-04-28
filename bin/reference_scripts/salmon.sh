#!/bin/bash

#SBATCH --time=4:00:00                # Maximum runtime
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=16                   # Number of threads per task
#SBATCH --mem=64gb                    # Memory per task
#SBATCH --job-name=Salmon_alignment   # Job name
#SBATCH --account=st-bstefans-1
#SBATCH --mail-user=mtello@student.ubc.ca

################################################################################


source activate fastqc_env

# Directories
GenomeDir="/scratch/st-bstefans-1/mtello/StewardCollab/mouse_transcriptome"
ReadsDir="/scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads"
OutputDir="/scratch/st-bstefans-1/mtello/StewardCollab/Salmon_aligned"

# Create transcriptome index
# cd "$GenomeDir"
# salmon index -t Mus_musculus.GRCm39.cdna.all.fa.gz -i mmusculus_index

# Create output directory if it doesn't exist
mkdir -p "$OutputDir"
cd "$OutputDir"

# Create a space-separated list of the matching files
# Reads1=$(ls "$ReadsDir"/*R1_trimmed.fastq | tr '\n' ' ')
# Reads2=$(ls "$ReadsDir"/*R2_trimmed.fastq | tr '\n' ' ')

salmon quant -i $GenomeDir/mmusculus_index -l IU \
         -1 /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF2-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF3-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF4-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF5-2-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF8-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF9-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S10-SJ-301123_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S3-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S4-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S5-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S6-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S8-SJ-010824_L1-S2_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S8-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V1-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V10-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V2-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V3-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V4-SJ-010824_L1_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V7-SJ-010824_L1-S3_R1_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V7-SJ-010824_L1_R1_trimmed.fastq \
         -2 /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF2-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF3-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF4-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF5-2-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF8-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/PF9-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S10-SJ-301123_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S3-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S4-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S5-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S6-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S8-SJ-010824_L1-S2_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/S8-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V1-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V10-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V2-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V3-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V4-SJ-010824_L1_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V7-SJ-010824_L1-S3_R2_trimmed.fastq /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/V7-SJ-010824_L1_R2_trimmed.fastq \
         -p 16 --validateMappings -o mmusculus_quant_trimmed

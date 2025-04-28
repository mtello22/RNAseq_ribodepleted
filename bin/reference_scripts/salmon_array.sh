#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --job-name=Salmon_alignment_array
#SBATCH --account=st-bstefans-1
#SBATCH --mail-user=mtello@student.ubc.ca
#SBATCH --array=1-20  # Adjust to the number of samples

# 1) Activate environment or load Salmon module
source activate fastqc_env

# 2) Define directories
GenomeDir="/scratch/st-bstefans-1/mtello/StewardCollab/mouse_transcriptome"
ReadsDir="/scratch/st-bstefans-1/mtello/StewardCollab/reads"
OutputDir="/scratch/st-bstefans-1/mtello/StewardCollab/Salmon_aligned"

# 3) Create list of sample prefixes (Option A: Hard-coded array)
samples=(
    "PF2-SJ-010824_S25_L001"
    "PF3-SJ-010824_S21_L001"
    "PF4-SJ-010824_S23_L001"
    "PF5-2-SJ-010824_S26_L001"
    "PF8-SJ-010824_S22_L001"
    "PF9-SJ-010824_S24_L001"
    "S10-SJ-301123_S12_L001"
    "S3-SJ-010824_S20_L001"
    "S4-SJ-010824_S16_L001"
    "S5-SJ-010824_S18_L001"
    "S6-SJ-010824_S17_L001"
    "S8-SJ-010824_S19_L001"
    "S8-SJ-010824_S2_L001"
    "V1-SJ-010824_S28_L001"
    "V10-SJ-010824_S15_L001"
    "V2-SJ-010824_S13_L001"
    "V3-SJ-010824_S29_L001"
    "V4-SJ-010824_S14_L001"
    "V7-SJ-010824_S27_L001"
    "V7-SJ-010824_S3_L001"
)

# 4) Grab the appropriate sample for this Slurm array task
sample="${samples[$SLURM_ARRAY_TASK_ID-1]}"

# 5) Build paths for R1 and R2
R1="${ReadsDir}/${sample}_R1_001.fastq"
R2="${ReadsDir}/${sample}_R2_001.fastq"

# (Optional) Check that the files exist
if [[ ! -f "$R1" ]]; then
  echo "Error: Cannot find R1 at $R1"
  exit 1
fi
if [[ ! -f "$R2" ]]; then
  echo "Error: Cannot find R2 at $R2"
  exit 1
fi

# 6) Make sure output directory exists
mkdir -p "$OutputDir"

# 7) Run Salmon
salmon quant \
    -i "${GenomeDir}/mmusculus_index" \
    -l IU \
    -1 "$R1" \
    -2 "$R2" \
    -p 16 \
    --validateMappings \
    -o "${OutputDir}/${sample}_quant"

echo "Finished processing $sample!"

#!/bin/bash
 
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --job-name=fastqc
#SBATCH --account=st-bstefans-1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtello@student.ubc.ca
 
################################################################################

#conda init
source activate fastqc_env


# Base directory containing the subdirectories with FASTQ files
BASE_DIR="/arc/project/st-bstefans-1/datasets/StewardCollab/R576-409024225/20240207-Illumina-mRNA-717883172"

# Iterate over each subdirectory
for dir in "$BASE_DIR"/*; do
  # Check if it's a directory
  if [ -d "$dir" ]; then
    # Find the R1 and R2 FASTQ files
    R1_FILE=$(find "$dir" -name "*_R1_001.fastq.gz")
    R2_FILE=$(find "$dir" -name "*_R2_001.fastq.gz")

    # Ensure both R1 and R2 files exist
    if [ -n "$R1_FILE" ] && [ -n "$R2_FILE" ]; then
      # Extract sample name from directory name
      SAMPLE_NAME=$(basename "$dir")

      # Define output file names
      OUT_R1="/scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/${SAMPLE_NAME}_R1_trimmed.fastq.gz"
      OUT_R2="/scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/${SAMPLE_NAME}_R2_trimmed.fastq.gz"

      # Run fastp
      echo "Processing $SAMPLE_NAME..."
      fastp -i "$R1_FILE" -I "$R2_FILE" -o "$OUT_R1" -O "$OUT_R2" \
        --disable_quality_filtering --detect_adapter_for_pe --low_complexity_filter --correction --trim_poly_g \
        --overrepresentation_analysis --dedup --dup_calc_accuracy 5

      echo "Finished processing $SAMPLE_NAME."
    else
      echo "Warning: Missing R1 or R2 files in $dir. Skipping..."
    fi
  fi
done
#!/bin/bash

#SBATCH --time=4:00:00                # Maximum runtime
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=16                   # Number of threads per task
#SBATCH --mem=64gb                    # Memory per task
#SBATCH --job-name=STAR_alignment     # Job name
#SBATCH --array=1-20                  # Job array index (one per sample)
#SBATCH --output=STAR_%A_%a.log       # Output log files per task
#SBATCH --account=st-bstefans-1
#SBATCH --mail-user=mtello@student.ubc.ca

################################################################################

source fastqc_env

# Directories
GenomeDir="/scratch/st-bstefans-1/mtello/StewardCollab/mouse_hm_genome"
TrimmedReadsDir="/scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads"
OutputDir="/scratch/st-bstefans-1/mtello/StewardCollab/STAR_aligned"

# Create output directory if it doesn't exist
mkdir -p "$OutputDir"

# List all R1 input files and store them in an array
R1_files=($TrimmedReadsDir/*_R1_trimmed.fastq)

# Get the file for the current job array task
R1=${R1_files[$SLURM_ARRAY_TASK_ID-1]}  # SLURM_ARRAY_TASK_ID starts at 1, array index starts at 0
R2=${R1/_R1_/_R2_}                     # Automatically get the corresponding R2 file
SampleName=$(basename "$R1" _R1_trimmed.fastq)  # Extract sample name

# Print the files being processed (for debugging)
echo "Processing sample: $SampleName"
echo "R1: $R1"
echo "R2: $R2"

# Run STAR alignment for the current sample
STAR --runThreadN 16 \
     --genomeDir "$GenomeDir" \
     --readFilesIn "$R1" "$R2" \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outFileNamePrefix "$OutputDir/${SampleName}"

echo "Finished processing sample: $SampleName"

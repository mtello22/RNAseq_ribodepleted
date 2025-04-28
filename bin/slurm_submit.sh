#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64gb
#SBATCH --job-name=RiboSeq_NF
#SBATCH --account=st-bstefans-1
#SBATCH --mail-user=mtello@student.ubc.ca
#SBATCH --mail-type=ALL

################################################################################

# Load modules needed for Nextflow
module load gcc nextflow

# Optional: activate conda base if needed for local tools
source ~/.bashrc
conda activate rnaseq_riboD_env

export NXF_HOME="/scratch/st-bstefans-1/mtello/.nextflow"
export NXF_PLUGINS_DIR="$NXF_HOME/plugins"

# Disable outbound telemetry ping
export NXF_DISABLE_TELEMETRY=true

# Set project directory
BaseDir="/scratch/st-bstefans-1/mtello/ribodepleted_rnaseq"

# Change into project folder
cd "$BaseDir"

# Run pipeline with Sockeye profile
nextflow run main.nf -profile sockeye -resume

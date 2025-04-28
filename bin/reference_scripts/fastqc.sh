#!/bin/bash
 
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --job-name=fastqc
#SBATCH --account=st-bstefans-1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtello@student.ubc.ca
 
################################################################################

module load perl 
module load openjdk

source activate fastqc_env

PERL5LIB=$PATH:~/perl5/lib/perl5/
export PERL5LIB

# INPUT_FILES=$(ls /arc/project/st-bstefans-1/datasets/AndrewCollab/R576-409024225/*/*/*.fastq.gz)
INPUT_FILES=$(ls /scratch/st-bstefans-1/mtello/StewardCollab/trimmed_reads/*.fastq.gz)

OUTPUT_DIR="/scratch/st-bstefans-1/mtello/StewardCollab/FASTQC"

fastqc $INPUT_FILES --noextract --outdir $OUTPUT_DIR

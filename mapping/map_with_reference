#!/bin/bash
#SBATCH --job-name=bowtie_align   # job name
#SBATCH --ntasks=1                # number of tasks
#SBATCH --cpus-per-task=8         # number of CPUs per task
#SBATCH --mem=16G                 # memory per node
#SBATCH --time=24:00:00           # time limit
#SBATCH --output=bowtie_align.%j.out  # output file name
#SBATCH --error=bowtie_align.%j.err   # error file name

#Purge any previous modules
module purge


# load Bowtie and SAMtools modules (modify as needed for your system)
module load bioinfo/samtools/1.21
module load bioinfo/bowtie/2.5.4


# set variables

reference_genome=/work/user/mcastanier/mapping/scripts/refDB.fasta
input_dir=/work/user/mcastanier/mapping/Seq_raw_read
output_dir=/work/user/mcastanier/mapping/Seq_raw_read/outputBow


BOWTIEBUILDBIN=bowtie2-build;
BOWTIEBIN=bowtie2;

# Build bowtie index
mkdir BWT_INDEX;
cd BWT_INDEX;
$BOWTIEBUILDBIN --quiet $reference_genome refDB;
cd ..;

# loop over input files
for input_file in $input_dir/*_R1_clean.fastq.gz; do
    # get input file name without extension
    input_base=$(basename $input_file _R1_clean.fastq.gz)
    # set output file prefix
    output_prefix=$input_base.aligned_consensus


# Run bowtie
$BOWTIEBIN --very-sensitive --no-unal -I 50 -X 700 --fr --threads 3 -x BWT_INDEX/refDB -1 $input_dir/${input_base}_R1_clean.fastq.gz -2 $input_dir/${input_base}_R2_clean.fastq.gz -S $output_dir/${input_base}_refDB.sam &> log_${input_base}_bowtieMap.log;

done

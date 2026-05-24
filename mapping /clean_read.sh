#!/bin/bash
#SBATCH --job-name=clean_reads
#SBATCH --output=clean_reads.out
#SBATCH --error=clean_reads.err
#SBATCH --partition=workq
#SBATCH --time=24:00:00         
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16

module purge
module load bioinfo/fastp/0.23.2

#################################
# Loop over all read files
#################################

for R1_IN in *_R1_001.fastq.gz
do
    # Extract sample prefix using basename
    READSPFXONE=$(basename "$R1_IN" _R1_001.fastq.gz)

    R2_IN=${READSPFXONE}_R2_001.fastq.gz

    echo "Processing sample: $READSPFXONE"

    # Output cleaned reads
    R1_CLEAN=${READSPFXONE}_R1_clean.fastq.gz
    R2_CLEAN=${READSPFXONE}_R2_clean.fastq.gz

    # fastp reports
    FASTP_HTML=${READSPFXONE}_fastp.html
    FASTP_JSON=${READSPFXONE}_fastp.json

    fastp \
        -i "$R1_IN" \
        -I "$R2_IN" \
        -o "$R1_CLEAN" \
        -O "$R2_CLEAN" \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --n_base_limit 5 \
        --detect_adapter_for_pe \
        --thread $SLURM_CPUS_PER_TASK \
        --html "$FASTP_HTML" \
        --json "$FASTP_JSON"

    echo "Finished: $READSPFXONE"
    echo "----------------------------------"

done

echo "All samples processed."

#!/bin/bash
#SBATCH  -p workq
#SBATCH --job-name=Robust_consensus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=RobustCons.%j.out
#SBATCH --error=RobustCons.%j.err

module load bioinfo/samtools/1.21
module load bioinfo/Bcftools/1.21

############################################
# CONVERT SAM TO BAM 
############################################

for SAM in *.sam;
do
if [ -f "$SAM" ]; then
BASENAME=$(basename ${SAM} .sam)
echo "Converting ${SAM} to BAM..."

samtools view -@ ${SLURM_CPUS_PER_TASK} -b \
${SAM} > ${BASENAME}.bam
fi
done


############################################
# PARAMETERS
############################################

REF=refDB.fasta
MIN_DP=5                     # position-level masking
CONTIG_MIN_MEAN_DP=1         # contig-level absence threshold
CONSENSUS_DIR=Robust_consensus2

mkdir -p ${CONSENSUS_DIR}

# Index reference if needed
if [ ! -f ${REF}.fai ]; then
    samtools faidx ${REF}
fi

############################################
# LOOP OVER BAM FILES
############################################

for FILE in *.bam;
do
    BASENAME=$(basename ${FILE} .bam)
    echo "=============================="
    echo "Processing ${BASENAME}"
    echo "=============================="

    ############################################
    # Sort + index BAM
    ############################################
    samtools sort -@ ${SLURM_CPUS_PER_TASK} \
        -o ${BASENAME}.sorted.bam \
        ${FILE}

    samtools index ${BASENAME}.sorted.bam

    ############################################
    # 2 Call ALL sites 
    ############################################
    bcftools mpileup \
        -Ou \
        -f ${REF} \
        -a FORMAT/DP \
        -q 15 -Q 15 \
        --max-depth 5000 \
        ${BASENAME}.sorted.bam | \
    bcftools call \
        -m --ploidy 1 -A \
        -Oz \
        -o ${BASENAME}.all.vcf.gz

    bcftools index ${BASENAME}.all.vcf.gz

    ############################################
    # 3 Compute mean depth per contig
    ############################################
    samtools depth -aa ${BASENAME}.sorted.bam | \
    awk '{
        sum[$1]+=$3
        count[$1]++
    }
    END {
        for (c in sum)
            print c"\t"sum[c]/count[c]
    }' > ${BASENAME}_contig_mean_depth.txt

    ############################################
    # Build mask
    ############################################
    > ${BASENAME}_mask.bed

    while read CONTIG MEAN_DP; do

        # If contig nearly absent ? mask entire contig
        if (( $(echo "$MEAN_DP < ${CONTIG_MIN_MEAN_DP}" | bc -l) )); then

            echo "Masking entire contig ${CONTIG} (mean DP = ${MEAN_DP})"

            LENGTH=$(samtools faidx ${REF} ${CONTIG} | \
                     grep -v ">" | tr -d '\n' | wc -c)

            echo -e "${CONTIG}\t0\t${LENGTH}" \
                >> ${BASENAME}_mask.bed

        else
            # Mask only low coverage positions
            samtools depth -aa -r ${CONTIG} ${BASENAME}.sorted.bam | \
            awk -v MIN_DP=${MIN_DP} \
                '$3 < MIN_DP {print $1"\t"$2-1"\t"$2}' \
                >> ${BASENAME}_mask.bed
        fi

    done < ${BASENAME}_contig_mean_depth.txt

    echo "Masked positions:"
    wc -l ${BASENAME}_mask.bed

    ############################################
    # Build consensus from  VCF
    ############################################
    bcftools consensus \
        -f ${REF} \
        -m ${BASENAME}_mask.bed \
        ${BASENAME}.all.vcf.gz \
        > ${CONSENSUS_DIR}/${BASENAME}_consensus.fa

done

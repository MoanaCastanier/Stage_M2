#!/bin/bash

# Script pour lancer Snippy sur tous les assemblages avec --rgid
# Usage: ./snippy_batch.sh dossier_fasta reference.fasta

FASTA_DIR="$1"
REFERENCE="$2"
OUTPUT_DIR="snippy_results_$FASTA_DIR"

REFERENCE_ABS=$(readlink -f "$REFERENCE")
mkdir -p "$OUTPUT_DIR"

# Boucle pour soumettre un job par échantillon
for fasta in "$FASTA_DIR"/*.{fasta,fa,fna}; do
    [ -f "$fasta" ] || continue
    
    fasta_abs=$(readlink -f "$fasta")
    sample=$(basename "$fasta" | sed 's/\.\(fasta\|fa\|fna\)$//')
    
    sbatch <<EOF
#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00 
#SBATCH --mem=16G   
#SBATCH --cpus-per-task=8
#SBATCH -J snippy_${sample}
#SBATCH -o ${OUTPUT_DIR}/${sample}_snippy.log

module load devel/Miniconda/Miniconda3
module load bioinfo/Snippy/4.6.0

snippy --cpus 8 \
       --outdir ${OUTPUT_DIR}/${sample}_snippy \
       --ref ${REFERENCE_ABS} \
       --ctgs ${fasta_abs} \
       --rgid ${sample}
EOF

done

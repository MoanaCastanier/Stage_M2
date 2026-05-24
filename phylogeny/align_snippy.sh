#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00 
#SBATCH --mem=32G   
#SBATCH --cpus-per-task=4
#SBATCH -J snippy_core
#SBATCH -o snippy_core.log

module load devel/Miniconda/Miniconda3
module load bioinfo/Snippy/4.6.0

SNIPPY_DIR="$1"
REFERENCE="$2"

snippy-core --ref "$REFERENCE" \
            --prefix "core_$SNIPPY_DIR" \
            $SNIPPY_DIR/*_snippy

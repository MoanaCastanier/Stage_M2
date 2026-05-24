#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00 
#SBATCH --mem=64G   
#SBATCH --cpus-per-task=16
#SBATCH -J phylogeny
#SBATCH -o phylogeny.log

module load devel/Miniconda/Miniconda3
#Need gcc-12.2.0
module load compilers/gcc/12.2.0
module load bioinfo/IQ-TREE/3.0.1

ALIGNMENT="$1"

# Construire l'arbre phylogénétique
iqtree -s "$ALIGNMENT" -m MFP -bb 1000 -nt AUTO


Script allowing the reconstruction of the phylogenies of the symbionts *Buchnera* and *Serratia* in 3 steps. All scripts are written in bash and made to be launched on a cluster thanks to the first lines (#SBATCH). Can be modified according to the needs and demands of the cluster.

Step 1: snippy_SNP_indiv.sh, SNP calling carried out using snippy for each individual in a file. SNPs are relative to a reference sequence provided to the script. 

Step 2: align_snippy.sh, snippy_core command that uses the SNP calling from step 1 to align all genomes with the reference sequence

Step 3: Reconstruction of phylogenetic trees by IQtree based on the complete alignment of align_snippy.sh 

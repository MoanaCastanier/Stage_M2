All scripts are written in bash and made to be launched on a cluster thanks to the first lines (#SBATCH). Can be modified according to the needs and demands of the cluster.

**Map with reference** : Uses Bowtie to map raw reads along selected lookup sequences.

**Make consensus** : uses Bcftools and Samtools to create a consensus sequence for each individual. Uses the output of map_with_reference.sh (BAM file). 

**Extract Serratia/Buchnera** : Use the output from make_consensus.sh to extract the consensus sequences for each of *Serratia* or *Buchnera*.

L'ensemble des scripts est écrit en bash et fait pour être lancer sur un cluster grâce aux premieres lignes (#SBATCH). Peut être modifier en fonction des besoins et des demandes du cluster. 

**Map with reference** : utilise Bowtie pour mapper les raw read le long de séquences de références choisie. 

**Extract Serratia/Buchnera** : utilise les sortie de map_with_reference.sh (fichier BAM) et extrait les sortie appartenant a Serratia et Buchnera.

**Make consensus** : utilise Bcftools et Samtools pour créer une séquence consensus pour chaque individu. Utilise les sortie de extract_Buchnera.sh et extract_Serratia.sh (fichier BAM specifique à un symbiote). 



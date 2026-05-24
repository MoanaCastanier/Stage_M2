#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=ExtrCons.%j.out
#SBATCH --error=ExtrCons.%j.err

mkdir SerratiaBcf

for fa in *.fa; do
    out="SerratiaBcf/Serratia_cons_${fa}"

    awk '
        /^>/ { keep = ($0 ~ /Serratia_chr/) }
        keep
    ' "$fa" > "$out"

done

#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=ExtrCons.%j.out
#SBATCH --error=ExtrCons.%j.err

mkdir BuchneraBcf

for fa in *.fa; do
    out="BuchneraBcf/Buchnera_cons_${fa}"

    awk '
        /^>/ { keep = ($0 ~ /Buchnera_aphidicola_str_5A/) }
        keep
    ' "$fa" > "$out"

done

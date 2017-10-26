#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p normal
#SBATCH -J imcrblst
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user kyclark@email.arizona.edu

OUT_DIR="$SCRATCH/imicrobe-blast/test"

[[ -d "$OUT_DIR" ]] && rm -rf "$OUT_DIR"

run.sh -q "$WORK/data/blast/pov.fa" -q "$WORK/data/blast/dolphin.fa" -o $OUT_DIR

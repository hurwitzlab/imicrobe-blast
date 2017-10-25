#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -p development
#SBATCH -J imcrblst
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user kyclark@email.arizona.edu

OUT_DIR="$SCRATCH/imicrobe-blast/test"

[[ -d "$OUT_DIR" ]] && rm -rf "$OUT_DIR"

run.sh -q "$WORK/data/blast/test1.fa" -q "$WORK/data/blast/test2.fa" -o $OUT_DIR

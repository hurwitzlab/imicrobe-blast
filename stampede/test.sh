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

if [[ -d $OUT_DIR ]]; then
  rm -rf $OUT_DIR
fi

run.sh -q "$SCRATCH/imicrobe-blast/test.fa" -o $OUT_DIR

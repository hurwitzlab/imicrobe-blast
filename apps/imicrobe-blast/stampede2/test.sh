#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p skx-normal
#SBATCH -J imcrblst
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user jklynch@email.arizona.edu

OUT_DIR="$SCRATCH/imicrobe-blast/test"

[[ -d "$OUT_DIR" ]] && rm -rf "$OUT_DIR"

DATA_DIR=/work/05066/imicrobe/iplantc.org/data/imicrobe/projects/37/samples

# these data files are under 1MB

sh run.sh -q "$DATA_DIR/739/ICOMM_SMPL_LTERSTATION1_REPLICATION1.fa" -q "$DATA_DIR/740/ICOMM_SMPL_LTERSTATION1_REPLICATION2.fa" -o $OUT_DIR

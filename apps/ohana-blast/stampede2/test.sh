#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -p skx-normal
#SBATCH -J test-ohana-blast
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user jklynch@email.arizona.edu

OUT_DIR="$SCRATCH/ohana-blast/test"

[[ -d "$OUT_DIR" ]] && rm -rf "$OUT_DIR"

# these files are under 1MB
DATA_DIR=/work/05066/imicrobe/iplantc.org/data/imicrobe/projects/37/samples
sh run.sh -q "$DATA_DIR/739/ICOMM_SMPL_LTERSTATION1_REPLICATION1.fa" -q "$DATA_DIR/740/ICOMM_SMPL_LTERSTATION1_REPLICATION2.fa" -o $OUT_DIR

# these files are about 31MB and 22MB
#DATA_DIR=/work/05066/imicrobe/iplantc.org/data/imicrobe/projects/94/samples
#sh run.sh -q "$DATA_DIR/1570/CAM_SMPL_002229.fa" -q "$DATA_DIR/1569/CAM_SMPL_002228.fa" -o $OUT_DIR

# these files are over 1.6GB
#DATA_DIR=/work/05066/imicrobe/iplantc.org/data/imicrobe/projects/263/samples/
#sh run.sh -q "$DATA_DIR/5199/HOT227_1_0500m-readpool.fa" -q "$DATA_DIR/5197/HOT227_1_0075m-readpool.fa" -o "$OUT_DIR/large"

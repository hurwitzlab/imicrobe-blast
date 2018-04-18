#!/usr/bin/env bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -p skx-normal
#SBATCH -J pack-files-ohana-contigs
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user jklynch@email.arizona.edu


module load tacc-singularity

IMG=$WORK/jklynch/project/imicrobe/apps/imicrobe-blast/singularity/imicrobe-blast-0.0.5.img
OHANA_HOT_DATA=/work/05066/imicrobe/iplantc.org/data/ohana/HOT

export BLAST_DB_OUTPUT=/work/05066/imicrobe/iplantc.org/data/ohana/blast

singularity exec $IMG pack_file_paths \
    --split-count 24 \
    --db-uri sqlite:///${BLAST_DB_OUTPUT}/ohana_contigs_seq_db.sqlite \
    --prefix ohana-contigs- \
    --max-workers 48

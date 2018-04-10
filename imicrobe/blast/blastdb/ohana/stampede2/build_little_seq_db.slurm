#!/usr/bin/env bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -p skx-normal
#SBATCH -J make-ohana-seq-db
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user jklynch@email.arizona.edu


module load tacc-singularity

IMG=$WORK/jklynch/project/imicrobe/apps/imicrobe-blast/singularity/imicrobe-blast-0.0.5.img
OHANA_HOT_DATA=/work/05066/imicrobe/iplantc.org/data/ohana/HOT

export BLAST_DB_OUTPUT=/work/05066/imicrobe/iplantc.org/data/ohana/blast

# delete the old database if there is one
#rm -f ${BLAST_DB_OUTPUT}/ohana_little_seq_db.sqlite

singularity exec $IMG build_seq_db \
    --fasta-globs ${OHANA_HOT_DATA}/HOT_224_*/**/contigs.fa \
    --db-uri sqlite:///${BLAST_DB_OUTPUT}/ohana_little_seq_db.sqlite \
    --valid-files-fp ohana_valid_fasta_files.txt \
    --invalid-files-fp ohana_invalid_fasta_files.txt \
    --max-workers 48

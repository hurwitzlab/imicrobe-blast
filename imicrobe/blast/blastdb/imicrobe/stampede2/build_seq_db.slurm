#!/usr/bin/env bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -p skx-normal
#SBATCH -J make-imicrobe-seq-db
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user joshua.kevin.lynch@gmail.com

#
# This script is intended to be run more than once since it has been taking
# more than 48 hours to build the sequence database for iMicrobe.
# No changes are required. Just let it run and then run it again until it
# has nothing to do.
#

module load tacc-singularity

IMG=$WORK/jklynch/project/imicrobe/apps/imicrobe-blast/singularity/imicrobe-blast-0.0.5.img
IMICROBE_DATA=/work/05066/imicrobe/iplantc.org/data/imicrobe/projects

export BLAST_DB_OUTPUT=/work/05066/imicrobe/iplantc.org/data/imicrobe/blast

rm ${BLAST_DB_OUTPUT}/imicrobe_seq_db.sqlite
rm imicrobe_valid_fasta_files.txt
rm imicrobe_invalid_fasta_files.txt

# running out of memory with 48 workers
# running out of memory with 36 workers
singularity exec $IMG build_seq_db \
    --fasta-globs ${IMICROBE_DATA}/**/*.fa,${IMICROBE_DATA}/**/*.fna,${IMICROBE_DATA}/**/*.fasta \
    --db-uri sqlite:///${BLAST_DB_OUTPUT}/imicrobe_seq_db.sqlite \
    --valid-files-fp imicrobe_valid_fasta_files.txt \
    --invalid-files-fp imicrobe_invalid_fasta_files.txt \
    --max-workers 24

#!/usr/bin/env bash
#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -p skx-normal
#SBATCH -J make-blast-db
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user jklynch@email.arizona.edu

module load blast
module load tacc-singularity

module use /scratch/01255/siliu/modulefiles
module load launcher/3.2

IMG=../../../../../singularity/imicrobe-blast-0.0.5.img
IMICROBE_PROJECTS=/work/05066/imicrobe/iplantc.org/data/imicrobe/projects
FILE_GROUPS=24

export BLAST_DB_OUTPUT=/work/05066/imicrobe/iplantc.org/data/blast/db24

mkdir -p $BLAST_DB_OUTPUT

find ${IMICROBE_PROJECTS} \
    -type f \
    -regextype posix-egrep \
    -regex ".+\.(fa|fna|fasta)$" \
    -size +0c \
    | sort \
    | uniq \
    | singularity exec ${IMG} filter_lines bad-fasta-files.txt \
    | singularity exec ${IMG} pack_file_paths -n ${FILE_GROUPS} --prefix imicrobe-

#
# go!
#

# no need to change these
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM

# must define LAUNCHER_WORKDIR somewhere
export LAUNCHER_WORKDIR=`pwd`
echo "LAUNCHER_WORKDIR: $LAUNCHER_WORKDIR"

# by default LAUNCHER_PPN is 8 (?)
# use LAUNCHER_PPN to control the number of jobs on each node
export LAUNCHER_PPN=${FILE_GROUPS}

export LAUNCHER_SCHEDULER=dynamic

export LAUNCHER_JOB_FILE=jobs${FILE_GROUPS}.launcher
$LAUNCHER_DIR/paramrun

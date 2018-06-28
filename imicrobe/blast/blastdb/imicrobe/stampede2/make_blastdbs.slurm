#!/usr/bin/env bash
#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -p skx-normal
#SBATCH -J make-imicrobe-blast-db
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user joshua.kevin.lynch@gmail.com

module load blast
module load tacc-singularity

IMG=/work/05066/imicrobe/singularity/imicrobe-blast-0.0.5.img
IMICROBE_PROJECTS=/work/05066/imicrobe/iplantc.org/data/imicrobe/projects
FILE_GROUPS=24

export BLAST_DB_OUTPUT=/work/05066/imicrobe/iplantc.org/data/blast/db24

mkdir -p $BLAST_DB_OUTPUT

#find ${IMICROBE_PROJECTS} \
#    -type f \
#    -regextype posix-egrep \
#    -regex ".+\.(fa|fna|fasta)$" \
#    -size +0c \
#    | sort \
#    | uniq \
#    | singularity exec ${IMG} filter_lines bad-fasta-files.txt \
#    | singularity exec ${IMG} pack_file_paths -n ${FILE_GROUPS} --prefix imicrobe-

cat imicrobe_valid_fasta_files.txt | singularity exec ${IMG} pack_file_paths -n ${FILE_GROUPS} --prefix imicrobe-

#
# go!
#

# assume parallel is already on $PATH
# define the env_parallel function (important!)
. $(which env_parallel.bash)

# nodelist
scontrol show hostnames $SLURM_JOB_NODELIST > nodelist

# $JOBFILE is the list of jobs
JOBFILE=makeblastdb_24.jobs

# PPN is processes per node
PPN=24

# env_parallel sends the shell environment to each process it starts

env_parallel --max-procs ${PPN} --joblog paralleljobs.log  --sshloginfile nodelist --group < ${JOBFILE}

echo " "
echo "GNU Parallel Job Complete"
echo " "